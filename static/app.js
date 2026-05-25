// PGCP webapp — frontend.
//
// Loads the preset list, lets the user pick a region, displays sensors on a
// Leaflet/OSM map, and talks to the Julia backend for live R recomputation.

const state = {
  regions: [],
  current: null,
  pk: 0.9,
  disabled: new Set(),
  map: null,
  layers: {
    polygon: null,
    sensors: null,
  },
  sensorMarkers: new Map(),
  sensorDisks: new Map(),
  evaluatePending: false,
  evaluateQueued: false,

  editSensors: [],
  originalSensors: [],
  dirty: false,
  nextNewIndex: 1,
  replanInFlight: false,
  newRadius: 0.089,

  // Popup-open guard: prevents an accidental "add sensor" when the user
  // interacts with a popup that overlaps the map (slider, delete button, etc).
  popupOpenCount: 0,
  popupClosedAt: 0,
};

const RADIUS_SMALL = 0.089;
const RADIUS_LARGE = 0.178;

// Sensor dot radius: scaled to current Birnbaum importance.
// At p_k = 0.9 the maximum possible importance is R/p_k ≈ 0.497 for an
// essential sensor; using 0.5 as the normalization keeps the scale stable
// across the slider range. Inactive (low-importance) sensors shrink to a
// readable floor of 3 px.
const DOT_MIN_PX = 3;
const DOT_MAX_PX = 12;
const IMPORTANCE_REF = 0.5;
let currentImportance = new Map(); // sensor name -> Birnbaum importance
let currentR = 0;                  // last R value returned by /evaluate
let currentEssential = new Set();  // sensors whose disable would alone tank R now

// A sensor is "currently essential" iff its Birnbaum importance equals R/p_k,
// i.e. disabling it alone reduces R to 0 from the current state. The
// equivalent characterization (size-1 min-cut of the residual system) is what
// the orange ring should highlight. We use a relative tolerance because the
// underlying probabilities are products of p_k values and accumulate FP error.
function computeCurrentEssential(R, pk) {
  const out = new Set();
  if (R < 1e-12 || pk < 1e-12) return out;
  const target = R / pk;
  for (const [name, imp] of currentImportance) {
    if (state.disabled.has(name)) continue;
    if (target > 0 && Math.abs(imp - target) / target < 1e-6) {
      out.add(name);
    }
  }
  return out;
}

// ---------------- API helpers ----------------

async function api(path, body) {
  const opts = { headers: { "Content-Type": "application/json" } };
  if (body !== undefined) {
    opts.method = "POST";
    opts.body = JSON.stringify(body);
  }
  const resp = await fetch(path, opts);
  if (!resp.ok) {
    throw new Error(`${path} -> ${resp.status} ${resp.statusText}`);
  }
  return await resp.json();
}

// ---------------- Coordinate helpers ----------------

// Convert normalized (x, y) in [0,1] x [0,h] to (lat, lon) via affine bbox map.
function normToLatLon(x, y, region) {
  const bn = region.bbox_normalized; // [x0, y0, x1, y1]
  const bl = region.bbox_latlon;     // {min_lat, max_lat, min_lon, max_lon}
  const u = (x - bn[0]) / (bn[2] - bn[0]);
  const v = (y - bn[1]) / (bn[3] - bn[1]);
  const lon = bl.min_lon + u * (bl.max_lon - bl.min_lon);
  const lat = bl.min_lat + v * (bl.max_lat - bl.min_lat);
  return [lat, lon];
}

function radiusNormToMeters(rNorm, region) {
  return rNorm * region.side_meters;
}

// ---------------- Rendering ----------------

function clearLayers() {
  if (state.layers.polygon) { state.map.removeLayer(state.layers.polygon); state.layers.polygon = null; }
  if (state.layers.sensors) { state.map.removeLayer(state.layers.sensors); state.layers.sensors = null; }
  state.sensorMarkers.clear();
  state.sensorDisks.clear();
}

// ---------------- Coordinate inverse ----------------

function latLonToNorm(lat, lon, region) {
  const bn = region.bbox_normalized;
  const bl = region.bbox_latlon;
  const u = (lon - bl.min_lon) / (bl.max_lon - bl.min_lon);
  const v = (lat - bl.min_lat) / (bl.max_lat - bl.min_lat);
  return [
    bn[0] + u * (bn[2] - bn[0]),
    bn[1] + v * (bn[3] - bn[1]),
  ];
}

function renderRegion(region) {
  clearLayers();
  const bl = region.bbox_latlon;
  const bounds = [[bl.min_lat, bl.min_lon], [bl.max_lat, bl.max_lon]];
  state.map.fitBounds(bounds, { padding: [20, 20] });

  // Polygon outline.
  state.layers.polygon = L.polygon(region.polygon_latlon, {
    color: "#333",
    weight: 2,
    fillColor: "#3a3",
    fillOpacity: 0.06,
    interactive: false,
  }).addTo(state.map);

  // Sensors: a faint coverage disk + a clickable marker.
  state.layers.sensors = L.layerGroup().addTo(state.map);
  const essential = new Set(region.precomputed.essential_sensors || []);
  const importanceByName = new Map(
    (region.precomputed.importance || []).map(r => [r.name, r.birnbaum])
  );
  const maxImp = Math.max(0.0001, ...importanceByName.values());

  // Seed importance with the precomputed (default-state) values; they will be
  // overwritten by the first /evaluate response.
  currentImportance = new Map(
    (region.precomputed.importance || []).map(r => [r.name, r.birnbaum])
  );

  for (const c of region.sensors) {
    addSensorToMap(c, region);
  }
  refreshSensorStyles();
}

function addSensorToMap(c, region) {
  const [lat, lon] = normToLatLon(c.x, c.y, region);
  const radiusMeters = radiusNormToMeters(c.radius, region);

  const disk = L.circle([lat, lon], {
    radius: radiusMeters,
    color: "#1f6feb",
    weight: 1,
    opacity: 0.35,
    fillColor: "#1f6feb",
    fillOpacity: 0.05,
    interactive: false,
  });
  disk.addTo(state.layers.sensors);
  state.sensorDisks.set(c.name, disk);

  const marker = L.circleMarker([lat, lon], {
    radius: DOT_MIN_PX,
    color: "#1f6feb",
    weight: 2,
    fillColor: "#1f6feb",
    fillOpacity: 0.85,
  });
  marker.bindTooltip(buildTooltip(c, region, false), { sticky: true });
  marker.on("click", (e) => {
    // L.CircleMarker does NOT auto-stop map click propagation (unlike L.Marker).
    // We use three layered defenses:
    //   1. native stopPropagation so the click doesn't bubble through the DOM
    //   2. Leaflet's internal `_stopped` flag so its `_fireDOMEvent` exits early
    //   3. (in onMapClick) ignoring clicks whose target is `.leaflet-interactive`
    if (e.originalEvent) {
      e.originalEvent.stopPropagation();
      e.originalEvent._stopped = true;
    }
    onSensorClick(c.name, e);
  });
  attachDragHandlers(marker, c.name);
  marker.addTo(state.layers.sensors);
  state.sensorMarkers.set(c.name, marker);
}

function removeSensorFromMap(name) {
  const m = state.sensorMarkers.get(name);
  const d = state.sensorDisks.get(name);
  if (m) { state.layers.sensors.removeLayer(m); state.sensorMarkers.delete(name); }
  if (d) { state.layers.sensors.removeLayer(d); state.sensorDisks.delete(name); }
}

function buildTooltip(c, region, isEssential) {
  const imp = currentImportance.get(c.name) || 0;
  return `<div class="sensor-tooltip">` +
    `<strong>センサ ${c.name}</strong><br>` +
    `半径 ${(c.radius * region.side_meters).toFixed(0)} m<br>` +
    `Birnbaum 重要度: ${imp.toFixed(4)}` +
    (isEssential ? "<br><em>このセンサ単独で被覆崩壊</em>" : "") +
  `</div>`;
}

function importanceToPixels(imp) {
  const t = Math.max(0, Math.min(1, imp / IMPORTANCE_REF));
  return DOT_MIN_PX + (DOT_MAX_PX - DOT_MIN_PX) * t;
}

function refreshSensorStyles() {
  if (!state.current) return;
  for (const c of state.current.sensors) {
    const marker = state.sensorMarkers.get(c.name);
    const disk = state.sensorDisks.get(c.name);
    if (!marker) continue;
    const isDisabled = state.disabled.has(c.name);
    const isEssential = currentEssential.has(c.name);
    const imp = currentImportance.get(c.name) || 0;
    const pixR = isDisabled ? DOT_MIN_PX : importanceToPixels(imp);
    marker.setRadius(pixR);
    marker.setStyle({
      color: isDisabled ? "#777" : (isEssential ? "#d99c00" : "#1f6feb"),
      fillColor: isDisabled ? "#bbb" : "#1f6feb",
      fillOpacity: isDisabled ? 0.4 : 0.85,
      weight: isEssential ? 3 : 2,
    });
    marker.setTooltipContent(buildTooltip(c, state.current, isEssential));
    if (disk) {
      disk.setStyle({
        opacity: isDisabled ? 0.08 : 0.35,
        fillOpacity: isDisabled ? 0.0 : 0.05,
      });
    }
  }
}

function renderImportanceFromMap() {
  const container = document.getElementById("importance-bars");
  container.innerHTML = "";
  const items = Array.from(currentImportance.entries())
    .map(([name, v]) => ({ name, birnbaum: v }))
    .filter(it => it.birnbaum > 0)
    .sort((a, b) => b.birnbaum - a.birnbaum)
    .slice(0, 10);
  if (items.length === 0) {
    container.innerHTML = `<div class="muted small">この故障状態では、どのセンサも信頼性に寄与しません (R = 0)。</div>`;
    return;
  }
  const maxV = Math.max(0.0001, ...items.map(i => i.birnbaum));
  for (const it of items) {
    const row = document.createElement("div");
    row.className = "imp-row";
    row.innerHTML =
      `<span class="imp-name">#${it.name}</span>` +
      `<span class="imp-bar"><span class="imp-bar-fill" style="width: ${(it.birnbaum / maxV * 100).toFixed(1)}%"></span></span>` +
      `<span class="imp-val">${it.birnbaum.toFixed(3)}</span>`;
    container.appendChild(row);
  }
}

function renderFailureList() {
  const div = document.getElementById("failure-list");
  if (state.disabled.size === 0) {
    div.textContent = "センサをクリックすると故障させられます。";
    return;
  }
  const names = Array.from(state.disabled).sort((a, b) => Number(a) - Number(b));
  div.innerHTML = `故障中: <strong>${names.length}</strong> 個 ― ${names.map(n => `#${n}`).join(", ")}`;
}

// ---------------- Interactions ----------------

function onSensorClick(name, event) {
  if (state.replanInFlight) return;
  // Shift+click opens the per-sensor popup (radius slider, delete).
  // A bare click just toggles failure — the primary demo action.
  if (event?.originalEvent?.shiftKey) {
    openSensorPopup(name, event);
    return;
  }
  toggleFailure(name);
}

// Radius bounds for the per-sensor slider. The smallest is well below the
// PDS minimum used in the paper (200 m), and the largest is roughly the
// 400-m sensor that already appears in the campus dataset; this range is
// chosen to be expressive without leading to a triangulation that the
// solver cannot finish quickly.
const RADIUS_MIN_NORM = 0.010;
const RADIUS_MAX_NORM = 0.300;

function openSensorPopup(name, event) {
  const sensor = state.editSensors.find(s => s.name === name);
  if (!sensor) return;
  const marker = state.sensorMarkers.get(name);
  if (!marker) return;

  const sideM = state.current.side_meters;
  const html = `
    <div class="sensor-popup">
      <div><strong>センサ ${sensor.name}</strong></div>
      <div class="muted small">
        半径 <span data-role="meters">${Math.round(sensor.radius * sideM)}</span> m
        <span class="muted">(正規化 <span data-role="norm">${sensor.radius.toFixed(3)}</span>)</span>
      </div>
      <input type="range" class="radius-slider"
             data-role="radius"
             min="${RADIUS_MIN_NORM}" max="${RADIUS_MAX_NORM}"
             step="0.001" value="${sensor.radius}">
      <div class="sensor-popup-actions">
        <button class="danger" data-action="delete">削除</button>
      </div>
    </div>
  `;

  marker.unbindPopup();
  marker.bindPopup(html, { autoClose: true, closeButton: true });
  marker.openPopup();

  setTimeout(() => {
    const popupEl = marker.getPopup()?.getElement();
    if (!popupEl) return;

    // Stop click / mousedown / mouseup / dblclick from bubbling into the map.
    // Leaflet's default disableClickPropagation covers mousedown/touchstart but
    // NOT click, so without this, a click on the slider would fire onMapClick
    // and add a duplicate sensor next to the popup. (Bug observed in v0.1.)
    L.DomEvent.on(
      popupEl,
      "click dblclick mousedown mouseup touchstart touchend",
      L.DomEvent.stopPropagation
    );

    const slider = popupEl.querySelector("[data-role='radius']");
    const metersEl = popupEl.querySelector("[data-role='meters']");
    const normEl = popupEl.querySelector("[data-role='norm']");
    const disk = state.sensorDisks.get(name);

    slider?.addEventListener("input", () => {
      const r = parseFloat(slider.value);
      sensor.radius = r;
      if (metersEl) metersEl.textContent = Math.round(r * sideM);
      if (normEl) normEl.textContent = r.toFixed(3);
      if (disk) disk.setRadius(radiusNormToMeters(r, state.current));
      markDirty();
    });

    popupEl.querySelector("[data-action='delete']")?.addEventListener("click", () => {
      const idx = state.editSensors.findIndex(s => s.name === name);
      if (idx >= 0) {
        state.editSensors.splice(idx, 1);
        removeSensorFromMap(name);
        markDirty();
        marker.closePopup();
      }
    });
  }, 0);
}

function toggleFailure(name) {
  if (state.disabled.has(name)) state.disabled.delete(name);
  else state.disabled.add(name);
  refreshSensorStyles();
  renderFailureList();
  scheduleEvaluate();
}

// ---- Placement mode: drag-to-move via manual mouse tracking on CircleMarker.

function attachDragHandlers(marker, name) {
  let dragging = false;
  let downAt = null;
  marker.on("mousedown", (e) => {
    if (state.replanInFlight) return;
    // Only treat as a drag if the user actually moves the cursor; a bare
    // mousedown-up at the same spot should fall through to the click handler
    // (toggle failure). We arm the listener now but skip committing in onUp
    // when no motion was seen.
    dragging = true;
    downAt = e.containerPoint;
    state.map.dragging.disable();
    state.map.on("mousemove", onMove);
    document.addEventListener("mouseup", onUp, { once: true });
  });
  let moved = false;
  function onMove(e) {
    if (!dragging) return;
    moved = true;
    marker.setLatLng(e.latlng);
    const disk = state.sensorDisks.get(name);
    if (disk) disk.setLatLng(e.latlng);
  }
  function onUp() {
    if (!dragging) return;
    dragging = false;
    state.map.off("mousemove", onMove);
    state.map.dragging.enable();
    if (!moved) { moved = false; return; }  // bare click — leave for click handler
    moved = false;
    const ll = marker.getLatLng();
    const [nx, ny] = latLonToNorm(ll.lat, ll.lng, state.current);
    const sensor = state.editSensors.find(s => s.name === name);
    if (sensor) { sensor.x = nx; sensor.y = ny; }
    markDirty();
  }
}

function onMapClick(e) {
  if (state.replanInFlight) return;
  // (1) If any popup is open, the user is interacting with it — never add.
  if (state.popupOpenCount > 0) return;
  // (2) If a popup was just closed (within the last 250 ms), this click is
  // the close action, not a deliberate empty-map click — skip.
  if (Date.now() - state.popupClosedAt < 250) return;
  // (3) Belt and suspenders: target inside any popup means we should skip.
  if (e.originalEvent?.target?.closest?.(".leaflet-popup")) return;
  // (4) Click landed on an existing sensor marker. The polygon outline and
  // coverage disks are created with `interactive: false`, so the only
  // `.leaflet-interactive` SVG path under the cursor is a sensor.
  if (e.originalEvent?.target?.classList?.contains("leaflet-interactive")) return;
  const [nx, ny] = latLonToNorm(e.latlng.lat, e.latlng.lng, state.current);
  if (nx < 0 || nx > 1 || ny < 0 || ny > 1) return;
  const name = `+${state.nextNewIndex++}`;
  const sensor = { name, x: nx, y: ny, radius: state.newRadius };
  state.editSensors.push(sensor);
  addSensorToMap(sensor, state.current);
  refreshSensorStyles();
  markDirty();
}

function clearAllSensors() {
  if (state.replanInFlight) return;
  for (const s of state.editSensors) removeSensorFromMap(s.name);
  state.editSensors = [];
  markDirty();
}

function markDirty() {
  state.dirty = true;
  document.body.classList.add("dirty");
  const el = document.getElementById("replan-status");
  el.classList.remove("running");
  el.classList.add("dirty");
  el.textContent = `配置を変更しました。「再計算」を押してください。`;
  refreshSensorCount();
}

function clearDirty() {
  state.dirty = false;
  document.body.classList.remove("dirty");
  refreshSensorCount();
}

function refreshSensorCount() {
  const el = document.getElementById("sensor-count");
  if (!el) return;
  el.textContent = state.editSensors.length;
}

async function doReplan(sensors) {
  state.replanInFlight = true;
  const status = document.getElementById("replan-status");
  status.classList.remove("dirty");
  status.classList.add("running");
  status.textContent = "BDD を再構築中…";
  try {
    const res = await api(`/api/regions/${state.current.id}/replan`, {
      sensors,
      pk: state.pk,
      disabled: Array.from(state.disabled).filter(n => sensors.some(s => s.name === n)),
    });
    // Merge new state into state.current.
    state.current.sensors = res.sensors;
    state.current.triangulation = res.triangulation;
    state.current.solver = { convergence: res.convergence, phi_nodes: null };
    state.current.n_sensors = res.n_sensors;
    // Replace stale precomputed min-path / min-cut so the quick-action buttons
    // in analysis mode refer to the NEW BDD's minimal sets.
    if (res.min_cut_sets)  state.current.precomputed.min_cut_sets  = res.min_cut_sets;
    if (res.min_path_sets) state.current.precomputed.min_path_sets = res.min_path_sets;
    // Rebuild the map and re-evaluate.
    state.editSensors = sensors.map(s => ({ ...s }));
    state.disabled = new Set(Array.from(state.disabled).filter(n => sensors.some(s => s.name === n)));
    renderRegion(state.current);
    currentImportance = new Map(Object.entries(res.importance));
    currentR = res.R;
    currentEssential = computeCurrentEssential(res.R, res.pk);
    document.getElementById("R-value").textContent = formatR(res.R);
    document.getElementById("R-hint").textContent =
      `(${res.n_sensors} センサ, ${res.elapsed_ms.toFixed(0)} ms で再計算)`;
    refreshSensorStyles();
    renderImportanceFromMap();
    renderFailureList();
    clearDirty();
    status.classList.remove("running");
    status.textContent = `再計算完了。R = ${formatR(res.R)} (${res.convergence ? "収束" : "未収束"})`;
  } catch (err) {
    console.error(err);
    status.classList.remove("running");
    status.classList.add("dirty");
    status.textContent = "再計算に失敗しました: " + err.message;
  } finally {
    state.replanInFlight = false;
  }
}

async function scheduleEvaluate() {
  if (!state.current) return;
  if (state.evaluatePending) { state.evaluateQueued = true; return; }
  state.evaluatePending = true;
  try {
    const res = await api(`/api/regions/${state.current.id}/evaluate`, {
      pk: state.pk,
      disabled: Array.from(state.disabled),
    });
    document.getElementById("R-value").textContent = formatR(res.R);
    document.getElementById("R-hint").textContent =
      `(p_k = ${res.pk.toFixed(3)}, 故障 ${res.disabled.length} 個, 評価 ${res.elapsed_ms.toFixed(1)} ms)`;
    if (res.importance) {
      currentImportance = new Map(Object.entries(res.importance));
      currentR = res.R;
      currentEssential = computeCurrentEssential(res.R, res.pk);
      refreshSensorStyles();
      renderImportanceFromMap();
    }
  } catch (err) {
    console.error(err);
    document.getElementById("R-value").textContent = "?";
    document.getElementById("R-hint").textContent = "サーバ通信に失敗しました。";
  } finally {
    state.evaluatePending = false;
    if (state.evaluateQueued) {
      state.evaluateQueued = false;
      scheduleEvaluate();
    }
  }
}

function formatR(R) {
  if (R >= 0.01) return R.toFixed(4);
  if (R > 0)     return R.toExponential(2);
  return "0";
}

function applyRegionPayload(region) {
  state.current = region;
  state.disabled = new Set();
  state.pk = region.default_pk;
  state.originalSensors = (region.original_sensors || region.sensors).map(s => ({ ...s }));
  state.editSensors = region.sensors.map(s => ({ ...s }));
  clearDirty();
  state.nextNewIndex = 1;
  refreshSensorCount();
  // Sync the "new sensor radius" meters readout to this region's scale.
  const meters = Math.round(state.newRadius * region.side_meters);
  const mEl = document.getElementById("new-radius-meters");
  if (mEl) mEl.textContent = meters;
  currentR = region.precomputed.R;
  currentEssential = new Set(region.precomputed.essential_sensors || []);
  document.getElementById("pk-slider").value = String(region.default_pk);
  document.getElementById("pk-value").textContent = region.default_pk.toFixed(3);
  document.getElementById("region-info").textContent =
    `${region.subtitle || ""} — ${region.n_sensors} センサ / 収束: ${region.solver.convergence ? "✓" : "未収束"}`;
  document.getElementById("R-value").textContent = formatR(region.precomputed.R);
  document.getElementById("R-hint").textContent =
    `(初期値 p_k = ${region.default_pk.toFixed(2)})`;
  renderRegion(region);
  renderImportanceFromMap();
  renderFailureList();
  scheduleEvaluate();
}

// ---------------- Map setup ----------------

function initMap() {
  state.map = L.map("map", { zoomControl: true });
  L.tileLayer("https://tile.openstreetmap.org/{z}/{x}/{y}.png", {
    attribution: '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a>',
    maxZoom: 19,
  }).addTo(state.map);
  state.map.setView([34.4, 132.7], 14);
  state.map.on("click", onMapClick);
  // Track popup open/close so onMapClick can ignore clicks happening through
  // an open popup or just after one closes.
  state.map.on("popupopen",  () => { state.popupOpenCount++; });
  state.map.on("popupclose", () => {
    state.popupOpenCount = Math.max(0, state.popupOpenCount - 1);
    state.popupClosedAt = Date.now();
  });
}

// ---------------- Boot ----------------

async function selectRegion(id) {
  const region = await api(`/api/regions/${encodeURIComponent(id)}`);
  applyRegionPayload(region);
}

async function boot() {
  initMap();
  const regions = await api("/api/regions");
  state.regions = regions;
  const sel = document.getElementById("region-select");
  sel.innerHTML = "";
  for (const r of regions) {
    const opt = document.createElement("option");
    opt.value = r.id;
    opt.textContent = r.display_name;
    sel.appendChild(opt);
  }
  sel.addEventListener("change", () => selectRegion(sel.value));

  document.getElementById("pk-slider").addEventListener("input", (e) => {
    state.pk = parseFloat(e.target.value);
    document.getElementById("pk-value").textContent = state.pk.toFixed(3);
    scheduleEvaluate();
  });

  document.getElementById("reset-failures").addEventListener("click", () => {
    state.disabled.clear();
    refreshSensorStyles();
    renderFailureList();
    scheduleEvaluate();
  });

  // Min cut: pick a small min-cut and DISABLE those sensors → R drops (often
  // all the way to 0 if the chosen cut is a singleton essential).
  document.getElementById("quick-cut").addEventListener("click", () => {
    if (!state.current) return;
    const smallest = state.current.precomputed.min_cut_sets?.smallest || [];
    if (!smallest.length) return;
    const pick = smallest[Math.floor(Math.random() * smallest.length)];
    state.disabled = new Set(pick);
    refreshSensorStyles();
    renderFailureList();
    scheduleEvaluate();
  });

  // Min path: pick a min-path set (sensors that JOINTLY suffice for coverage)
  // and disable every sensor NOT in that set → only the chosen path is alive.
  // R then equals p_k ^ |path| under independence, which the slider can
  // demonstrate live by varying p_k.
  document.getElementById("quick-path").addEventListener("click", () => {
    if (!state.current) return;
    const smallest = state.current.precomputed.min_path_sets?.smallest || [];
    if (!smallest.length) return;
    const pick = smallest[Math.floor(Math.random() * smallest.length)];
    const pickSet = new Set(pick);
    state.disabled = new Set();
    for (const s of state.current.sensors) {
      if (!pickSet.has(s.name)) state.disabled.add(s.name);
    }
    refreshSensorStyles();
    renderFailureList();
    scheduleEvaluate();
  });

  document.getElementById("replan-button").addEventListener("click", () => {
    if (state.replanInFlight) return;
    if (state.editSensors.length === 0) {
      const el = document.getElementById("replan-status");
      el.classList.add("dirty");
      el.textContent = "センサが 0 個では計算できません。1 つ以上配置してください。";
      return;
    }
    doReplan(state.editSensors.map(s => ({ ...s })));
  });
  document.getElementById("reset-button").addEventListener("click", () => {
    if (state.replanInFlight) return;
    doReplan(state.originalSensors.map(s => ({ ...s })));
  });
  document.getElementById("clear-all-button").addEventListener("click", clearAllSensors);

  const defaultRadiusSlider = document.getElementById("new-radius-slider");
  const defaultRadiusMeters = document.getElementById("new-radius-meters");
  defaultRadiusSlider.addEventListener("input", () => {
    state.newRadius = parseFloat(defaultRadiusSlider.value);
    if (state.current) {
      defaultRadiusMeters.textContent = Math.round(state.newRadius * state.current.side_meters);
    } else {
      defaultRadiusMeters.textContent = (state.newRadius * 2252).toFixed(0);
    }
  });

  document.getElementById("tiles-toggle").addEventListener("change", (e) => {
    document.body.classList.toggle("no-tiles", !e.target.checked);
  });

  if (regions.length) {
    sel.value = regions[0].id;
    await selectRegion(regions[0].id);
  } else {
    document.getElementById("R-hint").textContent = "プリセットが見つかりません。";
  }
}

window.addEventListener("DOMContentLoaded", () => {
  boot().catch(err => {
    console.error(err);
    alert("初期化に失敗しました: " + err.message);
  });
});
