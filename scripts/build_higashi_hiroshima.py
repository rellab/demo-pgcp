#!/usr/bin/env python3
"""Build the higashi_hiroshima preset from existing experiment5 artifacts.

Reads:
  - pgcp-experiments/experiment5/processed/campus_normalized.json
  - pgcp-experiments/experiment5/export.geojson
  - pgcp-experiments/experiment5/data/campus_sensors.json
  - pgcp-experiments/experiment5/results/campus_initial_triangulation.json
  - pgcp-experiments/experiment5/results/campus_importance.csv
  - pgcp-experiments/experiment5/results/campus_minsets.json
  - pgcp-experiments/experiment5/results/campus_summary.json
  - pgcp-experiments/experiment5/results/campus_sweep_pk.csv

Writes:
  - webapp/presets/higashi_hiroshima/meta.json
  - webapp/presets/higashi_hiroshima/polygon.json
  - webapp/presets/higashi_hiroshima/sensors.json
  - webapp/presets/higashi_hiroshima/triangulation.json
  - webapp/presets/higashi_hiroshima/precomputed.json
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent.parent
SRC = REPO / "pgcp-experiments" / "experiment5"
DST = REPO / "webapp" / "presets" / "higashi_hiroshima"


def load_json(path: Path):
    with path.open("r", encoding="utf-8") as fh:
        return json.load(fh)


def main() -> None:
    DST.mkdir(parents=True, exist_ok=True)

    norm = load_json(SRC / "processed" / "campus_normalized.json")
    geo = load_json(SRC / "export.geojson")
    sensors_in = load_json(SRC / "data" / "campus_sensors.json")
    mesh = load_json(SRC / "results" / "campus_initial_triangulation.json")
    minsets = load_json(SRC / "results" / "campus_minsets.json")
    summary = load_json(SRC / "results" / "campus_summary.json")

    # Extract lat/lon polygon (first ring of first Polygon feature).
    feat = geo["features"][0]
    geom = feat["geometry"]
    if geom["type"] == "Polygon":
        latlon_ring = geom["coordinates"][0]
    elif geom["type"] == "MultiPolygon":
        # Pick the largest ring by vertex count as a heuristic.
        latlon_ring = max(
            (ring for poly in geom["coordinates"] for ring in poly[:1]),
            key=len,
        )
    else:
        raise RuntimeError(f"unexpected geometry type: {geom['type']}")

    # GeoJSON coordinates are [lon, lat]. Convert to [lat, lon] tuples for clarity.
    polygon_latlon = [[lat, lon] for lon, lat in latlon_ring]

    # Compute lat/lon bbox.
    lats = [p[0] for p in polygon_latlon]
    lons = [p[1] for p in polygon_latlon]
    min_lat, max_lat = min(lats), max(lats)
    min_lon, max_lon = min(lons), max(lons)

    # Normalized bbox (always (0,0)-(1, h_ratio) per campus_normalized.json).
    bbox_n = norm["bbox_normalized"]  # [x0, y0, x1, y1]
    assert (bbox_n[0], bbox_n[1]) == (0, 0), "expected normalized bbox to start at (0,0)"

    meta = {
        "id": "higashi_hiroshima",
        "display_name": "広島大学 東広島キャンパス",
        "subtitle": "Higashi-Hiroshima Campus, Hiroshima University",
        "side_meters": norm["side_meters"],
        "bbox_normalized": bbox_n,
        "bbox_latlon": {
            "min_lat": min_lat,
            "max_lat": max_lat,
            "min_lon": min_lon,
            "max_lon": max_lon,
        },
        "n_sensors": len(sensors_in["circles"]),
        "default_pk": sensors_in.get("reliability", 0.9),
        # Prefer the converged maxLevel recorded in campus_summary.json over
        # the configuration knob in campus_sensors.json (which is the search
        # cap, not the converged depth). For this preset they happen to differ
        # (8 vs 11), and using the converged value makes Φ₁ ≡ Φ₂ hold at load
        # time, so the frontend's "convergence" indicator is meaningful.
        "maxlevel": summary.get("maxlevel", sensors_in.get("maxlevel", 11)),
    }

    # polygon.json carries both representations.
    polygon = {
        "polygon_normalized": norm["polygon"],
        "polygon_latlon": polygon_latlon,
    }

    # sensors.json (normalized only; lat/lon resolved in the frontend via meta bbox).
    sensors = {
        "circles": sensors_in["circles"],
        "reliability": sensors_in.get("reliability", 0.9),
        "maxlevel": sensors_in.get("maxlevel", 11),
    }

    # triangulation.json — strip the _meta block but keep vertex / triangle arrays
    # and the sensor->vertex index map (handy for the frontend).
    tri = {
        "vertices": mesh["vertices"],
        "triangles": mesh["triangles"],
    }
    if "sensor_vertex_index" in mesh:
        tri["sensor_vertex_index"] = mesh["sensor_vertex_index"]

    # Read importance CSV and reshape into a per-sensor map keyed by name.
    importance = []
    with (SRC / "results" / "campus_importance.csv").open("r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            importance.append(
                {
                    "name": row["name"],
                    "x": float(row["x"]),
                    "y": float(row["y"]),
                    "radius": float(row["radius"]),
                    "birnbaum": float(row["birnbaum"]),
                    "criticality1": float(row["criticality1"]),
                    "criticality0": float(row["criticality0"]),
                    "structure": float(row["structure"]),
                }
            )

    # Read pk_sweep CSV.
    pk_sweep = []
    with (SRC / "results" / "campus_sweep_pk.csv").open("r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            pk_sweep.append(
                {
                    "pk": float(row["p_k"]),
                    "R_lower": float(row["R_lower"]),
                    "R_upper": float(row["R_upper"]),
                }
            )

    # essential = sensors that appear as singleton min-cuts (R drops to 0 if disabled).
    essential = [
        sset[0] for sset in minsets["min_cut_sets"]["smallest"] if len(sset) == 1
    ]

    precomputed = {
        "R": summary["R_lower"],  # converged: lower == upper
        "default_pk": summary["reliability_p"],
        "peak_nodes": summary.get("peak_nodes"),
        "n_sensors": summary["n_sensors"],
        "convergence": summary.get("convergence", True),
        "phi_nodes": minsets.get("phi_nodes"),
        "importance": importance,
        "essential_sensors": essential,
        "min_cut_sets": minsets["min_cut_sets"],
        "min_path_sets": minsets["min_path_sets"],
        "pk_sweep": pk_sweep,
    }

    (DST / "meta.json").write_text(
        json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    (DST / "polygon.json").write_text(
        json.dumps(polygon, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    (DST / "sensors.json").write_text(
        json.dumps(sensors, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    (DST / "triangulation.json").write_text(
        json.dumps(tri, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    (DST / "precomputed.json").write_text(
        json.dumps(precomputed, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    print(f"wrote 5 files to {DST}")
    print(f"  n_sensors        = {meta['n_sensors']}")
    print(f"  bbox_latlon      = {meta['bbox_latlon']}")
    print(f"  R                = {precomputed['R']:.10f}")
    print(f"  essential        = {essential}")
    print(f"  n_min_cut_sets   = {minsets['min_cut_sets'].get('count_enumerated')}")
    print(f"  n_min_path_sets  = {minsets['min_path_sets'].get('count_enumerated')}")


if __name__ == "__main__":
    main()
