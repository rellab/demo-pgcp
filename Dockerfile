# PGCP webapp container — self-contained build.
#
# Builds on stock julia:1.10 (no dependency on the cudd-julia base image of
# pgcp-experiments). Installs MiniCUDD (built from source from upstream git),
# DelaunayTriangulation, Oxygen, and the few JSON / HTTP / GeometryBasics /
# ProgressMeter packages the server and the solver need. Skip Pkg.test for
# MiniCUDD here — the upstream test suite is unrelated to our use.
#
# Build (from the webapp/ directory):
#   docker compose up -d --build
#
# Or by hand:
#   docker build -t pgcp-webapp -f Dockerfile .
#   docker run --rm -p 8080:8080 pgcp-webapp

FROM julia:1.10

ENV DEBIAN_FRONTEND=noninteractive \
    LD_LIBRARY_PATH="/usr/local/lib:/usr/local/lib64:/usr/lib:/usr/lib/x86_64-linux-gnu"

# Build toolchain for MiniCUDD (Pkg.build invokes ./configure && make).
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    curl \
    ca-certificates \
    pkg-config \
    file \
    git \
    automake \
    autoconf \
    libtool \
    m4 \
    perl \
    && rm -rf /var/lib/apt/lists/*

# MiniCUDD: not in the General registry, built from source.
RUN julia --color=yes -e 'using Pkg; \
    Pkg.add(url="https://github.com/JuliaReliab/MiniCUDD.jl.git"); \
    Pkg.build("MiniCUDD"); \
    Pkg.precompile()'

# Everything else from the General registry. Kept minimal — we deliberately
# do not pull in Plots / Makie / etc to keep the image small.
RUN julia --color=yes -e 'using Pkg; \
    Pkg.add(["JSON", "JSON3", "HTTP", "GeometryBasics", "ProgressMeter", \
             "DelaunayTriangulation", "Oxygen"]); \
    Pkg.precompile()'

WORKDIR /app
COPY . /app

# Pre-instantiate the active project so the first request does not pay the cost.
RUN julia --project=/app/server -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

EXPOSE 8080
CMD ["julia", "--project=/app/server", "/app/server/server.jl"]
