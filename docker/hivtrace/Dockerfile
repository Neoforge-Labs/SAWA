# Pinned tn93 image for HIV-DR-Watch clustering stage.
#
# `tn93` (github.com/veg/tn93) is the CDC HIV-TRACE distance tool.
# It computes Tamura-Nei 1993 pairwise distances on aligned sequences —
# the standard for HIV molecular epidemiology since ~2010.
#
# Build:
#   docker build -t hivdrwatch/tn93:1.0.13 docker/hivtrace
#
# Usage (from cluster/distance.py):
#   docker run --rm \
#     -v $PWD/data/interim/curate:/work \
#     hivdrwatch/tn93:1.0.13 \
#     tn93 -t 0.015 -o /work/distances.csv /work/aligned.fasta
#
# Output CSV: ID1,ID2,Distance — one row per pair BELOW the threshold.
# We don't emit above-threshold pairs (there are billions of them and
# they're by definition not clusters).

FROM debian:12-slim

ARG TN93_VERSION=1.0.13

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential cmake git ca-certificates \
    && git clone --depth 1 --branch v${TN93_VERSION} https://github.com/veg/tn93.git /tmp/tn93 \
    && cmake -B /tmp/tn93/build -S /tmp/tn93 \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
    && cmake --build /tmp/tn93/build --parallel \
    && cmake --install /tmp/tn93/build \
    && rm -rf /tmp/tn93 \
    && apt-get purge -y build-essential cmake git \
    && apt-get autoremove -y \
    # Keep libgomp1 (OpenMP runtime, needed by tn93 at runtime) and
    # ca-certificates (in case downstream tools need TLS). Purging
    # build-essential also purges these transitively, so re-install
    # them explicitly here after the purge.
    && apt-get install -y --no-install-recommends libgomp1 ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Sanity check at build time — this used to be silent success but the
# missing libgomp broke it at runtime. Now if the runtime lib is
# missing this build layer fails loudly.
RUN tn93 --help 2>&1 | head -3 \
    && ldd $(which tn93) | grep -q "not found" && exit 1 || true

WORKDIR /work
ENTRYPOINT ["tn93"]
CMD ["--help"]
