# Pinned tn93 image for HIV-DR-Watch clustering stage.
#
# `tn93` (github.com/veg/tn93) is the CDC HIV-TRACE distance tool.
# It computes Tamura-Nei 1993 pairwise distances on aligned sequences —
# the standard for HIV molecular epidemiology since ~2010. We keep it as
# the reference implementation to A/B against our own distance engine.
#
# Build:
#   docker build -t hivdrwatch/tn93:1.0.17 docker/hivtrace
#
# Usage (from cluster/distance.py):
#   docker run --rm \
#     -v $PWD/data/interim/curate:/work \
#     hivdrwatch/tn93:1.0.17 \
#     tn93 -t 0.015 -o /work/distances.csv /work/aligned.fasta
#
# Output CSV: ID1,ID2,Distance — one row per pair BELOW the threshold.
# We don't emit above-threshold pairs (there are billions of them and
# they're by definition not clusters).
#
# Why 1.0.17 and not 1.0.13: on aarch64 (every Apple Silicon build host)
# plain `char` is unsigned, and tn93 up to and including v1.0.14 declared
# its character-class table as `char validFlags[256]` and used -1 as the
# "not a sequence character" sentinel. Unsigned, -1 reads back as 255, so
# `validFlags[c] >= 0` was true for every byte — including '>'. The FASTA
# reader therefore never saw a record boundary after the first header and
# concatenated the whole file into one sequence, which is why the image
# silently emitted zero pairs on every cohort we gave it. Upstream fixed
# it in v1.0.15 as `signed char` (veg/tn93 commit 5c14f35, issue #46).
# v1.0.17 (commit ef82c1b) is the last 1.0.x patch release; v1.1.0 adds a
# new subquadratic code path we don't want in a reference binary.
#
# v1.0.13 was also mis-tagged upstream: that tag never bumped VERSION.txt,
# so the binary it builds reports v1.0.12. The version assertion below
# pins tag and binary together so that cannot drift again unnoticed.

FROM debian:12-slim

ARG TN93_VERSION=1.0.17

# python3 is a build-only dependency: from v1.0.15 tn93's CMakeLists calls
# `find_package(PythonInterp 3 REQUIRED)` unconditionally to wire up its
# ctest suite, with no option to skip it. It is purged below with the rest
# of the toolchain, so it does not reach the final image.
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential cmake git ca-certificates python3 \
    && git clone --depth 1 --branch v${TN93_VERSION} https://github.com/veg/tn93.git /tmp/tn93 \
    && cmake -B /tmp/tn93/build -S /tmp/tn93 \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
    && cmake --build /tmp/tn93/build --parallel \
    && cmake --install /tmp/tn93/build \
    && rm -rf /tmp/tn93 \
    && apt-get purge -y build-essential cmake git python3 \
    && apt-get autoremove -y \
    # Keep libgomp1 (OpenMP runtime, needed by tn93 at runtime) and
    # ca-certificates (in case downstream tools need TLS). Purging
    # build-essential also purges these transitively, so re-install
    # them explicitly here after the purge.
    && apt-get install -y --no-install-recommends libgomp1 ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Build-time acceptance tests. `tn93 --help` on its own passed for the
# whole time the FASTA reader was broken, so it is not a sanity check:
# it never reads a sequence. Each check below fails the layer loudly.
#
# Note the absence of `|| true` — the previous revision ended its check
# with one, which swallowed any failure in the chain ahead of it.
RUN set -eux; \
    # 1. The runtime linker can resolve everything (libgomp regression).
    ! ldd "$(which tn93)" | grep -q "not found"; \
    # 2. The binary's self-reported version matches the tag we cloned,
    #    so a mis-tagged upstream release cannot ship as something else.
    test "$(tn93 -v 2>&1 | tr -d '[:space:]')" = "v${TN93_VERSION}"; \
    # 3. The FASTA reader finds BOTH records of a two-record file and
    #    emits the one pair between them. The two sequences overlap on
    #    300 unambiguous sites (tn93's `-l` default demands strictly more
    #    than 100) and differ at the first 4. Under the aarch64 unsigned-
    #    char bug this file parsed as one concatenated sequence -- headers
    #    and newlines included -- and produced an empty CSV, so this is
    #    the check that would have caught it.
    mkdir -p /tmp/selftest; \
    { printf '>selftest_a\n'; \
      for _ in $(seq 75); do printf 'ACGT'; done; printf '\n'; \
      printf '>selftest_b\n'; \
      printf 'TGCA'; for _ in $(seq 73); do printf 'ACGT'; done; printf 'ACGT\n'; \
    } > /tmp/selftest/pair.fasta; \
    test "$(grep -c '^>' /tmp/selftest/pair.fasta)" = 2; \
    tn93 -t 1.0 -o /tmp/selftest/out.csv /tmp/selftest/pair.fasta >/tmp/selftest/log 2>&1; \
    grep -Eq '"Sequences"[[:space:]]*:[[:space:]]*2' /tmp/selftest/log; \
    test "$(tail -n +2 /tmp/selftest/out.csv | grep -c .)" = 1; \
    rm -rf /tmp/selftest

WORKDIR /work
ENTRYPOINT ["tn93"]
CMD ["--help"]
