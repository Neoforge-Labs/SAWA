# Multi-arch MAFFT image for HIV-DR-Watch alignment stage.
#
# Build (any arch):
#   docker build -t hivdrwatch/mafft:latest docker/mafft
#
# Two arch paths:
#   * amd64: pin MAFFT 7.526 from cbrc.jp .deb (matches everyone else's images)
#   * arm64: install from Debian 12's mafft package (7.505 as of bookworm).
#     Version parity across arch is imperfect but the differences don't
#     affect our --add --keeplength usage — alignment output is identical
#     to within noise. If you need bit-exact reproducibility, force
#     `--platform=linux/amd64` at build time (Rosetta on Apple Silicon).
#
# Usage (from curate/align.py):
#   docker run --rm \
#     -v $PWD/data/reference:/ref:ro \
#     -v $PWD/data/interim/curate:/work \
#     hivdrwatch/mafft:latest \
#     mafft --add /work/new.fasta --keeplength --reorder /ref/hxb2_pol.fasta \
#     > /work/aligned.fasta

FROM debian:12-slim

# TARGETARCH is populated by BuildKit ("amd64" / "arm64" / ...).
ARG TARGETARCH
ARG MAFFT_VERSION=7.526

RUN set -eux; \
    apt-get update; \
    if [ "$TARGETARCH" = "amd64" ]; then \
        apt-get install -y --no-install-recommends wget ca-certificates; \
        wget -q "https://mafft.cbrc.jp/alignment/software/mafft_${MAFFT_VERSION}-1_amd64.deb"; \
        dpkg -i "mafft_${MAFFT_VERSION}-1_amd64.deb"; \
        rm "mafft_${MAFFT_VERSION}-1_amd64.deb"; \
        apt-get purge -y wget ca-certificates; \
        apt-get autoremove -y; \
    else \
        # arm64 (Apple Silicon, ARM servers) — apt-provided mafft.
        # Version will differ from MAFFT_VERSION; recorded in /mafft.version below.
        apt-get install -y --no-install-recommends mafft; \
    fi; \
    rm -rf /var/lib/apt/lists/*

# Sanity check at build time — records the exact version baked in so
# `docker inspect` and mounted runs can see it.
RUN mafft --version 2>&1 | tee /mafft.version

WORKDIR /work
ENTRYPOINT ["mafft"]
CMD ["--help"]
