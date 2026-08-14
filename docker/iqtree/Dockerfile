# Multi-arch IQ-TREE 2 + TreeTime image for the phylo/ module.
#
# Bundles both tools in one image via bioconda — reviewers cloning the
# repo don't want to juggle two containers when the tools always run
# together (ML tree → time-scaled tree).
#
# Build (any arch):
#   docker build -t hivdrwatch/phylo:latest docker/iqtree
#
# Usage from phylo/iqtree.py:
#   docker run --rm \
#     -v $PWD/data/interim/curated:/in:ro \
#     -v $PWD/data/interim/trees/<subtype>:/out \
#     hivdrwatch/phylo:latest \
#     iqtree2 -s /in/<subtype>.fasta -m GTR+F+I+G4 --seed 42 --prefix /out/<subtype>
#
# Usage for TreeTime:
#   docker run --rm \
#     -v $PWD/data/interim/trees/<subtype>:/work \
#     hivdrwatch/phylo:latest \
#     treetime --tree /work/<subtype>.treefile --aln /work/aln.fasta \
#              --dates /work/dates.csv --outdir /work/timed
#
# Pinned versions land in /versions.txt so `docker inspect` and
# in-container calls can verify what's baked in — matches MODULES.md's
# "every claim has a benchmark" rule.

FROM mambaorg/micromamba:1.5.10-bookworm-slim

# Root for install; micromamba's default user is `mambauser` (uid 1000)
USER root

ARG IQTREE_VERSION=2.3.6
ARG TREETIME_VERSION=0.11.4

# bioconda channel provides both tools. Pin exact versions so reruns are
# bit-identical for the "manuscript numbers must be reproducible" rule.
RUN micromamba install -n base -y \
        -c conda-forge -c bioconda \
        "iqtree=${IQTREE_VERSION}" \
        "treetime=${TREETIME_VERSION}" \
        "python=3.11" \
    && micromamba clean --all --yes

# Make base env's binaries + libs visible to non-login shells.
ENV PATH=/opt/conda/bin:$PATH \
    LD_LIBRARY_PATH=/opt/conda/lib:${LD_LIBRARY_PATH}

# Sanity + version pin file. Read from Python via /versions.txt.
RUN iqtree2 --version 2>&1 | head -1 > /versions.txt \
    && treetime --version 2>&1 | head -1 >> /versions.txt \
    && cat /versions.txt

WORKDIR /work
# Default entrypoint is iqtree2; the caller overrides with treetime as needed.
ENTRYPOINT []
CMD ["iqtree2", "--help"]
