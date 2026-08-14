# HIV-TRACE clustering image.
#
# Wraps `hivnetworkcsv` from the `hivclustering` package — the CDC
# HIV-TRACE reference clusterer. Given a TN93 distances CSV
# (id1,id2,distance), it emits a JSON with nodes + edges from which we
# extract connected-component cluster membership.
#
# We deliberately DON'T bundle tn93 here — the tn93 image
# (hivdrwatch/tn93) already produces the CSV in the shape hivnetworkcsv
# expects, and having a single tn93 producer keeps our clusterer and
# HIV-TRACE working from bit-identical distances.
#
# Build:
#   docker build -t hivdrwatch/hivtrace:0.1 docker/hivtrace_full
#
# Usage (from cluster/hivtrace.py):
#   docker run --rm \
#     -v $PWD/work:/work \
#     hivdrwatch/hivtrace:0.1 \
#     hivnetworkcsv -i /work/distances.csv -t 0.015 -f plain -J /work/network.json

FROM python:3.12-slim

# hivclustering pulls scipy/numpy/networkx; grab a compiler layer just
# for the wheel build, then drop it. Keep libgomp1 for scipy runtime.
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential ca-certificates \
    # hivclustering is versioned 1.9.x on PyPI; pin to a stable minor.
    # (1.7.1 exists in some historical release notes but never landed on PyPI.)
    && pip install --no-cache-dir 'hivclustering>=1.9.5,<2.0' \
    && apt-get purge -y build-essential \
    && apt-get autoremove -y \
    && apt-get install -y --no-install-recommends libgomp1 ca-certificates \
    && rm -rf /var/lib/apt/lists/* /root/.cache

# Sanity check at build time.
RUN hivnetworkcsv --help 2>&1 | head -3

WORKDIR /work
ENTRYPOINT ["hivnetworkcsv"]
CMD ["--help"]
