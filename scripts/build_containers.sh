#!/usr/bin/env bash
# Build every container image on the EC2 host (x86_64). Tags are recorded in versions.yml.
set -euo pipefail
cd "$(dirname "$0")/../containers"
P=longread-platform-benchmark
docker build -q -t $P/minimap2:2.30     -f Dockerfile.minimap2 .
docker build -q -t $P/sniffles:2.8.0    -f Dockerfile.sniffles .
docker build -q -t $P/cutesv:2.1.4      -f Dockerfile.cutesv .
docker build -q -t $P/truvari:5.4.0     -f Dockerfile.truvari .
docker build -q -t $P/clair3:2.0.2      -f Dockerfile.variantcall .
docker build -q -t $P/happy:0.3.15      -f Dockerfile.happy .
docker build -q -t $P/whatshap:2.8      -f Dockerfile.phasing .
docker build -q -t $P/modkit:0.6.4      -f Dockerfile.modkit .
docker build -q -t $P/analysis:1.0      -f Dockerfile.analysis .
docker images | grep "^$P"
