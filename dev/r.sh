#!/usr/bin/env bash
# Run Rscript inside the FHI R container (gunnar-r) with this project mounted
# at /work and the project-local library .Rlib/ (odin2, dust2, monty, and this
# package once installed) first on the library path.
#
#   dev/r.sh -e 'devtools::test()'
#   dev/r.sh examples/01_single_population.R
#
# Override the image with R_IMAGE=... dev/r.sh ...
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
mkdir -p "$root/.Rlib"
exec docker run --rm \
  --user "$(id -u):$(id -g)" -e HOME=/tmp \
  -e R_LIBS=/work/.Rlib \
  -v "$root":/work -w /work \
  "${R_IMAGE:-gunnar-r:latest}" Rscript "$@"
