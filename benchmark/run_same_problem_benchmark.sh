#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${ROOT_DIR}/build"

module purge
module load gcc cmake lib/openblas

cmake -S "${ROOT_DIR}" -B "${BUILD_DIR}" -DCMAKE_BUILD_TYPE=Release >/dev/null
cmake --build "${BUILD_DIR}" -j >/dev/null

echo "== Fortran benchmark =="
timeout 120s "${BUILD_DIR}/bench-fortran-same-problem"

echo "== Julia benchmark =="
timeout 120s julia --startup-file=no --project="${ROOT_DIR}/julia" "${ROOT_DIR}/julia/benchmark/bench_julia_same_problem.jl"
