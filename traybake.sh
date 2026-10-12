#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ENV_PREFIX="$ROOT/.conda-env"

if [[ ! -x "$ENV_PREFIX/bin/python" ]]; then
    echo "TRAYBAKE environment not found; run ./setup_environment.sh first." >&2
    exit 1
fi

export PYTHONPATH="$ROOT/src${PYTHONPATH:+:$PYTHONPATH}"
export PATH="$ENV_PREFIX/bin:$PATH"
export OMP_NUM_THREADS=1
export OMP_THREAD_LIMIT=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
exec "$ENV_PREFIX/bin/python" -m traybake "$@"
