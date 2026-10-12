#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ENV_PREFIX="$ROOT/.conda-env"
QUPKAKE_DIR="$ROOT/dependencies/QupKake"
QUPKAKE_REPOSITORY="https://github.com/Shualdon/QupKake.git"
QUPKAKE_COMMIT="f7a294e91929e4e4f0a7100ad5cf5a64a34a7d4f"
PATCH_FILE="$ROOT/patches/qupkake-traybake-compat.patch"

if ! command -v conda >/dev/null 2>&1; then
    echo "conda is required to create the TRAYBAKE environment." >&2
    exit 1
fi

mkdir -p "$ROOT/dependencies"
if [[ ! -f "$QUPKAKE_DIR/pyproject.toml" ]]; then
    if [[ -e "$QUPKAKE_DIR" ]]; then
        echo "dependencies/QupKake exists but is not a valid QupKake checkout." >&2
        exit 1
    fi
    git clone --no-checkout "$QUPKAKE_REPOSITORY" "$QUPKAKE_DIR"
    git -C "$QUPKAKE_DIR" checkout --detach "$QUPKAKE_COMMIT"
fi

if ! git -C "$QUPKAKE_DIR" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    echo "dependencies/QupKake is not a Git checkout." >&2
    exit 1
fi
QUPKAKE_HEAD="$(git -C "$QUPKAKE_DIR" rev-parse HEAD)"
if [[ "$QUPKAKE_HEAD" != "$QUPKAKE_COMMIT" ]]; then
    echo "QupKake must be checked out at $QUPKAKE_COMMIT (found $QUPKAKE_HEAD)." >&2
    echo "Preserve any local work, then check out the pinned commit and rerun setup." >&2
    exit 1
fi

if [[ -d "$ENV_PREFIX" ]]; then
    echo "Environment already exists at $ENV_PREFIX; remove it explicitly to recreate it." >&2
    exit 1
fi

if git -C "$ROOT/dependencies/QupKake" apply --reverse --check "$PATCH_FILE" >/dev/null 2>&1; then
    echo "TRAYBAKE's QupKake compatibility patch is already applied."
elif git -C "$ROOT/dependencies/QupKake" apply --check "$PATCH_FILE" >/dev/null 2>&1; then
    git -C "$ROOT/dependencies/QupKake" apply "$PATCH_FILE"
else
    echo "QupKake checkout does not match the expected compatibility patch." >&2
    exit 1
fi
conda env create --prefix "$ENV_PREFIX" --file "$ROOT/environment.yml"
conda run --prefix "$ENV_PREFIX" python -m pip install --no-deps "$QUPKAKE_DIR"

echo "TRAYBAKE environment is ready. Run ./traybake.sh or conda activate $ENV_PREFIX"
