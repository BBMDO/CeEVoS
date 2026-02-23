#!/bin/bash
set -e

FOLDX_BIN=$(command -v FoldX || true)
if [ -z "$FOLDX_BIN" ]; then
  echo "ERRO: FoldX não está no PATH. Rode: which FoldX"
  exit 1
fi

PDB_DIR="$(pwd)/07_functional_inference/foldx/repaired"
OUT_DIR="$(pwd)/07_functional_inference/foldx/results"
mkdir -p "$OUT_DIR"

for pdb in "$PDB_DIR"/*.pdb; do
  base=$(basename "$pdb" .pdb)
  echo "Repair: $base"
  "$FOLDX_BIN" --command=RepairPDB --pdb="$pdb" > "$OUT_DIR/$base.repair.log" 2>&1 || echo "Repair falhou: $pdb"
done
