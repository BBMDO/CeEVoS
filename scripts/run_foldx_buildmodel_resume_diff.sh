#!/usr/bin/env bash
set -euo pipefail

PROJECT="$(pwd)"
REPAIRED="$PROJECT/07_functional_inference/foldx/repaired_renum"
MUTLIST="$PROJECT/07_functional_inference/foldx/nonempty_mutfiles_diff_chainA.list"
OUTROOT="$PROJECT/07_functional_inference/foldx/buildmodel_diff_chainA"
ROTABASE_SRC="/home/danilo/Danilo/rotabase.txt"

mkdir -p "$OUTROOT"

[[ -s "$MUTLIST" ]] || { echo "[ERROR] MUTLIST vazio: $MUTLIST" >&2; exit 1; }
[[ -f "$ROTABASE_SRC" ]] || { echo "[ERROR] rotabase.txt não encontrado: $ROTABASE_SRC" >&2; exit 1; }

while IFS= read -r mutfile; do
  [[ -z "$mutfile" ]] && continue
  [[ -f "$mutfile" ]] || { echo "[WARN] mutfile não existe: $mutfile"; continue; }

  base="$(basename "$mutfile" .txt)"  # OGxxxx__cereus_..._to_...__diff
  og="${base%%__*}"

  tmp="${base#*__}"                # cereus_..._to_...__diff
  fromtok="${tmp%%_to_*}"          # cereus_...

  pdb="$REPAIRED/${og}__${fromtok}_Repair.pdb"
  [[ -f "$pdb" ]] || { echo "[WARN] PDB não encontrado: $pdb"; continue; }

  rundir="$OUTROOT/$base/${og}__${fromtok}_Repair"
  mkdir -p "$rundir"

  if ls "$rundir"/Raw_*.fxout >/dev/null 2>&1; then
    echo "SKIP (já tem Raw): $base"
    continue
  fi

  echo "RUN BuildModel: $base | $(basename "$pdb")"

  cp "$pdb" "$rundir/$(basename "$pdb")"
  cp "$mutfile" "$rundir/individual_list.txt"
  ln -sf "$ROTABASE_SRC" "$rundir/rotabase.txt"

  (
    cd "$rundir"
    FoldX --command=BuildModel --pdb="$(basename "$pdb")" --mutant-file="individual_list.txt" > buildmodel.log 2>&1
  ) || {
    echo "[WARN] BuildModel falhou: $rundir"
    continue
  }

done < "$MUTLIST"

echo "DONE"
