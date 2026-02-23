#!/usr/bin/env bash
set -euo pipefail

PROJECT="$(pwd)"
REPAIRED="$PROJECT/07_functional_inference/foldx/repaired_renum"
MUTLIST="$PROJECT/07_functional_inference/foldx/nonempty_mutfiles_chainA.list"
OUTROOT="$PROJECT/07_functional_inference/foldx/buildmodel_nonempty"
ROTABASE_SRC="/home/danilo/Danilo/rotabase.txt"

mkdir -p "$OUTROOT"

if [[ ! -s "$MUTLIST" ]]; then
  echo "[ERROR] Lista de mutfiles não existe ou está vazia: $MUTLIST" >&2
  exit 1
fi
if [[ ! -f "$ROTABASE_SRC" ]]; then
  echo "[ERROR] rotabase.txt não encontrado em: $ROTABASE_SRC" >&2
  exit 1
fi

while IFS= read -r mutfile; do
  [[ -z "$mutfile" ]] && continue
  [[ ! -f "$mutfile" ]] && { echo "[WARN] mutfile não existe: $mutfile"; continue; }

  base="$(basename "$mutfile" .txt)"     # ex: OG0028976__cereus_19790_to_fffe_19967__pp0.95
  og="${base%%__*}"                      # OG0028976
  tmp="${base#*__}"                      # cereus_19790_to_fffe_19967__pp0.95
  fromtok="${tmp%%_to_*}"                # cereus_19790

  pdb="$REPAIRED/${og}__${fromtok}_Repair.pdb"
  if [[ ! -f "$pdb" ]]; then
    echo "[WARN] PDB repaired_renum não encontrado para $base: $pdb"
    continue
  fi

  rundir="$OUTROOT/$base/${og}__${fromtok}_Repair"
  mkdir -p "$rundir"

  # condição de "já finalizado"
  if ls "$rundir"/Dif_*.fxout >/dev/null 2>&1; then
    echo "SKIP (já finalizado): $base"
    continue
  fi

  echo "RUN BuildModel: $base | $(basename "$pdb")"

  # prepara inputs no rundir
  cp "$pdb" "$rundir/$(basename "$pdb")"
  cp "$mutfile" "$rundir/individual_list.txt"
  ln -sf "$ROTABASE_SRC" "$rundir/rotabase.txt"

  # roda FoldX
  (
    cd "$rundir"
    FoldX --command=BuildModel --pdb="$(basename "$pdb")" --mutant-file="individual_list.txt" > buildmodel.log 2>&1
  ) || {
    echo "[WARN] BuildModel falhou: $rundir"
    # mantém log para inspeção
    continue
  }

done < "$MUTLIST"

echo "DONE"
