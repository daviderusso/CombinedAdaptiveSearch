#!/usr/bin/env bash
set -euo pipefail

# ===================== CONFIG =====================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INPUT_DIR="${SCRIPT_DIR}/data"
BASE_OUTDIR="${SCRIPT_DIR}/runs"

# Numero di ripetizioni per ogni coppia (input, runall)
REPS=1

# Elenco coppie (solo nomi file, senza path)
inputs=(
   "W.xy"
#   "USA.txt"
#   "CRT.txt"
#   "E.txt"
)

runalls=(
  "runall_W.txt"
#  "runall_USA.txt"
#  "runall_CRT.txt"
#  "runall_E.txt"
)

# Parametri comuni
TL="60.0"
REDH="1"
RUNALL_FLAG="1"

NIT="10"
#NIT="20"

#PERC_RED="0.0"
#PERC_RED="0.01"
#PERC_RED="0.05"
PERC_RED="0.1"

# ===================== CHECKS =====================

if (( ${#inputs[@]} == 0 )); then
  echo "No inputs configured." >&2
  exit 1
fi

if (( ${#inputs[@]} != ${#runalls[@]} )); then
  echo "List size mismatch: inputs=${#inputs[@]}, runalls=${#runalls[@]}" >&2
  exit 1
fi

mkdir -p "$BASE_OUTDIR"

# ===================== RUNS =====================

for i in "${!inputs[@]}"; do
  in_file="${inputs[$i]}"
  runall_file="${runalls[$i]}"

  # Costruzione path completi
  in_path="${INPUT_DIR}/${in_file}"
  runall_path="${INPUT_DIR}/${runall_file}"

  if [[ ! -f "$in_path" ]]; then
    echo "Missing input file: $in_path" >&2
    exit 1
  fi

  if [[ ! -f "$runall_path" ]]; then
    echo "Missing runall file: $runall_path" >&2
    exit 1
  fi

  # Nome base senza estensione
  in_stem="${in_file%.*}"

  for ((rep=1; rep<=REPS; rep++)); do
    outdir="${BASE_OUTDIR}/${in_stem}_rep$(printf "%02d" "$rep")_nit${NIT}_perc${PERC_RED}"
    mkdir -p "$outdir"

    echo "==> Running: input=$in_path | rep=$rep/$REPS"

    ./compute_paths \
      --input "$in_path" \
      --tl "$TL" \
      --redh "$REDH" \
      --nit "$NIT" \
      --runall "$RUNALL_FLAG" \
      --inputrunall "$runall_path" \
      --outdir "$outdir" \
      --perc_red "$PERC_RED"

  done
done
