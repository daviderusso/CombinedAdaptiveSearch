#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$ROOT_DIR"

CC="gcc"
CFLAGS="-O2 -std=c11 -Wall -Wextra -Wno-unused-parameter"
LDFLAGS=""
OUT_DIR="build"

# Parametri run: modifica qui per scegliere dataset/risultati
GRAPH_FILE="data/W.xy"
#INSTANCES_FILE="data/160WInstancesHeur_AP.csv"
INSTANCES_FILE="data/1WInstancesHeur_AP.csv"
RES_DIR="res2"
TIME_LIMIT="60"
BIN="$OUT_DIR/bch"

mkdir -p "$OUT_DIR"

$CC $CFLAGS main.c -o "$BIN" -lm $LDFLAGS

"$BIN" --graph "$GRAPH_FILE" --instances "$INSTANCES_FILE" --timelimit "$TIME_LIMIT" --res-dir "$RES_DIR"
