#!/usr/bin/env bash
set -euo pipefail

ALIGNMENT="${1:-MATK_dory.aligned.fasta}"
THREADS="${THREADS:-AUTO}"
SEED="${SEED:-12345}"
PREFIX="${PREFIX:-MATK_dory_ML}"

if [[ ! -f "$ALIGNMENT" ]]; then
  echo "Error: alignment file not found: $ALIGNMENT" >&2
  echo "Usage: $0 [alignment.fasta]" >&2
  exit 1
fi

if command -v iqtree2 >/dev/null 2>&1; then
  IQTREE_BIN="iqtree2"
elif command -v iqtree >/dev/null 2>&1; then
  IQTREE_BIN="iqtree"
else
  echo "Error: iqtree2/iqtree not found in PATH." >&2
  echo "Install IQ-TREE 2 and rerun." >&2
  exit 1
fi

echo "Running IQ-TREE with ModelFinder + 10,000 rapid bootstrap replicates..."
echo "Alignment: $ALIGNMENT"
echo "Prefix:    $PREFIX"
echo "Seed:      $SEED"
echo "Threads:   $THREADS"

"$IQTREE_BIN" \
  -s "$ALIGNMENT" \
  -m MFP \
  -bb 10000 \
  -bnni \
  -nt "$THREADS" \
  -seed "$SEED" \
  -pre "$PREFIX"

echo
echo "Done. Key outputs:"
echo "  ${PREFIX}.treefile   (best ML tree)"
echo "  ${PREFIX}.iqtree     (run summary, selected model)"
echo "  ${PREFIX}.log        (run log)"
echo "  ${PREFIX}.contree    (consensus bootstrap tree)"
