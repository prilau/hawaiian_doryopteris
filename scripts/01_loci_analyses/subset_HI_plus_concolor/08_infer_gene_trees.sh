#!/bin/bash
# Infer an unrooted gene tree for every alignment in data_reduced/ with IQ-TREE:
#   -m MFP     ModelFinder Plus picks the best model per gene
#   -bb 10000  10,000 ultrafast bootstraps
#   -nt AUTO   autodetect threads
# IQ-TREE files land in iqtree_runs/<locus>/; the ML tree is copied to
# gene_trees_unrooted/<locus>.treefile
#
# IQ-TREE's own output is sent to iqtree_runs/<locus>/<locus>.log so the screen
# stays clean for the progress bar.
#
# Dirs can be overridden as positional args:
#   ./infer_gene_trees.sh data_reduced_vcf gene_trees_unrooted_vcf iqtree_runs_vcf
set -uo pipefail

IN=${1:-data_reduced}
OUT=${2:-gene_trees_unrooted}
WORK=${3:-iqtree_runs}

mkdir -p "$OUT" "$WORK"

TOTAL=$(ls "$IN"/*.fasta 2>/dev/null | wc -l | tr -d ' ')
[ "$TOTAL" -eq 0 ] && { echo "No .fasta files in $IN"; exit 1; }

fmt () { printf '%d:%02d:%02d' $(( $1 / 3600 )) $(( ($1 % 3600) / 60 )) $(( $1 % 60 )); }

draw () {   # draw <done> <label>
    local done=$1 label=$2 width=40
    local filled=$(( done * width / TOTAL ))
    local pct=$(( done * 100 / TOTAL ))
    local eta=0
    [ "$done" -gt 0 ] && eta=$(( SECONDS * (TOTAL - done) / done ))
    local bar
    bar=$(printf '%*s' "$filled" '' | tr ' ' '#')
    bar+=$(printf '%*s' "$(( width - filled ))" '' | tr ' ' '.')
    printf '\r[%s] %d/%d (%d%%)  elapsed %s  eta %s  %-20s' \
        "$bar" "$done" "$TOTAL" "$pct" "$(fmt $SECONDS)" "$(fmt $eta)" "$label"
}

DONE=0 FAIL=0
for f in "$IN"/*.fasta; do
    locus=$(basename "$f" .fasta)

    if [ -s "$OUT/$locus.treefile" ]; then
        DONE=$(( DONE + 1 )); draw "$DONE" "$locus (cached)"; continue
    fi
    if [ "$(grep -c '^>' "$f")" -lt 4 ]; then
        DONE=$(( DONE + 1 )); draw "$DONE" "$locus (skip <4 taxa)"; continue
    fi

    draw "$DONE" "$locus ..."
    mkdir -p "$WORK/$locus"
    if iqtree -s "$f" -m MFP -bb 10000 -nt AUTO -pre "$WORK/$locus/$locus" -redo \
            > "$WORK/$locus/$locus.screen" 2>&1; then
        cp "$WORK/$locus/$locus.treefile" "$OUT/$locus.treefile"
    else
        FAIL=$(( FAIL + 1 ))
        printf '\r%*s\r' 100 ''            # clear the bar line
        echo "!! $locus failed -- see $WORK/$locus/$locus.log"
    fi
    DONE=$(( DONE + 1 ))
    draw "$DONE" "$locus done"
done

draw "$DONE" "complete"
printf '\n\nFinished: %d trees in %s/  (%d failed)\n' "$(( DONE - FAIL ))" "$OUT" "$FAIL"
