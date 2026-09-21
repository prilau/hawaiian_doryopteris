#!/bin/bash
# Bundle the per-locus gene trees into ASTRAL input files:
#   astral_input/gene_trees_all.tre    all trees, unchanged (one per line)
#   astral_input/gene_trees_bs10.tre   same trees, internal branches with
#                                      UFBoot support < 10 collapsed to polytomies
set -uo pipefail

IN=gene_trees_unrooted
OUT=astral_input
THRESH=10

mkdir -p "$OUT"

# --- version 1: as-is, one newick per line --------------------------------
: > "$OUT/gene_trees_all.tre"
for f in "$IN"/*.treefile; do
    tr -d '\n' < "$f"; echo
done >> "$OUT/gene_trees_all.tre"
