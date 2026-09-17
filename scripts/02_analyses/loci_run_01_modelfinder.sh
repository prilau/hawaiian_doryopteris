#!/bin/bash

for filename in *.fasta; do
	./bin/iqtree2 -m TEST -s $filename
done

grep "chosen" loci/*.log > model_results.log
