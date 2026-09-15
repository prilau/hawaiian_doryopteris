#!/bin/bash
for filename in *.fasta; do
	iqtree -bb 1000 -nt 4 -s $filename
done

mv loci/*.treefile gt_output

rename -e 's/.treefile/.tre/'  *.treefile
