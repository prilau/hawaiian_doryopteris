#!/bin/zsh

for filename in *.realign.fasta; do
	trimal -in $filename -out "strictplus_$filename" -fasta -strictplus
done

rename -e 's/.realign.fasta//'  strictplus_*.realign.fasta
