#!/bin/zsh

for filename in *.fasta; do
    /Users/priscillalau/Downloads/muscle3.8 -in $filename -out "$filename.realign.fasta"
done

for filename in *.realign.fasta; do
	trimal -in $filename -out "strictplus_$filename" -fasta -strictplus
    #alternatively
    #trimal -in $filename -out "custom_$filename" -fasta -gt 0.9 -cons 60 -w 3
    #trimal -in $filename -out "custom_$filename" -fasta -gappy out
done

#rename -e 's/.realign.fasta//'  strictplus_*.realign.fasta
