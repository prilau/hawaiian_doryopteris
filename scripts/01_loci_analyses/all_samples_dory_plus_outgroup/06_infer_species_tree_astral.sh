# for astral

# make input file
cat gt_output/*.tre > in.trees

# collapse low support branches in gene trees (need newick utils)
nw_ed  in.trees 'i & b<=10' o > in_BS10.tre


# run astral
java -jar /Applications/Phylogeny_Programs/Astral/astral.5.7.7.jar -i in_BS10.tre -o output/dory_HI_BS10.tre 2> output/dory_HI_BS10.log
