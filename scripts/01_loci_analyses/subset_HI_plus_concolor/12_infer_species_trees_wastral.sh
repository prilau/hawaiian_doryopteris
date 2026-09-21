# wastral; subset of samples; nodes are not collapsed regardless of support in the gene trees
## Hawaiian specimens excluding D. subdecipiens; D. concolor as outgroup
wastral -i astral_input/gene_trees_50.tre -o astral_output_trees/wastral_50.tre
## Hawaiian specimens including D. subdecipiens; D. concolor as outgroup
wastral -i astral_input/gene_trees_50_with_subdecip.tre -o astral_output_trees/wastral_50_with_subdecip.tre
