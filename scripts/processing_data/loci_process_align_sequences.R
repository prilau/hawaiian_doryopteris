setwd("~/Box/project_Doryopteris/analyses/phyparts/")

t_iq <- ape::read.tree("~/Box/project_Doryopteris/analyses/preliminary_ml/output_normal_168hr.tre")
tips_to_keep <- t_iq$tip.label

files <- list.files("original_loci/")

#get_fastas <- function(f) {return(length(unlist(strsplit(f, split = "\\."))) == 2)}
#r <- unlist(lapply(files, get_fastas))
#files_fasta <- files[r]

for (i in 250:length(files)) {
  aln <- ape::read.dna(paste0("original_loci/",files[i]), format="fasta", as.matrix=TRUE)
  #rownames(aln) <- gsub("\t","", rownames(aln), fixed = T)
  aln_subbed <- aln[rownames(aln) %in% tips_to_keep,]
  aln_subbed <- ips::mafft(aln_subbed)
  ape::write.FASTA(aln_subbed, file = paste0("edited_loci/", files[i]))
}
