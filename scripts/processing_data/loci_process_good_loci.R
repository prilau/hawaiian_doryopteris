setwd("/Users/priscillalau/Box/project_Doryopteris/data/dory_combined/probe")

# compare number of duplicated taxa by loci and by taxon
library(ape)
library(tidyverse)

file_path <- "./"

# read in alignments 
alns_all <- list()
files <- list.files(path = file_path)
files_all <- files[grep("afterMerge.mafft.oneline.fasta", files)]

for (i in 1:length(files_all)) {
  a <- read.dna(file = paste0(file_path, files_all[i]), format = "fasta")
  if (is.null(a) == FALSE) {
    alns_all[[i]] <- as.alignment(a)
    names(alns_all)[i] <- strsplit(files_all[i],split = "\\.")[[1]][[1]]
  } else {
    alns_all[[i]] <- NA
    names(alns_all)[i] <- strsplit(files_all[i],split = "\\.")[[1]][[1]]
  }
}

# make data frame of taxon and number of seqs 

#get_extr_num <- function(names) {
#  n1 <- matrix(unlist(strsplit(names, split = "\\.")), ncol = 3, byrow = T)[,2]
#  gsub("_sp", "", x = n1)
#}
get_extr_num <- function(names) {
  matrix(unlist(strsplit(names, split = "_")), ncol = 2, byrow = T)[,1]
}


identify_dups <- function(tip, loci) {
  a <- alns_all[[loci]]
  if (class(a) == "alignment") {
    num_seqs <- sum(grepl(tip, a$nam))
    return((num_seqs))
  } else {return(NA)}
}

name_dict <- read.csv("../../../metadata/names_to_extr_num.csv")
all_names <- name_dict$ex

res <- as.data.frame(matrix(nrow = length(all_names), ncol = length(alns_all)))
colnames(res) <- names(alns_all)
rownames(res) <- all_names

for (i in 1:nrow(res)) {
  for (j in 1:ncol(res)) {
    tip_name <- all_names[i]
    locus_name <- names(alns_all)[j]
    res[tip_name,locus_name] <- identify_dups(tip = tip_name, loci = locus_name)
  }
}

res_subset <- res[,which(is.na(alns_all) == FALSE)]


perc_good_row <- function(row) {
  sum(row > 1)/sum(row == 1)
}
perc_good_col <- function(col) {
  sum(col > 1)/sum(col == 1)
}

v_row <- apply(res_subset, 1, FUN = perc_good_row)
res_subset_rowsort <- res_subset[order(v_row,decreasing=T),]

v_col <- apply(res_subset, 2, FUN = perc_good_col)
res_subset_colsort <- res_subset[,order(v_col,decreasing=T)]

good_loci <- v_col[which(v_col <= 0.1)]

### copy good loci into a new directory
for (i in 1:length(good_loci)) {
  from = paste0("nodups/",names(good_loci)[i],".nodups.fasta")
  to = paste0("good_loci/",names(good_loci)[i],".nodups.fasta")
  file.copy(from = from, to = to)
}


files <- list.files(path = "good_loci/")
files_all <- files[grep(".fasta", files)]

rename_alignment <- function(a, new_names) {
  for (i in 1:length(names(a))){
    old <- names(a)[i]
    if (grepl("_sp",old)){
      old <- unlist(strsplit(old, "_"))[1]
    }
    new <- new_names$name[which(new_names$ex == old)]
    paste0("Old: ", old, ", New: ", new)
    names(a)[i] <- new
  }
  return(a)
}

name_dict <- read.csv("~/Box/project_Doryopteris/metadata/names_to_extr_num.csv")

for (i in 1:length(files_all)) {
  a <- ape::read.FASTA(files_all[i])
  a_renamed <- rename_alignment(a,name_dict)
  new_file <- strsplit(files_all[i],split = "\\.")[[1]][[1]]
  ape::write.FASTA(a_renamed,file = paste0("~/Box/project_Doryopteris/data/dory_combined/probe/good_loci/good_loci_renamed/",new_file,".fasta"))
}





##### identify good loci and concat nodups of 
##### good loci for prelim analysis

good_loci <- v_col[which(v_col <= 0.1)]
#dir.create(paste0(file_path, "/good_loci"))
loci_to_concat <- list() 
for (i in 1:length(good_loci)) {
  from = paste0("data/dory_combined/probe/nodups/",names(good_loci)[i],".nodups.fasta")
  #to = paste0(file_path, "bad_loci/",names(good_loci)[i],".afterMerge.mafft.oneline.phy")
  if (file.exists(from)) {
    loci_to_concat[[i]] <- read.dna(from, format = "fasta") 
  } else {loci_to_concat[[i]] <- NA}
  #write.dna(f, to, format = "sequential")
  #file.copy(from = from, to = to)
}
names(loci_to_concat) <- NULL
loci_to_concat[[length(good_loci) + 1]] <- TRUE
names(loci_to_concat)[length(good_loci) + 1] <- "check.names"
loci_to_concat[[length(good_loci) + 2]] <- TRUE
names(loci_to_concat)[length(good_loci) + 2] <- "fill.with.gaps"
loci_to_concat <- loci_to_concat[!is.na(loci_to_concat)]
concat <- do.call(cbind.DNAbin, loci_to_concat)
write.dna(concat, file = "data/dory_combined/probe/nodups/concat.phy", format = "sequential")
write.FASTA(concat, file = "data/dory_combined/probe/nodups/concat.fasta")
