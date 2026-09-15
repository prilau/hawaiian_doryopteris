# compare number of duplicated taxa by loci and by taxon
library(ape)
library(tidyverse)

setwd("~/Box/carrie2018/project_Doryopteris/")
file_path <- "data/dory_combined/probe/"
#file_path <- "data/dory_set1/Doryopteris_ProbeOnly/"

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

name_dict <- read.csv("metadata/names_to_extr_num.csv")
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

bad_taxa <- v_row[which(v_row > .1)]

bad_loci <- v_col[which(v_col > 0.1)]

# plotting
 cols <- rev(c("#d73027",
               "#d73027",
               "#d73027",
               "#d73027",
               "#d73027",
               "#d73027",
               "#d73027",
               "#d73027",
               "#4575b4",
               "#FFFFFF"))
 
 # plot histogram of v_col 
 hist(v_col, breaks = 100, col = "white", 
      main = "v_col metric", 
      xlab = "v_col", 
      ylab = "Number of columns")
 # log version
 hist(log(v_col), breaks = 100, col = "white", 
      main = "v_col metric", 
      xlab = "log of v_col", 
      ylab = "Number of columns")
 abline(v = log(0.1), col = "blue")

 
 pdf("figures/combined_probe_dup_taxa_unsorted.pdf", width = 20, height = 20)
 image(t(res_subset[rev(1:nrow(res_subset)),]), 
       col = cols,axes = FALSE)
 dev.off()
 
 pdf("figures/combined_probe_dup_taxa_colsort.pdf", width = 20, height = 20)
 image(t(res_subset_colsort[rev(1:nrow(res_subset_colsort)),]), 
       col = cols,axes = FALSE)
 # add stars for weird small Hawaii clade
 rows <- which(rownames(res_subset) %in% c("XZ038", "XZ021", "XZ043","XZ034","XZ026"))
 text(x = rep(.99, times = length(rows)), 
      y = 1-rows/nrow(res_subset), 
      labels = rep("*", times = length(rows)),
      col = "red")
 #rect(xleft = 0, xright = length(bad_loci)/ncol(res_subset_colsort), 
 #     ybottom = 0, ytop = 1, 
 #     border = "black", lwd = 5)
 dev.off()
 
 pdf("figures/combined_probe_dup_taxa_rowsort.pdf", width = 20, height = 20)
 image(t(res_subset_rowsort[rev(1:nrow(res_subset_rowsort)),]), 
       col = cols,axes = FALSE)
 #rect(xleft = 0, xright = 1, 
 #     ybottom = 1-length(bad_taxa)/nrow(res_subset_rowsort), ytop = 1, 
 #     border = "black", lwd = 5)
 dev.off()



### copy bad loci into a new directory for building gene trees
dir.create(paste0(file_path, "/bad_loci"))
for (i in 1:length(bad_loci)) {
  from = paste0(file_path,names(bad_loci)[i],".afterMerge.mafft.oneline.fasta")
  to = paste0(file_path, "bad_loci/",names(bad_loci)[i],".afterMerge.mafft.oneline.phy")
  f <- read.dna(from, format = "fasta")
  write.dna(f, to, format = "sequential")
  #file.copy(from = from, to = to)
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
