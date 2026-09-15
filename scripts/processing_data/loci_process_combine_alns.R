setwd("~/Box/carrie2018/project_Doryopteris/")
library(rMSA)

# function to rename seqs to extraction number and copy number
rename <- function(s) {paste0("XZ",unlist(strsplit(strsplit(s, split = "XZ")[[1]][[2]], split = "_sp")))}

### arguments 
#merge_what <- "full"
merge_what <- "probe"
#file_stem <- ".afterMerge.mafft.oneline.fasta"
file_stem <- ".nodups.phy"

if (merge_what == "probe") { 
  f1_dir = "Doryopteris_ProbeOnly/"
  f2_dir = "Doryopteris2_Probe/"
} else if (merge_what == "full") {
  f1_dir = "Doryopteris_Full/"
  f2_dir = "Doryopteris2_Full/"
  }

outdir <- paste0("data/dory_combined/", merge_what, "/nodups")
dir.create(outdir)

loci <- paste0("L", 1:450)

for (L in loci) {
  
  # read in alignment from first set 
  f1_name <- paste0("data/dory_set1/", 
                    f1_dir,
                    L, 
                    file_stem)
  f1_exists <- file.exists(f1_name)
  if (f1_exists) {
    if (grepl("fasta", file_stem)) {
      f1 <- readDNAStringSet(f1_name)
    } else { 
      f1 <- readDNAMultipleAlignment(f1_name, format = "phylip")
      f1 <- as(f1,"DNAStringSet")
      }
    if (length(f1) == 0) {
      rm(f1)
      f1_exists <- FALSE
      }
    }
  # read in alignment from second set 
  f2_name <- paste0("data/dory_set2/", 
                    f2_dir, 
                    L, 
                    file_stem)
  f2_exists <- file.exists(f2_name)
  if (f2_exists) {
    if (grepl("fasta", file_stem)) {
      f2 <- readDNAStringSet(f2_name)
    } else { 
      f2 <- readDNAMultipleAlignment(f2_name, format = "phylip")
      f2 <- as(f2,"DNAStringSet")
    }
    if (length(f2) == 0) {
      rm(f2)
      f1_exists <- FALSE
    }
  }
  
  # combine alignments, realign, and save
  if (f1_exists & f2_exists) {
    #print(paste0(L, ": combining both alignments"))
    combined <- c(f1,f2)
    combined@ranges@NAMES <- unlist(lapply(combined@ranges@NAMES, FUN = rename))
    # keep mafft quiet 
    combined_aln <- as(rMSA::mafft(x = combined), "DNAStringSet") 
    writeXStringSet(combined_aln, 
                    filepath = paste0(outdir,
                                      "/", 
                                      L, 
                                      ".nodups.fasta"), 
                    append=FALSE,
                    compress=FALSE,
                    compression_level=NA, 
                    format="fasta") 
  }
  
  # save just existing file
  if (f1_exists & !f2_exists) {
    print(paste0(L, ": only f1"))
    f1@ranges@NAMES <- unlist(lapply(f1@ranges@NAMES, FUN = rename))
    writeXStringSet(combined_aln, 
                    filepath = paste0(outdir,
                                      "/", 
                                      L, 
                                      ".nodups.fasta"), 
                    append=FALSE,
                    compress=FALSE,
                    compression_level=NA, 
                    format="fasta")  
  }
  
  # save just existing file
  if (!f1_exists & f2_exists) {
    print(paste0(L, ": only f2"))
    f2@ranges@NAMES <- unlist(lapply(f2@ranges@NAMES, FUN = rename))
    writeXStringSet(combined_aln, 
                    filepath = paste0(outdir,
                                      "/", 
                                      L, 
                                      ".nodups.fasta"), 
                    append=FALSE,
                    compress=FALSE,
                    compression_level=NA, 
                    format="fasta") 
  }
  
  # do nothing if neither exists
  if (!f1_exists & !f2_exists) {
    print(paste0(L, ": skipping, neither exists"))
  }
}
