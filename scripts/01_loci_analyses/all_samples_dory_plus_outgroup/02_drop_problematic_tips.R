library(ape)

setwd("/Users/priscillalau/Desktop/ASSIM/rLab/Doryopteris_project/data/fullphased/")

files <- list.files()
files_phy <- files[grep(".phy", files)]
files_phy <- files_phy[grep(".prune", files_phy, invert = T)]
files_phy <- files_phy[grep(".pl", files_phy, invert = T)]
alns_phy <- list()

#exclude blank files
empties <- c()
for (i in 1:length(files_phy)) {
  if (length(alns_phy[[i]]$seq) == 0) {
    empties <- c(empties, i)
  }
}
files_phy <- files_phy[-empties]

#read in phy alignments
for (i in 1:length(files_phy)) {
  alns <- read.dna(file = files_phy[i], format = "sequential")
  alns_phy[[i]] <- as.alignment(alns)
  names(alns_phy)[i] <- strsplit(files_phy[i],split = "\\.")[[1]][[1]]
}

#change $nam to extr_num and remove XZ050
for (i in 1:length(alns_phy)) {
  for (j in 1:length(alns_phy[[i]]$nam)) {
    if (is.na(alns_phy[[i]]$nam[j]) == TRUE) {
      name <- alns_phy[[i]]$nam[j]
      print(paste0("name: ", alns_phy[[i]]$nam[j],
                   "\ni: ", i,
                   "\nj: ", j))
    } else {
      name <- paste0("XZ", strsplit(alns_phy[[i]]$nam[j],
                                    split = "_XZ")[[1]][[2]])
      alns_phy[[i]]$nam[j] <- strsplit(name,
                                       split = "\t")[[1]][[1]]
    }
    if (grepl("XZ050", name) == TRUE | grepl("XZ066", name) == TRUE) {
      alns_phy[[i]]$nam[j] <- NA
      alns_phy[[i]]$seq[j] <- NA
      alns_phy[[i]]$nb <- alns_phy[[i]]$nb - 1
    }
  }
  alns_phy[[i]]$seq <- alns_phy[[i]]$seq[!is.na(alns_phy[[i]]$seq)]
  alns_phy[[i]]$nam <- alns_phy[[i]]$nam[!is.na(alns_phy[[i]]$nam)]
}

#exclude blank files again
empties2 <- c()
for (i in 1:length(alns_phy)) {
  if (length(alns_phy[[i]]$seq) == 0) {
    empties2 <- c(empties2, i)
  }
}
alns_phy <- alns_phy[-empties2]

#save as .phy
for (i in 1:length(alns_phy)) {
    bin <- as.DNAbin(alns_phy[[i]])
    # fix here - you need names(alns_phy)[i] not just alns_phy
    write.dna(bin, file = paste0(names(alns_phy)[i], ".new.phy"), format = "sequential")
}
