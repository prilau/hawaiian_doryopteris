setwd("/Users/priscillalau/Desktop/ASSIM/rLab/Doryopteris_project/data/trimmed/trees/")
library(ape)


tre <- c()
tree_all <- list.files()
tree_all <- tree_all[grep("\\.tre", tree_all)]

for (i in 1:length(tree_all)){
  tre[[i]] <- read.tree(tree_all[i])
}

outgroup <- c("L_crenulans_Brazil",
              "L_lomariacea_Brazil",
              "O_gleichenioides_Brazil",
              "O_pinnata_Brazil",
              "O_riedelii_Brazil",
              "L_quinquelobatum_Brazil",
              "L_paradoxa_Brazil")

needs_rooting <- numeric()
for (i in 1:length(tre)){
  has_outgroup <- any(outgroup %in% tre[[i]]$tip.label) == T
  if (has_outgroup) {
         outgroup_custom <- outgroup[outgroup %in% tre[[i]]$tip.label]
         if (length(outgroup_custom) == 1) {
           num <- which(tre[[i]]$tip.label == outgroup_custom)
         } else {
           num <- getMRCA(tre[[i]], outgroup_custom)
         }
         midpoint <- 0.5 * tre[[i]]$edge.length[which(tre[[i]]$edge[, 2] == num)]
         t_rooted <- phytools::reroot(tre[[i]], num, position = midpoint)
         write.tree(t_rooted, file = paste0("rooted/", "rooted_", tree_all[i]))
  } else {
    needs_rooting <- append(needs_rooting, i)
  }
}

needs_rooting <- tre[needs_rooting]



#rooting manually
tree_all <- list.files()
tree_all <- tree_all[grep("L397a", tree_all)]
tre <- read.tree(tree_all)
outgroup <- c("D_concolor_Brazil_1",
              "D_concolor_Bolivia")

num <- which(tre$tip.label == outgroup)
#or
num <- getMRCA(tre, outgroup)

midpoint <- 0.5 * tre$edge.length[which(tre$edge[, 2] == num)]
t_rooted <- phytools::reroot(tre, num, position = midpoint)
write.tree(t_rooted, file = paste0("rooted/", "rooted_", tree_all))
