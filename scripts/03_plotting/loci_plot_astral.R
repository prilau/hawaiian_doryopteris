# plot full astral tree
library(ape)
library(ggnewscale)
library(ggtree)
library(RevGadgets)
library(treeio)

##### read in and edit associated data ##### 
dat <- read.csv("metadata/modified_for_publication/combined_sample_for_github.csv") |> 
  filter(coll_num != "KR Wood 9071")

# edit localities
for (i in 1:length(dat$locality)) {
  if (dat$locality[i] == "Kahoolawe") {
    dat$locality[i] <- "Maui Nui"
  } else if (dat$locality[i] == "Maui") {
    dat$locality[i] <- "Maui Nui"
  } else if (dat$locality[i] == "Oahu ") {
    dat$locality[i] <- "Oahu"
  } 
}

# fix spellings
dat$locality[which(dat$locality == "Kauai")] <- "Kaua`i"
dat$locality[which(dat$locality == "Oahu")] <- "O`ahu"
dat$locality[which(dat$locality == "Hawai'i Island")] <- "Hawai`i Island"

# keep old names that match with tree tip labels
dat$old_name <- dat$name

dat$new_ID <- gsub("_", ". ", dat$new_ID)

#dat$index <- dat$name
# replace "Hawaii" in dat$name with collector number
for (i in 1:length(dat$name)){
  dat$name[i] <- paste0(dat$new_ID[i], " ", dat$coll_num[i])
  
  while(grepl("  ", dat$name[i])){
    dat$name[i] <- gsub("  ", " ", dat$name[i])
  }
  
}

# assign clade
dat$clade <- ifelse(dat$locality == "Kaua`i" | 
                       dat$locality == "O`ahu" | 
                       dat$locality == "Maui Nui" | 
                       dat$locality == "Hawai`i Island",
                     "hawaii",
                    ifelse(dat$new_ID == "D. concolor", "concolor",
                           ifelse(grepl("D.", dat$new_ID), "other dory", "not dory")))

##### read tree #####
tree <- read.tree("output/dory_HI_BS10.tre")

# reroot tree
outgroup <- c("L_crenulans_Brazil",
              "L_lomariacea_Brazil",
              "O_gleichenioides_Brazil",
              "O_pinnata_Brazil",
              "O_riedelii_Brazil",
              "L_quinquelobatum_Brazil",
              "L_paradoxa_Brazil")

outgroup <- outgroup[outgroup %in% tree$tip.label]
if (length(outgroup) == 1) {
  num <- which(tree$tip.label == outgroup)
} else {
  num <- getMRCA(tree, outgroup)
}
midpoint <- 0.5 * tree$edge.length[which(tree$edge[, 2] == num)]
t <- phytools::reroot(tree, num, position = midpoint)

# rescale tree height to 1
t$edge.length <- t$edge.length / max(node.depth.edgelength(t))
write.tree(t, file = "rooted_dory_BS10_BL1.tre")

t <- drop.tip(t, tip="D_decipiens_Hawaii_17")

# change/add locality to tip labels as well
for (i in 1:length(t$tip.label)){
  t$tip.label[i] <- dat$name[which(dat$old_name == t$tip.label[i])]
}

##### combine tree and data, edit tree #####
# get support values
support_values <- t$node.label
support_values <- c(rep(NA, times = length(t$tip.label)), support_values)
support_values <- as.numeric(support_values)
support_values[support_values < .6] <- NA

# make treedata object
dat <- dplyr::left_join(as_tibble(t), 
                        dat, 
                        by = c("label" = "name"))
t_dat <- as.treedata(dat)

# add support values 
t_dat@data$support <- support_values

# set line type to make terminal branches dashed and bold well-supported branches
support <- t_dat@data$support
t_dat@data$linewidth <- ifelse(is.na(support), "bad", "good")
ntips <- length(t_dat@phylo$tip.label)
t_dat@data$linetype <- c(rep("dashed", ntips),
                         rep("solid", nrow(t_dat@data) - ntips))






##### plot tree #####
colors <- c("hawaii" = "#FCAB10",
            "concolor" = "#acc8fa",
            "other dory" = "#214478",
            "Not Dory" = "#636363")

t_plot <- ggtree(t_dat, aes(linetype = linetype, size = linewidth)) +
  geom_tiplab(aes(color = clade), 
              size = 4, 
              fontface = "bold",
              nudge_x = .25) +
  geom_text2(aes(label=support), nudge_x = 0.3, nudge_y = 0, size = 2, hjust = 1, color="grey30") +
  scale_color_manual(values = colors) +
  scale_size_discrete("level", range=c(.5,1)) + 
  scale_linetype_identity() +
  scale_shape_manual(values = c(NA,23)) +
  xlim(0, 15) +
  theme(legend.position = 'none',
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent")) # get rid of legend panel bg)


t_plot

ggsave("~/figures/supp_fig_astral2_tree_full.pdf", t_plot, width = 8, height = 10, units = "in")

