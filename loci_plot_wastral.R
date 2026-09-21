library(RevGadgets)
library(ggtree)
library(treeio)



tree1 <- treeio::read.newick("~/Downloads/wastral_50.island.tre")
tree1 <- drop.tip(tree1, c("D_concolor_J_Game_94_001_Rarotonga",
                           "D_concolor_Ranker_1939_Moorea",             
                           "D_concolor_LY_Kuo_195_Taiwan",                
                           "D_concolor_KR_Wood_10787_Ua_Huka"))
tree1$edge.length <- ifelse(is.nan(tree1$edge.length), 1.0, tree1$edge.length)
tree1 <- reroot(tree1, node.number = 24)
tree1$node.label <- ifelse(as.numeric(tree1$node.label) > 0.6, round(x = as.numeric(tree1$node.label), digits = 2), "")


tree2 <- treeio::read.newick("~/Downloads/wastral_50_with_subdecip.island.tre")
tree2 <- drop.tip(tree2, c("D_concolor_J_Game_94_001_Rarotonga",
                           "D_concolor_Ranker_1939_Moorea",             
                           "D_concolor_LY_Kuo_195_Taiwan",                
                           "D_concolor_KR_Wood_10787_Ua_Huka"))
tree2$edge.length <- ifelse(is.nan(tree2$edge.length), 1.0, tree2$edge.length)
tree2 <- reroot(tree2, node.number = 31)
tree2$node.label <- ifelse(as.numeric(tree2$node.label) > 0.6, round(x = as.numeric(tree2$node.label), digits = 2), "")

p1 <- ggtree(tree1) + geom_tiplab(size=2, offset = 1, color="#90909080") +
  coord_cartesian(xlim=c(0,8)) +
  #geom_text(aes(label=node), color="red", size=2) +
  geom_nodelab(size=2, color="#505050", nudge_x=0.15)
p2 <- ggtree(tree2)  + geom_tiplab(size=2, offset = -5, color="#90909080") +
  scale_x_reverse() +
  #geom_text(aes(label=node), color="red", size=2) +
  geom_nodelab(size=2, color="#505050", , nudge_x=-0.15)

p2a <- flip(p2, 34, 47) |> 
  flip(35, 40) |> 
  flip(1, 2) |> 
  flip(4, 5) |> 
  flip(7, 42) |> 
  flip(56, 59) |> 
  flip(57, 58) |> 
  flip(17, 18) |> 
  flip(15, 48)


p3 <- p1|p2a
ggsave("~/Desktop/dory_two_trees.pdf", p3, height = 6, width = 16, units = "in")


t <- reroot(t, node.number = 124)

p <- ggtree(t) + geom_tiplab(size=2, color="#90909080") +
  coord_cartesian(xlim=c(0,8)) +
  #geom_text(aes(label=node), color="red", size=2) +
  geom_nodelab(size=2, color="#505050") +
  coord_cartesian(xlim=c(0,10))
ggsave("~/Desktop/dory_compare.pdf", p, height = 6, width = 10, units = "in")
