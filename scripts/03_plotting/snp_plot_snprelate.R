library(dplyr)

################
# plot HI only
################

tab <- read_rds("output/snprelate_HI_only_tab.rds")

table(tab$species)

kauai_col <- "#14d2dc"
maui_col <- "#fe938c"
oahu_col <- "#aa0a3c"
nothawaii_col <- "#acc8fa"

# one fill-able plotting symbol per species
sp_lv <- levels(tab$species)
sp_pch <- c(angelica = 21, decipiens = 22, decora = 23, subdecipiens = 24, takeuchii = 25)
stopifnot(all(sp_lv %in% names(sp_pch)))

# island (pop) palette + readable labels, keyed by the pop code
pop_col <- c("1" = oahu_col, "2" = kauai_col, "3" = maui_col)
pop_lab <- c("1" = "O'ahu",  "2" = "Kaua'i", "3" = "Maui Nui")
stopifnot(all(levels(tab$pop) %in% names(pop_col)))
pop_lv <- levels(tab$pop)

#label select samples in the PCA plot

# do this once, interactively (quartz() gives a clickable window; VS Code's viewer does not)
#quartz()
#plot(tab$EV1, tab$EV2, pch = 21, bg = pop_col[as.character(tab$pop)])
#picked <- identify(tab$EV1, tab$EV2, labels = tab$sample.id)
#tab$sample.id[picked]
#dev.off()
to_label <- c("D_decipiens_CM_Tribble_135", "D_decipiens_CM_Tribble_143",
               "D_decora_DD_Palmer_585", "D_subdecipiens_KR_Wood_18028")

label_points <- function(xvar, yvar, ids, ...) {
  idx <- match(ids, tab$sample.id)
  if (anyNA(idx)) warning("not in tab: ", paste(ids[is.na(idx)], collapse = ", "))
  text(tab[[xvar]][idx], tab[[yvar]][idx],
       labels = tab$sample.id[idx], xpd = TRUE, ...)
}

#plot PC 1 vs PC 2
pdf("PCA_HI_only_PC1vsPC2.pdf",10,10)
par(mar=c(6,6,1.5,2))
plot(tab$EV1,
     tab$EV2,
     col = "black",
     bg = alpha(pop_col[as.character(tab$pop)], 0.7),
     cex = 3,
     pch = sp_pch[as.character(tab$species)],
     cex.lab = 2,
     cex.main = 2,
     cex.axis = 1.5,
     main = "Hawaiian Clade",
     xlab = paste0("eigenvector 1 (",rounded[1],"%)"),
     ylab = paste0("eigenvector 2 (",rounded[2],"%)"))
label_points("EV1", "EV2", to_label, pos = 4, cex = 1.4)
dev.off()


#plot PC 3 vs PC 4
pdf("PCA_HI_only_PC3vsPC4.pdf",10,10)
par(mar=c(6,6,1.5,2))
plot(tab$EV3,
     tab$EV4,
     col = "black",
     bg = alpha(pop_col[as.character(tab$pop)], 0.7),
     cex = 3,
     pch = sp_pch[as.character(tab$species)],
     cex.lab = 2,
     cex.main = 2,
     cex.axis = 1.5,
     main = "",
     xlab=paste0("eigenvector 3 (",rounded[3],"%)"),
     ylab=paste0("eigenvector 4 (",rounded[4],"%)"))
label_points("EV3", "EV4", to_label, pos = 4, cex = 1.4)
dev.off()



###################
# plot HI+outgroup
###################
tab <- read_rds("output/snprelate_HI_plus_outgroup_tab.rds")

table(tab$species)

kauai_col <- "#14d2dc"
maui_col <- "#fe938c"
oahu_col <- "#aa0a3c"
nothawaii_col <- "#acc8fa"

# shapes match the Hawaii-only figures; concolor recycles the diamond (pch 23)
sp_lv <- levels(tab$species)
sp_pch <- c(angelica = 21, decipiens = 22, decora = 23, subdecipiens = 24, takeuchii = 25,
            concolor = 23)
stopifnot(all(sp_lv %in% names(sp_pch)))

# island (pop) palette + readable labels, keyed by the pop code
pop_col <- c("1" = oahu_col, "2" = kauai_col, "3" = maui_col, "5" = nothawaii_col)
pop_lab <- c("1" = "O'ahu",  "2" = "Kaua'i", "3" = "Maui Nui", "5" = "Not Hawai'i")
stopifnot(all(levels(tab$pop) %in% names(pop_col)))
pop_lv <- levels(tab$pop)


#label select samples in the PCA plot

# do this once
# quartz()
# plot(tab$EV1, tab$EV2, pch = 21, bg = pop_col[as.character(tab$pop)])
# picked <- identify(tab$EV1, tab$EV2, labels = tab$sample.id)
# tab$sample.id[picked]
# dev.off()

to_label <- c("D_concolor_LY_Kuo_195", "D_decipiens_CM_Tribble_135",
               "D_decipiens_CM_Tribble_143")

label_points <- function(xvar, yvar, ids, ...) {
  idx <- match(ids, tab$sample.id)
  if (anyNA(idx)) warning("not in tab: ", paste(ids[is.na(idx)], collapse = ", "))
  text(tab[[xvar]][idx], tab[[yvar]][idx],
       labels = tab$sample.id[idx], xpd = TRUE, ...)
}

#plot PC 1 vs PC 2
pdf("PCA_HI_plus_outgroup_PC1vsPC2.pdf",10,10)
par(mar=c(6,6,1.5,2))
plot(tab$EV1,
     tab$EV2,
     col = "black",
     bg = alpha(pop_col[as.character(tab$pop)], 0.7),
     cex = 3,
     pch = sp_pch[as.character(tab$species)],
     cex.lab = 2,
     cex.main = 2,
     cex.axis = 1.5,
     main = "Hawaiian Clade + Outgroup",
     xlab = paste0("eigenvector 1 (",rounded[1],"%)"),
     ylab = paste0("eigenvector 2 (",rounded[2],"%)"))
# combined legend: Island block on top, Species block directly below; matched symbol + font size
# leg_cex   <- 1.5
# leg_ptcex <- 2
# isl <- legend("topright", title = "Island", bty = "n",
#               legend = pop_lab[pop_lv],
#               pch = 22, col = "black", pt.bg = pop_col[pop_lv],
#               cex = leg_cex, pt.cex = leg_ptcex)
# legend(x = isl$rect$left, y = isl$rect$top - isl$rect$h, title = "Species", bty = "n",
#        legend = parse(text = paste0("italic('D. ", sp_lv, "')")),
#        pch = sp_pch[sp_lv], col = "black", pt.bg = "grey70",
#        cex = leg_cex, pt.cex = leg_ptcex)
label_points("EV1", "EV2", to_label, pos = 4, cex = 1.4)
dev.off()

#plot PC 3 vs PC 4
pdf("PCA_HI_plus_outgroup_PC3vsPC4.pdf",10,10)
par(mar=c(6,6,1.5,2))
plot(tab$EV3,
     tab$EV4,
     col = "black",
     bg = alpha(pop_col[as.character(tab$pop)], 0.7),
     cex = 3,
     pch = sp_pch[as.character(tab$species)],
     cex.lab = 2,
     cex.main = 2,
     cex.axis = 1.5,
     main = "",
     xlab=paste0("eigenvector 3 (",rounded[3],"%)"),
     ylab=paste0("eigenvector 4 (",rounded[4],"%)"))
label_points("EV3", "EV4", to_label, pos = 4, cex = 1.4)
dev.off()
