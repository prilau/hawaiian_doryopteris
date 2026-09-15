setwd("PCA_HI")

#load necessary packages
library(gdsfmt)
library(SNPRelate)
library(VariantAnnotation)
library(scales)


#call VCF file
vcf.fn <- "HI_only_linked_SNPs.recode.names_fixed.vcf"

#convert to GDS format for SNPRelate analysis
showfile.gds(closeall = TRUE)
snpgdsVCF2GDS(vcf.fn, "test.gds", method="biallelic.only", ignore.chr.prefix = "L")

#summarize data
snpgdsSummary("test.gds")
genofile <- snpgdsOpen("test.gds", readonly=FALSE)

#LD prune
set.seed(1000)
snpset <- snpgdsLDpruning(genofile, ld.threshold=10, maf=NaN, missing.rate=NaN, slide.max.n=5000,start.pos="random",autosome.only=FALSE)
names(snpset)
snpset.id <- unlist(snpset)

#run initial PCA
pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=2)


# determine the variance proportion (%) explained
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))
rounded <- head(round(pc.percent, 2))

# make a data.frame
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],    # the third eigenvector
                  EV4 = pca$eigenvect[,4],    # the fourth eigenvector
                  stringsAsFactors = FALSE)
head(tab)

#population code for color coding
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
popmap <- read.table("popmap.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
pop_code <- popmap$pop[match(pca$sample.id, popmap$sample.id)]

# sanity check: every sample got a code, counts look right
stopifnot(!anyNA(pop_code))
table(pop_code)

# species is the 2nd underscore-delimited token of the sample id (D_<species>_<collector>_<number>)
species_code <- sub("^D_([^_]+)_.*$", "\\1", pca$sample.id)

tab <- data.frame(sample.id = pca$sample.id,
                  pop = factor(pop_code),          # island code, kept for later use
                  species = factor(species_code),
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],    # the third eigenvector
                  EV4 = pca$eigenvect[,4],    # the fourth eigenvector
stringsAsFactors = FALSE)
head(tab)
table(tab$species)


kauai_col <- "#2ed1db"
maui_col <- "#fc948e"
oahu_col <- "#a80f3d"
nothawaii_col <- "#b9cff9"

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
pdf("PCA_unique_colors_PC1vsPC2.pdf",10,10)
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

#plot PC 2 vs PC 3
pdf("PCA_unique_colors_PC2vsPC3.pdf",10,10)
par(mar=c(6,6,1.5,2))
plot(tab$EV2,
     tab$EV3,
     col = "black",
     bg = alpha(pop_col[as.character(tab$pop)], 0.7),
     cex = 3,
     pch = sp_pch[as.character(tab$species)],
     cex.lab = 2,
     cex.main = 2,
     cex.axis = 1.5,
     main = "",
     xlab=paste0("eigenvector 2 (",rounded[2],"%)"),
     ylab=paste0("eigenvector 3 (",rounded[3],"%)"))
label_points("EV2", "EV3", to_label, pos = 4, cex = 1.4)
dev.off()

#plot PC 3 vs PC 4
pdf("PCA_unique_colors_PC3vsPC4.pdf",10,10)
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
