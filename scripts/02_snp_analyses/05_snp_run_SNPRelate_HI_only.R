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


write_rds(x = tab, file = "output/snprelate_HI_only_tab.rds")
