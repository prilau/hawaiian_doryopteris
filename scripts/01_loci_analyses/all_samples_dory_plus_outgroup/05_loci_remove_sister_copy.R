library(ape)
setwd("~/Box/carrie2018/project_Doryopteris/Doryopteris_FullPhased/strictplus_genetree/")

files <- list.files()
csv <- read.csv("~/Box/carrie2018/project_Doryopteris/metadata/names_to_extr_num.csv")
names <- csv$ex
num_copies1 <- matrix(nrow = length(names), ncol = length(files))
num_copies2 <- matrix(nrow = length(names), ncol = length(files))

for (i in 1: length(files)) {
  t <- read.tree(file = files[i])
  for (j in 1:length(names)) {
    n <- grep(names[j], t$tip.label)
    if (length(n) == 0) {
      num_copies1[j, i] <- 0
      num_copies2[j, i] <- 0
    } else if (length(n) ==1) {
      num_copies1[j, i] <- 1
      num_copies2[j, i] <- 1
    } else if (length(n) > 1) {
      mono <- is.monophyletic(phy = t, tips = n)
      if (mono == TRUE) {
        num_copies1[j, i] <- 1
        num_copies2[j, i] <- 2
      } else {
        num_copies1[j, i] <- 2
        num_copies2[j, i] <- 2
      }
    }
  }
}

rownames(num_copies1) <- names
rownames(num_copies2) <- names
colnames(num_copies1) <- unlist(strsplit(unlist(strsplit(files, split = "strictplus_")), split = ".fasta.tre"))
colnames(num_copies2) <- unlist(strsplit(unlist(strsplit(files, split = "strictplus_")), split = ".fasta.tre"))

# pdf(file = "~/Desktop/plot.pdf", height = 10, width = 5)
# par(mfrow = c(2, 1))
# image(t(num_copies1), col = c("white", "green", "red"), axes = TRUE)
# image(t(num_copies2), col = c("white", "green", "red"), axes = TRUE)
# dev.off()

cols <- c("grey","blue","orange")


pdf("~/Desktop/heatmap.pdf", width = 40, height = 80)
tt <- t(num_copies1)
par(oma=c(0,0,0,0)+1, mfrow = c(2,1))
image(x = 1:dim(num_copies1)[2], y=1:dim(num_copies1)[1], z = tt, xlab=NA, ylab=NA, xaxt="n", yaxt="n", col=cols)
axis(2, at=dim(num_copies1)[1]:1, labels = rownames(num_copies1), las=2, lwd.tick=0.1, lwd=0)
axis(3, at=1:dim(num_copies1)[2], labels = colnames(num_copies1), las=2, lwd.tick=0.1, lwd=0)
legend("bottomleft", legend=0:2, fill=cols, cex=3, horiz=TRUE, bty="n", inset=c(0,1.01), xpd=NA)

tt <- t(num_copies2)
image(x = 1:dim(num_copies2)[2], y=1:dim(num_copies2)[1], z = tt, xlab=NA, ylab=NA, xaxt="n", yaxt="n", col=cols)
axis(2, at=dim(num_copies2)[1]:1, labels = rownames(num_copies2), las=2, lwd.tick=0.1, lwd=0)
axis(3, at=1:dim(num_copies2)[2], labels = colnames(num_copies2), las=2, lwd.tick=0.1, lwd=0)
legend("bottomleft", legend=0:2, fill=cols, cex=3, horiz=TRUE, bty="n", inset=c(0,1.01), xpd=NA)
dev.off()



