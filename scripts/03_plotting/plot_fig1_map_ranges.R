# map Alstroemeriaceae localities
library(tidyverse)
library(rgbif)
library(sp)
library(maptools)
library(maps)
setwd("~/Box/project_Doryopteris/figures/")


# get and clean occurrence records for all Alstroemeriaceae
occ_dory <- occ_data(scientificName = "Doryopteris", country = "US", limit=200000)

# selecting Hawaiian species
occ_hawaii <- occ_dory$data[c(which(occ_dory$data$species == "Doryopteris angelica"),
                              which(occ_dory$data$species == "Doryopteris decipiens"),
                              which(occ_dory$data$species == "Doryopteris decora"),
                              which(occ_dory$data$species == "Doryopteris subdecipiens"),
                              which(occ_dory$data$species == "Doryopteris takeuchii"))
                             ,]

# plot all occurrences on map, color by genus

mycolors <- c("#ddbea9",
            "#a5a58d",
            "#bc6c25",
            "#34623f",
            "#550527")

names(mycolors) <- c("Doryopteris angelica",
                     "Doryopteris decipiens",
                     "Doryopteris decora",
                     "Doryopteris subdecipiens",
                     "Doryopteris takeuchii")
occ_hawaii$colors <- occ_hawaii$species

occ_hawaii$colors <- mycolors[occ_hawaii$colors]
occ_hawaii <- occ_hawaii[!is.na(occ_hawaii$decimalLongitude), ]
#coordinates(occ_hawaii) <- ~decimalLongitude + decimalLatitude


data(world2MapEnv)
png(file = "~/Desktop/map.png",
    width = 12, height = 8, units = "in", res = 300)
map('world', wrap=c(0,360))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "aliceblue")
map('world', wrap=c(0,360), fill = TRUE, add = T, col = "ivory")
points(occ_records, pch=16, col = paste0(occ_records$colors, "50"), cex = 0.7)
#legend(x = 200, y = 10, legend = names(mycolors),
#       col = mycolors, cex = 1.1, pch = 16, text.font = rep(3, times = 4))
dev.off()