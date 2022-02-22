#Plot heatmap of clustered CODEX data

library(purrr)
library(gplots)

clustering <- read.table("codex76_normalized_clustering_with_size.csv", sep=",", header=TRUE)
clustering <- clustering[,c(1,7:length(clustering[1,]))]
library(dplyr)
tmp <- clustering %>%  group_by(ClusterID) %>%  summarise_at(vars(colnames(clustering)[2:length(colnames(clustering))]), list(name = mean))
colnames(tmp) <- as.character(map(strsplit(colnames(tmp), split = "\\."), 1))
tmp <- as.matrix(tmp[,2:length(tmp[1,])])
rownames(tmp) <- sort(unique(clustering$ClusterID))

rownames(tmp) <- as.numeric(rownames(tmp))-424 #Renumber clusters starting with 1
par(cex.main=1, cex.lab=1, cex.axis=1)
heatmap.2(as.matrix(tmp), density.info="none", trace="none", scale="row", col=rev(heat.colors(12)))
