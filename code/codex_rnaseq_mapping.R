##Map each CODEX cluster to the scRNAseq/snRNAseq cluster with maximum Pearson correlation

library(plyr)
library(dplyr)
library(purrr)
library(Seurat)
library(gplots)
library(ggplot2)
library(RColorBrewer)


#Read in CODEX clustering, average expression of each protein by cluster
codex_clustering <- read.table("codex76_normalized_clustering_with_size.csv", sep=",", header=TRUE)
codex_clustering <- codex_clustering[,c(1,7:length(codex_clustering[1,]))]
codex_cluster_ids <- codex_clustering$ClusterID
codex_clustering <- codex_clustering %>%  group_by(ClusterID) %>%  summarise_at(vars(colnames(codex_clustering)[2:length(colnames(codex_clustering))]), list(name = mean))
colnames(codex_clustering) <- as.character(map(strsplit(colnames(codex_clustering), split = "\\."), 1))
codex_clustering <- as.matrix(codex_clustering[,2:length(codex_clustering[1,])])
rownames(codex_clustering) <- sort(unique(codex_cluster_ids))
rownames(codex_clustering) <- as.numeric(rownames(codex_clustering)) - 424 #Renumber VorTex cluster ids

#Read in integrated clustering
load("scsn.integrated.013122.Robj")
all.genes <- rownames(tis.integrated)
tis.integrated <- ScaleData(tis.integrated, features = all.genes)

#Restrict CODEX data to genes shared with scRNAseq/snRNAseq data
codex_clustering <- dplyr::select(as.data.frame(codex_clustering), -matches(c("size", "COLIV",  "pan")))

#Restrict integrated scRNAseq/snRNAseq data to shared genes
genes_int <- c("AQP1","CD34", "UMOD", "KRT7", "ITGB4", "CD7", "PDPN", "FCGR3A", "CD38", "PDCD1", "CD4",
              "CD8A", "ICOS", "CD3D",  "AQP2", "ITGAX", "CR2", "CD19","SDC1", "PECAM1",
             "PTPRC", "MKI67", "SYNPO","THY1", "HLA-DRA", "CD68", "CD9" ,"CALB1" )
int_expr <- tis.integrated[["RNA"]]@scale.data[rownames(tis.integrated[["RNA"]]@scale.data) %in% genes_int,]
int_expr <- t(int_expr[genes_int,])
int_expr <- as.data.frame(int_expr)

#Average expression of each gene by cluster
int_avg_expr <- int_expr %>%  group_by(Idents(tis.integrated)) %>%  summarise_at(vars(colnames(int_expr)[1:length(colnames(int_expr))]), list(name = mean))
int_avg_expr <- as.matrix(int_avg_expr)
clusts <- int_avg_expr[,1]
rownames(int_avg_expr) <- clusts
int_avg_expr <- int_avg_expr[,2:length(colnames(int_avg_expr))]

#Convert average CODEX and scRNA/snRNAseq to matrices
int_avg_expr <- as.matrix(int_avg_expr)
genes <- colnames(int_avg_expr)
clusters <- rownames(int_avg_expr)
codex_clustering <- as.matrix(codex_clustering)
codex_clustering <- t(codex_clustering)
int_avg_expr <- matrix(as.numeric(int_avg_expr),    # Convert to numeric matrix
                  ncol = ncol(int_avg_expr))
colnames(int_avg_expr) <- genes
rownames(int_avg_expr) <- clusters
int_avg_expr <- t(int_avg_expr)

#Find max Pearson correlation between each CODEX cluster and each scRNAseq/snRNAseq cluster
max_correlations <- matrix(0,nrow=dim(codex_clustering)[2], ncol=dim(int_avg_expr)[2])

for (j in 1:dim(codex_clustering)[2]) {
tmp <- vector(length=dim(int_avg_expr)[2])
for (i in 1:dim(int_avg_expr)[2]) {tmp[i] <- cor(int_avg_expr[,i], codex_clustering[,j])}
max_correlations[j,which.max(tmp)] <- 1
}

rownames(max_correlations) <- colnames(codex_clustering)
colnames(max_correlations) <- colnames(int_avg_expr)

col.order = c("EC-PTC", "EC-AEA", "EC-GC", "VSMC/P", "MC", "POD", "PT/PEC", "PT", "DTL", "ATL/TAL", "TAL", "DCT", "CNT", "PC", "tPC-IC", "IC", "FIB", "T", "B", "MYL", "NKT/NKC", "MAC")

col.side.colors = c("#DE7F7D", "#DE7F7D", "#DE7F7D", "#DE7F7D", "#429898", "#429898","#F19A48", "#F19A48", "#3481B9", "#3481B9", "#3481B9", "#F19A48", "#9498FA", "#9498FA", "#9498FA", "#9498FA", "#82878A", "#867E57", "#867E57", "#867E57", "#867E57", "#867E57")

col.side.colors = col.side.colors[col.order %in% colnames(max_correlations)]
col.order = col.order[col.order %in% colnames(max_correlations)]

row.order = c("436",  "449", "451", "452", "448", "425",  "437","427", "430",  "438" ,"443", "447","455", "435", "428" ,"431" ,"441","444","440","454","434","426", "453","445", "429", "450" ,"432", "446" , 
 "433" ,"439","442")

row.order = as.numeric(row.order) - 424

max_correlations = max_correlations[-c(9,15,18),] #remove background clusters
row.order = row.order[as.character(row.order) %in% rownames(max_correlations)]
row.order <- as.character(row.order)

max_correlations = max_correlations[row.order, col.order]

pdf("220211_int_pearcorr.pdf")
heatmap.2(max_correlations, scale="none", trace="none", dendrogram="none", colsep=c(1:30), rowsep=c(1:30), sepwidth=c(0.01,0.01), adjCol = c(1,0.5), col=brewer.pal(n=8, name="Reds"), Rowv=FALSE, Colv=FALSE, ColSideColors=col.side.colors)
dev.off()

####Plot heatmap with side-by-side codex and single-cell data
for (i in rownames(max_correlations)) {
    filename=paste("heatmaps/", i, ".pdf", sep="")
pdf(filename)
    my_matrix <- rbind(codex_clustering[,i], int_avg_expr[,colnames(max_correlations)[which.max(max_correlations[i,])]])
    rownames(my_matrix) <- c("CODEX", "scRNA-seq")
        my_data_frame <- as.data.frame.table(my_matrix)
    names(my_data_frame) <- c("Tech", "gene", "expr")
print(ggplot(data=my_data_frame, aes(Tech, gene, fill=expr)) + geom_tile() + coord_fixed() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=14)) +
     theme(axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
          axis.title.x=element_blank()) +
scale_fill_stepsn(colors=rev(heat.colors(n=12))) +
guides(fill=guide_legend(title="expression")))
dev.off()
    print(i)
}