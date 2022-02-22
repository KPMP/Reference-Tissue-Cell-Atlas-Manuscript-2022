#To integrate single-cell sequencing, single-nucleus sequencing, and LMD datasets, we first identify genes measured both in the LMD dataset and in the corresponding single-cell transcriptomic dataset. From this set of shared genes, we restrict to the subset of genes showing variable expression in the single-cell dataset.  We then compute the Pearson correlation between each cell in the scaled single-cell dataset and each the average ratio LCM expression profile for each subsegment. We assign each cell to the LCM subsegment with the highest correlation value. To evaluate the assignments, we examine the normalized distribution of cells assigned to each LCM subsegment within each single-cell cluster. We find that there is strong concordance across the datasets, with the majority of cells from each single-cell cluster assigned to the corresponding LCM segment (for example, proximal tubule cells are assigned to the proximal tubule subsegment, while podocytes are assigned to the glomerular subsegment).

#R version 4.1.0
library(Seurat)
library(dplyr)
library(tidyr) 
load("scsn.integrated.013122.Robj") #clustered integrated sc/sn dataset

#Find variable features, scale
tis.integrated <- FindVariableFeatures(tis.integrated, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(tis.integrated)
tis.integrated <- ScaleData(tis.integrated, features = all.genes)


#Read in and reformat LCM dataset
lcm <- read.table("Sub-segmenal_RNASeq_analysis_ttest.txt", header=TRUE)
lcm_table <- cbind(as.character(lcm$Segment), as.character(lcm$Gene), lcm$Log2_fc)
lcm_table <- as.data.frame(lcm_table)
colnames(lcm_table) <- c("Segment", "Gene", "Log2_fc")
lcm_table$Log2_fc <- as.numeric(as.character((lcm_table$Log2_fc)))
lcm_table$Gene <- as.character(lcm_table$Gene)
lcm_table$Segment <- as.character(lcm_table$Segment)
lcm_table2 <- pivot_wider(lcm_table, names_from=Segment, values_from=Log2_fc)
lcm_table2 <- as.data.frame(lcm_table2)
rownames(lcm_table2) <- lcm_table2$Gene
lcm_table2 <- select(lcm_table2, -matches("Gene"))

#Restrict to genes that are variable in sc/sn data and also found in LCM dataset
var_shared <- intersect(tis.integrated[["RNA"]]@var.features, lcm$Gene)

#Find Pearson correlation between each cell in integrated dataset and each LCM segment
d <- cor(tis.integrated[["RNA"]]@scale.data[var_shared,],lcm_table2[var_shared,])

#Assign each cell to segment with maximum Pearson correlation
assignments <- colnames(d)[apply(d,1,which.max)]
assignments <- cbind(rownames(d), assignments)
assignments <- cbind(assignments, as.character(Idents(tis.integrated)))
colnames(assignments) = c("cell_id", "assignments", "celltype")

#Count number of cells in each cluster that are assigned to each LCM segment (both raw counts and scaled across LCM segments for each cell cluster)
b <- table(assignments[,"assignments"], assignments[,"celltype"])
b <- cbind(rownames(b), b)
colnames(b) <- c("lcm_segment", colnames(b)[2:length(colnames(b))])
tmp_table_long <- pivot_longer(as.data.frame(b), -lcm_segment, names_to = "clust", values_to = "count")
scaled_table <- scale(table(assignments[,"assignments"], assignments[,"celltype"]))
tmp_rownames <- rownames(scaled_table)
scaled_table <- cbind(tmp_rownames, scaled_table)
scaled_table_long <- pivot_longer(as.data.frame(scaled_table), -tmp_rownames, names_to = "clust", values_to = "count_scaled")
tmp_table_long <- cbind(tmp_table_long, scaled_table_long$count_scaled)
names(tmp_table_long) <- c("lcm_segment", "clust", "count", "count_scaled")
tmp_table_long$count_scaled <- as.numeric(as.character(tmp_table_long$count_scaled))
tmp_table_long$count <- as.numeric(as.character(tmp_table_long$count))

#Reorder segments and sc/sn clusters
tmp_table_long$lcm_segment <- factor(tmp_table_long$lcm_segment, levels=c("Glom", "ProxTub", "TAL",  "DCT", "CD", "INT"))
tmp_table_long$int_clust <- factor(tmp_table_long$clust, levels=c("EC-PTC", "EC-AEA", "EC-GC",   "MC", "POD", "PT/PEC", "PT", "DTL", "ATL/TAL", "TAL", "DCT", "CNT", "PC", "tPC-IC", "IC", "FIB", "VSMC/P","T", "B", "MYL", "NKT/NKC", "MAC"))

#Plot
col2 <- c("#67001F", "#B2182B", "#D6604D", "#F4A582",
     		     		"#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
					   	      		 "#4393C3", "#2166AC", "#053061")
library(ggplot2)
p <- ggplot(tmp_table_long, aes(int_clust, lcm_segment)) + geom_tile(aes(fill=count_scaled)) + geom_text(aes(label=count)) + theme(text = element_text(size=36), 
axis.text.x = element_text(angle = 35, hjust = 1)) + scale_fill_gradientn(colors=rev(col2))

ggsave("integrated_mappings_flipped.pdf", p, height = 6, width = 24)
