unlink(".RData")
rm(list = ls());
library(SingleCellExperiment)
library(Seurat)
library(ggplot2)
library(sctransform)
library(mclust)
library(dplyr)
library(beeswarm)
library(clustree)
library(Matrix)
library(gridExtra)
library(reticulate)
library(harmony)
use_python("c://ProgramData//Anaconda3//")
library(umap)
library(stringr)
library(doParallel)
library(scales)

# Change to data directory
data_dir <- "~/Data"
start_dir <- NULL
if(dir.exists(data_dir))
{
  abs_data_dir <- normalizePath(data_dir)
  start_dir <- getwd()
  if(abs_data_dir!=start_dir) setwd(abs_data_dir)
  else start_dir <- NULL
}

####################### Specify directories - BEGIN

base_directory = paste(getwd(),"/",sep='');#"D:/KPMP_reference_atlas_code/"
base_directory = "D:/KPMP_reference_atlas_code/"

data_subdirectory_name = "WU_PMID29980650//"
baseFileName_groups = list("Control" = c(paste("MTS.kidney.dge",sep='')))

input_directory = paste(base_directory, "Experimental_data/SingleCell_datasets/", sep='')
input_fileName = "Wu_MTS.kidney_PMID29980650.dge.txt"
input_complete_fileName = paste(input_directory,input_fileName,sep='')

results_directory = paste(base_directory,"Results/WU_PMID29980650_singleCellAnalysis/",sep='')
dir.create(results_directory)
results_documentation_directory = paste(results_directory,"Documentation/",sep='')
dir.create(results_documentation_directory)

####################### Specify directories - END

{#Begin - Read essential cell type genes library
  directory = paste(base_directory,"Experimental_data//Sample_metadata_additional_datasets//",sep='')
  library_fileName = "Essential_genes.txt"
  complete_library_fileName = paste(directory,library_fileName,sep='')
  gene_library = read.csv(file=complete_library_fileName,header=TRUE,sep='\t')
}#End - Read essential cell type genes library

###Start single cell analysis
baseAnalysisFileName = "Wu_PMID33850129"

minimum_feature_count_per_cell = 200;
maximum_feature_count_per_cell = 5000
top_considered_features = 2000;
pc_dimensions = 1:30;
minimum_cells = 3;
Generate_figures_for_visualization = TRUE
mitochondrial_label = "^MT-";
max_percentage_mitochondrial_genes = 20;
logfc_thresholds = 0;
additional_vars_to_regress = "percent.mt"

options(future.globals.maxSize= 2000000*1024^2)

Raw_data <- read.table(file=input_complete_fileName,header=TRUE,stringsAsFactors = FALSE,sep='\t')

rowsums = rowSums(Raw_data)
indexNonZeroGeneCounts = which(rowsums!=0)
Raw_data = Raw_data[indexNonZeroGeneCounts,]

seurat_object = CreateSeuratObject(Raw_data,min.cells=minimum_cells,min.features=0)

###################################################################################################################################
#Quality control: Removal of cell with too few and too many nFeatures and too high mitochondrial feature counts

seurat_object <- subset(x = seurat_object, subset = nFeature_RNA >= 100)
seurat_object[["percent.mt"]] <- PercentageFeatureSet(object = seurat_object, pattern = mitochondrial_label)

if (Generate_figures_for_visualization)
{#Begin - Generate violin plots
  vln_plots = list()

  # Visualize QC metrics as a violin plot
  vlnplot_nFeature_RNA = VlnPlot(object = seurat_object, features = c("nFeature_RNA"), ncol = 1)
  vlnplot_nFeature_RNA = vlnplot_nFeature_RNA + geom_hline(yintercept=minimum_feature_count_per_cell,col="gray60",size=2)
  vlnplot_nFeature_RNA = vlnplot_nFeature_RNA + geom_hline(yintercept=maximum_feature_count_per_cell,col="gray60",size=2)
  vlnplot_nCount_RNA = VlnPlot(object = seurat_object, features = c("nCount_RNA"), ncol = 1)
  vlnplot_percent_mt = VlnPlot(object = seurat_object, features = c("percent.mt"), ncol = 1)
  vlnplot_percent_mt = vlnplot_percent_mt + geom_hline(yintercept=max_percentage_mitochondrial_genes,col="gray60",size=2)
  vln_plots[[length(vln_plots)+1]] = vlnplot_nFeature_RNA
  vln_plots[[length(vln_plots)+1]] = vlnplot_nCount_RNA
  vln_plots[[length(vln_plots)+1]] = vlnplot_percent_mt

  min_cols_count = 3;
  complete_png_fileName = paste(results_documentation_directory,"ViolinPlot_qualityControl.png",sep='')
  cols_count = min(min_cols_count,length(vln_plots))
  rows_count = ceiling(length(vln_plots)/cols_count)
  png(complete_png_fileName,width=800*cols_count,height=500*rows_count,res=75);
  do.call("grid.arrange",c(vln_plots,nrow=rows_count,ncol=cols_count))
  dev.off()
}#End - Generate violin plots

seurat_object <- subset(x = seurat_object, subset = nFeature_RNA <= maximum_feature_count_per_cell & nFeature_RNA >= minimum_feature_count_per_cell & percent.mt<=max_percentage_mitochondrial_genes)
seurat_object <- SCTransform(seurat_object, verbose=TRUE,variable.features.n = top_considered_features, vars.to.regress = c(additional_vars_to_regress))
seurat_object <- RunPCA(seurat_object, npcs = 30, verbose = FALSE)
seurat_object <- RunUMAP(seurat_object, dims = pc_dimensions, verbose = FALSE)
seurat_object <- RunTSNE(seurat_object, npcs = 30, verbose = FALSE)
seurat_object <- FindNeighbors(object = seurat_object, dims = pc_dimensions)

selected_resolution = 0.76
seurat_object = FindClusters(object = seurat_object, resolution = selected_resolution)

if (Generate_figures_for_visualization)
{#Begin - Visualize cell neighbors
  group = "seurat_clusters"
  neighborhood_plots=list()
  umap_neighborhood_plots=list()
  pt_size = 1;
  tsne_plot_cluster = DimPlot(object = seurat_object, dims=c(1,2), reduction = "tsne", pt.size=pt_size,label=TRUE,group.by = group)
  tsne_plot_cluster = tsne_plot_cluster + ggtitle(paste("TSNE - ",group," - Res: ",selected_resolution,sep=''))
  neighborhood_plots[[length(neighborhood_plots)+1]] = tsne_plot_cluster
  pca_plot_cluster = DimPlot(object = seurat_object, dims=c(1,2), reduction = "pca", pt.size=pt_size,label=TRUE,group.by = group)
  pca_plot_cluster = pca_plot_cluster + ggtitle(paste("PCA - ",group," - Res: ",selected_resolution,sep=''))
  neighborhood_plots[[length(neighborhood_plots)+1]] = pca_plot_cluster
  umap_plot_cluster = DimPlot(object = seurat_object, dims=c(1,2), reduction = "umap", pt.size=pt_size,label=TRUE,group.by = group)
  umap_plot_cluster = umap_plot_cluster + ggtitle(paste("UMAP - ",group," - Res: ",selected_resolution,sep=''))
  neighborhood_plots[[length(neighborhood_plots)+1]] = umap_plot_cluster
  umap_neighborhood_plots[[length(umap_neighborhood_plots)+1]] = umap_plot_cluster

  complete_png_fileName = paste(results_documentation_directory,"CellMaps.png",sep='')
  cols_count = min(min_cols_count,length(neighborhood_plots))
  rows_count = ceiling(length(neighborhood_plots)/cols_count)
  png(complete_png_fileName,width=600*cols_count,height=500*rows_count,res=75);
  do.call("grid.arrange",c(neighborhood_plots,nrow=rows_count,ncol=cols_count))
  dev.off()

  complete_png_fileName = paste(results_documentation_directory,"UMAP_cellMaps.png",sep='')
  cols_count = 1
  rows_count = ceiling(length(umap_neighborhood_plots)/cols_count)
  png(complete_png_fileName,width=600*cols_count,height=500*rows_count,res=75);
  do.call("grid.arrange",c(umap_neighborhood_plots,nrow=rows_count,ncol=cols_count))
  dev.off()
}#End - Visualize cell neighbors

{#Begin - Identify and write barcodes of each cluster
  col_names = c("Barcode","Cluster","Resolution")
  col_length = length(col_names);
  row_names = 1:length(seurat_object$seurat_clusters)
  row_length = length(row_names)
  barcode_cluster_associations = as.data.frame(array(NA,c(row_length,col_length),dimnames=list(row_names,col_names)))
  barcode_cluster_associations$Barcode = names(seurat_object$seurat_clusters)
  barcode_cluster_associations$Cluster = seurat_object$seurat_clusters
  barcode_cluster_associations$Resolution = selected_resolution;

  clusterBarcodes_fileName = "Barcodes_setaCC13_generated_from_WU_PMID29980650.txt"
  complete_clusterBarcodes_fileName = paste(results_directory,clusterBarcodes_fileName,sep='')
  write.table(barcode_cluster_associations,file=complete_clusterBarcodes_fileName,quote=FALSE,row.names=FALSE,sep='\t')
}#End - Identify and write barcodes of each cluster

{#Begin - Identify cluster marker genes and subject them to enrichment analysis of essential cell type genes
  cluster.markers = FindAllMarkers(object = seurat_object,
                                   logfc.threshold = 0,
                                   only.pos = TRUE,
                                   assay="SCT",
                                   slot="data")
  all_clusters = unique(cluster.markers$cluster)
  cluster.markers$Rank = -1;
  indexAll=1
  for (indexAll in 1:length(all_clusters))
  {#Begin
    current_cluster = all_clusters[indexAll]
    indexCurrentCluster = which(cluster.markers$cluster==current_cluster)
    current_markers = cluster.markers[indexCurrentCluster,]
    current_markers = current_markers[order(abs(current_markers$avg_log2FC),decreasing = TRUE),]
    current_markers = current_markers[order(abs(current_markers$p_val),decreasing = TRUE),]
    current_markers$Rank = rank(current_markers$p_val)
    cluster.markers[indexCurrentCluster,] = current_markers
  }#End
  indexMinusOne =  which(cluster.markers$Rank==-1)
  if (length(indexMinusOne)>0) { rm(cluster.markers)}

  background_genes = rownames(Raw_data)
  indexSigAdjPvalue = which(cluster.markers$p_val_adj<=0.05)
  indexBelowMaxRank = which(cluster.markers$Rank<=300)
  indexKeep = indexSigAdjPvalue[indexSigAdjPvalue %in% indexBelowMaxRank]
  cluster.markers = cluster.markers[indexKeep,]

  Col_names = c("Cluster","Cell_type","Pvalue","MinusLog10_pvalue","Rank","Overlap_count")
  Col_length = length(Col_names)
  Row_names = 1
  Row_length = length(Row_names)
  enrichment_results_base_line = array(NA,c(Row_length,Col_length),dimnames = list(Row_names,Col_names))
  enrichment_results_base_line = as.data.frame(enrichment_results_base_line)
  enrichment_results = c()

  clusters = unique(cluster.markers$cluster)
  indexKeep = which(gene_library$Essential_gene %in% background_genes)
  gene_library = gene_library[indexKeep,]
  cellTypes = unique(gene_library$Cell_type)
  length_background_genes = length(background_genes)
  for (indexCluster in 1:length(clusters))
  {#Beginn
    current_cluster = clusters[indexCluster]
    indexCurrentCluster = which(cluster.markers$cluster==current_cluster)
    current_markers = cluster.markers[indexCurrentCluster,]
    current_marker_genes = current_markers$gene
    for (indexCellType in 1:length(cellTypes))
    {#Begin
      cellType = cellTypes[indexCellType]
      indexCurrentCellType = which(gene_library$Cell_type==cellType)
      current_cellType_library = gene_library[indexCurrentCellType,]
      current_cellType_genes = current_cellType_library$Essential_gene

      overlap_length = length(current_cellType_genes[current_cellType_genes %in% current_marker_genes])
      a = overlap_length
      b = length(current_cellType_genes) - a
      c = length(current_marker_genes) - a
      d = length_background_genes - a - b - c

      contingency_table = array(NA,c(2,2))
      contingency_table[1,1] = a
      contingency_table[2,1] = b
      contingency_table[1,2] = c
      contingency_table[2,2] = d

      fisher_results = fisher.test(contingency_table,alternative = "greater")
      pvalue = fisher_results$p.value

      new_enrichments_results = enrichment_results_base_line
      new_enrichments_results$Cluster = current_cluster
      new_enrichments_results$Cell_type = cellType
      new_enrichments_results$Pvalue = pvalue
      new_enrichments_results$MinusLog10_pvalue = -log10(pvalue)
      new_enrichments_results$Overlap_count = a
      if (length(enrichment_results)==0) { enrichment_results = new_enrichments_results }
      else { enrichment_results = rbind(enrichment_results,new_enrichments_results) }
    }#End
  }#End

  for (indexCluster in 1:length(clusters))
  {#Beginn
    current_cluster = clusters[indexCluster]
    indexCurrentCluster = which(enrichment_results$Cluster==current_cluster)
    currentCluster_enrichment = enrichment_results[indexCurrentCluster,]
    currentCluster_enrichment$Rank = rank(currentCluster_enrichment$Pvalue)
    enrichment_results[indexCurrentCluster,] = currentCluster_enrichment
  }#End

  indexAtLeast3OverlapCount = which(enrichment_results$Overlap_count>2)
  enrichment_results = enrichment_results[indexAtLeast3OverlapCount,]
  enrichment_results = enrichment_results[order(enrichment_results$Pvalue),]
  enrichment_results = enrichment_results[order(enrichment_results$Cluster),]


  complete_enrichment_results_fileName = paste(results_directory,"Cell_type_enrichment_results.txt")
  write.table(enrichment_results,file=complete_enrichment_results_fileName,col.names=TRUE,row.names=FALSE,quote=FALSE,sep='\t')
}#End - Identify cluster marker genes and subject them to enrichment analysis of essential cell type genes


# Change back to start-up directory
if(!is.null(start_dir)) setwd(start_dir)
