unlink(".RData")
rm(list = ls(all.names=TRUE));
gc()
library(gridExtra)
library(grid)
library(ggplot2)


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

base_directory = paste(getwd(),"/",sep='');
base_directory = "D://KPMP_reference_atlas_code//"

data_directory = paste(base_directory,"Results/Marker_genes_and_proteins/",sep='')
data_fileName = "KPMP_combined_data_maxP_maxAdjP0.05_top2000MarkerGenesProteins.txt"
data_complete_fileName = paste(data_directory,data_fileName,sep='')

top_considered_genes = 300

data = read.csv(file=data_complete_fileName,header=TRUE,stringsAsFactors = FALSE, sep='\t')

#if (min(data$Value_1st)<1.3) { rm(data)  }

indexKeep = which(data$Fractional_rank_based_on_value_types<=top_considered_genes)
data = data[indexKeep,]

metabolism_onto_fileName = paste(base_directory,"MBCO_datasets/Sphingolipid_metabolism.txt",sep='')
#metabolism_onto_fileName = "D:/C-sharp-files/EnrichR_datasets/Self/Custom_metabolism_Homo_sapiens_2018July17.txt"
metabolism = read.csv(file=metabolism_onto_fileName,header=TRUE,stringsAsFactors = FALSE,sep='\t')

sl_pathways = c("Serine incorporation into ER membrane",
                "Ceramide synthesis",
                "Sphingomyelin synthesis from ceramide",
                "Sphingosin metabolism"
                )

indexKeep=which(metabolism$Scp %in% sl_pathways)
metabolism = metabolism[indexKeep,]

metabolism$Scp_factor = factor(metabolism$Scp,levels=sl_pathways)
metabolism = metabolism[order(metabolism$Target_gene_symbol),]
metabolism = metabolism[order(metabolism$Scp_factor),]

scp_genes = unique(metabolism$Target_gene_symbol)
genes_in_correct_order = rev(scp_genes)

indexKeep = which(data$Gene_symbol %in% scp_genes)
data = data[indexKeep,]
data$Cell_segment_dataset = paste(data$Cell_segment,"-",data$Dataset,sep='')

remove_cell_segment_datasets = c("Glom - Mesangial-LMD RNASeq","Glom - EC-LMD RNASeq","Glom - EPC-LMD RNASeq",
                                 "Glom - Mesangial-LMD Proteomics","Glom - EC-LMD Proteomics","Glom - EPC-LMD Proteomics",
                                 "TI - Interstitium-LMD Proteomics","TI - IC-LMD Proteomics","TI - MAC Proteomics","TI - PC Proteomics","TI - TAL Proteomics",
                                 "TI - MAC-LMD Proteomics","TI - PC-LMD Proteomics","TI - TAL-LMD Proteomics",
                                 "Glom - Mesangial-NSC Proteomics","Glom - EC-NSC Proteomics","Glom - EPC-NSC Proteomics",
                                 "CD - IC-LMD RNASeq")
indexRemove = which(data$Cell_segment_dataset %in% remove_cell_segment_datasets)
indexKeep = 1:length(data[,1])
indexKeep = indexKeep[!indexKeep %in% indexRemove]
data = data[indexKeep,]

rename_cell_segment_dataset = c("POD-SC RNASeq I" = "POD-SC RNASeq",
                                "Glom - POD-LMD RNASeq" = "Glom-LMD RNASeq",
                                "Glom - POD-LMD Proteomics" = "Glom-LMD Proteomics",
                                "Glom - POD-NSC Proteomics" = "Glom-NSC Proteomics",
                                "CD - PC-LMD RNASeq" = "CD-LMD RNASeq",
                                "TI - PT-LMD Proteomics"="TI-LMD Proteomics"
)

for (indexRename in 1:length(rename_cell_segment_dataset))
{#Begin
   old_name = names(rename_cell_segment_dataset)[indexRename]
   new_name = rename_cell_segment_dataset[indexRename]
   indexOld = which(data$Cell_segment_dataset==old_name)
   data$Cell_segment_dataset[indexOld] = new_name
}#End

Plot_list = list()

cell_segment_dataset_in_correct_order = c("Glom-LMD RNASeq",
                                          "Glom-LMD Proteomics",
                                          "Glom-NSC Proteomics",
                                          "POD-SN RNASeq","POD-SC RNASeq",
                                          "MC-SN RNASeq",
                                          "PEC/LOH-SC RNASeq",
                                          "EPC-SN RNASeq",
                                          "TI-LMD Proteomics",
                                          "PT-NSC Proteomics",
                                          "ProxTub-LMD RNASeq",
                                          "PT-1-SN RNASeq",
                                          "PT-2-SN RNASeq",
                                          "PT-4-SN RNASeq",
                                          "PT-5-SN RNASeq",
                                          "PT-1-SC RNASeq",
                                          "PT-2-SC RNASeq",
                                          "PT-3-SC RNASeq",
                                          "PT-4-SC RNASeq",
                                          "PT-5-SC RNASeq",
                                          "PT-6-SC RNASeq",
                                          "PT-7-SC RNASeq",
                                          "DTL-SN RNASeq",
                                          "DTL-SC RNASeq",
                                          "ATL-1-SN RNASeq","ATL-2-SN RNASeq","ATL-3-SN RNASeq",
                                          "TAL-LMD RNASeq",
                                          "DCT-LMD RNASeq",
                                          "CNT-SN RNASeq",
                                          "CD-LMD RNASeq",
                                          "PC-1-SN RNASeq",
                                          "PC-3-SN RNASeq",
                                          "PC-SC RNASeq",
                                          "IC-A1-SN RNASeq",
                                          "IC-A2-SN RNASeq",
                                          "IC-B-SN RNASeq",
                                          "IC-A-SC RNASeq",
                                          "IC-B-SC RNASeq",
                                          "INT-SN RNASeq",
                                          "MAC-SC RNASeq"
)

data$Gene_symbol_factor = factor(data$Gene_symbol,levels=genes_in_correct_order)
data$Cell_segment_dataset_factor = factor(data$Cell_segment_dataset,levels=cell_segment_dataset_in_correct_order)
data$Value_1st = as.numeric(data$Value_1st)
Plot = ggplot(data,aes(x=Cell_segment_dataset_factor, y=Gene_symbol_factor, fill=Fractional_rank_based_on_value_types))
Plot = Plot + geom_tile() + guides(fill=guide_legend(title='-log10(p)'))
Plot = Plot + scale_fill_gradient(low="cadetblue1" , high="dodgerblue3",limits=c(1,top_considered_genes))
#Plot = Plot + ggtitle(paste(gsub("_"," ",ontology),": ",cellType,sep=''))
Plot = Plot + geom_text(aes(label=Fractional_rank_based_on_value_types))
Plot = Plot + theme(plot.title = element_text(size=50,face=2))
Plot = Plot + guides(fill=guide_legend(title="Rank"))
Plot = Plot + theme(legend.title = element_text(size=20))
Plot = Plot + theme(axis.title.x = element_text(size=20))
Plot = Plot + theme(axis.text.x = element_text(size=20,angle=90,hjust=0,vjust=3))
Plot = Plot + scale_x_discrete(position = "top")
Plot = Plot + theme(axis.text.y = element_text(size=25))
#Plot = Plot + theme(legend.position="none")
Plot = Plot + xlab("")
Plot = Plot + ylab("")
Plot_list[[length(Plot_list)+1]] = Plot

length_plots = length(Plot_list)
max_plot_per_figure = 4;
figures_count = ceiling(length_plots / max_plot_per_figure)
for (indexF in 1:figures_count)
{#Begin
  indexStart = (indexF-1)* max_plot_per_figure + 1
  indexEnd = min(c((indexF) * max_plot_per_figure,length_plots))
  current_plots = list()
  for (indexCurrentPlot in indexStart:indexEnd)
  {#Begin
    current_plots[[length(current_plots)+1]] = Plot_list[[indexCurrentPlot]]
  }#End
  fileName = paste("Sphingomyeline_top",top_considered_genes,"genes_","figNo",indexF,".png",sep='')
  complete_png_fileName = paste(data_directory,fileName,sep='')
  cols_count = min(2,length(current_plots))
  rows_count = ceiling(length(current_plots)/cols_count)
  png(complete_png_fileName,width=640*cols_count,height=600*rows_count,res=75);
  do.call("grid.arrange",c(current_plots,nrow=rows_count,ncol=cols_count))
  dev.off()
}#End



