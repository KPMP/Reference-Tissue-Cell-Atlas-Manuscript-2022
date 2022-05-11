#Copyright 2022.The code was written by Jens Hansen working for the Ravi Iyengar Lab
#and the Kidney Precision Medicine Project (KPMP)
#It is made availabe under the Apache License, Version 2.0 (the "License");
#you may not use this file except in compliance with the License.
#You may obtain a copy of the License at
#
#http://www.apache.org/licenses/LICENSE-2.0
#
#Unless required by applicable law or agreed to in writing, software
#distributed under the License is distributed on an "AS IS" BASIS,
#WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#See the License for the specific language governing permissions and
#limitations under the License.

rm(list = ls());
library(Rtsne);
library(ClassDiscovery);
library(colorspace);
library(dendextend);
library(circlize);
library(ape)
library(gplots)
library(ggplot2)
library(gridExtra)
library(grid)

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

###################### Function Definitions ######################


singleNucleus_color = "darkred";
singleCellMich_color = "red3";
singleCellUCSF_color= "firebrick1"
lcm_rnaSeq_color = "plum"
lcm_prot_color = "lightskyblue"
nsc_prot_color = "blue"

base_directory = paste(getwd(),"/",sep='');#"D:/KPMP_reference_atlas_code/"
base_directory = "D:/KPMP_reference_atlas_code/"
directory = paste(base_directory,"/Results/GlomVsProxTub/",sep='');
fileName= "Standardized_dataset_collapsed_on_method_sharedGenesProteins.txt";
fileName_singleMethods = "StandardizedData_log2Ratios.txt";
complete_fileName = paste(directory,fileName,sep='')
complete_fileName_singleMethods = paste(directory,fileName_singleMethods,sep='')


pdf_complete_fileName = paste(directory,"Standardized_dataset_collapsed_on_method_sharedGenesProteins",".pdf",sep='');
pdf(pdf_complete_fileName, width=8.5, height=11);

Data = read.table(file=complete_fileName,stringsAsFactors = FALSE, sep="\t", quote="",header=TRUE)

integration_terms = unique(Data$KPMP_data_integration_term)[1]
max_median = max(abs(Data$ScoreOfInterest_median))
max_average = max(abs(Data$ScoreOfInterest_average))



#Plot rnaSeq vs proteomics
method1 = "Proteomics"
method2 = "RNASeq"

method1_label = paste("Combined\n",method1)
method2_label = paste("Combined\n",method2)
color1 = "blue"
color2 = "red"

cool = rainbow(41, start=rgb2hsv(col2rgb('cyan'))[1], end=rgb2hsv(col2rgb('blue'))[1])
warm = rainbow(41, start=rgb2hsv(col2rgb('darkred'))[1], end=rgb2hsv(col2rgb('orange1'))[1])

warm = gray.colors(41,start=0.1, end=0.6)

cols = c(rev(warm))
Color_map <- colorRampPalette(cols)(41)


integration_term = integration_terms[1];
indexCurrent = which(Data$KPMP_data_integration_term==integration_term)
current_data = Data[indexCurrent,]
current_data = current_data[order(current_data$Gene_symbol),]
indexMethod1 = which(current_data$Method==method1)
indexMethod2 = which(current_data$Method==method2)
current_data_method1 = current_data[indexMethod1,]
current_data_method2 = current_data[indexMethod2,]
indexUnequal = which(current_data_method1$Gene_symbol!=current_data_method2$Gene_symbol)
if (length(indexUnequal)!=0) { break; }
indexUnequal = which(current_data_method1$Gene_symbol!=current_data_method2$Gene_symbol)
if (length(indexUnequal)!=0) { break; }
plot_data = current_data_method1;

plot_data$ScoreOfInterest_median_method1 = current_data_method1$ScoreOfInterest_median;
plot_data$ScoreOfInterest_average_method1 = current_data_method1$ScoreOfInterest_average;
plot_data$ScoreOfInterest_populationSD_method1 = current_data_method1$ScoreOfInterest_populationSD;

plot_data$ScoreOfInterest_median_method2 = current_data_method2$ScoreOfInterest_median;
plot_data$ScoreOfInterest_average_method2 = current_data_method2$ScoreOfInterest_average;
plot_data$ScoreOfInterest_populationSD_method2 = current_data_method2$ScoreOfInterest_populationSD;

plot_data$Plot_text = ""
indexAboveThreshold = which( (plot_data$ScoreOfInterest_average_method1>3)
                            &(plot_data$ScoreOfInterest_average_method2>1))
indexBelowThreshold = which( (plot_data$ScoreOfInterest_average_method1<(0))
                            &(plot_data$ScoreOfInterest_average_method2<(0)))
indexLabel = c(indexAboveThreshold,indexBelowThreshold)
label_genes = c("PODXL","CLIC5","NES","NPHS2","AIF1","PTPRO","NES","CUBN",
                "VIL1","UGT2B7","GPX3","MIOX")  #"PDZK1","LRP2","FTL",
indexLabel = which(plot_data$Gene_symbol %in% label_genes)


plot_data$Plot_text[indexLabel] = plot_data$Gene_symbol[indexLabel]
plot_data$Point_color = "black"
plot_data$Point_color[indexLabel] = "white"

correlation = round(100*cor(plot_data$ScoreOfInterest_average_method1,plot_data$ScoreOfInterest_average_method2,method="pearson"))/100
title_color = Color_map[round((correlation-0.35)*100)]
plot_title = paste("cor: ",correlation,"",sep='')
graph = ggplot(plot_data,aes(x=ScoreOfInterest_average_method1,y=ScoreOfInterest_average_method2,label=Plot_text))
graph = graph + geom_point(size=0.5,shape=21,aes(color=Point_color,fill=Point_color))
graph = graph + theme(panel.background = element_rect(fill="white",color="lightgray"))
graph = graph + scale_fill_manual(values=unique(plot_data$Point_color))
graph = graph + scale_color_manual(values=unique(plot_data$Point_color))
graph = graph + geom_hline(yintercept=0,col="lightgray")
graph = graph + geom_vline(xintercept=0,col="lightgray")
graph = graph + geom_text(size=3.0)
graph = graph + theme(axis.title=element_text(size=30,hjust=0.5,face="bold"))
graph = graph + theme(axis.text = element_text(size=30))
graph = graph + theme(axis.title.x = element_text(size=30, color = color1))
graph = graph + theme(axis.title.y = element_text(size=30, color = color2))
graph = graph + theme(plot.title = element_text(size=30,hjust=0.5,face="bold",color=title_color))
graph = graph + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_blank(), axis.line = element_line(colour = "black"))
graph = graph + ylim(-max_average,max_average) + xlim(-max_average,max_average)
graph = graph + ggtitle(plot_title) + xlab(method1_label) + ylab(method2_label)
graph = graph + theme(aspect.ratio=1)
average_graph = graph;

grid.arrange(average_graph,ncol=1,nrow=1)


#Plot individual methods
indexPremiere = which(Data$Method=="Single cell Premiere")
Data$Method[indexPremiere] = "SC RNASeq Premiere"

unique_methods = unique(Data$Method)
unique_methods = unique_methods[!unique_methods %in% c("RNASeq","Proteomics")]

indexProt = c(grep("proteomics",unique_methods),grep("Proteomics",unique_methods))
indexSingleCell_premiere = which(unique_methods == "SC RNASeq")
indexSingleCell_ucsf = which(unique_methods == "SC RNASeq II")
indexSingleNucleus = which(unique_methods == "SN RNASeq")
indexLCMRNASeq = which(unique_methods == "LMD RNASeq")
indexLCMProteomics = which(unique_methods == "LMD Proteomics")
indexNSCProteomics = which(unique_methods == "NSC Proteomics")

indexRightOrder = c(indexSingleNucleus,indexSingleCell_premiere,indexSingleCell_ucsf,indexLCMRNASeq,indexLCMProteomics,indexNSCProteomics)

unique_methods = unique_methods[indexRightOrder]
unique_methods_labels = c("SN RNASeq","SC RNASeq","SC RNASeq II","LMD RNASeq","LMD Proteomics","NSC proteomics")
unique_colors = c(singleNucleus_color,singleCellMich_color,singleCellUCSF_color,lcm_rnaSeq_color,lcm_prot_color,nsc_prot_color)
length_unique_methods = length(unique_methods)

compare_methods1 = c()
compare_methods2 = c()
for (indexM1 in 1:(length(unique_methods)-1))
{#Begin
   method2 = unique_methods[indexM1]
   for (indexM2 in (indexM1+1):(length(unique_methods)))
   {#Begin
      method1 = unique_methods[indexM2]
      {#Begin
         compare_methods1 = c(compare_methods1,method1)
         compare_methods2 = c(compare_methods2,method2)
      }#End
   }#End
}#End

methods1 = c(compare_methods1)
methods2 = c(compare_methods2)

indexRow = 0;
indexCol = 0;

layout_matrix = array(NA,c(length_unique_methods,length_unique_methods),dimnames = list(unique_methods,unique_methods))

average_plots = vector('list',length(methods1))

for (indexMethod in 1:length(methods1))
{#Begin
   method1 = methods1[indexMethod]
   method2 = methods2[indexMethod]

   layout_index1 = which(unique_methods==method1)
   layout_index2 = which(unique_methods==method2)
   method1_label = unique_methods_labels[layout_index1]
   method2_label = unique_methods_labels[layout_index2]
   color1 = unique_colors[layout_index1]
   color2 = unique_colors[layout_index2]
   if (length(layout_index1)!=1) { break }
   if (length(layout_index2)!=1) { break }
   layout_matrix[layout_index1,layout_index2] = indexMethod;
   layout_matrix[layout_index2,layout_index1] = indexMethod;

   integration_term = integration_terms[1];
   indexCurrent = which(Data$KPMP_data_integration_term==integration_term)
   current_data = Data[indexCurrent,]
   current_data = current_data[order(current_data$Gene_symbol),]
   indexMethod1 = which(current_data$Method==method1)
   indexMethod2 = which(current_data$Method==method2)
   current_data_method1 = current_data[indexMethod1,]
   current_data_method2 = current_data[indexMethod2,]
   indexUnequal = which(current_data_method1$Gene_symbol!=current_data_method2$Gene_symbol)
   if (length(indexUnequal)!=0) { break; }
   plot_data = current_data_method1;

   plot_data$ScoreOfInterest_median_method1 = current_data_method1$ScoreOfInterest_median;
   plot_data$ScoreOfInterest_average_method1 = current_data_method1$ScoreOfInterest_average;
   plot_data$ScoreOfInterest_populationSD_method1 = current_data_method1$ScoreOfInterest_populationSD;

   plot_data$ScoreOfInterest_median_method2 = current_data_method2$ScoreOfInterest_median;
   plot_data$ScoreOfInterest_average_method2 = current_data_method2$ScoreOfInterest_average;
   plot_data$ScoreOfInterest_populationSD_method2 = current_data_method2$ScoreOfInterest_populationSD;

   correlation = round(100*cor(plot_data$ScoreOfInterest_average_method1,plot_data$ScoreOfInterest_average_method2,method="pearson"))/100
   title_color = Color_map[round((correlation-0.35)*100)]

   plot_title = paste("cor: ",correlation,"",sep='')
   graph = ggplot(plot_data,aes(x=ScoreOfInterest_average_method1,y=ScoreOfInterest_average_method2,label=Gene_symbol))
   graph = graph + geom_hline(yintercept=0,col="lightgray")
   graph = graph + geom_vline(xintercept=0,col="lightgray")
   #graph = graph + geom_text(size=2.0)
   graph = graph + geom_point(size=0.5)
   graph = graph + theme(axis.title.x = element_text(color = color1,size=10))
   graph = graph + theme(axis.title.y = element_text(color = color2,size=10))

   graph = graph + theme(plot.title = element_text(size=15,hjust=0.5,face="bold",color=title_color))
   graph = graph + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         panel.background = element_blank(), axis.line = element_line(colour = "black"))
   graph = graph + ylim(-max_average,max_average) + xlim(-max_average,max_average)
   graph = graph + ggtitle(plot_title) + xlab("") + ylab("")
   #graph = graph + ylab(method2_label)
   #graph = graph + xlab(method1_label)
   graph = graph + theme(aspect.ratio=1)
   average_graph = graph;

   average_plots[[indexMethod]] = average_graph;
}#End


index12 = layout_matrix[1,2]
index13 = layout_matrix[1,3]
index14 = layout_matrix[1,4]
index15 = layout_matrix[1,5]
index16 = layout_matrix[1,6]
index23 = layout_matrix[2,3]
index24 = layout_matrix[2,4]
index25 = layout_matrix[2,5]
index26 = layout_matrix[2,6]
index34 = layout_matrix[3,4]
index35 = layout_matrix[3,5]
index36 = layout_matrix[3,6]
index45 = layout_matrix[4,5]
index46 = layout_matrix[4,6]
index56 = layout_matrix[5,6]


blank <- grid.rect(gp=gpar(col="white"))


grid.arrange(average_plots[[index16]],average_plots[[index15]],average_plots[[index14]],average_plots[[index13]],average_plots[[index12]],
             average_plots[[index26]],average_plots[[index25]],average_plots[[index24]],average_plots[[index23]],blank,
             average_plots[[index36]],average_plots[[index35]],average_plots[[index34]],blank,                   blank,
             average_plots[[index46]],average_plots[[index45]],blank,                   blank,                   blank,
             average_plots[[index56]],blank,                   blank,                   blank,                   blank,
             ncol=length(layout_matrix[,1])-1,nrow=length(layout_matrix[1,])-1)

dev.off()


# Change back to start-up directory
if(!is.null(start_dir)) setwd(start_dir)
