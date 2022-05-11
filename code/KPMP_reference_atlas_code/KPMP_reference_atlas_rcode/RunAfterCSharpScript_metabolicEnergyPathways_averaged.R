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

unlink(".RData")
rm(list = ls(all.names=TRUE));
gc()
library(ggplot2)
library(gridExtra)
library(stringr)
library(ggridges)

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
base_directory = "D:/KPMP_reference_atlas_code/"

directory = paste(base_directory,"Results/Metabolic_activities/",sep='')
enrichment_fileName = "Energy_scps_rulesBasedAveraged_maxP_maxAdjP0.05_top500DEGsDEPs.txt";
complete_fileName = paste(directory,enrichment_fileName,sep='')

{#Begin - Define colors copy pasted from other script
  scp_color_list = list(   "pO2" = "firebrick"
                           ,"Renal elimination of drugs and toxins" = "cornflowerblue"
                           ,"Drug and toxin export via membrane transport proteins" = "orange"
                           ,"Citric acid cycle" = "red3"
                           ,"Peroxisomal beta oxidation" = "orange4"
                           ,"Peroxisomal beta oxidation specific enzymes" = "orange4"
                           ,"Mitochondrial and peroxisomal beta oxidation shared enzymes" = "navajowhite3"
                           ,"Mitochondrial beta oxidation specific enzymes" = "orange3"
                           ,"Mitochondrial beta oxidation" = "orange3"
                           ,"Lactate dehydrogenase" = "dodgerblue3"
                           ,"HMG-CoA pathway of ketone body formation specific enzymes" = "navajowhite4"
                           ,"Ketone body metabolism" = "tomato3"
                           ,"Ketone body catabolism" = "orange2"
                           ,"Ketone body catabolism specific enzymes" = "orange2"
                           ,"Ketone body formation and catabolism shared enzymes" = "bisque3"
                           ,"Pyruvate dehydrogenase" = "orangered2"
                           ,"Aerobic glycolysis" = "orangered2"
                           ,"Anaerobic glycolysis" = "dodgerblue3"
                           ,"Glycolysis" = "wheat2"
                           ,"Glycolysis specific enzymes" = "orange"
                           ,"Glycolysis and gluconeogenesis shared enzymes" = "bisque3"
                           ,"Glycolysis and gluconeogenesis" = "gray40"
                           ,"Gluconeogenesis specific enzymes" = "navajowhite4"
                           ,"Fructose and mannose metabolism" = "cornflowerblue"
                           ,"Gluatamate and glutamine metabolism" = "organge"
                           ,"Oxidative phosphorylation" = "darkorchid"
                           ,"MCT2" = "blue"
                           ,"MCT4" = "firebrick1"
                           ,"NADH-coenzyme Q oxidoreductase" = "orange1"
                           ,"Succinate-Q oxidoreductase" = "orange3"
                           ,"Cytochrome c oxidase" = "orangered1"
                           ,"Electron-transferring-flavoprotein dehydrogenase" = "orangered3"
                           ,"ATP synthase" = "red1"
                           ,"Coenzyme Q10" = "red3")

}#End - Define colors copy pasted from other script

centers = c("NSC Proteomics","LMD Proteomics","LMD RNASeq","SC RNASeq I","SN RNASeq")
centers = c("TCA","Aerobic","Anaerobic")#,"pO2","S3 pO2","S2 pO2","S1 pO2")

enrichment = read.csv(file=complete_fileName,stringsAsFactors = FALSE, sep='\t')

enrichment$Center = "error"
indexAerobic = which(enrichment$Scp %in% c("Aerobic glycolysis","Mitochondrial beta oxidation","Peroxisomal beta oxidation","Ketone body catabolism"))
enrichment$Center[indexAerobic] = "Aerobic"
indexCitric = which(enrichment$Scp=="Citric acid cycle")
enrichment$Center[indexCitric] = "TCA"
indexCitric = which(enrichment$Scp %in% c("Anaerobic glycolysis","Glycolysis","Glycolysis and gluconeogenesis","Ketone body metabolism"))
enrichment$Center[indexCitric] = "Anaerobic"
indexKeep = which(enrichment$Center %in% centers)
enrichment = enrichment[indexKeep,]
enrichment$KPMP_data_integration_term = gsub("_allPatients","",enrichment$KPMP_data_integration_term)
integration_terms = unique(enrichment$KPMP_data_integration_term)
integration_terms = c("Proximal_tubule","Descending limb","Descending_limb_medulla","Thin_ascending_limb","Thin_ascending_limb_medulla","Thick_ascending_limb","Thick_ascending_limb_medulla","Distal Convoluted Tubule","Connecting tubule","Principal cells","Principal_cell_medulla","Intercalated cells","Intercalated_cell_medulla")
Stacked_bardiagrams = list()

enrichment$KPMP_data_integration_term = gsub("-allPatients","",enrichment$KPMP_data_integration_term)

indexI=1
stacked_legend_scps_count=0
stacked_legend = ggplot() + theme_void()
for (indexI in 1:length(integration_terms))
{#Begin
  integration_term = integration_terms[indexI]
  indexCurrentIntegration_term = which(enrichment$KPMP_data_integration_term==integration_term)
  integration_enrichment = enrichment[indexCurrentIntegration_term,]
  plot_data = integration_enrichment
  indexMissing = which(!centers %in% plot_data$Center)
  if (length(indexMissing)>0)
  {#Begin
    for (indexIndex in 1:length(indexMissing))
    {#Begin
      add_plot_data = plot_data[1,]
      add_plot_data$Center = centers[indexMissing[indexIndex]]
      add_plot_data$Averaged_minusLog10Pvalue = 0;
      add_plot_data$Percentage_averaged_minusLog10Pvalue = 0;
      plot_data = rbind(plot_data,add_plot_data)
    }#End
  }#End

  plot_data$Scp = factor(plot_data$Scp,levels=names(scp_color_list))
  plot_data = plot_data[order(plot_data$Scp),]
  plot_data$Center = factor(plot_data$Center,levels=centers)
  current_plot_data_scps = unique(plot_data$Scp)
  stacked_colors = replicate(length(current_plot_data_scps),"gray")
  for (indexCG in 1:length(current_plot_data_scps))
  {#Begin
    current_group_scp = current_plot_data_scps[indexCG]
    if (current_group_scp %in% names(scp_color_list))
    {#Begin
      stacked_colors[indexCG] = scp_color_list[[current_group_scp]]
    }#End
  }#End

  headline = paste(integration_term,"\n(",unique(plot_data$Averaged_entities_count)," assays)",sep='')

  #Stacked_bardiagram = ggplot(data = plot_data,aes(y=Center,x=Percentage_averaged_minusLog10Pvalue,group=Scp,fill=Scp))
  #Stacked_bardiagram = Stacked_bardiagram + xlab("-log10(p) [%]") + ylab("")
  Stacked_bardiagram = ggplot(data = plot_data,aes(y=Center,x=Averaged_minusLog10Pvalue,group=Scp,fill=Scp))
  Stacked_bardiagram = Stacked_bardiagram + xlab("-log10(p)") + ylab("")
  Stacked_bardiagram = Stacked_bardiagram + geom_bar(stat="identity",color="black")
  Stacked_bardiagram = Stacked_bardiagram + xlim(c(0,11))
  Stacked_bardiagram = Stacked_bardiagram + ggtitle(headline) + theme(title = element_text(size=20,face=2,vjust=0.5))
  Stacked_bardiagram = Stacked_bardiagram + theme(axis.text.x = element_blank(),axis.ticks=element_blank())
  Stacked_bardiagram = Stacked_bardiagram + scale_fill_manual(values=stacked_colors)
  Stacked_bardiagram = Stacked_bardiagram + theme(axis.text.x = element_text(angle=0,size=15,face=2,vjust=0.5,hjust=0.5))
  Stacked_bardiagram = Stacked_bardiagram + theme(axis.text.y = element_text(size=17,face=2))

  if (length(unique(plot_data$Scp))>stacked_legend_scps_count)
  {#Begin
    current_stacked_legend = cowplot::get_legend(Stacked_bardiagram)
    stacked_legend = current_stacked_legend
    stacked_legend_scps_count = length(unique(plot_data$Scp))
  }#End
  Stacked_bardiagram = Stacked_bardiagram + theme(legend.position = "none")
  Stacked_bardiagram = Stacked_bardiagram + guides(fill=FALSE)

  Stacked_bardiagrams[[length(Stacked_bardiagrams)+1]] = Stacked_bardiagram
}#End
Stacked_bardiagrams[[length(Stacked_bardiagrams)+1]] = stacked_legend

png_fileName = str_replace(enrichment_fileName,".txt",".png")
complete_png_fileName = paste(directory,png_fileName,sep='')
cols_count = 1
rows_count = ceiling(length(Stacked_bardiagrams)+1)/cols_count
height_per_figure = 60+30*6
png(complete_png_fileName,width=600*cols_count,height=height_per_figure*rows_count,res=100);
do.call("grid.arrange",c(Stacked_bardiagrams,nrow=rows_count,ncol=cols_count))
dev.off()


# Change back to start-up directory
if(!is.null(start_dir)) setwd(start_dir)
