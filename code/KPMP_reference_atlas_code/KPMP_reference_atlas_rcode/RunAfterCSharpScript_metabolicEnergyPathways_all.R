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

base_directory = paste(getwd(),"/",sep='');#
base_directory = "D:/KPMP_reference_atlas_code/"
directory = paste(base_directory,"Results/Metabolic_activities/",sep='')

{#Begin - Read and add custom metabolism and transmembrane transport enrichment results
   enrichment = c()
   input_directory = paste(base_directory,"Results//Metabolic_activities//",sep='')
   fileName_metabolism = "Energy_scps_all_cellSubTypes_segments_datasets.txt"
   fileName_metabolism_png = "Energy_scps_all_cellSubTypes_segments_datasets.png"
   complete_fileName_metabolism = paste(input_directory,fileName_metabolism,sep='')
   enrichment_metabolism = read.csv(file=complete_fileName_metabolism,sep="\t",stringsAsFactors = FALSE,header=TRUE)
   enrichment = rbind(enrichment_metabolism)#,enrichment_tmtransport,enrichment_mbcoL3)
}#End - Read and add custom metabolism and transmembrane transport enrichment results

results_directory = input_directory


{#Begin - Adjust cell types in sample names: PC-CNT -> CNT-PC
  enrichment$Dataset.and.cell..sub.type.or.segment = gsub("PC-CNT","CNT-PC",enrichment$Dataset.and.cell..sub.type.or.segment)
  enrichment$Dataset.and.cell..sub.type.or.segment = gsub("ProxTub","PT",enrichment$Dataset.and.cell..sub.type.or.segment)
  enrichment$Dataset.and.cell..sub.type.or.segment = gsub("CD - PC","CD",enrichment$Dataset.and.cell..sub.type.or.segment)
  enrichment$Dataset.and.cell..sub.type.or.segment = gsub("Glom - POD","Glom",enrichment$Dataset.and.cell..sub.type.or.segment)
  enrichment$Dataset.and.cell..sub.type.or.segment = gsub("INT - FIB","Int",enrichment$Dataset.and.cell..sub.type.or.segment)
  enrichment$Dataset.and.cell..sub.type.or.segment = gsub("TI-PT","TI",enrichment$Dataset.and.cell..sub.type.or.segment)
  enrichment$Dataset.and.cell..sub.type.or.segment = gsub("SC RNASeq I","SC RNASeq",enrichment$Dataset.and.cell..sub.type.or.segment)
  indexRemove1 = grep("CD - IC",enrichment$Dataset.and.cell..sub.type.or.segment)
  indexRemove2 = grep("INT - MAC",enrichment$Dataset.and.cell..sub.type.or.segment)
  indexRemove3 = grep("Glom - EC",enrichment$Dataset.and.cell..sub.type.or.segment)
  indexRemove4 = grep("Glom - EPC",enrichment$Dataset.and.cell..sub.type.or.segment)
  indexRemove5 = grep("Glom - Mesangial",enrichment$Dataset.and.cell..sub.type.or.segment)
  indexRemove = unique(c(indexRemove1,indexRemove2,indexRemove3,indexRemove4,indexRemove5))
  indexKeep = 1:length(enrichment[,1])
  indexKeep = indexKeep[!indexKeep %in% indexRemove]
  enrichment = enrichment[indexKeep,]
}#End - Adjust cell types in sample names: PC-CNT -> CNT-PC

SCP_groups = list(  # "Ion transport" = c("Sodium reabsorption","Chloride reabsorption","Phosphate reabsortpion","Hydrogen secretion","Bicarbonate reabsorption")#"Chloride transporter","Magnesium transporter","Phosphate transporter")
                    #,"Ion channels" = c("Sodium channel","Potassium channel","Calcium channel","Chloride channel","Magnesium channel","Phosphate channel")
                    #,"Vesicle exocytosis" = c("Exosome secretion","Exomer-mediated vesicle traffic from TGN to plasma membrane","Microtubule-based secretory vesicle traffic","Vesicle fusion with plasma membrane","Vesicle tethering at plasma membrane","Clathrin-mediated vesicle traffic from TGN to endosomal lysosomal system","Retrograde vesicle traffic from endosome to TGN")
                    #,"Vesicle exocytosis II" = c("Vesicle exocytosis","Vesicle traffic between TGN and endosom")
                    #,"Small molecule transport" = c("Amino acid transporter","Glucose transporter")
                    # "NaKATPase" = c("Sodium potassium ATPase alpha subunit","Sodium potassium ATPase beta subunit","Regulation of sodium potassium ATPase")
                    # "Renal elimination of drugs and toxins" = c("Renal elimination of drugs and toxins","Drug and toxin export via membrane transport proteins")
                    #,"Overall energy I" = c("Mitochondrial beta oxidation","Peroxisomal beta oxidation","Citric acid cycle","Glycolysis","Pyruvate dehydrogenase","Lactate dehydrogenase","Keton body metabolism")
                     "Beta oxidation" = c("Mitochondrial beta oxidation specific enzymes","Mitochondrial and peroxisomal beta oxidation shared enzymes","Peroxisomal beta oxidation specific enzymes")
                    ,"Keton body metabolism" = c("HMG-CoA pathway of ketone body formation specific enzymes","Ketone body catabolism specific enzymes","Ketone body formation and catabolism shared enzymes")
                    #,"Beta oxidation II" = c("Mitochondrial beta oxidation","Peroxisomal beta oxidation")
                    ,"Glucose metabolism" = c("Gluconeogenesis specific enzymes","Glycolysis and gluconeogenesis shared enzymes","Glycolysis specific enzymes","Pyruvate dehydrogenase","Lactate dehydrogenase")#"Fructose-2-6-bisphosphate metabolism",
                    #,"Cellular detoxification" = c("Cellular detoxification")
                    ,"Citric acid cycle" = c("Citric acid cycle")
                    ,"Oxidative phosphorylation" = c("NADH-coenzyme Q oxidoreductase","Succinate-Q oxidoreductase","Cytochrome c oxidase","Electron-transferring-flavoprotein dehydrogenase","ATP synthase","Coenzyme Q10")
                    ,"Metabolic summary II" = c("Mitochondrial beta oxidation","Peroxisomal beta oxidation","Aerobic glycolysis","Anaerobic glycolysis","Ketone body catabolism")
                    ,"Metabolic summary" = c("Mitochondrial beta oxidation","Peroxisomal beta oxidation","Glycolysis","Pyruvate dehydrogenase","Lactate dehydrogenase","Ketone body catabolism")
                    #,"Mitochondrial energy" = c("Oxidative phosphorylation","Citric acid cycle")
                    #,"Gluconeogenesis substrates" = c("Fructose and mannose metabolism","Glutamate and glutamine metabolism","MCT1","MCT2","MCT3","MCT4")
                    #,"Urea cycle" = c("Urea cycle")
                    #,"GPI" = c("GPI-anchor biosynthesis at ER","GPI-anchor attachment to protein","Attached GPI-anchor maturation")
                    #,"Phosphate" = c("Phosphate transmembrane transport","Parathyroid hormone receptor signaling","Phosphate reabsorption")
)

{#Begin - Remove selected cell subtypes
   remove_ucsd_cell_subtypes = c("PT-3","PC-2","Unk")
   indexUCSD = grep("SN RNASeq UCSDWU",enrichment$Dataset.and.cell..sub.type.or.segment)
   indexRemove_ucsd_subtypes = c();
   for (indexRemove in 1:length(remove_ucsd_cell_subtypes))
   {#Begin
      indexRemove_ucsd_subtypes = c(indexRemove_ucsd_subtypes,grep(remove_ucsd_cell_subtypes[indexRemove],enrichment$Dataset.and.cell..sub.type.or.segment))
   }#End
   indexRemove_ucsd_subtypes = indexRemove_ucsd_subtypes[indexRemove_ucsd_subtypes %in% indexUCSD]
   remove_samples = c("LMD Proteomics-TI - IC","LMD Proteomics-TI - TAL","LMD Proteomics-TI - MAC","LMD Proteomics-TI - PC")
   indexRemove_samples = c();
   for (indexRemove in 1:length(remove_samples))
   {#Begin
     indexRemove_samples = c(indexRemove_samples,grep(remove_samples[indexRemove],enrichment$Dataset.and.cell..sub.type.or.segment))
   }#End

   indexKeep = 1:length(enrichment[,1])
   indexKeep = indexKeep[!indexKeep %in% indexRemove_ucsd_subtypes]
   indexKeep = indexKeep[!indexKeep %in% indexRemove_samples]
   enrichment = enrichment[indexKeep,]
}#End - Remove selected cell subtypes

{#Begin - Define cell subtype colors
   podGlo_color = "aquamarine4";
   PTTI_color = "orange2"
   PT4_color = "orange2"
   PT5_color = "orange2"
   cnt_pc_color = "indianred2"
   cortex_PCCD_color = "mediumorchid2";
   medulla_PCCD_color = "slateblue4";
   cortex_henle_color = "dodgerblue1";
   medulla_henle_color = "navy";
   cellType_other_color = "black"
   interstitium_color = "gray40"
   endothelial_color = "cyan3"
   immune_color = "yellow"

   center_cell_subtype_colors_list = list()
   premiere_cell_subtype_colors = list()
   premiere_cell_subtype_colors[["POD"]] = podGlo_color
   premiere_cell_subtype_colors[["vSMC/MC"]] = podGlo_color
   premiere_cell_subtype_colors[["vSMC/P"]] = podGlo_color
   premiere_cell_subtype_colors[["GC-EC"]] = podGlo_color
   premiere_cell_subtype_colors[["PT-1"]] = PTTI_color
   premiere_cell_subtype_colors[["PT-2"]] = PTTI_color
   premiere_cell_subtype_colors[["PT-3"]] = PTTI_color
   premiere_cell_subtype_colors[["PT-4"]] = PTTI_color
   premiere_cell_subtype_colors[["PT-5"]] = PTTI_color
   premiere_cell_subtype_colors[["PT-6"]] = PTTI_color
   premiere_cell_subtype_colors[["PT-7"]] = PTTI_color
   premiere_cell_subtype_colors[["DTL"]] = cortex_henle_color
   premiere_cell_subtype_colors[["ATL"]] = cortex_henle_color
   premiere_cell_subtype_colors[["TAL"]] = cortex_henle_color
   premiere_cell_subtype_colors[["DCT"]] = PTTI_color
   premiere_cell_subtype_colors[["CNT"]] = PTTI_color
   premiere_cell_subtype_colors[["CNT-PC"]] = cnt_pc_color
   premiere_cell_subtype_colors[["PC"]] = cortex_PCCD_color
   premiere_cell_subtype_colors[["tPC-IC"]] = cortex_PCCD_color
   premiere_cell_subtype_colors[["IC-A"]] = cortex_PCCD_color
   premiere_cell_subtype_colors[["IC-B"]] = cortex_PCCD_color
   premiere_cell_subtype_colors[["FIB"]] = interstitium_color
   premiere_cell_subtype_colors[["EC-AEA"]] = endothelial_color
   premiere_cell_subtype_colors[["EC-PT"]] = endothelial_color
   premiere_cell_subtype_colors[["MAC"]] = immune_color
   premiere_cell_subtype_colors[["MON"]] = immune_color
   premiere_cell_subtype_colors[["NKC"]] = immune_color
   premiere_cell_subtype_colors[["B cell"]] = immune_color
   premiere_cell_subtype_colors[["T Cell"]] = immune_color
   premiere_cell_subtype_colors[["T-Cyt"]] = immune_color
   premiere_cell_subtype_colors[["T-Cyt-MEM"]] = immune_color
   center_cell_subtype_colors_list[["SC RNASeq"]] = premiere_cell_subtype_colors

   ucsd_cell_subtype_colors = list()
   ucsd_cell_subtype_colors[["MC"]] = podGlo_color
   ucsd_cell_subtype_colors[["POD"]] = podGlo_color
   ucsd_cell_subtype_colors[["EPC"]] = podGlo_color
   ucsd_cell_subtype_colors[["vSMC/P"]] = podGlo_color
   ucsd_cell_subtype_colors[["PT-1"]] = PTTI_color
   ucsd_cell_subtype_colors[["PT-2"]] = PTTI_color
   ucsd_cell_subtype_colors[["PT-3"]] = PTTI_color
   ucsd_cell_subtype_colors[["PT-4"]] = PTTI_color
   ucsd_cell_subtype_colors[["PT-5"]] = PTTI_color
   ucsd_cell_subtype_colors[["Unk"]] = PTTI_color
   ucsd_cell_subtype_colors[["DTL"]] = medulla_henle_color
   ucsd_cell_subtype_colors[["ATL-1"]] = medulla_henle_color
   ucsd_cell_subtype_colors[["ATL-2"]] = medulla_henle_color
   ucsd_cell_subtype_colors[["ATL-3"]] = medulla_henle_color
   ucsd_cell_subtype_colors[["TAL-1"]] = medulla_henle_color
   ucsd_cell_subtype_colors[["TAL-2"]] = cortex_henle_color
   ucsd_cell_subtype_colors[["DCT"]] = PTTI_color
   ucsd_cell_subtype_colors[["CNT"]] = PTTI_color
   ucsd_cell_subtype_colors[["PC-1"]] = cortex_PCCD_color
   ucsd_cell_subtype_colors[["PC-2"]] = cortex_PCCD_color
   ucsd_cell_subtype_colors[["PC-3"]] = medulla_PCCD_color
   ucsd_cell_subtype_colors[["IC-A1"]] = cortex_PCCD_color
   ucsd_cell_subtype_colors[["IC-A2"]] = medulla_PCCD_color
   ucsd_cell_subtype_colors[["IC-B"]] = cortex_PCCD_color
   ucsd_cell_subtype_colors[["INT"]] = interstitium_color
   ucsd_cell_subtype_colors[["EC-1"]] = endothelial_color
   ucsd_cell_subtype_colors[["EC-2"]] = endothelial_color
   ucsd_cell_subtype_colors[["EC-3"]] = endothelial_color
   ucsd_cell_subtype_colors[["EC-4"]] = endothelial_color
   ucsd_cell_subtype_colors[["IMM"]] = immune_color
   center_cell_subtype_colors_list[["SN RNASeq"]] = ucsd_cell_subtype_colors

   liao_cell_subtype_colors = list()
   liao_cell_subtype_colors[["PT1"]] = PTTI_color
   liao_cell_subtype_colors[["PT2"]] = PTTI_color
   liao_cell_subtype_colors[["PT3"]] = PTTI_color
   liao_cell_subtype_colors[["PT4"]] = PTTI_color
   liao_cell_subtype_colors[["PT5"]] = PTTI_color
   liao_cell_subtype_colors[["PT6"]] = PTTI_color
   liao_cell_subtype_colors[["PT7"]] = PTTI_color
   liao_cell_subtype_colors[["PC"]] = cortex_PCCD_color
   liao_cell_subtype_colors[["LOH/DCT/IC"]] = cortex_henle_color
   center_cell_subtype_colors_list[["SC RNASeq Liao et al"]] = liao_cell_subtype_colors

   lmd_cell_subtype_colors = list()
   lmd_cell_subtype_colors[["Glom - POD"]] = podGlo_color
   lmd_cell_subtype_colors[["PT"]] = PTTI_color
   lmd_cell_subtype_colors[["TAL"]] = cortex_henle_color
   lmd_cell_subtype_colors[["DCT"]] = PTTI_color
   lmd_cell_subtype_colors[["CD"]] = cortex_PCCD_color
   center_cell_subtype_colors_list[["LMD RNASeq"]] = lmd_cell_subtype_colors
}#End - Define cell subtype colors

{#Begin - Keep only enrichment lines of selected centers
  selected_technologies = c( "SC RNASeq"
                            ,"SN RNASeq"
                            ,"LMD RNASeq"
                            ,"LMD Proteomics"
                            ,"NSC Proteomics"
                            )
  enrichment$Center = "error"
  indexKeep = c()
  for (indexST in 1:length(selected_technologies))
  {#Begin
    indexCurrentCenter=grep(selected_technologies[indexST],enrichment$Dataset.and.cell..sub.type.or.segment)
    indexKeep = c(indexKeep,indexCurrentCenter)
    enrichment$Center[indexCurrentCenter] = selected_technologies[indexST]
  }#End
  enrichment = enrichment[indexKeep,]
}#End - Keep only enrichment lines of selected centers

function_name = "Set cell types and subtypes to enrichment lines"
{#Begin - Set cell subtypes to enrichment lines
  enrichment$Cell_subtype = "error"
  enrichment$Cell_type = "error"
  sampleNames = unique(enrichment$Dataset.and.cell..sub.type.or.segment)
  for (indexS in 1:length(sampleNames))
  {#Begin
    sampleName = sampleNames[indexS]
    splitStrings = str_split(sampleName,"-")[[1]]

    if (length(splitStrings)==2)
    {#Begin
       cell_subtype = splitStrings[2]
    }#End
    if (length(splitStrings)>2)
    {#Begin
       cell_subtype = paste(splitStrings[2:length(splitStrings)],collapse="-")
    }#End
    if (!is.na(as.numeric(splitStrings[3]))) { cell_type = splitStrings[2] }
    if (is.na(as.numeric(splitStrings[3])))
    {#Begin
       if (splitStrings[2]=="IC") { cell_type = "IC"}
       if (splitStrings[2]=="EC") { cell_type = "EC"}
       if (splitStrings[2]=="TI") { cell_type = "TI"}
       if ((splitStrings[2]!="TI")&(splitStrings[2]!="EC")&(splitStrings[2]!="IC")) { cell_type = cell_subtype }
    }#End

    indexCurrentSampleName = which(enrichment$Dataset.and.cell..sub.type.or.segment==sampleName)
    enrichment$Cell_subtype[indexCurrentSampleName] = cell_subtype
    enrichment$Cell_type[indexCurrentSampleName] = cell_type
  }#End

  indexError = which(enrichment$Cell_subtype=="error")
  if (length(indexError)>0) { rm(enrichment); error_message = function_name }
  indexError = which(enrichment$Cell_type=="error")
  if (length(indexError)>0) { rm(enrichment); error_message = function_name }
}#End  Set cell subtypes to enrichment lines

{#Begin - Remove duplicated cells from LMD RNASeq and adjust LMD RNASeq cell subtype names
  indexRemove1 = which(enrichment$Cell_subtype=="Glom - EC")
  indexRemove2 = which(enrichment$Cell_subtype=="Glom - Mesangial")
  indexRemove3 = which(enrichment$Cell_subtype=="CD - IC")
  indexLMD = which(enrichment$Center=="LMD RNASeq OSUIU")
  indexRemove = c(indexRemove1,indexRemove2,indexRemove3)
  indexRemove = indexRemove[indexRemove %in% indexLMD]
  indexKeep = 1:length(enrichment[,1])
  indexKeep = indexKeep[!indexKeep %in% indexRemove]
  enrichment = enrichment[indexKeep,]

  unique(enrichment$Cell_subtype[indexLMD])
  indexLMD = which(enrichment$Center=="LMD RNASeq OSUIU")
  indexGlom = which(enrichment$Cell_subtype=="Glom - POD")
  indexGlom = indexGlom[indexGlom %in% indexLMD]
  indexCD = which(enrichment$Cell_subtype=="CD - PC")
  indexCD = indexCD[indexCD %in% indexLMD]
  enrichment$Cell_subtype[indexGlom] = "Glom"
  enrichment$Cell_subtype[indexCD] = "CD"
}#End - Remove duplicated cells from LMD RNASeq and adjust LMD RNASeq cell subtype names

{#Begin - Specify cortex and medulla cell types
  indexTAL1 = which(enrichment$Cell_subtype=="TAL-1")
  indexTAL2 = which(enrichment$Cell_subtype=="TAL-2")
  indexTAL = which(enrichment$Cell_subtype=="TAL")
  indexPC1 = which(enrichment$Cell_subtype=="PC-1")
  indexPC2 = which(enrichment$Cell_subtype=="PC-2")
  indexPC3 = which(enrichment$Cell_subtype=="PC-3")
  indexPC = which((enrichment$Cell_subtype=="PC"))
  indexCD = which((enrichment$Cell_subtype=="CD"))
  indexCNT_PC = which(enrichment$Cell_subtype=="CNT-PC")
  indexIC = which((enrichment$Cell_subtype=="IC"))
  indexICA1 = which(enrichment$Cell_subtype=="IC-A1")
  indexICA2 = which(enrichment$Cell_subtype=="IC-A2")
  indexICB_medulla = which((enrichment$Cell_subtype=="IC-B")&(enrichment$Center=="SN RNASeq"))
  indexICB_cortex = which((enrichment$Cell_subtype=="IC-B")&(enrichment$Center=="SC RNASeq"))
  indexICA = which(enrichment$Cell_subtype=="IC-A")
  indextPC_IC = which(enrichment$Cell_subtype=="tPC-IC")
  indexICB = which(enrichment$Cell_subtype=="IC-B")
  enrichment$Cell_type[indexTAL] = "TAL_medulla"
  enrichment$Cell_type[indexTAL1] = "TAL_medulla"
  enrichment$Cell_type[indexTAL2] = "TAL_cortex"
  enrichment$Cell_type[indexPC] = "PC cortex"
  enrichment$Cell_type[indexCNT_PC] = "PC cortex"
  enrichment$Cell_type[indextPC_IC] = "PC cortex"
  enrichment$Cell_type[indexCD] = "PC cortex"
  enrichment$Cell_type[indexPC1] = "PC cortex"
  enrichment$Cell_type[indexPC2] = "PC cortex"
  enrichment$Cell_type[indexPC3] = "PC medulla"
  enrichment$Cell_type[indexIC] = "IC cells cortex"
  enrichment$Cell_type[indexICA] = "IC cortex"
  enrichment$Cell_type[indexICA1] = "IC cortex"
  enrichment$Cell_type[indexICA2] = "IC medulla"
  enrichment$Cell_type[indexICB_medulla] = "IC medulla"
  enrichment$Cell_type[indexICB_cortex] = "IC medulla"
}#End - Specify cortex and medulla cell types

function_name = "Add cell types, segments and figure groups to enrichment lines"
{#Begin - Add cell types, segments and figure groups to enrichment lines
  orderedCellType_figureGroup_list = c( "MC" = "Glomerulus"
                                        ,"Glom" = "Glomerulus"
                                        ,"GC-EC" = "Glomerulus"
                                        ,"TI" = "Nephron"
                                        ,"vSMC/MC" = "Glomerulus"
                                       ,"GC" = "Endothelial"
                                       ,"EC" = "Endothelial"
                                       ,"PT" = "Nephron"
                                       ,"DTL" = "Nephron"
                                       ,"ATL" = "Nephron"
                                       ,"TAL" = "Nephron"
                                       ,"TAL_medulla" = "Nephron"
                                       ,"TAL_cortex" = "Nephron"
                                       ,"DCT" = "Nephron"
                                       ,"CNT" = "Nephron"
                                       ,"PC cortex" = "Nephron"
                                       ,"IC cortex" = "Nephron"
                                       ,"PC medulla" = "Nephron"
                                       ,"PC" = "Nephron"
                                       ,"IC medulla" = "Nephron"
                                       ,"IC" = "Nephron"
                                       ,"B Cell" = "Immune and interstitium"
                                       ,"MON" = "Immune and interstitium"
                                       ,"NKC" = "Immune and interstitium"
                                       ,"PEC/LOH" = "Epithelial"
                                       ,"POD" = "Glomerulus"
                                       ,"T Cell" = "Immune and interstitium"
                                       ,"T-CYT" = "Immune and interstitium"
                                       ,"T-CYT-MEM" = "Immune and interstitium"
                                       ,"EPC" = "Epithelial"
                                       ,"IMM" = "Immune and interstitium"
                                       ,"vSMC/P" = "Immune and interstitium"
                                       ,"INT" = "Immune and interstitium"
                                       ,"FIB" = "Immune and interstitium"
                                       ,"MAC" = "Immune and interstitium"
                                       ,"vSMC" = "Vascular"
  )

  orderedCellType_segment_list = c("MC" = "Glomerulus"
                                   ,"Glom" = "Glomerulus"
                                   ,"GC-EC" = "Glomerulus"
                                   ,"TI" = "Proximal_tubule"
                                   ,"vSMC/MC" = "Glomerulus"
                                  ,"GC" = "Endothelial"
                                  ,"EC" = "Endothelial"
                                  ,"PT" = "Proximal_tubule"
                                  ,"DTL" = "Loop of Henle"
                                  ,"ATL" = "Loop of Henle"
                                  ,"TAL" = "Loop of Henle"
                                  ,"DCT" = "Distal Convoluted Tubule"
                                  ,"TAL_medulla" = "Loop of Henle"
                                  ,"TAL_cortex" = "Loop of Henle"
                                  ,"CNT" = "Connecting tubule"
                                  ,"PC cortex" = "Collecting duct"
                                  ,"IC cortex" = "Collecting duct"
                                  ,"IC" = "Collecting duct"
                                  ,"PC medulla" = "Collecting duct"
                                  ,"PC" = "Collecting duct"
                                  ,"IC medulla" = "Collecting duct"
                                  ,"B Cell" = "Immune and interstitium"
                                  ,"MON" = "Immune and interstitium"
                                  ,"NKC" = "Immune and interstitium"
                                  ,"PEC/LOH" = "Epithelial"
                                  ,"POD" = "Glomerulus"
                                  ,"T Cell" = "Immune and interstitium"
                                  ,"T-CYT" = "Immune and interstitium"
                                  ,"T-CYT-MEM" = "Immune and interstitium"
                                  ,"EPC" = "Epithelial"
                                  ,"IMM" = "Immune and interstitium"
                                  ,"vSMC/P" = "Vascular"
                                  ,"INT" = "Interstitium"
                                  ,"FIB" = "Interstitium"
                                  ,"MAC" = "Immune system"
                                  ,"vSMC" = "Vascular"
                                  )
  cellTypes_for_segments = names(orderedCellType_segment_list)
  cellTypes_for_figureGroup = names(orderedCellType_figureGroup_list)
  cellTypes_in_enrichment = unique(enrichment$Cell_type)

  if ( (length(which(!cellTypes_for_segments %in% cellTypes_for_figureGroup))!=0)
      |(length(which(!cellTypes_for_figureGroup %in% cellTypes_for_segments))!=0))
  { rm(enrichment); error_message = paste(function_name," 1",sep='') }

  if (length(which(!cellTypes_in_enrichment %in% cellTypes_in_enrichment))!=0)
  { rm(enrichment); error_message = paste(function_name," 2",sep='') }

  cellTypes = cellTypes_for_segments

  enrichment$Subsegment = "error"
  enrichment$Figure_group = "error"

  for (indexST in 1:length(cellTypes))
  {#Begin
     cellType = cellTypes[indexST]
     segment = orderedCellType_segment_list[[cellType]]
     figureGroup = orderedCellType_figureGroup_list[[cellType]]
     indexCurrentCellType = which(enrichment$Cell_type==cellType)
     if (length(indexCurrentCellType)>0)
     {#Begin
        enrichment$Subsegment[indexCurrentCellType] = segment
        enrichment$Figure_group[indexCurrentCellType] = figureGroup
    }#End
  }#End

  indexError = which(enrichment$Subsegment=="error")
  unique(enrichment$Cell_type[indexError])
  if (length(indexError)>0)  { rm(enrichment); error_message = paste(function_name," 3",sep='') }
}#End - Add cell types, segments and figure groups to enrichment lines

function_name = "Set center cell subtype order based on cell type and decreasing sodium transport activity"
{#Begin - Set center cell subtype order based on cell type and decreasing sodium transport activity
   center_cell_subtype_orders = list()
   center_cell_subtype_orders[["LMD Proteomics"]]=c("Glom","TI")
   center_cell_subtype_orders[["NSC Proteomics"]]=c("Glom","PT")
   center_cell_subtype_orders[["LMD RNASeq"]]=c("Glom","PT","TAL","DCT","CD","Int")
   center_cell_subtype_orders[["SC RNASeq"]]=c("GC-EC","POD","vSMC/MC","PEC/LOH","PT-5","PT-1","PT-2","PT-3","PT-6","PT-7","PT-4","DTL","ATL","TAL","DCT","CNT","CNT-PC","PC","tPC-IC","IC-A","IC-B","MAC","MON","NKC","B Cell","T-CYT","T-CYT-MEM","T Cell","EC-AEA","EC-PT","FIB")
   center_cell_subtype_orders[["SN RNASeq"]]=c("POD","MC","EPC","PT-1","PT-4","PT-5","PT-2","DTL","ATL-2","ATL-1","ATL-3","TAL-1","TAL-2","DCT","CNT","PC-1","IC-A1","IC-B","PC-3","IC-A2","EC-1","EC-2","EC-3","EC-4","vSMC/P","IMM","INT")
}#End - Set center cell subtype order based on cell type and decreasing sodium transport activity

########################### Group reabsorption profiles
{#Begin - Plot stacked bardiagrams
  remove_legend_from_individual_plots=TRUE


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
                         ,"Ketone body catabolism" = "orange2"
                         ,"Ketone body catabolism specific enzymes" = "orange2"
                         ,"Ketone body formation and catabolism shared enzymes" = "bisque3"
                         ,"Pyruvate dehydrogenase" = "orangered2"
                         ,"Aerobic glycolysis" = "orangered2"
                         ,"Anaerobic glycolysis" = "dodgerblue3"
                         ,"Glycolysis" = "wheat2"
                         ,"Glycolysis specific enzymes" = "orange"
                         ,"Glycolysis and gluconeogenesis shared enzymes" = "bisque3"
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

stacked_colors = list( #"Beta oxidation" = c("darkorchid1","orange4","orange1","cornflowerblue")
                      #,"Beta oxidation II" = c("orange","cornflowerblue","darkorchid","blue")
                       "Vesicle exocytosis" = c("orange","cornflowerblue","darkorchid","blue","firebrick2","green","plum")
                      ,"Vesicle exocytosis II" = c("orange","cornflowerblue","darkorchid","blue","firebrick2","green")
                      #,"Ion transport" = c("orange","blue")
                      ,"Ion channels" = c("darkorchid1","skyblue","orange","green","cornflowerblue")
                      ,"Overall energy I" = c("navy","blue","lightskyblue","orange","firebrick2","mediumorchid","pink")
                      ,"Overall energy II" = c("navy","blue","orange","firebrick2","mediumorchid","pink")
                      ,"NaKATPase" = c("cornflowerblue","orange","mediumorchid")
                      ,"Small molecule transport" = c("cornflowerblue","orange")
                      ,"Mitochondrial energy" = c("cornflowerblue","orange")
                      ,"Urea cycle" = c("cornflowerblue")
                      ,"Glucose metabolism" = c("navy","blue","lightskyblue","orange","mediumorchid","cornflowerblue","firebrick2")
                      ,"Gluconeogenesis substrates" = c("cornflowerblue","orange","blue","firebrick","darkorchid","yellow","purple")
                      ,"Renal elimination of drugs and toxins" = c("cornflowerblue","orange","blue")
                      ,"GPI" = c("cornflowerblue","orange","blue")
)

scp_dependencies_list = list(  "Anaerobic glycolysis" = c("Glycolysis only enzymes","Lactate dehydrogenase")
                         ,"Aerobic glycolysis" = c("Glycolysis only enzymes","Pyruvate dehydrogenase")
                         ,"Peroxisomal beta oxidation" = c("Beta oxidation only peroxisomal enzymes")
                         ,"Mitochondrial beta oxidation" = c("Beta oxidation only mitochondrial enzymes")
                         ,"Glycolysis" = c("Glycolysis only enzymes")
                         ,"Keton body catabolism" = c("Keton body catabolism specific enzymes")
                       )
depending_scps = names(scp_dependencies_list)
scp_exclusion_list = list("Glycolysis" = c("Lactate dehydrogenase","Pyruvate dehydrogenase"))
scps_for_exclusion = names(scp_exclusion_list)

figure_groups = "Nephron"#unique(enrichment$Figure_group)
indexF = 1
for (indexF in 1:length(figure_groups))
{#Begin
  figure_group = figure_groups[indexF]
  indexCurrentFigureGroup = which(enrichment$Figure_group==figure_group)
  current_enrichment = enrichment[indexCurrentFigureGroup,]

  Stacked_bardiagrams = list()

  center_group_max_minusLog10Pvalue = list()
  if (figure_group == "HNephron")
  {#Begin - Set maximum values
     premiere_group_max_minusLog10Pvalue = list("Beta oxidation II" = 13,
                                                "Beta oxidation" = 13,
                                                "Glucose metabolism" = 13,
                                                "Oxidative phosphorylation" = 13,
                                                "Substrates for gluconeogenesis" = 13,
                                                "Urea" = 13
                                                )
     ucsd_group_max_minusLog10Pvalue = list("Beta oxidation" = 10,
                                            "Beta oxidation II" = 10,
                                            "Glucose metabolism" = 10,
                                            "Oxidative phosphorylation" = 10,
                                            "Substrates for gluconeogenesis" = 10,
                                            "Urea" = 10
                                           )
     lmd_group_max_minusLog10Pvalue = list("Beta oxidation II" = 2.5,
                                           "Beta oxidation" = 2.5,
                                           "Glucose metabolism" = 2.5,
                                           "Oxidative phosphorylation" = 2.5,
                                           "Substrates for gluconeogenesis" = 2.5,
                                           "Urea" = 2.5
                                          )

     center_group_max_minusLog10Pvalue = list("SC RNASeq" = premiere_group_max_minusLog10Pvalue,
                                              "SN RNASeq" = ucsd_group_max_minusLog10Pvalue,
                                              "LMD RNASeq" = lmd_group_max_minusLog10Pvalue
     )
  }#End - Set maximum values

  indexGroup = 5
  for (indexGroup in 1:length(SCP_groups))
  {#Begin
     current_enrichment_local = current_enrichment
     current_group = names(SCP_groups)[indexGroup]
     current_group_scps = SCP_groups[[current_group]]
     legend_scps_count = 0

     if (current_group %in% c("Metabolic summary","Metabolic summary II"))
     {#Begin - Remove pathways that depend on child pathways
       centers = unique(current_enrichment_local$Center)
       indexCenter=4
       new_current_enrichment_local = c()
       for (indexCenter in 1:length(centers))
       {#Begin
         center = centers[indexCenter]
         indexCurrentCenter = which(current_enrichment_local$Center==center)
         center_enrichment = current_enrichment_local[indexCurrentCenter,]
         cell_subtypes = unique(center_enrichment$Cell_subtype)
         indexCS=11
         for (indexCS in 1:length(cell_subtypes))
         {#Begin
            cell_subtype = cell_subtypes[indexCS]
            indexCurrentCellSubtype = which(center_enrichment$Cell_subtype==cell_subtype)
            cell_subtype_enrichment = center_enrichment[indexCurrentCellSubtype,]
            indexRemove=c()
            indexDS=1
            for (indexDS in 1:length(depending_scps))
            {#Begin
              depending_scp = depending_scps[indexDS]
              scps_that_need_to_be_there = scp_dependencies_list[[depending_scp]]
              indexScps_need_to_be_there = which(cell_subtype_enrichment$Subcellular.process..SCP %in% scps_that_need_to_be_there)
              if (length(cell_subtype_enrichment$Subcellular.process..SCP[indexScps_need_to_be_there])!=length(unique(cell_subtype_enrichment$Subcellular.process..SCP[indexScps_need_to_be_there])))
              { cell_subtype_enrichment= c() }
              if (min(cell_subtype_enrichment$Minus_log10_pvalue[indexScps_need_to_be_there])==0)
              {#Begin
                indexRemove_add = which(cell_subtype_enrichment$Subcellular.process..SCP==depending_scp)
                indexRemove = c(indexRemove,indexRemove_add)
              }#End
            }#End
            indexEx=1
            for (indexEx in 1:length(scps_for_exclusion))
            {#Begin
              scp_for_exclusion = scps_for_exclusion[indexEx]
              scps_that_must_not_be_there = scp_exclusion_list[[scp_for_exclusion]]
              indexScps_must_not_be_there = which(cell_subtype_enrichment$Subcellular.process..SCP %in% scps_that_must_not_be_there)
              if (length(cell_subtype_enrichment$Subcellular.process..SCP[indexScps_must_not_be_there])!=length(unique(cell_subtype_enrichment$Subcellular.process..SCP[indexScps_must_not_be_there])))
              { cell_subtype_enrichment= c() }
              if (sum(cell_subtype_enrichment$Minus_log10_pvalue[indexScps_must_not_be_there])>0)
              {#Begin
                indexRemove_add = which(cell_subtype_enrichment$Subcellular.process..SCP==scp_for_exclusion)
                indexRemove = c(indexRemove,indexRemove_add)
              }#End
            }#End
            new_cell_subtype_enrichment = cell_subtype_enrichment
            new_cell_subtype_enrichment$Minus_log10_pvalue[indexRemove] = 0;
            if (length(new_current_enrichment_local)==0) { new_current_enrichment_local = new_cell_subtype_enrichment }
            else { new_current_enrichment_local = rbind(new_current_enrichment_local,new_cell_subtype_enrichment) }
         }#End
       }#End
       current_enrichment_local = new_current_enrichment_local
     }#End - Remove pathways that depend on child pathways

     indexCurrentGroup = which(current_enrichment_local$Subcellular.process..SCP %in% current_group_scps)
     current_group_enrichment = current_enrichment_local[indexCurrentGroup,]

     ontology = current_group_enrichment$Ontology[1]
     current_centers = selected_technologies#c("SC RNASeq","SN RNASeq","LMD RNASeq") #,"SC RNASeq Liao et al"
     current_centers = current_centers[current_centers %in% unique(current_enrichment$Center)]
     indexCenter = 5
     stacked_legend = ggplot() + theme_void()
     stacked_legend_scps_count = 0
     for (indexCenter in 1:length(current_centers))
     {#Begin
        current_center = current_centers[indexCenter]
        center_cell_subtype_order = c(center_cell_subtype_orders[[current_center]],"empty")
        center_cell_subtype_colors = center_cell_subtype_colors_list[[current_center]]
        indexCurrentCenter_enrichment = which(current_group_enrichment$Center==current_center)
        indexNonZero = which(current_group_enrichment$Minus_log10_pvalue!=0)
        indexCurrentCenter_enrichment_nonZero = indexCurrentCenter_enrichment[indexCurrentCenter_enrichment %in% indexNonZero]
        if (length(indexCurrentCenter_enrichment_nonZero)==0) { Stacked_bardiagrams[[length(Stacked_bardiagrams)+1]] = ggplot() + theme_void() }
        if (length(indexCurrentCenter_enrichment_nonZero)>0)
        {#Begin
           center_currentGroup_collapsed_enrichment = current_group_enrichment[indexCurrentCenter_enrichment,]
           center_currentGroup_collapsed_enrichment$Cell_subtype = factor(center_currentGroup_collapsed_enrichment$Cell_subtype,levels=center_cell_subtype_order)
           plot_data = center_currentGroup_collapsed_enrichment;

           plot_data$Scp_factor = factor(plot_data$Subcellular.process..SCP,levels=current_group_scps)
           plot_data = plot_data[order(plot_data$Scp_factor),]
           current_plot_data_scps = unique(plot_data$Subcellular.process..SCP)
           stacked_colors = replicate(length(current_plot_data_scps),"gray")
           for (indexCG in 1:length(current_plot_data_scps))
           {#Begin
              current_group_scp = current_plot_data_scps[indexCG]
              if (current_group_scp %in% names(scp_color_list))
              {#Begin
                 stacked_colors[indexCG] = scp_color_list[[current_group_scp]]
              }#End
           }#End

           cell_subtypes_for_color = center_cell_subtype_order[center_cell_subtype_order %in% plot_data$Cell_subtype]

           text_colors=c()
           for (indexDataForColor in 1:length(cell_subtypes_for_color))
           {#Begin
              text_colors = c(text_colors,center_cell_subtype_colors[[cell_subtypes_for_color[indexDataForColor]]])
           }#End

           headline = paste(current_group,"\n",current_center,sep='')

           Stacked_bardiagram = ggplot(data=plot_data,aes(x=Cell_subtype,y=Minus_log10_pvalue,group=Scp_factor,fill=Scp_factor))
           Stacked_bardiagram = Stacked_bardiagram + geom_bar(stat="identity",color="black")
           Stacked_bardiagram = Stacked_bardiagram + ggtitle(headline) + theme(title = element_text(size=20,face=2,vjust=0.5))
           Stacked_bardiagram = Stacked_bardiagram + theme(axis.text.x = element_blank(),axis.ticks=element_blank())
           Stacked_bardiagram = Stacked_bardiagram + xlab("") + ylab("")
           Stacked_bardiagram = Stacked_bardiagram + scale_fill_manual(values=stacked_colors)
           Stacked_bardiagram = Stacked_bardiagram + theme(axis.text.x = element_text(angle=90,size=15,face=2,vjust=0.5,hjust=1,color=text_colors))
           Stacked_bardiagram = Stacked_bardiagram + theme(axis.text.y = element_text(size=17,face=2))

           if (length(unique(plot_data$Subcellular.process..SCP))>stacked_legend_scps_count)
           {#Begin
             current_stacked_legend = cowplot::get_legend(Stacked_bardiagram)
             stacked_legend = current_stacked_legend
             stacked_legend_scps_count = length(unique(plot_data$Subcellular.process..SCP))
           }#End

           if (remove_legend_from_individual_plots)
           {#Begin
              Stacked_bardiagram = Stacked_bardiagram + theme(legend.position = "none")
              Stacked_bardiagram = Stacked_bardiagram + guides(fill=FALSE)
           }#End
           if (current_center %in% names(center_group_max_minusLog10Pvalue))
           {#Begin
              group_max_minusLog10Pvalue = center_group_max_minusLog10Pvalue[[current_center]]
              if (current_group %in% names(group_max_minusLog10Pvalue))
              { Stacked_bardiagram = Stacked_bardiagram + ylim(0, group_max_minusLog10Pvalue[[current_group]]) }
           }#End
           Stacked_bardiagrams[[length(Stacked_bardiagrams)+1]] = Stacked_bardiagram
        }#End
     }#End
     #Stacked_bardiagrams[[length(Stacked_bardiagrams)+1]] = stacked_legend
  }#End

  Layout_matrix_base = rbind(c( 1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,4,4,4,5,5,5,6,6,6,6,6,6,6,6),
                             c( 1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,4,4,4,5,5,5,6,6,6,6,6,6,6,6),
                             c( 1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,4,4,4,5,5,5,6,6,6,6,6,6,6,6))
  Layout_matrix_base = rbind(c( 1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,4,4,4,5,5,5),
                             c( 1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,4,4,4,5,5,5),
                             c( 1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,4,4,4,5,5,5))
  Layout_matrix = c()

  for (indexSCP in 1:length(names(SCP_groups)))
  {#Begin
    add_layout_matrix = Layout_matrix_base + (indexSCP-1) * max(Layout_matrix_base)
    if (length(Layout_matrix)==0) { Layout_matrix = add_layout_matrix }
    else { Layout_matrix = rbind(Layout_matrix,add_layout_matrix) }
  }#End

  complete_png_fileName = paste(results_directory,fileName_metabolism_png,sep='')
  cols_count = min(length(current_centers)+1,length(Stacked_bardiagrams))
  cols_count = min(length(current_centers),length(Stacked_bardiagrams))
  rows_count = ceiling(length(Stacked_bardiagrams)/cols_count)
  if (remove_legend_from_individual_plots)  { width_per_plot = 350 }
  if (!remove_legend_from_individual_plots) { width_per_plot = 900 }

  png(complete_png_fileName,width=width_per_plot*cols_count,height=400*rows_count,res=100);
  #grid.arrange(Stacked_bardiagrams[[1]],Stacked_bardiagrams[[2]],Stacked_bardiagrams[[3]],Stacked_bardiagrams[[4]],Stacked_bardiagrams[[5]],Stacked_bardiagrams[[6]],
  #             Stacked_bardiagrams[[7]],Stacked_bardiagrams[[8]],Stacked_bardiagrams[[9]],Stacked_bardiagrams[[10]],Stacked_bardiagrams[[11]],Stacked_bardiagrams[[12]],
  #             Stacked_bardiagrams[[13]],Stacked_bardiagrams[[14]],Stacked_bardiagrams[[15]],Stacked_bardiagrams[[16]],Stacked_bardiagrams[[17]],Stacked_bardiagrams[[18]],
  #             Stacked_bardiagrams[[19]],Stacked_bardiagrams[[20]],Stacked_bardiagrams[[21]],Stacked_bardiagrams[[22]],Stacked_bardiagrams[[23]],Stacked_bardiagrams[[24]],
  #             Stacked_bardiagrams[[25]],Stacked_bardiagrams[[26]],Stacked_bardiagrams[[27]],Stacked_bardiagrams[[28]],Stacked_bardiagrams[[29]],Stacked_bardiagrams[[30]],
  #             layout_matrix = Layout_matrix)
  grid.arrange(Stacked_bardiagrams[[1]],Stacked_bardiagrams[[2]],Stacked_bardiagrams[[3]],Stacked_bardiagrams[[4]],Stacked_bardiagrams[[5]],
               Stacked_bardiagrams[[6]],Stacked_bardiagrams[[7]],Stacked_bardiagrams[[8]],Stacked_bardiagrams[[9]],Stacked_bardiagrams[[10]],
               Stacked_bardiagrams[[11]],Stacked_bardiagrams[[12]],Stacked_bardiagrams[[13]],Stacked_bardiagrams[[14]],Stacked_bardiagrams[[15]],
               Stacked_bardiagrams[[16]],Stacked_bardiagrams[[17]],Stacked_bardiagrams[[18]],Stacked_bardiagrams[[19]],Stacked_bardiagrams[[20]],
               Stacked_bardiagrams[[21]],Stacked_bardiagrams[[22]],Stacked_bardiagrams[[23]],Stacked_bardiagrams[[24]],Stacked_bardiagrams[[25]],
               layout_matrix = Layout_matrix)
  dev.off()

}#End - Plot grouped reabsorption profiles

}#End - Plot stacked bardiagrams
########################### Group reabsorption profiles


# Change back to start-up directory
if(!is.null(start_dir)) setwd(start_dir)
