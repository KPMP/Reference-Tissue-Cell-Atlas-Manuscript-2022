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

{#Begin - Open libraries and document versions
Col_names = c("Library","Version")
Col_length = length(Col_names)
Row_names = 1
Row_length= length(Row_names)
version_documentation_line = array(NA,c(Row_length,Col_length),dimnames=list(Row_names,Col_names))
version_documentation_line = as.data.frame(version_documentation_line)
version_documentations = c()

libraries = c("SingleCellExperiment","Seurat","ggplot2","sctransform","mclust","dplyr","beeswarm","clustree","Matrix","gridExtra","reticulate","umap","doParallel")
for (indexL in 1:length(libraries))
{#Begin
  current_library = libraries[indexL]
  library(current_library,character.only=TRUE)
  new_version_documentation_line = version_documentation_line
  new_version_documentation_line$Library = current_library
  new_version_documentation_line$Version = packageVersion(current_library)
  if (length(version_documentations)==0)
  {#Begin
    version_documentations = new_version_documentation_line
  }#End
  else
  {#Begin
    version_documentations = rbind(version_documentations,new_version_documentation_line)
  }#End
}#End

r_sessionInfo = sessionInfo()
python_sessionInfo = py_config()
}#End - Open libraries and document versions


#Region specify files
base_directory = paste(getwd(),"/",sep='');
base_directory = "D:/KPMP_reference_atlas_code/"
output_directory = paste(base_directory,"Experimental_data//IntegratedCluster_AvgExpression_RNAcounts//",sep='')
dir.create(output_directory)

{#Begin - Read seurat object
  seurat_directory = paste(base_directory,"Experimental_data//SingleCell_datasets//",sep='')
  seurat_complete_fileName = paste(seurat_directory,"scsn.integrated.013122.Robj",sep='')
  load(seurat_complete_fileName,verbose = TRUE)
}#End - Read seurat object

{#Begin - Read metadata
  metadata_directory = paste(base_directory,"//Experimental_data//Sample_metadata_additional_datasets//",sep='')
  premiere_metadata_fileName = "Single_cell_RNASeq_PREMIERE_metadata.txt";
  ucsd_metadata_fileName = "Single_nucleus_RNASeq_UCSDWU_metadata.txt";
  premiere_complete_metadata_fileName = paste(metadata_directory,premiere_metadata_fileName,sep='')
  ucsd_complete_metadata_fileName = paste(metadata_directory,ucsd_metadata_fileName,sep='')
  ucsf_metadata_fileName = "Single_cell_RNASeq_UCSF_metadata.txt"
  ucsf_complete_metadata_fileName = paste(metadata_directory,ucsf_metadata_fileName,sep='')

  keep_colnames = c("Library","Patient_id","Center")

  premiere_metadata = read.table(file=premiere_complete_metadata_fileName,header = TRUE,stringsAsFactors = FALSE,sep='\t')
  premiere_center_label = "Premiere-"
  premiere_metadata$Library = paste(premiere_center_label,premiere_metadata$Library,sep='')
  premiere_metadata$Patient_id = paste(premiere_center_label,premiere_metadata$Patient_id,sep='')
  premiere_metadata$Center = "SC RNASeq I" #Has to match the dataset names in C# script defined in the class KPMP_dataset_name_class

  indexColKeep = which(colnames(premiere_metadata) %in% keep_colnames)
  premiere_metadata = premiere_metadata[,indexColKeep]
  premiere_metadata$Library = gsub("[.]1","-1",premiere_metadata$Library)
  ucsd_metadata = read.table(file=ucsd_complete_metadata_fileName,header = TRUE,stringsAsFactors = FALSE,sep='\t')
  ucsd_center_label = "UCSD-"
  ucsd_metadata$Library = paste(ucsd_center_label,ucsd_metadata$Library,sep='')
  ucsd_metadata$Patient_id = paste(ucsd_center_label,ucsd_metadata$Patient_id,sep='')
  ucsd_metadata$Center = "SN RNASeq" #Has to match the dataset names in C# script defined in the class KPMP_dataset_name_class
  indexColKeep = which(colnames(ucsd_metadata) %in% keep_colnames)
  ucsd_metadata = ucsd_metadata[,indexColKeep]
  ucsf_metadata = read.table(file=ucsf_complete_metadata_fileName,header = TRUE,stringsAsFactors = FALSE,sep='\t')
  ucsf_center_label = "UCSF-"
  ucsf_metadata$Library = paste(ucsf_center_label,ucsf_metadata$Library,sep='')
  ucsf_metadata$Patient_id = paste(ucsf_center_label,ucsf_metadata$Patient_id,sep='')
  ucsf_metadata$Center = "SC RNASeq II" #Has to match the dataset names in C# script defined in the class KPMP_dataset_name_class
  indexColKeep = which(colnames(ucsf_metadata) %in% keep_colnames)
  ucsf_metadata = ucsf_metadata[,indexColKeep]
  ucsf_metadata$Library = gsub("UCSF-18-139","UCSF-18-139-1",ucsf_metadata$Library)
  ucsf_metadata$Library = gsub("UCSF-18-142","UCSF-18-142-6",ucsf_metadata$Library)
  ucsf_metadata$Library = gsub("UCSF-18-162","UCSF-18-162-6",ucsf_metadata$Library)
  ucsf_metadata$Library = gsub("UCSF-18-312","UCSF-18-312-6",ucsf_metadata$Library)
  ucsf_metadata$Library = gsub("UCSF-18-342","UCSF-18-342-6",ucsf_metadata$Library)
  metadata = rbind(premiere_metadata,ucsf_metadata,ucsd_metadata)
}#End - Read metadata

{#Begin - Add participants and centers to seurat object
  libraries = unique(tis.integrated$orig.ident)
  tis.integrated$Participant = "error"
  tis.integrated$Center = "error"
  for (indexL in 1:length(libraries))
  {#Begin
     current_library = libraries[indexL]
     indexMetadata = which(metadata$Library==current_library)
     indexTis = which(tis.integrated$orig.ident==current_library)
     tis.integrated$Participant[indexTis] = metadata$Patient_id[indexMetadata]
     tis.integrated$Center[indexTis] = metadata$Center[indexMetadata]
  }#End
  indexError = which(tis.integrated$Participant=="error")
  if (length(indexError)>0) { rm(tis.integrated) }
}#End - Add participants and centers to seurat object

{#Begin - Calculate and write averages

tis.integrated = subset(tis.integrated, subset = celltype %in% c("PT","POD"))

avg_slot = "counts"
avg_data = "RNA"


participants = unique(tis.integrated$Participant)
centers = unique(tis.integrated$Center)
indexP=1
combined_average_expression = c()
for (indexP in 1:length(participants))
{#Begin
   current_participant = participants[indexP]
   current_participant_tis_integrated = subset(x = tis.integrated, subset = (Participant == current_participant))
   centers = unique(current_participant_tis_integrated$Center)
   for (indexCenter in 1:length(centers))
   {#Begin
      current_center = centers[indexCenter]
      current_participant_center_tis_integrated = subset(x=current_participant_tis_integrated, subset = (Center == current_center))
      cellTypes = unique(current_participant_center_tis_integrated$celltype)
      for (indexCellType in 1:length(cellTypes))
      {#Begin
         current_cellType = cellTypes[indexCellType]
         current_cellType_tis_integrated = subset(x=current_participant_center_tis_integrated, subset = (celltype == current_cellType))

         current_average_expression = AverageExpression(current_cellType_tis_integrated,slot=avg_slot,assays=avg_data)
         current_participant_document = current_participant
         current_average_expression = as.data.frame(current_average_expression$RNA)
         current_participant_document = gsub(premiere_center_label,"",current_participant_document)
         current_participant_document = gsub(ucsd_center_label,"",current_participant_document)
         current_participant_document = gsub(ucsf_center_label,"",current_participant_document)
         current_average_expression$Cell_type = current_cellType
         current_average_expression$Patient = current_participant_document
         current_average_expression$Dataset = current_center
         current_average_expression$Seurat_data = avg_data
         current_average_expression$Seurat_slot = avg_slot
         current_average_expression$Symbol = rownames(current_average_expression)

         if (length(combined_average_expression)==0) { combined_average_expression = current_average_expression }
         else { combined_average_expression = rbind(combined_average_expression,current_average_expression) }
      }#End
   }#End
}#End

indexValue = which(colnames(combined_average_expression)=="all")
colnames(combined_average_expression)[indexValue] = "Expression_value"

combined_participant_center_complete_fileName = paste(output_directory,"AverageRNAExpression_",avg_data,"_",avg_slot,".txt",sep='')
write.table(combined_average_expression,file=combined_participant_center_complete_fileName,row.names=FALSE,col.names=TRUE,sep='\t',quote = FALSE)

}#End - Calculate and write averages


# Change back to start-up directory
if(!is.null(start_dir)) setwd(start_dir)
