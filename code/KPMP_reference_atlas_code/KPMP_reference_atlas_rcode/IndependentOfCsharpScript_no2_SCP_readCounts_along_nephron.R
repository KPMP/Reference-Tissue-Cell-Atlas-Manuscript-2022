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

##############################################################################################################################################

#Here specify the location of the main directory
main_directory = paste(getwd(),"/",sep='');#"D:/KPMP_reference_atlas_code/"
main_directory = "D:/KPMP_reference_atlas_code/"
#Here specify the location of the main directory

##############################################################################################################################################

{#Begin - Open libraries and document versions
Col_names = c("Library","Version")
Col_length = length(Col_names)
Row_names = 1
Row_length= length(Row_names)
version_documentation_line = array(NA,c(Row_length,Col_length),dimnames=list(Row_names,Col_names))
version_documentation_line = as.data.frame(version_documentation_line)
version_documentations = c()

#libraries = c("ggVennDiagram","cowplot","SingleCellExperiment","scales","RColorBrewer","Seurat","SeuratDisk","biomaRt","harmony","gridExtra","ggplot2","ggridges","sctransform","mclust","dplyr","beeswarm","clustree","Matrix","gridExtra","reticulate","umap","stringr","rhdf5","progress")
libraries = c("ggVennDiagram","cowplot","SingleCellExperiment","scales","RColorBrewer","Seurat","biomaRt","harmony","gridExtra","ggplot2","ggridges","sctransform","mclust","dplyr","beeswarm","clustree","Matrix","gridExtra","reticulate","umap","stringr","rhdf5","progress")
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

#use_python("c://ProgramData//Anaconda3//")
r_sessionInfo = sessionInfo()
#python_sessionInfo = py_config()
}#End - Open libraries and document versions

singleCellData_directory = paste(main_directory,"Experimental_data//SingleCell_datasets//",sep='')
results_base_directory = paste(main_directory,"Results//Sodium_glucose_reabsorption//",sep='')
dir.create(results_base_directory)

library_of_interest = c("NaAndGluTMTransport")

read_premiere = TRUE
if (read_premiere)
{#Begin - Read and process premiere data
   rawData_fileName = "Premiere_Raw.data_062220.txt"

   complete_rawDataFileName = paste(singleCellData_directory,rawData_fileName,sep='')
   list.files(complete_rawDataFileName)
   Premiere_raw_data_matrix_input = read.table(file=complete_rawDataFileName, header=TRUE, check.names=FALSE)
   rownames(Premiere_raw_data_matrix_input) = gsub("_","-",rownames(Premiere_raw_data_matrix_input))

   annotation_fileName = "lcm_cell_assignment_barcodes_premiere_2021April22_withXs.txt";
   complete_annotationFileName = paste(singleCellData_directory,annotation_fileName,sep='')
   premiere_annotations = read.table(file=complete_annotationFileName, header=TRUE, stringsAsFactors = FALSE,check.names=FALSE,sep='\t')
   indexPCIC = which(premiere_annotations$Cluster=="tPC-IC")
   premiere_annotations$Cluster[indexPCIC] = "IC"
   indexPCCNT = which(premiere_annotations$Cluster=="PC-CNT")
   premiere_annotations$Cluster[indexPCCNT] = "PC"
   premiere_annotations$Library = "error"
   for (indexPA in 1:length(premiere_annotations[,1]))
   {#Begin
      splitStrings = str_split(premiere_annotations$Barcode[indexPA],"_")[[1]]
      premiere_annotations$Library[indexPA] = splitStrings[1]
      premiere_annotations$Library[indexPA] = gsub("X","",premiere_annotations$Library[indexPA])
   }#End
   premiere_annotations$Barcode = gsub("-",".",premiere_annotations$Barcode)

   metadata_fileName = "Single_cell_RNASeq_PREMIERE_metadata.txt"
   complete_metadataFileName = paste(singleCellData_directory,metadata_fileName,sep='')
   premiere_metadata = read.table(file=complete_metadataFileName, header=TRUE, check.names=FALSE,stringsAsFactors = FALSE, sep='\t')
   premiere_metadata$Center = "SC RNASeq PREMIERE"

   indexMissing = which(!premiere_annotations$Library %in% premiere_metadata$Library)
   libraries = unique(premiere_annotations$Library)
   premiere_annotations$Participant = "error"
   premiere_annotations$XIST = "error"
   premiere_annotations$Tissue_type = "eror"
   premiere_annotations$Tissue_collection = "error"
   premiere_annotations$Sex_metadata = "unknown"
   premiere_annotations$Age = "unknown"
   for (indexL in 1:length(libraries))
   {#Begin
      library = libraries[indexL]
      indexLibraryAnnotations = which(premiere_annotations$Library==library)
      indexLibraryMetadata = which(premiere_metadata$Library==library)
      premiere_annotations$Participant[indexLibraryAnnotations] = premiere_metadata$Patient_id[indexLibraryMetadata]
      premiere_annotations$Tissue_type[indexLibraryAnnotations] = premiere_metadata$Tissue_type[indexLibraryMetadata]
      premiere_annotations$XIST[indexLibraryAnnotations] = premiere_metadata$XIST[indexLibraryMetadata]
      premiere_annotations$Tissue_collection[indexLibraryAnnotations] = premiere_metadata$Tissue_collection[indexLibraryMetadata]
      premiere_annotations$Sex_metadata[indexLibraryAnnotations] = premiere_metadata$Sex[indexLibraryMetadata]
      premiere_annotations$Age[indexLibraryAnnotations] = premiere_metadata$Age[indexLibraryMetadata]
   }#End

   indexKeepRow_premiere = which(rowSums(Premiere_raw_data_matrix_input)!=0)
   Premiere_raw_data_matrix_input = Premiere_raw_data_matrix_input[indexKeepRow_premiere,]

   premiere_pubmedId = "PMID: 32107344"
}#End - Read and process premiere data

read_ucsd = TRUE
if (read_ucsd)
{#Begin - Read and process ucsd data
  rawData_fileName = "GSE121862_UCSD-WU_Single_Nuclei_Cluster_Annotated_Raw_UMI_Matrix.tsv"

  complete_rawDataFileName = paste(singleCellData_directory,rawData_fileName,sep='')
  UCSD_raw_data_matrix_input = read.table(file=complete_rawDataFileName, header=TRUE, check.names=FALSE)
  rownames(UCSD_raw_data_matrix_input) = gsub("_","-",rownames(UCSD_raw_data_matrix_input))

  annotation_fileName = "lcm_cell_assignment_barcodes_ucsd_2021April22.txt";
  complete_annotationFileName = paste(singleCellData_directory,annotation_fileName,sep='')
  ucsd_annotations = read.table(file=complete_annotationFileName, header=TRUE, stringsAsFactors = FALSE,check.names=FALSE,sep='\t')
  indexPCCNT = which(ucsd_annotations$Cluster=="DL")
  ucsd_annotations$Cluster[indexPCCNT] = "DTL"

  metadata_fileName = "Single_nucleus_RNASeq_UCSDWU_metadata.txt"
  complete_metadataFileName = paste(singleCellData_directory,metadata_fileName,sep='')
  ucsd_metadata = read.table(file=complete_metadataFileName, header=TRUE, check.names=FALSE,stringsAsFactors = FALSE, sep='\t')
  ucsd_metadata$Center = "SN RNASeq UCSDWU"
  ucsd_annotations$Participant = "error"
  ucsd_annotations$Tissue_type = "error"
  ucsd_annotations$Tissue_collection = "error"

  for (indexPA in 1:length(ucsd_annotations[,1]))
  {#Begin
    splitStrings = str_split(ucsd_annotations$Barcode[indexPA],"_")[[1]]
    library = splitStrings[2]
    ucsd_annotations$Library[indexPA] = library
    indexMetadata_library = which(ucsd_metadata$Library==library)
    if (length(indexMetadata_library)!=1) { error }
    ucsd_annotations$Participant[indexPA] = ucsd_metadata$Patient_id[indexMetadata_library]
    ucsd_annotations$Tissue_type[indexPA] = ucsd_metadata$Tissue_type[indexMetadata_library]
    ucsd_annotations$Tissue_collection[indexPA] = ucsd_metadata$Tissue_collection[indexMetadata_library]
  }#End
  indexMissing = which(!ucsd_annotations$Library %in% ucsd_metadata$Library)
  indexMissing = which(ucsd_annotations$Participant =="error")
  indexMissing = which(ucsd_annotations$Tissue_type =="error")

  {#Begin - Remove selected cell subtypes
    remove_ucsd_cell_subtypes = c("PT-3","PC-2","Unk")
    indexRemove = which(ucsd_annotations$Cluster  %in% remove_ucsd_cell_subtypes)
    indexKeep = 1:length(ucsd_annotations[,1])
    indexKeep = indexKeep[!indexKeep %in% indexRemove]
    ucsd_annotations = ucsd_annotations[indexKeep,]
    barcodes = ucsd_annotations$Barcode
    indexKeep = which(colnames(UCSD_raw_data_matrix_input) %in% barcodes)
    UCSD_raw_data_matrix_input = UCSD_raw_data_matrix_input[,indexKeep]
    libraries = unique(ucsd_annotations$Library)
    indexKeep = which(ucsd_metadata$Library %in% libraries)
    ucsd_metadata = ucsd_metadata[indexKeep,]
  }#End - Remove selected cell subtypes

  indexKeepRow_ucsd = which(rowSums(UCSD_raw_data_matrix_input)!=0)
  UCSD_raw_data_matrix_input = UCSD_raw_data_matrix_input[indexKeepRow_ucsd,]

  ucsd_pubmedId = "PMID: 31249312"
}#End - Read and process ucsd data

read_wu_data = TRUE
if (read_wu_data)
{#Begin - Read and process wu data - PMID29980650

  fileName = "Wu_MTS.kidney_PMID29980650.dge.txt"
  completeFileName = paste(singleCellData_directory,fileName,sep='')

  wu_raw_data_matrix_input = read.table(file=completeFileName,header=TRUE,stringsAsFactors = FALSE,sep='\t')

  barcode_fileName = "Barcodes_setaCC13_generated_from_WU_PMID29980650.txt"
  complete_barcode_fileName = paste(main_directory,"Results/WU_PMID29980650_singleCellAnalysis/",barcode_fileName,sep='')
  wu_annotations = read.table(file=complete_barcode_fileName,header=TRUE,stringsAsFactors = FALSE,sep='\t')

  cluster_cellType_list = list("0" = "PT",
                               "1" = "TAL",
                               "2" = "DTL",
                               "3" = "PC",
                               "4" = "PT",
                               "5" = "Unknown",
                               "6" = "PT",
                               "7" = "POD",
                               "8" = "IC",
                               "9" = "DCT",
                               "10" = "CNT",
                               "11" = "EC",
                               "12" = "Unknown")


  wu_annotations$Tissue_type = "Unknown"
  wu_annotations$Tissue_collection = "Unknown"
  wu_annotations$Sex = "Unknown"
  wu_annotations$Library = "AdultKidney"
  wu_annotations$Participant = "AdultKidney"


  clusters = names(cluster_cellType_list)
  for (indexClust in 1:length(clusters))
  {#Begin
     cluster = clusters[indexClust]
     indexCurrentCluster = which(wu_annotations$Cluster==cluster)
     wu_annotations$Cluster[indexCurrentCluster] = cluster_cellType_list[[cluster]]
  }#End

  wu_pubmedId = "PMID: 29980650"
}#End - Read and process wu data - PMID29980650

read_muto_data = TRUE
{#Begin - Read and process Muto - PMID: 33850129
   seurat_object_file = "Muto_singleNucleus_seurat_PMID33850129.rds"
   complete_seurate_object_file = paste(singleCellData_directory,seurat_object_file,sep='')
   muto_object <- readRDS(file = complete_seurate_object_file)

   muto_raw_data_matrix_input = as.matrix(muto_object@assays$RNA@counts)

   {#Begin - Convert ensembl ids to gene symbols and remove all rows that do not refer to genes in library of interest
     ensembl_ids = rownames(muto_raw_data_matrix_input)

     human_ensembl_mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
     gene_symbols <- getBM(filters="ensembl_gene_id",attributes= c("ensembl_gene_id","hgnc_symbol"),values=ensembl_ids,mart= human_ensembl_mart)
     gene_symbols = gene_symbols[order(gene_symbols$ensembl_gene_id),]
     indexesDuplicated = which(duplicated(gene_symbols$ensembl_gene_id))
     indexUnique = which(!duplicated(gene_symbols$ensembl_gene_id))
     for (indexIndex in 1:length(indexesDuplicated))
     {#Begin
        indexDuplicated = indexesDuplicated[indexIndex]
        if (gene_symbols$ensembl_gene_id[indexDuplicated-1]!=gene_symbols$ensembl_gene_id[indexDuplicated]) { rm(gene_symbols) }
        gene_symbols$hgnc_symbol[indexDuplicated-1] = paste(gene_symbols$hgnc_symbol[indexDuplicated-1],gene_symbols$hgnc_symbol[indexDuplicated],sep=';')
     }#End
     gene_symbols = gene_symbols[indexUnique,]
     muto_raw_data_matrix_input = as.data.frame(muto_raw_data_matrix_input)
     muto_raw_data_matrix_input$Ensembl_id = rownames(muto_raw_data_matrix_input)
     indexKeep = which(muto_raw_data_matrix_input$Ensembl_id %in% gene_symbols$ensembl_gene_id)
     muto_raw_data_matrix_input = muto_raw_data_matrix_input[indexKeep,]
     muto_raw_data_matrix_input = muto_raw_data_matrix_input[order(muto_raw_data_matrix_input$Ensembl_id),]
     indexMismatch = which(muto_raw_data_matrix_input$Ensembl_id!=gene_symbols$ensembl_gene_id)
     if (length(indexMismatch)>0) { rm(muto_raw_data_matrix_input)}
     muto_raw_data_matrix_input$Symbol = gene_symbols$hgnc_symbol

     {#Begin - Remove all rows that do not contain genes in library of interest (in this case no need for merging multiple rows with different ensembl IDs mapping to the same gene symbol)
       directory_libryaries = paste(main_directory,"MBCO_datasets//",sep='')
       complete_library_fileName = paste(directory_libryaries,library_of_interest,"_scpGeneAssociations.txt",sep='')
       library = read.csv(file=complete_library_fileName,stringsAsFactors = FALSE, sep='\t')
       all_genes_in_library = library$Target_gene_symbol
       indexKeep = which(muto_raw_data_matrix_input$Symbol %in% all_genes_in_library)
       muto_raw_data_matrix_input = muto_raw_data_matrix_input[indexKeep,]
       muto_raw_data_matrix_input = muto_raw_data_matrix_input[order(muto_raw_data_matrix_input$Symbol),]
       indexesDuplicated = which(duplicated(muto_raw_data_matrix_input$Symbol))
       if (length(indexesDuplicated)>0) { rm(muto_raw_data_matrix_input) }
       rm(directory_libryaries)
       rm(complete_library_fileName)
       rm(all_genes_in_library)
       rm(library)
     }#End - Remove all rows that do not contain genes in library of interest (in this case no need for merging multiple rows with different ensembl IDs mapping to the same gene symbol)

     indexColumnsWithNumbers = which(!colnames(muto_raw_data_matrix_input) %in% c("Symbol","Ensembl_id"))
     rownames(muto_raw_data_matrix_input) = muto_raw_data_matrix_input$Symbol

     muto_raw_data_matrix_input = as.matrix(muto_raw_data_matrix_input[,indexColumnsWithNumbers])
   }#End - Convert ensembl ids to gene symbols and remove all rows that do not refer to genes in library of interest

   Col_names = colnames(premiere_annotations)
   Col_length = length(Col_names)
   Row_names = 1:dim(muto_raw_data_matrix_input)[2]
   Row_length = length(Row_names)
   muto_annotations = as.data.frame(array(NA,c(Row_length,Col_length),dimnames=list(Row_names,Col_names)))

   muto_annotations$Barcode = names(muto_object$cell_type)
   muto_annotations$Cluster = as.character(muto_object$author_cell_type)
   muto_annotations$Library = as.character(muto_object$library_uuid)
   muto_annotations$Participant = as.character(muto_object$donor_uuid)
   muto_annotations$Tissue_type = as.character(muto_object$tissue)
   muto_annotations$Tissue_collection = as.character(muto_object$sample_preservation_method)
   muto_annotations$Sex = as.character(muto_object$sex)

   indexDCT = which((muto_annotations$Cluster=="DCT1")|(muto_annotations$Cluster=="DCT2"))
   muto_annotations$Cluster[indexDCT] = gsub("DCT","DCT-",muto_annotations$Cluster[indexDCT])
   indexIC = which((muto_annotations$Cluster=="ICA")|(muto_annotations$Cluster=="ICB"))
   muto_annotations$Cluster[indexIC] = gsub("IC","IC-",muto_annotations$Cluster[indexIC])
   indexPOD = which(muto_annotations$Cluster=="PODO")
   muto_annotations$Cluster[indexPOD] = "POD"
   indexPT_VCAM = which(muto_annotations$Cluster=="PT_VCAM1")
   muto_annotations$Cluster[indexPT_VCAM] = "PT-VCAM1"
   indexEndo = which(muto_annotations$Cluster=="ENDO")
   muto_annotations$Cluster[indexEndo] = "EC"
   indexMC = which(muto_annotations$Cluster=="MES")
   muto_annotations$Cluster[indexMC] = "MC"
   muto_pubmedId = "PMID: 33850129"
}#End - Read and process Muto - PMID: 33850129

loop_datasets = c("SingleNucleus","SingleCellAndNucleus")
indexDatasetLoop=1

for (indexDatasetLoop in 1:length(loop_datasets))
{#Begin - Calculate reabsorption profiles for selected datasets

dataset = loop_datasets[indexDatasetLoop]
algorithm_of_interest = "Total_sum"
running_plot_number = 0
error_message = ""

results_directory = paste(results_base_directory,dataset,"_",library_of_interest,"//",sep='')
dir.create(results_directory)

scpGroupSpecifi_ylims = list()

if (dataset=="SingleCellAndNucleus")
{#Begin
  annotations = list( "SC RNASeq PREMIERE" = premiere_annotations
                     ,"SN RNASeq UCSDWU" = ucsd_annotations
                     ,"SN RNASeq Muto" = muto_annotations
                     ,"SN RNASeq Wu" = wu_annotations
  )
  Raw_datas = list( "SC RNASeq PREMIERE" = Premiere_raw_data_matrix_input
                   ,"SN RNASeq UCSDWU" = UCSD_raw_data_matrix_input
                   ,"SN RNASeq Muto" = muto_raw_data_matrix_input
                   ,"SN RNASeq Wu" = wu_raw_data_matrix_input
  )
  centers = c( "SC RNASeq PREMIERE"
              ,"SN RNASeq UCSDWU"
              ,"SN RNASeq Muto"
              ,"SN RNASeq Wu"
  )

  pubMedIds = c("SC RNASeq PREMIERE" = premiere_pubmedId
                ,"SN RNASeq UCSDWU" = ucsd_pubmedId
                ,"SN RNASeq Muto" = muto_pubmedId
                ,"SN RNASeq Wu" = wu_pubmedId
  )
  keep_centers_for_each_analysis=centers
}#End

if (dataset=="SingleNucleus")
{#Begin
  annotations = list( "SN RNASeq UCSDWU" = ucsd_annotations
                      ,"SN RNASeq Muto" = muto_annotations
                      ,"SN RNASeq Wu" = wu_annotations
  )
  Raw_datas = list( "SN RNASeq UCSDWU" = UCSD_raw_data_matrix_input
                    ,"SN RNASeq Muto" = muto_raw_data_matrix_input
                    ,"SN RNASeq Wu" = wu_raw_data_matrix_input
  )
  centers = c( "SN RNASeq UCSDWU"
               ,"SN RNASeq Muto"
               ,"SN RNASeq Wu"
  )

  pubMedIds = c("SN RNASeq UCSDWU" = ucsd_pubmedId
                ,"SN RNASeq Muto" = muto_pubmedId
                ,"SN RNASeq Wu" = wu_pubmedId
  )
  keep_centers_for_each_analysis=centers
}#End

#######################################################################

{#Begin - Generate and write metadata spreadsheets from annotations
  {#Begin - Define metadata spreadsheet
    Col_names = c("Analysis","Dataset","Participant","Cell_count","PMID","Library_count")
    Col_length = length(Col_names)
    Row_names = 1
    Row_length = length(Row_names)
    metadata_base_line = array(NA,c(Row_length,Col_length),dimnames = list(Row_names,Col_names))
    metadata_base_line = as.data.frame(metadata_base_line)
    overall_metadata = c()
  }#End - Define metadata spreadsheet

  indexA=1
  for (indexA in 1:length(annotations))
  {#Begin
     current_annotation = annotations[[indexA]]
     current_dataset = names(annotations)[indexA]
     participants = unique(current_annotation$Participant)
     indexS=1
     for (indexS in 1:length(participants))
     {#Begin
        participant = participants[indexS]
        indexCurrentParticipant = which(current_annotation$Participant==participant)
        new_metadata_line = metadata_base_line
        new_metadata_line$Analysis = dataset
        new_metadata_line$Dataset = current_dataset
        new_metadata_line$Participant = participant
        new_metadata_line$Cell_count = length(current_annotation$Barcode[indexCurrentParticipant])
        new_metadata_line$Library_count = length(unique(current_annotation$Library[indexCurrentParticipant]))
        new_metadata_line$PMID = pubMedIds[[current_dataset]]
        if (length(overall_metadata)==0) { overall_metadata = new_metadata_line }
        else { overall_metadata = rbind(overall_metadata,new_metadata_line)}
     }#End
  }#End
  fileName = "Metadata_summary.txt"
  complete_fileName = paste(results_directory,fileName,sep='')
  write.table(overall_metadata,file=complete_fileName,sep='\t',quote=FALSE,row.names=FALSE,col.names = TRUE)
}#End - Generate and write metadata spreadsheets from annotations

#######################################################################

{#Begin - Add cell type and physio segment to annotations
  cellType_physioSegment = list("PT" = "ProxTub",
                                "DTL" = "Henle",
                                "ATL" = "Henle",
                                "TAL" = "Henle",
                                "DCT" = "DistTub",
                                "CNT" = "CD",
                                "PC" = "CD",
                                "IC" = "CD",
                                "B cell" = "noPhysio",
                                "EC" = "noPhysio",
                                "EPC" = "noPhysio",
                                "FIB" = "noPhysio",
                                "GC" = "noPhysio",
                                "IMM" = "noPhysio",
                                "INT" = "noPhysio",
                                "MAC" = "noPhysio",
                                "MC" = "noPhysio",
                                "MON" = "noPhysio",
                                "LEUK" = "noPhysio",
                                "NKC" = "noPhysio",
                                "PEC" = "noPhysio",
                                "POD" = "noPhysio",
                                "Unknown" = "noPhysio",
                                "T" = "noPhysio",
                                "vSMC/MC" = "noPhysio",
                                "vSMC/P" = "noPhysio"
                                )

  indexAnno=3
  for (indexAnno in 1:length(annotations))
  {#Begin
    current_annotations = annotations[[indexAnno]]
    current_annotations$Cell_type = "error"
    current_annotations$Physio_segment = "error"
    cell_subtypes = unique(current_annotations$Cluster)
    indexCellSubtype=1
    for (indexCellSubtype in 1:length(cell_subtypes))
    {#Begin
      cell_subtype = cell_subtypes[indexCellSubtype]
      indexCurrent_cell_subytpe = which(current_annotations$Cluster==cell_subtype)
      cell_type = strsplit(cell_subtype,"-")[[1]][1]
      if (cell_type=="") { current_annotations = c(); if (error_message=="") { error_message = paste("missing cell type: ",cell_type,s)}}
      current_annotations$Cell_type[indexCurrent_cell_subytpe] = cell_type
      physio_segment = cellType_physioSegment[[cell_type]]
      if (!is.null(physio_segment))
      { current_annotations$Physio_segment[indexCurrent_cell_subytpe] = physio_segment }
    }#End
    indexError = which(current_annotations$Cell_type=="error")
    if (length(indexError)) { current_annotations = c() }
    indexError = which(current_annotations$Physio_segment=="error")
    if (length(indexError)) { current_annotations = c() }
    annotations[[indexAnno]] = current_annotations
  }#End
}#End - Add cell type and physio segment to annotations

keep_only_nephron_cells = TRUE
if (keep_only_nephron_cells)
{#Begin - Keep only nephron cells in datasets (influences normalization before merging of datasets)
  indexAn = 2
  for (indexAn in 1:length(annotations))
  {#Begin
     current_annotation = annotations[[indexAn]]
     current_raw_data_matrix_input = Raw_datas[[indexAn]]
     keep_physio_segments = c("DistTub","CD","Henle","ProxTub")
     indexKeep = which(current_annotation$Physio_segment %in% keep_physio_segments)
     current_annotation = current_annotation[indexKeep,]
     barcodes = current_annotation$Barcode
     indexColKeep = which(colnames(current_raw_data_matrix_input) %in% barcodes)
     current_raw_data_matrix_input = current_raw_data_matrix_input[,indexColKeep]
     annotations[[indexAn]] = current_annotation
     Raw_datas[[indexAn]] = current_raw_data_matrix_input
  }#End
}#End - Keep only nephron cells in datasets (influences normalization before merging of datasets)

#######################################################################

{#Begin - Read and process libraries
  directory_scp_hierarchy = paste(main_directory,"MBCO_datasets//",sep='')
  hierarchy_fileNames = paste(library_of_interest,"_scpHierarchy.txt",sep='')
  complete_hierarchy_fileNames = paste(directory_scp_hierarchy,hierarchy_fileNames,sep='')
  scp_hierarchy = c()
  for (indexH in 1:length(complete_hierarchy_fileNames))
  {#Begin
     add_hierarchy = read.csv(file=complete_hierarchy_fileNames[indexH],header=TRUE,sep='\t')
     scp_hierarchy = rbind(scp_hierarchy,add_hierarchy)
  }#End

  indexParent = which(scp_hierarchy$Parent_scp=="Sodium large neutral amino acid symporter B2L")
  indexChild = which(scp_hierarchy$Child_scp=="Cationic amino acid vs large neutral amino acid sodium antiporter")
  indexPrioritarize = indexParent[indexParent %in% indexChild]
  scp_hierarchy$Parent_priotirization[indexPrioritarize] = 1
  
  indexParent = which(scp_hierarchy$Parent_scp=="Sodium potassium symporter")
  indexChild = which(scp_hierarchy$Child_scp=="Sodium potassium chloride symporter")
  indexPrioritarize = indexParent[indexParent %in% indexChild]
  scp_hierarchy$Parent_priotirization[indexPrioritarize] = 1
  
  {#Begin - Split parent_scps arrays and add as new lines
    new_hierarchy = c()
    for (indexH in 1:length(scp_hierarchy[,1]))
    {#Begin
       current_line = scp_hierarchy[indexH,]
       stringSplits = strsplit(current_line$Parent_scp,";")[[1]]
       for (indexS in 1:length(stringSplits))
       {#Begin
          new_line = current_line
          new_line$Parent_scp = stringSplits[indexS]
          if (length(new_hierarchy)==0) { new_hierarchy = new_line }
          else { new_hierarchy = rbind(new_hierarchy,new_line) }
       }#End
    }#End
    scp_hierarchy = new_hierarchy
  }#End - Split parent_scps arrays and add as new lines

  directory_libryaries = directory_scp_hierarchy
  library_fileNames = paste(library_of_interest,"_scpGeneAssociations.txt",sep='')
  combined_library = c()
  if (length(library_fileNames)>0)
  {#Begin
     for (indexLibrary in 1:length(library_fileNames))
     {#Begin
       complete_library_fileName = paste(directory_libryaries,library_fileNames[indexLibrary],sep='')
       add_library = read.csv(file=complete_library_fileName,stringsAsFactors = FALSE, sep='\t')
       if (length(combined_library)==0) { combined_library = add_library }
       else { combined_library = rbind(combined_library,add_library) }
     }#End
  }#End

  combined_library = combined_library[order(combined_library$Target_gene_symbol),]
  combined_library = combined_library[order(combined_library$Scp),]
  indexDuplicated = which(duplicated(paste(combined_library$Scp,combined_library$Target_gene_symbol,sep='')))
  if (length(indexDuplicated)>0) { combined_library = c();error_message = "duplicated combined libary"}
  combined_library$Gene_weight_factor = 1
}#End - Read and process libraries

if (length(combined_library)==0) { if (error_message=="") { error_message = "Read and process libraries"} }

add_to_fileNames = ""
only_selected_compartment_genes = TRUE
selected_compartment = "Plasma membrane"
if (only_selected_compartment_genes)
{#Begin - Remove all non-plasma membrane genes from libraries except for non-plasma membrane scps
   nonPlasmaMembraneParentSCPs = c("Energy")
   nonPlasmaMembraneSCPs = c()
   current_scps = nonPlasmaMembraneParentSCPs
   while (length(current_scps) > 0)
   {#Begin
      indexParent = which(scp_hierarchy$Parent_scp %in% current_scps)
      nonPlasmaMembraneSCPs = c(nonPlasmaMembraneSCPs,scp_hierarchy$Child_scp[indexParent])
      current_scps = scp_hierarchy$Child_scp[indexParent]
   }#End
   minium_confidence=4
   jensenlab_fileName = "human_compartment_integrated_full_2020February21.tsv"
   jensenlab_directory = paste(main_directory,"jensenLab_compartments//",sep='')
   jensenlab_completeFileName = paste(jensenlab_directory,jensenlab_fileName,sep='')
   jensenlab = read.csv(file=jensenlab_completeFileName,sep='\t',stringsAsFactors = FALSE)
   indexPM = which(jensenlab$Compartment==selected_compartment)
   indexAboveCutoff=which(jensenlab$Confidence>=minium_confidence)
   indexPMGenes = indexPM[indexPM %in% indexAboveCutoff]
   pm_genes = unique(jensenlab$Gene_symbol[indexPMGenes])
   indexKeep1 = which(combined_library$Target_gene_symbol %in% pm_genes)
   indexKeep2 = which(combined_library$Scp %in% nonPlasmaMembraneSCPs)
   indexKeep = unique(c(indexKeep1,indexKeep2))
   indexRemove = which(!combined_library$Target_gene_symbol %in% pm_genes)
   removed_library = combined_library[indexRemove,]
   combined_library =combined_library[indexKeep,]
   add_to_fileNames = paste(add_to_fileNames,"_",selected_compartment,sep='')
}#End - Remove all non-plasma membrane genes from libraries except for non-plasma membrane scps

if (length(combined_library)==0) { if (error_message=="") { error_message = "Remove all non-plasma membrane genes from libraries"} }

{#Begin - Plot # NCBI gene symbols per SCP
  ncbi_symbol_group_scps = list("Sodium"  = c("Sodium L2B","Sodium B2L"),
                                "Glucose" = c("Glucose L2B","Glucose B2L"))


  raw_data_names = names(Raw_datas)
  library = combined_library
  library_name = paste("only ",selected_compartment,sep='')
  Scp_genes_plots = list();
  scp_groups = names(ncbi_symbol_group_scps)
  for (indexScpGroup in 1:length(scp_groups))
  {#Begin
      scp_group = scp_groups[indexScpGroup]
      scps = ncbi_symbol_group_scps[[scp_group]]
      indexCurrent_scps = which(combined_library$Scp %in% scps)
      if (length(indexCurrent_scps)>0)
      {#Begin
         currentSCPs_library = combined_library[indexCurrent_scps,]
         table_currentSCPs_library = as.data.frame(table(currentSCPs_library$Scp))
         table_currentSCPs_library$Var1 = gsub(" ","\n",table_currentSCPs_library$Var1)
         Scp_genes_plot = ggplot(table_currentSCPs_library,aes(x=Var1,y=Freq))
         Scp_genes_plot = Scp_genes_plot + geom_bar(stat="identity")
         Scp_genes_plot = Scp_genes_plot + ylab(paste("# NCBI gene symbols",sep='')) + xlab("")
         Scp_genes_plot = Scp_genes_plot + ggtitle(paste(scp_group," - ",library_name,sep=''))
         Scp_genes_plot = Scp_genes_plot + theme(axis.text.x = element_text(angle=90,size=40,vjust=0.5,hjust=1))
         Scp_genes_plot = Scp_genes_plot + theme(axis.text.y = element_text(angle=0,size=40))
         Scp_genes_plot = Scp_genes_plot + theme(axis.title.y = element_text(angle=90,size=40))
         Scp_genes_plot = Scp_genes_plot + theme(plot.title = element_text(angle=0,size=40,hjust=0.5))
         Scp_genes_plots[[length(Scp_genes_plots)+1]] = Scp_genes_plot
      }#End
  }#End
  complete_png_fileName = paste(results_directory,"SCP_gene_symbols_count",add_to_fileNames,".png",sep='')
  cols_count = min(length(Scp_genes_plots),2)
  rows_count = ceiling(length(Scp_genes_plots)/cols_count)
  png(complete_png_fileName,width=2000*cols_count,height=1000*rows_count,res=100);
  do.call("grid.arrange",c(Scp_genes_plots,nrow=rows_count,ncol=cols_count))
  dev.off()
}#End - Plot # NCBI gene symbols per SCP

selected_SCPs_for_barplot=list()
if (library_of_interest=="NaAndGluTMTransport")
{#Begin - Define SCPs for nephron function in barplot L2B profiles
  selected_SCPs_for_barplot = list(   "Sodium transporter only" = c("Sodium L2B","Sodium physio")
                                     ,"Cell counts" = c("Cell counts","Sodium_noParaPT_claudin2 physio")
                                     ,"Sodium net L2B" = c("Sodium net L2B","Sodium physio")
                                     ,"Sodium_PTm37_henlem30 net L2B" = c("Sodium net L2B","Sodium_PTm37_henlem30 physio")
                                     ,"Glucose net L2B" = c("Glucose L2B","Glucose physio")
                                  )
}#End - Define SCPs for nephron function in barplot L2B profiles

if (length(combined_library)==0) { if (error_message=="") { error_message = "Define SCPs for nephron function in barplot L2B profiles"} }

selected_SCPs_for_scp_subgroups=list()#######################################################################
{#Begin - Define SCPs for sub scp distribution accross cell types - for transmembrane transport
  subscp_label = " sub "
  rm(selected_SCPs_for_scp_subgroups)
  if (library_of_interest=="NaAndGluTMTransport")
  {#Begin
    selected_SCPs_for_scp_subgroups = list(  "Sodium L2B M" = c( "Sodium L2B"
                                           #,"NKCC1 (SLC12A2)"
                                            ,"NKCC2 (SLC12A1)"
                                           #,"Sodium potassium chloride symporter"
                                           #,"Sodium potassium symporter"
                                           #,"Sodium chloride symporter"
                                            ,"NCC (SLC12A3)"
                                            ,"CNT1 (SLC28A1)"
                                            ,"OCTN2 (SLC22A5)"
                                            ,"CT1 (SLC6A8)"
                                           #,"Sodium chloride choline symporter"
                                           #,"Sodium chloride glycine symporter"
                                           #,"TauT"
                                           #,"Sodium chloride betaine symporter"
                                           #,"Sodium chloride creatine symporter"
                                           #,"Sodium potassium chloride symporter"
                                           #,"GAT1","GAT3"
                                           #,"Sodium chloride GABA symporter"
                                            ,"NBCn1 (SLC4A7)","EAAT1 (SLC1A3)","NPT4 (SLC17A3)"
                                            ,"PiT-2 (SLC20A2)"
                                            ,"PiT-1 (SLC20A1)"
                                            ,"Sodium chloride symporter"
                                           #,"Sodium chloride beta-alanine symporter"
                                           #,"Sodium chloride betaine symporter"
                                           #,"Sodium chloride cationic amino acid symporter"
                                           #,"Sodium chloride neutral amino acid symporter"
                                            ,"NBCe1 (SLC4A4)"
                                           #,"Sodium bicarbonate symporter"
                                           #,"Sodium phosphate symporter, Type II"
                                           #,"SLC17A2"
                                            ,"NPT4 (SLC17A3)","NPT1 (SLC17A1)"#,"SLC17A4"
                                            ,"NPT2a (SLC34A1)"#,"NPT2b","NPT2c","PiT-1","PiT-3"
                                            ,"Sodium phosphate symporter, Type II"
                                            ,"Sodium phosphate symporter"
                                            ,"NaS1 (SLC13A1)"
                                           #,"NDCBE (SLC4A8)"
                                           #,"Sodium sulfate symporter"
                                           #,"Sodium hydrogen symporter"
                                            ,"Sodium amino acid symporter"
                                           #,"SGLT1 (SLC5A1)"
                                            ,"SGLT2 (SLC5A2)"
                                           #,"Sodium sugar symporter"
                                            ,"SMCT2 (SLC5A12)"
                                           #,"NaDC1"
                                           #,"Sodium mono-, di- or tricarboxylate symporter"
                                           #,"Sodium carnitine, choline or betaine symporter"
                                           #,"Sodium small organic molecule symporter"
                                           #,"Sodium organic molecule symporter"
                                           #,"Sodium biologically active small organic molecule symporter"
                                           #,"Sodium biochemical cofactor symporter"
                                            ,"Sodium L2B by symporter"
                                           #,"NHE1","NHE2","NHE3","NHE4","NHE8"
                                            ,"NHE3 (SLC9A3)","NHE8 (SLC9A8)"
                                           #,"Sodium vs hydrogen antiporter"
                                           #,"Sodium hydrogen amino acid vs potassium antiporter"
                                            ,"EAAT3 (SLC1A1)"
                                           #,"Sodium vs potassium antiporter"
                                           #,"Sodium vs organic cation antiporter"
                                           #,"Sodium L2B by antiporter"
                                           #,"Sodium small organic molecule symporter"
                                           #,"Sodium creatine symporter"
                                            ,"Sodium L2B by antiporter"
                        ),
      "Sodium B2L M" = c(  "Sodium B2L"
                           ,"Sodium B2L by symporter"
                           ,"Sodium B2L by antiporter"
                           ,"NaDC3 (SLC13A3)"
                           ,"GAT2 (SLC6A13)"
                           ,"gamma+LAT1 (SLC7A7)"
                           ,"SNAT1 (SLC38A1)"
                           #,"SNAT2 (SLC38A2)"
                           ,"SNAT3 (SLC38A3)"
                           #,"SNAT4 (SLC38A4)"
                           #,"SNAT5 (SLC38A5)"
                           #,"SNAT7 (SLC38A7)"
                           #,"Sodium B2L by symporter" #SCP is missing
                           #,"Sodium B2L by unimolecular transport"
                           ,"NCX1 (SLC8A1)"
                           #,"NCX2 (SLC8A2)"
                           #,"NCX3 (SLC8A3)"
                           #,"NCKX3 (SLC24A3)"
                           #,"Calcium vs sodium antiporter"
                           #,"Calcium potassium vs sodium antiporter"
                        ),
      "Glucose L2B M" = c(   "Glucose L2B"
                             ,"SGLT1 (SLC5A1)"
                             ,"SGLT2 (SLC5A2)"
                             #,"SGLT5 (SLC5A10)"
                             ,"Sodium glucose symporter"
      )
    )
  }#End

  scp_subgroups = names(selected_SCPs_for_scp_subgroups)
  missing_scps = c()

  all_add_combined_library = c()
  error_messages = c()
  missing_scps = c()
  indexScpSubGroup=2
  for (indexScpSubGroup in 1:length(scp_subgroups))
  {#Begin - Add scps foreach scp sub group
     current_hierarchy = scp_hierarchy
     scp_subgroup = scp_subgroups[indexScpSubGroup]
     final_parent_scp = gsub(" M","",scp_subgroup)
     full_subscp_label = paste(subscp_label,scp_subgroup,sep='')
     scps = unique(selected_SCPs_for_scp_subgroups[[scp_subgroup]])
     scp_ancestors_list = list()
     indexScp=7
     for (indexScp in 1:length(scps))
     {#Begin - Fill scp_ancestors_list and check for missing scps
        current_scp = scps[indexScp]
        current_ancestor_scps = c()
        if (!current_scp %in% scp_hierarchy$Child_scp) { missing_scps = c(missing_scps,current_scp) }
        else
        {#Begin - Fill scp_ancestors_list
           search_scps = current_scp
           while (length(search_scps)>0)
           {#Begin
              indexParentScp = which(tolower(scp_hierarchy$Child_scp) %in% tolower(search_scps))
              search_scp_strings = scp_hierarchy$Parent_scp[indexParentScp]
              search_scps = c()
              if (length(search_scp_strings)>0)
              {#Begin
                 for (indexSearch in 1:length(search_scp_strings))
                 {#Begin
                    search_scps = c(search_scps,strsplit(search_scp_strings[indexSearch],";")[[1]])
                 }#End
                 current_ancestor_scps = c(current_ancestor_scps,search_scps)
              }#End
           }#End
           if (  (!final_parent_scp %in% current_ancestor_scps)
                &(final_parent_scp!=current_scp) )
           { error_messages = c(error_messages,paste(current_scp," is no child of ",final_parent_scp,sep='')) }
           scp_ancestors_list[[current_scp]] = current_ancestor_scps
        }#End - Fill scp_ancestors_list
     }#End - Fill scp_ancestors_list and check for missing scps

     indexCurrent_scps = which(combined_library$Scp %in% scps)
     add_combined_library = combined_library[indexCurrent_scps,]
     add_unique_genes_length = length(unique(add_combined_library$Target_gene_symbol))

     indexScp = 1
     for (indexScp in 1:length(scps))
     {#Begin - Remove all scp genes from offspring scps
        current_scp = scps[indexScp]
        ancestor_scps = scp_ancestors_list[[current_scp]]
        indexCurrentAncestorScps = which(add_combined_library$Scp %in% ancestor_scps)
        indexCurrentScp = which(add_combined_library$Scp == current_scp)
        current_scp_genes = add_combined_library$Target_gene_symbol[indexCurrentScp]
        indexCurrentGenes = which(add_combined_library$Target_gene_symbol %in% current_scp_genes)
        indexCurrentScpGenes_in_ancestor = indexCurrentAncestorScps[indexCurrentAncestorScps %in% indexCurrentGenes]
        indexKeep = 1:length(add_combined_library[,1])
        indexKeep = indexKeep[!indexKeep %in% indexCurrentScpGenes_in_ancestor]
        add_combined_library = add_combined_library[indexKeep,]
     }#End - Remove all scp genes from offspring scps

     indexDuplicated = which(duplicated(add_combined_library$Target_gene_symbol))
     if (length(indexDuplicated)>0)
     {#Begin
       duplicated_genes = unique(add_combined_library$Target_gene_symbol[indexDuplicated])
       for (indexDuplicated in 1:length(duplicated_genes))
       {#Begin
          duplicated_gene = duplicated_genes[indexDuplicated]
          indexCurrentDuplicated = which(add_combined_library$Target_gene_symbol==duplicated_gene)
          add_combined_library$Gene_weight_factor[indexCurrentDuplicated] = 1 / length(indexCurrentDuplicated)
       }#End
     }#End

     if (sum(add_combined_library$Gene_weight_factor) != add_unique_genes_length)
     { error_messages = c(error_messages,"Remove all scp genes from offspring scps") }

     add_combined_library$Scp = paste(add_combined_library$Scp,full_subscp_label,sep='')
     scps = paste(scps,full_subscp_label,sep='')
     scps = scps[!scps %in% missing_scps]

     #add_combined_library$Scp = gsub("L2B by symporter","other symporter",add_combined_library$Scp)
     #add_combined_library$Scp = gsub("L2B by antiporter","other antiporter",add_combined_library$Scp)
     #scps = gsub("L2B by symporter","other symporter",scps)
     #scps = gsub("L2B by antiporter","other antiporter",scps)

     selected_SCPs_for_scp_subgroups[[scp_subgroup]] = scps

     if (length(all_add_combined_library)==0) { all_add_combined_library = add_combined_library }
     else { all_add_combined_library = rbind(all_add_combined_library,add_combined_library)}
  }#End - Add scps foreach scp sub group

  combined_library = rbind(combined_library,all_add_combined_library)

  if (length(combined_library)==0) { if (error_message=="") { error_message = "Define SCPs for SCP cross-comparison"} }
  if (length(error_messages)>0) { if (error_message=="") { error_message = "Check error messages"; combined_library = c() }}
  if (length(missing_scps)>0) { if (error_message=="") { error_message = "Check missing scps"; combined_library = c() }}

}#End - Define SCPs for sub scp distribution accross cell types - for transmembrane transport ####here

if (length(combined_library)==0) { if (error_message=="") { error_message = "Define SCPs for x=Cell subtype, group=SCP and add SCPs to combined library"} }

reserveScps_segment = list("Sodium amino acid symporter" = c("Henle","DistTub","CD"))

identify_scps_with_high_read_counts = FALSE
if (identify_scps_with_high_read_counts)  { selected_SCPs = unique(combined_library$Scp) }
scp_subgroup_scps = gsub(" M","",names(selected_SCPs_for_scp_subgroups))
if (!identify_scps_with_high_read_counts) { selected_SCPs = unlist(c(selected_SCPs_for_scp_subgroups,selected_SCPs_for_barplot,scp_subgroup_scps)) }
selected_SCPs = unique(c(selected_SCPs,names(reserveScps_segment)))

{#Begin - Keep only selected SCPs in library and check for duplicates within each selected SCP
   combined_library2 = combined_library
   combined_library = combined_library2
   unique_selected_SCPs = unique(c(unlist(selected_SCPs)))
   unique_selected_SCPs = unique_selected_SCPs[order(unique_selected_SCPs)]
   indexKeep = which(combined_library$Scp %in% unique_selected_SCPs)
   combined_library = combined_library[indexKeep,]
   for (indexSelected in 1:length(unique_selected_SCPs))
   {#Begin
      selected_scp = unique_selected_SCPs[indexSelected]
      indexCurrentSCP = which(combined_library$Scp==selected_scp)
      if (length(combined_library$Target_gene_symbol[indexCurrentSCP])!=length(unique(combined_library$Target_gene_symbol[indexCurrentSCP]))) { combined_library = c();if (error_message=="") { error_message = paste("Keep only selected SCPs in library and check for duplicates within each selected SCP-index: ",indexSelected,": ",selected_scp,sep='') };break }
   }#End
}#End - Keep only selected SCPs in library and check for duplicates within each selected SCP

if (length(combined_library)==0) { if (error_message=="") { error_message = "Keep only selected SCPs in library and check for duplicates within each selected SCP"} }

{#Begin - Define read_counts_cell_subtypes
  Col_names = c("Cell_type","Cell_type_factor","Cell_subtype","Physio_segment","Cell_subtype_counts","Center","SCP","Read_counts","Average_read_counts","Data_type")
  Col_length = length(Col_names)
  Row_names = 1;
  Row_length = length(Row_names)
  read_count_cell_subtype_base_line = array(NA,c(Row_length,Col_length),dimnames=list(Row_names,Col_names))
  read_count_cell_subtype_base_line = as.data.frame(read_count_cell_subtype_base_line)
  read_counts_cell_subtypes = c()
}#End - Define read_counts_cell_subtypes

if (algorithm_of_interest == "Total_sum")
{#Begin - Calculate read counts per cell subtype and scp
  scps = unique(c(unlist(selected_SCPs)))
  scps = scps[order(scps)]
  indexPhysio = grep("physio",scps)
  indexKeep = 1:length(scps)
  indexKeep = indexKeep[!indexKeep %in% indexPhysio]
  scps = scps[indexKeep]
  indexScp=2
  scp_progress = progress_bar$new(total=length(scps))
  for (indexScp in 1:length(scps))
  {#Begin
    scp = scps[indexScp]
    scp_progress$tick()
    indexCurrentSCP = which(combined_library$Scp==scp)
    if (length(indexCurrentSCP)>0)
    {#Begin - if (length(indexCurrentSCP)>0)
    current_genes = combined_library$Target_gene_symbol[indexCurrentSCP]
    weight_factors = combined_library$Gene_weight_factor[indexCurrentSCP]
    names(weight_factors) = current_genes
    annotation_centers = names(annotations)
    indexCenter = 1
    for (indexCenter in 1:length(annotation_centers))
    {#Begin
      center = annotation_centers[indexCenter]
      current_annotation = annotations[[center]]
      current_rawData = Raw_datas[[center]]
      cellSubTypes = unique(current_annotation$Cluster)
      indexCellSubType = 1
      for (indexCellSubType in 1:length(cellSubTypes))
      {#Begin
        cellSubType = cellSubTypes[indexCellSubType]
        indexCurrentCellSubtype = which(current_annotation$Cluster==cellSubType)
        current_barcodes = current_annotation$Barcode[indexCurrentCellSubtype]
        current_cellTypes = unique(current_annotation$Cell_type[indexCurrentCellSubtype])
        if (length(current_cellTypes)==1) { cellType = current_cellTypes }
        else { cellType = paste(current_cellTypes,collapse=";") }
        current_physioSegments = unique(current_annotation$Physio_segment[indexCurrentCellSubtype])
        if (length(current_physioSegments)==1) { physio_segment = current_physioSegments }
        else { physio_segment = paste(current_physioSegments,collapse=";") }
        indexCurrentCols = which(colnames(current_rawData)%in%current_barcodes)

        unique_weight_factors = unique(weight_factors)
        total_read_counts_sum = 0
        indexWeightFactor = 1
        for (indexWeightFactor in 1:length(unique_weight_factors))
        {#Begin
           weight_factor = unique_weight_factors[indexWeightFactor]
           indexCurrentWeightFactor = which(weight_factors==weight_factor)
           current_weight_factors = weight_factors[indexCurrentWeightFactor]
           indexCurrentRows = which(rownames(current_rawData)%in%names(current_weight_factors))
           if (length(indexCurrentRows)>0)
           { total_read_counts_sum = total_read_counts_sum + sum(current_rawData[indexCurrentRows,indexCurrentCols]) * weight_factor }
        }#End

        indexCurrentRowGenes = which(rownames(current_rawData)[indexCurrentRows] %in% current_genes)

        current_read_count_cell_subtype_base_line = read_count_cell_subtype_base_line
        current_read_count_cell_subtype_base_line$Center = center
        current_read_count_cell_subtype_base_line$Cell_type = cellType
        current_read_count_cell_subtype_base_line$Cell_subtype = cellSubType
        current_read_count_cell_subtype_base_line$Physio_segment = physio_segment
        current_read_count_cell_subtype_base_line$Data_type = "1_Read counts"
        current_read_count_cell_subtype_base_line$Cell_subtype_counts = length(indexCurrentCols)
        current_read_count_cell_subtype_base_line$SCP = scp
        if ((length(indexCurrentCols)>0)&(total_read_counts_sum>0))
        {#Begin
          current_read_count_cell_subtype_base_line$Read_counts = total_read_counts_sum
          current_read_count_cell_subtype_base_line$Average_read_counts = current_read_count_cell_subtype_base_line$Read_counts/current_read_count_cell_subtype_base_line$Cell_subtype_counts
        }#End
        else
        {#Begin
          current_read_count_cell_subtype_base_line$Read_counts = 0
          current_read_count_cell_subtype_base_line$Average_read_counts = 0
        }#End

        read_count_cell_subtype_line = current_read_count_cell_subtype_base_line
        if (length(read_counts_cell_subtypes)==0) { read_counts_cell_subtypes = read_count_cell_subtype_line }
        else { read_counts_cell_subtypes = rbind(read_counts_cell_subtypes, read_count_cell_subtype_line) }
      }#End
    }#End
    }#End - if (length(indexCurrentSCP)>0)
  }#End
  scp_progress$terminate()
}#End - Calculate read counts per cell subtype and scp

if (identify_scps_with_high_read_counts)
{#Begin - Screen for scps with high read counts
  Get_all_children_scps = function(search_term)
  {#Begin
    indexAminoAcid = grep(search_term,scp_hierarchy$Parent_scp)
    scps = scp_hierarchy$Parent_scp[indexAminoAcid]
    child_scps = unique(scp_hierarchy$Child_scp[indexAminoAcid])
    scps = unique(c(scps,child_scps))
    while (length(child_scps)>0)
    {#Begin
      parent_scps = child_scps
      indexParentSCPs = which(scp_hierarchy$Parent_scp %in% parent_scps)
      child_scps = scp_hierarchy$Child_scp[indexParentSCPs]
      scps = unique(c(scps,child_scps))
    }#End
    return (scps)
  }#End


  allSCPs_read_counts = read_counts_cell_subtypes
  indexSubScps = grep(subscp_label,allSCPs_read_counts$SCP)
  indexKeep = 1:length(allSCPs_read_counts[,1])
  indexKeep = indexKeep[!indexKeep %in% indexSubScps]
  allSCPs_read_counts = allSCPs_read_counts[indexKeep,]
  allSCPs_read_counts = allSCPs_read_counts[order(allSCPs_read_counts$Read_counts,decreasing=TRUE),]

  indexAntiporter = grep("antiporter",allSCPs_read_counts$SCP)
  indexSymporter = grep("ymporter",allSCPs_read_counts$SCP)
  indexChannel = grep("hannel",allSCPs_read_counts$SCP)
  indexB2L = grep("B2L",allSCPs_read_counts$SCP)
  indexL2B = grep("L2B",allSCPs_read_counts$SCP)

  sodium_scps = Get_all_children_scps("odium")
  indexSodium = which(allSCPs_read_counts$SCP %in% sodium_scps)
  sodium_read_counts = allSCPs_read_counts[indexSodium,]
  sodium_read_counts = sodium_read_counts[order(sodium_read_counts$Read_counts,decreasing=TRUE),]
  sodium_scps = unique(sodium_read_counts$SCP)

  indexSodium_B2L = indexSodium[indexSodium %in% indexB2L]
  sodiumB2L_read_counts = allSCPs_read_counts[indexSodium_B2L,]
  sodiumB2L_read_counts = sodiumB2L_read_counts[order(sodiumB2L_read_counts$Read_counts,decreasing=TRUE),]
  sodium_B2L_scps = unique(sodiumB2L_read_counts$SCP)

  indexSodiumAntiporter = indexSodium[indexSodium %in% indexAntiporter]
  sodiumAntiporter_read_counts = allSCPs_read_counts[indexSodiumAntiporter,]
  sodiumAntiporter_read_counts = sodiumAntiporter_read_counts[order(sodiumAntiporter_read_counts$Read_counts,decreasing=TRUE),]
  sodium_antiporter_scps = unique(sodiumAntiporter_read_counts$SCP)

  indexSodiumSymporter = indexSymporter[indexSymporter %in% indexSodium]
  sodiumSymporter_read_counts = allSCPs_read_counts[indexSodiumSymporter,]
  sodiumSymporter_read_counts = sodiumSymporter_read_counts[order(sodiumSymporter_read_counts$Read_counts,decreasing=TRUE),]
  sodium_symporter_scps = unique(sodiumSymporter_read_counts$SCP)

  indexSodiumChannel = indexSodium[indexSodium %in% indexChannel]
  sodiumChannel_read_counts = allSCPs_read_counts[indexSodiumChannel,]
  sodiumChannel_read_counts = sodiumChannel_read_counts[order(sodiumChannel_read_counts$Read_counts,decreasing=TRUE),]
  sodium_channel_scps = unique(sodiumChannel_read_counts$SCP)

  hydrogen_scps = Get_all_children_scps("ydrogen")
  indexHydrogen = which(allSCPs_read_counts$SCP %in% hydrogen_scps)
  indexHydrogenAntiporter = indexHydrogen[indexHydrogen %in% indexAntiporter]
  hydrogenAntiporter_read_counts = allSCPs_read_counts[indexHydrogenAntiporter,]
  hydrogenAntiporter_read_counts = hydrogenAntiporter_read_counts[order(hydrogenAntiporter_read_counts$Read_counts,decreasing=TRUE),]
  hydrogen_antiporter_scps = unique(hydrogenAntiporter_read_counts$SCP)

  indexHydrogenSymporter = indexSymporter[indexSymporter %in% indexHydrogen]
  hydrogenSymporter_read_counts = allSCPs_read_counts[indexHydrogenSymporter,]
  hydrogenSymporter_read_counts = hydrogenSymporter_read_counts[order(hydrogenSymporter_read_counts$Read_counts,decreasing=TRUE),]
  hydrogen_symporter_scps = unique(hydrogenSymporter_read_counts$SCP)

  chloride_scps = Get_all_children_scps("hloride")
  indexChloride = which(allSCPs_read_counts$SCP %in% chloride_scps)
  indexChlorideAntiporter = indexChloride[indexChloride %in% indexAntiporter]
  chlorideAntiporter_read_counts = allSCPs_read_counts[indexChlorideAntiporter,]
  chlorideAntiporter_read_counts = chlorideAntiporter_read_counts[order(chlorideAntiporter_read_counts$Read_counts,decreasing=TRUE),]
  chloride_antiporter_scps = unique(chlorideAntiporter_read_counts$SCP)

  indexChlorideSymporter = indexSymporter[indexSymporter %in% indexChloride]
  chlorideSymporter_read_counts = allSCPs_read_counts[indexChlorideSymporter,]
  chlorideSymporter_read_counts = chlorideSymporter_read_counts[order(chlorideSymporter_read_counts$Read_counts,decreasing=TRUE),]
  chloride_symporter_scps = unique(chlorideSymporter_read_counts$SCP)

  indexChlorideChannel = indexChannel[indexChannel %in% indexChloride]
  chlorideChannel_read_counts = allSCPs_read_counts[indexChlorideChannel,]
  chlorideChannel_read_counts = chlorideChannel_read_counts[order(chlorideChannel_read_counts$Read_counts,decreasing=TRUE),]
  chloride_channel_scps = unique(chlorideChannel_read_counts$SCP)

  potassium_scps = Get_all_children_scps("otassium")
  indexPotassium = which(allSCPs_read_counts$SCP %in% potassium_scps)
  indexPotassiumAntiporter = indexPotassium[indexPotassium %in% indexAntiporter]
  potassiumAntiporter_read_counts = allSCPs_read_counts[indexPotassiumAntiporter,]
  potassiumAntiporter_read_counts = potassiumAntiporter_read_counts[order(potassiumAntiporter_read_counts$Read_counts,decreasing=TRUE),]
  potassium_antiporter_scps = unique(potassiumAntiporter_read_counts$SCP)

  indexPotassiumChannel = indexPotassium[indexPotassium %in% indexChannel]
  potassiumChannel_read_counts = allSCPs_read_counts[indexPotassiumChannel,]
  potassiumChannel_read_counts = potassiumChannel_read_counts[order(potassiumChannel_read_counts$Read_counts,decreasing=TRUE),]
  potassium_channelr_scps = unique(potassiumChannel_read_counts$SCP)

  bicarbonate_scps = Get_all_children_scps("icarbonate")
  indexBicarbonate = which(allSCPs_read_counts$SCP %in% bicarbonate_scps)
  indexBicarbonateSymporter = indexBicarbonate[indexBicarbonate %in% indexSymporter]
  bicarbonateSymporter_read_counts = allSCPs_read_counts[indexBicarbonateSymporter,]
  bicarbonateSymporter_read_counts = bicarbonateSymporter_read_counts[order(bicarbonateSymporter_read_counts$Read_counts,decreasing=TRUE),]
  bicarbonate_symporter_scps = unique(bicarbonateSymporter_read_counts$SCP)

  indexBicarbonateAntiporter = indexBicarbonate[indexBicarbonate %in% indexAntiporter]
  bicarbonateAntiporter_read_counts = allSCPs_read_counts[indexBicarbonateAntiporter,]
  bicarbonateAntiporter_read_counts = bicarbonateAntiporter_read_counts[order(bicarbonateAntiporter_read_counts$Read_counts,decreasing=TRUE),]
  bicarbonate_antiporter_scps = unique(bicarbonateAntiporter_read_counts$SCP)

  phosphate_scps = Get_all_children_scps("hosphate")
  indexPhosphate = which(allSCPs_read_counts$SCP %in% phosphate_scps)
  indexPhosphateAntiporter = indexPhosphate[indexPhosphate %in% indexAntiporter]
  phosphateAntiporter_read_counts = allSCPs_read_counts[indexPhosphateAntiporter,]
  phosphateAntiporter_read_counts = phosphateAntiporter_read_counts[order(phosphateAntiporter_read_counts$Read_counts,decreasing=TRUE),]
  phosphate_antiporter_scps = unique(phosphateAntiporter_read_counts$SCP)

  indexPhosphateSymporter = indexPhosphate[indexPhosphate %in% indexSymporter]
  phosphateSymporter_read_counts = allSCPs_read_counts[indexPhosphateSymporter,]
  phosphateSymporter_read_counts = phosphateSymporter_read_counts[order(phosphateSymporter_read_counts$Read_counts,decreasing=TRUE),]
  phosphate_symporter_scps = unique(phosphateSymporter_read_counts$SCP)

  glucose_scps = Get_all_children_scps("lucose")
  indexGlucose = which(allSCPs_read_counts$SCP %in% glucose_scps)
  indexGlucoseSymporter = indexGlucose[indexGlucose %in% indexSymporter]
  glucoseSymporter_read_counts = allSCPs_read_counts[indexGlucoseSymporter,]
  glucoseSymporter_read_counts = glucoseSymporter_read_counts[order(glucoseSymporter_read_counts$Read_counts,decreasing=TRUE),]
  glucose_symporter_scps = unique(glucoseSymporter_read_counts$SCP)

  calcium_scps = Get_all_children_scps("alcium")
  indexCalcium = which(allSCPs_read_counts$SCP %in% calcium_scps)
  indexCalcium_channel = indexCalcium[indexCalcium %in% indexChannel]
  calciumChannel_read_counts = allSCPs_read_counts[indexCalcium_channel,]
  calciumChannel_read_counts = calciumChannel_read_counts[order(calciumChannel_read_counts$Read_counts,decreasing=TRUE),]
  calcium_channel_scps = unique(calciumChannel_read_counts$SCP)

  aa_scps = Get_all_children_scps("mino acid")
  indexAA = which(allSCPs_read_counts$SCP %in% aa_scps)
  indexAA_na = indexAA[indexAA %in% c(indexSodiumSymporter,indexSodiumAntiporter)]
  aa_read_counts = allSCPs_read_counts[indexAA_na,]
  aa_read_counts = aa_read_counts[order(aa_read_counts$Read_counts,decreasing = FALSE),]
  aa_scps = unique(aa_read_counts$SCP)

}#End - Screen for scps with high read counts

if (library_of_interest=="NaAndGluTMTransport")
{#Begin - Generate net L2B as a duplicate of L2B
  net_L2B_label = " net L2B"
  indexL2B = grep("L2B",read_counts_cell_subtypes$SCP)
  indexSubScp = grep(subscp_label,read_counts_cell_subtypes$SCP)
  indexL2B = indexL2B[!indexL2B %in% indexSubScp]
  net_read_counts = read_counts_cell_subtypes[indexL2B,]
  net_read_counts$SCP = gsub(" L2B",net_L2B_label,net_read_counts$SCP)
  read_counts_cell_subtypes = rbind(read_counts_cell_subtypes,net_read_counts)
}#End - Generate net L2B as a duplicate of L2B

remove_B2L_scps_from_net_L2B=TRUE
if ((remove_B2L_scps_from_net_L2B)&(library_of_interest=="NaAndGluTMTransport"))
{#Begin - Remove B2L scps from net L2B
  indexL2B = grep("net L2B",read_counts_cell_subtypes$SCP)
  indexSub = grep(subscp_label,read_counts_cell_subtypes$SCP)
  indexL2B = indexL2B[!indexL2B %in% indexSub]
  L2B_scps = unique(read_counts_cell_subtypes$SCP[indexL2B])
  read_counts_cell_subtypes = read_counts_cell_subtypes[order(read_counts_cell_subtypes$SCP),]
  indexReab = 1
  for (indexReab in 1:length(L2B_scps))
  {#Begin
     L2B_scp = L2B_scps[indexReab]
     splitStrings = strsplit(L2B_scp," ")[[1]]
     current_ion = splitStrings[1]
     if (length(splitStrings)-2>1)
     {#Begin
        for (indexSplit in 2:(length(splitStrings)-2))
        {#Begin
           current_ion = paste(current_ion," ",splitStrings[indexSplit],sep='')
        }#End
     }#End
     B2L_scp = paste(current_ion," B2L",sep='')
     indexCurrentB2L = which(read_counts_cell_subtypes$SCP==B2L_scp)
     if (length(indexCurrentB2L)>0)
     {#Begin
        indexCurrentL2B = which(read_counts_cell_subtypes$SCP==L2B_scp)
        if (length(indexCurrentL2B)==length(indexCurrentB2L))
        {#Begin
           if (length(which(read_counts_cell_subtypes$Cell_subtype[indexCurrentL2B]!=read_counts_cell_subtypes$Cell_subtype[indexCurrentB2L]))!=0)
           {#Begin
             read_counts_cell_subtypes = c(); if (error_message=="") { error_message = paste("Remove B2L scps from net L2B: ",L2B_scp,sep='') }
           }#End
           read_counts_cell_subtypes$Read_counts[indexCurrentL2B] = read_counts_cell_subtypes$Read_counts[indexCurrentL2B] - read_counts_cell_subtypes$Read_counts[indexCurrentB2L]
        }#End
        if (length(indexCurrentL2B)!=length(indexCurrentB2L))
        {#Begin
           read_counts_cell_subtypes = c(); if (error_message=="") { error_message = paste("Remove B2L scps from net L2B: ",L2B_scp,sep='') }
        }#End
     }#End
  }#End
  add_to_fileNames = paste(add_to_fileNames,"_B2LRemoved",sep='')
}#End - Remove B2L scps from net L2B

if ((remove_B2L_scps_from_net_L2B)&(library_of_interest=="NaAndGluTMTransport"))
{#Begin - Generate new SCPs with cell counts and sodium net L2B divided by cell counts as read counts
  indexSodiumNetL2B = which(read_counts_cell_subtypes$SCP=="Sodium net L2B")
  new_read_counts = read_counts_cell_subtypes[indexSodiumNetL2B,]
  new_read_counts$SCP = "Cell counts"
  new_read_counts$Read_counts = new_read_counts$Cell_subtype_counts
  new_read_counts$Average_read_counts = 1

  new_read_counts2 = read_counts_cell_subtypes[indexSodiumNetL2B,]
  new_read_counts2$SCP = "Sodium net L2B per cell counts"
  new_read_counts2$Read_counts = new_read_counts$Read_counts / new_read_counts$Cell_subtype_counts
  new_read_counts2$Average_read_counts = 1

  read_counts_cell_subtypes = rbind(read_counts_cell_subtypes,new_read_counts)
  read_counts_cell_subtypes = rbind(read_counts_cell_subtypes,new_read_counts2)
}#End - Generate new SCPs with cell counts and sodium net L2B divided by cell counts as read counts

readCounts_fileName = "Read_counts.txt"
complete_readCounts_fileName = paste(results_directory,readCounts_fileName,sep='')
write.table(read_counts_cell_subtypes,file=complete_readCounts_fileName,sep='\t',quote = FALSE)

{#Begin - Define cell subtype colors
  podGlo_color = "aquamarine4";
  PTTI_color = "orange2"
  PT4_color = "orange2"
  PT5_color = "orange2"
  cnt_pc_color = "indianred2"
  cortex_PCCD_color = "darkorchid1";
  medulla_PCCD_color = "darkorchid4";
  cortex_henle_color = "dodgerblue1";
  medulla_henle_color = "dodgerblue4";
  cellType_other_color = "black"
  interstitium_color = "gray40"
  endothelial_color = "cyan3"
  immune_color = "yellow"
  ureter_color = "white"

  center_cell_subtype_colors_list = list()
  premiere_cell_subtype_colors = list()
  premiere_cell_subtype_colors[["POD"]] = podGlo_color
  premiere_cell_subtype_colors[["vSMC_MC"]] = podGlo_color
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
  premiere_cell_subtype_colors[["CNT"]] = cortex_PCCD_color
  premiere_cell_subtype_colors[["CNT_PC"]] = cortex_PCCD_color
  premiere_cell_subtype_colors[["PC"]] = cortex_PCCD_color
  premiere_cell_subtype_colors[["tPC_IC"]] = cortex_PCCD_color
  premiere_cell_subtype_colors[["IC-A"]] = cortex_PCCD_color
  premiere_cell_subtype_colors[["IC-B"]] = cortex_PCCD_color
  premiere_cell_subtype_colors[["IC"]] = cortex_PCCD_color
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
  premiere_cell_subtype_colors[["Ureter"]] = ureter_color
  center_cell_subtype_colors_list[["SC RNASeq PREMIERE"]] = premiere_cell_subtype_colors

  ucsd_cell_subtype_colors = list()
  ucsd_cell_subtype_colors[["MC"]] = podGlo_color
  ucsd_cell_subtype_colors[["POD"]] = podGlo_color
  ucsd_cell_subtype_colors[["EPC"]] = podGlo_color
  ucsd_cell_subtype_colors[["vSMC_P"]] = podGlo_color
  ucsd_cell_subtype_colors[["PT-1"]] = PTTI_color
  ucsd_cell_subtype_colors[["PT-2"]] = PTTI_color
  ucsd_cell_subtype_colors[["PT-3"]] = PTTI_color
  ucsd_cell_subtype_colors[["PT-4"]] = PTTI_color
  ucsd_cell_subtype_colors[["PT-5"]] = PTTI_color
  ucsd_cell_subtype_colors[["Unk"]] = PTTI_color
  ucsd_cell_subtype_colors[["DTL"]] = medulla_henle_color
  ucsd_cell_subtype_colors[["DL"]] = medulla_henle_color
  ucsd_cell_subtype_colors[["ATL-1"]] = medulla_henle_color
  ucsd_cell_subtype_colors[["ATL-2"]] = medulla_henle_color
  ucsd_cell_subtype_colors[["ATL-3"]] = medulla_henle_color
  ucsd_cell_subtype_colors[["TAL-1"]] = medulla_henle_color
  ucsd_cell_subtype_colors[["TAL-2"]] = cortex_henle_color
  ucsd_cell_subtype_colors[["DCT"]] = PTTI_color
  ucsd_cell_subtype_colors[["CNT"]] = cortex_PCCD_color
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
  ucsd_cell_subtype_colors[["IMM"]] = immune_color
  ucsd_cell_subtype_colors[["Ureter"]] = ureter_color
  center_cell_subtype_colors_list[["SN RNASeq UCSDWU"]] = ucsd_cell_subtype_colors

  muto_cell_subtype_colors = list()
  muto_cell_subtype_colors[["POD"]] = podGlo_color
  muto_cell_subtype_colors[["PT"]] = PTTI_color
  muto_cell_subtype_colors[["PT-VCAM1"]] = PTTI_color
  muto_cell_subtype_colors[["DTL"]] = cortex_henle_color
  muto_cell_subtype_colors[["TAL"]] = cortex_henle_color
  muto_cell_subtype_colors[["DCT-1"]] = PTTI_color
  muto_cell_subtype_colors[["DCT-2"]] = PTTI_color
  muto_cell_subtype_colors[["CNT"]] = cortex_PCCD_color
  muto_cell_subtype_colors[["PC"]] = cortex_PCCD_color
  muto_cell_subtype_colors[["IC-A"]] = cortex_PCCD_color
  muto_cell_subtype_colors[["IC-B"]] = cortex_PCCD_color
  muto_cell_subtype_colors[["Ureter"]] = ureter_color
  center_cell_subtype_colors_list[["SN RNASeq Muto"]] = muto_cell_subtype_colors

  wu_cell_subtype_colors = list()
  wu_cell_subtype_colors[["POD"]] = podGlo_color
  wu_cell_subtype_colors[["PT"]] = PTTI_color
  wu_cell_subtype_colors[["DTL"]] = cortex_henle_color
  wu_cell_subtype_colors[["ATL"]] = cortex_henle_color
  wu_cell_subtype_colors[["TAL"]] = cortex_henle_color
  wu_cell_subtype_colors[["DCT"]] = PTTI_color
  wu_cell_subtype_colors[["CNT"]] = cortex_PCCD_color
  wu_cell_subtype_colors[["PC"]] = cortex_PCCD_color
  wu_cell_subtype_colors[["IC"]] = cortex_PCCD_color
  wu_cell_subtype_colors[["Ureter"]] = ureter_color
  center_cell_subtype_colors_list[["SN RNASeq Wu"]] = wu_cell_subtype_colors

  merge_premiere_subtype_colors = premiere_cell_subtype_colors[names(premiere_cell_subtype_colors)!="Ureter"]
  merge_ucsd_subtype_colors = ucsd_cell_subtype_colors[names(ucsd_cell_subtype_colors)!="Ureter"]
  names(merge_ucsd_subtype_colors) = paste(names(merge_ucsd_subtype_colors),"-U",sep='')
  names(merge_premiere_subtype_colors) = paste(names(merge_premiere_subtype_colors),"-P",sep='')
  merge_subtype_colors = c(merge_premiere_subtype_colors,merge_ucsd_subtype_colors)
  merge_subtype_colors[["Ureter"]] = ureter_color

  center_cell_subtype_colors_list[["SC SN RNASeq merged"]] = merge_subtype_colors
  center_cell_subtype_colors_list[["SC SN RNASeq merged male"]] = merge_subtype_colors
  center_cell_subtype_colors_list[["SC SN RNASeq merged female"]] = merge_subtype_colors

  physiological_cell_subtype_colors = list()
  physiological_cell_subtype_colors[["ProxTub"]] = PTTI_color
  physiological_cell_subtype_colors[["Henle"]] = cortex_henle_color
  physiological_cell_subtype_colors[["CD"]] = cortex_PCCD_color
  physiological_cell_subtype_colors[["Ureter"]] = ureter_color
  center_cell_subtype_colors_list[["Physiological"]] = physiological_cell_subtype_colors
}#End - Define cell subtype colors

if (library_of_interest=="NaAndGluTMTransport")
{#Begin - Set center cell subtype order based on cell type and decreasing sodium transport activity and set first cell subtype of physio segments
  scps_for_sorting = list()
  scps_for_sorting[["NaAndGluTMTransport"]] = "Sodium L2B"
  scp_for_sorting = scps_for_sorting[[library_of_interest]]

  cell_type_order = c("POD","PT","DL","DTL","ATL","TAL","TAL_medulla","TAL_cortex","DCT","CNT","CNT_PC","CD","CD_cortex","PC","PC_cortex","tPC_IC","IC","IC_A","IC_B","IC_cortex","PC_medulla","IC_medulla")

  center_physioSegment_firstCellSubType_array = list()
  center_cell_subtype_orders = list()

  indexRemove = which(!read_counts_cell_subtypes$Cell_type %in% cell_type_order)
  if (length(indexRemove)!=0) { read_counts_cell_subtypes=c(); if (error_message=="") { error_message = paste("Center cell type order contains not all cells: ",unique(read_counts_cell_subtypes$Cell_type[indexRemove]),sep='') } }
  #read_counts_cell_subtypes = read_counts_cell_subtypes[indexKeep,]
  centers = unique(read_counts_cell_subtypes$Center)
  indexCenter=1
  for (indexCenter in 1:length(centers))
  {#Begin
    center = centers[indexCenter]
    indexCurrentCenter = which(read_counts_cell_subtypes$Center == center)
    indexCurrentSCP = which(read_counts_cell_subtypes$SCP == scp_for_sorting)
    indexCurrent = indexCurrentCenter[indexCurrentCenter %in% indexCurrentSCP]
    centerScp_read_counts_cell_subtypes = read_counts_cell_subtypes[indexCurrent,]
    centerScp_read_counts_cell_subtypes = centerScp_read_counts_cell_subtypes[order(centerScp_read_counts_cell_subtypes$Average_read_counts,decreasing=TRUE),]
    centerScp_read_counts_cell_subtypes$Cell_type_factor = factor(centerScp_read_counts_cell_subtypes$Cell_type,levels=cell_type_order)
    centerScp_read_counts_cell_subtypes = centerScp_read_counts_cell_subtypes[order(centerScp_read_counts_cell_subtypes$Cell_type_factor),]
    center_cell_subtype_order = unique(centerScp_read_counts_cell_subtypes$Cell_subtype)
    read_counts_cell_subtypes[indexCurrent,] = centerScp_read_counts_cell_subtypes
    center_cell_subtype_orders[[center]] = center_cell_subtype_order

    physioScp_read_counts_cell_subtypes = centerScp_read_counts_cell_subtypes
    indexKeep = which(!duplicated(physioScp_read_counts_cell_subtypes$Physio_segment))
    physioScp_read_counts_cell_subtypes = physioScp_read_counts_cell_subtypes[indexKeep,]
    center_physioSegment_firstCellSubType = list()
    for (indexPCT in 1:length(physioScp_read_counts_cell_subtypes[,1]))
    {#Begin
      physio_cell_type = physioScp_read_counts_cell_subtypes$Physio_segment[indexPCT]
      center_physioSegment_firstCellSubType[[physio_cell_type]] = physioScp_read_counts_cell_subtypes$Cell_subtype[indexPCT]
    }#End
    center_physioSegment_firstCellSubType[["Ureter"]] = "Ureter"
    center_physioSegment_firstCellSubType_array[[center]] = center_physioSegment_firstCellSubType
  }#End
}#End - Set center cell subtype order based on cell type and decreasing sodium transport activity and set first cell subtype of physio cell types

if (library_of_interest=="NaAndGluTMTransport")
{#Begin - Add physiological measurements
  physio_segments              =                    c("ProxTub","Henle","DistTub","CD","Ureter")
  costanzo_L2B_values = list("Sodium"    = c(    67,     25,       5,     3,    0),
                                      "Potassium" = c(    67,     20,       0,     0,   13),
                                      "Phosphate" = c(    85,      0,       0,     0,   15),
                                      "Calcium"   = c(    67,     25,       8,     0,    0),
                                      "Magnesium" = c(    30,     60,       5,     0,    5))

  physio_segments           =                      c("ProxTub","Henle","DistTub","CD","Ureter")
  berne_L2B_values = list("Sodium"      = c(    67,     20,       7,     5,    1),
                                   "Potassium"   = c(    67,     20,      10,     3,    0),#
                                   "Phosphate"   = c(    80,     -1,      10,    -1,   10),#
                                   "Calcium"     = c(    70,     20,       5,     5,    0),#
                                   "Magnesium"   = c(    30,     65,       1,     1,    3),#
                                   "Bicarbonate" = c(    85,     10,      -1,     5,    0))#

  physio_segments           =                        c("ProxTub","Henle","DistTub","CD","Ureter")
  schmidt_L2B_values = list("Sodium"      = c(       70,     20,        8,   1,   1),
                                     "Chloride"    = c(       70,     20,        8,   1,   1),#
                                     "Water"       = c(       70,     10,     12.5, 6.5,   1),#
                                     "Potassium"   = c(       65,     20,        3,   2,  10),#
                                     "Calcium"     = c(       55,     35,        5,   4,   1),#
                                     "Magnesium"   = c(       30,     50,       10,   0,  10))

  physio_segments                =                       c("ProxTub","Henle","DistTub","CD","Ureter")
  silbernagl_L2B_values = list("Sodium"       = c(    60,     32,       4,    -1,    4),
                                        "Phosphate"    = c(    65,     15,       7,    -1,   13),
                                        "Calcium"      = c(    57,     25,      15,    -1,    3),
                                        "Magnesium"    = c(    22,     60,       5,    -1,   13),
                                        "Water"        = c(    60,     26,      11,    -1,    3),  #Silbernagel
                                        "Valin"        = c(    90,     10,       0,     0,    0),
                                        "Glucose"      = c(    90,     10,       0,     0,    0),
                                        "Glycin"       = c(    90,     10,       0,     0,    0),
                                        "Valin Glycin" = c(    90,     10,       0,     0,    0))


  physio_scp_label = " physio"

  L2B_values_list = list("Costanzo" = costanzo_L2B_values,
                                  "Berne" = berne_L2B_values,
                                  "Silbernagl" = silbernagl_L2B_values,
                                  "Schmidt" = schmidt_L2B_values)
  L2B_value_sources = names(L2B_values_list)

  paracellular_L2B_sodium = list("ProxTub","Henle")
  paracellular_L2B = list("Sodium" = paracellular_L2B_sodium)

  {#Begin - Generate profiles without paracellular sodium L2B
    for (indexReab in 1:length(L2B_value_sources))
    {#Begin
      L2B_value_source = L2B_value_sources[indexReab]
      L2B_values = L2B_values_list[[L2B_value_source]]
      sodium_L2B_noParaPT = L2B_values$Sodium
      sodium_L2B_noParaPT[1] = sodium_L2B_noParaPT[1] * 0.63  #PMID: 28452575: 37%, PMID: 20385797, PMID: 4056034/Chapter 27: 2/3 Cl para, PMID: 6736248/Chapter 27: rabbit: PT: 2/3 <-> 1/3, PMID: 28617689 (32-64% for NaCl), PMID: 20385797: Claudin-2 deficient mice: 37% less paracellular Na+
      sodium_L2B_noParaPT[2] = sodium_L2B_noParaPT[2] * 0.70
      factor = 100 / sum(sodium_L2B_noParaPT)
      sodium_L2B_noParaPT = sodium_L2B_noParaPT * factor
      L2B_values[["Sodium_PTm37_henlem30"]] = sodium_L2B_noParaPT
      L2B_values_list[[L2B_value_source]] = L2B_values
    }#End
  }#End - Generate profiles without paracellular sodium L2B

  physio_reab_ions = names(L2B_values_list)
  L2B_not_100 = c()
  for (indexRList in 1:length(L2B_values_list))
  {#Begin
     L2B_values = L2B_values_list[[indexRList]]
     for (indexRV in 1:length(L2B_values))
     {#Begin
        currentIon_L2B_values = unlist(L2B_values[indexRV])
        indexNegative = which(currentIon_L2B_values < 0)
        currentIon_L2B_values[indexNegative] = 0;
        Sum = sum(unlist(currentIon_L2B_values))
        if (abs(round(Sum)-100)>2) { L2B_not_100=c(L2B_not_100,paste(names(L2B_values_list)[indexRList],"-",names(L2B_values)[indexRV],": ",Sum,sep='')) }
     }#End
  }#End

  physiological_read_counts_cell_subtypes = c()
  for (indexReab in 1:length(L2B_value_sources))
  {#Begin
     L2B_value_source = L2B_value_sources[indexReab]
     L2B_values = L2B_values_list[[L2B_value_source]]
     ions = names(L2B_values)
     if (length(ions)>0)
     {#Begin
        for (indexIon in 1:length(ions))
        {#Begin
          ion =ions[indexIon]
          current_values = L2B_values[[ion]]
          if (length(current_values)>0)
          {#Begin
            for (indexCurrentValues in 1:length(current_values))
            {#Begin
               current_value = current_values[indexCurrentValues]
               if (current_value==-1) { current_value = 0}
               if (current_value!=-1)
               {#Begin
                  physio_segment = physio_segments[indexCurrentValues]
                  #physio_cell_subtype = center_physioSegment_firstCellSubType[[physio_segment]]
                  read_count_cell_subtype_line = read_count_cell_subtype_base_line
                  read_count_cell_subtype_line$Cell_type = physio_segment
                  read_count_cell_subtype_line$Cell_subtype = physio_segment#physio_cell_subtype
                  read_count_cell_subtype_line$Center = L2B_value_source
                  read_count_cell_subtype_line$Physio_segment = physio_segment
                  read_count_cell_subtype_line$Data_type = "2_Molecules"
                  read_count_cell_subtype_line$Cell_subtype_counts = 0
                  read_count_cell_subtype_line$SCP = paste(ion,physio_scp_label,sep='')
                  read_count_cell_subtype_line$Read_counts = current_value
                  read_count_cell_subtype_line$Average_read_counts = current_value
                  if (length(physiological_read_counts_cell_subtypes)==0) { physiological_read_counts_cell_subtypes = read_count_cell_subtype_line }
                  else { physiological_read_counts_cell_subtypes = rbind(physiological_read_counts_cell_subtypes,read_count_cell_subtype_line)}
               }#End
            }#End
          }#End
        }#End
     }#End
  }#End
  read_counts_cell_subtypes = rbind(read_counts_cell_subtypes,physiological_read_counts_cell_subtypes)
  #center_cell_subtype_orders[["Physiological"]] = physiological_cell_subtypes
}#End - Add physiological measurements

if (library_of_interest=="NaAndGluTMTransport")
{#Begin - Plot cell counts
   centers = names(annotations)
   Plots = list()
   Tissue_type_plots = list()
   Tissue_collection_plots = list()
   Participant_plots = list()
   indexC=1
   for (indexC in 1:length(centers))
   {#Begin
      center = centers[indexC]
      center_cell_subtype_order = center_cell_subtype_orders[[center]]
      center_cell_subtype_colors = center_cell_subtype_colors_list[[center]]

      center_annotations = annotations[[center]]
      indexKeep = which(center_annotations$Cluster %in% center_cell_subtype_order)
      center_annotations = center_annotations[indexKeep,]
      center_annotations$Plot_value = 1
      center_annotations$Cell_subtype_factor = factor(center_annotations$Cluster,levels=center_cell_subtype_order)

      libaries_count = length(unique(center_annotations$Library))
      cells_count = length(center_annotations$Barcode)
      participants_count = length(unique(center_annotations$Participant))

      indexCurrentCenter = which(read_counts_cell_subtypes$Center==center)
      center_read_counts = read_counts_cell_subtypes[indexCurrentCenter,]


      center_read_counts$Cell_subtype_factor = factor(center_read_counts$Cell_subtype,levels=center_cell_subtype_order)
      center_read_counts = center_read_counts[order(center_read_counts$Cell_subtype_counts,decreasing=TRUE),]
      center_read_counts = center_read_counts[order(center_read_counts$Cell_subtype_factor),]
      indexUnique = which(!duplicated(center_read_counts$Cell_subtype_factor))
      center_read_counts = center_read_counts[indexUnique,]
      center_read_counts$Cell_subtype_counts_percent = 100*center_read_counts$Cell_subtype_counts/sum(center_read_counts$Cell_subtype_counts)

      indexNotureter = which(center_read_counts$Cell_subtype!="Ureter")
      center_read_counts = center_read_counts[indexNotureter,]

      cell_subtypes_for_color = unique(center_read_counts$Cell_subtype)
      text_colors = replicate(length(cell_subtypes_for_color),"gray77")
      for (indexDataForColor in 1:length(cell_subtypes_for_color))
      {#Begin
        if (length(center_cell_subtype_colors[[cell_subtypes_for_color[indexDataForColor]]])>0)
        { text_colors[indexDataForColor] = center_cell_subtype_colors[[cell_subtypes_for_color[indexDataForColor]]] }
      }#End

      cellsubtype_counts_fileName = paste("Cell_counts_",center,".txt",sep='')
      cellsubtype_counts_completeFileName = paste(results_directory,cellsubtype_counts_fileName,sep='');
      write.table(center_read_counts,file=cellsubtype_counts_completeFileName,quote=FALSE,col.names=TRUE,row.names=TRUE,sep='\t')

      Plot = ggplot(data=center_read_counts,aes(x=Cell_subtype_factor,y=Cell_subtype_counts))
      Plot = Plot + geom_bar(stat = "identity")
      Plot = Plot + theme(axis.text.x = element_text(angle=90,size=15,face=2,vjust=0.5,hjust=1,color=text_colors))
      Plot = Plot + theme(axis.text.y = element_text(size=15,face=2))
      Plot = Plot + ggtitle(paste("Cell counts: ",center, " (",libaries_count," libraries, ",participants_count," participants, ",cells_count," cells)",sep='')) + theme(title = element_text(size=10,face=2,hjust=1))
      Plots[[length(Plots)+1]] = Plot

      Plot = ggplot(data=center_read_counts,aes(x=Cell_subtype_factor,y=Cell_subtype_counts_percent))
      Plot = Plot + geom_bar(stat = "identity")
      Plot = Plot + theme(axis.text.x = element_text(angle=90,size=15,face=2,vjust=0.5,hjust=1,color=text_colors))
      Plot = Plot + theme(axis.text.y = element_text(size=15,face=2))
      Plot = Plot + ggtitle(paste("Cell counts: ",center, " (",libaries_count," libraries, ",participants_count," participants, ",cells_count," cells)",sep='')) + theme(title = element_text(size=10,face=2,hjust=1))
      Plots[[length(Plots)+1]] = Plot

      if (length(is.na(center_annotations$Tissue_type))==0)
      {#Begin
         Tissue_type_plot = ggplot(data=center_annotations,aes(x=Cell_subtype_factor,y=Plot_value,group=Tissue_type,color=Tissue_type))
         Tissue_type_plot = Tissue_type_plot + geom_bar(stat = "identity")
         Tissue_type_plot = Tissue_type_plot + theme(axis.text.x = element_text(angle=90,size=15,face=2,vjust=0.5,hjust=1,color=text_colors))
         Tissue_type_plot = Tissue_type_plot + theme(axis.text.y = element_text(size=15,face=2))
         Tissue_type_plot = Tissue_type_plot + ggtitle(paste("Cell counts: ",center, " (",libaries_count," libraries, ",participants_count," participants, ",cells_count," cells)",sep='')) + theme(title = element_text(size=10,face=2,hjust=1))
         Tissue_type_plots[[length(Tissue_type_plots)+1]] = Tissue_type_plot
      }#End

      Tissue_collection_plot = ggplot(data=center_annotations,aes(x=Cell_subtype_factor,y=Plot_value,group=Tissue_collection,color=Tissue_collection))
      Tissue_collection_plot = Tissue_collection_plot + geom_bar(stat = "identity")
      Tissue_collection_plot = Tissue_collection_plot + theme(axis.text.x = element_text(angle=90,size=15,face=2,vjust=0.5,hjust=1,color=text_colors))
      Tissue_collection_plot = Tissue_collection_plot + theme(axis.text.y = element_text(size=15,face=2))
      Tissue_collection_plot = Tissue_collection_plot + ggtitle(paste("Cell counts: ",center, " (",libaries_count," libraries, ",participants_count," participants, ",cells_count," cells)",sep='')) + theme(title = element_text(size=10,face=2,hjust=1))

      Tissue_collection_plots[[length(Tissue_collection_plots)+1]] = Tissue_collection_plot

      Participant_plot = ggplot(data=center_annotations,aes(x=Cell_subtype_factor,y=Plot_value,group=Participant,color=Participant))
      Participant_plot = Participant_plot + geom_bar(stat = "identity")
      Participant_plot = Participant_plot + theme(axis.text.x = element_text(angle=90,size=15,face=2,vjust=0.5,hjust=1,color=text_colors))
      Participant_plot = Participant_plot + theme(axis.text.y = element_text(size=15,face=2))
      Participant_plot = Participant_plot + ggtitle(paste("Cell counts: ",center, " (",libaries_count," libraries, ",participants_count," participants, ",cells_count," cells)",sep='')) + theme(title = element_text(size=10,face=2,hjust=1))

      Participant_plots[[length(Participant_plots)+1]] = Participant_plot

   }#End

   max_plots_per_figure = 8;
   figures_count = ceiling(length(Plots)/max_plots_per_figure)
   for (indexF in 1:figures_count)
   {#Begin
     startPlot = (indexF-1)*max_plots_per_figure+1
     endPlot = min(indexF*max_plots_per_figure,length(Plots))
     current_plots = Plots[startPlot:endPlot]
     rows_count = min(length(centers),length(current_plots))
     cols_count = ceiling(length(current_plots)/rows_count)
     complete_png_fileName = paste(results_directory,"Nephron_cell_counts_",L2B_value_source,"_no",indexF,".png",sep='')
     png(complete_png_fileName,width=700*cols_count,height=400*rows_count,res=100);
     do.call("grid.arrange",c(current_plots,nrow=rows_count,ncol=cols_count))
     dev.off()
   }#End

   max_plots_per_figure = 8;
   figures_count = ceiling(length(Tissue_type_plots)/max_plots_per_figure)
   if (figures_count>0)
   {#Begin
     for (indexF in 1:figures_count)
     {#Begin
       startPlot = (indexF-1)*max_plots_per_figure+1
       endPlot = min(indexF*max_plots_per_figure,length(Tissue_type_plots))
       current_plots = Tissue_type_plots[startPlot:endPlot]
       rows_count = min(length(centers),length(current_plots))
       cols_count = ceiling(length(current_plots)/rows_count)
       complete_png_fileName = paste(results_directory,"Nephron_cell_counts_tissue_type_",L2B_value_source,"_no",indexF,".png",sep='')
       png(complete_png_fileName,width=700*cols_count,height=400*rows_count,res=100);
       do.call("grid.arrange",c(current_plots,nrow=rows_count,ncol=cols_count))
       dev.off()
     }#End
   }#End

   max_plots_per_figure =8;
   figures_count = ceiling(length(Tissue_collection_plots)/max_plots_per_figure)
   for (indexF in 1:figures_count)
   {#Begin
     startPlot = (indexF-1)*max_plots_per_figure+1
     endPlot = min(indexF*max_plots_per_figure,length(Tissue_collection_plots))
     current_plots = Tissue_collection_plots[startPlot:endPlot]
     cols_count = min(length(centers),length(current_plots))
     rows_count = ceiling(length(current_plots)/cols_count)
     complete_png_fileName = paste(results_directory,"Nephron_cell_counts_tissue_collection_",L2B_value_source,"_no",indexF,".png",sep='')
     png(complete_png_fileName,width=700*cols_count,height=400*rows_count,res=100);
     do.call("grid.arrange",c(current_plots,nrow=rows_count,ncol=cols_count))
     dev.off()
   }#End

   max_plots_per_figure = 8;
   figures_count = ceiling(length(Participant_plots)/max_plots_per_figure)
   for (indexF in 1:figures_count)
   {#Begin
     startPlot = (indexF-1)*max_plots_per_figure+1
     endPlot = min(indexF*max_plots_per_figure,length(Participant_plots))
     current_plots = Participant_plots[startPlot:endPlot]
     cols_count = min(length(centers),length(current_plots))
     rows_count = ceiling(length(current_plots)/cols_count)
     complete_png_fileName = paste(results_directory,"Nephron_cell_counts_participants_",L2B_value_source,"_no",indexF,".png",sep='')
     png(complete_png_fileName,width=700*cols_count,height=400*rows_count,res=100);
     do.call("grid.arrange",c(current_plots,nrow=rows_count,ncol=cols_count))
     dev.off()
   }#End
}#End - Plot cell counts

{#Begin - Define plot colors
  alternating_colors_proxTub = c("orange","darkorange2")
  alternating_colors_henle = c("dodgerblue2","cyan2","blue")
  alternating_colors_distTub = c("orange","darkorange2")
  alternating_colors_cd = c("mediumorchid1","darkorchid2")

  frameColor_proxTub = "darkorange4"
  frameColor_henle = "blue4"
  frameColor_distTub = "darkorange4"
  frameColor_cd = "darkorchid4"
  frameColor_ureter = "gray30"

  {#Begin - Define color sets
    color_sets = list()
    color_sets[[1]] = c("blue")
    color_sets[[2]] = c("blue","orange")
    color_sets[[3]] = c("blue","orange","darkorchid")
    color_sets[[4]] = c("blue","orange","darkorchid","firebrick")
    color_sets[[5]] = c("blue","orange","darkorchid","firebrick","dodgerblue")
    color_sets[[6]] = c("blue","orange","darkorchid","firebrick","dodgerblue","orange3")
    color_sets[[7]] = c("blue","orange","darkorchid","firebrick","dodgerblue","orange3","navy")
    color_sets[[8]] = c("blue","orange","darkorchid","firebrick","dodgerblue","orange3","navy","yellow")
    color_sets[[9]] = c("blue","orange","darkorchid","firebrick","dodgerblue","orange3","navy","yellow","dodgerblue4")
    color_sets[[10]] = c("blue","orange","darkorchid","firebrick","dodgerblue","orange3","navy","yellow","dodgerblue4","plum")
    color_sets[[11]] = c("blue","orange","darkorchid","firebrick","dodgerblue","orange3","navy","yellow","dodgerblue4","plum","mediumorchid")
    color_sets[[12]] = c("orange","dodgerblue","mediumorchid","firebrick","dodgerblue","orange3","navy","yellow","dodgerblue4","plum","mediumorchid","gold","blue")
    color_sets[[13]] = c("orange","dodgerblue","mediumorchid","firebrick","dodgerblue","orange3","navy","yellow","dodgerblue4","plum","mediumorchid","gold","blue","mediumorchid1")
    color_sets[[14]] = c("orange","dodgerblue","mediumorchid","firebrick","dodgerblue","orange3","navy","yellow","dodgerblue4","plum","mediumorchid","gold","blue","mediumorchid1","cyan3")
    color_sets[[15]] = c("orange","dodgerblue","mediumorchid","firebrick","dodgerblue","orange3","navy","yellow","dodgerblue4","plum","mediumorchid","gold","blue","mediumorchid1","cyan3","pink")
    color_sets[[16]] = c("orange","dodgerblue","mediumorchid","firebrick","dodgerblue","orange3","navy","yellow","dodgerblue4","plum","mediumorchid","gold","blue","mediumorchid1","cyan3","pink","blueviolet")
    color_sets[[17]] = c("orange","dodgerblue","mediumorchid","firebrick","dodgerblue","orange3","navy","yellow","dodgerblue4","plum","mediumorchid","gold","blue","mediumorchid1","cyan3","pink","blueviolet","chocolate3")
    color_sets[[18]] = c("orange","lightskyblue","mediumorchid","navy","dodgerblue","firebrick","lightskyblue4","orange3","cyan","mediumorchid4","gold","blue","plum","darkred","cyan3","pink","blueviolet","chocolate3","yellow3")
    color_sets[[19]] = c("orange","dodgerblue","mediumorchid","firebrick","dodgerblue","orange3","navy","yellow","dodgerblue4","plum","mediumorchid","gold","blue","mediumorchid1","cyan3","pink","blueviolet","chocolate3","yellow3","aquamarine")
    color_sets[[20]] = c("orange","dodgerblue","mediumorchid","firebrick","dodgerblue","orange3","navy","yellow","dodgerblue4","plum","mediumorchid","gold","blue","mediumorchid1","cyan3","pink","blueviolet","chocolate3","yellow3","aquamarine","orange3")
    color_sets[[21]] = c("orange","dodgerblue","mediumorchid","firebrick","dodgerblue","orange3","navy","yellow","dodgerblue4","plum","mediumorchid","gold","blue","mediumorchid1","cyan3","pink","blueviolet","chocolate3","yellow3","aquamarine","orange3","cyan")
    color_sets[[22]] = c("orange","dodgerblue","mediumorchid","firebrick","dodgerblue","orange3","navy","yellow","dodgerblue4","plum","mediumorchid","gold","blue","mediumorchid1","cyan3","pink","blueviolet","chocolate3","yellow3","aquamarine","orange3","cyan","mediumorchid1")
    color_sets[[23]] = c("orange","dodgerblue","mediumorchid","firebrick","dodgerblue","orange3","navy","yellow","dodgerblue4","plum","mediumorchid","gold","blue","mediumorchid1","cyan3","pink","blueviolet","chocolate3","yellow3","aquamarine","orange3","cyan","darkblue","mediumorchid1")
    color_sets[[24]] = c("orange","dodgerblue","mediumorchid","firebrick","dodgerblue","orange3","navy","yellow","dodgerblue4","plum","mediumorchid","gold","blue","mediumorchid1","cyan3","pink","blueviolet","chocolate3","yellow3","aquamarine","orange3","cyan","darkblue","mediumorchid1","plum")
    color_sets[[25]] = c("orange","dodgerblue","mediumorchid","firebrick","dodgerblue","orange3","navy","yellow","dodgerblue4","plum","mediumorchid","gold","blue","mediumorchid1","cyan3","pink","blueviolet","chocolate3","yellow3","aquamarine","orange3","cyan","darkblue","mediumorchid1","plum","darkorchid4")
    color_sets[[26]] = c("orange","dodgerblue","mediumorchid","firebrick","dodgerblue","orange3","navy","yellow","dodgerblue4","plum","mediumorchid","gold","blue","mediumorchid1","cyan3","pink","blueviolet","chocolate3","yellow3","aquamarine","orange3","cyan","darkblue","mediumorchid1","plum","darkorchid4","yellow3")
    color_sets[[27]] = c("orange","dodgerblue","mediumorchid","firebrick","dodgerblue","orange3","navy","yellow","dodgerblue4","plum","mediumorchid","gold","blue","mediumorchid1","cyan3","pink","blueviolet","chocolate3","yellow3","aquamarine","orange3","cyan","darkblue","mediumorchid1","plum","darkorchid4","yellow3","navy")
    color_sets[[28]] = c("orange","dodgerblue","mediumorchid","firebrick","dodgerblue","orange3","navy","yellow","dodgerblue4","plum","mediumorchid","gold","blue","mediumorchid1","cyan3","pink","blueviolet","chocolate3","yellow3","aquamarine","orange3","cyan","darkblue","mediumorchid1","plum","darkorchid4","yellow3","navy","firebrick2")
    color_sets[[29]] = c("orange","dodgerblue","mediumorchid","firebrick","orange3","darkorchid","yellow","navy","dodgerblue4","plum","mediumorchid","gold","blue","cyan","mediumorchid1","pink","blueviolet","chocolate3","yellow3","aquamarine3","orange3","cyan3","darkblue","mediumorchid1","plum","darkorchid4","yellow3","navy","firebrick2","darkorchid3")
    color_sets[[30]] = c("orange","dodgerblue","mediumorchid","firebrick","dodgerblue","orange3","navy","yellow","dodgerblue4","plum","mediumorchid","gold","blue","mediumorchid1","cyan3","pink","blueviolet","chocolate3","yellow3","aquamarine","orange3","cyan","darkblue","mediumorchid1","plum","darkorchid4","yellow3","navy","firebrick2","darkorchid3","blue")
    color_sets[[31]] = c("orange","dodgerblue","mediumorchid","firebrick","dodgerblue","orange3","navy","yellow","dodgerblue4","plum","mediumorchid","gold","blue","mediumorchid1","cyan3","pink","blueviolet","chocolate3","yellow3","aquamarine","orange3","cyan","darkblue","mediumorchid1","plum","darkorchid4","yellow3","navy","firebrick2","darkorchid3","blue","plum3")
    color_sets[[32]] = c("orange","dodgerblue","mediumorchid","firebrick","dodgerblue","orange3","navy","yellow","dodgerblue4","plum","mediumorchid","gold","blue","mediumorchid1","cyan3","pink","blueviolet","chocolate3","yellow3","aquamarine","orange3","cyan","darkblue","mediumorchid1","plum","darkorchid4","yellow3","navy","firebrick2","darkorchid3","plum3","blue","orangered")
    color_sets[[33]] = c("orange","dodgerblue","mediumorchid","firebrick","dodgerblue","orange3","navy","yellow","dodgerblue4","plum","mediumorchid","gold","blue","mediumorchid1","cyan3","pink","blueviolet","chocolate3","yellow3","aquamarine","orange3","cyan","darkblue","mediumorchid1","plum","darkorchid4","yellow3","navy","firebrick2","darkorchid3","plum3","blue","orangered","skyblue")
    color_sets[[34]] = c("orange","dodgerblue","mediumorchid","firebrick","dodgerblue","orange3","navy","yellow","dodgerblue4","plum","mediumorchid","gold","blue","mediumorchid1","cyan3","pink","blueviolet","chocolate3","yellow3","aquamarine","orange3","cyan","darkblue","mediumorchid1","plum","darkorchid4","yellow3","navy","firebrick2","darkorchid3","plum3","blue","orangered","skyblue","gray")
    color_sets[[35]] = c("orange","dodgerblue","mediumorchid","firebrick","dodgerblue","orange3","navy","yellow","dodgerblue4","plum","mediumorchid","gold","blue","mediumorchid1","cyan3","pink","blueviolet","chocolate3","yellow3","aquamarine","orange3","cyan","darkblue","mediumorchid1","plum","darkorchid4","yellow3","navy","firebrick2","darkorchid3","plum3","blue","orangered","skyblue","gray","limegreen")
  }#End - Define color sets
}#End - Define plot colors



###############################################################################
##################### Segment specific net L2B

calculate_percent_read_counts_compared_to_reference_scp = TRUE
if ((algorithm_of_interest == "Total_sum") & (calculate_percent_read_counts_compared_to_reference_scp))
{#Begin - Calculate read counts in percent of net L2B for all scps within each scp_group over nephron
  nephron_percent_read_counts = c()
  merged_scp_groups = c(selected_SCPs_for_barplot,selected_SCPs_for_scp_subgroups)
  scp_groups = names(merged_scp_groups)
  for (indexScpGroup in 1:length(scp_groups))
  {#Begin
    scp_group = scp_groups[indexScpGroup]
    splitStrings = strsplit(scp_group," ")[[1]]
    ion = splitStrings[1]
    if (length(splitStrings)-2>1)
    {#Begin
       for (indexSplit in 2:(length(splitStrings)-2))
       {#Begin
          ion = paste(ion," ",splitStrings[indexSplit],sep='')
       }#End
    }#End
    reference_scp = paste(ion,net_L2B_label,sep='')
    indexCurrent_reference_scp = which(read_counts_cell_subtypes$SCP==reference_scp)
    reference_read_counts = read_counts_cell_subtypes[indexCurrent_reference_scp,]
    scps = merged_scp_groups[[scp_group]]
    indexCurrentGroup = which(read_counts_cell_subtypes$SCP %in% scps)
    scpGroup_read_counts = read_counts_cell_subtypes[indexCurrentGroup,]
    centers = unique(scpGroup_read_counts$Center)
    for (indexCenter in 1:length(centers))
    {#Begin
      center = centers[indexCenter]
      indexCurrentCenter = which(scpGroup_read_counts$Center==center)
      center_read_counts = scpGroup_read_counts[indexCurrentCenter,]
      indexCurrentCenter_reference = which(reference_read_counts$Center==center)
      center_reference_read_counts = reference_read_counts[indexCurrentCenter_reference,]
      dataTypes = unique(center_read_counts$Data_type)
      for (indexDataType in 1:length(dataTypes))
      {#Begin
        dataType = dataTypes[indexDataType]
        indexCurrentDataType = which(center_read_counts$Data_type==dataType)
        indexCurrentDataType_reference = which(center_reference_read_counts$Data_type==dataType)
        if (length(indexCurrentDataType)>0)
        {#Begin
          datatype_read_counts = center_read_counts[indexCurrentDataType,]
          if (length(indexCurrentDataType_reference)>0)
          {#Begin
             sum_reference_readCounts = abs(sum(center_reference_read_counts$Read_counts[indexCurrentDataType_reference]))
             datatype_read_counts$Reference = "Nephron net L2B"
          }#End
          else
          {#Begin
             sum_reference_readCounts = sum(datatype_read_counts$Read_counts)
             datatype_read_counts$Reference = "Nephron same scps"
          }#End
          if (sum_reference_readCounts!=0)
          { datatype_read_counts$Percent_read_counts_nephron = 100*datatype_read_counts$Read_counts/sum_reference_readCounts }
          else { datatype_read_counts$Percent_read_counts_nephron = 0 }
          datatype_read_counts$Scp_group = scp_group
          if (length(scps)==1) { datatype_read_counts$Group_scps = scps }
          else { datatype_read_counts$Group_scps = paste(scps,collapse=";") }
          if (length(nephron_percent_read_counts)==0) { nephron_percent_read_counts = datatype_read_counts }
          else { nephron_percent_read_counts = rbind(nephron_percent_read_counts,datatype_read_counts) }
        }#Begin
      }#End
    }#End
  }#End
}#End - Calculate read counts in percent of net L2B for all scps within each scp_group over nephron

if (algorithm_of_interest == "Total_sum")
{#Begin - Generate dataset with read counts of each scp (of each scp_group) combined for each segment

Col_names = c("Physio_segment","Center","SCP_group","Group_scps","SCP","Data_type","Read_counts","Reference","Percent_read_counts_nephron")
Col_length = length(Col_names)
Row_names = 1
Row_length = length(Row_names)
segment_scpGroup_read_count_base_line = array(NA,c(Row_length,Col_length),dimnames=list(Row_names,Col_names))
segment_scpGroup_read_count_base_line = as.data.frame(segment_scpGroup_read_count_base_line)
segment_scpGroup_read_counts = c()

centers = unique(nephron_percent_read_counts$Center)
for (indexC in 1:length(centers))
{#Begin
   center = centers[indexC]
   indexCurrentCenter = which(nephron_percent_read_counts$Center==center)
   center_read_counts = nephron_percent_read_counts[indexCurrentCenter,]
   dataTypes = unique(center_read_counts$Data_type)
   for (indexD in 1:length(dataTypes))
   {#Begin
      dataType = dataTypes[indexD]
      indexCurrentDataType = which(center_read_counts$Data_type==dataType)
      datatype_read_counts = center_read_counts[indexCurrentDataType,]
      scpGroups = unique(datatype_read_counts$Scp_group)
      for (indexScpGroup in 1:length(scpGroups))
      {#Begin
         scpGroup = scpGroups[indexScpGroup]
         indexCurrentScpGroup = which(datatype_read_counts$Scp_group==scpGroup)
         scpGroup_read_counts = datatype_read_counts[indexCurrentScpGroup,]
         scps = unique(scpGroup_read_counts$SCP)
         for (indexScp in 1:length(scps))
         {#Begin
            scp = scps[indexScp]
            indexCurrentScp = which(scpGroup_read_counts$SCP==scp)
            scp_read_counts = scpGroup_read_counts[indexCurrentScp,]
            segments = unique(scp_read_counts$Physio_segment)
            for (indexS in 1:length(segments))
            {#Begin
               segment = segments[indexS]
               indexCurrentSegment = which(scp_read_counts$Physio_segment==segment)
               segment_read_counts = scp_read_counts[indexCurrentSegment,]

               segment_scpGroup_read_count_line = segment_scpGroup_read_count_base_line
               segment_scpGroup_read_count_line$Reference = "Nephron"
               segment_scpGroup_read_count_line$Physio_segment = segment
               segment_scpGroup_read_count_line$Center = center
               segment_scpGroup_read_count_line$SCP = scp
               segment_scpGroup_read_count_line$Data_type = dataType
               segment_scpGroup_read_count_line$Read_counts = sum(segment_read_counts$Read_counts)
               segment_scpGroup_read_count_line$Percent_read_counts_nephron = sum(segment_read_counts$Percent_read_counts_nephron)
               current_scp_group = unique(segment_read_counts$Scp_group)
               if (length(current_scp_group)==1) { segment_scpGroup_read_count_line$SCP_group = current_scp_group }
               else {segment_scpGroup_read_count_line$SCP_group = paste(current_scp_group,collapse=";")}
               current_groupScps = unique(segment_read_counts$Group_scps)
               if (length(current_groupScps)==1) { segment_scpGroup_read_count_line$Group_scps = current_groupScps }
               else {segment_scpGroup_read_count_line$Group_scps = paste(current_groupScps,collapse=";")}

               if (length(segment_scpGroup_read_counts)==0) { segment_scpGroup_read_counts = segment_scpGroup_read_count_line }
               else { segment_scpGroup_read_counts = rbind(segment_scpGroup_read_counts,segment_scpGroup_read_count_line)}
            }#End
         }#End
      }#End
   }#End
}#End

}#End - Generate dataset with read counts of each scp (of each scp_group) combined for each segment

if (algorithm_of_interest == "Total_sum")
{#Begin - Calculate mean, sd error and eventually p-values of segment read counts for each SCP in each Scp_group
  segment_scpGroup_read_counts_for_mean = segment_scpGroup_read_counts
  length_centers = length(keep_centers_for_each_analysis)
  keep_centers = c(keep_centers_for_each_analysis,L2B_value_sources)

  indexKeep = which(segment_scpGroup_read_counts_for_mean$Center %in% keep_centers)
  segment_scpGroup_read_counts_for_mean = segment_scpGroup_read_counts_for_mean[indexKeep,]
  if (length(which(!keep_centers %in% segment_scpGroup_read_counts_for_mean$Center))>0) { segment_scpGroup_read_counts_for_mean = c() }
  datatype_considered_centers_count = list()
  dataTypes = unique(segment_scpGroup_read_counts_for_mean$Data_type)
  for (indexDT in 1:length(dataTypes))
  {#Begin
     dataType = dataTypes[indexDT]
     indexCurrentDataType = which(segment_scpGroup_read_counts$Data_type==dataType)
     datatype_considered_centers_count[[dataType]] = length(unique(segment_scpGroup_read_counts$Center[indexCurrentDataType]))
  }#End

  Col_names = c("Physio_segment","Centers","SCP","SCP_group","Group_scps","Data_type","Reference","Mean_percent_read_counts","SD_percent_read_counts","SD_error_percent_read_counts","Median_percent_read_counts","Considered_centers_count","Pvalue","Averaged_centers")
  Col_length = length(Col_names)
  Row_names = 1
  Row_length = length(Row_names)
  segment_scpGroup_average_read_counts_base_line = array(NA,c(Row_length,Col_length),dimnames=list(Row_names,Col_names))
  segment_scpGroup_average_read_counts_base_line = as.data.frame(segment_scpGroup_average_read_counts_base_line)
  segment_scpGroup_average_read_counts = c()

  segments = unique(segment_scpGroup_read_counts_for_mean$Physio_segment)
  indexS=1
  for (indexS in 1:length(segments))
  {#Begin
     segment = segments[indexS]
     indexCurrentSegment = which(segment_scpGroup_read_counts_for_mean$Physio_segment==segment)
     segment_read_counts = segment_scpGroup_read_counts_for_mean[indexCurrentSegment,]
     dataTypes = unique(segment_read_counts$Data_type)
     indexD=1
     for (indexD in 1:length(dataTypes))
     {#Begin
        dataType = dataTypes[indexD]
        indexCurrentDataType = which(segment_read_counts$Data_type==dataType)
        datatype_read_counts = segment_read_counts[indexCurrentDataType,]
        scp_groups = unique(datatype_read_counts$SCP_group)
        if (length(scp_groups)!=length(unique(scp_groups))) { datatype_read_counts = NULL }
        indexScpGroup=3
        for (indexScpGroup in 1:length(scp_groups))
        {#Begin
           scp_group = scp_groups[indexScpGroup]
           indexCurrentScpGroups = which(datatype_read_counts$SCP_group==scp_group)
           scpGroup_read_counts = datatype_read_counts[indexCurrentScpGroups,]
           if (dataType=="1_Read counts") { current_centers_count = datatype_considered_centers_count[[dataType]] }
           if (dataType=="2_Molecules") { current_centers_count = length(unique(scpGroup_read_counts$Center)) }
           scps = unique(scpGroup_read_counts$SCP)
           indexScp=1
           for (indexScp in 1:length(scps))
           {#Begin
              scp = scps[indexScp]
              indexCurrentScp = which(scpGroup_read_counts$SCP==scp)
              scp_read_counts = scpGroup_read_counts[indexCurrentScp,]
              percent_read_counts = scp_read_counts$Percent_read_counts
              if (length(percent_read_counts) > current_centers_count) { percent_read_counts = replicate(current_centers_count,-999999999)}
              while (length(percent_read_counts) < current_centers_count) { percent_read_counts = c(percent_read_counts,0) }
              segment_scpGroup_average_read_counts_line = segment_scpGroup_average_read_counts_base_line
              segment_scpGroup_average_read_counts_line$Physio_segment = segment
              segment_scpGroup_average_read_counts_line$Centers = paste(scp_read_counts$Center,collapse=";")
              segment_scpGroup_average_read_counts_line$SCP = scp

              segment_scpGroup_average_read_counts_line$SCP_group = scp_group
              current_groupScps = unique(scp_read_counts$Group_scps)
              if (length(current_groupScps)==1) { segment_scpGroup_average_read_counts_line$Group_scps = current_groupScps }
              else { segment_scpGroup_average_read_counts_line$Group_scps = paste(current_groupScps,collapse=";") }
              segment_scpGroup_average_read_counts_line$Data_type = dataType
              current_reference = unique(scp_read_counts$Reference)
              if (length(current_reference)==1) { segment_scpGroup_average_read_counts_line$Reference = current_reference }
              else { segment_scpGroup_average_read_counts_line$Reference = paste(current_reference,collapse=";") }
              segment_scpGroup_average_read_counts_line$Mean_percent_read_counts = mean(percent_read_counts)
              segment_scpGroup_average_read_counts_line$Averaged_centers = length(percent_read_counts)
              if (length(scp_read_counts[1,])>1)
              {#Begin
                 segment_scpGroup_average_read_counts_line$SD_percent_read_counts = sd(percent_read_counts)
                 segment_scpGroup_average_read_counts_line$SD_error_percent_read_counts = segment_scpGroup_average_read_counts_line$SD_percent_read_counts / sqrt(length(percent_read_counts))
              }#End
              else
              {#Begin
                segment_scpGroup_average_read_counts_line$SD_percent_read_counts = 0
                segment_scpGroup_average_read_counts_line$SD_error_percent_read_counts = 0
              }#End
              segment_scpGroup_average_read_counts_line$Median_percent_read_counts = median(percent_read_counts)
              segment_scpGroup_average_read_counts_line$Considered_centers_count = current_centers_count
              segment_scpGroup_average_read_counts_line$Pvalue = -1

              if (grepl(net_L2B_label,scp))
              {#Begin
                 ion = str_replace(scp_group,net_L2B_label,"")
                 ion_physio_scp = paste(ion,physio_scp_label,sep='')
                 indexPhysioScp = which(segment_read_counts$SCP==ion_physio_scp)
                 indexGroup = which(segment_read_counts$SCP_group==scp_group)
                 indexPhysioScpGroup = indexPhysioScp[indexPhysioScp %in% indexGroup] #for physio differences in n in the different segments are accepted
                 if ((length(indexPhysioScpGroup)>1)&(length(percent_read_counts)>1))
                 {#Begin
                    current_physio_read_counts = segment_read_counts[indexPhysioScpGroup,]
                    ttest = t.test(x=scp_read_counts$Percent_read_counts,y=percent_read_counts,alternative="two.sided",var.equal = FALSE,paired=FALSE)
                    segment_scpGroup_average_read_counts_line$Pvalue = ttest$p.value
                 }#End
              }#End
              if (length(segment_scpGroup_average_read_counts)==0) { segment_scpGroup_average_read_counts = segment_scpGroup_average_read_counts_line}
              else { segment_scpGroup_average_read_counts = rbind(segment_scpGroup_average_read_counts,segment_scpGroup_average_read_counts_line)}
           }#End
        }#End
     }#End
  }#End

  indexReadCounts = which(segment_scpGroup_average_read_counts$Data_type=="1_Read counts")
  indexMolecules = which(segment_scpGroup_average_read_counts$Data_type=="2_Molecules")
  if (length(unique(segment_scpGroup_average_read_counts$Averaged_centers[indexReadCounts]))!=1) { segment_scpGroup_average_read_counts= c() }
  if (unique(segment_scpGroup_average_read_counts$Averaged_centers[indexReadCounts])!=length_centers) { segment_scpGroup_average_read_counts= c() }
  if (length(unique(segment_scpGroup_average_read_counts$Averaged_centers[indexMolecules]))==1) { segment_scpGroup_average_read_counts= c() }
  segment_scpGroup_average_read_counts$SD_percent_read_counts = NULL
}#End - Calculate mean and sd error of segment read counts for each SCP in each Scp_group

segmentMean_fileName = "Read_counts_mean_segment.txt"
complete_segmentMean_fileName = paste(results_directory,segmentMean_fileName,sep='')
write.table(segment_scpGroup_average_read_counts,file=complete_segmentMean_fileName,row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)

pdf_width = 8   #8.27 = DIN A4
pdf_height = 6 #11.69 = DIN A4


segment_fixed_max_y = 100
#segment_fixed_min_y = -40

if (algorithm_of_interest == "Total_sum")
{#Begin - Visualize mean and sd of segment percent read counts
  scp_groups = unique(segment_scpGroup_average_read_counts$SCP_group)
  scp_groups = names(selected_SCPs_for_barplot)
  scp_groups = scp_groups[order(scp_groups)]
  equal_length_text = "                                                            "

  if (library_of_interest=="NaAndGluTMTransport")
  {#Begin
  scp_groups = c( "Cell counts","Sodium net L2B"
                 ,"Sodium_PTm37_henlem30 net L2B","Glucose net L2B"
                 )
  }#End

  segments_in_correct_order = c("ProxTub","Henle","DistTub","CD","Ureter")
  Barplots = list()
  indexScpGroup=3
  for (indexScpGroup in 1:length(scp_groups))
  {#Begin
     scp_group = scp_groups[indexScpGroup]
     indexCurrentGroup = which(segment_scpGroup_average_read_counts$SCP_group==scp_group)
     if (length(indexCurrentGroup)==0) { Barplots[[length(Barplots)+1]] = ggplot() + geom_blank() }
     else
     {#Begin
        scpGroup_read_counts = segment_scpGroup_average_read_counts[indexCurrentGroup,]
        scpGroup_read_counts$Physio_segment_dataType = paste(scpGroup_read_counts$Physio_segment,"\n",scpGroup_read_counts$Data_type,"\n",scpGroup_read_counts$SCP,sep='')
        scpGroup_read_counts$Physio_segment_dataType = paste(equal_length_text,scpGroup_read_counts$Physio_segment_dataType,equal_length_text,sep='\n')
        scpGroup_read_counts$Physio_segment_factor = factor(scpGroup_read_counts$Physio_segment,levels=segments_in_correct_order)
        scpGroup_read_counts = scpGroup_read_counts[order(scpGroup_read_counts$Data_type,decreasing=TRUE),]
        scpGroup_read_counts = scpGroup_read_counts[order(scpGroup_read_counts$Physio_segment_factor),]
        scpGroup_read_counts$Physio_segment_dataType = factor(scpGroup_read_counts$Physio_segment_dataType,levels=scpGroup_read_counts$Physio_segment_dataType)

        Colors = replicate(length(scpGroup_read_counts[,1]),"gray")
        Frame_colors = replicate(length(scpGroup_read_counts[,1]),"gray")
        for (indexCurrent in 1:length(scpGroup_read_counts[,1]))
        {#Begin - set alternating colors
           if (grepl("ProxTub",scpGroup_read_counts$Physio_segment_dataType[indexCurrent]))
           {#Begin
              Colors[indexCurrent] = "orange"
              Frame_colors[indexCurrent] = "orange4"
           }#End
           if (grepl("Henle",scpGroup_read_counts$Physio_segment_dataType[indexCurrent]))
           {#Begin
              Colors[indexCurrent] = "dodgerblue2"
              Frame_colors[indexCurrent] = "dodgerblue4"
           }#End
           if (grepl("DistTub",scpGroup_read_counts$Physio_segment_dataType[indexCurrent]))
           {#Begin
              Colors[indexCurrent] = "orange"
              Frame_colors[indexCurrent] = "orange4"
           }#End
           if (grepl("CD",scpGroup_read_counts$Physio_segment_dataType[indexCurrent]))
           {#Begin
             Colors[indexCurrent] = "mediumorchid"
             Frame_colors[indexCurrent] = "mediumorchid4"
           }#End
           if (grepl("Ureter",scpGroup_read_counts$Physio_segment_dataType[indexCurrent]))
           {#Begin
             Frame_colors[indexCurrent] = "gray30"
           }#End
           if (grepl("_Molecules",scpGroup_read_counts$Physio_segment_dataType[indexCurrent]))
           {#Begin
              Colors[indexCurrent] = "gray"
           }#End
        }#End

        centers = strsplit(paste(unique(scpGroup_read_counts$Centers),collapse=";"),";")[[1]]
        centers = unique(centers)
        centers = gsub("SN RNASeq ","",centers)
        centers = gsub("SC RNASeq ","",centers)
        centers = gsub("LMD RNASeq ","",centers)
        center_string = ""
        max_centers_in_one_line = round(length(centers)/2)
        centers_in_one_line = 0
        for (indexCenter in 1:length(centers))
        {#Begin
           current_center = centers[indexCenter]
           if (indexCenter!=1) { center_string = paste(center_string,", ",sep='')}
           if (centers_in_one_line==max_centers_in_one_line)
           {#Begin
              center_string = paste(center_string,"\n",current_center,sep='');
              centers_in_one_line = 1
           }#End
           else
           {#Begin
              center_string = paste(center_string,current_center,sep='')
              centers_in_one_line = centers_in_one_line + 1
           }#End
        }#End

        scpGroup_read_counts = scpGroup_read_counts[order(scpGroup_read_counts$Physio_segment_dataType),]
        indexNonNegativePvalues = which(scpGroup_read_counts$Pvalue!=-1)
        pvalues_string = ""
        if (length(indexNonNegativePvalues)>0)
        {#Begin
          pvalues_string = paste("pvalues:",round(scpGroup_read_counts$Pvalue[indexNonNegativePvalues]*100)/100,sep='',collapse=' ')
        }#End


        plot_title = paste(scp_group,"\n",center_string,"\n",pvalues_string,sep='')

        means = scpGroup_read_counts$Mean_percent_read_counts
        indexPostiveMeans = which(means>0)
        indexNegativeMeans = which(means<0)
        scpGroup_read_counts$SD_error_percent_read_counts_lower_value = scpGroup_read_counts$Mean_percent_read_counts;
        scpGroup_read_counts$SD_error_percent_read_counts_upper_value = scpGroup_read_counts$Mean_percent_read_counts + scpGroup_read_counts$SD_error_percent_read_counts;
        scpGroup_read_counts$SD_error_percent_read_counts_lower_value[indexNegativeMeans] = 0;
        scpGroup_read_counts$SD_error_percent_read_counts_upper_value[indexNegativeMeans] = scpGroup_read_counts$SD_error_percent_read_counts[indexNegativeMeans];

        #scpGroup_read_counts$SD_percent_read_counts[indexNegativeMeans] = -scpGroup_read_counts$SD_percent_read_counts[indexNegativeMeans]

        sd_errors = scpGroup_read_counts$SD_error_percent_read_counts_upper_value
        sd_errors[is.na(sd_errors)] = 0

        if (length(indexPostiveMeans)>0)
        { max_read_counts = max(c(100,round(10*scpGroup_read_counts$SD_error_percent_read_counts_upper_value[indexPostiveMeans])/10)) }
        if (length(indexPostiveMeans)==0) { max_read_counts = 0 }

        #if (length(indexNegativeMeans)>0)
        #{ min_read_counts = means[indexNegativeMeans]+sds[indexNegativeMeans]*1.1 }
        #if (length(indexNegativeMeans)==0) { min_read_counts = 0}

        if (length(indexNegativeMeans)>0)
        { min_read_counts = means[indexNegativeMeans]*1.1 }
        if (length(indexNegativeMeans)==0) { min_read_counts = 0}

        if (exists("segment_fixed_max_y")) { max_read_counts = segment_fixed_max_y }
        if (exists("segment_fixed_min_y")) { min_read_counts = segment_fixed_min_y }

        Barplot = ggplot(scpGroup_read_counts,aes(x=Physio_segment_dataType,y=Mean_percent_read_counts,fill=Physio_segment_dataType))
        Barplot = Barplot + geom_bar(stat="identity",color=Frame_colors,size=1)
        Barplot = Barplot + geom_errorbar(aes(ymin=Mean_percent_read_counts,ymax=Mean_percent_read_counts+SD_error_percent_read_counts),size=1,width=0.3,linetype=1,color=Frame_colors)
        Barplot = Barplot + ggtitle(plot_title)
        Barplot = Barplot + theme(legend.position="none")
        Barplot = Barplot + theme(axis.text.y = element_text(size=10))
        Barplot = Barplot + scale_fill_manual(values = Colors)
        Barplot = Barplot + xlab("") + ylab("molecule L2B [%]\nUMI counts [%]")
        Barplot = Barplot + theme(axis.title.y = element_text(size=17,hjust=0.5))
        Barplot = Barplot + theme(axis.text.y = element_text(size=17,hjust=0.5,face=2))
        Barplot = Barplot + theme(axis.text.x = element_text(angle=90,size=7,hjust=1,vjust=0.5,face=2))#,color=text_colors))
        Barplot = Barplot + theme(title = element_text(size=10,face=2,hjust=0.5))
        Barplot = Barplot + scale_y_continuous(limits=c(min_read_counts,max_read_counts),breaks=c(-40,-20,0,20,40,60,80,100),oob = rescale_none)
        Barplots[[length(Barplots)+1]] = Barplot
     }#End
  }#End

  considered_centers_count = max(segment_scpGroup_average_read_counts$Considered_centers_count)
  if (length(considered_centers_count)==1) { n_string = considered_centers_count }
  else { n_string = paste(considered_centers_count,collapse=";") }

  max_plots_per_figure = 1
  figures_count = ceiling(length(Barplots)/max_plots_per_figure)
  for (indexF in 1:figures_count)
  {#Begin
    startPlot = (indexF-1)*max_plots_per_figure+1
    endPlot = min(indexF*max_plots_per_figure,length(Barplots))
    current_plots = Barplots[startPlot:endPlot]
    cols_count = min(2,length(current_plots))
    rows_count = ceiling(length(current_plots)/cols_count)
    complete_pdf_fileName = paste(results_directory,"R_meanSdReads_per_segment_no",indexF,"_",n_string,"ns",add_to_fileNames,".pdf",sep='')
    pdf(complete_pdf_fileName,width=pdf_width,height=pdf_height);
    do.call("grid.arrange",c(current_plots,nrow=rows_count,ncol=cols_count))
    dev.off()
  }#End
}#End - Visualize mean and sd of segment percent read counts

##################### Segment specific net L2B
###############################################################################
##################### Segment specific net L2B
###############################################################################
##################### Cell type specific scp distribution

{#Begin - Generate dataset read counts collapsed on each Cell_type/Center/SCP/Data_type
  hr = read_counts_cell_subtypes
  indexKeep = which(hr$Data_type=="1_Read counts")
  hr = hr[indexKeep,]

  cellType_collapsed_read_counts = c()
  hr$Combine = paste(hr$Cell_type,hr$Center,hr$SCP,hr$Data_type,sep="_")
  combines = unique(hr$Combine)
  for (indexC in 1:length(combines))
  {#Begin
     combine = combines[indexC]
     indexCurrentCombine = which(hr$Combine==combine)
     current_hr = hr[indexCurrentCombine,]
     current_hr$Read_counts[1] = sum(current_hr$Read_counts)
     current_hr$Cell_type_counts[1] = sum(current_hr$Cell_subtype_counts)
     if (length(cellType_collapsed_read_counts)==0) { cellType_collapsed_read_counts = current_hr[1,] }
     else { cellType_collapsed_read_counts = rbind(cellType_collapsed_read_counts,current_hr[1,]) }
  }#End
  cellType_collapsed_read_counts$Combine = NULL
  cellType_collapsed_read_counts$Cell_subtype = NULL
  cellType_collapsed_read_counts$Cell_subtype_counts = NULL
  cellType_collapsed_read_counts$Average_read_counts = NULL
}#End - Generate dataset read counts collapsed on each Cell_type/Center/SCP/Data_type

if ((algorithm_of_interest == "Total_sum") & (calculate_percent_read_counts_compared_to_reference_scp))
{#Begin - Calculate percent read counts for all scps within each scp_group over CELL TYPE
  cellType_percent_read_counts = c()
  merged_scp_groups = c(selected_SCPs_for_barplot,selected_SCPs_for_scp_subgroups)
  scp_groups = names(selected_SCPs_for_scp_subgroups)
  scp_groups = scp_groups[order(scp_groups)]
  for (indexScpGroup in 1:length(scp_groups))
  {#Begin
    scp_group = scp_groups[indexScpGroup]
    splitStrings = strsplit(scp_group," ")[[1]]
    ion = splitStrings[1]
    if (length(splitStrings)-2>1)
    {#Begin
      for (indexSplit in 2:(length(splitStrings)-2))
      {#Begin
        ion = paste(ion," ",splitStrings[indexSplit],sep='')
      }#End
    }#End
    reference_scp = paste(ion,net_L2B_label,sep='')
    indexReferenceScp = which(cellType_collapsed_read_counts$SCP==reference_scp)
    ref_sr = cellType_collapsed_read_counts[indexReferenceScp,]
    scps = selected_SCPs_for_scp_subgroups[[scp_group]]
    indexCurrentGroup = which(cellType_collapsed_read_counts$SCP %in% scps)
    if (length(indexCurrentGroup)>0)
    {#Begin
      sr = cellType_collapsed_read_counts[indexCurrentGroup,]
      sr$Combine = paste(sr$Center,sr$Data_type,sep='_')#sr$Cell_type,
      ref_sr$Combine = paste(ref_sr$Center,ref_sr$Data_type,sep='_')#sr$Cell_type,
      combines = unique(sr$Combine)
      for (indexC in 1:length(combines))
      {#Begin
        combine = combines[indexC]
        indexCurrentCombine = which(sr$Combine==combine)
        current_sr = sr[indexCurrentCombine,]
        indexCurrentCombine_reference = which(ref_sr$Combine==combine)
        if (length(indexCurrentCombine_reference)>0)
        {#Begin
          reference_sum_read_counts = abs(sum(ref_sr$Read_counts[indexCurrentCombine_reference]))
          current_sr$Reference_for_percent = reference_scp
        }#End
        if (length(indexCurrentCombine_reference)==0)
        {#Begin
          reference_sum_read_counts = sum(current_sr$Read_counts)
          current_sr$Reference_for_percent = scp_group
        }#End
        if (reference_sum_read_counts>0)
        { current_sr$Percent_read_count = 100*current_sr$Read_counts / reference_sum_read_counts }
        else { current_sr$Percent_read_count = 0 }
        current_sr$SCP_group = scp_group
        if (length(cellType_percent_read_counts)==0) { cellType_percent_read_counts=current_sr}
        else { cellType_percent_read_counts = rbind(cellType_percent_read_counts,current_sr) }
      }#End
    }#End
  }#End
  cellType_percent_read_counts$Read_counts=NULL
  cellType_percent_read_counts$Average_read_counts=NULL
  cellType_percent_read_counts$Cum_read_counts=NULL
  cellType_percent_read_counts$Combine=NULL
}#End - Calculate percent read counts for all scps within each scp_group over CELL TYPE

visualize_B2L_as_negative_L2B = TRUE
if ((algorithm_of_interest == "Total_sum")&(visualize_B2L_as_negative_L2B)&(calculate_percent_read_counts_compared_to_reference_scp))
{#Begin - Set B2L read_counts and percent_read_counts_as_negative
  indexB2LSCPGroups = grep("B2L",cellType_percent_read_counts$SCP_group)
  indexL2BSCPGroups = grep("L2B",cellType_percent_read_counts$SCP_group)
  if (length(indexB2LSCPGroups)>0)
  {#Begin
    cellType_percent_read_counts$Percent_read_count[indexB2LSCPGroups] = -cellType_percent_read_counts$Percent_read_count[indexB2LSCPGroups]
    cellType_percent_read_counts$SCP_group[indexB2LSCPGroups] = gsub(" B2L "," L2B ",cellType_percent_read_counts$SCP_group[indexB2LSCPGroups])
  }#End
  #cellType_percent_read_counts$SCP[indexB2LSCPGroups] = gsub(" B2L "," L2B ",cellType_percent_read_counts$SCP[indexB2LSCPGroups])
}#End - Set B2L read_counts and percent_read_counts_as_negative

{#Begin - Calculate mean and sd error of cell type read counts for each SCP in each Scp_group
  stP = cellType_percent_read_counts
  keep_centers = keep_centers_for_each_analysis

  indexKeep = which(stP$Center %in% keep_centers)
  stP = stP[indexKeep,]
  if (length(which(!keep_centers %in% stP$Center))>0) { stP = c()}
  length_centers = length(keep_centers)

  stP$Combine = paste(stP$Cell_type,stP$SCP_group,stP$SCP,sep="_")
  combines = unique(stP$Combine)
  cellType_mean_percent_read_counts = c()
  double_scps = c()
  combine_with_double_scps = c()
  indexCombine = 1
  for (indexCombine in 1:length(combines))
  {#Begin
     current_combine = combines[indexCombine]
     indexCurrentCombine = which(stP$Combine==current_combine)
     current_stP = stP[indexCurrentCombine,]
     current_scps = unique(current_stP$SCP)
     if (length(unique(sign(current_stP$Percent_read_count[current_stP$Percent_read_count!=0])))>1)
     {#Begin - Check for double SCPs
        double_scps = c(double_scps,current_scps)
        combine_with_double_scps = c(combine_with_double_scps,current_combine)
     }#End - Check for double SCPs

     current_percent_read_counts = current_stP$Percent_read_count
     if (length(current_percent_read_counts)>length_centers) { current_percent_read_counts = replicate(length_centers, -99999999) }
     while (length(current_percent_read_counts)<length_centers) { current_percent_read_counts = c(current_percent_read_counts,0) }

     mean_current_stP = current_stP[1,]
     mean_current_stP$Mean_percent_read_counts = mean(current_percent_read_counts)
     mean_current_stP$Number_of_averaged_centers = length(current_percent_read_counts)
     if (length(current_stP$Percent_read_count)==1)
     {#Begin
         mean_current_stP$SD_percent_read_counts = 0
         mean_current_stP$SD_error_percent_read_counts = 0
     }#End
     if (length(current_stP$Percent_read_count)>1)
     {#Begin
         mean_current_stP$SD_percent_read_counts = sd(current_percent_read_counts)
         mean_current_stP$SD_error_percent_read_counts = mean_current_stP$SD_percent_read_counts / sqrt(length(current_percent_read_counts))
     }#End
     if (length(cellType_mean_percent_read_counts)==0) { cellType_mean_percent_read_counts = mean_current_stP }
     else { cellType_mean_percent_read_counts = rbind(cellType_mean_percent_read_counts,mean_current_stP) }
  }#End

  double_scps = unique(double_scps)
  if (length(double_scps)>0) { cellType_mean_percent_read_counts = c(); if (error_message=="") { error_message = "Calculate mean and sd of cell type read counts for each SCP in each Scp_group: double assignments" }}

  if (length(unique(cellType_mean_percent_read_counts$Number_of_averaged_centers))!=1) { cellType_mean_percent_read_counts = c() }
  if (unique(cellType_mean_percent_read_counts$Number_of_averaged_centers) != length_centers) { cellType_mean_percent_read_counts = c() }

  cellType_mean_percent_read_counts$SD_percent_read_counts=NULL
  cellType_mean_percent_read_counts$Combine=NULL
  cellType_mean_percent_read_counts$Center=NULL
  cellType_mean_percent_read_counts$Percent_cum_read_counts=NULL
  cellType_mean_percent_read_counts$Read_counts = NULL
  cellType_mean_percent_read_counts$Keep_even_if_no_read_counts = FALSE
}#End - Calculate mean and sd error of cell type read counts for each SCP in each Scp_group

cellType_mean_percent_read_counts$SCP_original = cellType_mean_percent_read_counts$SCP

if (algorithm_of_interest == "Total_sum")
{#Begin - Add trimmed hierarchy organization to cellType_mean_percent_read_counts
  Add_hierarchy_to_readCounts_recursive = function(scpGroup_mean_read_counts,trimmed_scp_hierarchy,current_parent,consecutive_no,add_to_name,subscp_and_category_label,characters_to_add)
  {#Begin
     add_to_name = paste(characters_to_add,add_to_name,sep='')
     parent_has_children_in_read_counts = FALSE
     indexCurrentHierarchy = which(trimmed_scp_hierarchy$Parent_scp==current_parent)
     if (length(indexCurrentHierarchy)>0)
     {#Begin
        current_children = unique(trimmed_scp_hierarchy$Child_scp[indexCurrentHierarchy])
        current_children = current_children[order(current_children)]
        indexChild=3
        for (indexChild in 1:length(current_children))
        {#Begin
           current_child = current_children[indexChild]
           indexCurrentReadCounts = which(scpGroup_mean_read_counts$SCP==paste(current_child,subscp_and_category_label,sep=''))
           if (length(indexCurrentReadCounts)>0)
           {#Begin
              parent_has_children_in_read_counts = TRUE
           }#End
           return_list = Add_hierarchy_to_readCounts_recursive(scpGroup_mean_read_counts,trimmed_scp_hierarchy,current_child,consecutive_no,add_to_name,subscp_and_category_label,characters_to_add)
           scpGroup_mean_read_counts = return_list[[1]]
           consecutive_no = return_list[[2]]
           parent_has_children_in_read_counts_internal = return_list[[3]]
           if (parent_has_children_in_read_counts_internal) { parent_has_children_in_read_counts = TRUE }
        }#End
     }#End
     indexParent_in_readCounts = which(scpGroup_mean_read_counts$SCP==paste(current_parent,subscp_and_category_label,sep=''))
     if (length(indexParent_in_readCounts)>0)
     {#Begin
       if (parent_has_children_in_read_counts)
       {#Begin
          scpGroup_mean_read_counts$Add_in_front_of_SCP[indexParent_in_readCounts] = add_to_name
       }#End
       if (!parent_has_children_in_read_counts)
       {#Begin
          scpGroup_mean_read_counts$Add_in_front_of_SCP[indexParent_in_readCounts] = add_to_name
       }#End
       consecutive_no = consecutive_no - 1
       scpGroup_mean_read_counts$Consecutive_no[indexParent_in_readCounts] = consecutive_no
     }#End
     if (parent_has_children_in_read_counts)
     {#Begin
       consecutive_no = consecutive_no - 1
       if (length(indexParent_in_readCounts)>0)
       {#Begin
          scpGroup_mean_read_counts$Keep_even_if_no_read_counts[indexParent_in_readCounts] = TRUE
       }#End
       if (length(indexParent_in_readCounts)==0)
       {#Begin
           add_read_counts = scpGroup_mean_read_counts[1,]
           add_read_counts$SCP = current_parent
           add_read_counts$Add_in_front_of_SCP = add_to_name
           add_read_counts$Consecutive_no = 100
           add_read_counts$Percent_read_count = 0
           add_read_counts$Mean_percent_read_counts = 0;
           add_read_counts$SD_error_percent_read_counts = 0
           add_read_counts$Consecutive_no = consecutive_no
           add_read_counts$Keep_even_if_no_read_counts = TRUE
           scpGroup_mean_read_counts = rbind(scpGroup_mean_read_counts,add_read_counts)
       }#End
     }#End
     return_list = list(scpGroup_mean_read_counts,consecutive_no,parent_has_children_in_read_counts)
     return (return_list)
  }#End

  scp_groups = unique(cellType_mean_percent_read_counts$SCP_group)
  cellType_mean_percent_read_counts$Consecutive_no = 0
  indexScpGroup = 1
  new_cellType_mean_percent_read_counts = c()
  cellType_mean_percent_read_counts$Add_in_front_of_SCP = ""
  indexScpGroup = 1
  characters_to_add = "   "
  for (indexScpGroup in 1:length(scp_groups))
  {#Begin
     scp_group = scp_groups[indexScpGroup]
     splitStrings = str_split(scp_group," ")[[1]]
     ion = splitStrings[1]
     length_ion_name = length(splitStrings)-2
     if (length_ion_name>1)
     {#Begin
        for (indexI in 2:length_ion_name)
        {#Begin
           ion = paste(ion,splitStrings[indexI],sep=' ')
        }#End
     }#End
     indexCurrentGroup = which(cellType_mean_percent_read_counts$SCP_group==scp_group)
     scpGroup_mean_read_counts = cellType_mean_percent_read_counts[indexCurrentGroup,]
     all_scps = unique(scpGroup_mean_read_counts$SCP)
     final_parent_scps = c(paste(ion," L2B",sep=''),paste(ion," B2L",sep=''),paste(ion," channel",sep=''))
     indexFinalParent = 1
     for (indexFinalParent in 1:length(final_parent_scps))
     {#Begin
        final_parent = final_parent_scps[indexFinalParent]
        {#Begin - Generate trimmed hierarchy
           indexCurrentChildren = which(scp_hierarchy$Parent_scp==final_parent)
           trimmed_scp_hierarchy = scp_hierarchy[indexCurrentChildren,]
           while (length(indexCurrentChildren)>0)
           {#Begin
              next_parents = scp_hierarchy$Child_scp[indexCurrentChildren]
              indexCurrentChildren = which(scp_hierarchy$Parent_scp %in% next_parents)
              trimmed_scp_hierarchy = rbind(trimmed_scp_hierarchy,scp_hierarchy[indexCurrentChildren,])
           }#End
           trimmed_scp_hierarchy = trimmed_scp_hierarchy[order(trimmed_scp_hierarchy$Child_scp,trimmed_scp_hierarchy$Parent_priotirization),]
           indexUnique = which(!duplicated(trimmed_scp_hierarchy$Child_scp,sep=''))
           trimmed_scp_hierarchy = trimmed_scp_hierarchy[indexUnique,]
        }#End - Generate trimmed hierarchy

        consecutive_no = (indexFinalParent) * 100
        subscp_and_category_label = paste(subscp_label,final_parent," M",sep='')
        add_to_name = ""
        return_list = Add_hierarchy_to_readCounts_recursive(scpGroup_mean_read_counts, trimmed_scp_hierarchy, final_parent, consecutive_no, add_to_name, subscp_and_category_label,characters_to_add)
        scpGroup_mean_read_counts = return_list[[1]]
     }#End
     if (length(new_cellType_mean_percent_read_counts)==0) { new_cellType_mean_percent_read_counts = scpGroup_mean_read_counts }
     else { new_cellType_mean_percent_read_counts = rbind(new_cellType_mean_percent_read_counts,scpGroup_mean_read_counts) }
  }#End
  for (indexNew in 1:length(new_cellType_mean_percent_read_counts[,1]))
  {#Begin
     new_cellType_mean_percent_read_counts$Add_in_front_of_SCP[indexNew] = substr(new_cellType_mean_percent_read_counts$Add_in_front_of_SCP[indexNew],nchar(characters_to_add)+1,nchar(new_cellType_mean_percent_read_counts$Add_in_front_of_SCP[indexNew]))
  }#End
  new_cellType_mean_percent_read_counts$SCP = paste(new_cellType_mean_percent_read_counts$Add_in_front_of_SCP,new_cellType_mean_percent_read_counts$SCP,sep='')
  cellType_mean_percent_read_counts = new_cellType_mean_percent_read_counts
}#End - Add trimmed hierarchy organization to cellType_mean_percent_read_counts

indexEmpty = which(cellType_mean_percent_read_counts$SCP=="")

Check_for_matching_read_counts=TRUE
if (Check_for_matching_read_counts)
{#Begin - Check if cellType means agree with segment means
  if (library_of_interest=="NaAndGluTMTransport")
  {#Begin
    cellTypeScpGroup_map_to_segmentScp = list("Sodium L2B M" = "Sodium net L2B")
  }#End

   indexReadCounts = which(segment_scpGroup_average_read_counts$Data_type=="1_Read counts")
   segment_readCounts_scpGroup_average_read_counts = segment_scpGroup_average_read_counts[indexReadCounts,]
   for (indexMap in 1:length(cellTypeScpGroup_map_to_segmentScp))
   {#Begin
      cellTypeScpGroup = names(cellTypeScpGroup_map_to_segmentScp)[indexMap]
      segmentScpGroup = cellTypeScpGroup_map_to_segmentScp[[cellTypeScpGroup]]
      indexCurrentScpGroup_cellType = which(cellType_mean_percent_read_counts$SCP_group==cellTypeScpGroup)
      indexCurrentScpGroup_segment = which(segment_readCounts_scpGroup_average_read_counts$SCP_group==segmentScpGroup)
      currentScpGroup_cellType_mean_percent_read_counts = cellType_mean_percent_read_counts[indexCurrentScpGroup_cellType,]
      currentScpGroup_segment_mean_percent_read_counts = segment_readCounts_scpGroup_average_read_counts[indexCurrentScpGroup_segment,]
      cellType_physio_segments = unique(currentScpGroup_cellType_mean_percent_read_counts$Physio_segment)
      segment_physio_segments = unique(currentScpGroup_segment_mean_percent_read_counts$Physio_segment)
      if (length(which(!cellType_physio_segments %in% segment_physio_segments))>0) { cellType_mean_percent_read_counts = c() }
      if (length(which(!segment_physio_segments %in% cellType_physio_segments))>0) { cellType_mean_percent_read_counts = c() }
      indexPhysio=3
      for (indexPhysio in 1:length(cellType_physio_segments))
      {#Begin
         cellType_physio_segment = cellType_physio_segments[indexPhysio]
         indexCurrentPhysio_cellType = which(currentScpGroup_cellType_mean_percent_read_counts$Physio_segment==cellType_physio_segment)
         indexCurrentPhysio_segment = which(currentScpGroup_segment_mean_percent_read_counts$Physio_segment==cellType_physio_segment)
         if (round(100000*sum(currentScpGroup_segment_mean_percent_read_counts$Mean_percent_read_counts[indexCurrentPhysio_segment]))!=round(100000*sum(currentScpGroup_cellType_mean_percent_read_counts$Mean_percent_read_counts[indexCurrentPhysio_cellType]))) {  if (error_message=="") { error_message = paste("mean read count mismatch: ",cellType_physio_segment," ",cellTypeScpGroup,sep='');cellType_mean_percent_read_counts = c() }}
      }#End
   }#End
}#End - Check if cellType means agree with segment means

cellTypeMean_fileName = "Read_counts_mean_cellType.txt"
complete_cellTypeMean_fileName = paste(results_directory,cellTypeMean_fileName,sep='')
write.table(cellType_mean_percent_read_counts,file=complete_cellTypeMean_fileName,row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)

#cellType_fixed_max_y = 90
#cellType_fixed_min_y = 0

legend_columns_count=1
{#Begin - Visualize mean and sd error of cell type percent read counts
  legend_key_size_in_cm = 0.18
  pdf_width = 8.27   #8.27 = DIN A4
  pdf_height = 5 #11.69 = DIN A4
  legend_key_text_size = 5.45

  scp_groups = unique(cellType_mean_percent_read_counts$SCP_group)
  Barplots = list()

  defined_scpGroups = names(selected_SCPs_for_scp_subgroups)
  indexScpGroup=2
  for (indexScpGroup in 1:length(scp_groups))
  {#Begin
     scp_group = scp_groups[indexScpGroup]
     current_ion = strsplit(scp_group," ")[[1]]
     full_L2B_subscp_label = paste(subscp_label,scp_group,sep='')
     full_B2L_subscp_label = gsub("L2B","B2L",full_L2B_subscp_label)
     indexCurrentGroup = which(cellType_mean_percent_read_counts$SCP_group==scp_group)
     scpGroup_read_counts = cellType_mean_percent_read_counts[indexCurrentGroup,]
     indexKeep = unique(c(which(scpGroup_read_counts$Mean_percent_read_counts!=0),which(scpGroup_read_counts$Keep_even_if_no_read_counts==TRUE)))
     if (length(indexKeep)>0)
     {#Begin
       scpGroup_read_counts = scpGroup_read_counts[indexKeep,]
       cellTypes_in_correct_order = c("PT","DTL","ATL","TAL","DCT","CNT","PC","tPC_IC","IC")
       scpGroup_read_counts$Cell_type_factor = factor(scpGroup_read_counts$Cell_type,levels=cellTypes_in_correct_order)

       {#Begin - Remove L2B and B2L labels and adjust secreation SCP names, if also in L2B
          indexB2L = grep(full_B2L_subscp_label,scpGroup_read_counts$SCP)
          indexL2B = grep(full_L2B_subscp_label,scpGroup_read_counts$SCP)
          scpGroup_read_counts$SCP = gsub(full_L2B_subscp_label,"",scpGroup_read_counts$SCP)
          scpGroup_read_counts$SCP = gsub(full_B2L_subscp_label,"",scpGroup_read_counts$SCP)
          indexB2LSCP_in_L2B = which(scpGroup_read_counts$SCP[indexB2L] %in% scpGroup_read_counts$SCP[indexL2B])
          if (length(indexB2LSCP_in_L2B)>0)
          {#Begin
             indexAddSpace = indexB2L[indexB2LSCP_in_L2B]
             scpGroup_read_counts$SCP[indexAddSpace] = paste(scpGroup_read_counts$SCP[indexAddSpace]," ",sep='')
          }#End
       }#End - Remove L2B and B2L labels and adjust secreation SCP names, if also in L2B

       if (library_of_interest=="NaAndGluTMTransport")
       {#Begin - Set plot_scps_in_correct_order
          indexCurrentIon = grep(current_ion,names(selected_SCPs_for_scp_subgroups))
          current_names = names(selected_SCPs_for_scp_subgroups)[indexCurrentIon]
          if (  (length(grep("L2B",current_names))>0)
               &(length(grep("B2L",current_names))>0)
               &(length(grep("channel",scp_group))==0)
               &(length(grep("Sodium potassium ATPase",scp_group))==0))
          {#Begin
             indexL2B = grep("L2B",names(selected_SCPs_for_scp_subgroups))
             indexB2L = grep("B2L",names(selected_SCPs_for_scp_subgroups))
             indexReabSecr = c(indexL2B,indexB2L)
             indexCurrentIon = indexCurrentIon[indexCurrentIon %in% indexReabSecr]
             if (length(indexCurrentIon)==2)
             {#Begin
                add_read_counts = scpGroup_read_counts[1,]
                add_read_counts$SCP = "";
                add_read_counts$Consecutive_no = 100
                add_read_counts$Percent_read_count=0
                add_read_counts$Mean_percent_read_counts=0;
                add_read_counts$SD_error_percent_read_counts=0
                scpGroup_read_counts = rbind(scpGroup_read_counts,add_read_counts)
             }#End
          }#End
          scpGroup_read_counts = scpGroup_read_counts[order(scpGroup_read_counts$Consecutive_no),]
       }#End - Set plot_scps_in_correct_order
       plot_scps_in_correct_order = unique(scpGroup_read_counts$SCP)

       scpGroup_read_counts$SCP_factor = factor(scpGroup_read_counts$SCP,levels=plot_scps_in_correct_order)

       {#Begin - Initialize error bar visualizations
         scpGroup_read_counts$Error_bar_y_position = -1000
         scpGroup_read_counts$Error_bar_x_position = -1000
         scpGroup_read_counts = scpGroup_read_counts[order(scpGroup_read_counts$SCP_factor,decreasing=TRUE),]
         scpGroup_read_counts = scpGroup_read_counts[order(scpGroup_read_counts$Cell_type_factor,decreasing=FALSE),]
         error_bar_shifts = c(replicate(length(unique(scpGroup_read_counts$SCP)),c(0)))
         cellTypes = unique(scpGroup_read_counts$Cell_type)
         for (indexCT in 1:length(cellTypes))
         {#Begin - Set postion of error bars
           cellType = cellTypes[indexCT]
           indexCurrentCelltype = which(scpGroup_read_counts$Cell_type==cellType)
           for (indexSign in 1:2)
           {#Begin
             if (indexSign==1) { indexCurrentSign = which(scpGroup_read_counts$Mean_percent_read_counts>=0) }
             if (indexSign==2) { indexCurrentSign = which(scpGroup_read_counts$Mean_percent_read_counts<0) }
             indexCurrentCelltype_sign = indexCurrentCelltype[indexCurrentCelltype %in% indexCurrentSign]
             scpGroup_read_counts$Error_bar_y_position[indexCurrentCelltype_sign] = cumsum(scpGroup_read_counts$Mean_percent_read_counts[indexCurrentCelltype_sign])
             current_x_shifts = error_bar_shifts[1:length(indexCurrentCelltype_sign)]
             scpGroup_read_counts$Error_bar_x_position[indexCurrentCelltype_sign] = indexCT + current_x_shifts
           }#End
         }#End - Initialize error bar visualizations

         error_bar_width = 0.4
         error_bar_size = 0.3

         indexNegative = which(scpGroup_read_counts$Mean_percent_read_counts<0)
         scpGroup_read_counts$SD_error_percent_read_counts_upper_value = scpGroup_read_counts$Error_bar_y_position + scpGroup_read_counts$SD_error_percent_read_counts
         scpGroup_read_counts$SD_error_percent_read_counts_lower_value = scpGroup_read_counts$Error_bar_y_position
         if (length(indexNegative)>0)
         {#Begin
           scpGroup_read_counts$SD_error_percent_read_counts_upper_value[indexNegative] = scpGroup_read_counts$Error_bar_y_position[indexNegative] - scpGroup_read_counts$Mean_percent_read_counts[indexNegative] + scpGroup_read_counts$SD_error_percent_read_counts[indexNegative]
           scpGroup_read_counts$SD_error_percent_read_counts_lower_value[indexNegative] = scpGroup_read_counts$Error_bar_y_position[indexNegative] - scpGroup_read_counts$Mean_percent_read_counts[indexNegative]
         }#End
       }#End - Initialize error bars

       scpGroup_read_counts = scpGroup_read_counts[order(scpGroup_read_counts$SCP_factor),]
       indexNotZero = which(scpGroup_read_counts$Percent_read_count!=0)
       scps_for_colors = unique(scpGroup_read_counts$SCP[indexNotZero])
       scps_all = unique(scpGroup_read_counts$SCP)

       zero_line_color = "black"#"gray60"
       Colors_for_colored_scps = color_sets[[length(scps_for_colors)]]
       Colors = replicate(length(scps_all),"white")
       for (indexSCP in 1:length(scps_all))
       {#Begin - Set fill colors
         current_scp = scps_all[indexSCP]
         indexCurrentSCP_in_scps_for_colors = which(scps_for_colors==current_scp)
         if (length(indexCurrentSCP_in_scps_for_colors)>0)
         {#Begin
            Frame_colors[indexSCP] = Colors_for_colored_scps[indexCurrentSCP_in_scps_for_colors]
            Colors[indexSCP] = Colors_for_colored_scps[indexCurrentSCP_in_scps_for_colors]
         }#End
       }#End - Set fill colors

       {#Begin - Set frame colors
          scpGroup_read_counts = scpGroup_read_counts[order(scpGroup_read_counts$SCP_factor,decreasing=FALSE),]
          scpGroup_read_counts = scpGroup_read_counts[order(scpGroup_read_counts$Cell_type_factor,decreasing=FALSE),]

          indexNotZero = which(scpGroup_read_counts$Percent_read_count!=0)
          Frame_colors = replicate(length(scpGroup_read_counts[,1]),zero_line_color)
          for (indexIndex in 1:length(indexNotZero))
          {#Begin
             current_scp = scpGroup_read_counts$SCP[indexNotZero[indexIndex]]
             indexColor = which(scps_for_colors==current_scp)
             if (length(indexColor)!=0)
             {#Begin
                Frame_colors[indexNotZero[indexIndex]] = Colors_for_colored_scps[indexColor]
             }#End
          }#End
       }#End - Set frame colors

       {#Begin - Set y_breaks
          if (!scp_group %in% names(scpGroupSpecifi_ylims))
          {#Begin - Identify max positive, negative and total heights
            cell_types = unique(scpGroup_read_counts$Cell_type)
            max_positive_height = 0
            max_negative_height = 0
            max_total_height = 0
            indexPositive = which(scpGroup_read_counts$Mean_percent_read_counts>0)
            indexNegative = which(scpGroup_read_counts$Mean_percent_read_counts<0)
            for (indexCellType in 1:length(cell_types))
            {#Begin
              cell_type = cell_types[indexCellType]
              indexCurrentCellType = which(scpGroup_read_counts$Cell_type==cell_type)
              indexCurrentPostive = indexCurrentCellType[indexCurrentCellType %in% indexPositive]
              indexCurrentNegative = indexCurrentCellType[indexCurrentCellType %in% indexNegative]
              current_postive_height = sum(abs(scpGroup_read_counts$Mean_percent_read_counts[indexCurrentPostive]))
              current_negative_height = sum(abs(scpGroup_read_counts$Mean_percent_read_counts[indexCurrentNegative]))
              current_total_height = sum(abs(scpGroup_read_counts$Mean_percent_read_counts[indexCurrentCellType]))
              if (current_postive_height>max_positive_height) { max_positive_height = current_postive_height }
              if (current_negative_height>max_negative_height) { max_negative_height = current_negative_height }
              if (current_total_height>max_total_height) { max_total_height = current_total_height }
            }#End
            max_total_height = round(max_total_height * 1.2)
            max_positive_height = round(max_positive_height * 1.2)
            max_negative_height = round(max_negative_height * 1.2)
          }#End - Identify max positive, negative and total heights
          if (scp_group %in% names(scpGroupSpecifi_ylims))
          {#Begin
            current_ylims = scpGroupSpecifi_ylims[[scp_group]]
            max_positive_height = current_ylims[["Upper_limit"]]
            max_negative_height = abs(current_ylims[["Lower_limit"]])
            max_total_height = max_positive_height + max_negative_height
          }#End

          number_of_marks = 5

          height_per_mark = max_total_height / number_of_marks
          considered_height_distances = c(2.5,5,10,20,25)
          for (indexCons in 1:length(considered_height_distances))
          {#Begin
             considered_height_distance = considered_height_distances[indexCons]
             if ((round(height_per_mark/considered_height_distance)) <= 1)
             {#Begin
                height_per_mark = considered_height_distance
                break;
             }#End
          }#End

          positive_marks = ceiling(max_positive_height / height_per_mark)
          positive_marks = positive_marks[positive_marks !=0]
          negative_marks = ceiling(max_negative_height / height_per_mark)
          negative_marks = negative_marks[negative_marks !=0]
          y_breaks = c()
          if (length(negative_marks)>0)
          { y_breaks = c(y_breaks,-(1:negative_marks)*height_per_mark) }
          y_breaks = c(y_breaks,0)
          if (length(positive_marks)>0)
          { y_breaks = c(y_breaks,(1:positive_marks)*height_per_mark) }
       }#End - Set y_breaks

       plot_title = scp_group
       plot_title = gsub("L2B","transporter",plot_title)
       plot_title = gsub(" M","",plot_title)
       max_sd_error = max(scpGroup_mean_read_counts$SD_error_percent_read_counts)

       max_y = max(scpGroup_read_counts$SD_error_percent_read_counts_upper_value)*1.1
       min_y = min(0, scpGroup_read_counts$Error_bar_y_position)*1.1

       if (exists("cellType_fixed_max_y")) { max_y = cellType_fixed_max_y }
       if (exists("cellType_fixed_min_y")) { min_y = cellType_fixed_min_y }

       Barplot = ggplot(scpGroup_read_counts,aes(x=Cell_type_factor,y=Mean_percent_read_counts,fill=SCP_factor,color=Colors))
       if (max_sd_error==0)
       {#Begin
          Barplot = Barplot + geom_bar(stat="identity",size=0)
       }#End
       if (max_sd_error > 0)
       {#Begin
          Barplot = Barplot + geom_bar(stat="identity",color=Frame_colors,size=error_bar_size)
          Barplot = Barplot + geom_errorbar(aes(x=Error_bar_x_position,ymin=SD_error_percent_read_counts_lower_value,ymax=SD_error_percent_read_counts_upper_value),size=error_bar_size,width=error_bar_width,linetype=1,color=Frame_colors)
       }#End
       if (scp_group %in% names(scpGroupSpecifi_ylims))
       {#Begin
          Barplot = Barplot + scale_y_continuous(breaks=y_breaks,limits=c(-max_negative_height,max_positive_height))
       }#End
       if (!scp_group %in% names(scpGroupSpecifi_ylims))
       {#Begin
          Barplot = Barplot + scale_y_continuous(breaks=y_breaks)
       }#End
       #Barplot = Barplot + geom_errorbar(aes(x=Error_bar_x_position,ymin=Error_bar_y_position-SD_percent_read_counts,ymax=Error_bar_y_position+SD_percent_read_counts),size=1,width=0.05,linetype=1,color="black")#,color=Frame_colors
       Barplot = Barplot + ggtitle(plot_title)
       Barplot = Barplot + scale_fill_manual(values = Colors)
       Barplot = Barplot + guides(fill=guide_legend(ncol=legend_columns_count,))
       Barplot = Barplot + theme(legend.key.size = unit(legend_key_size_in_cm, 'cm'))
       Barplot = Barplot + theme(legend.key = element_rect(fill="white",color=NA))
       #Barplot = Barplot + theme(legend.key = element_rect())
       Barplot = Barplot + geom_hline(yintercept=0,color=zero_line_color,size=1.5)
       ylabel = "error"
       ylabel = "mRNA levels [%]"
       Barplot = Barplot + xlab("") + ylab(ylabel)
       Barplot = Barplot + theme(axis.title.y = element_text(size=12,color="black",face=2,hjust=0.5))
       Barplot = Barplot + theme(axis.text.y = element_text(size=12,color="black",face=2,hjust=0.5))
       Barplot = Barplot + theme(axis.text.x = element_text(angle=0,size=12,color="black",hjust=0.5,vjust=0.5,face=2))#,color=text_colors))
       Barplot = Barplot + theme(legend.text = element_text(size=legend_key_text_size,color="black",face=2))
       Barplot = Barplot + theme(legend.title = element_blank())
       Barplot = Barplot + scale_y_continuous(limits=c(min_y,max_y),breaks=c(-40,-20,0,20,40,60,80,100),oob = rescale_none)
       #Barplot = Barplot + ylim(c(-40,100))

       Barplot = Barplot + theme(plot.title = element_text(size=14,face=2,hjust=0.5))
       Legend = get_legend(Barplot)
       Barplot = Barplot + theme(legend.position = "none")

       Barplots[[length(Barplots)+1]] = Barplot
       Barplots[[length(Barplots)+1]] = Legend
     }#End
     else
     {#Begin
       Barplots[[length(Barplots)+1]] = ggplot() + theme_void()
     }#End
  }#End

  max_plots_per_figure = 2
  figures_count = ceiling(length(Barplots)/max_plots_per_figure)
  for (indexF in 1:figures_count)
  {#Begin
    startPlot = (indexF-1)*max_plots_per_figure+1
    endPlot = min(indexF*max_plots_per_figure,length(Barplots))
    current_plots = Barplots[startPlot:endPlot]
    cols_count = min(2,length(current_plots))
    rows_count = ceiling(length(current_plots)/cols_count)
    complete_pdf_fileName = paste(results_directory,"R_meanSdReads_per_cellType_no",indexF,"_",n_string,"ns",add_to_fileNames,".pdf",sep='')
    pdf(complete_pdf_fileName,width=pdf_width,height=pdf_height);
    do.call("grid.arrange",c(current_plots,nrow=rows_count,ncol=cols_count))
    dev.off()
  }#End
}#End - Visualize mean and sd error of cell type percent read counts

##################### Cell type specific scp distribution
###############################################################################
###############################################################################

}#End - Calculate reabsorption profiles for selected datasets


# Change back to start-up directory
if(!is.null(start_dir)) setwd(start_dir)
