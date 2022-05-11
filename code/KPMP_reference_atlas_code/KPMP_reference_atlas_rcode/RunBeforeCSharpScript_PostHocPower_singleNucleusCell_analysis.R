#Copyright 2022.The code was written by Jens Hansen (if not stated otherwise) working for the Ravi Iyengar Lab
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


{#Begin - Extract command line input arguments - written by Yuguang Xiong
   # Determine script invocation method.
   getInvokeMethod <- function()
   {
      # Get all command arguments including R system arguments.
      cmd_args <- commandArgs(FALSE)

      # Get the value of "--file" option.
      matched_pattern <- regexpr("(?<=^--file=).+", cmd_args, perl=TRUE)
      prog_cmd <- regmatches(cmd_args, matched_pattern)
      n_prog_cmd <- length(prog_cmd)

      # Get the availability of "--args" option.
      args_opt_avail <- "--args" %in% cmd_args

      # Determine invocation method based on n_prog_cmd and args_opt_avail.
      invoke_method <- NULL
      if(n_prog_cmd == 0) invoke_method <- "R"
      else if(n_prog_cmd == 1) invoke_method <- "Rscript"
      else invoke_method <- "unknown"

      # Return invocation method.
      return(invoke_method)
   }

   # Obtain the path of this program in Rscript command line.
   getProgPath <- function(invoke_method)
   {
      prog_path <- NULL

      # Determine the function/config path.
      if(invoke_method == "Rscript")
      {
         # Retrieve function path from the path of this main program.
         cmd_args <- commandArgs(FALSE)
         # Get the value of "--file" option.
         matched_pattern <- regexpr("(?<=^--file=).+", cmd_args, perl=TRUE)
         prog_path <- regmatches(cmd_args, matched_pattern)
      }

      # Return the function/config path.
      return(prog_path)
   }

   # Check the command line of script invocation by Rscript.
   checkInvokeCommandLine <- function(invoke_method, n_args_min)
   {
      # Initialize checking status.
      cmd_status <- TRUE

      # Only check command line for Rscript invocation.
      if(invoke_method == "Rscript")
      {
         # Get the program name.
         prog_name <- basename(getProgPath(invoke_method))
         # Get the command line of script invocation by Rscript.
         cmd_args <- commandArgs(TRUE)
         if(length(cmd_args) < n_args_min)
         {
            cat(paste("Usage:", prog_name, "[tis_center] [cores_count]\n"))
            cmd_status <- FALSE
         }
      }
      else if(invoke_method == "R")
      {
         # Set cmd_status to TRUE to allow running this script in an interactive R console.
         cmd_status <- TRUE
      }
      else
      {
         warning("This script needs to be launched in either interactive R console or Rscript command line!")
         cmd_status <- FALSE
      }

      # Return checking status.
      return(cmd_status)
   }

   # Obtain the value of a specified input arguments in Rscript command line or use the default.
   getArg <- function(invoke_method, arg_name, arg_pos, default_value=NULL, verbose=FALSE)
   {
      # invoke_method: "R" or "Rscript".

      arg <- NULL

      # Determine the function/config path.
      if(invoke_method == "Rscript")
      {
         # Type 2: for the invocation by Rscript at system terminal.
         cmd_args <- commandArgs(TRUE)
         if(length(cmd_args) > arg_pos-1) arg <- cmd_args[arg_pos]
         else
         {
            if(verbose) warning(paste("No", arg_name, "is provided and use the default value!"))
            arg <- default_value
         }
         if(is.null(arg) && verbose) warning(paste(arg_name, "is NULL when invoke_method is \"Rscript\"!"))
      }
      else if(invoke_method == "R")
      {
         # Type 1: for the invocation by "R" at R terminal.
         arg <- default_value
         if(is.null(arg) && verbose) warning(paste(arg_name, "is NULL when invoke_method is \"R\"!"))
      }
      else
      {
         if(verbose) warning("invoke_method must be either \"R\" or \"Rscript\"!")
      }

      # Return the function/config path.
      return(arg)
   }

   ################### Main Program ###################

   # Determine script invocation method.
   invoke_method <- getInvokeMethod()
   if(invoke_method == "unknown") stop("Script must be invoked by either R or Rscript!")

   # Obtain invoked R script command.
   stopifnot(checkInvokeCommandLine(invoke_method,1))

   # Determine the first input argument.
   tis_center <- getArg(invoke_method=invoke_method, arg_name="tis_center", arg_pos=1, default_value="Premiere")
   if(tis_center!="Premiere" && tis_center!="UCSD") stop("tis_center must be either Premiere or UCSD!")
   cores_count <- getArg(invoke_method=invoke_method, arg_name="cores_count", arg_pos=2, default_value=1)
   cores_count <- as.integer(cores_count)
   if(!is.na(cores_count)) {if(cores_count<=0) stop("cores_count must be a positive integer!")}
   else stop("cores_count must be a numeric number!")
}#End - Extract command line input arguments - written by Yuguang Xiong

tis_center = "UCSD"
tis_center = "Premiere"
base_directory = paste(getwd(),"/",sep='');
base_directory = "D:/KPMP_reference_atlas_code/"

set = "2019June_revisedScAdv_2022January31"

{#Begin - Open libraries and document versions - BEGIN
   Col_names = c("Library","Version")
   Col_length = length(Col_names)
   Row_names = 1
   Row_length= length(Row_names)
   version_documentation_line = array(NA,c(Row_length,Col_length),dimnames=list(Row_names,Col_names))
   version_documentation_line = as.data.frame(version_documentation_line)
   version_documentations = c()

   libraries = c("SingleCellExperiment","Seurat","ggplot2","sctransform","mclust","dplyr","beeswarm","clustree","Matrix","gridExtra","reticulate","umap","doParallel")
   #use_python("c://ProgramData//Anaconda3//")
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
   #python_sessionInfo = py_config()
}#End - Open libraries and document versions - END

{#Begin - Specify file names and generate directories
    directory = paste(base_directory,"Experimental_data/SingleCell_datasets/",sep='');
    metadata_directory = paste(base_directory,"/Experimental_data/Sample_metadata_additional_datasets/",sep='');

    if (tis_center=="UCSD")
    {#Begin
       rawData_fileName = paste("GSE121862_UCSD-WU_Single_Nuclei_Cluster_Annotated_Raw_UMI_Matrix.tsv",sep='')
       randomization_fileName = paste("Randomizations_SN RNASeq UCSDWU_2020June22.txt",sep='')
       metadata_fileName = "Single_nucleus_RNASeq_UCSDWU_metadata.txt"
    }#End
    if (tis_center=="Premiere")
    {#Begin
       rawData_fileName = paste("Premiere_Raw.data_062220.txt",sep='')
       randomization_fileName = paste("Randomizations_SC RNASeq PREMIERE_2020June22.txt",sep='')
       metadata_fileName = "Single_cell_RNASeq_PREMIERE_metadata.txt"
    }#End
    complete_fileName_essential_genes = paste(base_directory,"/Experimental_data/Sample_metadata_additional_datasets/Essential_genes.txt",sep='')

    ###Generate results directory - BEGIN
    baseFileName = paste(tis_center,"_",set,sep='')
    results_directory = paste(base_directory,"Results/PostHocPower_seurat/",baseFileName,"/",sep='');
    dir.create(results_directory,recursive = TRUE)
    results_documentation_directory = paste(results_directory,"Documentations/",sep='')
    dir.create(results_documentation_directory)
    ###Generate results directory - END
}#End - Specify file names and generate directories

{#Begin - Write version documentations
   version_control_fileName = "AA_version_documentation.txt"
   r_session_fileName = "AA_R_session.txt"
   python_configuration_fileName = "AA_python_config.txt"
   complete_version_control_fileName = paste(results_documentation_directory,version_control_fileName,sep='')
   complete_r_session_fileName = paste(results_documentation_directory,r_session_fileName,sep='')
   #complete_python_configuration_fileName = paste(results_documentation_directory,python_configuration_fileName,sep='')

   sink(file=complete_r_session_fileName)
   r_sessionInfo
   sink()

   #sink(file=complete_python_configuration_fileName)
   #python_sessionInfo
   #sink()

   write.table(version_documentations,file=complete_version_control_fileName,quote=FALSE,row.names=FALSE,col.names=TRUE,sep='\t')
}#End - Write version documentations

{#Begin - Read data
   #Region Read randomziations
   complete_randomization_fileName = paste(metadata_directory, randomization_fileName,sep='')
   Randomizations = read.table(file=complete_randomization_fileName,header = TRUE,stringsAsFactors = FALSE,sep='\t')
   indexKeep = which(Randomizations$Considered_patients_count>1)
   Randomizations = Randomizations[indexKeep,]
   #endregion

   #Region Read rawdata
   complete_rawDataFileName = paste(directory,rawData_fileName,sep='')
   list.files(directory)
   Raw_data_matrix_input = read.table(file=complete_rawDataFileName, header=TRUE, check.names=FALSE)
   rownames(Raw_data_matrix_input) = gsub("_","-",rownames(Raw_data_matrix_input))
   #endregion

   #Region Read metadata
   complete_metadata_fileName = paste(metadata_directory,metadata_fileName,sep='')
   metadata = read.table(complete_metadata_fileName,header = TRUE,stringsAsFactors = FALSE,sep='\t')
   metadata$Library = gsub("-",".",metadata$Library)
   if (tis_center == "UCSD") { metadata$Library = paste(metadata$Library,"_",sep='') }
   #endregion

   #Region Read essential genes
   cellType_essential_genes = read.table(file=complete_fileName_essential_genes,header=TRUE,sep='\t',stringsAsFactors = FALSE)
   #endregion
}#End- Read data

{#Begin - Specify parameter
   if (tis_center=="Premiere")
   {#Begin
      minimum_feature_count_per_cell = 500;
   }#End
   if (tis_center=="UCSD")
   {#Begin
      minimum_feature_count_per_cell = 400;
   }#End

   maximum_feature_count_per_cell = 5000;
   minimum_cells = 3;
   max_percentage_mitochondrial_genes = 50;
   top_considered_features = 2000;
   max_pc_dimensions = c(30);
   cluster_selected_resolutions = (5:20)/10
   additional_vars_to_regress = c("");#c("percent.mt");#percent.mt","S.Score","G2M.Score");
   only_upregulated_genes = TRUE;
   mitochondrial_gene_label = "^MT-"#
   alpha_cluster_DEGs_adj_pvalues = 0.05;
   keep_top_ranked_marker_for_cellType_detection=300;
   minimum_distance_between_minLog10Pvalues_of_first_and_second_predicted_cellType = -log10(0.05);
   findMarkersAssays = c("SCT")#Currently, in this version, only the first assay will be used
   findMarkersSlots = c("data")#Currently, in this version, only the first slot will be used
   if (tis_center == "UCSD")
   {#Begin
      cellTypes_to_be_detected = unique(cellType_essential_genes$Cell_type)
      cellTypes_to_be_detected = cellTypes_to_be_detected[!cellTypes_to_be_detected %in% c("MAC","MON","Tcells","EC-AEA-DVR")]
   }#End
   if (tis_center == "Premiere")
   {#Begin
      cellTypes_to_be_detected = unique(cellType_essential_genes$Cell_type)
   }#End
}#End - Specify parameter

{#Begin - Document parameters
   Col_names = c("Parameter","Value")
   Col_length = length(Col_names);
   Row_names = 1;
   Row_length= length(Row_names)
   parameter_base_line = array(NA,c(Row_length,Col_length),dimnames=list(Row_names, Col_names))
   parameter_base_line = as.data.frame(parameter_base_line)
   parameters = c()

   parameter_line = parameter_base_line
   parameter_line$Parameter = "Base file name";
   parameter_line$Value = baseFileName;
   parameters = parameter_line

   parameter_line = parameter_base_line
   parameter_line$Parameter = "Dataset";
   parameter_line$Value = set;
   parameters = rbind(parameters,parameter_line)

   parameter_line = parameter_base_line
   parameter_line$Parameter = "randomization_fileName";
   parameter_line$Value = randomization_fileName;
   parameters = rbind(parameters,parameter_line)

   parameter_line = parameter_base_line
   parameter_line$Parameter = "rawData_fileName";
   parameter_line$Value = rawData_fileName;
   parameters = rbind(parameters,parameter_line)

   parameter_line = parameter_base_line
   parameter_line$Parameter = "Minimum feature count per cell";
   parameter_line$Value = minimum_feature_count_per_cell;
   parameters = rbind(parameters,parameter_line)

   parameter_line = parameter_base_line
   parameter_line$Parameter = "maximum_feature_count_per_cell";
   parameter_line$Value = maximum_feature_count_per_cell;
   parameters = rbind(parameters,parameter_line)

   parameter_line = parameter_base_line
   parameter_line$Parameter = "minimum_cells";
   parameter_line$Value = minimum_cells;
   parameters = rbind(parameters,parameter_line)

   parameter_line = parameter_base_line
   parameter_line$Parameter = "max_percentage_mitochondrial_genes";
   parameter_line$Value = max_percentage_mitochondrial_genes;
   parameters = rbind(parameters,parameter_line)

   parameter_line = parameter_base_line
   parameter_line$Parameter = "mitochondrial_gene_label ";
   parameter_line$Value = mitochondrial_gene_label;
   parameters = rbind(parameters,parameter_line)

   parameter_line = parameter_base_line
   parameter_line$Parameter = "top_considered_features";
   parameter_line$Value = top_considered_features;
   parameters = rbind(parameters,parameter_line)

   parameter_line = parameter_base_line
   parameter_line$Parameter = "additional_vars_to_regress";
   parameter_line$Value = paste(additional_vars_to_regress,sep='',collapse=", ");
   parameters = rbind(parameters,parameter_line)

   parameter_line = parameter_base_line
   parameter_line$Parameter = "cluster_selected_resolutions";
   parameter_line$Value = paste(cluster_selected_resolutions,sep='',collapse=", ");
   parameters = rbind(parameters,parameter_line)

   parameter_line = parameter_base_line
   parameter_line$Parameter = "Find markers assays";
   parameter_line$Value = paste(findMarkersAssays,sep='',collapse=", ");
   parameters = rbind(parameters,parameter_line)

   parameter_line = parameter_base_line
   parameter_line$Parameter = "Find markers slots";
   parameter_line$Value = paste(findMarkersSlots,sep='',collapse=", ");
   parameters = rbind(parameters,parameter_line)

   parameter_line = parameter_base_line
   parameter_line$Parameter = "Cell types to be detected";
   parameter_line$Value = paste(cellTypes_to_be_detected,collapse=';')
   parameters = rbind(parameters,parameter_line)

   parameter_line = parameter_base_line
   parameter_line$Parameter = "only upregulated genes";
   parameter_line$Value = only_upregulated_genes;
   parameters = rbind(parameters,parameter_line)

   parameter_line = parameter_base_line
   parameter_line$Parameter = "adj pvalue alpha for cluster DEGs";
   parameter_line$Value = alpha_cluster_DEGs_adj_pvalues;
   parameters = rbind(parameters,parameter_line)

   parameter_line = parameter_base_line
   parameter_line$Parameter = "keep_top_ranked_marker_for_cellType_detection";
   parameter_line$Value = keep_top_ranked_marker_for_cellType_detection;
   parameters = rbind(parameters,parameter_line)

   parameter_line = parameter_base_line
   parameter_line$Parameter = "Essential genes file name";
   parameter_line$Value = basename(complete_fileName_essential_genes)
   parameters = rbind(parameters,parameter_line)

   parameter_line = parameter_base_line
   parameter_line$Parameter = "minimum_distance_between_minLog10Pvalues_of_first_and_second_predicted_cellType"
   parameter_line$Value = minimum_distance_between_minLog10Pvalues_of_first_and_second_predicted_cellType;
   parameters = rbind(parameters,parameter_line)

   parameter_line = parameter_base_line
   parameter_line$Parameter = "parallel cores count";
   parameter_line$Value = cores_count;
   parameters = rbind(parameters,parameter_line)

   system_time = Sys.time()
   system_time = gsub(":","_",system_time)
   system_time = gsub("-","_",system_time)
   system_time = gsub(" ","_",system_time)

   png_completeFileName = paste(results_documentation_directory,"Parameter_table_",system_time,".png",sep='');
   png(png_completeFileName, width=1000, height=1000);
   grid.table(parameters, rows=NULL);#, show.rownames=FALSE)
   dev.off()

   complete_paramter_fileName = paste(results_documentation_directory,"Parameters_table_",system_time,".txt",sep='');
   write.table(parameters,complete_paramter_fileName,row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
}#End - Document parameters

#currentRando_nos = c(4441,4443,4445,3766,4447,3568,3608,3610,3691,3615,3658,3659,60,61,62,66,69,73,74,76,78,3683,3682)
#indexKeep = which(Randomizations$Randomization_no %in% c(2))
#Randomizations = Randomizations[indexKeep,]

{#Begin - Register and fill parrallel cores
  #registerDoParallel(cores=cores_count)
  length_randomizations = length(Randomizations[,1])
  randomizations_per_core = length_randomizations/cores_count

  #Randomizations = Randomizations[(length(Randomizations[,1])-3):length(Randomizations[,1]),];

  if (length_randomizations<cores_count) { cores_count = length_randomizations }

  parallel_clusters = makeCluster(cores_count)
  clusterEvalQ(parallel_clusters, {library(SingleCellExperiment)
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
               use_python("c://ProgramData//Anaconda3//")
               library(umap)
               library(doParallel)}
  );

  clusterExport(cl=parallel_clusters, varlist= list("set","minimum_feature_count_per_cell","maximum_feature_count_per_cell","minimum_cells","max_percentage_mitochondrial_genes","metadata",
                                           "top_considered_features","max_pc_dimensions","cluster_selected_resolutions","additional_vars_to_regress","only_upregulated_genes","mitochondrial_gene_label","cellType_essential_genes",
                                           "findMarkersAssays","findMarkersSlots","alpha_cluster_DEGs_adj_pvalues","keep_top_ranked_marker_for_cellType_detection","complete_fileName_essential_genes","cellTypes_to_be_detected",
                                           "minimum_distance_between_minLog10Pvalues_of_first_and_second_predicted_cellType",
                                           "Raw_data_matrix_input","baseFileName","directory","results_directory","tis_center","results_documentation_directory"),
                envir=.GlobalEnv);

  randomizations_list = vector(mode = "list", length = cores_count)
  indexCore = 0;
  Randomizations$Core = -1;

  for (indexRando in 1:length_randomizations)
  {#Begin
    indexCore = indexCore+1;
    if (indexCore>cores_count) {indexCore=1}
    randomizations_list[[indexCore]] = c(randomizations_list[[indexCore]],indexRando)
    Randomizations$Core = indexCore
  }#End

  all_randos = c()
  for (indexCore in 1:cores_count)
  {#Begin
    all_randos = c(all_randos,randomizations_list[[indexCore]])
  }#End

  all_randos = all_randos[order(all_randos)]
  length(all_randos)
  length(unique(all_randos))

  for (indexCore in 1:round(cores_count/2))
  {#Begin
    randomizations_list[[indexCore]] = rev(randomizations_list[[indexCore]])
  }#End

  for (indexCore in 1:cores_count)
  {#Begin
    current_indexes = randomizations_list[[indexCore]]
    #startIndex = min(floor((indexCore-1) * randomizations_per_core+1),length_randomizations);
    #endIndex = min(floor(indexCore * randomizations_per_core),length_randomizations)
    core_randomizations = Randomizations[current_indexes,]
    core_randomizations$Cluster_no = indexCore
    core_randomizations_length = length(core_randomizations[,1])
    clusterCall(parallel_clusters[indexCore], function(d) {assign('core_randomizations', d, pos=.GlobalEnv)}, core_randomizations)
    clusterCall(parallel_clusters[indexCore], function(d) {assign('indexCore', d, pos=.GlobalEnv)}, indexCore)
  }#End

  cluster_generation_correct = TRUE;

  combined_core_randomizations <- do.call('rbind', clusterEvalQ(parallel_clusters, core_randomizations))
  combined_core_randomizations = combined_core_randomizations[order(combined_core_randomizations$Randomization_no),]
  if (length(combined_core_randomizations[,1]) != length(Randomizations[,1])) { cluster_generation_correct = FALSE }
  if (cluster_generation_correct)
  {#Begin
    for (indexR in 1:length(combined_core_randomizations[,1]))
    {#Begin
      if (combined_core_randomizations$Randomization_no[indexR] != Randomizations$Randomization_no[indexR]) { cluster_generation_correct=FALSE;}
      if (combined_core_randomizations$Considered_patients_count[indexR] != Randomizations$Considered_patients_count[indexR]) { cluster_generation_correct=FALSE;}
      if (combined_core_randomizations$ReadWrite_considered_patients[indexR] != Randomizations$ReadWrite_considered_patients[indexR]) { cluster_generation_correct=FALSE;}
      #if (combined_core_randomizations$ReadWrite_considered_centers[indexR] != Randomizations$ReadWrite_considered_centers[indexR]) { cluster_generation_correct=FALSE;}
      if (combined_core_randomizations$Finished[indexR] != Randomizations$Finished[indexR]) { cluster_generation_correct=FALSE;}
    }#End
  }#End

  if (!cluster_generation_correct) { stopCluster(parallel_clusters) }

  complete_randomization_docu_fileName = paste(results_documentation_directory,"Randomizations_assigned_",system_time,".txt",sep='')
  write.table(Randomizations,file=complete_randomization_docu_fileName,quote=FALSE,row.names=FALSE,sep='\t')
}#End - Register and fill parallel cores

clusterEvalQ(parallel_clusters,
{#Begin - Parallel foreach

   Generate_and_process_seurat_object = function(indexCore, set, Raw_data_matrix, current_randomization_string, minimum_cells, results_directory, results_documentation_directory, tis_center, metadata, minimum_feature_count_per_cell, maximum_feature_count_per_cell, mitochondrial_gene_label, max_percentage_mitochondrial_genes, additional_vars_to_regress, top_considered_features, Generate_figures_for_visualization)
   {#Begin - Generate_and_process_seurat_object
      Continue_function = TRUE
      if (Continue_function)
      {#Begin - Create seurat object
         seuset = tryCatch( { CreateSeuratObject(Raw_data_matrix,project = current_randomization_string , min.cells=minimum_cells,min.features=0) },
                            error=function(cond) { return (NA) },
                            #warning=function(w) {},
                            finally = {});
         if ((!exists("seuset"))||(typeof(seuset)!="S4"))
         {#Begin
            complete_analysisFinished_fileName = paste(results_directory,"AA_Analysis_interrupted_CreateSeuratObject_set",set,"_core",indexCore,"_",current_randomization_string,".txt",sep='')
            write.table(paste(current_randomization_string," interrupted: CreateSeuratObject",sep=''),file=complete_analysisFinished_fileName,quote=FALSE,row.names=FALSE,sep='')
            Continue_function = FALSE;
         }#End
      }#End - Create seurat object

      if (Continue_function)
      {#Begin - Label patients, tissue collection, tissue type
         barcodes = colnames(seuset$RNA);
         patients = replicate(length(barcodes),"error");
         centers = replicate(length(barcodes),"error");
         tissue_collections = replicate(length(barcodes),"error")
         tissue_types = replicate(length(barcodes),"error")
         length_patients = length(patients)
         for (indexP in 1:length_patients)
         {#Begin
            barcode = barcodes[indexP]
            splitStrings = strsplit(barcode,"_")[[1]]
            if (tis_center=="UCSD") { patient = splitStrings[2] }
            if (tis_center=="Pr") { patient = splitStrings[1] }
            patients[indexP] = patient;
            centers[indexP] = tis_center;
         }#End

         all_indexBarcodes = c()
         for (indexMetadata in 1:length(metadata[,1]))
         {#Begin
            current_metadata_line = metadata[indexMetadata,]
            indexBarcodes = grep(current_metadata_line$Library,barcodes)
            all_indexBarcodes = c(all_indexBarcodes,indexBarcodes)
            patients[indexBarcodes] = current_metadata_line$Patient_id
            tissue_collections[indexBarcodes] = current_metadata_line$Tissue_collection
            tissue_types[indexBarcodes] = current_metadata_line$Tissue_type
         }#End

         length(unique(all_indexBarcodes)) == length(all_indexBarcodes)
         indexError = which(tissue_types=="error")
         if (length(indexError)>0) { rm(seuset) }

         seuset$Patient = patients
         seuset$Tissue_collection = tissue_collections
         seuset$Tissue_type = tissue_types
      }#End - Label patients, tissue collection, tissue type

      if (Continue_function)
      {#Begin - Quality control: Removal of cell with too few and too many nFeatures and too high mitochondrial feature counts

         seuset[["percent.mt"]] <- PercentageFeatureSet(object = seuset, pattern = mitochondrial_gene_label)

         if (Generate_figures_for_visualization)
         {#Begin - Plot QC violin plots
            # Visualize QC metrics as a violin plot
            group = "Patient"
            vlnplot_nFeature_RNA = VlnPlot(object = seuset, features = c("nFeature_RNA"), ncol = 1, group.by=group)
            vlnplot_nFeature_RNA = vlnplot_nFeature_RNA + geom_hline(yintercept=minimum_feature_count_per_cell,col="green",size=1)
            vlnplot_nFeature_RNA = vlnplot_nFeature_RNA + geom_hline(yintercept=maximum_feature_count_per_cell,col="green",size=1)
            vlnplot_nCount_RNA = VlnPlot(object = seuset, features = c("nCount_RNA"), ncol = 1, group.by=group)
            vlnplot_percent_mt = VlnPlot(object = seuset, features = c("percent.mt"), ncol = 1, group.by=group)
            vlnplot_percent_mt = vlnplot_percent_mt + geom_hline(yintercept=max_percentage_mitochondrial_genes,col="green",size=1)

            complete_png_fileName = paste(results_documentation_directory,current_randomization_string,"_data_quality_plots.png",sep='')
            png(complete_png_fileName,width=9000,height=9000,res=350);
            grid.arrange(vlnplot_nFeature_RNA,vlnplot_nCount_RNA,vlnplot_percent_mt,nrow=3,ncol=1)
            dev.off()
         }#End - Plot QC violin plots

         if (length(unique(seuset$percent.mt))==1)
         {#Begin - Remove cells with too few or many features
            seuset = tryCatch( { subset(x = seuset, subset = (nFeature_RNA <= maximum_feature_count_per_cell) & (nFeature_RNA >= minimum_feature_count_per_cell))},
                               error=function(cond) { return (NA) },
                               #warning=function(w) {},
                               finally = {});
         }#End - Remove cells with too few or many features

         if (length(unique(seuset$percent.mt))>1)
         {#Begin - Remove cells with too few, too many features and too much % mitochondrial RNA
            seuset = tryCatch( { subset(x = seuset, subset = (nFeature_RNA <= maximum_feature_count_per_cell) & (nFeature_RNA >= minimum_feature_count_per_cell) & (percent.mt<=max_percentage_mitochondrial_genes))},
                               error=function(cond) { return (NA) },
                               #warning=function(w) {},
                               finally = {});
         }#End - Remove cells with too few, too many features and too much % mitochondrial RNA

         if ((!exists("seuset"))||(typeof(seuset)!="S4"))
         {#Begin
            complete_analysisFinished_fileName = paste(results_directory,"AA_Analysis_interrupted_subset_set",set,"_core",indexCore,"_",current_randomization_string,".txt",sep='')
            write.table(paste(current_randomization_string," interrupted: subset",sep=''),file=complete_analysisFinished_fileName,quote=FALSE,row.names=FALSE,sep='')
            Continue_function = FALSE;
         }#End

      }#End - Quality control: Removal of cell with too few and too many nFeatures and too high mitochondrial feature counts

      if (Continue_function)
      {#Begin - Add cell cycle genes, if among additional_vars_to_regress

         if (("S.Score" %in% additional_vars_to_regress)||("G2M.Score" %in% additional_vars_to_regress))
         {#Begin
            s.genes <- cc.genes$s.genes;
            g2m.genes <- cc.genes$g2m.genes;
            seuset = tryCatch( { CellCycleScoring(seuset,s.features=s.genes, g2m.features=g2m.genes, set.ident = TRUE)},
                               error=function(cond) { return (NA) },
                               #warning=function(w) {},
                               finally = {});
            if ((!exists("seuset"))||(typeof(seuset)!="S4"))
            {#Begin
               complete_analysisFinished_fileName = paste(results_directory,"AA_Analysis_interrupted_cellcycle_set",set,"_core",indexCore,"_randoNo",current_randomization_no,".txt",sep='')
               write.table(paste("RandoNo",current_randomization_no," interrupted: cellcycle",sep=''),file=complete_analysisFinished_fileName,quote=FALSE,row.names=FALSE,sep='')
               Continue_function = FALSE;
            }#End
         }#End

      }#End - Add cell cycle genes, if among additional_vars_to_regress

      if (Continue_function)
      {#Begin - SC transform data
         seuset = tryCatch( { SCTransform(seuset, verbose=TRUE,variable.features.n = top_considered_features) },
                            error=function(cond) { return (NA) })
         #warning=function(w) {},
         #finally = {});
         if ((!exists("seuset"))||(typeof(seuset)!="S4"))
         {#Begin
            complete_analysisFinished_fileName = paste(results_directory,"AA_Analysis_interrupted_SCTransform_set",set,"_core",indexCore,"_",current_randomization_string,".txt",sep='')
            write.table(paste(current_randomization_string," interrupted: SCTransform",sep=''),file=complete_analysisFinished_fileName,quote=FALSE,row.names=FALSE,sep='')
            Continue_function = FALSE;
         }#End
      }#End - SC transform data

      if (Continue_function)
      {#Begin - Perform linear dimensional reduction
         seuset = tryCatch( { RunPCA(object = seuset, features = VariableFeatures(object = seuset)) },
                            error=function(cond) { return (NA) },
                            #warning=function(w) {},
                            finally = {});
         if ((!exists("seuset"))||(typeof(seuset)!="S4"))
         {#Begin
            complete_analysisFinished_fileName = paste(results_directory,"AA_Analysis_interrupted_RunPCA_set",set,"_core",indexCore,"_",current_randomization_string,".txt",sep='')
            write.table(paste(current_randomization_string," interrupted: RunPCA",sep=''),file=complete_analysisFinished_fileName,quote=FALSE,row.names=FALSE,sep='')
            Continue_function=FALSE;
         }#End

         if (Continue_function)
         {#Begin
            pca_genes = seuset[["pca"]]
            pca_genes = as.data.frame(pca_genes[,])

            pca_fileName = paste(results_directory,current_randomization_string,"_PCA_genes",".txt",sep='')
            pca_genes$Gene_symbol = rownames(pca_genes);
            write.table(pca_genes,file=pca_fileName,row.names=FALSE,quote=FALSE,sep='\t')
         }#End
      }#End - Perform linear dimensional reduction
      return (seuset)
   }#End - Generate_and_process_seurat_object

   Identify_neighbors_and_generate_umaps_and_tsne_if_generate_figures = function(indexCore, set, seuset, results_directory, current_randomization_string, pc_dimensions, Generate_figures_for_visualization)
   {#Begin - Identify_neighbors_and_generate_umaps_and_tsne_if_generate_figures
      Continue_function = TRUE
      seuset = tryCatch( { FindNeighbors(object = seuset, dims = pc_dimensions) },
                         error=function(cond) { return (NA) },
                         #warning=function(w) {},
                         finally = {});
      if ((!exists("seuset"))||(typeof(seuset)!="S4"))
      {#Begin
         complete_analysisFinished_fileName = paste(results_directory,"AA_Analysis_interrupted_FindNeighbors_set",set,"_core",indexCore,"_",current_randomization_string,".txt",sep='')
         write.table(paste(current_randomization_string," interrupted: FindNeighbors",sep=''),file=complete_analysisFinished_fileName,quote=FALSE,row.names=FALSE,sep='')
         Continue_function = FALSE
      }#End

      if ((Generate_figures_for_visualization)&(Continue_function))
      {#Begin - Do tsne and umap dimensionality reductions
         seuset <- RunUMAP(seuset, reduction="pca", dims = pc_dimensions, verbose = FALSE)
         seuset <- RunTSNE(seuset, npcs = 30, reduction="pca", verbose = FALSE)
      }#End - Do tsne and umap dimensionality reductions
      return (seuset)
   }#End - Identify_neighbors_and_generate_umaps_and_tsne_if_generate_figures

    core_randomizations_length = length(core_randomizations[,1])
    bg_genes_raw_data_matrix_input = unique(rownames(Raw_data_matrix_input))
    #complete_cluster_rawData_matrix_input_file = paste(results_directory,"Raw_data_input_core",indexCore,".txt");
    #write.table(Raw_data_matrix_input,file=complete_cluster_rawData_matrix_input_file,row.names = TRUE, col.names = TRUE,sep='\t',quote = FALSE)

    for (indexRando in 1:core_randomizations_length)
    {#Begin
       continue_randomizations_loop = TRUE
       #Region Read rawdata
       #Raw_data_matrix_input2 = read.table(file=complete_cluster_rawData_matrix_input_file, header=TRUE, check.names=FALSE)
       #endregion

       startTime_of_current_randomizationNo = Sys.time()

       Generate_figures_for_visualization = FALSE;

       current_randomization_no = core_randomizations$Randomization_no[indexRando];

       complete_analysisStarted_fileName = paste(results_directory,"AAA_Analysis_started_set",set,"_core",indexCore,"_randoNo",current_randomization_no,".txt",sep='')
       write.table(paste("RandoNo",current_randomization_no," finished",sep=''),file=complete_analysisStarted_fileName,quote=FALSE,row.names=FALSE,sep='')
       if (current_randomization_no==0) { Generate_figures_for_visualization = TRUE; }
       current_randomization_add_string = paste("_randoNo",current_randomization_no,sep='')
       current_randomization_string = paste("RandoNo",current_randomization_no,sep='')

       {#Begin - Generate Raw_data_matrix by removing libraries from Raw_data_matrix_input
            if (!exists("Raw_data_matrix_input"))
            {#Begin
               Raw_data_matrix_input = cbind(Raw_data_matrix,Raw_data_matrix_unused)
            }#End
            current_patients = strsplit(core_randomizations$ReadWrite_considered_patients[indexRando],';')[[1]]
            barcodes = colnames(Raw_data_matrix_input);
            patients_of_barcodes = barcodes;
            length_patients = length(patients_of_barcodes)
            for (indexP in 1:length_patients)
            {#Begin
               patient = patients_of_barcodes[indexP]
               splitStrings = strsplit(patient,"_")[[1]]
               if (tis_center=="UCSD") { patient = splitStrings[2] }
               if (tis_center=="Premiere") { patient = splitStrings[1] }
               patients_of_barcodes[indexP] = patient;
            }#End
            indexKeepBarcodes = which(patients_of_barcodes %in% current_patients)
            indexRemoveBarcodes = 1:length(patients_of_barcodes)
            indexRemoveBarcodes = indexRemoveBarcodes[!indexRemoveBarcodes %in% indexKeepBarcodes]
            Raw_data_matrix = Raw_data_matrix_input[,indexKeepBarcodes]
            Raw_data_matrix_unused = Raw_data_matrix_input[,indexRemoveBarcodes]
            rm(Raw_data_matrix_input)
       }#End - Generate Raw_data_matrix by removing libraries from Raw_data_matrix_input

       {#Begin - Create seurat object
          seuset = Generate_and_process_seurat_object(indexCore, set, Raw_data_matrix, current_randomization_string, minimum_cells, results_directory, results_documentation_directory, tis_center, metadata, minimum_feature_count_per_cell, maximum_feature_count_per_cell, mitochondrial_gene_label, max_percentage_mitochondrial_genes, additional_vars_to_regress, top_considered_features, Generate_figures_for_visualization)
          if ((!exists("seuset"))||(typeof(seuset)!="S4")) { continue_max_pc_dimensionts_loop = FALSE}
       }#End - Create seurat object

       if (continue_randomizations_loop)
       {#Begin - if (continue_randomizations_loop)
        for (indexMaxPcDimension in 1:length(max_pc_dimensions))
        {#Begin - indexMaxPcDimension
          continue_max_pc_dimensionts_loop = TRUE;
          max_pc_dimension = max_pc_dimensions[indexMaxPcDimension]
          max_pc_dimension_string = paste("maxPC",max_pc_dimension,sep='');
          pc_dimensions = 1:max_pc_dimension

          if ((!exists("seuset"))||(typeof(seuset)!="S4"))
          {#Begin - Rebuild Seuset, if missing
             seuset = Generate_and_process_seurat_object(indexCore, set, Raw_data_matrix, current_randomization_string, minimum_cells, results_directory, results_documentation_directory, tis_center, metadata, minimum_feature_count_per_cell, maximum_feature_count_per_cell, mitochondrial_gene_label, max_percentage_mitochondrial_genes, additional_vars_to_regress, top_considered_features, Generate_figures_for_visualization)
             if ((!exists("seuset"))||(typeof(seuset)!="S4")) { continue_max_pc_dimensionts_loop = FALSE}
          }#End - Rebuild Seuset, if missing

          if (continue_max_pc_dimensionts_loop)
          {#Begin - Identify neighbors
             seuset = Identify_neighbors_and_generate_umaps_and_tsne_if_generate_figures(indexCore, set, seuset, results_directory, current_randomization_string, pc_dimensions, Generate_figures_for_visualization)
          }#End - Identify neighbors

          if (continue_max_pc_dimensionts_loop)
          {#Begin - Reset cellType_was_detected_list
             cellType_was_detected_list = list()
             cellTypes_to_be_detected_length = length(cellTypes_to_be_detected)
             for (indexCellType in 1:cellTypes_to_be_detected_length)
             {#Begin
                cellType = cellTypes_to_be_detected[indexCellType]
                cellType_was_detected_list[[cellType]] = FALSE
             }#End
          }#End - Reset cellType_was_detected_list
          max_count_of_identified_cell_types_at_same_resolution = -1
          max_percent_of_cells_assigned_to_cell_types_count=-1

          if (continue_max_pc_dimensionts_loop)
          {#Begin - Define resolution per cell type lines
             Col_names = c("First_resolution_that_identified_cell_type","Previous_resolution","Cell_type_1st","Cell_type_2nd","Cell_type_1st_minLog10Pvalue","Cell_type_2nd_minLog10Pvalue","Diff_first_second_minLog10Pvalues","Minimum_diff_minLog10Pvalues","Cells_count","Cell_cluster_count","Cluster_no","Total_cluster_count")
             Col_length = length(Col_names)
             Row_names = 1
             Row_length = length(Row_names)
             resolution_documentation_base_line = array(NA,c(Row_length,Col_length),dimnames = list(Row_names,Col_names))
             resolution_documentation_base_line = as.data.frame(resolution_documentation_base_line)
             resolution_documentations = c();
          }#End - Define resolution per cell type lines

          if (continue_max_pc_dimensionts_loop)
          {#Begin - Define resolutions per 1st and 2nd cell type differences
             Col_names = c("Selected","Resolution","Cluster_count","Cell_types_count","Percent_of_cells_assigned_to_cell_types_count","Cells_assigned_to_cell_types_count","Cells_count","New_cell_types","Disappearing_cell_types","Quantile02_minusLog10Differences","Quantile05_minusLog10Differences","Quantile10_minusLog10Differences","Average_minusLog10Differences","Median_minusLog10Differences","SD_minusLog10Differences","Valid")

             Col_length = length(Col_names)
             Row_names = 1
             Row_length = length(Row_names)
             resolution_minusLog10Differences_documentation_base_line = array(NA,c(Row_length,Col_length),dimnames = list(Row_names,Col_names))
             resolution_minusLog10Differences_documentation_base_line = as.data.frame(resolution_minusLog10Differences_documentation_base_line)
             resolution_minusLog10Differences_documentations = c();
          }#End - Define resolutions per 1st and 2nd cell type differences

          if (continue_max_pc_dimensionts_loop)
          {#Begin - Identify resolution for real analysis
             indexResolution=1
             length_cluster_selected_resolutions = length(cluster_selected_resolutions)
             previous_cluster_selected_resolution = -1;
             previous_cluster_count = -1;
             all_cell_types_detected = FALSE
             indexFindMarkersAssay = 1
             for (indexResolution in 1:length_cluster_selected_resolutions)
             {#Begin - indexResolution
                if (!all_cell_types_detected)
                {#Begin - if (!all_cell_types_detected)

                   continue_resolution_loop = TRUE
                   if ((!exists("seuset"))||(typeof(seuset)!="S4"))
                   {#Begin - Rebuild Seuset and identify neighbors
                      seuset = Generate_and_process_seurat_object(indexCore, set, Raw_data_matrix, current_randomization_string, minimum_cells, results_directory, results_documentation_directory, tis_center, metadata, minimum_feature_count_per_cell, maximum_feature_count_per_cell, mitochondrial_gene_label, max_percentage_mitochondrial_genes, additional_vars_to_regress, top_considered_features, Generate_figures_for_visualization)
                      if ((!exists("seuset"))||(typeof(seuset)!="S4")) { continue_resolution_loop = FALSE}
                      if (continue_resolution_loop)
                      {#Begin
                         seuset = Identify_neighbors_and_generate_umaps_and_tsne_if_generate_figures(indexCore, set, seuset, results_directory, current_randomization_string, pc_dimensions, Generate_figures_for_visualization)
                         if ((!exists("seuset"))||(typeof(seuset)!="S4")) { continue_resolution_loop = FALSE }
                      }#End
                   }#End - Rebuild Seuset and identify neighbors

                   if (continue_resolution_loop)
                   {#Begin - Define cluster selected resolution and strings for file names
                      cluster_selected_resolution = cluster_selected_resolutions[indexResolution];
                      cluster_selected_resolution_string = paste("res",cluster_selected_resolution,sep='')
                      current_randomization_plus_maxPcDimension = paste(current_randomization_string,"_",max_pc_dimension_string,sep='')
                   }#End - Define cluster selected resolution and strings for file names

                   if (continue_resolution_loop)
                   {#Begin - Identify clusters and set cluster_count
                      seuset = tryCatch( { FindClusters(object = seuset, resolution = cluster_selected_resolution) },
                                           error=function(cond) { return (NA) },
                                           finally = {});

                      if ((!exists("seuset"))||(typeof(seuset)!="S4"))
                      { continue_resolution_loop = FALSE; }

                      if (continue_resolution_loop)
                      {#Begin
                         seuset$Cluster <- as.numeric(as.character(Idents(seuset)));
                         cluster_count = length(table(seuset$Cluster))
                      }#End
                   }#End - Identify clusters and set cluster_count

                   if ((continue_resolution_loop)&(cluster_count>previous_cluster_count))
                   {#Begin - if (cluster_count>previous_cluster_count)
                       previous_cluster_count = cluster_count

                       if (continue_resolution_loop)
                       {#Begin - Identify markers of each cluster versus all other clusters
                          Idents(seuset) = seuset$Cluster

                          findMarkersAssay = findMarkersAssays[indexFindMarkersAssay]
                          findMarkersSlot = findMarkersSlots[indexFindMarkersAssay]

                          if (continue_resolution_loop)
                          {#Begin - Find markers using FindAllMarkers
                              clusterNOs = unique(seuset$Cluster)
                              combined_clusterDEGs <- tryCatch( { FindAllMarkers(object = seuset, min.pct = 0.1, only.pos = only_upregulated_genes, assay=findMarkersAssay, slot=findMarkersSlot) },
                                                                  error=function(cond) { return (NA) },
                                                                  finally = {});
                              if ((!exists("combined_clusterDEGs"))||(is.null(combined_clusterDEGs)))
                              { continue_resolution_loop = FALSE }

                              if (continue_resolution_loop)
                              {#Begin - Finalize combined_clusterDEGs
                                 combined_clusterDEGs$cluster = as.numeric(as.character(combined_clusterDEGs$cluster))
                                 missing_clusterNOs = which(!clusterNOs %in% combined_clusterDEGs$cluster)

                                 if (length(missing_clusterNOs)>0)
                                 {#Begin - Add noDEGs lines for missing clusters
                                    Col_names = c("p_val","avg_log2FC","pct.1","pct.2","p_val_adj","gene","cluster")
                                    Col_length = length(Col_names)
                                    Row_names = missing_clusterNOs
                                    Row_length = length(Row_names)
                                    add_cluster.markers = as.data.frame(array(NA,c(Row_length,Col_length),dimnames=list(Row_names,Col_names)))
                                    add_cluster.markers$cluster = missing_clusterNOs
                                    add_cluster.markers$p_val = 1;
                                    add_cluster.markers$avg_log2FC = 0;
                                    add_cluster.markers$pct.1 = 0;
                                    add_cluster.markers$pct.2 = 0;
                                    add_cluster.markers$p_val_adj = 1;
                                    add_cluster.markers$gene = "noDEGs"
                                    combined_clusterDEGs = rbind(combined_clusterDEGs,add_cluster.markers)
                                 }#End - Add noDEGs lines for missing clusters

                                 combined_clusterDEGs$Resolution = cluster_selected_resolution;
                                 combined_clusterDEGs$Assay = findMarkersAssay;
                                 combined_clusterDEGs$Slot = findMarkersSlot;
                                 combined_clusterDEGs$Cell_type_1st = ""
                                 combined_clusterDEGs$Cell_type_1st_pvalue = -1;
                                 combined_clusterDEGs$Cell_type_2nd = ""
                                 combined_clusterDEGs$Cell_type_2nd_pvalue = -1;
                                 combined_clusterDEGs = as.data.frame(combined_clusterDEGs)
                                 indexNa = which(is.na(combined_clusterDEGs$p_val))
                                 if (length(indexNa)>0) { combined_clusterDEGs$p_val[indexNa] = 1; }
                                 indexNa = which(is.na(combined_clusterDEGs$p_val_adj))
                                 if (length(indexNa)>0) { combined_clusterDEGs$p_val_adj[indexNa] = 1; }
                                 indexNa = which(is.na(combined_clusterDEGs$avg_logFC))
                                 if (length(indexNa)>0) { combined_clusterDEGs$avg_logFC[indexNa] = 0; }
                                 combined_clusterDEGs = as.data.frame(combined_clusterDEGs)
                              }#End - Finalize combined_clusterDEGs
                          }#End - Find markers using FindAllMarkers
                       }#End - Identify markers of each cluster versus all other clusters

                       if (continue_resolution_loop)
                       {#Begin - Assign cell types to clusters and update cellType_was_detected_list
                          cellTypes_to_be_detected_length = length(cellTypes_to_be_detected)
                          cellTypes = names(cellType_was_detected_list)
                          bg_genes = bg_genes_raw_data_matrix_input
                          length_bg_genes = length(bg_genes)
                          indexEgInBg = which(cellType_essential_genes$Essential_gene %in% bg_genes)
                          bg_essential_genes = cellType_essential_genes[indexEgInBg,]
                          eg_cellTypes = unique(bg_essential_genes$Cell_type)
                          length_eg_cellTypes = length(eg_cellTypes)
                          cellType = names(cellType_was_detected_list)[indexCellType]
                          indexCurrentClusterNo = 1
                          current_cluster_nos = unique(combined_clusterDEGs$cluster)
                          length_cluster = length(current_cluster_nos)
                          seuset$Cell_type_1st = "Not assigned"; #String is used below, if change here, change below as well
                          seuset$Cell_type_1st_minLog10Pvalue = -1;
                          seuset$Cell_type_2nd = "Not assigned";
                          seuset$Cell_type_2nd_minLog10Pvalue = -1;
                          detectedCellTypes_clusterCount_list = list();
                          detectedCellTypes_minLog10Pvalue_list = list();
                          detectedCellTypes_secondCellType_list = list();
                          detectedCellTypes_secondMinLog10Pvalue_list = list();
                          detectedCellTypes_cellsCount_list = list()
                          detectedCellTypes_clusterNo_list = list()
                          all_minusLog10Pvalue_distances = c()
                          for (indexCurrentClusterNo in 1:length_cluster)
                          {#Begin - Identify cell types for each cluster and change cell type detected and set lastIndexResolution_of_new_cellType_detected
                              current_cluster_no = current_cluster_nos[indexCurrentClusterNo]
                              indexCurrentDegs = which(combined_clusterDEGs$cluster == current_cluster_no)
                              current_cluster_degs = combined_clusterDEGs[indexCurrentDegs,]
                              current_cluster_degs = current_cluster_degs[order(current_cluster_degs$p_val_adj),]
                              current_cluster_degs = current_cluster_degs[1:keep_top_ranked_marker_for_cellType_detection,]
                              indexSig = which(current_cluster_degs$p_val_adj <= alpha_cluster_DEGs_adj_pvalues)
                              current_cluster_symbols = unique(current_cluster_degs$gene[indexSig])
                              current_cluster_symbols_count = length(current_cluster_symbols)
                              lowest_pvalue = 99;
                              second_lowest_pvalue = 99;
                              lowest_pvalue_cell_type = "";
                              second_lowest_pvalue_cell_type = "";
                              indexEG_cellTypes=1
                              for (indexEG_cellTypes in 1:length_eg_cellTypes)
                              {#Begin - Identify most and 2nd most significant cell Type using Fisher's Exact Test
                                  eg_cellType = eg_cellTypes[indexEG_cellTypes]
                                  indexCurrentEssential = which(bg_essential_genes$Cell_type==eg_cellType)
                                  current_essential_gene_symbols = unique(bg_essential_genes$Essential_gene[indexCurrentEssential])
                                  overlap_genes = current_cluster_symbols[current_cluster_symbols %in% current_essential_gene_symbols]
                                  a = length(overlap_genes)
                                  b = current_cluster_symbols_count - a
                                  c = length(current_essential_gene_symbols) - a
                                  d = length_bg_genes - a - b - c

                                  confidence_table = array(NA,c(2,2))
                                  confidence_table[1,1] = a
                                  confidence_table[2,1] = b
                                  confidence_table[1,2] = c
                                  confidence_table[2,2] = d

                                  fisher_test = fisher.test(confidence_table,alternative="greater")
                                  current_pvalue = fisher_test$p.value
                                  if (current_pvalue <= lowest_pvalue)
                                  {#Begin
                                      if (lowest_pvalue <= second_lowest_pvalue)
                                      {#Begin
                                          second_lowest_pvalue = lowest_pvalue;
                                          second_lowest_pvalue_cell_type = lowest_pvalue_cell_type;
                                      }#End
                                      lowest_pvalue = current_pvalue;
                                      lowest_pvalue_cell_type = eg_cellType
                                  }#End
                                  if ((current_pvalue > lowest_pvalue) & (current_pvalue <= second_lowest_pvalue))
                                  {#Begin
                                      second_lowest_pvalue = current_pvalue;
                                      second_lowest_pvalue_cell_type = eg_cellType
                                  }#End
                              }#End - Identify most and 2nd most significant cell Type using Fisher's Exact Test
                              minusLog10_first = 0
                              if (lowest_pvalue <= 1)
                              { minusLog10_first = -log10(lowest_pvalue) }
                              minusLog10_second = 0
                              if (second_lowest_pvalue <= 1)
                              { minusLog10_second = -log10(second_lowest_pvalue) }
                              currentCluster_cell_counts = length(which(seuset$Cluster == current_cluster_no))
                              if ((minusLog10_first - minusLog10_second) >= minimum_distance_between_minLog10Pvalues_of_first_and_second_predicted_cellType)
                              {#Begin
                                  indexSeuset_currentClusterNo = which(seuset$Cluster == current_cluster_no)
                                  all_minusLog10Pvalue_distances = c(all_minusLog10Pvalue_distances,(minusLog10_first - minusLog10_second))
                                  seuset$Cell_type_1st[indexSeuset_currentClusterNo] = lowest_pvalue_cell_type;
                                  seuset$Cell_type_1st_minLog10Pvalue[indexSeuset_currentClusterNo] = minusLog10_first;
                                  seuset$Cell_type_2nd[indexSeuset_currentClusterNo] = second_lowest_pvalue_cell_type;
                                  seuset$Cell_type_2nd_minLog10Pvalue[indexSeuset_currentClusterNo] = minusLog10_second;
                                  if (lowest_pvalue_cell_type %in% names(detectedCellTypes_clusterCount_list))
                                  {#Begin
                                      detectedCellTypes_clusterCount_list[[lowest_pvalue_cell_type]] = detectedCellTypes_clusterCount_list[[lowest_pvalue_cell_type]] + 1
                                      detectedCellTypes_minLog10Pvalue_list[[lowest_pvalue_cell_type]] = c(detectedCellTypes_minLog10Pvalue_list[[lowest_pvalue_cell_type]],minusLog10_first);
                                      detectedCellTypes_secondCellType_list[[lowest_pvalue_cell_type]] = c(detectedCellTypes_secondCellType_list[[lowest_pvalue_cell_type]],second_lowest_pvalue_cell_type);
                                      detectedCellTypes_secondMinLog10Pvalue_list[[lowest_pvalue_cell_type]] = c(detectedCellTypes_secondMinLog10Pvalue_list[[lowest_pvalue_cell_type]],minusLog10_second);
                                      detectedCellTypes_cellsCount_list[[lowest_pvalue_cell_type]] = c(detectedCellTypes_cellsCount_list[[lowest_pvalue_cell_type]],currentCluster_cell_counts)
                                      detectedCellTypes_clusterNo_list[[lowest_pvalue_cell_type]] = c(detectedCellTypes_clusterNo_list[[lowest_pvalue_cell_type]],current_cluster_no)
                                  }#End
                                  else
                                  {#Begin
                                      detectedCellTypes_clusterCount_list[[lowest_pvalue_cell_type]] = 1;
                                      detectedCellTypes_minLog10Pvalue_list[[lowest_pvalue_cell_type]] = minusLog10_first;
                                      detectedCellTypes_secondCellType_list[[lowest_pvalue_cell_type]] = second_lowest_pvalue_cell_type;
                                      detectedCellTypes_secondMinLog10Pvalue_list[[lowest_pvalue_cell_type]] = minusLog10_second;
                                      detectedCellTypes_cellsCount_list[[lowest_pvalue_cell_type]] = currentCluster_cell_counts
                                      detectedCellTypes_clusterNo_list[[lowest_pvalue_cell_type]] = current_cluster_no
                                  }#End
                              }#End
                          }#End - Identify cell types for each cluster and change cell type detected and set lastIndexResolution_of_new_cellType_detected

                          allCelltypes = unique(seuset$Cell_type_1st)
                          allCelltypes = allCelltypes[!allCelltypes %in% "Not assigned"]
                          identified_cellTypes_count_at_current_resolution = length(allCelltypes)

                          indexDetectedInPreviousResolution = which(cellType_was_detected_list==TRUE)
                          cellTypesDetectedInPreviousResolution = names(cellType_was_detected_list)[indexDetectedInPreviousResolution]
                          indexDisappearingCellTypes = which(!cellTypesDetectedInPreviousResolution %in% allCelltypes)
                          disappearingCellTypes = cellTypesDetectedInPreviousResolution[indexDisappearingCellTypes]
                          indexNewCellTypes = which(!allCelltypes %in% cellTypesDetectedInPreviousResolution)
                          new_cellTypes = allCelltypes[indexNewCellTypes]

                          {#Begin - Document overall assignment success
                              resolution_minusLog10Differences_documentation_line = resolution_minusLog10Differences_documentation_base_line
                              resolution_minusLog10Differences_documentation_line$Resolution = cluster_selected_resolution
                              resolution_minusLog10Differences_documentation_line$Average_minusLog10Differences = mean(all_minusLog10Pvalue_distances)
                              resolution_minusLog10Differences_documentation_line$Median_minusLog10Differences = median(all_minusLog10Pvalue_distances)
                              resolution_minusLog10Differences_documentation_line$SD_minusLog10Differences = sd(all_minusLog10Pvalue_distances)
                              resolution_minusLog10Differences_documentation_line$Quantile10_minusLog10Differences = quantile(all_minusLog10Pvalue_distances,probs=0.1)
                              resolution_minusLog10Differences_documentation_line$Quantile05_minusLog10Differences = quantile(all_minusLog10Pvalue_distances,probs=0.05)
                              resolution_minusLog10Differences_documentation_line$Quantile02_minusLog10Differences = quantile(all_minusLog10Pvalue_distances,probs=0.02)
                              resolution_minusLog10Differences_documentation_line$Cluster_count = length_cluster
                              resolution_minusLog10Differences_documentation_line$Cell_types_count = identified_cellTypes_count_at_current_resolution
                              if (length(new_cellTypes)==0)
                              { resolution_minusLog10Differences_documentation_line$New_cell_types = "" }
                              if (length(new_cellTypes)==1)
                              { resolution_minusLog10Differences_documentation_line$New_cell_types = new_cellTypes }
                              if (length(new_cellTypes)>1)
                              { resolution_minusLog10Differences_documentation_line$New_cell_types = paste(new_cellTypes,collapse=';') }
                              if (length(disappearingCellTypes)==0)
                              { resolution_minusLog10Differences_documentation_line$Disappearing_cell_types = "" }
                              if (length(disappearingCellTypes)==1)
                              { resolution_minusLog10Differences_documentation_line$Disappearing_cell_types = disappearingCellTypes; }
                              if (length(disappearingCellTypes)>1)
                              { resolution_minusLog10Differences_documentation_line$Disappearing_cell_types = paste(disappearingCellTypes,collapse=';'); }
                              if (length(disappearingCellTypes)==0)
                              { resolution_minusLog10Differences_documentation_line$Valid = TRUE }
                              if (length(disappearingCellTypes)>0)
                              { resolution_minusLog10Differences_documentation_line$Valid = FALSE }
                              indexCellsNotAssigned = which(seuset$Cell_type_1st=="Not assigned")
                              resolution_minusLog10Differences_documentation_line$Cells_count = length(seuset$Cell_type_1st)
                              resolution_minusLog10Differences_documentation_line$Cells_assigned_to_cell_types_count = resolution_minusLog10Differences_documentation_line$Cells_count - length(indexCellsNotAssigned)
                              resolution_minusLog10Differences_documentation_line$Percent_of_cells_assigned_to_cell_types_count = 100 * resolution_minusLog10Differences_documentation_line$Cells_assigned_to_cell_types_count / resolution_minusLog10Differences_documentation_line$Cells_count
                              resolution_minusLog10Differences_documentation_line$Selected = FALSE
                              if (length(resolution_minusLog10Differences_documentations)==0) { resolution_minusLog10Differences_documentations = resolution_minusLog10Differences_documentation_line }
                              else { resolution_minusLog10Differences_documentations = rbind(resolution_minusLog10Differences_documentations,resolution_minusLog10Differences_documentation_line) }
                          }#End - Document overall assignment success

                          percent_of_cells_assigned_to_cell_types_count = resolution_minusLog10Differences_documentation_line$Percent_of_cells_assigned_to_cell_types_count

                          if (  ( (identified_cellTypes_count_at_current_resolution > max_count_of_identified_cell_types_at_same_resolution)
                                 |((identified_cellTypes_count_at_current_resolution == max_count_of_identified_cell_types_at_same_resolution)&(percent_of_cells_assigned_to_cell_types_count > max_percent_of_cells_assigned_to_cell_types_count)))
                                 &(length(disappearingCellTypes)==0))
                          {#Begin - Document detection of new cell types, update max count of identified cell types and write cluster DEGs
                              resolution_minusLog10Differences_documentations$Selected = FALSE
                              indexCurrentResolution = which(resolution_minusLog10Differences_documentations$Resolution==cluster_selected_resolution)
                              resolution_minusLog10Differences_documentations$Selected[indexCurrentResolution] = TRUE
                              if (length(new_cellTypes)>0)
                              {#Begin - Document resolution of new cell types and add new cell type to cellType detected list
                                  for (indexNew in 1:length(new_cellTypes))
                                  {#Begin - Document resolution of new cell types
                                      new_cell_type = new_cellTypes[indexNew]
                                      if (new_cell_type %in% names(cellType_was_detected_list))
                                      {#Begin
                                         if (cellType_was_detected_list[[new_cell_type]]==FALSE)
                                         {#Begin
                                            this_cellType_clusters_count = detectedCellTypes_clusterCount_list[[new_cell_type]]
                                            for (indexThisCellTypeCluster in 1:this_cellType_clusters_count)
                                            {#Begin
                                               resolution_documentation_line = resolution_documentation_base_line
                                               resolution_documentation_line$Cell_type_1st = new_cell_type
                                               resolution_documentation_line$Cell_type_1st_minLog10Pvalue = detectedCellTypes_minLog10Pvalue_list[[new_cell_type]][indexThisCellTypeCluster]
                                               resolution_documentation_line$Cell_type_2nd = detectedCellTypes_secondCellType_list[[new_cell_type]][indexThisCellTypeCluster]
                                               resolution_documentation_line$Cell_type_2nd_minLog10Pvalue = detectedCellTypes_secondMinLog10Pvalue_list[[new_cell_type]][indexThisCellTypeCluster]
                                               resolution_documentation_line$Diff_first_second_minLog10Pvalues = resolution_documentation_line$Cell_type_1st_minLog10Pvalue - resolution_documentation_line$Cell_type_2nd_minLog10Pvalue
                                               resolution_documentation_line$Minimum_diff_minLog10Pvalues = minimum_distance_between_minLog10Pvalues_of_first_and_second_predicted_cellType
                                               resolution_documentation_line$First_resolution_that_identified_cell_type = cluster_selected_resolution
                                               resolution_documentation_line$Previous_resolution = previous_cluster_selected_resolution
                                               resolution_documentation_line$Cell_cluster_count = detectedCellTypes_clusterCount_list[[new_cell_type]];
                                               resolution_documentation_line$Cells_count = detectedCellTypes_cellsCount_list[[new_cell_type]][indexThisCellTypeCluster];
                                               resolution_documentation_line$Cluster_no = detectedCellTypes_clusterNo_list[[new_cell_type]][indexThisCellTypeCluster];
                                               resolution_documentation_line$Total_cluster_count = length_cluster
                                               if (length(resolution_documentations)==0) { resolution_documentations = resolution_documentation_line; }
                                               else { resolution_documentations = rbind(resolution_documentations,resolution_documentation_line) }
                                            }#End
                                         }#End
                                      }#End
                                  }#End - Document resolution of new cell types
                              }#End - Document resolution of new cell types and add new cell type to cellType detected list
                              if (length(new_cellTypes)>0)
                              {#Begin - if (length(new_cellTypes)>0) - Add new cell type to cellType_was_detected_list
                                 for (indexNew in 1:length(new_cellTypes))
                                 {#Begin - Add new cell type to cellType_was_detected_list, if all previous cell types were re-identified
                                     new_cell_type = new_cellTypes[indexNew]
                                     cellType_was_detected_list[[new_cell_type]] = TRUE;
                                 }#End - Add new cell type to cellType_was_detected_list, if all previous cell types were re-identified
                              }#End - if (length(new_cellTypes)>0) - Add new cell type to cellType_was_detected_list
                              {#Begin - Identify and override cell type barcode annotations
                                  col_names = c("Barcode","Cluster","Cell_type","Cell_type_1st_minLog10Pvalue","Cell_type_2nd","Cell_type_2nd_minLog10Pvalue","Patient","Resolution")
                                  col_length = length(col_names);
                                  row_names = 1:length(seuset$Cell_type_1st)
                                  row_length = length(row_names)
                                  barcode_cellType_associations = as.data.frame(array(NA,c(row_length,col_length),dimnames=list(row_names,col_names)))
                                  barcode_cellType_associations$Barcode = names(seuset$Cell_type_1st)
                                  barcode_cellType_associations$Cluster = seuset$Cluster
                                  barcode_cellType_associations$Patient = seuset$Patient;
                                  barcode_cellType_associations$Cell_type = seuset$Cell_type_1st
                                  barcode_cellType_associations$Cell_type_1st_minLog10Pvalue = seuset$Cell_type_1st_minLog10Pvalue
                                  barcode_cellType_associations$Cell_type_2nd = seuset$Cell_type_2nd
                                  barcode_cellType_associations$Cell_type_2nd_minLog10Pvalue = seuset$Cell_type_2nd_minLog10Pvalue
                                  barcode_cellType_associations$Resolution = cluster_selected_resolution;
                              }#End - Identify and override cell type barcode annotations
                              {#Begin - Write cluster DEGs per resolution
                                  clusterDEGs_fileName = paste("ClusterDEGs_",current_randomization_plus_maxPcDimension,"_",findMarkersAssay,"_",findMarkersSlot,"_",cluster_selected_resolution_string,".txt",sep='');
                                  complete_clusterDEGs_fileName = paste(results_directory,clusterDEGs_fileName,sep='')
                                  write.table(combined_clusterDEGs,file=complete_clusterDEGs_fileName,quote=FALSE,row.names=FALSE,sep='\t')
                              }#End - Write cluster DEGs per resolution
                              max_count_of_identified_cell_types_at_same_resolution = identified_cellTypes_count_at_current_resolution
                              max_percent_of_cells_assigned_to_cell_types_count = percent_of_cells_assigned_to_cell_types_count
                          }#End - Document detection of new cell types, update max count of identified cell types and write cluster DEGs
                       }#End - Assign cell types to clusters and update cellType_was_detected_list

                       if (continue_resolution_loop)
                       {#Begin - Analyse if all cell types detected
                           length_cellTypes = length(names(cellType_was_detected_list))
                           all_cell_types_detected = !FALSE %in% (cellType_was_detected_list);
                       }#End - Analyse if all cell types detected

                       seuset$Cluster = NULL
                       seuset$Cell_type_1st = NULL
                       seuset$Cell_type_2nd = NULL
                       seuset$Cell_type_1st_minLog10Pvalue = NULL
                       seuset$Cell_type_2nd_minLog10Pvalue = NULL
                   }#End - if (cluster_count>previous_cluster_count)
             }#End - if (!all_cell_types_detected)
             }#End - indexResolution
          }#End - Identify resolution for real analysis

          if (continue_max_pc_dimensionts_loop)
          {#Begin - Write resolution documentations
             resolutions_fileName = paste("CellType_detection_resolutions_",current_randomization_plus_maxPcDimension,"_",findMarkersAssay,"_",findMarkersSlot,".txt",sep='');
             complete_resolutions_fileName = paste(results_directory,resolutions_fileName,sep='')
             write.table(resolution_documentations,file=complete_resolutions_fileName,quote=FALSE,row.names=FALSE,sep='\t')

             minLog10Pvalues_fileName = paste("CellType_detection_minLog10PDiffs_resolutions_",current_randomization_plus_maxPcDimension,"_",findMarkersAssay,"_",findMarkersSlot,".txt",sep='');
             complete_minLog10Pvalues_fileName = paste(results_directory,minLog10Pvalues_fileName,sep='')
             write.table(resolution_minusLog10Differences_documentations,file=complete_minLog10Pvalues_fileName,quote=FALSE,row.names=FALSE,sep='\t')
          }#End - Write resolution documentations

          if (continue_max_pc_dimensionts_loop)
          {#Begin - Perform real analysis at selected resolution
              continue_real_analysis = TRUE

              indexRealAnalysis = which(resolution_minusLog10Differences_documentations$Selected==TRUE)
              if (length(indexRealAnalysis)!=1) { rm(Raw_data_matrix); rm(seuset);}
              cluster_selected_resolution = resolution_minusLog10Differences_documentations$Resolution[indexRealAnalysis]

              {#Begin - Write barcode_cellType_associations that were identified above
                 complete_clusterBarcodes_fileName = paste(results_directory,"Barcode_cellType_associations_",current_randomization_plus_maxPcDimension,".txt",sep='')
                 write.table(barcode_cellType_associations,file=complete_clusterBarcodes_fileName,quote=FALSE,row.names=FALSE,sep='\t')
              }#End - Write barcode_cellType_associations that were identified above

              if ((!exists("seuset"))||(typeof(seuset)!="S4"))
              {#Begin - Rebuild Seuset and identify neighbors
                 seuset = Generate_and_process_seurat_object(indexCore, set, Raw_data_matrix, current_randomization_string, minimum_cells, results_directory, results_documentation_directory, tis_center, metadata, minimum_feature_count_per_cell, maximum_feature_count_per_cell, mitochondrial_gene_label, max_percentage_mitochondrial_genes, additional_vars_to_regress, top_considered_features, Generate_figures_for_visualization)
                 if ((!exists("seuset"))||(typeof(seuset)!="S4")) { continue_real_analysis = FALSE}
                 if (continue_real_analysis)
                 {#Begin
                    seuset = Identify_neighbors_and_generate_umaps_and_tsne_if_generate_figures(indexCore, set, seuset, results_directory, current_randomization_string, pc_dimensions, Generate_figures_for_visualization)
                    if ((!exists("seuset"))||(typeof(seuset)!="S4")) { continue_real_analysis = FALSE }
                 }#End
              }#End - Rebuild Seuset and identify neighbors

              if (continue_real_analysis)
              {#Begin - Reassign cell types as saved in barcode cellType associations
                  seuset$Cluster = barcode_cellType_associations$Cluster
                  seuset$Cell_type_1st = barcode_cellType_associations$Cell_type
                  seuset$Cell_type_2nd = barcode_cellType_associations$Cell_type_2nd
                  seuset$Cell_type_1st_minLog10Pvalue = barcode_cellType_associations$Cell_type_1st_minLog10Pvalue
                  seuset$Cell_type_2nd_minLog10Pvalue = barcode_cellType_associations$Cell_type_2nd_minLog10Pvalue
                  Idents(seuset) = seuset$Cell_type_1st
              }#End - Reassign cell types as saved in barcode cellType associations

              if ((continue_real_analysis)&(Generate_figures_for_visualization))
              {#Begin - Plot cell neighborhoods - quality control
                   pt_size=3;
                   pt_size=0.01;

                   feature_plot_mt_tsn = FeaturePlot(seuset, features=c("percent.mt"),reduction="tsne") + ggtitle(paste("Percent mt TSNE - Res: ",cluster_selected_resolution,sep=''))
                   feature_plot_mt_pca = FeaturePlot(seuset, features=c("percent.mt"),reduction="pca") + ggtitle(paste("Percent mt PCA - Res: ",cluster_selected_resolution,sep=''))
                   feature_plot_mt_umap = FeaturePlot(seuset, features=c("percent.mt"),reduction="umap") + ggtitle(paste("Percent mt UMAP - Res: ",cluster_selected_resolution,sep=''))
                   feature_plot_nC_tsn = FeaturePlot(seuset, features=c("nCount_RNA"),reduction="tsne") + ggtitle(paste("nCount RNA TSNE - Res: ",cluster_selected_resolution,sep=''))
                   feature_plot_nC_pca = FeaturePlot(seuset, features=c("nCount_RNA"),reduction="pca") + ggtitle(paste("nCount RNA PCA - Res: ",cluster_selected_resolution,sep=''))
                   feature_plot_nC_umap = FeaturePlot(seuset, features=c("nCount_RNA"),reduction="umap") + ggtitle(paste("nCount RNA UMAP - Res: ",cluster_selected_resolution,sep=''))
                   feature_plot_nF_tsn = FeaturePlot(seuset, features=c("nFeature_RNA"),reduction="tsne") + ggtitle(paste("nFeature RNA TSNE - Res: ",cluster_selected_resolution,sep=''))
                   feature_plot_nF_pca = FeaturePlot(seuset, features=c("nFeature_RNA"),reduction="pca") + ggtitle(paste("nFeature RNA PCA - Res: ",cluster_selected_resolution,sep=''))
                   feature_plot_nF_umap = FeaturePlot(seuset, features=c("nFeature_RNA"),reduction="umap") + ggtitle(paste("nFeature RNA UMAP - Res: ",cluster_selected_resolution,sep=''))

                   complete_png_fileName = paste(results_documentation_directory,current_randomization_plus_maxPcDimension,"_cellMap_features.png",sep='')
                   png(complete_png_fileName,width=12000,height=9000,res=350);
                   grid.arrange(feature_plot_mt_pca, feature_plot_nC_pca,  feature_plot_nF_pca,
                                feature_plot_mt_tsn, feature_plot_nC_tsn,  feature_plot_nF_tsn,
                                feature_plot_mt_umap,feature_plot_nC_umap, feature_plot_nF_umap,nrow=3,ncol=3)
                   dev.off()
              }#End - Plot cell neighborhoods - quality control

              if ((continue_real_analysis)&(Generate_figures_for_visualization))
              {#Begin - Generate cell neighborhood plots by cluster and cell types
                   pt_size=0.01;
                   groups = c("Cluster","Cell_type_1st","Patient")
                   Plots = list();
                   for (indexG in 1:length(groups))
                   {#Begin
                      group = groups[indexG]
                      pca_plot = DimPlot(object = seuset, dims=c(1,2), reduction = "pca", pt.size=pt_size,label=TRUE,group.by = group)
                      pca_plot = pca_plot + ggtitle(paste("PCA - ",group," - Res: ",cluster_selected_resolution,sep=''))
                      Plots[[length(Plots)+1]] = pca_plot
                      tsne_plot = DimPlot(object = seuset, dims=c(1,2), reduction = "tsne", pt.size=pt_size,label=TRUE,group.by = group)
                      tsne_plot = tsne_plot + ggtitle(paste("TSNE - ",group," - Res: ",cluster_selected_resolution,sep=''))
                      Plots[[length(Plots)+1]] = tsne_plot
                      umap_plot = DimPlot(object = seuset, dims=c(1,2), reduction = "umap", pt.size=pt_size,label=TRUE,group.by = group)
                      umap_plot = umap_plot + ggtitle(paste("UMAP - ",group," - Res: ",cluster_selected_resolution,sep=''))
                      Plots[[length(Plots)+1]] = umap_plot
                   }#End
                   complete_png_fileName = paste(results_documentation_directory,current_randomization_plus_maxPcDimension,"_cellMaps.png",sep='')
                   png(complete_png_fileName,width=3000*3,height=3000*length(groups),res=350);
                   grid.arrange(Plots[[1]],Plots[[2]],Plots[[3]],Plots[[4]],Plots[[5]],Plots[[6]],Plots[[7]],Plots[[8]],Plots[[9]],nrow=3,ncol=length(groups));
                   dev.off()
              }#End - Generate cell neighborhood plots by cluster

              if (continue_real_analysis)
              {#Begin - Identify and Markers of each cellType versus all other cellType (almost copy paste from above, except cluster NO replaced by cell types)

                   findMarkersAssay = findMarkersAssays[indexFindMarkersAssay]
                   findMarkersSlot = findMarkersSlots[indexFindMarkersAssay]

                   {#Begin - Find markers using FindAllMarkers
                      combined_clusterDEGs <- tryCatch( { FindAllMarkers(object = seuset, min.pct = 0.1, only.pos = only_upregulated_genes, assay=findMarkersAssay, slot=findMarkersSlot) },
                                                        error=function(cond) { return (NA) },
                                                        finally = {});
                      if ((!exists("combined_clusterDEGs"))||(is.null(combined_clusterDEGs)))
                      { continue_real_analysis = FALSE }

                      if (continue_real_analysis)
                      {#Begin - Finalize combined_clusterDEGs
                         combined_clusterDEGs$cluster = as.character(combined_clusterDEGs$cluster)
                         missing_cellTypes = which(!cellTypes_to_be_detected %in% combined_clusterDEGs$cluster)

                         if (length(missing_cellTypes)>0)
                         {#Begin - Add noDEGs lines for missing clusters
                            Col_names = c("p_val","avg_log2FC","pct.1","pct.2","p_val_adj","gene","cluster")
                            Col_length = length(Col_names)
                            Row_names = missing_cellTypes
                            Row_length = length(Row_names)
                            add_cluster.markers = as.data.frame(array(NA,c(Row_length,Col_length),dimnames=list(Row_names,Col_names)))
                            add_cluster.markers$cluster = missing_cellTypes
                            add_cluster.markers$p_val = 1;
                            add_cluster.markers$avg_log2FC = 0;
                            add_cluster.markers$pct.1 = 0;
                            add_cluster.markers$pct.2 = 0;
                            add_cluster.markers$p_val_adj = 1;
                            add_cluster.markers$gene = "noDEGs"
                            combined_clusterDEGs = rbind(combined_clusterDEGs,add_cluster.markers)
                         }#End - Add noDEGs lines for missing clusters

                         combined_clusterDEGs$Cell_type = combined_clusterDEGs$cluster;
                         combined_clusterDEGs$Resolution = cluster_selected_resolution;
                         combined_clusterDEGs$Assay = findMarkersAssay;
                         combined_clusterDEGs$Slot = findMarkersSlot;
                         combined_clusterDEGs = as.data.frame(combined_clusterDEGs)
                         combined_clusterDEGs = as.data.frame(combined_clusterDEGs)
                         indexNa = which(is.na(combined_clusterDEGs$p_val))
                         if (length(indexNa)>0) { combined_clusterDEGs$p_val[indexNa] = 1; }
                         indexNa = which(is.na(combined_clusterDEGs$p_val_adj))
                         if (length(indexNa)>0) { combined_clusterDEGs$p_val_adj[indexNa] = 1; }
                         indexNa = which(is.na(combined_clusterDEGs$avg_logFC))
                         if (length(indexNa)>0) { combined_clusterDEGs$avg_logFC[indexNa] = 0; }
                      }#End - Finalize combined_clusterDEGs
                   }#End - Find markers using FindAllMarkers

              }#End - Identify DEGs of each cellType versus all other cellType (copy paste from above)

              if (continue_real_analysis)
              {#Begin - Write merged cluster DEGs
                   clusterDEGs_fileName = paste("CellTypeDEGsAfterMerging_",current_randomization_plus_maxPcDimension,"_",findMarkersAssay,"_",findMarkersSlot,".txt",sep='');
                   complete_clusterDEGs_fileName = paste(results_directory,clusterDEGs_fileName,sep='')
                   write.table(combined_clusterDEGs,file=complete_clusterDEGs_fileName,quote=FALSE,row.names=FALSE,sep='\t')
              }#End - Write merged cluster DEGs

              if (continue_real_analysis)
              {#Begin - Identify and write cellType patient counts and cellType_percentageCell_perPatient
                   cell_types = unique(barcode_cellType_associations$Cell_type)
                   patients = unique(seuset$Patient)
                   col_names = patients
                   col_length = length(col_names);
                   row_names = cell_types
                   row_length = length(row_names)
                   percentagePatient_perCellType = as.data.frame(array(0,c(row_length,col_length),dimnames=list(row_names,col_names)))

                   col_names = cell_types
                   col_length = length(col_names);
                   row_names = patients
                   row_length = length(row_names)
                   percentageCelltype_perPatient = as.data.frame(array(0,c(row_length,col_length),dimnames=list(row_names,col_names)))

                   col_names = c("CellType","Patient","Cells_count","Resolution")
                   col_length = length(col_names);
                   row_names = 1:(length(cell_types)*length(patients))
                   row_length = length(row_names)
                   cellCounts_per_CellType_and_patient = as.data.frame(array(0,c(row_length,col_length),dimnames=list(row_names,col_names)))

                   {#Begin - Calculate cell counts per patient and cell type combination
                      indexCellTypePatient = 0;
                      for (indexPatient in 1:length(patients))
                      {#Begin
                         patient = patients[indexPatient]
                         indexCurrentPatient = which(barcode_cellType_associations$Patient==patient)
                         for (indexCellType in 1:length(cell_types))
                         {#Begin
                            cell_type = cell_types[indexCellType]
                            indexCurrentCellType = which(barcode_cellType_associations$Cell_type==cell_type)
                            indexCurrentCellTypePatient = indexCurrentPatient[indexCurrentPatient %in% indexCurrentCellType]
                            indexCellTypePatient = indexCellTypePatient + 1;
                            cellCounts_per_CellType_and_patient$CellType[indexCellTypePatient] = cell_type;
                            cellCounts_per_CellType_and_patient$Patient[indexCellTypePatient] = patient;
                            cellCounts_per_CellType_and_patient$Cells_count[indexCellTypePatient] = length(indexCurrentCellTypePatient)
                         }#End
                      }#End
                      cellCounts_per_CellType_and_patient$Resolution = cluster_selected_resolution
                   }#End - Calculate cell counts per patient and cell type combination

                   {#Begin - Calculate percentagePatient per cell type
                      percentagePatient_perCellType_colnames = colnames(percentagePatient_perCellType)
                      for (indexCellType in 1:length(cell_types))
                      {#Begin
                         cell_type = cell_types[indexCellType];
                         indexCurrent = which(barcode_cellType_associations$Cell_type==cell_type);
                         current_barcode_cellType_associations = barcode_cellType_associations[indexCurrent,];
                         patient_table = table(current_barcode_cellType_associations$Patient);
                         current_patients = names(patient_table);
                         indexCurrentCellType = which(rownames(percentagePatient_perCellType)==cell_type)
                         for (indexPatient in 1:length(current_patients))
                         {#Begin
                             current_patient = current_patients[indexPatient];
                             indexCol = which(percentagePatient_perCellType_colnames==current_patient)
                             indexPatientTable = which(names(patient_table)==current_patient)
                             percentagePatient_perCellType[indexCurrentCellType,indexCol]=patient_table[indexPatientTable]
                         }#End
                         percentagePatient_perCellType[indexCurrentCellType,] = 100 * percentagePatient_perCellType[indexCurrentCellType,] / sum(percentagePatient_perCellType[indexCurrentCellType,])
                      }#End
                      percentagePatient_perCellType$RowSums = rowSums(percentagePatient_perCellType)
                      percentagePatient_perCellType$Resolution = cluster_selected_resolution;
                      percentagePatient_perCellType$CellType = rownames(percentagePatient_perCellType)
                      indexResolution = which(colnames(percentagePatient_perCellType)=="Resolution")
                      indexCellType = which(colnames(percentagePatient_perCellType)=="CellType")
                      indexRest = 1:length(percentagePatient_perCellType[1,])
                      indexRest = indexRest[!indexRest %in% c(indexCellType,indexResolution)]
                      percentagePatient_perCellType = percentagePatient_perCellType[,c(indexCellType,indexResolution,indexRest)]
                      percentagePatient_perCellType = percentagePatient_perCellType[order(percentagePatient_perCellType$CellType),]
                   }#End - Calculate percentagePatient per cell type

                   {#Begin - Calculate percentageCellType per patient
                      percentageCelltype_perPatient_colnames = colnames(percentageCelltype_perPatient)
                      for (indexPatient in 1:length(patients))
                      {#Begin
                         patient = patients[indexPatient]
                         indexCurrent = which(barcode_cellType_associations$Patient==patient)
                         current_barcode_cellType_associations = barcode_cellType_associations[indexCurrent,];
                         cellType_table = table(current_barcode_cellType_associations$Cell_type);
                         current_cellTypes = names(cellType_table);
                         indexCurrentPatient = which(rownames(percentageCelltype_perPatient)==patient)
                         for (indexCurrentCellType in 1:length(current_cellTypes))
                         {#Begin
                            current_cellType = current_cellTypes[indexCurrentCellType];
                            indexCol = which(percentageCelltype_perPatient_colnames==current_cellType)
                            indexCellTypeTable = which(names(cellType_table)==current_cellType)
                            if (length(indexCellTypeTable)>0)
                            {#Begin
                               percentageCelltype_perPatient[indexCurrentPatient,indexCol] = cellType_table[indexCellTypeTable]
                            }#End
                         }#End
                         percentageCelltype_perPatient[indexCurrentPatient,] = 100 * percentageCelltype_perPatient[indexCurrentPatient,] / sum(percentageCelltype_perPatient[indexCurrentPatient,])
                      }#End
                      percentageCelltype_perPatient$RowSums = rowSums(percentageCelltype_perPatient)
                      percentageCelltype_perPatient$Resolution = cluster_selected_resolution;
                      percentageCelltype_perPatient$Patient = rownames(percentageCelltype_perPatient)
                      indexResolution = which(colnames(percentageCelltype_perPatient)=="Resolution")
                      indexPatient = which(colnames(percentageCelltype_perPatient)=="Patient")
                      indexRest = 1:length(percentageCelltype_perPatient[1,])
                      indexRest = indexRest[!indexRest %in% c(indexPatient,indexResolution)]
                      percentageCelltype_perPatient = percentageCelltype_perPatient[,c(indexPatient,indexResolution,indexRest)]
                      percentageCelltype_perPatient = percentageCelltype_perPatient[order(percentageCelltype_perPatient$Patient),]
                   }#End - Calculate percentageCellType per patient

                   complete_percentagePatient_perCellType_fileName = paste(results_directory,"CellCounts_PercentPatient_perCellType_",current_randomization_plus_maxPcDimension,".txt",sep='')
                   write.table(percentagePatient_perCellType,file=complete_percentagePatient_perCellType_fileName,quote=FALSE,row.names=FALSE,sep='\t')

                   complete_percentageCelType_perPatient_fileName = paste(results_directory,"CellCounts_PercentCellType_perPatient_",current_randomization_plus_maxPcDimension,".txt",sep='')
                   write.table(percentageCelltype_perPatient,file=complete_percentageCelType_perPatient_fileName,quote=FALSE,row.names=FALSE,sep='\t')

                   complete_cellType_patientCounts = paste(results_directory,"CellCounts_per_cellType_and_patients_",current_randomization_plus_maxPcDimension,".txt",sep='')
                   write.table(cellCounts_per_CellType_and_patient,file=complete_cellType_patientCounts,quote=FALSE,row.names=FALSE,sep='\t')

                   if (Generate_figures_for_visualization)
                   {#Begin - plot cell type patient distributions
                      Stacked_bars = ggplot(cellCounts_per_CellType_and_patient,aes(x=CellType,y=Cells_count,fill=Patient))
                      Stacked_bars = Stacked_bars + geom_bar(stat="identity")
                      Stacked_bars = Stacked_bars + theme(plot.title = element_text(size = 40, face = "bold", hjust=0.5))
                      Stacked_bars = Stacked_bars + theme(axis.text = element_text(size = 25, face = "bold",angle=90))
                      Stacked_bars = Stacked_bars + theme(legend.text = element_text(size = 28, face = "bold",angle=0))
                      Stacked_bars = Stacked_bars + theme(axis.text.x = element_text(vjust=0.5))
                      Stacked_bars = Stacked_bars + theme(axis.title = element_text(size = 30, face = "bold"))

                      complete_png_fileName = paste(results_documentation_directory,current_randomization_plus_maxPcDimension,"_cellType_patient_cellCounts.png",sep='')
                      png(complete_png_fileName,width=9000,height=9000,res=350);
                      print(Stacked_bars)
                      dev.off()
                   }#End - plot cell type patient distributions

              }#End - Identify and write cluster cluster patient counts and cluster_percentageCell_perPatient

              if (continue_real_analysis)
              {#Begin - Identify and write cell type tissue collections
                   col_names = c("Cell_type","Tissue_collection","Cell_count","Per_cell_count","Resolution")
                   col_length = length(col_names);
                   row_names = 1
                   row_length = length(row_names)
                   cellType_tissue_collection_base_line = as.data.frame(array(0,c(row_length,col_length),dimnames=list(row_names,col_names)))
                   cellType_tissue_collections = c()

                   cellTypes = unique(seuset$Cell_type_1st)
                   for (indexCellTypes in 1:length(cellTypes))
                   {#Begin
                      cellType = cellTypes[indexCellTypes]
                      indexCurrent = which(seuset$Cell_type_1st==cellType);
                      current_tissue_collections = seuset$Tissue_collection[indexCurrent]
                      table_tc = table(current_tissue_collections)
                      total_cell_count = sum(table_tc)
                      for (indexT in 1:length(table_tc))
                      {#Begin
                         cellType_tissue_collection_line = cellType_tissue_collection_base_line
                         cellType_tissue_collection_line$Cell_type = cellType
                         cellType_tissue_collection_line$Tissue_collection = names(table_tc)[indexT]
                         cellType_tissue_collection_line$Cell_count = table_tc[indexT]
                         cellType_tissue_collection_line$Per_cell_count = 100 * cellType_tissue_collection_line$Cell_count / total_cell_count
                         cellType_tissue_collection_line$Resolution =  cluster_selected_resolution
                         if (length(cellType_tissue_collections)==0) { cellType_tissue_collections = cellType_tissue_collection_line }
                         else { cellType_tissue_collections = rbind(cellType_tissue_collections,cellType_tissue_collection_line) }
                      }#End
                   }#End

                   complete_cellType_tcs_fileName = paste(results_directory,"CellType_tissueCollection_",current_randomization_plus_maxPcDimension,".txt",sep='')
                   write.table(cellType_tissue_collections,file=complete_cellType_tcs_fileName,quote=FALSE,row.names=FALSE,sep='\t')
              }#End - Identify and write cluster tissue collections

              if (continue_real_analysis)
              {#Begin - Identify and write cellType tissue types
                   col_names = c("Cell_type","Tissue_type","Cell_count","Per_cell_count","Resolution")
                   col_length = length(col_names);
                   row_names = 1
                   row_length = length(row_names)
                   cellType_tissue_type_base_line = as.data.frame(array(0,c(row_length,col_length),dimnames=list(row_names,col_names)))
                   cellType_tissue_types = c()

                   cellTypes = unique(seuset$Cell_type_1st)
                   for (indexCellType in 1:length(cellTypes))
                   {#Begin
                      cellType = cellTypes[indexCellType]
                      indexCurrent = which(seuset$Cell_type_1st==cellType);
                      current_tissue_types = seuset$Tissue_type[indexCurrent]
                      table_tt = table(current_tissue_types)
                      total_cell_count = sum(table_tt)
                      for (indexT in 1:length(table_tt))
                      {#Begin
                         cellType_tissue_type_line = cellType_tissue_type_base_line
                         cellType_tissue_type_line$Cell_type = cellType
                         cellType_tissue_type_line$Tissue_type = names(table_tt)[indexT]
                         cellType_tissue_type_line$Cell_count = table_tt[indexT]
                         cellType_tissue_type_line$Per_cell_count = 100 * cellType_tissue_type_line$Cell_count / total_cell_count
                         cellType_tissue_type_line$Resolution = cluster_selected_resolution
                         if (length(cellType_tissue_types)==0) { cellType_tissue_types = cellType_tissue_type_line }
                         else { cellType_tissue_types = rbind(cellType_tissue_types,cellType_tissue_type_line) }
                      }#End
                   }#End

                   complete_cellType_tcs_fileName = paste(results_directory,"CellType_tissueType_",current_randomization_plus_maxPcDimension,".txt",sep='')
                   write.table(cellType_tissue_types,file=complete_cellType_tcs_fileName,quote=FALSE,row.names=FALSE,sep='\t')
              }#End - Identify and write cluster tissue types

              if (continue_real_analysis)
              {#Begin - Identify and write readCounts, mitochondrial gene expression, feature counts per cluster
                   unique_cellTypes = unique(seuset$Cell_type_1st)
                   length_cellTypes = length(unique_cellTypes)

                   col_names = c("Resolution","Cell_type","NUMI_RNA_average","NUMI_RNA_sd","NUMI_RNA_median","NFeature_RNA_average","NFeature_RNA_sd","NFeature_RNA_median","Percent_mt_average","Percent_mt_sd","Percent_mt_median","Cell_count")
                   col_length = length(col_names);
                   row_names = 1:length_cellTypes
                   row_length = length(row_names)
                   quality_control = as.data.frame(array(NA,c(row_length,col_length),dimnames=list(row_names,col_names)))

                   quality_control$Resolution = cluster_selected_resolution;

                   for (indexCellType in 1:length_cellTypes)
                   {#Begin
                      cell_type = unique_cellTypes[indexCellType]
                      indexCurrent = which(seuset$Cell_type_1st==cell_type)
                      current_nFeature_RNA = seuset$nFeature_RNA[indexCurrent]
                      current_nUMI_RNA = seuset$nCount_RNA[indexCurrent]
                      current_percent_mt = seuset$percent.mt[indexCurrent];
                      quality_control$Resolution = cluster_selected_resolution
                      quality_control$Cell_type[indexCellType] = as.character(cell_type);
                      quality_control$NUMI_RNA_average[indexCellType] = mean(current_nUMI_RNA);
                      quality_control$NUMI_RNA_median[indexCellType] = median(current_nUMI_RNA);
                      quality_control$NUMI_RNA_sd[indexCellType] = sd(current_nUMI_RNA);
                      quality_control$NFeature_RNA_average[indexCellType] = mean(current_nFeature_RNA);
                      quality_control$NFeature_RNA_median[indexCellType] = median(current_nFeature_RNA);
                      quality_control$NFeature_RNA_sd[indexCellType] = sd(current_nFeature_RNA);
                      quality_control$Percent_mt_average[indexCellType] = mean(current_percent_mt);
                      quality_control$Percent_mt_median[indexCellType] = median(current_percent_mt);
                      quality_control$Percent_mt_sd[indexCellType] = sd(current_percent_mt);
                      quality_control$Cell_count[indexCellType] = length(current_nUMI_RNA);
                   }#End

                   complete_quality_control_fileName = paste(results_directory,"Quality_control_",current_randomization_plus_maxPcDimension,".txt",sep='')
                   write.table(quality_control,file=complete_quality_control_fileName,quote=FALSE,row.names=FALSE,sep='\t')
              }#End - Identify and write readCounts, mitochondrial gene expression, feature counts per cluster

              if ((continue_real_analysis)&(Generate_figures_for_visualization))
              {#Begin - Generate QC metrics per cell type in violin plot
                  group = "Cell_type_1st"
                  vlnplot_nFeature_RNA = VlnPlot(object = seuset, features = c("nFeature_RNA"), ncol = 1, group.by=group)
                  vlnplot_nFeature_RNA = vlnplot_nFeature_RNA + geom_hline(yintercept=minimum_feature_count_per_cell,col="green",size=1)
                  vlnplot_nFeature_RNA = vlnplot_nFeature_RNA + geom_hline(yintercept=maximum_feature_count_per_cell,col="green",size=1)
                  vlnplot_nCount_RNA = VlnPlot(object = seuset, features = c("nCount_RNA"), ncol = 1, group.by=group)
                  vlnplot_percent_mt = VlnPlot(object = seuset, features = c("percent.mt"), ncol = 1, group.by=group)
                  vlnplot_percent_mt = vlnplot_percent_mt + geom_hline(yintercept=max_percentage_mitochondrial_genes,col="green",size=1)

                  complete_png_fileName = paste(results_documentation_directory,current_randomization_plus_maxPcDimension,"_cellType_quality_plots.png",sep='')
                  png(complete_png_fileName,width=9000,height=9000,res=350);
                  grid.arrange(vlnplot_nFeature_RNA,vlnplot_nCount_RNA,vlnplot_percent_mt,nrow=3,ncol=1)
                  dev.off()
              }#End - Generate QC metrics per cell type in violin plot

          }#End - Perform real analysis at last resolution that detected a new cell type

          previous_cluster_selected_resolution = cluster_selected_resolution

        }#End - indexMaxPcDimension
       }#End  - if (continue_randomizations_loop)

       if (exists("seuset")) { rm(seuset) }
       if (exists("barcode_cluster_associations")) { rm(barcode_cluster_associations) }
       if (exists("combined_clusterDEGs")) { rm(combined_clusterDEGs) }
       gc()

        endTime_of_current_randomizationNo = Sys.time()

        Col_names = c("Core","IndexRando","Core_randomizations","Started","Finished","Total_cells_count")
        Col_length = length(Col_names)
        Row_names = 1
        Row_length = length(Row_names)
        current_report = array(NA,c(Row_length,Col_length),dimnames = list(Row_names,Col_names))
        current_report = as.data.frame(current_report)
        current_report$Core = indexCore
        current_report$IndexRando = indexRando
        current_report$Started = startTime_of_current_randomizationNo
        current_report$Finished = endTime_of_current_randomizationNo
        current_report$Total_cells_count = length(Raw_data_matrix[1,])

        complete_analysisFinished_fileName = paste(results_directory,"AA_Analysis_finished_set",set,"_core",indexCore,"_randoNo",current_randomization_no,".txt",sep='')
        write.table(current_report,file=complete_analysisFinished_fileName,quote=FALSE,col.names=TRUE,row.names=FALSE,sep='\t')
        if (file.exists(complete_analysisStarted_fileName)) { file.remove(complete_analysisStarted_fileName) }
    }#End  for loop

    complete_coreFinished_fileName = paste(results_documentation_directory,"AAA_Analysis_finished_core",indexCore,".txt",sep='')
    write.table(paste("Core ",indexCore," DONE",sep=''),file=complete_coreFinished_fileName,quote=FALSE,row.names=FALSE,sep='')
}#End - Parallel foreach
)


registerDoSEQ()
unregister <- function() {
   env <- foreach:::.foreachGlobals
   rm(list=ls(name=env), pos=env)
}
unregister()
stopCluster(parallel_clusters)
rm(parallel_clusters)
gc()


# Change back to start-up directory
if(!is.null(start_dir)) setwd(start_dir)
