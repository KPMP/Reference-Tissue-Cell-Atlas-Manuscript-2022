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

#Used for paper
library(gridExtra)
library(grid)
library(ggplot2)
unlink(".RData")
rm(list = ls(all.names=TRUE));

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

options(future.globals.maxSize= 20000000*1024^2)

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
            cat(paste("Usage:", prog_name, "[tis_center]\n"))
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
   tis_center <- getArg(invoke_method=invoke_method, arg_name="tis_center", arg_pos=1, default_value="arg1_default")
   if(tis_center!="Premiere" && tis_center!="UCSD") stop("tis_center must be either Premiere or UCSD!")
}#End - Extract command line input arguments - written by Yuguang Xiong

{#Begin - Predefinitions

base_directory = paste(getwd(),"/",sep='');#"D:/KPMP_reference_atlas_code/"
base_directory = "D:/KPMP_reference_atlas_code/"
tis_center = "Premiere"
tis_center = "UCSD"

if (tis_center=="Premiere")
{#Begin - Analyze Premiere settings
directory = paste(base_directory,"Results/PostHocPower_Premiere_AAAS_2022February14_finalResults/",sep='')
baseName = "Premiere_AAAS_2022February14_mPC30_SCT_data_300DEGsCT_300DEGsSCP"
considered_cell_types = c("POD","EC-GLO","PT","DTL","ATL","TAL","DCT","PC","IC","INT","EC-AVR","EC-AEA-DVR","MAC","MON","Tcell")
}#End - Analyze Premiere settings

if (tis_center=="UCSD")
{#Begin - #Analyze UCSD settings
directory = paste(base_directory,"Results/PostHocPower_UCSD_AAAS_2022February14_finalResults/",sep='')
baseName = "UCSD_AAAS_2022February14_mPC30_SCT_data_300DEGsCT_300DEGsSCP";
considered_cell_types = c("POD","MC","EC-GLO","PT","PT/DTL","DTL","ATL","TAL","DCT","CNT","PC","IC","INT/vSMC/P","EC-AVR")
}#End - #Analyze UCSD settings

fileName = paste(baseName,".txt",sep='')
pdf_fileName = paste(baseName,".pdf",sep='')
pdf_suppl_fileName = paste(baseName,"_suppl.pdf",sep='')
complete_fileName = paste(directory,fileName,sep='')
complete_pdf_fileName = paste(directory,pdf_fileName,sep='')

LeavePatientsOut = read.csv(complete_fileName,header=TRUE, stringsAsFactors = FALSE, sep='\t')
LeavePatientsOut = LeavePatientsOut[LeavePatientsOut$Considered_patients>1,]

correlation_cutoff_line = 0.8;
cutoff_line_width = 1;
standardEnrichment_cutoff_line = 20;
dynamicEnrichment_cutoff_line = 14;
patient_count_cutoff_line1 = 99; #15;
patient_count_cutoff_line2 = 99; #25;
cutoff_line_type = "dashed"
cutoff_line_color = "black"
set_line_types = c("solid","dashed")
axis_color = "gray25"

entityClassNames = unique(LeavePatientsOut$EntityClassName)
indexFisher = grep("Fisher exact - ",entityClassNames)
cellTypes = entityClassNames[indexFisher]
cellTypes = gsub("Fisher exact - ","",cellTypes)
indexKeep = which(!cellTypes %in% c("Not assigned","NKC",1:30))
cellTypes = cellTypes[indexKeep]

cellTypes_abbr_list = list("POD" = "POD",
                           "MC" = "MC",
                           "MC/vSMC/P" = "MC/vSMC/P",
                           "EC-GLO" = "GLO EC",
                           "PT" = "PT",
                           "PT/DTL" = "PT/DTL",
                           "DTL" = "DTL",
                           "ATL" = "ATL",
                           "ATL/TAL" = "ATL/TAL",
                           "TAL" = "TAL",
                           "DCT" = "DCT",
                           "CNT" = "CNT",
                           "CNT/PC" = "CNT/PC",
                           "PC" = "PC",
                           "PC/IC" = "PC/IC",
                           "IC" = "IC",
                           "INT" = "FIB",
                           "EC-AVR" = "AVR EC",
                           "EC-AEA-DVR" = "AEA-DVR EC",
                           "vSMC/P" = "vSMC/P",
                           "MAC" = "MAC",
                           "MON" = "MON",
                           "Tcell" = "T cell",
                           "Bcell" = "B cell",
                           "NKC" = "NKC",
                           "Multiple assignments" = "Multiple assignments",
                           "CNT/EC-GLO" = "CNT/EC-GLO",
                           "CNT/IC" = "CNT/IC",
                           "DTL/DCT" = "DTL/DCT",
                           "DTL/TAL" = "DTL/TAL",
                           "EC-GLO/EC-AEA-DVR" = "EC-GLO/EC-AEA-DVR",
                           "MC/EC-GLO" = "MC/EC-GLO",
                           "PC/vSMC/P" = "PC/vSMC/P",
                           "PT/CNT" = "PT/CNT",
                           "PT/EC-GLO" = "PT/EC-GLO",
                           "TAL/DCT" = "TAL/DCT",
                           "DCT/CNT" = "DCT/CNT",
                           "DTL/CNT" = "DTL/CNT",
                           "DTL/INT" = "DTL/INT",
                           "EC-AEA-DVR/INT" = "EC-AEA-DVR/INT",
                           "EC-AVR/EC-AEA-DVR" = "EC-AVR/EC-AEA-DVR",
                           "EC-GLO/EC-AVR" = "EC-GLO/EC-AVR",
                           "PT/EC-AVR" = "PT/EC-AVR",
                           "PT/PC" = "PT/PC",
                           "PT/vSMC/P" = "PT/vSMC/P",
                           "ATL/CNT" = "ATL/CNT",
                           "ATL/PC" = "ATL/PC",
                           "DTL/PC" = "DTL/PC",
                           "INT/EC-GLO" = "INT/EC-GLO",
                           "INT/vSMC/P" = "INT/vSMC/P",
                           "PT/TAL" = "PT/TAL",
                           "POD/CNT" = "POD/CNT",
                           "PT/DCT" = "PT/DCT",
                           "PT/INT" = "PT/INT",
                           "DCT/EC-GLO" = "DCT/EC-GLO",
                           "POD/PT" = "POD/PT",
                           "TAL/PC" = "TAL/PC",
                           "CNT/vSMC/P" = "CNT/vSMC/P",
                           "POD/INT" = "POD/INT",
                           "TAL/CNT" = "TAL/CNT",
                           "TAL/vSMC/P" = "TAL/vSMC/P",
                           "DTL/EC-AEA-DVR" = "DTL/EC-AEA-DVR",
                           "ATL/vSMC/P" = "ATL/vSMC/P",
                           "CNT/INT" = "CNT/INT",
                           "DTL/vSMC/P" = "DTL/vSMC/P",
                           "POD/DCT" = "POD/DCT",
                           "PT/ATL" = "PT/ATL",
                           "DCT/IC" = "DCT/IC",
                           "TAL/MAC" = "TAL/MAC",
                           "MC/DTL" = "MC/DTL",
                           "DTL/IC" = "DTL/IC"
)
indexNotIdentical = which(names(cellTypes_abbr_list)!=cellTypes_abbr_list)
indexMissing = which(!considered_cell_types %in% names(cellTypes_abbr_list))
if (length(indexMissing)>0) { rm(LeavePatientsOut)}
indexKeep = which(names(cellTypes_abbr_list)%in%considered_cell_types)
cellTypes_abbr_list = cellTypes_abbr_list[indexKeep]

axis_lwd = 0.2;
main_datapoints_pch = 16;
second_main_datapoints_pch = 1
lwd=second_main_datapoints_lwd = 0.1
main_datapoints_cex = 0.6;


plot_margin_bottom=1.7#2;
plot_margin_left=4;
plot_margin_top=1.7;
plot_margin_right=0.7;

min_x = 0;
max_x = ceiling(max(LeavePatientsOut$Considered_patients) / 5) * 5

cd_color = "darkorchid3"
henle_color = "dodgerblue3"
glom_color = "aquamarine4"
tubule_color = "darkorange3"
immune_color = "lightgoldenrod4"
endothelial_color = "indianred3"
interstitial_color = "gray40"


cellType_color_list = list( "POD" = glom_color
                           ,"PT" = tubule_color
                           ,"PT/DTL" = tubule_color
                           ,"ATL" = henle_color
                           ,"ATL/TAL" = henle_color
                           ,"TAL" = henle_color
                           ,"CNT/PC" = cd_color
                           ,"PC" = cd_color
                           ,"MC" = glom_color
                           ,"MC/vSMC/P" = glom_color
                           ,"IC" = cd_color
                           ,"DCT" = tubule_color
                           ,"EC-GLO" = glom_color
                           ,"EC-AVR" = endothelial_color
                           ,"EC-AEA-DVR" = endothelial_color
                           ,"INT/vSMC/P" = endothelial_color
                           ,"vSMC/P" = endothelial_color
                           ,"MAC" = immune_color
                           ,"MON" = immune_color
                           ,"Bcell" = immune_color
                           ,"Bcell/Tcell" = immune_color
                           ,"Tcell" = immune_color
                           ,"INT" = immune_color
                           ,"DTL" = henle_color
                           ,"CNT" = cd_color
                           ,"Multiple assignments" = "black")

cellTypes = cellTypes[cellTypes != "Multiple assignments"]
cellTypes = names(cellTypes_abbr_list)

}#End - Predefinitions

pdf(complete_pdf_fileName, width=8.5, height=11);
par(mfrow=c(2,2))

base_array = c( 1, 2, 3, 4,
                1, 2, 3, 4,
                5, 6, 7, 8,
                5, 6, 7, 8,
                9,10,11,12,
                9,10,11,12,
               13,14,15,16,
               13,14,15,16,
               17,18,19,20,
               17,18,19,20)
matrix_array = c(base_array,base_array+max(base_array))
layout_matrix = matrix(matrix_array,nrow=2*10,ncol=4,byrow=TRUE)

figures_per_set = 16;
layout(layout_matrix);
par(mar=c(2,2,2,2))

indexFigurePage = 0;

for(indexCellType in 1:length(cellTypes))
{#Begin
   cellType = cellTypes[indexCellType]
   cellType_string1 = paste(" - ",cellType,"_",sep='')
   cellType_string2 = paste(" - ",cellType," vs",sep='')
   current_color = "gray"
   if (cellType %in% names(cellType_color_list))
   { current_color = cellType_color_list[[cellType]] }

   indexCurrentCellType1 = grep(cellType_string1,paste(LeavePatientsOut$EntityClassName,"_",sep=''))
   indexCurrentCellType2 = grep(cellType_string2,LeavePatientsOut$EntityClassName)
   indexCurrentCellType = unique(c(indexCurrentCellType1,indexCurrentCellType2))

   if (length(indexCurrentCellType)>0)
   {#Begin --- all plot of current cell type
      indexFigurePage = indexFigurePage+1;
      if (indexFigurePage>max(layout_matrix)) { indexFigurePage=1; }
      currentCellType_leavePatientsOut = LeavePatientsOut[indexCurrentCellType,]

      ##### Cell type detected - Start
      indexCurrentEntityClass = grep("Fisher exact - ",currentCellType_leavePatientsOut$EntityClassName)
      if (length(indexCurrentEntityClass)==0)
      {#Begin
         plot(indexFigurePage,xlim=c(0,1),ylim=c(0,1),ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n');
      }#End
      else
      {#Begin
         current_fisher = currentCellType_leavePatientsOut[indexCurrentEntityClass,]
         current_fisher_set0 = current_fisher[current_fisher$Set_no==0,]
         current_fisher_set1 = current_fisher[current_fisher$Set_no==1,]
         plot_data = current_fisher_set1;

         plot_data = plot_data[order(plot_data$Considered_patients,decreasing=TRUE),]

         indexesAbove95 = which(plot_data$Percent_randomizationSets_count>=95)
         indexAbove95 = 0
         if (length(indexesAbove95)>0) { indexAbove95 = max(indexesAbove95) }

         if (indexAbove95 > 0)
         { cellTypeDetected_cutoff = plot_data$Considered_patients[indexAbove95] }

         indexFigurePage = indexFigurePage + 1;
         min_y = 0;
         max_y = 100;
         main_title = cellTypes_abbr_list[[cellType]]

         par(mar=c(plot_margin_bottom,plot_margin_left,plot_margin_top,plot_margin_right))
         plot(indexFigurePage, type="n",ylim=c(min_y,max_y), xlim=c(min_x,max_x),frame=FALSE,axes=FALSE,ann=FALSE)
         title(main=main_title,cex.main=1.5)
         axis(1,lwd=0.3,col=axis_color,cex.axis=1.3)
         axis(2,lwd=0.3,col=axis_color,cex.axis=1.3)

         indexSet=1;
         points(plot_data$Considered_patients,plot_data$Percent_randomizationSets_count,type="p",col=current_color,pch=main_datapoints_pch,cex=0.8)

         if (indexAbove95>0)
         { points(c(cellTypeDetected_cutoff,cellTypeDetected_cutoff),c(0,max_y),type="l",col=cutoff_line_color,lty=cutoff_line_type,lwd=1.2) }
         points(c(0,max_x),c(95,95),type="l",col=cutoff_line_color,lty=cutoff_line_type,lwd=1.2)

         par(mar=c(plot_margin_bottom,plot_margin_left+0.5,plot_margin_top,plot_margin_right))
         title(ylab="[%]",cex.lab=1.3)
         points(c(patient_count_cutoff_line1,patient_count_cutoff_line1),c(min_y,max_y),type="l",col=axis_color,lwd=0.3,lty=cutoff_line_type)
         points(c(patient_count_cutoff_line2,patient_count_cutoff_line2),c(min_y,max_y),type="l",col=axis_color,lwd=0.3,lty=cutoff_line_type)
      }#End
      ##### Cell type detected - END
   }#End --- all plot of current cell type

}#End - Main figure

dev.off()


# Change back to start-up directory
if(!is.null(start_dir)) setwd(start_dir)
