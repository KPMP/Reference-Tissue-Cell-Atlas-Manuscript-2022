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

#Open libraries and document versions - BEGIN
Col_names = c("Library","Version")
Col_length = length(Col_names)
Row_names = 1
Row_length= length(Row_names)
version_documentation_line = array(NA,c(Row_length,Col_length),dimnames=list(Row_names,Col_names))
version_documentation_line = as.data.frame(version_documentation_line)
version_documentations = c()

libraries = c("Rtsne","ClassDiscovery","colorspace","dendextend","circlize","ape","fields","colormap","gplots")

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
#Open libraries and document versions - END

base_directory = paste(getwd(),"/",sep='');#
base_directory = "D:/KPMP_reference_atlas_code/"
directory = paste(base_directory,"/Results/GlomVsProxTub/",sep='');

cellType_podGlom_color_rnaSeq = "aquamarine3";
cellType_PTTI_color_rnaSeq = "darkorange3"
cellType_podGlom_color_proteomics = "aquamarine4";
cellType_PTTI_color_proteomics = "darkorange4"
cellType_PCCD_color_rnaSeq = "darkorchid3";
cellType_other_color_rnaSeq = "black"

correlation_methods = c("pearson");#"spearman"#

###Write version documentations - BEGIN
results_documentation_directory = paste(directory,"Documentation/",sep='')
dir.create(results_documentation_directory)

version_control_fileName = "AA_version_documentation.txt"
r_session_fileName = "AA_R_session.txt"
complete_version_control_fileName = paste(results_documentation_directory,version_control_fileName,sep='')
complete_r_session_fileName = paste(results_documentation_directory,r_session_fileName,sep='')

sink(file=complete_r_session_fileName)
r_sessionInfo
sink()

write.table(version_documentations,file=complete_version_control_fileName,quote=FALSE,row.names=FALSE,col.names=TRUE,sep='\t')
###Write version control - END


indexC=1
for (indexC in 1:length(correlation_methods))
{#Begin

  correlation_method = correlation_methods[indexC];

  datasets = c(
    "GlomVsProxTub_Single_value_all_degsdeps",
    "GlomVsProxTub_Single_value_all_shared_degsdeps",
    "GlomVsProxTub_Log2_ratio_singlepatient_all_degsdeps",
    "GlomVsProxTub_Log2_ratio_singlepatient_all_shared_degsdeps"
  );

  #pdf_complete_fileName = paste(directory,"KPMP_subsegmental_correlation_",correlation_method,".pdf",sep='');
  #pdf(pdf_complete_fileName, width=8.5, height=11);

  indexChange_larger_one=0;

  indexF=4
  for (indexF in 1:length(datasets))
  {#Begin
    dataset = datasets[indexF]
    fileName = paste(dataset,".txt",sep='');
    complete_fileName = paste(directory,fileName,sep='');

    Data = read.table(file=complete_fileName,stringsAsFactors = FALSE, sep="\t", quote="",header=TRUE)

    Colnames = colnames(Data)
    indexSymbol = which(Colnames=="NCBI_symbol")
    if (length(indexSymbol)!=0)
    {
      Colnames[indexSymbol] = "Symbol"
    }

    Colnames = gsub(" ",".",Colnames)
    Colnames = gsub("Rna.","",Colnames)
    Colnames = gsub(".vs.Glom","",Colnames)
    Colnames = gsub(".vs.Pod","",Colnames)
    Colnames = gsub(".vs.POD","",Colnames)
    Colnames = gsub(".vs.ProxTub","",Colnames)
    Colnames = gsub("counts","",Colnames)
    Colnames = gsub(".Log2.ratio.singlepatient.Podocyte.glomerulus","",Colnames)
    Colnames = gsub(".Log2.ratio.singlepatient.Proximal.tubule","",Colnames)
    Colnames = gsub("_Log2_ratio_singlepatient",".",Colnames)
    Colnames = gsub("[.]"," ",Colnames)
    Colnames = gsub("_"," ",Colnames)
    Colnames = gsub("LMD Proteomics OSUIU ProxTub","LMD Proteomics OSUIU TI",Colnames)
    Colnames = gsub(" POD"," Pod",Colnames)

    colnames(Data) = Colnames

    indexKeep = which(!is.na(Data$Symbol))
    Data = Data[indexKeep,]

    rownames(Data) = Data$Symbol;

    colnames_remove = c("Symbol","RowSum","RowMean","RowSampleSD","RowSampleCV","Description","Gene_ontology","Molecular Biology of the Cell","MGI","KEGG","Own","Transfac","Various","Human_symbols","NCB Idescription","NCBI geneID","SCPs","Human symbols","Gene ontology","RowSum","RowMean","RowSampleSD","RowSampleCV","Description","ReadWrite SCPs","ReadWrite human symbols")
    indexColKeep = 1:length(Data[1,])
    indexColRemove = which(colnames(Data) %in% colnames_remove)
    indexColKeep = indexColKeep[!indexColKeep %in% indexColRemove]
    Data = Data[,indexColKeep];
    Colnames = colnames(Data)
    Rownames = rownames(Data);

    add_to_file = "";
    if (grepl("GlomVsProxTub_Single_value",dataset))
    {#Begin
       if (min(Data)==0)
       {#Begin
           Data = Data + 1
       }#End
       Data = log10(Data)
       add_to_file = "log10_";
    }#End


    dist.function = function(x) as.dist((1-cor((x),method=correlation_method))/2)
    hclust.function = function(x) hclust(x,method="average")

    col_dist = dist.function(Data)
    col_hc = hclust.function(col_dist)
    col_dend = as.dendrogram(col_hc)
    indexOriginalData_col = order.dendrogram(col_dend);
    col_labels_in_dendrogram = Colnames[indexOriginalData_col]
    Data_reordered = Data[,indexOriginalData_col]

    row_dist = dist.function(t(Data) )
    row_hc = hclust.function(row_dist)
    row_dend = as.dendrogram(row_hc)
    indexOriginalData_row = order.dendrogram(row_dend);
    row_labels_in_dendrogram = Rownames[indexOriginalData_row]
    Data_reordered = Data_reordered[rev(indexOriginalData_row),] #Data_reordered = Data_reordered[rev(indexOriginalData_row),] this would be correct, so rotated dendrogram 180 degrees

    Colors = replicate(length(col_labels_in_dendrogram),"gray");

    indexPodGlom1 = grep("Pod",col_labels_in_dendrogram);
    indexProxTub1 = grep("ProxTub",col_labels_in_dendrogram);
    indexPodGlom2 = grep("Glom",col_labels_in_dendrogram);
    indexProxTub2 = grep(" TI ",col_labels_in_dendrogram);

    indexRNASeq = grep("RNASeq",col_labels_in_dendrogram)
    indexProteomics = grep("Proteomics",col_labels_in_dendrogram)
    indexProxTub = c(indexProxTub1,indexProxTub2)
    indexPodGlom = c(indexPodGlom1,indexPodGlom2)

    indexPodGlom_prot = indexPodGlom[indexPodGlom %in% indexProteomics]
    indexPodGlom_rnaSeq = indexPodGlom[indexPodGlom %in% indexRNASeq]

    indexProxTub_prot = indexProxTub[indexProxTub %in% indexProteomics]
    indexProxTub_rnaSeq = indexProxTub[indexProxTub %in% indexRNASeq]

    Colors[indexProxTub_rnaSeq] = cellType_PTTI_color_rnaSeq
    Colors[indexPodGlom_rnaSeq] = cellType_podGlom_color_rnaSeq
    Colors[indexProxTub_prot] = cellType_PTTI_color_proteomics
    Colors[indexPodGlom_prot] = cellType_podGlom_color_proteomics

    col_dend = color_labels(col_dend,col=Colors)

    complete_pdf_fileName = paste(directory,"KPMP_subsegmental_", "correlation_",correlation_method,"_",add_to_file,dataset,"_colDendogram.pdf",sep='')
    pdf(complete_pdf_fileName,height=11.7,width=5);


    par(cex=0.8,font=2);
    #circlize_dendrogram(dend)
    par(cex=0.2,font=1);
    par(mar=c(1,1,1,80)+0.1)


    col_dend <- set(col_dend, "labels_cex", 3)
    plot(col_dend,horiz=TRUE);
    title(main=dataset)

    dev.off();

    complete_pdf_fileName = paste(directory,"KPMP_subsegmental_","correlation_",correlation_method,"_",add_to_file,dataset,"_matrix.pdf",sep='')
    pdf(complete_pdf_fileName,height=11.7,width=8.3);

    cool = rainbow(50, start=rgb2hsv(col2rgb('cyan'))[1], end=rgb2hsv(col2rgb('blue'))[1])
    warm = rainbow(50, start=rgb2hsv(col2rgb('red'))[1], end=rgb2hsv(col2rgb('orange'))[1])
    cols = c(rev(cool), '#ffffff',rev(warm))
    if (min(Data_reordered)<0)
    {#Begin
      Color_map <- colorRampPalette(cols)(255)
    }#End
    if (min(Data_reordered)>=0)
    {#Begin
      Color_map <- colorRampPalette(c("#ffffff",rev(warm)))(255)
    }#End


    par(mar=c(7,4,2,5)+0.1)
    image((as.matrix(Data_reordered))[length(Data_reordered[,1]):1,],axes=FALSE, col= Color_map);
    image.plot(t(Data_reordered),legend.only=TRUE, col= Color_map);
    box(col="black",lty="solid")
    dimensions = dim(Data_reordered)
    if (length(row_labels_in_dendrogram)<=100) { grid(dimensions[1],dimensions[2],col="black",lty="solid",lwd=0.02) }
    dev.off();

    get_branches_heights=FALSE
    if (get_branches_heights)
    {#Begin - Get branches heights
      branches_heights = unique(get_branches_heights(col_dend,sort=TRUE))
      top_heights = rev(branches_heights[(length(branches_heights)-7):length(branches_heights)])
      indexTop=8
      plot(col_dend,horiz=TRUE,main=top_heights[indexTop])
      abline(v=top_heights[indexTop],col="red",lty=2)

    }#End - Get branches heights
  }#End
}#End


# Change back to start-up directory
if(!is.null(start_dir)) setwd(start_dir)
