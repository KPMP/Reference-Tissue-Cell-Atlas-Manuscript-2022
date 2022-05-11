All code files and datasets discussed below should be copied into the listed directories that are within the
base directory "KPMP_reference_atlas_code/".

Download datasets:

-----------------------------
KPMP experimental datasets:
https://www.kpmp.org/doi-collection/10-48698-z30t-0a62

Copy all subdirectories into the subdirectory "/Experimental_data/"
-----------------------------



-----------------------------
Molecular Biology of the Cell Ontology (www.mbc-ontology.org):
link to website: "https://github.com/SBCNY/Molecular-Biology-of-the-Cell/MBCO_windows_application/MBCO_datasets/"
Download the files:
"MBCO_v1.1_inferred_SCP_relationships.txt"
"MBCO_v1.1_gene-SCP_associations_human.txt"
Copy these files into subdirectory: "/MBCO_datasets/"
-----------------------------


-----------------------------
Gene Ontology Biological Processes:
link to website: "https://maayanlab.cloud/Enrichr/#libraries"
Download the file:
"GO_Biological_Process_2018"
(We downloaded this file on 2018 July 17.)
Copy these files into subdirectory: "/GeneOntology_datasets/"
Rename to "Geneontology_BiologicalProcess_2018_2018July17.txt"
-----------------------------


-----------------------------
Metabolic and transmembrane transport ontologies
link to website: "https://github.com/SBCNY/Molecular-Biology-of-the-Cell/Additional_MBCO_gene_libraries/"
Download the files:
"Metabolic_energy_generation_pathways_human.txt"
"NaAndGluTMTransport_scpGeneAssociations.txt"
"Parent_child_transmembrane_transport_scps.txt"                                        
"Sphingolipid_metabolism.txt"
Copy these files into subdirectory: "/MBCO_datasets/"
-----------------------------


-----------------------------
jensenLab_compartments
link to website: "https://compartments.jensenlab.org/Downloads"
Download "All channels integrated: human"
(We downloaded this file on 2020 February 21.)
Copy this file into subdirectory "/jensenLab_compartments/"
Rename downloaded file "human_compartment_integrated_full.tsv" to "human_compartment_integrated_full_2020February21.tsv"
-----------------------------


-----------------------------
Single nucleus dataset by Wu et al. PMID:29980650 
Download "GSE114156_Human.kidney.dge.txt.gz" from NCBI GEO website - GSE114156
(The last time we downloaded this file was 2022 March 17.)
Copy file into subdirectory "/Experimental_data/SingleCell_datasets/"
Rename file to "Wu_MTS.kidney_PMID29980650.dge.txt"
-----------------------------


-----------------------------
Single nucleus dataset by Muto et al. 2021 Apr 13;12(1):2190. doi: 10.1038/s41467-021-22368-w. PMID:33850129 
Link to website: https://cellxgene.cziscience.com/collections/9b02383a-9358-4f0f-9795-a891ec523bcc
Download the dataset "Single cell transcriptional and chromatin accessibility profiling redefine cellular heterogeneity in the adult human kidney - RNAseq" in ".rds" format
(Our last download of this file was on 2022 March 17.)
Copy file into subdirectory "/Experimental_data/SingleCell_datasets/"
Rename file to "Muto_singleNucleus_seurat_PMID33850129.rds"
-----------------------------

--------------------------------------------------------------------------------------------------------------------------------------------
Run R- and C#-scripts

Go to the directory "KPMP_reference_atlas_code/" that was generated after extracting the downloaded KPMP_reference_atlas_code.zip file.

-----------------------------
Set 1 - Run these R-scripts first
Go into the subdirectory = "KPMP_reference_atlas_rCode/"
Run r-scripts by using a suitable editor such as R-studoe (order does not matter):
RunBeforeCSharpScript_PostHocPower_singleNucleusCell_analysis.R
RunBeforeCSharpScript_calculate_meanExpressionValues_in_integrated_clusterAnalysis.R

Notes: 1) In each script the base directory needs to be specified.
          Default is: base_directory = "D:/KPMP_reference_atlas_code/"
       2) In the PostHocPower analysis script, the TIS center needs to be specified.
          Default is: "Premiere", alternative "UCSD". The following C# script will request the results
          obtained for both TISs
-----------------------------


-----------------------------
Set 2 - Run c-sharp script after running all scripts in set 1
Go into the subdirectory = "KPMP_reference_atlas_csharp_code\"
Open C# script using a suitable editor such as Microsoft Visual Studio
KPMP_reference_kidney_tissue_atlas.sln
Notes: 1) The script has to be within the directory "KPMP_reference_atlas_csharp_code\". Otherwise it won't find the input and output directories.
       2) A progress report will be written as a text-file ("Report.txt") that can be found in the subdirectory "Results\".
-----------------------------



-----------------------------
Set 3 - Run these R-scripts after running c-sharp script of set 2
Go into the subdirectory = "KPMP_reference_atlas_rCode/"
Run r-scripts by using a suitable editor such as R-studio (order does not matter):
RunAfterCSharpScript_metabolicEnergyPathways_all.R
RunAfterCSharpScript_metabolicEnergyPathways_averaged.R
RunAfterCSharpScript_PostHocPower_visualizeResults.R 
RunAfterCSharpScript_sphingomyosin_metabolism.R
RunAfterCSharpScript_crossOmicsComparison_correlationCoefficients.R
RunAfterCSharpScript_crossOmicsComparison_hierarchicalClustering.R

Notes: 1) In each script the base directory needs to be specified.
          Default is: base_directory = "D:/KPMP_reference_atlas_code/"
       2) Within the PostHocPower visualization scripts the center has to be specified: "Premiere" or "UCSD"
       3) The last script is time-consuming.
-----------------------------


-----------------------------
Set 4 - Run this R-scripts in the given order (they are independent of the scripts in set 1,2 and 3)
Go into the subdirectory = "KPMP_reference_atlas_rCode/"
Run r-scripts by using a suitable editor such as R-studio in the given order:
Rscript IndependentOfCsharpScript_no1_Wu_PMID29980650_Single_cell_analysis.R
Rscript IndependentOfCsharpScript_no2_SCP_readCounts_along_nephron.R

Notes: 1) In each script the base directory needs to be specified.
          Default is: base_directory = "D:/KPMP_reference_atlas_code/"
       2) The second r-script needs access to the internet, because it downloads ensembl gene symbol annotations during runtime
          (using the command 'useDataset' from the 'biomaRt' package). Our results are based on the ensembl annotation downloaded
          on 2022 March 25.
-----------------------------
