/*
Copyright 2022.The code was written by Jens Hansen working for the Ravi Iyengar Lab
and the Kidney Precision Medicine Project (KPMP)
It is made availabe under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

Please see http://www.mbc-ontology.org/ or https://github.com/SBCNY/Molecular-Biology-of-the-Cell for a windows application offering the same and additional functionalities.
The script below is extracted from that application.
Please acknowledge our publication by citing the following reference:


Please check out www.mbc-ontology.org for a windows application that allows analysis of custom data sets with the
Molecular Biology of the Cell Ontology. It includes the functionalities applied to the KPMP datasets in the reference atlas paper and additional functionalities.

*/



using System;
using Enumerations;
using Analysis;


class RNASeq_analysis
{
    public static void Write_task_report_comment(string comment)
    {
        System.IO.StreamWriter report = new System.IO.StreamWriter(Global_directory_class.Get_report_fileName(), true);
        report.WriteLine(comment);
        report.Close();
    }

    public static void Write_task_report_function_start(string functionName)
    {
        System.IO.StreamWriter report = new System.IO.StreamWriter(Global_directory_class.Get_report_fileName(), true);
        report.WriteLine("----------------------------------------------------------------------");
        report.WriteLine("Begin - " + functionName);
        report.WriteLine("");
        report.Close();
    }

    public static void Write_task_report_function_end(string functionName)
    {
        System.IO.StreamWriter report = new System.IO.StreamWriter(Global_directory_class.Get_report_fileName(), true);
        report.WriteLine("");
        report.WriteLine("End - " + functionName);
        report.WriteLine("----------------------------------------------------------------------");
        report.WriteLine();
        report.WriteLine();
        report.Close();
    }

    static void Main_kidney()
    {
        if (System.IO.File.Exists(Global_directory_class.Get_report_fileName()))
        { System.IO.File.Delete(Global_directory_class.Get_report_fileName()); }

        string report_function;

        #region Generate marker genes and proteins
        report_function = "Generate and write marker genes and proteins for enrichment and humanBase analysis";
        Write_task_report_function_start(report_function);
        Write_task_report_comment("This function depends on the experimental gene and protein lists in the folder experimental data");
        KPMP_analysis_reference_atlas_manuscript_class.Generate_and_write_marker_genes_and_proteins_for_enrichment_and_humanBase_analysis();
        Write_task_report_function_end(report_function);
        #endregion

        #region Enrichment analysis for the generation of pathway networks
        report_function = "Generate scp networks using Molecular Biology of the Cell Ontology and dynamic enrichment analysis";
        Write_task_report_function_start(report_function);
        Write_task_report_comment("This function depends on the marker genes and proteins generated above.");
        Write_task_report_comment("It also depends on the files \"MBCO_v1.1_gene-SCP_associations_human.txt\"");
        Write_task_report_comment("and \"MBCO_v1.1_inferred_SCP_relationships.txt\" that can be downloaded");
        Write_task_report_comment("from www.mbc-ontology.org or github.com/SBCNY/Molecular-Biology-of-the-Cell.");
        Write_task_report_comment("You can find both files on the website in the folder \"\\MBCO_windows_application\\MBCO_datsets\\\"");
        Write_task_report_comment("They need to be copied into the folder \"\\KPMP_reference_atlas_code\\MBCO_datasets\\\" on your hard drive");
        KPMP_analysis_reference_atlas_manuscript_class.Generate_scp_networks_using_MBCO_and_dynamic_enrichment_analysis();
        Write_task_report_comment("Generated networks can be visualized with a graph editor such as yED.");
        Write_task_report_comment("");
        Write_task_report_comment("The windows application that can be downloaded from this website allows analysis of custom data sets");
        Write_task_report_comment("with the Molecular Biology of the Cell Ontology and offers more functionalities than this script.");
        Write_task_report_function_end(report_function);
        #endregion

        #region Post hoc power analysis for single cell RNAseq
        report_function = "Process PostHoc Power analysis results for single cell RNAseq";
        Write_task_report_function_start(report_function);
        Write_task_report_comment("This function depends on the output of the R-script");
        Write_task_report_comment("\"RunBeforeCSharpScript_PostHocPower_singleNucleusCell_analysis.R\" \"Premiere\" 1");
        KPMP_analysis_reference_atlas_manuscript_class.Analyse_leavePatientsOutAnalysis("Premiere_AAAS_2022February14");
        Write_task_report_comment("The output of this function can be visualized with the R-script");
        Write_task_report_comment("\"RunAfterCSharpScript_PostHocPower_visualizeResults.R\" \"Premiere\"");
        Write_task_report_function_end(report_function);
        #endregion

        #region Post hoc power analysis for single nuclues RNAseq
        report_function = "Process PostHoc Power analysis results for single nucleus RNAseq";
        Write_task_report_function_start(report_function);
        Write_task_report_comment("This function depends on the output of the R-script");
        Write_task_report_comment("\"RunBeforeCSharpScript_PostHocPower_singleNucleusCell_analysis.R\" \"UCSD\" 1");
        KPMP_analysis_reference_atlas_manuscript_class.Analyse_leavePatientsOutAnalysis("UCSD_AAAS_2022February14");
        Write_task_report_comment("The output of this function can be visualized with the R-script");
        Write_task_report_comment("\"RunAfterCSharpScript_PostHocPower_visualizeResults.R\" \"UCSD\"");
        Write_task_report_function_end(report_function);
        #endregion

        #region Predicted energy metabolic pathways
        report_function = "Predict metabolic pathway activities involved in energy generation in a rules-based manner";
        Write_task_report_function_start(report_function);
        Write_task_report_comment("This function depends on the marker genes and proteins generated by a previous function");
        KPMP_analysis_reference_atlas_manuscript_class.Predict_metabolic_pathway_activities_in_a_rules_based_manner();
        Write_task_report_comment("The output of this function can be visualized with the R-scripts");
        Write_task_report_comment("\"RunAfterCSharpScript_metabolicEnergyPathways_averaged.R\"");
        Write_task_report_comment("and \"RunAfterCSharpScript_metabolicEnergyPathways_all.R\"");
        Write_task_report_function_end(report_function);
        #endregion

        #region Cross omic comparision - Figure 2C/D/E
        report_function = "Compare participant specific data across omic technologies";
        Write_task_report_function_start(report_function);
        Write_task_report_comment("This function is very time consumptive and depends on the output of the R-script");
        Write_task_report_comment("\"RunBeforeCSharpScript_calculate_meanExpressionValues_in_integrated_clusterAnalysis.R\"");
        KPMP_analysis_reference_atlas_manuscript_class.Compare_participant_specific_data_across_omic_technologies();
        Write_task_report_comment("Its results can be visualized with the R-scripts");
        Write_task_report_comment("\"RunAfterCSharpScript_crossOmicsComparison_hierarchicalClustering.R\"");
        Write_task_report_comment("and \"RunAfterCSharpScript_crossOmicsComparison_correlationCoefficients.R\"");
        Write_task_report_function_end(report_function);
        #endregion
    }

    static void Main()
    {
        switch (Global_class.Input_data)
        {
            case Input_data_enum.Kidney:
                Main_kidney();
                break;
            default:
                throw new Exception("input data is not considered");
        }
    }
}
