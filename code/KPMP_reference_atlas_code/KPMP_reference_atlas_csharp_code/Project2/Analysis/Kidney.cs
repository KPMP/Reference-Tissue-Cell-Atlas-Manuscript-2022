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

*/


using System;
using System.Drawing;
using System.Collections.Generic;
using Highthroughput_data;
using Enumerations;
using ReadWrite;
using Enrichment_2018;
using KPMP;
using Mbco_enrichment;
using MBCO_data;
using Mbco_ontology_enrichment;
using MBCO_network_and_nodes;
using MBCO_input_for_windows_form_application;

namespace Analysis
{
    class KPMP_analysis_reference_atlas_manuscript_class
    {
        public static void Generate_and_write_marker_genes_and_proteins_for_enrichment_and_humanBase_analysis()
        {
            string directory;
            string dataset;
            string results_directory = Global_directory_class.Marker_genes_proteins_directory;
            float max_pvalue_or_adj_pvalue = 0.05F;
            int keep_top_DEGsDEPs_for_mbco_enrichment = 300;
            int keep_top_DEGsDEPs_for_combined_dataset = 2000;
            KPMP_value_type_enum bulk_value_type_1st = KPMP_value_type_enum.Minus_log10_pvalue;
            KPMP_value_type_enum single_value_type_1st = KPMP_value_type_enum.Minus_log10_pvalue_adjusted;
            KPMP_value_type_enum bulk_value_type_2nd = KPMP_value_type_enum.Log2_ratioavg;
            KPMP_value_type_enum single_value_type_2nd = KPMP_value_type_enum.Log_ratioavg;
            string[] keep_patientIDs = new string[] { "allPatients", "AllPatients" };// Get_keep_patientIDs();

            KPMP_standardized_dataset_class combined_standardized_data = new KPMP_standardized_dataset_class();
            KPMP_standardized_dataset_class current_standardized_data;

            KPMP_singleRNASeqCluster_class single_cluster;
            KPMP_integration_paper_metadata_class add_dataset_patient;
            KPMP_integration_paper_metadata_class combined_dataset_patient = new KPMP_integration_paper_metadata_class();

            #region LMD RNASeq
            dataset = KPMP_dataset_name_class.Lmd_rnaseq_iuosu;
            string lmd_rnaSeq_parameter_string = dataset.ToString();
            KPMP_subsegmental_rnaSeq_class subsegmental_RNASeq = new KPMP_subsegmental_rnaSeq_class();
            subsegmental_RNASeq.Generate_by_reading_r_output("Integration_paper_SubSegmentalRNASeq_OSUIU_2019-05-21//TTest//", "Sub-segmenal_RNASeq_analysis_ttest.txt", dataset);
            add_dataset_patient = subsegmental_RNASeq.Get_deep_copy_of_dataset_patient_instance();
            combined_dataset_patient.Add_other(add_dataset_patient);
            subsegmental_RNASeq.Keep_only_lines_with_indicated_patientIDs(keep_patientIDs);
            subsegmental_RNASeq.Keep_only_indicated_valueTypes(new KPMP_value_type_enum[] { KPMP_value_type_enum.Log2_ratioavg, KPMP_value_type_enum.Minus_log10_pvalue });
            subsegmental_RNASeq.Keep_only_lines_with_minimum_selectedValueType_and_all_other_values_of_the_same_lines(-(float)Math.Log10(max_pvalue_or_adj_pvalue), KPMP_value_type_enum.Minus_log10_pvalue);
            current_standardized_data = subsegmental_RNASeq.Generate_standardized_dataset_instance(bulk_value_type_1st, bulk_value_type_2nd);
            current_standardized_data.Check_if_all_values_are_positive();
            combined_standardized_data.Add_deep_copy_of_other_after_checking_if_both_have_not_been_filtered_by_cutoff(current_standardized_data);

            current_standardized_data = new KPMP_standardized_dataset_class();
            #endregion

            #region LMD proteomics
            KPMP_subsegmental_proteomics_class subsegmental_prot = new KPMP_subsegmental_proteomics_class();
            dataset = KPMP_dataset_name_class.Lmd_proteomics_iuosu;
            string LMD_proteomics_parameter_string = dataset.ToString();
            subsegmental_prot.Generate_for_regular_analysis("Integration_paper_SubSegmentalProteomics_OSUIU_2019-06-14", dataset);
            add_dataset_patient = subsegmental_prot.Get_deep_copy_of_dataset_patient_instance();
            combined_dataset_patient.Add_other(add_dataset_patient);
            subsegmental_prot.Keep_only_lines_with_indicated_patientIDs(keep_patientIDs);
            subsegmental_prot.Keep_only_lines_with_minimum_selectedValueType_and_all_other_values_of_the_same_lines(-(float)Math.Log10(max_pvalue_or_adj_pvalue), KPMP_value_type_enum.Minus_log10_pvalue);
            current_standardized_data = subsegmental_prot.Generate_kpmp_standardized_dataset_instance(bulk_value_type_1st, bulk_value_type_2nd);
            current_standardized_data.Check_if_all_values_are_positive();
            combined_standardized_data.Add_deep_copy_of_other_after_checking_if_both_have_not_been_filtered_by_cutoff(current_standardized_data);
            current_standardized_data = new KPMP_standardized_dataset_class();
            #endregion

            #region Spatial metabolomics
            KPMP_metabolite_class metabolites = new KPMP_metabolite_class();
            metabolites.Generate(new string[] { "18-139", "18-142", "18-342" });
            DE_class de_metabolites = metabolites.Generate_de_instance_with_metabolite_hmdbIDs_filled_with_highest_rank_for_associated_molecular_formula("");
            de_metabolites.Write_file(results_directory, KPMP_dataset_name_class.Spatial_metabolomics + ".txt");
            #endregion

            #region Single nucleus UCSD
            single_cluster = new KPMP_singleRNASeqCluster_class();
            directory = Global_directory_class.Experimental_data_directory + "Integration_paper_SingleNucleusRNASeq_UCSD_DEGs_2019-06-10\\";
            dataset = KPMP_dataset_name_class.SingleNucleus_ucsd;
            single_cluster.Generate_and_override_current_array(directory, dataset);
            single_cluster.Remove_selected_cluster(new string[] { "PT-3", "PC-2", "Unk" });
            add_dataset_patient = single_cluster.Get_deep_copy_of_dataset_patient_instance();
            combined_dataset_patient.Add_other(add_dataset_patient);
            single_cluster.Keep_only_lines_with_indicated_patientIDs(keep_patientIDs);
            single_cluster.Keep_only_lines_with_adjusted_pvalue_below_alpha(max_pvalue_or_adj_pvalue);
            current_standardized_data = single_cluster.Generate_standardized_dataset_on_expression_values(single_value_type_1st, single_value_type_2nd);
            combined_standardized_data.Add_deep_copy_of_other_after_checking_if_both_have_not_been_filtered_by_cutoff(current_standardized_data);
            current_standardized_data = new KPMP_standardized_dataset_class();
            #endregion

            #region Near single cell proteomics
            KPMP_nearSingleCell_proteomics_UCSF_class nearSingleCell_prot = new KPMP_nearSingleCell_proteomics_UCSF_class();
            dataset = KPMP_dataset_name_class.NearSingleCell_proteomics_ucsf;
            string NSC_proteomics_parameter_string = dataset.ToString();
            nearSingleCell_prot.Generate_including_reversal_of_fold_changes_for_pt("Integration_paper_NearSingleCellProteomics_UCSF_2019-06-14", dataset);
            add_dataset_patient = nearSingleCell_prot.Get_deep_copy_of_dataset_patient_instance();
            combined_dataset_patient.Add_other(add_dataset_patient);
            nearSingleCell_prot.Keep_only_indicated_patientIDs(keep_patientIDs);
            nearSingleCell_prot.Keep_only_lines_with_minimum_selectedValueType_and_all_other_values_of_the_same_lines(-(float)Math.Log10(max_pvalue_or_adj_pvalue), KPMP_value_type_enum.Minus_log10_pvalue);
            current_standardized_data = nearSingleCell_prot.Generate_kpmp_standardized_dataset_instance(bulk_value_type_1st, bulk_value_type_2nd);
            current_standardized_data.Check_if_all_values_positive_or_zero();
            combined_standardized_data.Add_deep_copy_of_other_after_checking_if_both_have_not_been_filtered_by_cutoff(current_standardized_data);
            current_standardized_data = new KPMP_standardized_dataset_class();
            #endregion

            #region Single cell Premiere
            single_cluster = new KPMP_singleRNASeqCluster_class();
            directory = Global_directory_class.Experimental_data_directory + "Integration_paper_SingleCellRNASeq_PREMIERE_2019-06\\";
            dataset = KPMP_dataset_name_class.SingleCell_premiere;
            single_cluster.Generate_and_override_current_array(directory, dataset);
            add_dataset_patient = single_cluster.Get_deep_copy_of_dataset_patient_instance();
            combined_dataset_patient.Add_other(add_dataset_patient);
            single_cluster.Keep_only_lines_with_indicated_patientIDs(keep_patientIDs);
            single_cluster.Keep_only_lines_with_adjusted_pvalue_below_alpha(max_pvalue_or_adj_pvalue);
            current_standardized_data = single_cluster.Generate_standardized_dataset_on_expression_values(single_value_type_1st, single_value_type_2nd);
            current_standardized_data.Check_if_all_values_positive_or_zero();
            combined_standardized_data.Add_deep_copy_of_other_after_checking_if_both_have_not_been_filtered_by_cutoff(current_standardized_data);
            current_standardized_data = new KPMP_standardized_dataset_class();
            #endregion

            combined_standardized_data.Keep_only_lines_with_1stValue_aboveEqual_cutoff_and_check_if_1st_value_is_among_input_values(-Math.Log10(max_pvalue_or_adj_pvalue), KPMP_value_type_enum.Minus_log10_pvalue, KPMP_value_type_enum.Minus_log10_pvalue_adjusted);
            combined_standardized_data.Keep_genes_with_indicated_valueTypes1st2nd_if_below_maxFractionalRank_and_check_if_at_least_one_kept_value1stor2nd_is_not_infinity_keep_all_other_genes(keep_top_DEGsDEPs_for_combined_dataset, KPMP_value_type_enum.Minus_log10_pvalue, KPMP_value_type_enum.Log2_ratioavg);
            combined_standardized_data.Keep_genes_with_indicated_valueTypes1st2nd_if_below_maxFractionalRank_and_check_if_at_least_one_kept_value1stor2nd_is_not_infinity_keep_all_other_genes(keep_top_DEGsDEPs_for_combined_dataset, KPMP_value_type_enum.Minus_log10_pvalue_adjusted, KPMP_value_type_enum.Log_ratioavg);
            combined_standardized_data.Check_if_recalcuated_fractional_ranks_based_on_descending_values1st_and_values2nd_agree_with_existing();
            combined_standardized_data.Check_if_all_fractional_ranks_are_below_cutoff(keep_top_DEGsDEPs_for_combined_dataset);
            combined_standardized_data.Check_if_all_values_positive_or_zero();
            combined_standardized_data.Check_if_no_values_is_NaN();
            combined_standardized_data.Write(results_directory, "KPMP_combined_data_maxP_maxAdjP" + max_pvalue_or_adj_pvalue + "_top" + keep_top_DEGsDEPs_for_combined_dataset + "MarkerGenesProteins.txt", true);
            combined_standardized_data.Keep_genes_with_indicated_valueTypes1st2nd_if_below_maxFractionalRank_and_check_if_at_least_one_kept_value1stor2nd_is_not_infinity_keep_all_other_genes(keep_top_DEGsDEPs_for_mbco_enrichment, KPMP_value_type_enum.Minus_log10_pvalue, KPMP_value_type_enum.Log2_ratioavg);
            combined_standardized_data.Keep_genes_with_indicated_valueTypes1st2nd_if_below_maxFractionalRank_and_check_if_at_least_one_kept_value1stor2nd_is_not_infinity_keep_all_other_genes(keep_top_DEGsDEPs_for_mbco_enrichment, KPMP_value_type_enum.Minus_log10_pvalue_adjusted, KPMP_value_type_enum.Log_ratioavg);

            combined_standardized_data.Write_one_file_foreach_dataset(Global_directory_class.Marker_genes_proteins_directory + "MBCO_enrichment_input\\", max_pvalue_or_adj_pvalue + "_top" + keep_top_DEGsDEPs_for_mbco_enrichment + "MarkerGenesProteins.txt", true);

            Dictionary<string, Color> datasetName_color_dict = new Dictionary<string, Color>();
            datasetName_color_dict.Add(KPMP_dataset_name_class.Lmd_rnaseq_iuosu, Color.Pink);
            datasetName_color_dict.Add(KPMP_dataset_name_class.Lmd_proteomics_iuosu, Color.LightSkyBlue);
            datasetName_color_dict.Add(KPMP_dataset_name_class.NearSingleCell_proteomics_ucsf, Color.Blue);
            datasetName_color_dict.Add(KPMP_dataset_name_class.SingleNucleus_ucsd, Color.DarkRed);
            datasetName_color_dict.Add(KPMP_dataset_name_class.SingleCell_premiere, Color.Red);

            MBCO_windows_form_input_class mbco_data = new MBCO_windows_form_input_class();
            mbco_data.Generate_from_KPMP_standardized_dataset_instance_and_add_to_array(combined_standardized_data);
            mbco_data.Add_colors(datasetName_color_dict);
            mbco_data.Write_one_file_for_each_center(results_directory + "MBCO_input_for_windows_application\\");
        }

        #region Predict metabolic activities for nephron segments
        private static Enrichment2018_results_class Analyse_metabolic_pathways_of_combined_standardized_data_and_return_enrichment_results_disease(out string parameter_string)
        {
            string read_directory = Global_directory_class.Marker_genes_proteins_directory;
            string enrichment_results_subdirectory = Global_directory_class.MBCO_enrichment_results_directory;

            KPMP_standardized_dataset_class combined_standardized_data = new KPMP_standardized_dataset_class();
            float max_pvalue = 0.05F;
            combined_standardized_data.Read(read_directory, "KPMP_combined_data_maxP_maxAdjP" + max_pvalue + "_top2000MarkerGenesProteins.txt", true);
            combined_standardized_data.Check_if_all_pvalues_adjusted_pvalues_are_below_cutoff(max_pvalue);
            combined_standardized_data.Check_if_all_values_positive_or_zero();

            KPMP_value_type_enum value_type_1st_single = KPMP_value_type_enum.Minus_log10_pvalue_adjusted;
            KPMP_value_type_enum value_type_2nd_single = KPMP_value_type_enum.Log_ratioavg;
            KPMP_value_type_enum value_type_1st_bulk = KPMP_value_type_enum.Minus_log10_pvalue;
            KPMP_value_type_enum value_type_2nd_bulk = KPMP_value_type_enum.Log2_ratioavg;

            int top_x_genesProteins = 500;
            parameter_string = "maxP_maxAdjP" + max_pvalue + "_top" + top_x_genesProteins + "DEGsDEPs";

            combined_standardized_data.Keep_genes_with_indicated_valueTypes1st2nd_if_below_maxFractionalRank_and_check_if_at_least_one_kept_value1stor2nd_is_not_infinity_keep_all_other_genes(top_x_genesProteins, value_type_1st_single, value_type_2nd_single);
            combined_standardized_data.Keep_genes_with_indicated_valueTypes1st2nd_if_below_maxFractionalRank_and_check_if_at_least_one_kept_value1stor2nd_is_not_infinity_keep_all_other_genes(top_x_genesProteins, value_type_1st_bulk, value_type_2nd_bulk);
            combined_standardized_data.Check_if_all_fractional_ranks_are_below_cutoff(top_x_genesProteins);

            Dictionary<string, Dictionary<string, string>> center_cellSubType_setIntegrationTerm_dict = Get_center_cellSubType_setIntegrationTerm_dict_for_medullary_cell_types();
            combined_standardized_data.Reset_integrationTerm_for_indicated_center_cellTypes(center_cellSubType_setIntegrationTerm_dict);


            string[] conditions = combined_standardized_data.Get_all_unique_ordered_datasets();

            Downstream_analysis_2020_class downstream = new Downstream_analysis_2020_class();
            downstream.Options.Data_value_signs_of_interest = new Data_value_signs_of_interest_enum[] { Data_value_signs_of_interest_enum.Combined };
            downstream.Options.Ontologies = new Ontology_type_enum[] { Ontology_type_enum.Metabolic_energy_generation_pathways };
            downstream.Options.Top_top_x_predictions = 9999;
            downstream.Options.Only_keep_filtered_results = true;
            downstream.Options.Add_missing_scps_identified_in_other_conditions = true;
            downstream.Options.Write_results = false;

            KPMP_standardized_dataset_class combined_standardized_data_current_dataset;

            string[] datasets = combined_standardized_data.Get_all_unique_ordered_datasets();
            string dataset;
            int datasets_length = datasets.Length;
            string[] bg_genes;
            Enrichment2018_results_class combined_enrichment_results = new Enrichment2018_results_class();
            for (int indexD = 0; indexD < datasets_length; indexD++)
            {
                dataset = datasets[indexD];
                combined_standardized_data_current_dataset = combined_standardized_data.Deep_copy();
                combined_standardized_data_current_dataset.Keep_only_indicated_datasets(dataset);
                DE_class de = combined_standardized_data_current_dataset.Generate_de_instance_alternatively_and_fill_with_one_of_indicated_valueTypes_and_check_if_duplicated(value_type_2nd_bulk, value_type_2nd_single);
                bg_genes = combined_standardized_data_current_dataset.Get_bgGenes_of_indicated_dataset(dataset);
                downstream.Generate(bg_genes);
                Enrichment2018_results_class add_enrichment_results = downstream.Analyse_de_instance_and_return_unfiltered_enrichment_results(de, enrichment_results_subdirectory, "")[0];
                combined_enrichment_results.Add_other(add_enrichment_results);
            }
            return combined_enrichment_results;
        }

        private static Dictionary<string, Dictionary<string, string>> Get_center_cellSubType_setIntegrationTerm_dict_for_medullary_cell_types()
        {
            Dictionary<string, Dictionary<string, string>> center_cellSubType_setIntegrationTerm_dict = new Dictionary<string, Dictionary<string, string>>();
            Dictionary<string, string> ucsd_cellSubType_setIntegrationTerm = new Dictionary<string, string>();
            ucsd_cellSubType_setIntegrationTerm.Add("DTL", "Descending_limb_medulla-allPatients");
            ucsd_cellSubType_setIntegrationTerm.Add("ATL-1", "Thin_ascending_limb_medulla-allPatients");
            ucsd_cellSubType_setIntegrationTerm.Add("ATL-2", "Thin_ascending_limb_medulla-allPatients");
            ucsd_cellSubType_setIntegrationTerm.Add("ATL-3", "Thin_ascending_limb_medulla-allPatients");
            ucsd_cellSubType_setIntegrationTerm.Add("TAL-1", "Thick_ascending_limb_medulla-allPatients");
            ucsd_cellSubType_setIntegrationTerm.Add("PC-3", "Principal_cell_medulla-allPatients");
            ucsd_cellSubType_setIntegrationTerm.Add("IC-A2", "Intercalated_cell_medulla-allPatients");
            center_cellSubType_setIntegrationTerm_dict.Add(KPMP_dataset_name_class.SingleNucleus_ucsd, ucsd_cellSubType_setIntegrationTerm);
            return center_cellSubType_setIntegrationTerm_dict;
        }

        public static void Predict_metabolic_pathway_activities_in_a_rules_based_manner()
        {
            string directory = Global_directory_class.Metabolic_activities_results_directory;
            string parameter_string;
            Dictionary<string, Dictionary<string, string>> center_cellSubType_setIntegrationTerm_dict = Get_center_cellSubType_setIntegrationTerm_dict_for_medullary_cell_types();
            Enrichment2018_results_class combined_enrichment_results = Analyse_metabolic_pathways_of_combined_standardized_data_and_return_enrichment_results_disease(out parameter_string);

            Enrichment2018_results_class combined_enrichment_results_write = combined_enrichment_results.Deep_copy();
            combined_enrichment_results_write.Extract_integrationGroup_from_sampleName_and_modify_sampleName();
            combined_enrichment_results_write.Collapse_bulk_transcriptomics_proteomics_on_segments();
            //combined_enrichment_results_write.Keep_only_lines_with_non_zero_minusLog10Pvalues();
            combined_enrichment_results_write.Write(directory, "Energy_scps_all_cellSubTypes_segments_datasets.txt");

            KPMP_singleCell_annotation_class premiere_annotation = new KPMP_singleCell_annotation_class();
            premiere_annotation.Generate_by_reading("lcm_cell_assignment_barcodes_premiere_2021April22.txt");
            KPMP_singleCell_annotation_class ucsd_annotation = new KPMP_singleCell_annotation_class();
            ucsd_annotation.Generate_by_reading("lcm_cell_assignment_barcodes_ucsd_2021April22.txt");

            KPMP_singleCell_subtype_count_class subtype_counts = new KPMP_singleCell_subtype_count_class();
            subtype_counts.Generate_from_annotations_add_to_array(premiere_annotation, KPMP_dataset_name_class.SingleCell_premiere, "-allPatients");
            subtype_counts.Generate_from_annotations_add_to_array(ucsd_annotation, KPMP_dataset_name_class.SingleNucleus_ucsd, "-allPatients");
            subtype_counts.Reset_kpmpIntegration_term_as_indicated_in_dictionary(center_cellSubType_setIntegrationTerm_dict);
            subtype_counts.Write(directory, "Cell_subtype_counts.txt");

            Dictionary<string, string[]> first_addDefaultFirstSibling_if_no_siblings_dict = new Dictionary<string, string[]>();
            first_addDefaultFirstSibling_if_no_siblings_dict.Add("Glycolysis and gluconeogenesis shared enzymes", new string[] { "Glycolysis specific enzymes", "Gluconeogenesis specific enzymes" });
            first_addDefaultFirstSibling_if_no_siblings_dict.Add("Ketone body formation and catabolism shared enzymes", new string[] { "Ketone body catabolism specific enzymes", "HMG-CoA pathway of ketone body formation specific enzymes" });
            first_addDefaultFirstSibling_if_no_siblings_dict.Add("Mitochondrial and peroxisomal beta oxidation shared enzymes", new string[] { "Mitochondrial beta oxidation", "Peroxisomal beta oxidation" });

            Dictionary<string, string[]> second_addDefaultFirstSibling_if_no_siblings_dict = new Dictionary<string, string[]>();
            second_addDefaultFirstSibling_if_no_siblings_dict.Add("Glycolysis specific enzymes", new string[] { "Pyruvate dehydrogenase", "Lactate dehydrogenase" });
            second_addDefaultFirstSibling_if_no_siblings_dict.Add("Lactate dehydrogenase", new string[] { "Glycolysis specific enzymes", "Gluconeogenesis specific enzymes" });
            second_addDefaultFirstSibling_if_no_siblings_dict.Add("Pyruvate dehydrogenase", new string[] { "Glycolysis specific enzymes" });

            Dictionary<string, string[]> parentScp_childScpsRequirement_dict = new Dictionary<string, string[]>();
            parentScp_childScpsRequirement_dict.Add("Anaerobic glycolysis", new string[] { "Glycolysis specific enzymes", "Lactate dehydrogenase" });
            parentScp_childScpsRequirement_dict.Add("Aerobic glycolysis", new string[] { "Glycolysis specific enzymes", "Pyruvate dehydrogenase" });
            parentScp_childScpsRequirement_dict.Add("Peroxisomal beta oxidation", new string[] { "Peroxisomal beta oxidation specific enzymes" });
            parentScp_childScpsRequirement_dict.Add("Mitochondrial beta oxidation", new string[] { "Mitochondrial beta oxidation specific enzymes" });
            parentScp_childScpsRequirement_dict.Add("Glycolysis", new string[] { "Glycolysis specific enzymes" });
            parentScp_childScpsRequirement_dict.Add("Ketone body catabolism", new string[] { "Ketone body catabolism specific enzymes" });

            Dictionary<string, string[]> parentScp_childScpExclustion_dict = new Dictionary<string, string[]>();
            //parentScp_childScpExclustion_dict.Add("Glycolysis", new string[] { "Lactate dehydrogenase", "Pyruvate dehydrogenase" });
            parentScp_childScpExclustion_dict.Add("Glycolysis", new string[] { "Anaerobic glycolysis", "Aerobic glycolysis" });
            parentScp_childScpExclustion_dict.Add("Glycolysis and gluconeogenesis", new string[] { "Glycolysis specific enzymes", "Gluconeogenesis specific enzymes" });
            parentScp_childScpExclustion_dict.Add("Ketone body metabolism", new string[] { "Ketone body catabolism specific enzymes", "HMG-CoA pathway of ketone body formation specific enzymes" });

            KPMP_singleCell_averageMinusLog10Pvalue_class averageEnrichment_selected_scps = new KPMP_singleCell_averageMinusLog10Pvalue_class();

            Enrichment2018_results_class enrichment_results_selected_scps = combined_enrichment_results.Deep_copy();
            enrichment_results_selected_scps.Add_first_sibling_as_default_sibling_if_all_siblings_are_missing(first_addDefaultFirstSibling_if_no_siblings_dict);
            enrichment_results_selected_scps.Add_first_sibling_as_default_sibling_if_all_siblings_are_missing(second_addDefaultFirstSibling_if_no_siblings_dict);
            enrichment_results_selected_scps.Keep_scp_only_if_all_indicated_requiredScps_have_non_zero_minusLog10Pvalue(parentScp_childScpsRequirement_dict);
            enrichment_results_selected_scps.Remove_scp_only_if_at_least_one_indicated_exclusionScp_has_non_zero_minusLog10value(parentScp_childScpExclustion_dict);
            averageEnrichment_selected_scps.Generate_from_enrichment_results_and_cell_subtype_counts_and_add_to_array(enrichment_results_selected_scps, subtype_counts);

            averageEnrichment_selected_scps.Keep_only_indicated_centers(new string[] { KPMP_dataset_name_class.SingleCell_premiere, KPMP_dataset_name_class.SingleNucleus_ucsd, KPMP_dataset_name_class.Lmd_rnaseq_iuosu });
            averageEnrichment_selected_scps.Generate_new_enrichment_results_lines_by_averaging_average_values_over_all_assays_of_same_integration_term_and_add_to_array("RNAseq");
            averageEnrichment_selected_scps.Keep_only_indicated_centers("RNAseq");
            averageEnrichment_selected_scps.Write(directory, "Energy_scps_rulesBasedAveraged_" + parameter_string + ".txt");
        }
        #endregion

        private static Ontology_enrichment_class Set_kpmp_integrationGroup_and_setWithinIntegrationGroup(Ontology_enrichment_class onto_enrich)
        {
            Ontology_enrichment_line_class onto_enrich_line;
            List<Ontology_enrichment_line_class> keep = new List<Ontology_enrichment_line_class>();
            int onto_length = onto_enrich.Enrich.Length;
            string[] splitStrings;
            int splitStrings_length;
            List<string> final_splitStrings = new List<string>();
            for (int indexO = 0; indexO < onto_length; indexO++)
            {
                onto_enrich_line = onto_enrich.Enrich[indexO];
                splitStrings = onto_enrich_line.Sample_name.Split('$');
                splitStrings_length = splitStrings.Length;
                onto_enrich_line.Sample_name = splitStrings[1] + ": " + splitStrings[2];
                onto_enrich_line.Integration_group = (string)splitStrings[0].Replace("-allPatients","");
                onto_enrich_line.Set_within_integration_group = (string)splitStrings[1].Clone();
                if (!onto_enrich_line.Integration_group.Equals("no integration term"))
                {
                    keep.Add(onto_enrich_line);
                }
            }
            onto_enrich.Enrich = keep.ToArray();
            return onto_enrich;
        }

        public static void Generate_scp_networks_using_MBCO_and_dynamic_enrichment_analysis()
        {
            Mbc_enrichment_pipeline_class mbco_enrichment_pipeline = new Mbc_enrichment_pipeline_class( Ontology_type_enum.Mbco_level3 );
            mbco_enrichment_pipeline.Options.Maximum_pvalue_for_standardDynamicEnrichment = 1F;
            mbco_enrichment_pipeline.Options.Scp_levels_for_dynamicEnrichment = new int[] { 3 };
            mbco_enrichment_pipeline.Options.Kept_top_predictions_dynamicEnrichment_per_level = new int[] { -1,   //Position 0 will be ignored
                                                                                                            -1,   //Position 1 will be ignored
                                                                                                            -1,   //Position 2 will be ignored
                                                                                                             7,   //Keep top x level 3 SCPs/SCP-unions
                                                                                                            -1  };//Position 4 will be ignored
            mbco_enrichment_pipeline.Options.Numbers_of_merged_scps_for_dynamicEnrichment_per_level = new int[5][];

            #region KPMP reference atlas datasets
            string[] custom_dataset_names;
            string custom_dataset_name;
            int top_degs = 300;
            int top_deps = 300;
            float pvalue = 0.05F;
            float singleCellNucleus_adjPvalue = 0.05F;
            string topDEGsDEPs_string = "_maxAdjPvalue" + singleCellNucleus_adjPvalue + "_maxPvalue" + pvalue + "_top" + top_degs + "Genes_top" + top_deps + "Proteins";
            custom_dataset_names = new string[] { KPMP_dataset_name_class.SingleNucleus_ucsd + "_adjPvalue" + singleCellNucleus_adjPvalue + "_top" + top_degs + "MarkerGenesProteins.txt",
                                                  KPMP_dataset_name_class.SingleCell_premiere + "_adjPvalue" + singleCellNucleus_adjPvalue + "_top" + top_degs + "MarkerGenesProteins.txt",
                                                  KPMP_dataset_name_class.Lmd_rnaseq_iuosu + "_pvalue" + pvalue + "_top" + top_degs + "MarkerGenesProteins.txt",
                                                  KPMP_dataset_name_class.Lmd_proteomics_iuosu + "_pvalue" + pvalue + "_top" + top_deps + "MarkerGenesProteins.txt",
                                                  KPMP_dataset_name_class.NearSingleCell_proteomics_ucsf + "_pvalue" + pvalue + "_top" + top_deps + "MarkerGenesProteins.txt"
                                                };
            #endregion

            int datasets_length = custom_dataset_names.Length;
            string[] keep_sample_names = new string[datasets_length];
            string[] bg_genes_names = new string[datasets_length];
            for (int indexD = 0; indexD < datasets_length; indexD++)
            {
                custom_dataset_name = custom_dataset_names[indexD];
                bg_genes_names[indexD] = System.IO.Path.GetFileNameWithoutExtension(custom_dataset_name) + "_bgGenes.txt";
                keep_sample_names[indexD] = custom_dataset_name.Split('_')[0];
            }

            string bg_genes_name;

            string[] bg_genes;
            Data_class data;

            Ontology_enrichment_class current_dynamic_enrichment;
            Ontology_enrichment_class combined_dynamic_enrichments = new Ontology_enrichment_class();

            mbco_enrichment_pipeline.Generate_for_all_runs();

            for (int indexD = 0; indexD < datasets_length; indexD++)
            {
                custom_dataset_name = custom_dataset_names[indexD];
                bg_genes_name = bg_genes_names[indexD];

                Custom_data_readWriteOptions_class custom_data_readWriteOptions = new Custom_data_readWriteOptions_class(custom_dataset_name);
                Custom_data_class custom_data = new Custom_data_class();
                custom_data.Generate_custom_data_instance(custom_data_readWriteOptions);
                data = custom_data.Generate_new_data_instance();
                if (data.Data.Length == 0) { throw new System.Exception(); }
                bg_genes = ReadWriteClass.Read_string_array_and_remove_non_letters_from_beginning_and_end_of_each_line(Global_directory_class.Marker_genes_proteins_for_MBCO_directory + bg_genes_name);
                mbco_enrichment_pipeline.Reset_mbco_and_adjust_bg_genes(bg_genes);
                current_dynamic_enrichment = mbco_enrichment_pipeline.Analyse_data_instance_via_dynamic_enrichment_analysis(data);
                combined_dynamic_enrichments.Add_other(current_dynamic_enrichment);
                mbco_enrichment_pipeline.Check_for_duplicated_enrichment_lines(current_dynamic_enrichment.Enrich);
            }

            combined_dynamic_enrichments = Set_kpmp_integrationGroup_and_setWithinIntegrationGroup(combined_dynamic_enrichments);
            combined_dynamic_enrichments.Write(Global_directory_class.MBCO_enrichment_results_directory, "Dynamic_scp_level3_enrichment_results.txt");
            System.Drawing.Color[] set_colors = new System.Drawing.Color[] { System.Drawing.Color.DarkRed, System.Drawing.Color.Red, System.Drawing.Color.Pink, System.Drawing.Color.LightSkyBlue, System.Drawing.Color.Blue };
            string[] set_abbreviations = new string[] { "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z" };

            Dictionary<string, Network_class> integrationGroup_network_dict = new Dictionary<string, Network_class>();// add_in.Generate_addProcess_network_dictionary();
            string[] podocyte_metabolite_pathways = new string[] { "Arachidonic acid metabolism", "Linolenic acid metabolism", "Sphingolipid metabolism", "GPI-anchor synthesis and attachment to protein", "Phosphoglyceride biosynthesis" };
            string[] proxTub_metabolite_pathways = new string[] { "Galactose metabolism", "Purine metabolism", "Phosphoglyceride biosynthesis", "Glycolysis and Gluconeogenesis", "Fructose and mannose metabolism", "D-Arginine and D-ornithine metabolism", "Carnitine shuttle", "Carnitine biosynthesis and transport" };
            string[] tal_metabolite_pathways = new string[] { "Purine metabolism", "Glycolysis and Gluconeogenesis", "Phosphoglyceride biosynthesis" };
            string[] ic_metabolite_pathways = new string[] { "Purine metabolism", "Phosphoglyceride biosynthesis" };
            string[] pc_metabolite_pathways = new string[] { "Purine metabolism", "Glycolysis and Gluconeogenesis", "Phosphoglyceride biosynthesis" };
            Network_class podocyte_metabolite_network = new Network_class();
            podocyte_metabolite_network.Add_single_nodes(podocyte_metabolite_pathways);
            Network_class proxTub_metabolite_network = new Network_class();
            proxTub_metabolite_network.Add_single_nodes(proxTub_metabolite_pathways);
            Network_class tal_metabolite_network = new Network_class();
            tal_metabolite_network.Add_single_nodes(tal_metabolite_pathways);
            Network_class pc_metabolite_network = new Network_class();
            pc_metabolite_network.Add_single_nodes(pc_metabolite_pathways);
            Network_class ic_metabolite_network = new Network_class();
            ic_metabolite_network.Add_single_nodes(ic_metabolite_pathways);

            integrationGroup_network_dict.Add("Podocyte_glomerulus", podocyte_metabolite_network);
            integrationGroup_network_dict.Add("Proximal_tubule", proxTub_metabolite_network);
            integrationGroup_network_dict.Add("Thick_ascending_limb", tal_metabolite_network);
            integrationGroup_network_dict.Add("Principal cells", pc_metabolite_network);
            integrationGroup_network_dict.Add("Intercalated cells", ic_metabolite_network);

            string keep_sample_name;
            int keep_length = keep_sample_names.Length;

            Dictionary<string, System.Drawing.Color> setAbbreviation_Color_dict = new Dictionary<string, System.Drawing.Color>();
            Dictionary<string, string> set_setAbbreviation_dict = new Dictionary<string, string>();
            for (int indexK = 0; indexK < keep_length; indexK++)
            {
                keep_sample_name = keep_sample_names[indexK];
                set_setAbbreviation_dict.Add(keep_sample_name, set_abbreviations[indexK]);
                setAbbreviation_Color_dict.Add(set_abbreviations[indexK], set_colors[indexK]);
            }
            setAbbreviation_Color_dict.Add("X", System.Drawing.Color.Lime);

            Mbc_network_based_integration_class mbco_network_integration = new Mbc_network_based_integration_class(mbco_enrichment_pipeline.Options.Ontology);
            mbco_network_integration.Options.Top_quantile_probability_of_scp_interactions_for_dynamic_enrichment = mbco_enrichment_pipeline.Options.Top_quantile_of_scp_interactions_for_dynamicEnrichment_per_level;
            mbco_network_integration.Options.Add_edges_that_connect_scps_between_sets = false;
            mbco_network_integration.Generate();

            mbco_network_integration.Options.Add_edges_that_connect_scps_between_sets = true;
            mbco_network_integration.Generate_and_write_integrative_network_for_dynamic_enrichment_results_of_each_integrationGroupName_only_defined_sets(combined_dynamic_enrichments, setAbbreviation_Color_dict, set_setAbbreviation_dict, "Dynamic_level3SCPs_", integrationGroup_network_dict);
        }

        public static void Analyse_leavePatientsOutAnalysis(string baseFileName)
        {
            string[] findMarkers_assay_slots; 
            string findMarkers_assay_slot;
            int[] maxPrincipalComponents = new int[] { };
            int topDEGs_for_SCPs = 300;
            int[] top_DEGs_for_cellType_identification_array = new int[] { 300 };
            int top_DEGs_for_cellType_identification;
            int top_DEGs_for_cellType_identification_array_length = top_DEGs_for_cellType_identification_array.Length;
            string parameter_string;
            string pc_string;
            string subdirectory = baseFileName + "\\";
            string final_results_subdirectory = "PostHocPower_" + baseFileName + "_finalResults\\";
            ReadWriteClass.Create_subdirectory_in_results_directory_if_it_does_not_exist(final_results_subdirectory);
            KPMP_singleNucleusCell_randomization_set_class randomization_sets = new KPMP_singleNucleusCell_randomization_set_class();
            switch (baseFileName)
            {
                case "Premiere_AAAS_2022February14":
                    randomization_sets.Generate_by_reading_for_singleCellDataset("Randomizations_SC RNASeq PREMIERE_2020June22.txt");
                    maxPrincipalComponents = new int[] { 30 };
                    findMarkers_assay_slots = new string[] { "SCT_data" };
                    break;
                case "UCSD_AAAS_2022February14":
                    randomization_sets.Generate_by_reading_for_singleCellDataset("Randomizations_SN RNASeq UCSDWU_2020June22.txt");
                    maxPrincipalComponents = new int[] { 30 };
                    findMarkers_assay_slots = new string[] { "SCT_data" };
                    break;
                default:
                    throw new Exception();
            }
            int findMarkers_assay_slots_length = findMarkers_assay_slots.Length;
            int max_principal_components_length = maxPrincipalComponents.Length;
            int max_principal_component;
            for (int indexPC = 0; indexPC < max_principal_components_length; indexPC++)
            {
                max_principal_component = maxPrincipalComponents[indexPC];
                for (int indexCellTypeDEGs = 0; indexCellTypeDEGs < top_DEGs_for_cellType_identification_array_length; indexCellTypeDEGs++)
                {
                    top_DEGs_for_cellType_identification = top_DEGs_for_cellType_identification_array[indexCellTypeDEGs];
                    for (int indexFindMarkers = 0; indexFindMarkers < findMarkers_assay_slots_length; indexFindMarkers++)
                    {
                        findMarkers_assay_slot = findMarkers_assay_slots[indexFindMarkers];

                        parameter_string = Get_parameter_string(max_principal_component, findMarkers_assay_slot, top_DEGs_for_cellType_identification, topDEGs_for_SCPs);
                        pc_string = "maxPC" + max_principal_component;

                        KPMP_barcode_cellType_class barcodes = new KPMP_barcode_cellType_class();
                        barcodes.Generate(subdirectory, pc_string);

                        KPMP_leavePatientsOut_cellType_assignment_class cellType_assignments = new KPMP_leavePatientsOut_cellType_assignment_class();
                        cellType_assignments.Options.Top_degs_for_cell_type_assignment = top_DEGs_for_cellType_identification;
                        cellType_assignments.Generate_from_barcode_cellType_instance(barcodes, randomization_sets);

                        KPMP_leavePatientsOut_analysis_class leavePatients_out_analysis = new KPMP_leavePatientsOut_analysis_class();
                        leavePatients_out_analysis.Generate_fisherExact_from_cell_type_assignments(cellType_assignments);
                        leavePatients_out_analysis.Write(final_results_subdirectory, baseFileName + "_" + parameter_string + ".txt");
                    }
                }
            }
        }

        private static string Get_parameter_string(int max_principal_component, string findMarkers_assay_slot, int keepTopDEGs_for_cellType, int keepTopDEGs_for_scps)
        {
            return "mPC" + max_principal_component + "_" + findMarkers_assay_slot + "_" + keepTopDEGs_for_cellType + "DEGsCT" + "_" + keepTopDEGs_for_scps + "DEGsSCP";
        }
 
        public static void Compare_participant_specific_data_across_omic_technologies()
        {
            KPMP_standardized_dataset_class standardized_dataset = new KPMP_standardized_dataset_class();
            KPMP_integration_paper_metadata_class combined_dataset_patients = new KPMP_integration_paper_metadata_class();
            KPMP_standardized_dataset_class add_standardized_dataset;

            string patient_samples_subdirectory = "IntegratedCluster_AvgExpression_RNAcounts\\";
            string results_directory = Global_directory_class.Results_directory +  "GlomVsProxTub\\";
            string baseFileName = "GlomVsProxTub";

            KPMP_singleNucleusCell_averageExpression_class singleNucleusCell_rnaSeq = new KPMP_singleNucleusCell_averageExpression_class();
            singleNucleusCell_rnaSeq.Generate(); 
            combined_dataset_patients.Add_other(singleNucleusCell_rnaSeq.Get_deep_copy_of_dataset_patient());
            singleNucleusCell_rnaSeq.Merge_proxTub_cluster_by_keeping_gene_with_highest_expression();
            singleNucleusCell_rnaSeq.Replace_all_cellType0_by_cellType1("PT", "ProxTub");

            add_standardized_dataset = singleNucleusCell_rnaSeq.Generate_standardized_dataset_with_all_valueTypes_filled_in_value_1st();
            standardized_dataset.Add_deep_copy_of_other_after_checking_if_both_have_not_been_filtered_by_cutoff(add_standardized_dataset);

            KPMP_subsegmentalRNASeq_rawCounts_class subsegRNASeq = new KPMP_subsegmentalRNASeq_rawCounts_class();
            subsegRNASeq.Generate();
            combined_dataset_patients.Add_other(subsegRNASeq.Get_deep_copy_of_dataset_patient_instance());
            add_standardized_dataset = subsegRNASeq.Generate_standardized_dataset_instance_with_all_values_filled_in_value_1st();
            standardized_dataset.Add_deep_copy_of_other_after_checking_if_both_have_not_been_filtered_by_cutoff(add_standardized_dataset);

            KPMP_subsegmental_proteomics_class subsegmental_proteomics = new KPMP_subsegmental_proteomics_class();
            subsegmental_proteomics.Generate_for_single_patient_analysis("Integration_paper_SubSegmentalProteomics_OSUIU_2019-06-14", KPMP_dataset_name_class.Lmd_proteomics_iuosu);
            combined_dataset_patients.Add_other(subsegmental_proteomics.Get_deep_copy_of_dataset_patient_instance());
            subsegmental_proteomics.Replace_all_segment0s_by_segment1s("TI", "ProxTub");
            add_standardized_dataset = subsegmental_proteomics.Generate_standardized_dataset_instance_with_all_values_filled_in_value_1st();
            standardized_dataset.Add_deep_copy_of_other_after_checking_if_both_have_not_been_filtered_by_cutoff(add_standardized_dataset);


            KPMP_nearSingleCell_proteomics_UCSF_class nearSingleCell_proteomics = new KPMP_nearSingleCell_proteomics_UCSF_class();
            nearSingleCell_proteomics.Generate_for_single_patient_analysis("Integration_paper_NearSingleCellProteomics_UCSF_2019-06-14", KPMP_dataset_name_class.NearSingleCell_proteomics_ucsf);
            nearSingleCell_proteomics.Replace_all_segment0s_by_segment1s("PT", "ProxTub");
            combined_dataset_patients.Add_other(nearSingleCell_proteomics.Get_deep_copy_of_dataset_patient_instance());
            add_standardized_dataset = nearSingleCell_proteomics.Generate_standardized_dataset_instance_with_all_values_filled_in_value_1st();
            standardized_dataset.Add_deep_copy_of_other_after_checking_if_both_have_not_been_filtered_by_cutoff(add_standardized_dataset);

            standardized_dataset.Add_integration_terms();
            standardized_dataset.Check_for_duplicates();
            standardized_dataset.Write(results_directory, "Standardized_dataset_collected.txt",true);

            combined_dataset_patients.Add_analysis_to_all_lines("Crossplatform comparison");
            combined_dataset_patients.Write_into_results_directory(patient_samples_subdirectory, "KPMP_dataset_patients.txt");

            string glom_integration_term = KPMP_data_integration_class.Get_podocyte_glomerulus_integration_term();
            string proxTub_integration_term = KPMP_data_integration_class.Get_proximal_tubule_integration_term();
            standardized_dataset.Keep_only_indicated_integration_types_and_datasetPatients_with_all_terms(new string[] { glom_integration_term, proxTub_integration_term });
            standardized_dataset.Add_missing_genes_to_each_dataset_condition_if_valueType1st_is_Average_or_singleValue_ifNot_throwException();
            standardized_dataset.Calculate_ratios_between_indicated_integration_terms_and_add_assuming_only_singleValues1st_for_term0_and_multiple_for_term1(glom_integration_term, proxTub_integration_term);
            standardized_dataset.Calculate_log2_ratio_values_1st_and_add_as_new_value_1st();
            standardized_dataset.Write(results_directory, "Standardized_dataset_for_glomVsProxTub_comparison_all.txt", false);

            string[] shared_genes = standardized_dataset.Get_all_genes_that_have_a_non_zero_value_of_indicated_value_type_1st_in_all_datasets(KPMP_value_type_enum.Single_value);

            KPMP_standardized_dataset_collapseOnMethod_class collapse_on_method;

            KPMP_standardized_dataset_class standardized_data_log2Ratios = standardized_dataset.Deep_copy();
            standardized_data_log2Ratios.Keep_only_indicated_valueTypes_1st(KPMP_value_type_enum.Log2_ratio_singlepatient);
            standardized_data_log2Ratios.Set_valueType_1st_as_score_of_interest();
            standardized_data_log2Ratios.Write(results_directory, "StandardizedData_log2Ratios.txt", true);

            collapse_on_method = new KPMP_standardized_dataset_collapseOnMethod_class();
            collapse_on_method.Generate_from_kpmp_standardized_dataset(standardized_data_log2Ratios);
            collapse_on_method.Write(results_directory, "Standardized_dataset_collapsed_on_method.txt");

            standardized_data_log2Ratios.Keep_only_indicated_gene_symbols(shared_genes);
            standardized_data_log2Ratios.Write(results_directory, "StandardizedData_log2Ratios_sharedGenesProteins.txt", true);
            collapse_on_method.Generate_from_kpmp_standardized_dataset(standardized_data_log2Ratios);
            collapse_on_method.Write(results_directory, "Standardized_dataset_collapsed_on_method_sharedGenesProteins.txt");

            KPMP_standardized_dataset_class standardized_data_log2Ratios_only_subsegmental = standardized_dataset.Deep_copy();
            standardized_data_log2Ratios_only_subsegmental.Keep_only_indicated_datasets(new string[] { "LMD Proteomics IU-OSU", "LMD RNASeq IU-OSU", "NSC proteomics UCSF" });
            string[] subsegmental_shared_genes = standardized_data_log2Ratios_only_subsegmental.Get_all_genes_that_have_a_non_zero_value_of_indicated_value_type_1st_in_all_datasets(KPMP_value_type_enum.Single_value);
            standardized_data_log2Ratios_only_subsegmental.Set_valueType_1st_as_score_of_interest();
            standardized_data_log2Ratios_only_subsegmental.Write(results_directory, "Standardized_data_log2Ratios_onlySubsegmental.txt", true);

            collapse_on_method = new KPMP_standardized_dataset_collapseOnMethod_class();
            collapse_on_method.Generate_from_kpmp_standardized_dataset(standardized_data_log2Ratios_only_subsegmental);
            collapse_on_method.Write(results_directory, "Standardized_dataset_collapsed_on_method_onlySubsegmental.txt");

            standardized_data_log2Ratios_only_subsegmental.Keep_only_indicated_gene_symbols(subsegmental_shared_genes);
            standardized_data_log2Ratios_only_subsegmental.Write(results_directory, "Standardized_data_log2Ratios_onlySubsegmental_sharedGenesProteins.txt", true);

            collapse_on_method = new KPMP_standardized_dataset_collapseOnMethod_class();
            collapse_on_method.Generate_from_kpmp_standardized_dataset(standardized_data_log2Ratios_only_subsegmental);
            collapse_on_method.Write(results_directory, "Standardized_dataset_collapsed_on_method_onlySubsegmental_sharedGeneProteins.txt");

            KPMP_value_type_enum[] value_types = new KPMP_value_type_enum[] { KPMP_value_type_enum.Log2_ratio_singlepatient, KPMP_value_type_enum.Single_value };
            KPMP_standardized_dataset_class current_standardized_dataset;
            foreach (KPMP_value_type_enum value_type in value_types)
            {
                current_standardized_dataset = standardized_dataset.Deep_copy();
                current_standardized_dataset.Keep_only_indicated_valueTypes_1st(value_type);

                DE_class de_all = current_standardized_dataset.Generate_de_instance_alternatively_and_fill_with_one_of_indicated_valueTypes_and_check_if_duplicated(value_type);
                de_all.Write_file(results_directory, baseFileName + "_" + value_type + "_all_degsdeps.txt");

                DE_class de_all_shared = de_all.Deep_copy();
                de_all_shared.Keep_only_stated_symbols(shared_genes);
                de_all_shared.Write_file(results_directory, baseFileName + "_" + value_type + "_all_shared_degsdeps.txt");
            }
        }
    }
}
