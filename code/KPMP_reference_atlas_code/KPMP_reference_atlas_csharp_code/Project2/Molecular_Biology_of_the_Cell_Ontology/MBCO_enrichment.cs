/*
Copyright 2022.The code was written by Jens Hansen working for the Ravi Iyengar Lab
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
Please acknowledge the MBC Ontology in your publications by citing the following reference:
Jens Hansen, David Meretzky, Simeneh Woldesenbet, Gustavo Stolovitzky, Ravi Iyengar: 
A flexible ontology for inference of emergent whole cell function from relationships between subcellular processes
Sci Rep. 2017 Dec18th
*/

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Drawing;
using Enumerations;
using MBCO_network_and_nodes;
using Mbco_ontology_enrichment;
using MBCO_gene_scp_associations;
using MBCO_data;
using Statistic;
using ReadWrite;
    

namespace Mbco_enrichment
{
 
    class Mbc_enrichment_pipeline_options_class
    {
        #region Fields for general Options
        private Ontology_type_enum ontology;
        public bool Report { get; set; }
        public Ontology_type_enum Ontology
        {
            get { return ontology; }
            set
            {
                if (!ontology.Equals(Ontology_type_enum.E_m_p_t_y)) { throw new Exception(); }
                ontology = value;
            }
        }
        public int Max_columns_per_analysis { get; set; }
        #endregion

        #region Fields for enrichment analysis
        public float Maximum_pvalue_for_standardDynamicEnrichment { get; set; }
        #endregion

        #region Fields for dynamic enrichment analysis
        public float[] Top_quantile_of_scp_interactions_for_dynamicEnrichment_per_level { get; set; }
        public int[] Scp_levels_for_dynamicEnrichment { get; set; }
        public int[][] Numbers_of_merged_scps_for_dynamicEnrichment_per_level { get; set; }
        public int[] Kept_top_predictions_dynamicEnrichment_per_level { get; set; }
        #endregion

        public Mbc_enrichment_pipeline_options_class(Ontology_type_enum ontology)
        {
            #region Parameter for general options
            Report = true;
            Ontology = ontology;
            Max_columns_per_analysis = 100;
            #endregion

            #region Parameter for standard and dynamic enrichment analysis
            Maximum_pvalue_for_standardDynamicEnrichment = 1;
            #endregion

            #region Fields for dynamic enrichment analysis
            Top_quantile_of_scp_interactions_for_dynamicEnrichment_per_level = new float[] { -1, -1, -1, 0.25F, -1 };
            Scp_levels_for_dynamicEnrichment = new int[] { 3 };
            Numbers_of_merged_scps_for_dynamicEnrichment_per_level = new int[5][];
            Numbers_of_merged_scps_for_dynamicEnrichment_per_level[3] = new int[] { 2, 3 };
            Kept_top_predictions_dynamicEnrichment_per_level = new int[] { -1, -1, -1, 7 };
            #endregion
        }

        public Mbc_enrichment_pipeline_options_class Deep_copy()
        {
            Mbc_enrichment_pipeline_options_class copy = (Mbc_enrichment_pipeline_options_class)this.MemberwiseClone();
            copy.Kept_top_predictions_dynamicEnrichment_per_level = Array_class.Deep_copy_array(this.Kept_top_predictions_dynamicEnrichment_per_level);
            int numbers_of_merged_scps_length = this.Numbers_of_merged_scps_for_dynamicEnrichment_per_level.Length;
            copy.Numbers_of_merged_scps_for_dynamicEnrichment_per_level = new int[numbers_of_merged_scps_length][];
            for (int indexNumber = 0; indexNumber < numbers_of_merged_scps_length; indexNumber++)
            {
                copy.Numbers_of_merged_scps_for_dynamicEnrichment_per_level[indexNumber] = Array_class.Deep_copy_array(this.Numbers_of_merged_scps_for_dynamicEnrichment_per_level[indexNumber]);
            }
            copy.Scp_levels_for_dynamicEnrichment = Array_class.Deep_copy_array(this.Scp_levels_for_dynamicEnrichment);
            copy.Top_quantile_of_scp_interactions_for_dynamicEnrichment_per_level = Array_class.Deep_copy_array(this.Top_quantile_of_scp_interactions_for_dynamicEnrichment_per_level);
            return copy;
        }
    }

    class Mbc_enrichment_pipeline_class
    {
        #region Fields
        public MBCO_association_class MBCO_association_unmodified { get; set; }
        public MBCO_association_class MBCO_association { get; set; }
        public Leave_out_scp_scp_network_class Leave_out_scp_network_for_dynamicEnrichment_analysis { get; set; }
        public Mbc_enrichment_pipeline_options_class Options { get; set; }
        public Fisher_exact_test_class Fisher { get; set; }
        public string[] Exp_bg_genes { get; set; }
        public string[] Final_bg_genes { get; set; }
        #endregion

        public Mbc_enrichment_pipeline_class(Ontology_type_enum ontology)
        {
            Options = new Mbc_enrichment_pipeline_options_class(ontology);
            Leave_out_scp_network_for_dynamicEnrichment_analysis = new Leave_out_scp_scp_network_class(this.Options.Ontology);
        }

        #region Generate
        private void Generate_mbco_association()
        {
            MBCO_association_unmodified = new MBCO_association_class();
            MBCO_association_unmodified.Generate_by_reading_safed_file(this.Options.Ontology);
        }

        private void Generate_leave_out_scp_network()
        {
            Leave_out_class leave_out = new Leave_out_class(Options.Ontology);
            leave_out.Generate_by_reading_safed_file();
            Leave_out_scp_network_for_dynamicEnrichment_analysis.Options.Top_quantile_of_considered_SCP_interactions_per_level = this.Options.Top_quantile_of_scp_interactions_for_dynamicEnrichment_per_level;
            Leave_out_scp_network_for_dynamicEnrichment_analysis.Generate_scp_scp_network_from_leave_out(leave_out);
            Leave_out_scp_network_for_dynamicEnrichment_analysis.Scp_nw.Nodes.Set_processLevel_for_all_nodes_to_level3();
        }

        public void Generate_for_all_runs(params MBCO_association_class[] foreign_ontology)
        {
            Generate_mbco_association();
            Generate_leave_out_scp_network();
        }

        public void Reset_mbco_and_adjust_bg_genes(string[] bg_genes)
        {
            this.MBCO_association = this.MBCO_association_unmodified.Deep_copy();
            string[] all_mbco_genes = MBCO_association.Get_all_distinct_ordered_symbols();
            if (bg_genes.Length > 0)
            {
                this.Exp_bg_genes = Array_class.Deep_copy_string_array(bg_genes);
                this.Final_bg_genes = Overlap_class.Get_intersection(this.Exp_bg_genes, all_mbco_genes);
            }
            else
            {
                this.Exp_bg_genes = new string[0];
                this.Final_bg_genes = Array_class.Deep_copy_string_array(all_mbco_genes);
            }
            this.MBCO_association.Keep_only_bg_symbols(this.Final_bg_genes);
            this.MBCO_association.Remove_background_genes_scp();
            Fisher = new Fisher_exact_test_class(Final_bg_genes.Length, false);
        }
        #endregion

        #region Check
        private void Check_if_no_duplicates_and_in_bg_genes(string[] input_overlap_genes)
        {
            string[] overlap_genes = input_overlap_genes.Distinct().OrderBy(l => l).ToArray();
            this.Final_bg_genes = this.Final_bg_genes.OrderBy(l => l).ToArray();
            string bg_gene;
            int bg_genes_length = this.Final_bg_genes.Length;
            int indexBG = 0;
            int stringCompare2 = -2;
            if (overlap_genes.Length != input_overlap_genes.Length) { throw new Exception(); }
            foreach (string gene in overlap_genes)
            {
                stringCompare2 = -2;
                while (stringCompare2 < 0)
                {
                    bg_gene = this.Final_bg_genes[indexBG];
                    stringCompare2 = bg_gene.CompareTo(gene);
                    if (stringCompare2 < 0)
                    {
                        indexBG++;
                    }
                }
                if (stringCompare2 != 0)
                {
                    throw new Exception();
                }
            }
        }

        public void Check_for_duplicated_enrichment_lines(Ontology_enrichment_line_class[] onto_enrich_lines)
        {
            int onto_enrich_length = onto_enrich_lines.Length;
            onto_enrich_lines = onto_enrich_lines.OrderBy(l => l.Sample_name).ThenBy(l => l.Scp_name).ToArray();
            Ontology_enrichment_line_class previous_onto_enrich_line;
            Ontology_enrichment_line_class onto_enrich_line;
            for (int indexO = 1; indexO < onto_enrich_length; indexO++)
            {
                previous_onto_enrich_line = onto_enrich_lines[indexO - 1];
                onto_enrich_line = onto_enrich_lines[indexO];
                if ((previous_onto_enrich_line.Sample_name.Equals(onto_enrich_line.Sample_name))
                    && (previous_onto_enrich_line.Scp_name.Equals(onto_enrich_line.Scp_name)))
                {
                    throw new Exception();
                }
            }
        }
        #endregion

        #region Data analysis
        private Dictionary<string, List<string>[]> Generate_processName_colIndex_overlapSymbols_dictionary_and_count_column_entries_in_data(ref Data_class data)
        {
            int data_empty_value = 0;
            Dictionary<string, List<string>[]> processName_overlapSymbols_dict = new Dictionary<string, List<string>[]>();
            int col_count = data.ColChar.Columns.Length;
            Data_line_class[] data_lines = data.Data;
            data_lines = data_lines.OrderBy(l => l.NCBI_official_symbol.Length).ThenBy(l => l.NCBI_official_symbol).ToArray();
            Data_line_class data_line;
            int data_length = data.Data.Length;

            Data_line_class col_entries_count_line = new Data_line_class("Column entries count", col_count);

            int indexMbco = 0;
            MBCO_association_line_class[] mbco_association_lines = MBCO_association.MBCO_associations;
            mbco_association_lines = mbco_association_lines.OrderBy(l => l.Symbol.Length).ThenBy(l => l.Symbol).ToArray();
            MBCO_association_line_class mbco_association_line;
            int mbco_associations_length = mbco_association_lines.Length;
            int valueCompare;
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                data_line = data_lines[indexData];
                valueCompare = -2;
                while ((indexMbco < mbco_associations_length) && (valueCompare <= 0))
                {
                    mbco_association_line = mbco_association_lines[indexMbco];
                    valueCompare = mbco_association_line.Symbol.Length.CompareTo(data_line.NCBI_official_symbol.Length);
                    if (valueCompare == 0)
                    {
                        valueCompare = String.CompareOrdinal(mbco_association_line.Symbol, data_line.NCBI_official_symbol);
                    }
                    if (valueCompare < 0)
                    {
                        indexMbco++;
                    }
                    else if (valueCompare == 0)
                    {
                        if (!processName_overlapSymbols_dict.ContainsKey(mbco_association_line.ProcessName))
                        {
                            processName_overlapSymbols_dict.Add(mbco_association_line.ProcessName, new List<string>[col_count]);
                            for (int indexCol = 0; indexCol < col_count; indexCol++)
                            {
                                processName_overlapSymbols_dict[mbco_association_line.ProcessName][indexCol] = new List<string>();
                            }
                        }
                        for (int indexCol = 0; indexCol < col_count; indexCol++)
                        {
                            if (data_line.Columns[indexCol] != data_empty_value)
                            {
                                processName_overlapSymbols_dict[mbco_association_line.ProcessName][indexCol].Add(data_line.NCBI_official_symbol);
                            }
                        }
                        indexMbco++;
                    }
                }
                for (int indexCol = 0; indexCol < col_count; indexCol++)
                {
                    if (data_line.Columns[indexCol] != data_empty_value)
                    {
                        col_entries_count_line.Columns[indexCol]++;
                    }
                }
            }
            data.Column_entries_count_line = col_entries_count_line;
            return processName_overlapSymbols_dict;
        }

        private Ontology_enrichment_line_class[] Generate_dynamic_enrichment_lines(Dictionary<string, List<string>[]> processName_colIndex_overlapSymbols_dict, Data_class data)
        {
            List<Ontology_enrichment_line_class> enrichment_lines = new List<Ontology_enrichment_line_class>();
            string[] processNames = processName_colIndex_overlapSymbols_dict.Keys.ToArray();
            string processName;
            int processNames_length = processNames.Length;
            Ontology_enrichment_line_class new_enrichment_line;
            List<Ontology_enrichment_line_class> enrichment_lines_list = new List<Ontology_enrichment_line_class>();

            Colchar_column_line_class[] columns = data.ColChar.Columns;
            Colchar_column_line_class column_line;
            Data_line_class data_column_entries_line = data.Column_entries_count_line;
            int col_length = columns.Length;
            List<string>[] overlap_symbols_of_each_column;
            for (int indexPN = 0; indexPN < processNames_length; indexPN++)
            {
                processName = processNames[indexPN];
                overlap_symbols_of_each_column = processName_colIndex_overlapSymbols_dict[processName];
                for (int indexCol = 0; indexCol < col_length; indexCol++)
                {
                    if (overlap_symbols_of_each_column[indexCol].Count > 0)
                    {
                        column_line = columns[indexCol];
                        new_enrichment_line = new Ontology_enrichment_line_class();
                        new_enrichment_line.Ontology_type = this.Options.Ontology;
                        new_enrichment_line.Scp_name = (string)processName.Clone();
                        new_enrichment_line.Sample_name = column_line.SampleName;
                        new_enrichment_line.Overlap_symbols = overlap_symbols_of_each_column[indexCol].ToArray();
                        new_enrichment_line.Overlap_count = new_enrichment_line.Overlap_symbols.Length;
                        new_enrichment_line.Experimental_symbols_count = (int)data_column_entries_line.Columns[indexCol];
                        enrichment_lines_list.Add(new_enrichment_line);
                    }
                }
            }
            return enrichment_lines_list.ToArray();
        }

        private void Add_missing_process_information_and_backgroundGenes_count_for_standard_scps(ref Ontology_enrichment_line_class[] enrichment_lines)
        {
            int background_genes_length = this.Final_bg_genes.Length;
            enrichment_lines = enrichment_lines.OrderBy(l => l.Scp_name).ToArray();
            int enrich_length = enrichment_lines.Length;
            Ontology_enrichment_line_class enrich_line;
            MBCO_association.Order_by_processName_symbol();
            int onto_length = MBCO_association.MBCO_associations.Length;
            int indexOnto = 0;
            MBCO_association_line_class mbco_association_line;
            int stringCompare;
            int process_symbols_count = 1;
            for (int indexE = 0; indexE < enrich_length; indexE++)
            {
                enrich_line = enrichment_lines[indexE];
                stringCompare = -2;
                while ((indexOnto < onto_length) && (stringCompare < 0))
                {
                    mbco_association_line = MBCO_association.MBCO_associations[indexOnto];
                    stringCompare = mbco_association_line.ProcessName.CompareTo(enrich_line.Scp_name);
                    if (stringCompare < 0)
                    {
                        indexOnto++;
                        process_symbols_count = 1;
                    }
                    else if (stringCompare == 0)
                    {
                        while ((indexOnto < onto_length - 1) && (mbco_association_line.ProcessName.Equals(MBCO_association.MBCO_associations[indexOnto + 1].ProcessName)))
                        {
                            process_symbols_count++;
                            indexOnto++;
                        }
                        if (!string.IsNullOrEmpty(mbco_association_line.ProcessID)) { enrich_line.Scp_id = (string)mbco_association_line.ProcessID.Clone(); }
                        enrich_line.ProcessLevel = mbco_association_line.ProcessLevel;
                        if (!string.IsNullOrEmpty(mbco_association_line.Parent_processName)) { enrich_line.Parent_scp_name = (string)mbco_association_line.Parent_processName.Clone(); }
                        enrich_line.Process_symbols_count = process_symbols_count;
                        enrich_line.Bg_symbol_count = background_genes_length;
                    }
                }
            }
        }

        private void Add_missing_process_information_and_backgroundGenes_count_for_scp_unions(ref Ontology_enrichment_line_class[] enrichment_lines, Dictionary<string, List<string>> scp_scpUnion_dict)
        {
            MBCO_association.Order_by_symbol_processName();
            Dictionary<string, int> scpUnion_symbols_length_dict = new Dictionary<string, int>();
            MBCO_association_line_class mbco_association_line;
            int mbco_association_length = MBCO_association.MBCO_associations.Length;
            List<string> scp_unions_of_current_scp;
            List<string> add_one_to_gene_count_of_scp_unions = new List<string>();
            Dictionary<string, int> scpUnion_genes_count_dict = new Dictionary<string, int>();
            Dictionary<string, int> scpUnion_level_dict = new Dictionary<string, int>();
            for (int indexMBCO = 0; indexMBCO < mbco_association_length; indexMBCO++)
            {
                mbco_association_line = MBCO_association.MBCO_associations[indexMBCO];
                if ((indexMBCO == 0) || (!mbco_association_line.Symbol.Equals(MBCO_association.MBCO_associations[indexMBCO - 1].Symbol)))
                {
                    add_one_to_gene_count_of_scp_unions.Clear();
                }
                if (scp_scpUnion_dict.ContainsKey(mbco_association_line.ProcessName))
                {
                    scp_unions_of_current_scp = scp_scpUnion_dict[mbco_association_line.ProcessName];
                    add_one_to_gene_count_of_scp_unions.AddRange(scp_unions_of_current_scp);
                    foreach (string scp_union_of_curren_scp in scp_unions_of_current_scp)
                    {
                        if (!scpUnion_level_dict.ContainsKey(scp_union_of_curren_scp))
                        {
                            scpUnion_level_dict.Add(scp_union_of_curren_scp, mbco_association_line.ProcessLevel);
                        }
                    }
                }
                if ((indexMBCO == mbco_association_length - 1) || (!mbco_association_line.Symbol.Equals(MBCO_association.MBCO_associations[indexMBCO + 1].Symbol)))
                {
                    add_one_to_gene_count_of_scp_unions = add_one_to_gene_count_of_scp_unions.Distinct().ToList();
                    foreach (string scp_union in add_one_to_gene_count_of_scp_unions)
                    {
                        if (!scpUnion_genes_count_dict.ContainsKey(scp_union)) { scpUnion_genes_count_dict.Add(scp_union, 0); }
                        scpUnion_genes_count_dict[scp_union]++;
                    }
                }
            }

            int background_symbols_count = this.Final_bg_genes.Length;
            int enrichment_lines_length = enrichment_lines.Length;
            Ontology_enrichment_line_class enrichment_line;
            for (int indexE = 0; indexE < enrichment_lines_length; indexE++)
            {
                enrichment_line = enrichment_lines[indexE];
                if (scpUnion_genes_count_dict.ContainsKey(enrichment_line.Scp_name))
                {
                    enrichment_line.Process_symbols_count = scpUnion_genes_count_dict[enrichment_line.Scp_name];
                    enrichment_line.ProcessLevel = scpUnion_level_dict[enrichment_line.Scp_name];
                    enrichment_line.Bg_symbol_count = background_symbols_count;
                    enrichment_line.Parent_scp_name = "No annotated parent SCP";
                    enrichment_line.Scp_id = "Dynamic SCP union";
                }
            }
        }

        private Dictionary<string, List<string>> Generate_all_scp_unions_and_add_them_with_overlap_symbols_to_scp_colIndex_overlapSymbols_dict(ref Dictionary<string, List<string>[]> processName_colIndex_overlapSymbols_dict)
        {
            string[] experimental_processNames = processName_colIndex_overlapSymbols_dict.Keys.ToArray();
            string[][] array_of_scps_in_one_union = Leave_out_scp_network_for_dynamicEnrichment_analysis.Generate_array_of_scp_unions_between_any_combination_between_two_or_three_neighboring_selected_scps(experimental_processNames);

            Dictionary<string, List<string>> scp_scpUnion_dict = new Dictionary<string, List<string>>();

            string[] consideredSCP_names = processName_colIndex_overlapSymbols_dict.Keys.ToArray();
            int data_col_length = processName_colIndex_overlapSymbols_dict[consideredSCP_names[0]].Length;

            int unions_length = array_of_scps_in_one_union.Length;
            string[] scps_of_one_union;
            string scp;
            int scps_of_one_union_length;
            StringBuilder sb = new StringBuilder();
            string name_of_scp_union;
            List<string>[] unionScpOverlapSymbols_of_each_column = new List<string>[data_col_length];
            List<string>[] currentScpOverlapSymbols_of_each_column;
            bool[] consider_column = new bool[data_col_length];
            for (int indexUnion = 0; indexUnion < unions_length; indexUnion++)
            {
                scps_of_one_union = array_of_scps_in_one_union[indexUnion];
                scps_of_one_union_length = scps_of_one_union.Length;

                #region Clear stringbuilder and reset unionScpOverlapSymbols_of_each_column and consider scp
                sb.Clear();
                unionScpOverlapSymbols_of_each_column = new List<string>[data_col_length];
                for (int indexCol = 0; indexCol < data_col_length; indexCol++)
                {
                    unionScpOverlapSymbols_of_each_column[indexCol] = new List<string>();
                    consider_column[indexCol] = true;
                }
                #endregion

                for (int indexScp = 0; indexScp < scps_of_one_union_length; indexScp++)
                {
                    scp = scps_of_one_union[indexScp];
                    if (indexScp != 0) { sb.AppendFormat("$"); }
                    sb.AppendFormat(scp);
                    currentScpOverlapSymbols_of_each_column = processName_colIndex_overlapSymbols_dict[scp];
                    for (int indexCol = 0; indexCol < data_col_length; indexCol++)
                    {
                        if (currentScpOverlapSymbols_of_each_column[indexCol].Count == 0)
                        {
                            consider_column[indexCol] = false;
                            unionScpOverlapSymbols_of_each_column[indexCol].Clear();
                        }
                        else if (consider_column[indexCol])
                        {
                            unionScpOverlapSymbols_of_each_column[indexCol].AddRange(currentScpOverlapSymbols_of_each_column[indexCol]);
                        }
                    }
                }
                for (int indexCol = 0; indexCol < data_col_length; indexCol++)
                {
                    unionScpOverlapSymbols_of_each_column[indexCol] = unionScpOverlapSymbols_of_each_column[indexCol].Distinct().OrderBy(l => l).ToList();
                }
                name_of_scp_union = sb.ToString();
                for (int indexScp = 0; indexScp < scps_of_one_union_length; indexScp++)
                {
                    scp = scps_of_one_union[indexScp];
                    if (!scp_scpUnion_dict.ContainsKey(scp)) { scp_scpUnion_dict.Add(scp, new List<string>()); }
                    scp_scpUnion_dict[scp].Add(name_of_scp_union);
                }
                for (int indexCol = 0; indexCol < data_col_length; indexCol++)
                {
                    if (consider_column[indexCol])
                    {
                        processName_colIndex_overlapSymbols_dict.Add(name_of_scp_union, unionScpOverlapSymbols_of_each_column);
                        break;
                    }
                }
            }
            return scp_scpUnion_dict;
        }

        private void Calculate_missing_pvalues(ref Ontology_enrichment_line_class[] enrichment_lines)
        {
            int enrich_length = enrichment_lines.Length;
            Ontology_enrichment_line_class enrich_line;
            int a; int b; int c; int d;
            for (int indexE = 0; indexE < enrich_length; indexE++)
            {
                enrich_line = enrichment_lines[indexE];
                if (enrich_line.Pvalue == -1)
                {
                    a = enrich_line.Overlap_count;
                    b = enrich_line.Experimental_symbols_count - a;
                    c = enrich_line.Process_symbols_count - a;
                    d = enrich_line.Bg_symbol_count - a - b - c;
                    if ((a < 0) || (b < 0) || (c < 0) || (d < 0))
                    {
                        Check_if_no_duplicates_and_in_bg_genes(enrich_line.Overlap_symbols);
                        string[] genes = MBCO_association.Get_all_symbols_of_process_names(enrich_line.Scp_name);
                        Check_if_no_duplicates_and_in_bg_genes(genes);
                        //Check_if_no_duplicates_and_in_bg_genes(this.Exp_bg_genes);
                        throw new Exception();
                    }
                    enrich_line.Pvalue = Fisher.Get_rightTailed_p_value(a, b, c, d);
                    if (enrich_line.Pvalue > 1)
                    {
                        if (enrich_line.Pvalue < 1.0001)
                        {
                            enrich_line.Pvalue = 1;
                        }
                        else
                        {
                            throw new Exception();
                        }
                    }
                }
            }
        }

        private void Calculate_minusLog10pvalues(ref Ontology_enrichment_line_class[] enrichment_lines)
        {
            int enrich_length = enrichment_lines.Length;
            Ontology_enrichment_line_class enrich_line;
            for (int indexE = 0; indexE < enrich_length; indexE++)
            {
                enrich_line = enrichment_lines[indexE];
                enrich_line.Minus_log10_pvalue = -(float)Math.Log10(enrich_line.Pvalue);
            }
        }

        public Ontology_enrichment_class Analyse_data_instance_via_dynamic_enrichment_analysis(Data_class data_input)
        {
            int columns_length = data_input.ColChar.Columns_length;
            int runs = (int)Math.Ceiling((float)columns_length / (float)Options.Max_columns_per_analysis);
            int currentFirstColumn = -1;
            int currentLastColumn = -1;
            List<int> keep_columns = new List<int>();
            List<Ontology_enrichment_line_class> dynamic_enrichment_lines_list = new List<Ontology_enrichment_line_class>();
            for (int indexRun = 0; indexRun < runs; indexRun++)
            {
                currentFirstColumn = indexRun * Options.Max_columns_per_analysis;
                currentLastColumn = Math.Min((indexRun + 1) * Options.Max_columns_per_analysis - 1, columns_length - 1);
                keep_columns.Clear();
                for (int keepColumn = currentFirstColumn; keepColumn <= currentLastColumn; keepColumn++)
                {
                    keep_columns.Add(keepColumn);
                }
                Data_class data = data_input.Deep_copy();
                data.Keep_only_input_columns_and_remove_all_rows_that_are_left_over_with_only_zero_values(keep_columns.ToArray());
                data.Keep_only_input_rowNames(this.Final_bg_genes);
                Dictionary<string, List<string>[]> scp_colIndex_overlapSymbols_dict = Generate_processName_colIndex_overlapSymbols_dictionary_and_count_column_entries_in_data(ref data);
                Dictionary<string, List<string>> scp_scpUnion_dict = Generate_all_scp_unions_and_add_them_with_overlap_symbols_to_scp_colIndex_overlapSymbols_dict(ref scp_colIndex_overlapSymbols_dict);
                Ontology_enrichment_line_class[] current_dynamic_enrichment_lines = Generate_dynamic_enrichment_lines(scp_colIndex_overlapSymbols_dict, data);
                Check_for_duplicated_enrichment_lines(current_dynamic_enrichment_lines);
                Add_missing_process_information_and_backgroundGenes_count_for_standard_scps(ref current_dynamic_enrichment_lines);
                Check_for_duplicated_enrichment_lines(current_dynamic_enrichment_lines);
                Add_missing_process_information_and_backgroundGenes_count_for_scp_unions(ref current_dynamic_enrichment_lines, scp_scpUnion_dict);
                Check_for_duplicated_enrichment_lines(current_dynamic_enrichment_lines);
                Calculate_missing_pvalues(ref current_dynamic_enrichment_lines);
                Check_for_duplicated_enrichment_lines(current_dynamic_enrichment_lines);
                dynamic_enrichment_lines_list.AddRange(current_dynamic_enrichment_lines);
            }
            Ontology_enrichment_line_class[] dynamic_enrichment_lines = dynamic_enrichment_lines_list.ToArray();

            Calculate_minusLog10pvalues(ref dynamic_enrichment_lines);
            Check_for_duplicated_enrichment_lines(dynamic_enrichment_lines);
            Ontology_enrichment_class dynamic_enrichment_filtered = new Ontology_enrichment_class();
            dynamic_enrichment_filtered.Add_other_lines(dynamic_enrichment_lines);
            dynamic_enrichment_filtered.Keep_top_ranked_predictions_per_level_for_each_sample_after_calculation_of_fractional_rank(Options.Kept_top_predictions_dynamicEnrichment_per_level);
            return dynamic_enrichment_filtered;
        }
        #endregion
    }

    class MBCO_network_based_integration_options_class
    {
        public bool Report { get; set; }
        public bool Add_edges_that_connect_scps_between_sets { get; set; }
        public float[] Top_quantile_probability_of_scp_interactions_for_dynamic_enrichment { get; set; }
        public Ontology_type_enum Ontology { get; set; }

        public int[] Keep_top_prediction_for_randomWalkWithRestart_based_networks_per_level { get; set; }

        public MBCO_network_based_integration_options_class(Ontology_type_enum ontology)
        {
            Ontology = ontology;
            Report = true;
            Add_edges_that_connect_scps_between_sets = false;
            Keep_top_prediction_for_randomWalkWithRestart_based_networks_per_level = new int[] { -1, -1, 5, 5, -1 };
        }
    }

    class Mbc_network_based_integration_class
    {
        #region Fields
        public Leave_out_scp_scp_network_class Leave_out_scp_network_dynamic_scps { get; set; }
        public MBCO_network_based_integration_options_class Options { get; set; }
        #endregion

        public Mbc_network_based_integration_class(Ontology_type_enum ontology)
        {
            this.Options = new MBCO_network_based_integration_options_class(ontology);
            Leave_out_scp_network_dynamic_scps = new Leave_out_scp_scp_network_class(this.Options.Ontology);
        }

        #region Generate
        private void Generate_leave_out_scp_networks()
        {
            Leave_out_class leave_out = new Leave_out_class(Options.Ontology);
            leave_out.Generate_by_reading_safed_file();

            Leave_out_scp_network_dynamic_scps.Options.Top_quantile_of_considered_SCP_interactions_per_level = this.Options.Top_quantile_probability_of_scp_interactions_for_dynamic_enrichment;
            Leave_out_scp_network_dynamic_scps.Generate_scp_scp_network_from_leave_out(leave_out);
            Leave_out_scp_network_dynamic_scps.Scp_nw.Transform_into_undirected_single_network_and_set_all_widths_to_one();
            Leave_out_scp_network_dynamic_scps.Scp_nw.Nodes.Set_processLevel_for_all_nodes_to_level3();
        }

        public void Generate()
        {
            Generate_leave_out_scp_networks();
        }
        #endregion

        public void Generate_and_write_integrative_network_for_dynamic_enrichment_results_of_each_integrationGroupName_only_defined_sets(Ontology_enrichment_class dynamic_onto_enrich, Dictionary<string, Color> setAbbreviation_color_dict, Dictionary<string, string> set_setAbbreviation_dict, string baseFile_name, Dictionary<string, Network_class> integrationGroup_addNetwork_dict)
        {
            string results_directory = Global_directory_class.MBCO_scp_networks_results_directory;
            Ontology_enrichment_line_class onto_enrich_unfiltered_line = new Ontology_enrichment_line_class();
            ReadWriteClass.Create_directory_if_it_does_not_exist(results_directory);
            dynamic_onto_enrich.Enrich = dynamic_onto_enrich.Enrich.OrderBy(l => l.Integration_group).ThenBy(l => l.Sample_name).ToArray();
            Ontology_enrichment_class current_standard_onto_enrich = new Ontology_enrichment_class();
            List<Ontology_enrichment_line_class> sameEntryTypeTimepoint_onto_enrich_list = new List<Ontology_enrichment_line_class>();
            int standard_onto_enrich_length = dynamic_onto_enrich.Enrich.Length;
            Ontology_enrichment_line_class onto_enrich_line;
            string fileName;

            Leave_out_scp_scp_network_class current_scp_network = new Leave_out_scp_scp_network_class(this.Options.Ontology);
            Leave_out_scp_scp_network_class current_enrichmentLine_scp_network;
            string[] scps;
            List<string> scps_of_current_set = new List<string>();

            string[] scps_for_dictionary = dynamic_onto_enrich.Get_all_distinct_scps_after_spliting_scp_unions();
            Dictionary<string, int> processName_processLevel_dict = new Dictionary<string, int>();
            foreach (string scp in scps_for_dictionary)
            {
                processName_processLevel_dict.Add(scp, 3);
            }

            #region Variables for same integration group
            Network_class combined_scpNetwork_ofCurrentIntegrationGroup = new Network_class();
            Dictionary<string, List<string>> scpsAndLegends_setAbbreviations_dict_ofCurrentIntegrationGroup = new Dictionary<string, List<string>>();
            int set_count_ofCurrentIntegrationGroup = 0;
            int sample_count_within_currentSetOfCurrentIntegrationGroup = 0;
            List<string> sets_ofCurrentIntegrationGroup = new List<string>();
            #endregion

            string[] missing_scps_and_legend;
            string setAbbreviation = "";
            string[] setAbbrevations_of_current_scp;
            string legend_name = "";
            List<Color> colors_of_current_scp = new List<Color>();

            List<Ontology_enrichment_line_class> current_onto_enrich = new List<Ontology_enrichment_line_class>();
            Dictionary<string, int> setAbbreviation_processLevel_dict = new Dictionary<string, int>();

            for (int indexO = 0; indexO < standard_onto_enrich_length; indexO++)
            {
                onto_enrich_line = dynamic_onto_enrich.Enrich[indexO];
                if ((indexO == 0)
                    || (!onto_enrich_line.Integration_group.Equals(dynamic_onto_enrich.Enrich[indexO - 1].Integration_group)))
                {
                    set_count_ofCurrentIntegrationGroup = 0;
                    combined_scpNetwork_ofCurrentIntegrationGroup = new Network_class();
                    scpsAndLegends_setAbbreviations_dict_ofCurrentIntegrationGroup.Clear();
                    sets_ofCurrentIntegrationGroup.Clear();
                }
                if ((indexO == 0)
                    || (!onto_enrich_line.Set_within_integration_group.Equals(dynamic_onto_enrich.Enrich[indexO - 1].Set_within_integration_group))
                    || (!onto_enrich_line.Integration_group.Equals(dynamic_onto_enrich.Enrich[indexO - 1].Integration_group)))
                {
                    sample_count_within_currentSetOfCurrentIntegrationGroup = 0;
                }
                if ((indexO == 0)
                    || (!onto_enrich_line.Integration_group.Equals(dynamic_onto_enrich.Enrich[indexO - 1].Integration_group))
                    || (!onto_enrich_line.Sample_name.Equals(dynamic_onto_enrich.Enrich[indexO - 1].Sample_name)))
                {
                    current_scp_network = new Leave_out_scp_scp_network_class(this.Options.Ontology);
                    set_count_ofCurrentIntegrationGroup++;
                    sample_count_within_currentSetOfCurrentIntegrationGroup++;
                    scps_of_current_set.Clear();
                    legend_name = onto_enrich_line.Complete_sample_name;
                    if (sample_count_within_currentSetOfCurrentIntegrationGroup == 1)
                    {
                        setAbbreviation = set_setAbbreviation_dict[onto_enrich_line.Set_within_integration_group];
                    }
                    else
                    {
                        setAbbreviation = set_setAbbreviation_dict[onto_enrich_line.Set_within_integration_group] + sample_count_within_currentSetOfCurrentIntegrationGroup;
                    }
                    if (!scpsAndLegends_setAbbreviations_dict_ofCurrentIntegrationGroup.ContainsKey(legend_name))
                    {
                        scpsAndLegends_setAbbreviations_dict_ofCurrentIntegrationGroup.Add(legend_name, new List<string>());
                    }
                    scpsAndLegends_setAbbreviations_dict_ofCurrentIntegrationGroup[legend_name].Add(setAbbreviation);
                    sets_ofCurrentIntegrationGroup.Add(onto_enrich_line.Set_within_integration_group);
                }
                scps = onto_enrich_line.Scp_name.Split('$');
                scps = scps.OrderBy(l => l).ToArray();
                foreach (string scp in scps)
                {
                    if (!scpsAndLegends_setAbbreviations_dict_ofCurrentIntegrationGroup.ContainsKey(scp)) { scpsAndLegends_setAbbreviations_dict_ofCurrentIntegrationGroup.Add(scp, new List<string>()); }
                    scpsAndLegends_setAbbreviations_dict_ofCurrentIntegrationGroup[scp].Add(setAbbreviation);
                    scpsAndLegends_setAbbreviations_dict_ofCurrentIntegrationGroup[scp] = scpsAndLegends_setAbbreviations_dict_ofCurrentIntegrationGroup[scp].Distinct().ToList();
                }

                current_enrichmentLine_scp_network = Leave_out_scp_network_dynamic_scps.Deep_copy_scp_network();
                current_enrichmentLine_scp_network.Scp_nw.Keep_only_input_nodeNames(scps.ToArray());
                current_scp_network.Scp_nw.Merge_this_network_with_other_network(current_enrichmentLine_scp_network.Scp_nw);
                scps_of_current_set.AddRange(scps);

                if ((indexO == standard_onto_enrich_length - 1)
                    || (!onto_enrich_line.Integration_group.Equals(dynamic_onto_enrich.Enrich[indexO + 1].Integration_group))
                    || (!onto_enrich_line.Sample_name.Equals(dynamic_onto_enrich.Enrich[indexO + 1].Sample_name)))
                {
                    missing_scps_and_legend = Overlap_class.Get_part_of_list1_but_not_of_list2(scps_of_current_set.ToArray(), current_scp_network.Scp_nw.Nodes.Get_all_nodeNames());
                    if (missing_scps_and_legend.Length > 0)
                    {
                        current_scp_network.Scp_nw.Add_single_nodes(missing_scps_and_legend);
                    }

                    current_scp_network.Scp_nw.Nodes.Set_level_for_all_nodes(10);
                    current_scp_network.Scp_nw.Nodes.Set_processLevel_for_all_nodes_based_on_dictionary(processName_processLevel_dict);
                    current_scp_network.Scp_nw.Add_single_nodes(legend_name);

                    combined_scpNetwork_ofCurrentIntegrationGroup.Merge_this_network_with_other_network(current_scp_network.Scp_nw);
                }

                if ((indexO == standard_onto_enrich_length - 1)
                    || (!onto_enrich_line.Integration_group.Equals(dynamic_onto_enrich.Enrich[indexO + 1].Integration_group)))
                {
                    if (set_count_ofCurrentIntegrationGroup >= 1)
                    {
                        fileName = results_directory + baseFile_name + "_" + onto_enrich_line.Integration_group;

                        if (integrationGroup_addNetwork_dict.ContainsKey(onto_enrich_line.Integration_group))
                        {
                            Network_class outer_network = integrationGroup_addNetwork_dict[onto_enrich_line.Integration_group].Deep_copy();
                            string[] outer_scps = outer_network.Get_all_scps();
                            outer_network.Merge_this_network_with_other_network(combined_scpNetwork_ofCurrentIntegrationGroup);
                            outer_network.Add_single_nodes("Added network");

                            combined_scpNetwork_ofCurrentIntegrationGroup = outer_network;

                            setAbbreviation = "X";
                            scpsAndLegends_setAbbreviations_dict_ofCurrentIntegrationGroup.Add("Added network", new List<string>());
                            scpsAndLegends_setAbbreviations_dict_ofCurrentIntegrationGroup["Added network"].Add(setAbbreviation);
                            foreach (string scp in outer_scps)
                            {
                                if (!scpsAndLegends_setAbbreviations_dict_ofCurrentIntegrationGroup.ContainsKey(scp)) { scpsAndLegends_setAbbreviations_dict_ofCurrentIntegrationGroup.Add(scp, new List<string>()); }
                                scpsAndLegends_setAbbreviations_dict_ofCurrentIntegrationGroup[scp].Add(setAbbreviation);
                                scpsAndLegends_setAbbreviations_dict_ofCurrentIntegrationGroup[scp] = scpsAndLegends_setAbbreviations_dict_ofCurrentIntegrationGroup[scp].Distinct().ToList();
                                if (!processName_processLevel_dict.ContainsKey(scp))
                                {
                                    processName_processLevel_dict.Add(scp,3);
                                }
                            }
                        }

                        string[] all_scps = combined_scpNetwork_ofCurrentIntegrationGroup.Nodes.Get_all_nodeNames();

                        #region Add_edges_that_connect_scps_between_sets)
                        current_scp_network = Leave_out_scp_network_dynamic_scps.Deep_copy_scp_network();
                        current_scp_network.Scp_nw.Keep_only_input_nodeNames(all_scps);
                        current_scp_network.Scp_nw.Switch_all_edge_types_to_input_type(NWedge_type_enum.Dashed_line);
                        current_scp_network.Scp_nw.Merge_this_network_with_other_network(combined_scpNetwork_ofCurrentIntegrationGroup);
                        combined_scpNetwork_ofCurrentIntegrationGroup = current_scp_network.Scp_nw;
                        #endregion

                        combined_scpNetwork_ofCurrentIntegrationGroup.Nodes.Set_level_for_all_nodes(100);
                        combined_scpNetwork_ofCurrentIntegrationGroup.Nodes.Set_processLevel_for_all_nodes_based_on_dictionary(processName_processLevel_dict);
                        combined_scpNetwork_ofCurrentIntegrationGroup.Keep_only_nodes_with_indicated_levels(new int[] { 3, 10, 100 });
                        string[] scps_100 = combined_scpNetwork_ofCurrentIntegrationGroup.Nodes.Get_all_nodeNames_of_indicated_levels(100);


                        string[] all_nodes_in_combined_scp_network = combined_scpNetwork_ofCurrentIntegrationGroup.Nodes.Get_all_nodeNames();

                        Dictionary<string, Shape_enum> scp_gene_shape_dict = new Dictionary<string, Shape_enum>();
                        Dictionary<string, Color[]> scp_gene_colors_dict = new Dictionary<string, Color[]>();
                        string[] all_scps_legends = scpsAndLegends_setAbbreviations_dict_ofCurrentIntegrationGroup.Keys.OrderBy(l => l).ToArray();
                        string scp_plus_sets_name;
                        string[] sets_of_current_scp;
                        combined_scpNetwork_ofCurrentIntegrationGroup.Nodes.Order_by_name();
                        int nodes_length = combined_scpNetwork_ofCurrentIntegrationGroup.Nodes.Nodes_length;
                        int indexNode = 0;
                        int stringCompare;
                        NetworkNode_line_class node_line;
                        foreach (string current_scpOrLegend in all_scps_legends)
                        {
                            sets_of_current_scp = scpsAndLegends_setAbbreviations_dict_ofCurrentIntegrationGroup[current_scpOrLegend].OrderBy(l => l).ToArray();
                            int sets_length = sets_of_current_scp.Length;
                            if (sets_length == 1)
                            {
                                scp_plus_sets_name = current_scpOrLegend;
                            }
                            else
                            {
                                scp_plus_sets_name = current_scpOrLegend;
                            }

                            stringCompare = -2;
                            while ((indexNode < nodes_length) && (stringCompare < 0))
                            {
                                node_line = combined_scpNetwork_ofCurrentIntegrationGroup.Nodes.Nodes[indexNode];
                                stringCompare = node_line.Name.CompareTo(current_scpOrLegend);
                                if (stringCompare < 0)
                                {
                                    indexNode++;
                                }
                                else if (stringCompare == 0)
                                {
                                    node_line.Name = (string)scp_plus_sets_name.Clone();
                                    indexNode++;
                                }
                            }
                            if (stringCompare != 0) { throw new Exception(); }
                            setAbbrevations_of_current_scp = scpsAndLegends_setAbbreviations_dict_ofCurrentIntegrationGroup[current_scpOrLegend].Distinct().ToArray();
                            colors_of_current_scp.Clear();
                            foreach (string setAbbreviation2 in setAbbrevations_of_current_scp)
                            {
                                colors_of_current_scp.Add(setAbbreviation_color_dict[setAbbreviation2.Substring(0, 1)]);
                            }
                            switch (setAbbrevations_of_current_scp.Length)
                            {
                                case 0:
                                    throw new Exception();
                                case 1:
                                    scp_gene_colors_dict.Add(scp_plus_sets_name, colors_of_current_scp.ToArray());
                                    scp_gene_shape_dict.Add(scp_plus_sets_name, Shape_enum.Ellipse);
                                    break;
                                default:
                                    scp_gene_colors_dict.Add(scp_plus_sets_name, colors_of_current_scp.ToArray());
                                    scp_gene_shape_dict.Add(scp_plus_sets_name, Shape_enum.Ellipse);
                                    break;
                            }
                        }
                        combined_scpNetwork_ofCurrentIntegrationGroup.Write_yED_nw_in_results_directory_with_nodes_colored_by_set_and_sized_by_number_of_different_colors_and_sameLevel_processes_grouped(fileName, Shape_enum.Diamond, scp_gene_shape_dict, scp_gene_colors_dict);
                    }
                }
            }
        }

    }
}
