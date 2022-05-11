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

using Enumerations;
using Highthroughput_data;
using Statistic;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Enrichment_2018
{
    class Ontology_abbreviation_class
    {
        public static string Get_abbreviation_of_ontology(Ontology_type_enum ontology)
        {
            string delimiter = "";
            switch (ontology)
            {
                case Ontology_type_enum.Mbco_level3:
                    return "MBCO2021" + delimiter + "L3";
                default:
                    return ontology.ToString();
            }
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class Downstream_analysis_2020_options_class
    {
        public Ontology_type_enum[] Ontologies { get; set; }
        public Data_value_signs_of_interest_enum[] Data_value_signs_of_interest { get; set; }
        public int Top_top_x_predictions { get; set; }
        public Organism_enum Organism { get; set; }
        public bool Write_results { get; set; }
        public bool Only_keep_filtered_results { get; set; }
        public bool Add_missing_scps_identified_in_other_conditions { get; set; }
        public bool Use_ontology_abbreviations_for_directories_and_files { get; set; }

        public Downstream_analysis_2020_options_class()
        {
            Organism = Global_class.Organism;
            Ontologies = new Ontology_type_enum[] { Ontology_type_enum.Metabolic_energy_generation_pathways };
            Data_value_signs_of_interest = new Data_value_signs_of_interest_enum[] { Data_value_signs_of_interest_enum.Combined, Data_value_signs_of_interest_enum.Upregulated, Data_value_signs_of_interest_enum.Downregulated };
            Top_top_x_predictions = 15;
            Write_results = true;
            Only_keep_filtered_results = false;
            Use_ontology_abbreviations_for_directories_and_files = true;
        }
    }

    class Downstream_analysis_2020_class
    {
        public Fisher_exact_test_class Fisher { get; set; }
        public Ontology_library_class[] Ontology_libraries { get; set; }
        string[] Experimental_bg_genes { get; set; }
        public Downstream_analysis_2020_options_class Options { get; set; }

        public Downstream_analysis_2020_class()
        {
            Options = new Downstream_analysis_2020_options_class();
        }

        #region Generate for all runs
        private string[] Read_background_genes()
        {
            string complete_go_fileName = Global_directory_class.Go_dataset_directory + "Geneontology_BiologicalProcess_2018_2018July17.txt";
            System.IO.StreamReader go_reader = new System.IO.StreamReader(complete_go_fileName);
            List<string> allGenes = new List<string>();
            char delimiter = '\t';
            string input_line;
            string[] splitStrings;
            int splitStrings_length;
            while ((input_line = go_reader.ReadLine())!=null)
            {
                splitStrings = input_line.Split(delimiter);
                splitStrings_length = splitStrings.Length;
                for (int indexSplit=1; indexSplit<splitStrings_length;indexSplit++)
                {
                    allGenes.Add(splitStrings[indexSplit]);
                }
            }

            allGenes.AddRange(new string[] { "ERO1L", "ERO1LB", "FYB", "GBP4", "HSFY2", "IFNA13", "KCNE1L", "KIRREL", "KRT33A", "KRT33B", "LILRA6", "MRE11A", "MTCP1", "NOMO1", "NRD1", "NUPL1", "NUPR1L", "PAK7", "PYCRL", "SLX1B",
                              "SRPR", "TOMM70A" }); //These genes are part of additional pathways in the intially used KPMP ontology that were removed before publication of the metabolic gene library.
                                                    //They are added here, because the enrichment algorithm consideres all pathways genes (no matter, if they are labeled as background genes or
                                                    //annotated to a real pathway) as background genes.

            return allGenes.Distinct().OrderBy(l => l).ToArray();
        }


        private void Generate_ontology_libraries()
        {
            string[] ontology_background_genes = Read_background_genes();
            System.IO.StreamWriter write = new System.IO.StreamWriter(Global_directory_class.Get_report_fileName(),true);
            write.WriteLine("{0} background genes", ontology_background_genes.Length);

            List<Ontology_library_line_class> background_gene_lines = new List<Ontology_library_line_class>();
            Ontology_library_line_class new_library_line;

            Ontology_type_enum[] ontologies = this.Options.Ontologies.Distinct().ToArray();
            Ontology_type_enum ontology;
            int ontologies_length = ontologies.Length;
            Ontology_libraries = new Ontology_library_class[ontologies_length];
            Ontology_library_class ontology_library;
            Dictionary<string, bool> bg_genes_dict = new Dictionary<string, bool>();
            for (int indexO = 0; indexO < ontologies_length; indexO++)
            {
                ontology = ontologies[indexO];
                ontology_library = new Ontology_library_class();
                ontology_library.Generate_by_reading(ontology, Options.Organism, this.Experimental_bg_genes);
                background_gene_lines.Clear();
                foreach (string bgGene in ontology_background_genes)
                {
                    if (!bg_genes_dict.ContainsKey(bgGene))
                    {
                        bg_genes_dict.Add(bgGene, true);
                    }
                    else
                    {
                        write.WriteLine("Duplicated bg gene: {0}", bgGene);
                    }
                    new_library_line = new Ontology_library_line_class();
                    new_library_line.Ontology = ontology;
                    new_library_line.Scp = "Background genes";
                    new_library_line.Target_gene_symbol = (string)bgGene.Clone();
                    new_library_line.Additional_information = "";
                    new_library_line.Organism = Options.Organism;
                    new_library_line.Organism_string = "";
                    new_library_line.Level = -1;
                    background_gene_lines.Add(new_library_line);
                }
                ontology_library.Add_to_array(background_gene_lines.ToArray());
                if (this.Experimental_bg_genes.Length>0) { ontology_library.Keep_only_indicated_genes(this.Experimental_bg_genes); }
                Ontology_libraries[indexO] = ontology_library;
            }
            write.Close();
        }

        public void Generate(params string[] bg_genes)
        {
            if (bg_genes.Length > 0)
            {
                this.Experimental_bg_genes = Array_class.Deep_copy_string_array(bg_genes);
            }
            else
            {
                // throw new Exception();
                this.Experimental_bg_genes = new string[0];
            }
            Generate_ontology_libraries();

            Fisher = new Fisher_exact_test_class(50000, false);
        }
        #endregion

        private Enrichment2018_results_class Do_enrichment_analysis_with_single_ontology_on_all_columns(Ontology_library_class ontology_library, DE_class de_input)
        {
            Ontology_type_enum current_ontology = ontology_library.Get_current_ontology_and_check_if_only_one();

            DE_class de = de_input.Deep_copy();
            string[] all_ontology_genes = ontology_library.Get_all_ordered_unique_gene_symbols();
            de.Keep_only_stated_symbols(all_ontology_genes);
            Dictionary<string, DE_line_class> deSymbol_deLine_dict = new Dictionary<string, DE_line_class>();
            DE_line_class de_line;
            int de_length = de.DE.Count;
            int de_colCount = de.ColChar.Columns.Count;
            for (int indexDe = 0; indexDe < de_length; indexDe++)
            {
                de_line = de.DE[indexDe];
                deSymbol_deLine_dict.Add(de_line.Gene_symbol, de_line);
            }

            List<Enrichment2018_results_line_class> enrichment_results_list = new List<Enrichment2018_results_line_class>();
            Dictionary<string, List<string>[]> pathway_deColumnsOverlapSymbols_dict = new Dictionary<string, List<string>[]>();
            Dictionary<string, int> pathway_pathwaySymbolsCount_dict = new Dictionary<string, int>();
            int library_length = ontology_library.Library.Length;
            Ontology_library_line_class ontology_library_line;
            DE_line_class current_line;
            string current_scp;
            string current_symbol;

            DE_line_class de_line_copy;
            DE_class de_copy = new DE_class();
            Dictionary<string, bool> gene_already_considered_dict = new Dictionary<string, bool>();
            Dictionary<string, bool> ontology_bgSymbols_dict = new Dictionary<string, bool>();

            for (int indexL = 0; indexL < library_length; indexL++)
            {
                ontology_library_line = ontology_library.Library[indexL];
                current_scp = ontology_library_line.Scp;
                current_symbol = ontology_library_line.Target_gene_symbol;
                if (!ontology_bgSymbols_dict.ContainsKey(current_symbol))
                {
                    ontology_bgSymbols_dict.Add(current_symbol, true);
                }
                if (!pathway_pathwaySymbolsCount_dict.ContainsKey(current_scp))
                {
                    pathway_pathwaySymbolsCount_dict.Add(current_scp, 0);
                }
                pathway_pathwaySymbolsCount_dict[current_scp]++;
                if (deSymbol_deLine_dict.ContainsKey(current_symbol))
                {
                    current_line = deSymbol_deLine_dict[current_symbol];
                    de_line_copy = current_line.Deep_copy();
                    de_copy.DE.Add(de_line_copy);
                    if (!pathway_deColumnsOverlapSymbols_dict.ContainsKey(ontology_library_line.Scp))
                    {
                        pathway_deColumnsOverlapSymbols_dict.Add(current_scp, new List<string>[de_colCount]);
                        for (int indexCol = 0; indexCol < de_colCount; indexCol++)
                        {
                            pathway_deColumnsOverlapSymbols_dict[current_scp][indexCol] = new List<string>();
                        }
                    }
                    for (int indexCol = 0; indexCol < de_colCount; indexCol++)
                    {
                        if (current_line.Columns[indexCol].Value != 0)
                        {
                            pathway_deColumnsOverlapSymbols_dict[current_scp][indexCol].Add(current_symbol);
                        }
                    }
                }
            }

            int[] experimental_symbols_in_columns_count = de.Get_non_zero_counts_of_each_column_in_indexOrder();
            string[] scps = pathway_deColumnsOverlapSymbols_dict.Keys.ToArray();
            int scps_length = scps.Length;
            string scp;
            List<string>[] currentScp_overlap_symbols_of_columns;
            List<string> currentScpColumn_overlap_symbols;
            Enrichment2018_results_line_class new_enrichment_results_line;
            List<Enrichment2018_results_line_class> new_enrichment_results_list = new List<Enrichment2018_results_line_class>();
            List<DE_column_characterization_line_class> colChar_columns = de.ColChar.Columns;
            string[] complete_sampleNames_of_columns = new string[de_colCount];
            int bg_genes_count = ontology_bgSymbols_dict.Keys.ToArray().Length;
            for (int indexCol = 0; indexCol < de_colCount; indexCol++)
            {
                complete_sampleNames_of_columns[indexCol] = de.ColChar.Columns[indexCol].Get_full_column_name();
            }
            int a; int b; int c; int d; double p_value;
            for (int indexScp = 0; indexScp < scps_length; indexScp++)
            {
                scp = scps[indexScp];
                currentScp_overlap_symbols_of_columns = pathway_deColumnsOverlapSymbols_dict[scp];
                for (int indexCol = 0; indexCol < de_colCount; indexCol++)
                {
                    currentScpColumn_overlap_symbols = currentScp_overlap_symbols_of_columns[indexCol];
                    if (currentScpColumn_overlap_symbols.Count > 0)
                    {
                        new_enrichment_results_line = new Enrichment2018_results_line_class();
                        new_enrichment_results_line.Scp = (string)scp.Clone();
                        new_enrichment_results_line.Overlap_symbols = currentScpColumn_overlap_symbols.Distinct().OrderBy(l => l).ToArray();
                        if (new_enrichment_results_line.Overlap_symbols.Length != currentScpColumn_overlap_symbols.Count) 
                        {
                            string comment = new_enrichment_results_line.Scp;
                            System.IO.StreamWriter report = new System.IO.StreamWriter(Global_directory_class.Get_report_fileName(), true);
                            report.WriteLine(comment);
                            comment = new_enrichment_results_line.ReadWrite_overlap_symbols + " (" + new_enrichment_results_line.Overlap_symbols.Length + " vs " + currentScpColumn_overlap_symbols.Count + ")";
                            report.WriteLine(comment);
                            comment = "indexSCP = " + indexScp + " (" + scp + "); indexCol = " + indexCol;
                            report.WriteLine(comment);
                            report.Close();
                            throw new Exception(); 
                        }
                        new_enrichment_results_line.Overlap_count = new_enrichment_results_line.Overlap_symbols.Length;
                        new_enrichment_results_line.Ontology = current_ontology;
                        new_enrichment_results_line.Experimental_genes_count = experimental_symbols_in_columns_count[indexCol];
                        new_enrichment_results_line.Scp_genes_count = pathway_pathwaySymbolsCount_dict[scp];
                        new_enrichment_results_line.Sample_name = (string)complete_sampleNames_of_columns[indexCol].Clone();
                        new_enrichment_results_line.Bg_genes_count = bg_genes_count;

                        a = new_enrichment_results_line.Overlap_count;
                        b = new_enrichment_results_line.Scp_genes_count - a;
                        c = new_enrichment_results_line.Experimental_genes_count - a;
                        d = new_enrichment_results_line.Bg_genes_count - a - b - c;

                        if ((a < 0) || (b < 0) || (c < 0) || (d < 0)) 
                        {
                            string comment = "a=" + a + ", b=" + b + ", c=" + c + ", d=" + d;
                            System.IO.StreamWriter report = new System.IO.StreamWriter(Global_directory_class.Get_report_fileName(), true);
                            report.WriteLine(comment);
                            report.Close();
                            throw new Exception(); 
                        }

                        p_value = Fisher.Get_rightTailed_p_value(a, b, c, d);
                        new_enrichment_results_line.Pvalue = p_value;
                        new_enrichment_results_line.Minus_log10_pvalue = (float)-Math.Log(p_value, 10);
                        new_enrichment_results_list.Add(new_enrichment_results_line);
                    }
                }
            }
            Enrichment2018_results_class enrichment_results = new Enrichment2018_results_class();
            enrichment_results.Enrichment_results = new_enrichment_results_list.ToArray();
            return enrichment_results;
        }

        private Enrichment2018_results_class[] Generate_filtered_enrichment_results(Enrichment2018_results_class[] enrichment_results_array)
        {
            int enrichment2018_results_length = enrichment_results_array.Length;
            Enrichment2018_results_class current_enrichment_results;
            Enrichment2018_results_class filtered_enrichment_results;
            Enrichment2018_results_class[] filtered_enrichment_results_array = new Enrichment2018_results_class[enrichment2018_results_length];
            for (int indexE = 0; indexE < enrichment2018_results_length; indexE++)
            {
                current_enrichment_results = enrichment_results_array[indexE];
                filtered_enrichment_results = current_enrichment_results.Deep_copy();
                filtered_enrichment_results.Keep_only_top_x_ranked_scps_per_condition(this.Options.Top_top_x_predictions);
                filtered_enrichment_results_array[indexE] = filtered_enrichment_results;
            }
            return filtered_enrichment_results_array;
        }

        public Enrichment2018_results_class[] Analyse_de_instance_and_return_unfiltered_enrichment_results(DE_class de_input, string subdirectory, string add_results_file_name)
        {
            DE_class de = de_input.Deep_copy();
            if (this.Experimental_bg_genes.Length > 0)
            {
                de.Keep_only_stated_symbols(this.Experimental_bg_genes); //Non ontology genes will be removed for each ontology separately
            }
            Ontology_library_class current_ontology_library;
            int ontologies_length = Ontology_libraries.Length;
            Enrichment2018_results_class new_enrichment_results;
            Enrichment2018_results_class[] enrichment_results = new Enrichment2018_results_class[ontologies_length];
            for (int indexO = 0; indexO < ontologies_length; indexO++)
            {
                current_ontology_library = Ontology_libraries[indexO];
                new_enrichment_results = Do_enrichment_analysis_with_single_ontology_on_all_columns(current_ontology_library, de);
                new_enrichment_results.Calculate_fractional_ranks_for_scps_based_on_selected_valuetype(Enrichment2018_value_type_enum.Minuslog10pvalue);
                new_enrichment_results.Order_by_ontology_sampleName_pvalue();
                new_enrichment_results.Remove_indicated_scps("Background genes");
                if (Options.Add_missing_scps_identified_in_other_conditions)
                {
                    new_enrichment_results.Add_missing_scps_that_were_detected_at_least_once_with_pvalue_one_and_indicated_rank(99999, de);
                }
                enrichment_results[indexO] = new_enrichment_results;
            }
            if (Options.Write_results)
            {
                 Write_enrichment_results_for_R(enrichment_results, subdirectory, add_results_file_name + "_all");
                 Enrichment2018_results_class[] filtered_enrichment_results = Generate_filtered_enrichment_results(enrichment_results);
                 Write_enrichment_results_for_R(filtered_enrichment_results, subdirectory, add_results_file_name + "_filtered");
            }
            return enrichment_results;
        }

        public void Write_enrichment_results_for_R(Enrichment2018_results_class[] enrichment_results, string subdirectory, string addition_at_end_of_file)
        {
            ReadWrite.ReadWriteClass.Create_directory_if_it_does_not_exist(Global_directory_class.Results_directory + subdirectory);
            subdirectory = subdirectory + "//";
            Ontology_type_enum ontology;
            string ontology_string = "error";
            foreach (Enrichment2018_results_class enrichment_result in enrichment_results)
            {
                if (enrichment_result.Enrichment_results.Length > 0)
                {
                    ontology = enrichment_result.Get_ontology_and_check_if_only_one();
                    if (Options.Use_ontology_abbreviations_for_directories_and_files)
                    {
                        ontology_string = Ontology_abbreviation_class.Get_abbreviation_of_ontology(ontology);
                    }
                    else
                    {
                        ontology_string = ontology.ToString();
                    }
                    enrichment_result.Write_for_r(subdirectory, ontology_string + addition_at_end_of_file + ".txt");
                }
            }
        }
    }


}