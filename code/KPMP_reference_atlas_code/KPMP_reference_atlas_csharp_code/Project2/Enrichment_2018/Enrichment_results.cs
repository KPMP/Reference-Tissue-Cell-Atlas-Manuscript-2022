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
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Highthroughput_data;
using Enumerations;
using ReadWrite;
using Statistic;


namespace Enrichment_2018
{
    enum Enrichment2018_color_specification_enum { E_m_p_t_y, Regular, Highlight };
    enum Enrichment2018_value_type_enum {  E_m_p_t_y, Minuslog10pvalue, Fractional_rank, MinusLog10fdr }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////

    class De_condition_line_class
    {
        public string SampleName { get; set; }

        public De_condition_line_class Deep_copy()
        {
            De_condition_line_class copy = (De_condition_line_class)this.MemberwiseClone();
            copy.SampleName = (string)this.SampleName.Clone();
            return copy;
        }
    }

    class Enrichment2018_results_line_class
    {
        public Ontology_type_enum  Ontology { get; set; }
        public string Scp { get; set; }
        public string Sample_name { get; set; }
        public double Pvalue { get; set; }
        public string Integration_group { get; set; }
        public double FDR { get; set; }

        public float Fractional_rank { get; set; }
        public float Ordinal_rank { get; set; }
        public float Minus_log10_pvalue { get; set; }
        public float Minus_log10_fdr { get; set; }
        public int Overlap_count { get; set; }
        public int Scp_genes_count { get; set; }
        public int Experimental_genes_count { get; set; }
        public Enrichment2018_color_specification_enum Color_specification { get; set; }
        public int Bg_genes_count { get; set; }
        public string[] Overlap_symbols { get; set; }

        public static bool Check_for_correct_ordering {  get { return Global_class.Check_ordering; } }

        public string ReadWrite_overlap_symbols
        {
            get { return ReadWriteClass.Get_writeLine_from_array(this.Overlap_symbols, Enrichment2018_results_readWriteOptions_class.Array_delimiter); }
            set { this.Overlap_symbols = ReadWriteClass.Get_array_from_readLine<string>(value, Enrichment2018_results_readWriteOptions_class.Array_delimiter); }
        }

        #region Order
        public static Enrichment2018_results_line_class[] Order_by_pvalue(Enrichment2018_results_line_class[] lines)
        {
            Dictionary<double, List<Enrichment2018_results_line_class>> pvalue_dict = new Dictionary<double, List<Enrichment2018_results_line_class>>();
            int lines_length = lines.Length;
            Enrichment2018_results_line_class line;
            for (int indexL = 0; indexL < lines_length; indexL++)
            {
                line = lines[indexL];
                if (!pvalue_dict.ContainsKey(line.Pvalue))
                {
                    pvalue_dict.Add(line.Pvalue, new List<Enrichment2018_results_line_class>());
                }
                pvalue_dict[line.Pvalue].Add(line);
            }

            double[] pvalues = pvalue_dict.Keys.ToArray();
            double pvalue;
            int pvalues_length = pvalues.Length;
            List<Enrichment2018_results_line_class> ordered_enrichment_results = new List<Enrichment2018_results_line_class>();
            pvalues = pvalues.OrderBy(l => l).ToArray();
            for (int indexP = 0; indexP < pvalues_length; indexP++)
            {
                pvalue = pvalues[indexP];
                ordered_enrichment_results.AddRange(pvalue_dict[pvalue]);
            }

            if (Check_for_correct_ordering)
            {
                #region Check for correct ordering

                int ordered_length = ordered_enrichment_results.Count;
                Enrichment2018_results_line_class previous_line;
                Enrichment2018_results_line_class current_line;
                for (int indexO = 1; indexO < ordered_length; indexO++)
                {
                    previous_line = ordered_enrichment_results[indexO - 1];
                    current_line = ordered_enrichment_results[indexO];
                    if (current_line.Pvalue.CompareTo(previous_line.Pvalue) < 0) { throw new Exception(); }
                }
                #endregion
            }
            return ordered_enrichment_results.ToArray();
        }

        public static Enrichment2018_results_line_class[] Order_by_ontology_sampleName_pvalue(Enrichment2018_results_line_class[] lines)
        {
            Dictionary<Ontology_type_enum, Dictionary<string, Dictionary<double, List<Enrichment2018_results_line_class>>>> ontology_sampleName_pvalue_dict = new Dictionary<Ontology_type_enum, Dictionary<string, Dictionary<double, List<Enrichment2018_results_line_class>>>>();
            Dictionary<string, Dictionary<double, List<Enrichment2018_results_line_class>>> sampleName_pvalue_dict = new Dictionary<string, Dictionary<double, List<Enrichment2018_results_line_class>>>();
            Dictionary<double, List<Enrichment2018_results_line_class>> pvalue_dict = new Dictionary<double, List<Enrichment2018_results_line_class>>();
            int lines_length = lines.Length;
            Enrichment2018_results_line_class line;
            for (int indexL = 0; indexL < lines_length; indexL++)
            {
                line = lines[indexL];
                if (!ontology_sampleName_pvalue_dict.ContainsKey(line.Ontology))
                {
                    ontology_sampleName_pvalue_dict.Add(line.Ontology, new Dictionary<string, Dictionary<double, List<Enrichment2018_results_line_class>>>());
                }
                if (!ontology_sampleName_pvalue_dict[line.Ontology].ContainsKey(line.Sample_name))
                {
                    ontology_sampleName_pvalue_dict[line.Ontology].Add(line.Sample_name, new Dictionary<double, List<Enrichment2018_results_line_class>>());
                }
                if (!ontology_sampleName_pvalue_dict[line.Ontology][line.Sample_name].ContainsKey(line.Pvalue))
                {
                    ontology_sampleName_pvalue_dict[line.Ontology][line.Sample_name].Add(line.Pvalue, new List<Enrichment2018_results_line_class>());
                }
                ontology_sampleName_pvalue_dict[line.Ontology][line.Sample_name][line.Pvalue].Add(line);
            }

            Ontology_type_enum[] ontologies = ontology_sampleName_pvalue_dict.Keys.ToArray();
            Ontology_type_enum ontology;
            int ontologies_length = ontologies.Length;
            string[] sampleNames;
            string sampleName;
            int sampleNames_length;
            double[] pvalues;
            double pvalue;
            int pvalues_length;

            ontologies = ontologies.OrderBy(l => l).ToArray();
            List<Enrichment2018_results_line_class> ordered_enrichment_results = new List<Enrichment2018_results_line_class>();
            for (int indexO = 0; indexO < ontologies_length; indexO++)
            {
                ontology = ontologies[indexO];
                sampleName_pvalue_dict = ontology_sampleName_pvalue_dict[ontology];
                sampleNames = sampleName_pvalue_dict.Keys.ToArray();
                sampleNames_length = sampleNames.Length;
                sampleNames = sampleNames.OrderBy(l => l).ToArray();
                for (int indexSN = 0; indexSN < sampleNames_length; indexSN++)
                {
                    sampleName = sampleNames[indexSN];
                    pvalue_dict = sampleName_pvalue_dict[sampleName];
                    pvalues = pvalue_dict.Keys.ToArray();
                    pvalues_length = pvalues.Length;
                    pvalues = pvalues.OrderBy(l => l).ToArray();
                    for (int indexP = 0; indexP < pvalues_length; indexP++)
                    {
                        pvalue = pvalues[indexP];
                        ordered_enrichment_results.AddRange(pvalue_dict[pvalue]);
                    }
                }
            }
            if (Check_for_correct_ordering)
            {
                #region Check for correct ordering
                int ordered_length = ordered_enrichment_results.Count;
                if (ordered_length != lines_length) { throw new Exception(); }
                Enrichment2018_results_line_class previous_line;
                Enrichment2018_results_line_class current_line;
                for (int indexOrder = 1; indexOrder < ordered_length; indexOrder++)
                {
                    previous_line = ordered_enrichment_results[indexOrder - 1];
                    current_line = ordered_enrichment_results[indexOrder];
                    if (current_line.Ontology.CompareTo(previous_line.Ontology) < 0) { throw new Exception(); }
                    if ((current_line.Ontology.Equals(previous_line.Ontology))
                        && (current_line.Sample_name.CompareTo(previous_line.Sample_name) < 0)) { throw new Exception(); }
                    if (   (current_line.Ontology.Equals(previous_line.Ontology))
                        && (current_line.Sample_name.Equals(previous_line.Sample_name))
                        && (current_line.Pvalue.CompareTo(previous_line.Pvalue) < 0)) { throw new Exception(); }
                }
                #endregion
            }
            return ordered_enrichment_results.ToArray();
        }

        public static Enrichment2018_results_line_class[] Order_by_ontology_sampleName_descendingMinusLog10Pvalue(Enrichment2018_results_line_class[] lines)
        {
            Dictionary<Ontology_type_enum, Dictionary<string, Dictionary<double, List<Enrichment2018_results_line_class>>>> ontology_sampleName_minusLog10Pvalue_dict = new Dictionary<Ontology_type_enum, Dictionary<string, Dictionary<double, List<Enrichment2018_results_line_class>>>>();
            Dictionary<string, Dictionary<double, List<Enrichment2018_results_line_class>>> sampleName_minusLog10Pvalue_dict = new Dictionary<string, Dictionary<double, List<Enrichment2018_results_line_class>>>();
            Dictionary<double, List<Enrichment2018_results_line_class>> minusLog10Pvalue_dict = new Dictionary<double, List<Enrichment2018_results_line_class>>();
            int lines_length = lines.Length;
            Enrichment2018_results_line_class line;
            for (int indexL = 0; indexL < lines_length; indexL++)
            {
                line = lines[indexL];
                if (!ontology_sampleName_minusLog10Pvalue_dict.ContainsKey(line.Ontology))
                {
                    ontology_sampleName_minusLog10Pvalue_dict.Add(line.Ontology, new Dictionary<string, Dictionary<double, List<Enrichment2018_results_line_class>>>());
                }
                if (!ontology_sampleName_minusLog10Pvalue_dict[line.Ontology].ContainsKey(line.Sample_name))
                {
                    ontology_sampleName_minusLog10Pvalue_dict[line.Ontology].Add(line.Sample_name, new Dictionary<double, List<Enrichment2018_results_line_class>>());
                }
                if (!ontology_sampleName_minusLog10Pvalue_dict[line.Ontology][line.Sample_name].ContainsKey(line.Minus_log10_pvalue))
                {
                    ontology_sampleName_minusLog10Pvalue_dict[line.Ontology][line.Sample_name].Add(line.Minus_log10_pvalue, new List<Enrichment2018_results_line_class>());
                }
                ontology_sampleName_minusLog10Pvalue_dict[line.Ontology][line.Sample_name][line.Minus_log10_pvalue].Add(line);
            }

            Ontology_type_enum[] ontologies = ontology_sampleName_minusLog10Pvalue_dict.Keys.ToArray();
            Ontology_type_enum ontology;
            int ontologies_length = ontologies.Length;
            string[] sampleNames;
            string sampleName;
            int sampleNames_length;
            double[] minusLog10_pvalues;
            double minusLog10_pvalue;
            int minusLog10_pvalues_length;

            ontologies = ontologies.OrderBy(l => l).ToArray();
            List<Enrichment2018_results_line_class> ordered_enrichment_results = new List<Enrichment2018_results_line_class>();
            for (int indexO = 0; indexO < ontologies_length; indexO++)
            {
                ontology = ontologies[indexO];
                sampleName_minusLog10Pvalue_dict = ontology_sampleName_minusLog10Pvalue_dict[ontology];
                sampleNames = sampleName_minusLog10Pvalue_dict.Keys.ToArray();
                sampleNames_length = sampleNames.Length;
                sampleNames = sampleNames.OrderBy(l => l).ToArray();
                for (int indexSN = 0; indexSN < sampleNames_length; indexSN++)
                {
                    sampleName = sampleNames[indexSN];
                    minusLog10Pvalue_dict = sampleName_minusLog10Pvalue_dict[sampleName];
                    minusLog10_pvalues = minusLog10Pvalue_dict.Keys.ToArray();
                    minusLog10_pvalues_length = minusLog10_pvalues.Length;
                    minusLog10_pvalues = minusLog10_pvalues.OrderByDescending(l => l).ToArray();
                    for (int indexMinus = 0; indexMinus < minusLog10_pvalues_length; indexMinus++)
                    {
                        minusLog10_pvalue = minusLog10_pvalues[indexMinus];
                        ordered_enrichment_results.AddRange(minusLog10Pvalue_dict[minusLog10_pvalue]);
                    }
                }
            }
            if (Check_for_correct_ordering)
            {
                #region Check for correct ordering
                int ordered_length = ordered_enrichment_results.Count;
                if (ordered_length != lines_length) { throw new Exception(); }
                Enrichment2018_results_line_class previous_line;
                Enrichment2018_results_line_class current_line;
                for (int indexOrder = 1; indexOrder < ordered_length; indexOrder++)
                {
                    previous_line = ordered_enrichment_results[indexOrder - 1];
                    current_line = ordered_enrichment_results[indexOrder];
                    if (current_line.Ontology.CompareTo(previous_line.Ontology) < 0) { throw new Exception(); }
                    if ((current_line.Ontology.Equals(previous_line.Ontology))
                        && (current_line.Sample_name.CompareTo(previous_line.Sample_name) < 0)) { throw new Exception(); }
                    if ((current_line.Ontology.Equals(previous_line.Ontology))
                        && (current_line.Sample_name.Equals(previous_line.Sample_name))
                        && (current_line.Minus_log10_pvalue.CompareTo(previous_line.Minus_log10_pvalue) > 0)) { throw new Exception(); }
                }
                #endregion
            }
            return ordered_enrichment_results.ToArray();
        }

        public static Enrichment2018_results_line_class[] Order_by_ontology_scpName_sampleName(Enrichment2018_results_line_class[] lines)
        {
            Dictionary<Ontology_type_enum, Dictionary<string, Dictionary<string, List<Enrichment2018_results_line_class>>>> ontology_scpName_sampleName_dict = new Dictionary<Ontology_type_enum, Dictionary<string, Dictionary<string, List<Enrichment2018_results_line_class>>>>();
            Dictionary<string, Dictionary<string, List<Enrichment2018_results_line_class>>> scpName_sampleName_dict = new Dictionary<string, Dictionary<string, List<Enrichment2018_results_line_class>>>();
            Dictionary<string, List<Enrichment2018_results_line_class>> sampleName_dict = new Dictionary<string, List<Enrichment2018_results_line_class>>();

            int lines_length = lines.Length;
            Enrichment2018_results_line_class line;
            for (int indexL = 0; indexL < lines_length; indexL++)
            {
                line = lines[indexL];
                if (!ontology_scpName_sampleName_dict.ContainsKey(line.Ontology))
                {
                    ontology_scpName_sampleName_dict.Add(line.Ontology, new Dictionary<string, Dictionary<string, List<Enrichment2018_results_line_class>>>());
                }
                if (!ontology_scpName_sampleName_dict[line.Ontology].ContainsKey(line.Scp))
                {
                    ontology_scpName_sampleName_dict[line.Ontology].Add(line.Scp, new Dictionary<string, List<Enrichment2018_results_line_class>>());
                }
                if (!ontology_scpName_sampleName_dict[line.Ontology][line.Scp].ContainsKey(line.Sample_name))
                {
                    ontology_scpName_sampleName_dict[line.Ontology][line.Scp].Add(line.Sample_name, new List<Enrichment2018_results_line_class>());
                }
                ontology_scpName_sampleName_dict[line.Ontology][line.Scp][line.Sample_name].Add(line);
            }

            Ontology_type_enum[] ontologies = ontology_scpName_sampleName_dict.Keys.ToArray();
            Ontology_type_enum ontology;
            int ontologies_length = ontologies.Length;
            string[] scpNames;
            string scpName;
            int scpNames_length;
            string[] sampleNames;
            string sampleName;
            int sampleNames_length;

            ontologies = ontologies.OrderBy(l => l).ToArray();
            List<Enrichment2018_results_line_class> ordered_enrichment_results = new List<Enrichment2018_results_line_class>();
            for (int indexO = 0; indexO < ontologies_length; indexO++)
            {
                ontology = ontologies[indexO];
                scpName_sampleName_dict = ontology_scpName_sampleName_dict[ontology];
                scpNames = scpName_sampleName_dict.Keys.ToArray();
                scpNames_length = scpNames.Length;
                scpNames = scpNames.OrderBy(l => l).ToArray();
                for (int indexSCP = 0; indexSCP < scpNames_length; indexSCP++)
                {
                    scpName = scpNames[indexSCP];
                    sampleName_dict = scpName_sampleName_dict[scpName];
                    sampleNames = sampleName_dict.Keys.ToArray();
                    sampleNames_length = sampleNames.Length;
                    sampleNames = sampleNames.OrderBy(l => l).ToArray();
                    for (int indexSN = 0; indexSN < sampleNames_length; indexSN++)
                    {
                        sampleName = sampleNames[indexSN];
                        ordered_enrichment_results.AddRange(sampleName_dict[sampleName]);
                    }
                }
            }

            if (Check_for_correct_ordering)
            {
                #region Check for correct ordering
                int ordered_length = ordered_enrichment_results.Count;
                if (ordered_length != lines_length) { throw new Exception(); }
                Enrichment2018_results_line_class previous_line;
                Enrichment2018_results_line_class current_line;
                for (int indexOrder = 1; indexOrder < ordered_length; indexOrder++)
                {
                    previous_line = ordered_enrichment_results[indexOrder - 1];
                    current_line = ordered_enrichment_results[indexOrder];
                    if (current_line.Ontology.CompareTo(previous_line.Ontology) < 0) { throw new Exception(); }
                    if ((current_line.Ontology.Equals(previous_line.Ontology))
                        && (current_line.Scp.CompareTo(previous_line.Scp) < 0)) { throw new Exception(); }
                    if ((current_line.Ontology.Equals(previous_line.Ontology))
                        && (current_line.Scp.Equals(previous_line.Scp))
                        && (current_line.Sample_name.CompareTo(previous_line.Sample_name) < 0)) { throw new Exception(); }
                }
            }
            #endregion

            return ordered_enrichment_results.ToArray();
        }
        #endregion

        public Enrichment2018_results_line_class()
        {
            this.Ordinal_rank = -1;
            this.Fractional_rank = -1;
            Color_specification = Enrichment2018_color_specification_enum.Regular;
            this.Integration_group = "";
        }

        public Enrichment2018_results_line_class Deep_copy()
        {
            Enrichment2018_results_line_class copy = (Enrichment2018_results_line_class)this.MemberwiseClone();
            copy.Scp = (string)this.Scp.Clone();
            copy.Sample_name = (string)this.Sample_name.Clone();
            int symbols_length = this.Overlap_symbols.Length;
            copy.Integration_group = (string)this.Integration_group.Clone();
            copy.Overlap_symbols = new string[symbols_length];
            for (int indexS=0; indexS<symbols_length; indexS++)
            {
                copy.Overlap_symbols[indexS] = (string)this.Overlap_symbols[indexS].Clone();
            }
            return copy;
        }
    }

    class Enrichment2018_results_readWriteOptions_class : ReadWriteOptions_base
    {
        public static char Array_delimiter { get { return ','; } }

        public Enrichment2018_results_readWriteOptions_class(string directory, string fileName)
        {
            ReadWriteClass.Create_directory_if_it_does_not_exist(directory);
            this.File = directory + fileName;
            this.Key_propertyNames = new string[] { "Ontology", "Integration_group", "Sample_name",                           "Scp",                       "Experimental_genes_count",                                     "Scp_genes_count",                                     "Bg_genes_count",             "Overlap_count",               "Pvalue", "Minus_log10_pvalue", "ReadWrite_overlap_symbols" };
            this.Key_columnNames = new string[]   { "Ontology", "Integration group", "Dataset and cell (sub)type or segment", "Subcellular process (SCP)", "Number of experimental genes that are part of background set", "Number of SCP genes that are part of background set", "Number of background genes", "Number of overlapping genes", "Pvalue", "Minus_log10_pvalue", "Overlapping gene symbols" };
            this.HeadlineDelimiters = new char[] { Global_class.Tab };
            this.LineDelimiters = new char[] { Global_class.Tab };
            this.File_has_headline = true;
            this.Report = ReadWrite_report_enum.Report_main;
        }
    }

    class Enrichment2018_results_mbcoResults_readOptions_class : ReadWriteOptions_base
    {
        public static char Array_delimiter { get { return ','; } }

        public Enrichment2018_results_mbcoResults_readOptions_class(string completeDirectory, string fileName)
        {
            string directory = completeDirectory;
            this.File = directory + fileName;
            this.Key_propertyNames = new string[] { "Sample_entryType", "Sample_name", "Scp", "Pvalue", "Minus_log10_pvalue", "Fractional_rank", "Overlap_count", "Scp_genes_count", "Experimental_genes_count", "ReadWrite_overlap_symbols" };
            this.Key_columnNames = new string[] { "EntryType", "Sample_name", "Scp_name", "Pvalue", "Minus_log10_pvalue", "Fractional_rank", "Overlap_count", "Process_symbols_count", "Experimental_symbols_count", "ReadWrite_overlap_symbols" };
            this.HeadlineDelimiters = new char[] { Global_class.Tab };
            this.LineDelimiters = new char[] { Global_class.Tab };
            this.File_has_headline = true;
            this.Report = ReadWrite_report_enum.Report_main;
        }
    }

    class Enrichment2018_results_readOptions_for_r_class : ReadWriteOptions_base
    {
        public Enrichment2018_results_readOptions_for_r_class(string directory, string fileName)
        {
            ReadWriteClass.Create_directory_if_it_does_not_exist(directory);
            this.File = directory + fileName;
            this.Key_propertyNames = new string[] { "Ontology", "Sample_name", "Scp", "Experimental_genes_count", "Scp_genes_count", "Bg_genes_count", "Overlap_count", "Pvalue", "Minus_log10_pvalue", "FDR", "Minus_log10_fdr", "Fractional_rank", "Overlap_count", "Scp_genes_count", "Experimental_genes_count", "ReadWrite_overlap_symbols", "Rank_with_roman_numbers" };
            this.Key_propertyNames = new string[] { "Ontology", "Sample_name", "Scp", "Experimental_genes_count", "Scp_genes_count", "Bg_genes_count", "Overlap_count", "Pvalue", "Minus_log10_pvalue", "FDR", "Minus_log10_fdr", "Fractional_rank", "Overlap_count", "Scp_genes_count", "Experimental_genes_count", "ReadWrite_overlap_symbols", "Rank_with_roman_numbers" };
            this.Key_columnNames = this.Key_propertyNames;
            this.HeadlineDelimiters = new char[] { Global_class.Tab };
            this.LineDelimiters = new char[] { Global_class.Tab };
            this.File_has_headline = true;
            this.Report = ReadWrite_report_enum.Report_main;
        }
    }

    class Enrichment2018_results_class
    {
        public Enrichment2018_results_line_class[] Enrichment_results { get; set; }

        #region Order
        public void Order_by_pvalue()
        {
            this.Enrichment_results = Enrichment2018_results_line_class.Order_by_pvalue(this.Enrichment_results);
        }

        public void Order_by_ontology_sampleName_pvalue()
        {
            this.Enrichment_results = Enrichment2018_results_line_class.Order_by_ontology_sampleName_pvalue(this.Enrichment_results);
        }

        public void Order_by_ontology_scpName_sampleName()
        {
            this.Enrichment_results = Enrichment2018_results_line_class.Order_by_ontology_scpName_sampleName(this.Enrichment_results);
        }
        #endregion

        public Enrichment2018_results_class()
        {
            this.Enrichment_results = new Enrichment2018_results_line_class[0];
        }

        public void Extract_integrationGroup_from_sampleName_and_modify_sampleName()
        {
            string[] splitStrings;
            string splitString;
            int splitStrings_length;
            StringBuilder sb_sample_name = new StringBuilder();
            string integration_group;
            bool allPatients_reached = false;
            string lol = "LMD Proteomics-Glom - EC - allPatients - Minus_log10_pvalue - Endothelial glomerula - allPatients";
            foreach (Enrichment2018_results_line_class enrich_line in this.Enrichment_results)
            {
                splitStrings = enrich_line.Sample_name.Split('-');
                splitStrings_length = splitStrings.Length;
                allPatients_reached = false;
                sb_sample_name.Clear();
                integration_group = "";
                for (int indexSplit=0; indexSplit<splitStrings_length;indexSplit++)
                {
                    splitString = splitStrings[indexSplit];
                    if (splitString=="allPatients") { allPatients_reached = true; }
                    if (!allPatients_reached)
                    {
                        if (sb_sample_name.Length>0) { sb_sample_name.AppendFormat("-"); }
                        sb_sample_name.AppendFormat(Text2_class.Remove_space_comma_semicolon_colon_underline_from_end_and_beginning_of_text(splitString));
                    }
                    else if (  (!splitString.Equals("allPatients"))
                             && (!splitString.Equals("Minus_log10_pvalue"))
                             && (!splitString.Equals("Minus_log10_pvalue_adjusted")))
                    {
                        if (String.IsNullOrEmpty(integration_group)) { integration_group = splitString; }
                        else { throw new Exception(); }
                    }
                }
                enrich_line.Integration_group = (string)integration_group.Clone();
                enrich_line.Sample_name = sb_sample_name.ToString();
            }
        }

        public void Collapse_bulk_transcriptomics_proteomics_on_segments()
        {
            int enrich_length = this.Enrichment_results.Length;
            Enrichment2018_results_line_class enrichment_line;
            List<Enrichment2018_results_line_class> keep = new List<Enrichment2018_results_line_class>();
            for (int indexE=0; indexE<enrich_length;indexE++)
            {
                enrichment_line = this.Enrichment_results[indexE];
                if (enrichment_line.Sample_name.Equals("LMD Proteomics-Glom-EC")) { }
                else if (enrichment_line.Sample_name.Equals("LMD Proteomics-Glom-EPC")) { }
                else if (enrichment_line.Sample_name.Equals("LMD Proteomics-Glom-Mesangial")) { }
                else if (enrichment_line.Sample_name.Equals("LMD Proteomics-Glom-POD"))
                {
                    enrichment_line.Sample_name = "LMD Proteomics-Glom";
                    keep.Add(enrichment_line);
                }
                else if (enrichment_line.Sample_name.Equals("NSC Proteomics-Glom-EC")) { }
                else if (enrichment_line.Sample_name.Equals("NSC Proteomics-Glom-EPC")) { }
                else if (enrichment_line.Sample_name.Equals("NSC Proteomics-Glom-Mesangial")) { }
                else if (enrichment_line.Sample_name.Equals("NSC Proteomics-Glom-POD"))
                {
                    enrichment_line.Sample_name = "NSC Proteomics-Glom";
                    keep.Add(enrichment_line);
                }
                else if (enrichment_line.Sample_name.Equals("LMD RNASeq-CD-IC")) { }
                else if (enrichment_line.Sample_name.Equals("LMD RNASeq-CD-PC"))
                {
                    enrichment_line.Sample_name = "LMD RNASeq-CD";
                    keep.Add(enrichment_line);
                }
                else if (enrichment_line.Sample_name.Equals("LMD RNASeq-Glom-EC")) { }
                else if (enrichment_line.Sample_name.Equals("LMD RNASeq-Glom-EPC")) { }
                else if (enrichment_line.Sample_name.Equals("LMD RNASeq-Glom-Mesangial")) { }
                else if (enrichment_line.Sample_name.Equals("LMD RNASeq-Glom-POD"))
                {
                    enrichment_line.Sample_name = "LMD RNASeq-CD";
                    keep.Add(enrichment_line);
                }
                else if (enrichment_line.Sample_name.Equals("LMD RNASeq-INT-FIB")) { }
                else if (enrichment_line.Sample_name.Equals("LMD RNASeq-INT-MAC"))
                {
                    enrichment_line.Sample_name = "LMD RNASeq-INT";
                    keep.Add(enrichment_line);
                }
                else
                {
                    keep.Add(enrichment_line);
                }
            }
            this.Enrichment_results = keep.ToArray();
        }

        private void Add_to_array(Enrichment2018_results_line_class[] add_enrichment_results)
        {
            int this_enrichment_results_length = this.Enrichment_results.Length;
            int add_enrichment_results_length = add_enrichment_results.Length;
            int new_enrichment_results_length = this_enrichment_results_length + add_enrichment_results_length;
            Enrichment2018_results_line_class[] new_enrichment_results = new Enrichment2018_results_line_class[new_enrichment_results_length];
            int indexNew = -1;
            for (int indexThis = 0; indexThis < this_enrichment_results_length; indexThis++)
            {
                indexNew++;
                new_enrichment_results[indexNew] = this.Enrichment_results[indexThis];
            }
            for (int indexAdd = 0; indexAdd < add_enrichment_results_length; indexAdd++)
            {
                indexNew++;
                new_enrichment_results[indexNew] = add_enrichment_results[indexAdd];
            }
            this.Enrichment_results = new_enrichment_results;
        }


        public Ontology_type_enum Get_ontology_and_check_if_only_one()
        {
            Ontology_type_enum ontology = this.Enrichment_results[0].Ontology;
            foreach (Enrichment2018_results_line_class enrichment_results_line in Enrichment_results)
            {
                if (!ontology.Equals(enrichment_results_line.Ontology))
                {
                    throw new Exception();
                }
            }
            return ontology;
        }

        public void Calculate_ordinal_ranks_for_scps_based_on_pvalue()
        {
            int enrichment_results_length = Enrichment_results.Length;
            Enrichment2018_results_line_class results_line;
            this.Enrichment_results = Enrichment2018_results_line_class.Order_by_ontology_sampleName_pvalue(this.Enrichment_results);
            int ordinal_rank = 0;
            for (int indexE = 0; indexE < enrichment_results_length; indexE++)
            {
                results_line = this.Enrichment_results[indexE];
                if ((indexE == 0)
                    || (!results_line.Ontology.Equals(this.Enrichment_results[indexE - 1].Ontology))
                    || (!results_line.Sample_name.Equals(this.Enrichment_results[indexE - 1].Sample_name)))
                {
                    ordinal_rank = 0;
                }
                ordinal_rank++;
                results_line.Ordinal_rank = ordinal_rank;
            }
        }

        public void Calculate_fractional_ranks_for_scps_based_on_selected_valuetype(Enrichment2018_value_type_enum selected_valueType)
        {
            int enrichment_results_length = Enrichment_results.Length;
            Enrichment2018_results_line_class results_line;
            Enrichment2018_results_line_class inner_results_line;
            int firstIndex_sameCondition = -1;
            switch (selected_valueType)
            {
                case Enrichment2018_value_type_enum.Minuslog10pvalue:
                    this.Enrichment_results = Enrichment2018_results_line_class.Order_by_ontology_sampleName_descendingMinusLog10Pvalue(this.Enrichment_results);
                    break;
                default:
                    throw new Exception();
            }
            //this.Enrichment_results = this.Enrichment_results.OrderBy(l => l.Ontology).ThenBy(l => l.Sample_timepoint).ThenBy(l => l.Sample_entryType).ThenBy(l => l.Sample_name).ThenByDescending(l => l.Minus_log10_pvalue).ToArray();
            int ordinal_rank = 0;
            List<float> current_ordinal_ranks = new List<float>();
            float fractional_rank;
            bool do_calculate_rank = false;
            for (int indexE = 0; indexE < enrichment_results_length; indexE++)
            {
                results_line = this.Enrichment_results[indexE];
                if ((indexE == 0)
                    || (!results_line.Ontology.Equals(this.Enrichment_results[indexE - 1].Ontology))
                    || (!results_line.Sample_name.Equals(this.Enrichment_results[indexE - 1].Sample_name)))
                {
                    ordinal_rank = 0;
                }
                switch (selected_valueType)
                {
                    case Enrichment2018_value_type_enum.Minuslog10pvalue:
                        if ((indexE == 0)
                            || (!results_line.Ontology.Equals(this.Enrichment_results[indexE - 1].Ontology))
                            || (!results_line.Sample_name.Equals(this.Enrichment_results[indexE - 1].Sample_name))
                            || (!results_line.Minus_log10_pvalue.Equals(this.Enrichment_results[indexE - 1].Minus_log10_pvalue)))
                        {
                            current_ordinal_ranks.Clear();
                            firstIndex_sameCondition = indexE;
                        }
                        break;
                    case Enrichment2018_value_type_enum.MinusLog10fdr:
                        if ((indexE == 0)
                            || (!results_line.Ontology.Equals(this.Enrichment_results[indexE - 1].Ontology))
                            || (!results_line.Sample_name.Equals(this.Enrichment_results[indexE - 1].Sample_name))
                            || (!results_line.Minus_log10_fdr.Equals(this.Enrichment_results[indexE - 1].Minus_log10_fdr)))
                        {
                            current_ordinal_ranks.Clear();
                            firstIndex_sameCondition = indexE;
                        }
                        break;
                    default:
                        throw new Exception();
                }
                ordinal_rank++;
                current_ordinal_ranks.Add(ordinal_rank);
                do_calculate_rank = false;
                switch (selected_valueType)
                {
                    case Enrichment2018_value_type_enum.Minuslog10pvalue:
                        if ((indexE == enrichment_results_length - 1)
                            || (!results_line.Ontology.Equals(this.Enrichment_results[indexE + 1].Ontology))
                            || (!results_line.Sample_name.Equals(this.Enrichment_results[indexE + 1].Sample_name))
                            || (!results_line.Minus_log10_pvalue.Equals(this.Enrichment_results[indexE + 1].Minus_log10_pvalue)))
                        {
                            do_calculate_rank = true;
                        }
                        break;
                    case Enrichment2018_value_type_enum.MinusLog10fdr:
                        if ((indexE == enrichment_results_length - 1)
                            || (!results_line.Ontology.Equals(this.Enrichment_results[indexE + 1].Ontology))
                            || (!results_line.Sample_name.Equals(this.Enrichment_results[indexE + 1].Sample_name))
                            || (!results_line.Minus_log10_fdr.Equals(this.Enrichment_results[indexE + 1].Minus_log10_fdr)))
                        {
                            do_calculate_rank = true;
                        }
                        break;
                    default:
                        throw new Exception();
                }
                if (do_calculate_rank)
                { 
                    if (current_ordinal_ranks.Count > 1)
                    {
                        fractional_rank = Math_class.Get_average(current_ordinal_ranks.ToArray());
                        for (int indexInner = firstIndex_sameCondition; indexInner<=indexE;indexInner++)
                        {
                            inner_results_line = this.Enrichment_results[indexInner];
                            inner_results_line.Fractional_rank = fractional_rank;
                        }
                    }
                    else if (current_ordinal_ranks.Count==1)
                    {
                        if (firstIndex_sameCondition!=indexE) { throw new Exception(); }
                        fractional_rank = current_ordinal_ranks[0];
                        results_line.Fractional_rank = fractional_rank;
                    }
                    else { throw new Exception(); }
                }
            }
        }

        public void Add_missing_scps_that_were_detected_at_least_once_with_pvalue_one_and_indicated_rank(float rank, DE_class de)
        {
            Dictionary<string, bool> sampleName_dict = new Dictionary<string, bool>();
            int enrichment_results_length = this.Enrichment_results.Length;
            Enrichment2018_results_line_class enrichment_results_line;
            Enrichment2018_results_line_class add_enrichment_results_line;
            List<Enrichment2018_results_line_class> add_enrichment_result_lines = new List<Enrichment2018_results_line_class>();
            DE_column_characterization_class colChar = de.ColChar;
            DE_column_characterization_line_class colChar_line;
            int columns_length = colChar.Columns.Count;
            for (int indexC = 0; indexC < columns_length; indexC++)
            {
                colChar_line = colChar.Columns[indexC];
                if (!sampleName_dict.ContainsKey(colChar_line.Combined_names))
                { sampleName_dict.Add(colChar_line.Combined_names, true); }
            }

            for (int indexE = 0; indexE < enrichment_results_length; indexE++)
            {
                enrichment_results_line = this.Enrichment_results[indexE];
                if (!sampleName_dict.ContainsKey(enrichment_results_line.Sample_name))
                { sampleName_dict.Add(enrichment_results_line.Sample_name, true); }
            }

            string[] sampleNames;
            string sampleName;
            int sampleNames_length;
            for (int indexE = 0; indexE < enrichment_results_length; indexE++)
            {
                enrichment_results_line = this.Enrichment_results[indexE];
                if ((indexE == 0)
                    || (!enrichment_results_line.Scp.Equals(this.Enrichment_results[indexE - 1].Scp)))
                {
                    sampleNames = sampleName_dict.Keys.ToArray();
                    sampleNames_length = sampleNames.Length;
                    for (int indexSN = 0; indexSN < sampleNames_length; indexSN++)
                    {
                        sampleName = sampleNames[indexSN];
                        sampleName_dict[sampleName] = false;
                    }
                }
                if (sampleName_dict[enrichment_results_line.Sample_name] == true) { throw new Exception(); }
                sampleName_dict[enrichment_results_line.Sample_name] = true;
                if ((indexE == enrichment_results_length - 1)
                    || (!enrichment_results_line.Scp.Equals(this.Enrichment_results[indexE + 1].Scp)))
                {
                    sampleNames = sampleName_dict.Keys.ToArray();
                    sampleNames_length = sampleNames.Length;
                    for (int indexSN = 0; indexSN < sampleNames_length; indexSN++)
                    {
                        sampleName = sampleNames[indexSN];
                        if (sampleName_dict[sampleName] == false)
                        {
                            add_enrichment_results_line = new Enrichment2018_results_line_class();
                            add_enrichment_results_line.Ontology = enrichment_results_line.Ontology;
                            add_enrichment_results_line.Sample_name = (string)sampleName.Clone();
                            add_enrichment_results_line.Scp = (string)enrichment_results_line.Scp.Clone();
                            add_enrichment_results_line.Overlap_count = 0;
                            add_enrichment_results_line.Scp_genes_count = -1;
                            add_enrichment_results_line.Pvalue = 1;
                            add_enrichment_results_line.Minus_log10_pvalue = 0;
                            add_enrichment_results_line.FDR = 1;
                            add_enrichment_results_line.Minus_log10_fdr = 0;
                            add_enrichment_results_line.Ordinal_rank = rank;
                            add_enrichment_results_line.Fractional_rank = rank;
                            add_enrichment_results_line.Color_specification = Enrichment2018_color_specification_enum.Regular;
                            add_enrichment_results_line.Overlap_symbols = new string[0];
                            add_enrichment_result_lines.Add(add_enrichment_results_line);
                        }
                    }
                }
            }
            Add_to_array(add_enrichment_result_lines.ToArray());
        }

        #region Keep

        #region Keep/Remove SCPs in dependence of other SCPs
        private Enrichment2018_results_line_class[] Add_first_sibling_as_default_sibling_if_all_siblings_are_missing_in_sameGroup_enrichment_lines(Enrichment2018_results_line_class[] sameGroup_enrichment_lines, Dictionary<string, string[]> addDefaultFirstSibling_if_no_siblings_dict)
        {
            string sampleName = sameGroup_enrichment_lines[0].Sample_name;
            Dictionary<string, bool> scps_inCurrentSample_dict = new Dictionary<string, bool>();
            int enrichment_length = sameGroup_enrichment_lines.Length;
            Enrichment2018_results_line_class enrichment_line;
            for (int indexE = 0; indexE < enrichment_length; indexE++)
            {
                enrichment_line = sameGroup_enrichment_lines[indexE];
                if (!enrichment_line.Sample_name.Equals(sampleName)) { throw new Exception(); }
                if ((!scps_inCurrentSample_dict.ContainsKey(enrichment_line.Scp))
                    && (  (enrichment_line.Minus_log10_pvalue > 0)
                        ||(enrichment_line.Minus_log10_pvalue==-10)))
                { scps_inCurrentSample_dict.Add(enrichment_line.Scp, true); }
            }
            string[] sibling_scps;
            bool add_first_sibling_as_default_sibling;
            Enrichment2018_results_line_class default_sibling_line;
            List<Enrichment2018_results_line_class> new_lines = new List<Enrichment2018_results_line_class>();
            for (int indexE = 0; indexE < enrichment_length; indexE++)
            {
                enrichment_line = sameGroup_enrichment_lines[indexE];
                if (  (addDefaultFirstSibling_if_no_siblings_dict.ContainsKey(enrichment_line.Scp))
                    &&((enrichment_line.Minus_log10_pvalue>0) || (enrichment_line.Minus_log10_pvalue == -10)))
                {
                    sibling_scps = addDefaultFirstSibling_if_no_siblings_dict[enrichment_line.Scp];
                    add_first_sibling_as_default_sibling = true;
                    foreach (string sibling_scp in sibling_scps)
                    {
                        if (scps_inCurrentSample_dict.ContainsKey(sibling_scp)) { add_first_sibling_as_default_sibling = false; }
                    }
                    if (add_first_sibling_as_default_sibling)
                    {
                        default_sibling_line = new Enrichment2018_results_line_class();
                        default_sibling_line.Experimental_genes_count = enrichment_line.Experimental_genes_count;
                        default_sibling_line.Bg_genes_count = enrichment_line.Bg_genes_count;
                        default_sibling_line.FDR = 1;
                        default_sibling_line.Fractional_rank = 999999;
                        //default_sibling_line.Group = (string)enrichment_line.Group.Clone();
                        default_sibling_line.Minus_log10_fdr = 0;
                        default_sibling_line.Minus_log10_pvalue = -10;
                        default_sibling_line.Ontology = enrichment_line.Ontology;
                        default_sibling_line.Overlap_count = 0;
                        default_sibling_line.Overlap_symbols = new string[0];
                        default_sibling_line.Pvalue = 1;
                        default_sibling_line.Sample_name = (string)enrichment_line.Sample_name.Clone();
                        default_sibling_line.Scp = (string)sibling_scps[0].Clone();
                        default_sibling_line.Scp_genes_count = -1;
                        new_lines.Add(default_sibling_line);
                    }
                }
            }
            new_lines.AddRange(sameGroup_enrichment_lines);
            return new_lines.ToArray();
        }

        public void Add_first_sibling_as_default_sibling_if_all_siblings_are_missing(Dictionary<string, string[]> addDefaultFirstSibling_if_no_siblings_dict)
        {
            this.Enrichment_results = this.Enrichment_results.OrderBy(l => l.Sample_name).ToArray();
            Enrichment2018_results_line_class enrichment_line;
            Enrichment2018_results_line_class[] add_to_keep;
            List<Enrichment2018_results_line_class> currentGroup_lines = new List<Enrichment2018_results_line_class>();
            List<Enrichment2018_results_line_class> keep = new List<Enrichment2018_results_line_class>();
            int enrichment_length = this.Enrichment_results.Length;
            for (int indexE = 0; indexE < enrichment_length; indexE++)
            {
                enrichment_line = this.Enrichment_results[indexE];
                if ((indexE == 0)
                    || (!enrichment_line.Sample_name.Equals(this.Enrichment_results[indexE - 1].Sample_name)))
                {
                    currentGroup_lines.Clear(); ;
                }
                currentGroup_lines.Add(enrichment_line);
                if ((indexE == enrichment_length - 1)
                    || (!enrichment_line.Sample_name.Equals(this.Enrichment_results[indexE + 1].Sample_name)))
                {
                    add_to_keep = Add_first_sibling_as_default_sibling_if_all_siblings_are_missing_in_sameGroup_enrichment_lines(currentGroup_lines.ToArray(), addDefaultFirstSibling_if_no_siblings_dict);
                    keep.AddRange(add_to_keep);
                }
            }
            this.Enrichment_results = keep.ToArray();
        }

        private Enrichment2018_results_line_class[] Keep_scp_only_if_all_indicated_requiredScps_have_non_zero_minusLog10value_in_sameGroup_enrichment_lines(Enrichment2018_results_line_class[] sameGroup_enrichment_lines, Dictionary<string, string[]> parentScp_childScpsRequirement_dict)
        {
            string sampleName = sameGroup_enrichment_lines[0].Sample_name;
            Dictionary<string, bool> scps_inCurrentSample_dict = new Dictionary<string, bool>();
            int enrichment_length = sameGroup_enrichment_lines.Length;
            Enrichment2018_results_line_class enrichment_line;
            for (int indexE = 0; indexE < enrichment_length; indexE++)
            {
                enrichment_line = sameGroup_enrichment_lines[indexE];
                if (!enrichment_line.Sample_name.Equals(sampleName)) { throw new Exception(); }
                if ((!scps_inCurrentSample_dict.ContainsKey(enrichment_line.Scp))
                    && (  (enrichment_line.Minus_log10_pvalue > 0)
                        || (enrichment_line.Minus_log10_pvalue == -10))) //default pathways have -10 log10pvalue
                { scps_inCurrentSample_dict.Add(enrichment_line.Scp, true); }
            }
            string[] requirement_scps;
            bool keep_line;
            List<Enrichment2018_results_line_class> keep = new List<Enrichment2018_results_line_class>();
            for (int indexE = 0; indexE < enrichment_length; indexE++)
            {
                enrichment_line = sameGroup_enrichment_lines[indexE];
                keep_line = true;
                if (parentScp_childScpsRequirement_dict.ContainsKey(enrichment_line.Scp))
                {
                    requirement_scps = parentScp_childScpsRequirement_dict[enrichment_line.Scp];
                    foreach (string requirement_scp in requirement_scps)
                    {
                        if (!scps_inCurrentSample_dict.ContainsKey(requirement_scp)) { keep_line = false; }
                    }
                }
                if (keep_line) { keep.Add(enrichment_line); }
            }
            return keep.ToArray();
        }

        public void Keep_scp_only_if_all_indicated_requiredScps_have_non_zero_minusLog10Pvalue(Dictionary<string, string[]> parentScp_childScpsRequirement_dict)
        {
            this.Enrichment_results = this.Enrichment_results.OrderBy(l => l.Sample_name).ToArray();
            Enrichment2018_results_line_class enrichment_line;
            Enrichment2018_results_line_class[] add_to_keep;
            List<Enrichment2018_results_line_class> currentGroup_lines = new List<Enrichment2018_results_line_class>();
            List<Enrichment2018_results_line_class> keep = new List<Enrichment2018_results_line_class>();
            int enrichment_length = this.Enrichment_results.Length;
            for (int indexE = 0; indexE < enrichment_length; indexE++)
            {
                enrichment_line = this.Enrichment_results[indexE];
                if ((indexE == 0)
                    || (!enrichment_line.Sample_name.Equals(this.Enrichment_results[indexE - 1].Sample_name)))
                {
                    currentGroup_lines.Clear(); ;
                }
                currentGroup_lines.Add(enrichment_line);
                if ((indexE == enrichment_length - 1)
                    || (!enrichment_line.Sample_name.Equals(this.Enrichment_results[indexE + 1].Sample_name)))
                {
                    add_to_keep = Keep_scp_only_if_all_indicated_requiredScps_have_non_zero_minusLog10value_in_sameGroup_enrichment_lines(currentGroup_lines.ToArray(), parentScp_childScpsRequirement_dict);
                    keep.AddRange(add_to_keep);
                }
            }
            this.Enrichment_results = keep.ToArray();
        }

        private Enrichment2018_results_line_class[] Remove_scp_only_if_at_least_one_indicated_exclusionScp_has_non_zero_minusLog10value_in_sameGroup_enrichment_lines(Enrichment2018_results_line_class[] sameGroup_enrichment_lines, Dictionary<string, string[]> parentScp_childScpsExclusion_dict)
        {
            string sampleName = sameGroup_enrichment_lines[0].Sample_name;
            Dictionary<string, bool> scps_inCurrentSample_dict = new Dictionary<string, bool>();
            int enrichment_length = sameGroup_enrichment_lines.Length;
            Enrichment2018_results_line_class enrichment_line;
            for (int indexE = 0; indexE < enrichment_length; indexE++)
            {
                enrichment_line = sameGroup_enrichment_lines[indexE];
                if (!enrichment_line.Sample_name.Equals(sampleName)) { throw new Exception(); }
                if ((!scps_inCurrentSample_dict.ContainsKey(enrichment_line.Scp))
                    && (  (enrichment_line.Minus_log10_pvalue > 0)
                        ||(enrichment_line.Minus_log10_pvalue == -10)))
                { scps_inCurrentSample_dict.Add(enrichment_line.Scp, true); }
            }
            string[] exclusion_scps;
            bool keep_line;
            List<Enrichment2018_results_line_class> keep = new List<Enrichment2018_results_line_class>();
            for (int indexE = 0; indexE < enrichment_length; indexE++)
            {
                enrichment_line = sameGroup_enrichment_lines[indexE];
                keep_line = true;
                if (parentScp_childScpsExclusion_dict.ContainsKey(enrichment_line.Scp))
                {
                    exclusion_scps = parentScp_childScpsExclusion_dict[enrichment_line.Scp];
                    foreach (string exclusion_scp in exclusion_scps)
                    {
                        if (scps_inCurrentSample_dict.ContainsKey(exclusion_scp)) 
                        {
                            keep_line = false; 
                        }
                    }
                }
                if (keep_line) { keep.Add(enrichment_line); }
            }
            return keep.ToArray();
        }

        public void Remove_scp_only_if_at_least_one_indicated_exclusionScp_has_non_zero_minusLog10value(Dictionary<string, string[]> parentScp_childScpsExclusion_dict)
        {
            this.Enrichment_results = this.Enrichment_results.OrderBy(l => l.Sample_name).ToArray();
            Enrichment2018_results_line_class enrichment_line;
            Enrichment2018_results_line_class[] add_to_keep;
            List<Enrichment2018_results_line_class> currentGroup_lines = new List<Enrichment2018_results_line_class>();
            List<Enrichment2018_results_line_class> keep = new List<Enrichment2018_results_line_class>();
            int enrichment_length = this.Enrichment_results.Length;
            for (int indexE = 0; indexE < enrichment_length; indexE++)
            {
                enrichment_line = this.Enrichment_results[indexE];
                if ((indexE == 0)
                    || (!enrichment_line.Sample_name.Equals(this.Enrichment_results[indexE - 1].Sample_name)))
                {
                    currentGroup_lines.Clear(); ;
                }
                currentGroup_lines.Add(enrichment_line);
                if ((indexE == enrichment_length - 1)
                    || (!enrichment_line.Sample_name.Equals(this.Enrichment_results[indexE + 1].Sample_name)))
                {
                    add_to_keep = Remove_scp_only_if_at_least_one_indicated_exclusionScp_has_non_zero_minusLog10value_in_sameGroup_enrichment_lines(currentGroup_lines.ToArray(), parentScp_childScpsExclusion_dict);
                    keep.AddRange(add_to_keep);
                }
            }
            this.Enrichment_results = keep.ToArray();
        }
        #endregion

        public void Remove_indicated_scps(params string[] remove_scps)
        {
            remove_scps = remove_scps.Distinct().OrderBy(l => l).ToArray();
            string remove_scp;
            int remove_scps_length = remove_scps.Length;
            int indexRemove = 0;
            int stringCompare = -2;

            this.Enrichment_results = this.Enrichment_results.OrderBy(l => l.Scp).ToArray();

            int enrichment_length = Enrichment_results.Length;
            Enrichment2018_results_line_class enrichment_results_line;
            List<Enrichment2018_results_line_class> keep = new List<Enrichment2018_results_line_class>();
            for (int indexE = 0; indexE < enrichment_length; indexE++)
            {
                enrichment_results_line = this.Enrichment_results[indexE];
                stringCompare = -2;
                while ((indexRemove < remove_scps_length) && (stringCompare < 0))
                {
                    remove_scp = remove_scps[indexRemove];
                    stringCompare = remove_scp.CompareTo(enrichment_results_line.Scp);
                    if (stringCompare < 0)
                    {
                        indexRemove++;
                    }
                }
                if (stringCompare != 0)
                {
                    keep.Add(enrichment_results_line);
                }
            }
            this.Enrichment_results = keep.ToArray();
            //Check_if_not_empty();
        }

        public void Keep_only_top_x_ranked_scps_per_condition(int top_x)
        {
            this.Order_by_ontology_sampleName_pvalue();
            int results_length = this.Enrichment_results.Length;
            Enrichment2018_results_line_class results_line;
            List<Enrichment2018_results_line_class> keep_enrichment = new List<Enrichment2018_results_line_class>();
            int kept_per_condition = 0;
            for (int indexE = 0; indexE < results_length; indexE++)
            {
                results_line = this.Enrichment_results[indexE];
                if ((indexE == 0)
                    || (!results_line.Ontology.Equals(this.Enrichment_results[indexE - 1].Ontology))
                    || (!results_line.Sample_name.Equals(this.Enrichment_results[indexE - 1].Sample_name)))
                {
                    kept_per_condition = 0;
                }
                if (kept_per_condition < top_x)
                {
                    kept_per_condition++;
                    keep_enrichment.Add(results_line);
                }
            }
            this.Enrichment_results = keep_enrichment.ToArray();
        }
        #endregion

        #region Add other
        public void Add_other(Enrichment2018_results_class add_enrich)
        {
            this.Add_to_array(add_enrich.Enrichment_results);
        }
        #endregion

        #region Copy read write
        public Enrichment2018_results_class Deep_copy()
        {
            Enrichment2018_results_class copy = (Enrichment2018_results_class)this.MemberwiseClone();
            int enrichment_length = this.Enrichment_results.Length;
            copy.Enrichment_results = new Enrichment2018_results_line_class[enrichment_length];
            for (int indexE=0; indexE<enrichment_length; indexE++)
            {
                copy.Enrichment_results[indexE] = this.Enrichment_results[indexE].Deep_copy();
            }
            return copy;
        }
        public void Write(string subdirectory, string fileName)
        {
            Enrichment2018_results_readWriteOptions_class readWriteOptions = new Enrichment2018_results_readWriteOptions_class(subdirectory, fileName);
            ReadWriteClass.WriteData(Enrichment_results, readWriteOptions);
        }
        public void Write_for_r(string subdirectory, string fileName)
        {
            Enrichment2018_results_readOptions_for_r_class readWriteOptions = new Enrichment2018_results_readOptions_for_r_class(subdirectory, fileName);
            ReadWriteClass.WriteData(Enrichment_results, readWriteOptions);
        }
        #endregion
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
  
}
