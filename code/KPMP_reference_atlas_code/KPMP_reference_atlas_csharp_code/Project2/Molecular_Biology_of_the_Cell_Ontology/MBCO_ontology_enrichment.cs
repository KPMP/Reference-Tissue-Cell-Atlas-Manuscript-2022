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
using Enumerations;
using ReadWrite;
using Statistic;

namespace Mbco_ontology_enrichment
{
    class Ontology_enrichment_line_class
    {
        #region Fields for MBCO
        public string Scp_id { get; set; }
        public string Scp_name { get; set; }
        public string Parent_scp_name { get; set; }
        public Ontology_type_enum Ontology_type { get; set; }
        public int ProcessLevel { get; set; }
        #endregion

        #region Fields enrichment results
        public float Relative_visitation_frequency { get; set; }
        public float Minus_log10_pvalue { get; set; }
        public double Pvalue { get; set; }
        public double Qvalue { get; set; }
        public double FDR { get; set; }
        public float Fractional_rank { get; set; }
        public int Overlap_count { get; set; }
        public int Process_symbols_count { get; set; }
        public int Experimental_symbols_count { get; set; }
        public int Bg_symbol_count { get; set; }
        public string[] Overlap_symbols { get; set; }
        public string[] Overlap_symbols_with_other_conditions { get; set; }
        public bool Protected_against_greedy_symbol_removal { get; set; }
        public string ReadWrite_overlap_symbols
        {
            get { return ReadWriteClass.Get_writeLine_from_array(Overlap_symbols, Ontology_enrich_readWriteOptions_class.Delimiter); }
            set { Overlap_symbols = ReadWriteClass.Get_array_from_readLine<string>(value, Ontology_enrich_readWriteOptions_class.Delimiter); }
        }
        public string ReadWrite_overlap_symbols_with_other_conditions
        {
            get { return ReadWriteClass.Get_writeLine_from_array(Overlap_symbols_with_other_conditions, Ontology_enrich_readWriteOptions_class.Delimiter); }
            set { Overlap_symbols_with_other_conditions = ReadWriteClass.Get_array_from_readLine<string>(value, Ontology_enrich_readWriteOptions_class.Delimiter); }
        }
        #endregion

        #region Fields for sample
        public string Sample_name { get; set; }
        public string Complete_sample_name { get { return Sample_name; } }
        public string Integration_group { get; set; }
        public string Set_within_integration_group { get; set; }
        #endregion

        public Ontology_enrichment_line_class()
        {
            Scp_id = Global_class.Empty_entry;
            Scp_name = Global_class.Empty_entry;
            Parent_scp_name = Global_class.Empty_entry;
            Sample_name = Global_class.Empty_entry;
            Integration_group = Global_class.Empty_entry;
            Overlap_symbols = new string[0];
            Pvalue = -1;
        }

        public bool Is_mbco_ontology()
        {
            switch (this.Ontology_type)
            {
                case Ontology_type_enum.Mbco_level3:
                    return true;
                default:
                    return false;
            }
        }

        #region Standard way
        public static Ontology_enrichment_line_class[] Order_by_sample_and_scpName(Ontology_enrichment_line_class[] onto_enrich_array)
        {
            onto_enrich_array = onto_enrich_array.OrderBy(l => l.Ontology_type).ThenBy(l => l.Sample_name).ThenBy(l => l.Scp_name).ToArray();
            return onto_enrich_array;
        }

        public static Ontology_enrichment_line_class[] Order_by_sampleName_and_scpName(Ontology_enrichment_line_class[] onto_enrich_array)
        {
            onto_enrich_array = onto_enrich_array.OrderBy(l => l.Sample_name).ThenBy(l => l.Scp_name).ToArray();
            return onto_enrich_array;
        }

        public static Ontology_enrichment_line_class[] Order_by_complete_sample_pvalue(Ontology_enrichment_line_class[] onto_enrich_array)
        {
            onto_enrich_array = onto_enrich_array.OrderBy(l => l.Ontology_type).ThenBy(l => l.Sample_name).ThenBy(l => l.Pvalue).ToArray();
            return onto_enrich_array;
        }
        public static Ontology_enrichment_line_class[] Order_level_entryType_timepoint_sampleName(Ontology_enrichment_line_class[] onto_enrich_array)
        {
            onto_enrich_array = onto_enrich_array.OrderBy(l => l.ProcessLevel).ThenBy(l => l.Sample_name).ToArray();
            return onto_enrich_array;
        }

        public static Ontology_enrichment_line_class[] Order_entryType_timepoint_sampleName(Ontology_enrichment_line_class[] onto_enrich_array)
        {
            onto_enrich_array = onto_enrich_array.OrderBy(l => l.Sample_name).ToArray();
            return onto_enrich_array;
        }

        public static Ontology_enrichment_line_class[] Order_entryType_timepoint_sampleName_level_pvalue(Ontology_enrichment_line_class[] onto_enrich_array)
        {
            onto_enrich_array = onto_enrich_array.OrderBy(l => l.Sample_name).ThenBy(l => l.ProcessLevel).ThenBy(l => l.Pvalue).ToArray();
            return onto_enrich_array;
        }

        public bool Equal_complete_sample(Ontology_enrichment_line_class other)
        {
            bool equal = ((this.Ontology_type.Equals(other.Ontology_type))
                          && (this.Sample_name.Equals(other.Sample_name)));
            return equal;
        }
        #endregion

        public Ontology_enrichment_line_class Deep_copy()
        {
            Ontology_enrichment_line_class copy = (Ontology_enrichment_line_class)this.MemberwiseClone();
            copy.Scp_id = (string)this.Scp_id.Clone();
            copy.Scp_name = (string)this.Scp_name.Clone();
            copy.Sample_name = (string)this.Sample_name.Clone();
            copy.Parent_scp_name = (string)this.Parent_scp_name.Clone();
            copy.Integration_group = (string)this.Integration_group.Clone();
            int symbols_length = this.Overlap_symbols.Length;
            copy.Overlap_symbols = new string[symbols_length];
            for (int indexS = 0; indexS < symbols_length; indexS++)
            {
                copy.Overlap_symbols[indexS] = (string)this.Overlap_symbols[indexS].Clone();
            }
            return copy;
        }
    }

    class Ontology_enrich_readWriteOptions_class : ReadWriteOptions_base
    {
        public static char Delimiter { get { return ','; } }

        public Ontology_enrich_readWriteOptions_class(string directory, string fileName)
        {
            ReadWriteClass.Create_directory_if_it_does_not_exist(directory);
            this.File = directory + fileName;
            Key_propertyNames = new string[] { "Ontology_type", "Integration_group", "Sample_name",                           "Scp_name",                  "Experimental_symbols_count",                                   "Process_symbols_count",                               "Bg_symbol_count",            "Overlap_count",               "Pvalue",  "Minus_log10_pvalue",   "Fractional_rank", "ReadWrite_overlap_symbols" };
            Key_columnNames = new string[]   { "Ontology",      "Integration group", "Dataset and cell (sub)type or segment", "Subcellular process (SCP)", "Number of experimental genes that are part of background set", "Number of SCP genes that are part of background set", "Number of background genes", "Number of overlapping genes", "P-value", "Minus log10(p-value)", "Fractional rank", "Overlapping gene symbols" };
            File_has_headline = true;
            LineDelimiters = new char[] { Global_class.Tab };
            HeadlineDelimiters = new char[] { Global_class.Tab };
            Report = ReadWrite_report_enum.Report_main;
        }
    }

    class Ontology_enrichment_class
    {
        public Ontology_enrichment_line_class[] Enrich { get; set; }

        public Ontology_enrichment_class()
        {
            Enrich = new Ontology_enrichment_line_class[0];
        }

        #region Filter
        public void Keep_top_ranked_predictions_per_level_for_each_sample_after_calculation_of_fractional_rank(int[] max_fractional_ranks_per_level)
        {
            Calculate_fractional_ranks_for_SCPs_within_each_sampleName_processLevel();
            int enrich_length = this.Enrich.Length;
            Ontology_enrichment_line_class enrichment_line;
            List<Ontology_enrichment_line_class> keep_onto_list = new List<Ontology_enrichment_line_class>();
            this.Enrich = this.Enrich.OrderBy(l => l.Ontology_type).ThenBy(l => l.Sample_name).ThenBy(l => l.ProcessLevel).ThenByDescending(l => l.Minus_log10_pvalue).ToArray();
            int kept_lines_count = 0;
            for (int indexE = 0; indexE < enrich_length; indexE++)
            {
                enrichment_line = this.Enrich[indexE];
                if (!enrichment_line.Is_mbco_ontology())
                {
                    throw new Exception();
                }
                if (enrichment_line.Fractional_rank <= max_fractional_ranks_per_level[enrichment_line.ProcessLevel])
                {
                    kept_lines_count++;
                    keep_onto_list.Add(enrichment_line);
                }
            }
            this.Enrich = keep_onto_list.ToArray();
        }
        #endregion

        #region Get
        public string[] Get_all_distinct_scps_after_spliting_scp_unions()
        {
            int enrich_length = this.Enrich.Length;
            List<string> scpNames = new List<string>();
            Ontology_enrichment_line_class enrich_line;
            string[] scps;
            string scp;
            int scps_length;
            for (int indexE = 0; indexE < enrich_length; indexE++)
            {
                enrich_line = this.Enrich[indexE];
                scps = enrich_line.Scp_name.Split('$');
                scps_length = scps.Length;
                for (int indexS = 0; indexS < scps_length; indexS++)
                {
                    scp = scps[indexS];
                    scpNames.Add(scp);
                }
            }
            return scpNames.Distinct().OrderBy(l => l).ToArray();
        }
        #endregion

        #region Add
        public void Add_other_lines(Ontology_enrichment_line_class[] other_lines)
        {
            int this_enrich_length = this.Enrich.Length;
            int other_enrich_length = other_lines.Length;
            int new_enrich_length = this_enrich_length + other_enrich_length;
            Ontology_enrichment_line_class[] new_enrich = new Ontology_enrichment_line_class[new_enrich_length];
            int indexNew = -1;
            for (int indexThis = 0; indexThis < this_enrich_length; indexThis++)
            {
                indexNew++;
                new_enrich[indexNew] = this.Enrich[indexThis];
            }
            for (int indexOther = 0; indexOther < other_enrich_length; indexOther++)
            {
                indexNew++;
                new_enrich[indexNew] = other_lines[indexOther].Deep_copy();
            }
            this.Enrich = new_enrich;
        }

        public void Add_other(Ontology_enrichment_class other)
        {
            Add_other_lines(other.Enrich);
        }
        #endregion

        #region Seperate SCPs predicted by dynamic enrichment analysis
        public void Separate_scp_unions_into_single_scps_and_keep_line_defined_by_lowest_pvalue_for_each_scp_and_add_scp_specific_genes(Ontology_enrichment_class standard_unfiltered)
        {
            int enrich_length = this.Enrich.Length;
            Ontology_enrichment_line_class onto_enrich_line;
            Ontology_enrichment_line_class singleScp_onto_enrich_line;
            List<Ontology_enrichment_line_class> onto_enrich_keep = new List<Ontology_enrichment_line_class>();
            Dictionary<string, bool> considered_scps_of_current_condition = new Dictionary<string, bool>();
            this.Enrich = this.Enrich.OrderBy(l => l.Sample_name).ThenBy(l => l.Pvalue).ToArray();
            string[] scps;
            string scp;
            int scps_length;
            for (int indexO = 0; indexO < enrich_length; indexO++)
            {
                onto_enrich_line = this.Enrich[indexO];
                if ((indexO == 0)
                    || (!onto_enrich_line.Sample_name.Equals(this.Enrich[indexO - 1].Sample_name)))
                {
                    considered_scps_of_current_condition.Clear();
                }
                scps = onto_enrich_line.Scp_name.Split('$');
                scps_length = scps.Length;
                for (int indexScp = 0; indexScp < scps_length; indexScp++)
                {
                    scp = scps[indexScp];
                    if (!considered_scps_of_current_condition.ContainsKey(scp))
                    {
                        considered_scps_of_current_condition.Add(scp, true);
                        singleScp_onto_enrich_line = onto_enrich_line.Deep_copy();
                        singleScp_onto_enrich_line.Scp_name = (string)scp.Clone();
                        singleScp_onto_enrich_line.Scp_id = "broken up dynamic results from " + onto_enrich_line.Scp_name;
                        singleScp_onto_enrich_line.Overlap_symbols = new string[0];
                        onto_enrich_keep.Add(singleScp_onto_enrich_line);
                    }
                }
            }
            //   if (onto_enrich_keep.Count<=this.Enrich.Length) { throw new Exception(); }
            this.Enrich = onto_enrich_keep.ToArray();

            this.Enrich = this.Enrich.OrderBy(l => l.Sample_name).ThenBy(l => l.Scp_name).ToArray();
            standard_unfiltered.Enrich = standard_unfiltered.Enrich.OrderBy(l => l.Sample_name).ThenBy(l => l.Scp_name).ToArray();

            int indexStandard = 0;
            int this_length = this.Enrich.Length;
            int standard_length = standard_unfiltered.Enrich.Length;
            Ontology_enrichment_line_class this_onto_enrich_line;
            Ontology_enrichment_line_class standard_onto_enrich_line = new Ontology_enrichment_line_class();
            int stringCompare;
            for (int indexThis = 0; indexThis < this_length; indexThis++)
            {
                this_onto_enrich_line = this.Enrich[indexThis];
                stringCompare = -2;
                while ((indexStandard < standard_length) && (stringCompare < 0))
                {
                    standard_onto_enrich_line = standard_unfiltered.Enrich[indexStandard];
                    stringCompare = standard_onto_enrich_line.Sample_name.CompareTo(this_onto_enrich_line.Sample_name);
                    if (stringCompare == 0)
                    {
                        stringCompare = standard_onto_enrich_line.Scp_name.CompareTo(this_onto_enrich_line.Scp_name);
                    }
                    if (stringCompare < 0) { indexStandard++; }
                }
                if (stringCompare != 0) { throw new Exception(); }
                this_onto_enrich_line.Overlap_symbols = Array_class.Deep_copy_string_array(standard_onto_enrich_line.Overlap_symbols);
            }
        }
        #endregion

        #region Rank
        public void Calculate_fractional_ranks_for_SCPs_within_each_sampleName_processLevel()
        {
            this.Enrich = this.Enrich.OrderBy(l => l.Sample_name).ThenBy(l => l.ProcessLevel).ThenBy(l => l.Pvalue).ToArray();
            int enrich_length = this.Enrich.Length;
            Ontology_enrichment_line_class enrich_line;
            Ontology_enrichment_line_class inner_enrich_line;
            int ordinal_rank = 0;
            int firstIndexSamePvalue = -1;
            List<float> current_ordinal_ranks = new List<float>();
            float fractional_rank;
            for (int indexE = 0; indexE < enrich_length; indexE++)
            {
                enrich_line = this.Enrich[indexE];
                if ((indexE == 0)
                    || (!enrich_line.Sample_name.Equals(this.Enrich[indexE - 1].Sample_name))
                    || (!enrich_line.ProcessLevel.Equals(this.Enrich[indexE - 1].ProcessLevel)))
                {
                    ordinal_rank = 0;
                }
                if ((indexE == 0)
                    || (!enrich_line.Sample_name.Equals(this.Enrich[indexE - 1].Sample_name))
                    || (!enrich_line.ProcessLevel.Equals(this.Enrich[indexE - 1].ProcessLevel))
                    || (!enrich_line.Pvalue.Equals(this.Enrich[indexE - 1].Pvalue)))
                {
                    current_ordinal_ranks.Clear();
                    firstIndexSamePvalue = indexE;
                }
                ordinal_rank++;
                current_ordinal_ranks.Add(ordinal_rank);
                if ((indexE == enrich_length - 1)
                    || (!enrich_line.Sample_name.Equals(this.Enrich[indexE + 1].Sample_name))
                    || (!enrich_line.ProcessLevel.Equals(this.Enrich[indexE + 1].ProcessLevel))
                    || (!enrich_line.Pvalue.Equals(this.Enrich[indexE + 1].Pvalue)))
                {
                    if (current_ordinal_ranks.Count > 1)
                    {
                        fractional_rank = Math_class.Get_average(current_ordinal_ranks.ToArray());
                        for (int indexInner = firstIndexSamePvalue; indexInner <= indexE; indexInner++)
                        {
                            inner_enrich_line = this.Enrich[indexInner];
                            inner_enrich_line.Fractional_rank = fractional_rank;
                        }
                    }
                    else if (current_ordinal_ranks.Count == 1)
                    {
                        if (firstIndexSamePvalue != indexE) { throw new Exception(); }
                        enrich_line.Fractional_rank = current_ordinal_ranks[0];
                    }
                    else { throw new Exception(); }
                }
            }
        }
        #endregion

        #region Write copy
        public void Write(string directory, string file_name)
        {
            this.Enrich = this.Enrich.OrderBy(l => l.Sample_name).ThenBy(l => l.ProcessLevel).ThenBy(l => l.Pvalue).ToArray();
            Ontology_enrich_readWriteOptions_class readWriteOptions = new Ontology_enrich_readWriteOptions_class(directory, file_name);
            ReadWriteClass.WriteData(Enrich, readWriteOptions);
        }

        public Ontology_enrichment_class Deep_copy()
        {
            Ontology_enrichment_class copy = (Ontology_enrichment_class)this.MemberwiseClone();
            int enrich_length = this.Enrich.Length;
            copy.Enrich = new Ontology_enrichment_line_class[enrich_length];
            for (int indexE = 0; indexE < enrich_length; indexE++)
            {
                copy.Enrich[indexE] = this.Enrich[indexE].Deep_copy();
            }
            return copy;
        }
        #endregion
    }
}
