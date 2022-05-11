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
using ReadWrite;
using Enumerations;
using Highthroughput_data;
using Statistic;
using System.IO;

namespace KPMP
{
    enum Statistical_score_of_interest_type_enum { E_m_p_t_y, Fractional_rank, Value_normalized_by_max, Value };

    class KPMP_standardized_dataset_line_class
    {
        public string Gene_symbol { get; set; }
        public string PatientId { get; set; }
        public string Dataset { get; set; }
        public string Cell_segment { get; set; }
        public string KPMP_data_integration_term { get; set; }
        public double Value_1st { get; set; }
        public double Value_2nd { get; set; }
        public float Fractional_rank_based_on_value_types { get; set; }
        public double Score_of_interest { get; set; }
        public Statistical_score_of_interest_type_enum Score_of_interest_type { get; set; }
        public KPMP_value_type_enum Value_type_1st { get; set; }
        public KPMP_value_type_enum Value_type_2nd { get; set; }

        public KPMP_standardized_dataset_line_class()
        {
            Gene_symbol = "";
            PatientId = "";
            Dataset = "";
            Cell_segment = "";
            KPMP_data_integration_term = "";
            Fractional_rank_based_on_value_types = -1;
            Value_1st = -1;
            Value_2nd = -1;
        }

        public static KPMP_standardized_dataset_line_class[] Order_by_condition_and_geneSymbol(KPMP_standardized_dataset_line_class[] standardized_data)
        {
            standardized_data = standardized_data.OrderBy(l => l.Dataset).ThenBy(l => l.Cell_segment).ThenBy(l => l.PatientId).ThenBy(l => l.Value_type_1st).ThenBy(l=>l.Value_type_2nd).ThenBy(l => l.Gene_symbol).ToArray();
            return standardized_data;
        }

        public bool Is_equal_condition(KPMP_standardized_dataset_line_class other)
        {
            bool equal = (this.Dataset.Equals(other.Dataset))
                          && (this.Cell_segment.Equals(other.Cell_segment))
                          && (this.PatientId.Equals(other.PatientId))
                          && (this.Value_type_1st.Equals(other.Value_type_1st))
                          && (this.Value_type_2nd.Equals(other.Value_type_2nd));
            return equal;
        }

        public KPMP_standardized_dataset_line_class Deep_copy()
        {
            KPMP_standardized_dataset_line_class copy = (KPMP_standardized_dataset_line_class)this.MemberwiseClone();
            copy.Gene_symbol = (string)this.Gene_symbol.Clone();
            copy.PatientId = (string)this.PatientId.Clone();
            copy.Dataset = (string)this.Dataset.Clone();
            copy.Cell_segment = (string)this.Cell_segment.Clone();
            copy.KPMP_data_integration_term = (string)this.KPMP_data_integration_term.Clone();
            return copy;
        }
    }

    class KPMP_standardized_bgGeneProtein_line_class
    {
        public string Dataset { get; set; }
        public string BgGeneProtein { get; set; }
    }

    class KPMP_standardized_dataset_readWriteOptions_class : ReadWriteOptions_base
    {
        public KPMP_standardized_dataset_readWriteOptions_class(string directory, string fileName)
        {
            this.File = directory + fileName;
            this.Key_propertyNames = new string[] { "PatientId", "Dataset", "Cell_segment", "KPMP_data_integration_term", "Gene_symbol", "Value_1st", "Value_type_1st", "Value_2nd", "Value_type_2nd", "Fractional_rank_based_on_value_types" };//,"Fractional_rank","Score_of_interest","Score_of_interest_type" };
            this.Key_columnNames = this.Key_propertyNames;
            this.HeadlineDelimiters = new char[] { Global_class.Tab };
            this.LineDelimiters = new char[] { Global_class.Tab };
            this.File_has_headline = true;
            this.Report = ReadWrite_report_enum.Report_main;
        }
    }

    class KPMP_standardized_dataset_bgGenesProteins_readWriteOptions_class : ReadWriteOptions_base
    {
        public KPMP_standardized_dataset_bgGenesProteins_readWriteOptions_class(string directory, string fileName)
        {
            this.File = directory + fileName;
            this.Key_propertyNames = new string[] { "Dataset", "BgGeneProtein" };
            this.Key_columnNames = this.Key_propertyNames;
            this.HeadlineDelimiters = new char[] { Global_class.Tab };
            this.LineDelimiters = new char[] { Global_class.Tab };
            this.File_has_headline = true;
            this.Report = ReadWrite_report_enum.Report_main;
        }
    }

    class KPMP_standardized_dataset_class
    {
        public KPMP_standardized_dataset_line_class[] KPMP_data { get; set; }
        public Dictionary<string,string[]> Dataset_bgGenesProteins_dict { get; set; }
        public bool Data_filtered_by_cutoff { get; set; }

        public KPMP_standardized_dataset_class()
        {
            this.KPMP_data = new KPMP_standardized_dataset_line_class[0];
            Dataset_bgGenesProteins_dict = new Dictionary<string, string[]>();
            Data_filtered_by_cutoff = false;
        }

        #region Check
        public void Check_for_duplicates()
        {
            int this_length = this.KPMP_data.Length;
            this.KPMP_data = this.KPMP_data.OrderBy(l => l.Cell_segment).ThenBy(l => l.Dataset).ThenBy(l => l.Gene_symbol).ThenBy(l => l.PatientId).ToArray();
            KPMP_standardized_dataset_line_class this_line;
            for (int indexThis = 0; indexThis < this_length; indexThis++)
            {
                this_line = this.KPMP_data[indexThis];
                if (  (indexThis!=0)
                    && (this_line.Cell_segment.Equals(this.KPMP_data[indexThis - 1].Cell_segment))
                    && (this_line.Dataset.Equals(this.KPMP_data[indexThis - 1].Dataset))
                    && (this_line.Gene_symbol.Equals(this.KPMP_data[indexThis - 1].Gene_symbol))
                    && (this_line.PatientId.Equals(this.KPMP_data[indexThis - 1].PatientId)))
                {
                    throw new Exception();
                }
            }
        }

        public void Check_if_all_values_positive_or_zero()
        {
            int zero_values = 0;
            foreach (KPMP_standardized_dataset_line_class kpmp_data_line in this.KPMP_data)
            {
                if (kpmp_data_line.Value_1st<0) { throw new Exception(); }
                if (kpmp_data_line.Value_2nd<0) { throw new Exception(); }
                if (kpmp_data_line.Value_2nd==0) { zero_values++; }
            }
        }

        public void Check_if_all_fractional_ranks_are_below_cutoff(float fractional_rank_cutoff)
        {
            foreach (KPMP_standardized_dataset_line_class kpmp_data_line in this.KPMP_data)
            {
                if (kpmp_data_line.Fractional_rank_based_on_value_types > fractional_rank_cutoff)
                {
                    throw new Exception();
                }
            }
        }

        public void Check_if_all_pvalues_adjusted_pvalues_are_below_cutoff(float pvalue_adj_pvalue_cutoff)
        {
            foreach (KPMP_standardized_dataset_line_class kpmp_data_line in this.KPMP_data)
            {
                switch (kpmp_data_line.Value_type_1st)
                {
                    case KPMP_value_type_enum.Average:
                    case KPMP_value_type_enum.Fdr:
                    case KPMP_value_type_enum.Log2_ratioavg:
                    case KPMP_value_type_enum.Log2_ratio_singlepatient:
                    case KPMP_value_type_enum.Log2_single_value:
                    case KPMP_value_type_enum.Median:
                    case KPMP_value_type_enum.Ratioavg:
                    case KPMP_value_type_enum.Ratio_singlepatient:
                        break;
                    case KPMP_value_type_enum.Pvalue:
                    case KPMP_value_type_enum.Pvalue_adjusted:
                        if (kpmp_data_line.Value_1st > pvalue_adj_pvalue_cutoff) { throw new Exception(); }
                        break;
                    case KPMP_value_type_enum.Minus_log10_fdr:
                    case KPMP_value_type_enum.Minus_log10_pvalue:
                    case KPMP_value_type_enum.Minus_log10_pvalue_adjusted:
                        if (kpmp_data_line.Value_1st<(-Math.Log10(pvalue_adj_pvalue_cutoff))) { throw new Exception(); }
                        break;
                    default:
                        throw new Exception();

                }
            }
        }

        public void Check_if_all_values_are_positive()
        {
            foreach (KPMP_standardized_dataset_line_class kpmp_data_line in this.KPMP_data)
            {
                if (kpmp_data_line.Value_1st <= 0) { throw new Exception(); }
                if (kpmp_data_line.Value_2nd <= 0) { throw new Exception(); }
            }
        }

        public void Check_if_no_values_is_NaN()
        {
            foreach (KPMP_standardized_dataset_line_class kpmp_data_line in this.KPMP_data)
            {
                if (Double.IsNaN(kpmp_data_line.Value_1st)) { throw new Exception(); }
                if (Double.IsNaN(kpmp_data_line.Value_2nd)) { throw new Exception(); }
            }
        }
        #endregion

        #region Add
        public void Add_to_array(KPMP_standardized_dataset_line_class[] add_kpmp_data)
        {
            int this_length = this.KPMP_data.Length;
            int add_length = add_kpmp_data.Length;
            int new_length = this_length + add_length;
            KPMP_standardized_dataset_line_class[] new_kpmp_data = new KPMP_standardized_dataset_line_class[new_length];
            int indexNew = -1;
            for (int indexThis = 0; indexThis < this_length; indexThis++)
            {
                indexNew++;
                new_kpmp_data[indexNew] = this.KPMP_data[indexThis];
            }
            for (int indexAdd = 0; indexAdd < add_length; indexAdd++)
            {
                indexNew++;
                new_kpmp_data[indexNew] = add_kpmp_data[indexAdd];
            }
            this.KPMP_data = new_kpmp_data;
        }

        public void Add_deep_copy_to_array(KPMP_standardized_dataset_line_class[] add_kpmp_data)
        {
            int this_length = this.KPMP_data.Length;
            int add_length = add_kpmp_data.Length;
            int new_length = this_length + add_length;
            KPMP_standardized_dataset_line_class[] new_kpmp_data = new KPMP_standardized_dataset_line_class[new_length];
            int indexNew = -1;
            for (int indexThis = 0; indexThis < this_length; indexThis++)
            {
                indexNew++;
                new_kpmp_data[indexNew] = this.KPMP_data[indexThis];
            }
            for (int indexAdd = 0; indexAdd < add_length; indexAdd++)
            {
                indexNew++;
                new_kpmp_data[indexNew] = add_kpmp_data[indexAdd].Deep_copy();
            }
            this.KPMP_data = new_kpmp_data;
        }

        private void Add_deep_copy_to_dictionary(Dictionary<string,string[]> add_dataset_bgGenesProteins_dict)
        {
            string[] add_keys = add_dataset_bgGenesProteins_dict.Keys.ToArray();
            string add_key;
            int add_keys_length = add_keys.Length;
            for (int indexAdd=0; indexAdd<add_keys_length;indexAdd++)
            {
                add_key = (string)add_keys[indexAdd].Clone();
                if (!this.Dataset_bgGenesProteins_dict.ContainsKey(add_key))
                {
                    this.Dataset_bgGenesProteins_dict.Add(add_key, Array_class.Deep_copy_string_array(add_dataset_bgGenesProteins_dict[add_key]));
                }
            }
        }

        public void Add_to_existing_instances(KPMP_standardized_dataset_line_class[] standard_data_lines, Dictionary<string,string[]> dataset_bgGenesProteins_dict)
        {
            Add_deep_copy_to_dictionary(dataset_bgGenesProteins_dict);
            Add_to_array(standard_data_lines);
        }

        public void Add_deep_copy_of_other_after_checking_if_both_have_not_been_filtered_by_cutoff(KPMP_standardized_dataset_class other)
        {
            if ((this.Data_filtered_by_cutoff)||(other.Data_filtered_by_cutoff)) { throw new Exception(); }
            Add_deep_copy_to_array(other.KPMP_data);
            Add_deep_copy_to_dictionary(other.Dataset_bgGenesProteins_dict);
        }
        #endregion

        private Fill_de_line_class Generate_fill_de_line_based_on_one_of_indicated_value_types_check_if_duplicated(KPMP_standardized_dataset_line_class kpmp_data_line, KPMP_value_type_enum[] potential_valueTypes)
        {
            Fill_de_line_class fill_de_line = new Fill_de_line_class();
            fill_de_line.Symbols_for_de = new string[] { (string)kpmp_data_line.Gene_symbol.Clone() };
            fill_de_line.Names_for_de = KPMP_data_integration_class.Get_names_for_de_in_correct_order_separated_by_minus(kpmp_data_line.Cell_segment, kpmp_data_line.PatientId, kpmp_data_line.Dataset, kpmp_data_line.Value_type_1st, kpmp_data_line.KPMP_data_integration_term);
            bool value_filled = false;
            foreach (KPMP_value_type_enum potential_valueType in potential_valueTypes)
            {
                if (kpmp_data_line.Value_type_1st == potential_valueType)
                {
                    if (value_filled) { throw new Exception(); }
                    value_filled = true;
                    fill_de_line.Value_for_de = kpmp_data_line.Value_1st;
                }
                if (kpmp_data_line.Value_type_2nd == potential_valueType)
                {
                    if (value_filled) { throw new Exception(); }
                    value_filled = true;
                    fill_de_line.Value_for_de = kpmp_data_line.Value_2nd;
                }
            }
            if (!value_filled) { throw new Exception(); }
            return fill_de_line;
        }

        public DE_class Generate_de_instance_alternatively_and_fill_with_one_of_indicated_valueTypes_and_check_if_duplicated(params KPMP_value_type_enum[] valueTypes)
        {
            Fill_de_line_class fill_de_line;
            List<Fill_de_line_class> fill_de_list = new List<Fill_de_line_class>();
            foreach (KPMP_standardized_dataset_line_class kpmp_data_line in KPMP_data)
            {
                fill_de_line = Generate_fill_de_line_based_on_one_of_indicated_value_types_check_if_duplicated(kpmp_data_line, valueTypes);
                fill_de_list.Add(fill_de_line);
            }
            DE_class de = new DE_class();
            de.Fill_with_data_alternatively(fill_de_list.ToArray());
            return de;
        }

        public string[] Get_bgGenes_of_indicated_dataset(string dataset)
        {
            string[] bgGenes = new string[0];
            if (Dataset_bgGenesProteins_dict.ContainsKey(dataset)) { bgGenes = Array_class.Deep_copy_string_array(Dataset_bgGenesProteins_dict[dataset]); }
            return bgGenes;
        }

        public void Add_missing_genes_to_each_dataset_condition_if_valueType1st_is_Average_or_singleValue_ifNot_throwException()
        {
            string[] all_genes = Get_all_unique_ordered_genes();
            string current_gene;
            int all_genes_length = all_genes.Length;
            int indexCurrentGene = 0;
            int stringCompare = -2;
            int data_length = this.KPMP_data.Length;
            KPMP_standardized_dataset_line_class standardized_line;
            KPMP_standardized_dataset_line_class new_standardized_line;
            List<KPMP_standardized_dataset_line_class> new_lines = new List<KPMP_standardized_dataset_line_class>();
            this.KPMP_data = KPMP_standardized_dataset_line_class.Order_by_condition_and_geneSymbol(this.KPMP_data);
            int currentCondition_geneCount = 0;
            int allCondition_geneCount = -1;
            for (int indexData=0; indexData<data_length;indexData++)
            {
                standardized_line = this.KPMP_data[indexData];
                switch (standardized_line.Value_type_1st)
                {
                    case KPMP_value_type_enum.Average:
                    case KPMP_value_type_enum.Single_value:
                        break;
                    default:
                        throw new Exception();
                }
                stringCompare = -2;
                if (  (indexData==0)
                    ||(!standardized_line.Is_equal_condition(this.KPMP_data[indexData-1])))
                {
                    indexCurrentGene = 0;
                    currentCondition_geneCount = 0;
                }
                currentCondition_geneCount++;
                stringCompare = -2;
                while (stringCompare<0)
                {
                    current_gene = all_genes[indexCurrentGene];
                    stringCompare = current_gene.CompareTo(standardized_line.Gene_symbol);
                    if (stringCompare<0)
                    {
                        new_standardized_line = new KPMP_standardized_dataset_line_class();
                        new_standardized_line.Dataset = (string)standardized_line.Dataset.Clone();
                        new_standardized_line.Cell_segment = (string)standardized_line.Cell_segment.Clone();
                        new_standardized_line.Gene_symbol = (string)current_gene.Clone();
                        new_standardized_line.KPMP_data_integration_term = (string)standardized_line.KPMP_data_integration_term.Clone();
                        new_standardized_line.PatientId = (string)standardized_line.PatientId.Clone();
                        new_standardized_line.Value_type_1st = standardized_line.Value_type_1st;
                        new_standardized_line.Value_type_2nd = KPMP_value_type_enum.No_selection;
                        new_standardized_line.Value_1st = 0;
                        new_standardized_line.Value_2nd = 0;
                        new_lines.Add(new_standardized_line);
                        indexCurrentGene++;
                        currentCondition_geneCount++;
                    }
                    else if (stringCompare==0)
                    {
                        indexCurrentGene++;
                    }
                }
                if ((indexData == data_length-1)
                    || (!standardized_line.Is_equal_condition(this.KPMP_data[indexData + 1])))
                {
                    while (indexCurrentGene<all_genes_length)
                    {
                        current_gene = all_genes[indexCurrentGene];
                        new_standardized_line = new KPMP_standardized_dataset_line_class();
                        new_standardized_line.Dataset = (string)standardized_line.Dataset.Clone();
                        new_standardized_line.Cell_segment = (string)standardized_line.Cell_segment.Clone();
                        new_standardized_line.Gene_symbol = (string)current_gene.Clone();
                        new_standardized_line.KPMP_data_integration_term = (string)standardized_line.KPMP_data_integration_term.Clone();
                        new_standardized_line.PatientId = (string)standardized_line.PatientId.Clone();
                        new_standardized_line.Value_type_1st = standardized_line.Value_type_1st;
                        new_standardized_line.Value_type_2nd = KPMP_value_type_enum.No_selection;
                        new_standardized_line.Value_1st = 0;
                        new_standardized_line.Value_2nd = 0;
                        new_lines.Add(new_standardized_line);
                        indexCurrentGene++;
                        currentCondition_geneCount++;
                    }
                }
                if ((indexData == data_length - 1)
                    || (!standardized_line.Is_equal_condition(this.KPMP_data[indexData + 1])))
                {
                    if (allCondition_geneCount==-1)
                    {
                        allCondition_geneCount = currentCondition_geneCount;
                    }
                    else if (allCondition_geneCount != currentCondition_geneCount)
                    {
                        throw new Exception();
                    }
                }

            }
            this.Add_to_array(new_lines.ToArray());
        }

        public void Calculate_ratios_between_indicated_integration_terms_and_add_assuming_only_singleValues1st_for_term0_and_multiple_for_term1(string integration_term0, string integration_term1)
        {
            int kpmp_data_length = this.KPMP_data.Length;
            KPMP_standardized_dataset_line_class kpmp_data_line;
            this.KPMP_data = this.KPMP_data.OrderBy(l => l.PatientId).ThenBy(l => l.Dataset).ThenBy(l => l.Gene_symbol).ThenBy(l => l.Value_type_1st).ThenBy(l => l.KPMP_data_integration_term).ToArray();
            KPMP_standardized_dataset_line_class ratio_kpmp_data_line_base;
            KPMP_standardized_dataset_line_class ratio_kpmp_data_line;
            List<KPMP_standardized_dataset_line_class> ratio_kpmp_data_list = new List<KPMP_standardized_dataset_line_class>();
            double value_term0 = -1;
            List<double> value_term1s = new List<double>();
            List<string> cell_segment1s = new List<string>();
            double value_term1;
            string cell_segment1;
            string cell_segment0 = "error";
            int cell_segment1s_length;
            string ratio_cell_segment1_name;
            for (int indexData=0; indexData<kpmp_data_length; indexData++)
            {
                kpmp_data_line = this.KPMP_data[indexData];
                if (!kpmp_data_line.Value_type_2nd.Equals(KPMP_value_type_enum.No_selection)) { throw new Exception(); }
                if (!kpmp_data_line.Value_type_1st.Equals(KPMP_value_type_enum.Single_value)) { throw new Exception(); }
                if ((indexData==0)
                    || (!kpmp_data_line.PatientId.Equals(this.KPMP_data[indexData - 1].PatientId))
                    //|| (!kpmp_data_line.KPMP_data_integration_term.Equals(this.KPMP_data[indexData - 1].KPMP_data_integration_term))
                    || (!kpmp_data_line.Dataset.Equals(this.KPMP_data[indexData - 1].Dataset))
                    || (!kpmp_data_line.Gene_symbol.Equals(this.KPMP_data[indexData - 1].Gene_symbol))
                    || (!kpmp_data_line.Value_type_1st.Equals(this.KPMP_data[indexData - 1].Value_type_1st)))
                {
                    value_term0 = -1;
                    cell_segment0 = "error";
                    value_term1s.Clear();
                    cell_segment1s.Clear();
                }
                if (kpmp_data_line.KPMP_data_integration_term.Equals(integration_term0))
                {
                    if (value_term0 != -1) { throw new Exception(); }
                    value_term0 = kpmp_data_line.Value_1st;
                    cell_segment0 = (string)kpmp_data_line.Cell_segment.Clone();
                }
                if (kpmp_data_line.KPMP_data_integration_term.IndexOf(integration_term1)==0)
                {
                    value_term1s.Add(kpmp_data_line.Value_1st);
                    cell_segment1s.Add(kpmp_data_line.Cell_segment);
                }
                if ((indexData == kpmp_data_length - 1)
                    || (!kpmp_data_line.PatientId.Equals(this.KPMP_data[indexData + 1].PatientId))
                    // || (!kpmp_data_line.KPMP_data_integration_term.Equals(this.KPMP_data[indexData + 1].KPMP_data_integration_term))
                    || (!kpmp_data_line.Dataset.Equals(this.KPMP_data[indexData + 1].Dataset))
                    || (!kpmp_data_line.Gene_symbol.Equals(this.KPMP_data[indexData + 1].Gene_symbol))
                    || (!kpmp_data_line.Value_type_1st.Equals(this.KPMP_data[indexData + 1].Value_type_1st)))
                {
                    if ((value_term0 == -1) || (value_term1s.Count == 0))
                    {
                        throw new Exception();
                    }
                    if (cell_segment1s.Count != cell_segment1s.Distinct().ToList().Count)
                    {
                        throw new Exception();
                    }
                    cell_segment1s_length = cell_segment1s.Count;
                    for (int index1 = 0; index1 < cell_segment1s_length; index1++)
                    {
                        cell_segment1 = cell_segment1s[index1];
                        value_term1 = value_term1s[index1];
                        ratio_cell_segment1_name = cell_segment1;

                        ratio_kpmp_data_line_base = new KPMP_standardized_dataset_line_class();
                        ratio_kpmp_data_line_base.Dataset = (string)kpmp_data_line.Dataset.Clone();
                        ratio_kpmp_data_line_base.Gene_symbol = (string)kpmp_data_line.Gene_symbol.Clone();
                        ratio_kpmp_data_line_base.PatientId = (string)kpmp_data_line.PatientId.Clone();
                        ratio_kpmp_data_line_base.Value_2nd = -1;
                        ratio_kpmp_data_line_base.Value_type_2nd = KPMP_value_type_enum.No_selection;
                        switch (kpmp_data_line.Value_type_1st)
                        {
                            case KPMP_value_type_enum.Single_value:
                                ratio_kpmp_data_line_base.Value_type_1st = KPMP_value_type_enum.Ratio_singlepatient;
                                break;
                            case KPMP_value_type_enum.Average:
                                ratio_kpmp_data_line_base.Value_type_1st = KPMP_value_type_enum.Ratioavg;
                                break;
                            default:
                                throw new Exception();
                        }
                        value_term0 = value_term0 + 1;
                        value_term1 = value_term1 + 1;
                        ratio_kpmp_data_line = ratio_kpmp_data_line_base.Deep_copy();
                        ratio_kpmp_data_line.Value_1st = value_term0 / value_term1;
                        if ((double.IsNaN(ratio_kpmp_data_line.Value_1st))
                            || (double.IsInfinity(ratio_kpmp_data_line.Value_1st)))
                        {
                            throw new Exception();
                        }
                        ratio_kpmp_data_line.Cell_segment = (string)cell_segment0.Clone() + " vs " + (string)cell_segment1.Clone();
                        ratio_kpmp_data_line.KPMP_data_integration_term = (string)integration_term0.Clone();
                        ratio_kpmp_data_list.Add(ratio_kpmp_data_line);

                        ratio_kpmp_data_line = ratio_kpmp_data_line_base.Deep_copy();
                        ratio_kpmp_data_line.Value_1st = value_term1 / value_term0;
                        if ((double.IsNaN(ratio_kpmp_data_line.Value_1st))
                            || (double.IsInfinity(ratio_kpmp_data_line.Value_1st)))
                        {
                            throw new Exception();
                        }
                        ratio_kpmp_data_line.Cell_segment = (string)cell_segment1.Clone() + " vs " + (string)cell_segment0.Clone();
                        ratio_kpmp_data_line.KPMP_data_integration_term = (string)integration_term1.Clone();
                        ratio_kpmp_data_list.Add(ratio_kpmp_data_line);
                    }
                }
            }
            Add_to_array(ratio_kpmp_data_list.ToArray());
        }

        public void Check_if_recalcuated_fractional_ranks_based_on_descending_values1st_and_values2nd_agree_with_existing()
        {
            int data_length = this.KPMP_data.Length;
            KPMP_standardized_dataset_line_class data_line;
            KPMP_standardized_dataset_line_class inner_data_line;
            this.KPMP_data = this.KPMP_data.OrderBy(l => l.Dataset).ThenBy(l => l.Cell_segment).ThenBy(l => l.PatientId).ThenBy(l => l.KPMP_data_integration_term).ThenBy(l => l.Value_type_1st).ThenBy(l => l.Value_type_2nd).ThenByDescending(l => l.Value_1st).ThenByDescending(l => l.Value_2nd).ToArray();
            int firstIndexSameValue = 0;
            int currentOrdinalRank = 0;
            List<float> sameValues_1st2nd_ordinalRanks = new List<float>();
            float final_rank;

            for (int indexKPMP = 0; indexKPMP < data_length; indexKPMP++)
            {
                data_line = this.KPMP_data[indexKPMP];
                if ((data_line.Value_1st == -1)
                    || (data_line.Value_2nd == -1)
                    || (data_line.Value_type_1st.Equals(KPMP_value_type_enum.E_m_p_t_y))
                    || (data_line.Value_type_2nd.Equals(KPMP_value_type_enum.E_m_p_t_y)))
                {
                    throw new Exception();
                }
                if ((indexKPMP == 0)
                    || (!data_line.Dataset.Equals(this.KPMP_data[indexKPMP - 1].Dataset))
                    || (!data_line.Cell_segment.Equals(this.KPMP_data[indexKPMP - 1].Cell_segment))
                    || (!data_line.PatientId.Equals(this.KPMP_data[indexKPMP - 1].PatientId))
                    || (!data_line.KPMP_data_integration_term.Equals(this.KPMP_data[indexKPMP - 1].KPMP_data_integration_term))
                    || (!data_line.Value_type_1st.Equals(this.KPMP_data[indexKPMP - 1].Value_type_1st))
                    || (!data_line.Value_type_2nd.Equals(this.KPMP_data[indexKPMP - 1].Value_type_2nd)))
                {
                    currentOrdinalRank = 0;
                }
                if ((indexKPMP == 0)
                    || (!data_line.Dataset.Equals(this.KPMP_data[indexKPMP - 1].Dataset))
                    || (!data_line.Cell_segment.Equals(this.KPMP_data[indexKPMP - 1].Cell_segment))
                    || (!data_line.PatientId.Equals(this.KPMP_data[indexKPMP - 1].PatientId))
                    || (!data_line.KPMP_data_integration_term.Equals(this.KPMP_data[indexKPMP - 1].KPMP_data_integration_term))
                    || (!data_line.Value_type_1st.Equals(this.KPMP_data[indexKPMP - 1].Value_type_1st))
                    || (!data_line.Value_type_2nd.Equals(this.KPMP_data[indexKPMP - 1].Value_type_2nd))
                    || (!data_line.Value_1st.Equals(this.KPMP_data[indexKPMP - 1].Value_1st))
                    || (!data_line.Value_2nd.Equals(this.KPMP_data[indexKPMP - 1].Value_2nd)))
                {
                    firstIndexSameValue = indexKPMP;
                    sameValues_1st2nd_ordinalRanks.Clear();
                }
                currentOrdinalRank++;
                sameValues_1st2nd_ordinalRanks.Add(currentOrdinalRank);
                if ((indexKPMP == data_length - 1)
                    || (!data_line.Dataset.Equals(this.KPMP_data[indexKPMP + 1].Dataset))
                    || (!data_line.Cell_segment.Equals(this.KPMP_data[indexKPMP + 1].Cell_segment))
                    || (!data_line.PatientId.Equals(this.KPMP_data[indexKPMP + 1].PatientId))
                    || (!data_line.KPMP_data_integration_term.Equals(this.KPMP_data[indexKPMP + 1].KPMP_data_integration_term))
                    || (!data_line.Value_type_1st.Equals(this.KPMP_data[indexKPMP + 1].Value_type_1st))
                    || (!data_line.Value_type_2nd.Equals(this.KPMP_data[indexKPMP + 1].Value_type_2nd))
                    || (!data_line.Value_1st.Equals(this.KPMP_data[indexKPMP + 1].Value_1st))
                    || (!data_line.Value_2nd.Equals(this.KPMP_data[indexKPMP + 1].Value_2nd)))
                {
                    if (sameValues_1st2nd_ordinalRanks.Count > 1)
                    {
                        final_rank = Math_class.Get_average(sameValues_1st2nd_ordinalRanks.ToArray());
                        for (int indexInner = firstIndexSameValue; indexInner <= indexKPMP; indexInner++)
                        {
                            inner_data_line = this.KPMP_data[indexInner];
                            if (inner_data_line.Fractional_rank_based_on_value_types != final_rank) { throw new Exception(); }
                        }
                    }
                    else if (sameValues_1st2nd_ordinalRanks.Count == 1)
                    {
                        if (firstIndexSameValue != indexKPMP) { throw new Exception(); }
                        if (data_line.Fractional_rank_based_on_value_types != sameValues_1st2nd_ordinalRanks[0]) { throw new Exception(); }
                    }
                    else { throw new Exception(); }
                }
            }
        }

        public void Calculate_fractional_ranks_based_on_descending_values1st_and_values2nd()
        {
            int data_length = this.KPMP_data.Length;
            KPMP_standardized_dataset_line_class data_line;
            KPMP_standardized_dataset_line_class inner_data_line;
            this.KPMP_data = this.KPMP_data.OrderBy(l => l.Dataset).ThenBy(l => l.Cell_segment).ThenBy(l => l.PatientId).ThenBy(l => l.KPMP_data_integration_term).ThenBy(l => l.Value_type_1st).ThenBy(l => l.Value_type_2nd).ThenByDescending(l => l.Value_1st).ThenByDescending(l => l.Value_2nd).ToArray();
            int firstIndexSameValue = 0;
            int currentOrdinalRank = 0;
            List<float> sameValues_1st2nd_ordinalRanks = new List<float>();
            float final_rank;

            for (int indexKPMP = 0; indexKPMP < data_length; indexKPMP++)
            {
                data_line = this.KPMP_data[indexKPMP];
                if ((data_line.Value_1st == -1)
                    //|| (data_line.Value_2nd == -1)
                    || (data_line.Value_type_1st.Equals(KPMP_value_type_enum.E_m_p_t_y))
                    || (data_line.Value_type_2nd.Equals(KPMP_value_type_enum.E_m_p_t_y)))
                {
                    throw new Exception();
                }
                if ((indexKPMP == 0)
                    || (!data_line.Dataset.Equals(this.KPMP_data[indexKPMP - 1].Dataset))
                    || (!data_line.Cell_segment.Equals(this.KPMP_data[indexKPMP - 1].Cell_segment))
                    || (!data_line.PatientId.Equals(this.KPMP_data[indexKPMP - 1].PatientId))
                    || (!data_line.KPMP_data_integration_term.Equals(this.KPMP_data[indexKPMP - 1].KPMP_data_integration_term))
                    || (!data_line.Value_type_1st.Equals(this.KPMP_data[indexKPMP - 1].Value_type_1st))
                    || (!data_line.Value_type_2nd.Equals(this.KPMP_data[indexKPMP - 1].Value_type_2nd)))
                {
                    currentOrdinalRank = 0;
                }
                if ((indexKPMP == 0)
                    || (!data_line.Dataset.Equals(this.KPMP_data[indexKPMP - 1].Dataset))
                    || (!data_line.Cell_segment.Equals(this.KPMP_data[indexKPMP - 1].Cell_segment))
                    || (!data_line.PatientId.Equals(this.KPMP_data[indexKPMP - 1].PatientId))
                    || (!data_line.KPMP_data_integration_term.Equals(this.KPMP_data[indexKPMP - 1].KPMP_data_integration_term))
                    || (!data_line.Value_type_1st.Equals(this.KPMP_data[indexKPMP - 1].Value_type_1st))
                    || (!data_line.Value_type_2nd.Equals(this.KPMP_data[indexKPMP - 1].Value_type_2nd))
                    || (!data_line.Value_1st.Equals(this.KPMP_data[indexKPMP - 1].Value_1st))
                    || (!data_line.Value_2nd.Equals(this.KPMP_data[indexKPMP - 1].Value_2nd)))
                {
                    firstIndexSameValue = indexKPMP;
                    sameValues_1st2nd_ordinalRanks.Clear();
                }
                currentOrdinalRank++;
                sameValues_1st2nd_ordinalRanks.Add(currentOrdinalRank);
                if ((indexKPMP == data_length - 1)
                    || (!data_line.Dataset.Equals(this.KPMP_data[indexKPMP + 1].Dataset))
                    || (!data_line.Cell_segment.Equals(this.KPMP_data[indexKPMP + 1].Cell_segment))
                    || (!data_line.PatientId.Equals(this.KPMP_data[indexKPMP + 1].PatientId))
                    || (!data_line.KPMP_data_integration_term.Equals(this.KPMP_data[indexKPMP + 1].KPMP_data_integration_term))
                    || (!data_line.Value_type_1st.Equals(this.KPMP_data[indexKPMP + 1].Value_type_1st))
                    || (!data_line.Value_type_2nd.Equals(this.KPMP_data[indexKPMP + 1].Value_type_2nd))
                    || (!data_line.Value_1st.Equals(this.KPMP_data[indexKPMP + 1].Value_1st))
                    || (!data_line.Value_2nd.Equals(this.KPMP_data[indexKPMP + 1].Value_2nd)))
                {
                    if (sameValues_1st2nd_ordinalRanks.Count > 1)
                    {
                        final_rank = Math_class.Get_average(sameValues_1st2nd_ordinalRanks.ToArray());
                        for (int indexInner = firstIndexSameValue; indexInner <= indexKPMP; indexInner++)
                        {
                            inner_data_line = this.KPMP_data[indexInner];
                            inner_data_line.Fractional_rank_based_on_value_types = final_rank;
                        }
                    }
                    else if (sameValues_1st2nd_ordinalRanks.Count == 1)
                    {
                        if (firstIndexSameValue != indexKPMP) { throw new Exception(); }
                        data_line.Fractional_rank_based_on_value_types = sameValues_1st2nd_ordinalRanks[0];
                    }
                    else { throw new Exception(); }
                }
            }
        }

        public void Set_valueType_1st_as_score_of_interest()
        {
            foreach (KPMP_standardized_dataset_line_class data_line in this.KPMP_data)
            {
                data_line.Score_of_interest = data_line.Value_1st;
                data_line.Score_of_interest_type = Statistical_score_of_interest_type_enum.Value;
            }
        }

        public void Calculate_log2_ratio_values_1st_and_add_as_new_value_1st()
        {
            List<KPMP_standardized_dataset_line_class> add = new List<KPMP_standardized_dataset_line_class>();
            KPMP_standardized_dataset_line_class add_line;
            foreach (KPMP_standardized_dataset_line_class data_line in this.KPMP_data)
            {
                if (!data_line.Value_type_2nd.Equals(KPMP_value_type_enum.No_selection)) { throw new Exception(); }
                switch (data_line.Value_type_1st)
                {
                    case KPMP_value_type_enum.Ratio_singlepatient:
                        add_line = data_line.Deep_copy();
                        add_line.Value_1st = (float)Math.Log(add_line.Value_1st, 2);
                        if (  (double.IsNaN(add_line.Value_1st))
                            ||(double.IsInfinity(add_line.Value_1st)))
                        {
                            throw new Exception();
                        }
                        switch (data_line.Value_type_1st)
                        {
                            case KPMP_value_type_enum.Ratio_singlepatient:
                                add_line.Value_type_1st = KPMP_value_type_enum.Log2_ratio_singlepatient;
                                break;
                            case KPMP_value_type_enum.Single_value:
                                add_line.Value_type_1st = KPMP_value_type_enum.Log2_single_value;
                                break;
                            default:
                                throw new Exception();
                        }
                        add_line.Value_2nd = 0;
                        add_line.Value_type_2nd = KPMP_value_type_enum.No_selection;
                        add.Add(add_line);
                        break;
                    case KPMP_value_type_enum.Single_value:
                        break;
                    default:
                        throw new Exception();
                }
            }
            Add_to_array(add.ToArray());
        }

        public void Add_integration_terms()
        {
            foreach (KPMP_standardized_dataset_line_class standardized_data_line in this.KPMP_data)
            {
                standardized_data_line.KPMP_data_integration_term = KPMP_data_integration_class.Get_kpmp_integration_term(standardized_data_line.Cell_segment);
            }
        }

        public void Reset_integrationTerm_for_indicated_center_cellTypes(Dictionary<string, Dictionary<string, string>> center_cellSubType_setIntegrationTerm_dict)
        {
            int kpmp_data_length = this.KPMP_data.Length;
            KPMP_standardized_dataset_line_class standardized_line;
            Dictionary<string, string> cellSubType_setIntegrationTerm_dict;
            for (int indexData =0; indexData<kpmp_data_length;indexData++)
            {
                standardized_line = this.KPMP_data[indexData];
                if (center_cellSubType_setIntegrationTerm_dict.ContainsKey(standardized_line.Dataset))
                {
                    cellSubType_setIntegrationTerm_dict = center_cellSubType_setIntegrationTerm_dict[standardized_line.Dataset];
                    if (cellSubType_setIntegrationTerm_dict.ContainsKey(standardized_line.Cell_segment))
                    {
                        standardized_line.KPMP_data_integration_term = (string)cellSubType_setIntegrationTerm_dict[standardized_line.Cell_segment].Clone();
                    }
                }
            }
        }

        #region Keep
        public void Keep_genes_with_indicated_valueTypes1st2nd_if_below_maxFractionalRank_and_check_if_at_least_one_kept_value1stor2nd_is_not_infinity_keep_all_other_genes(float max_fractional_rank, KPMP_value_type_enum value_type_1st, KPMP_value_type_enum value_type_2nd)
        {
            Calculate_fractional_ranks_based_on_descending_values1st_and_values2nd();
            int data_length = this.KPMP_data.Length;
            KPMP_standardized_dataset_line_class kpmp_standardized_data_line;
            Dictionary<string, Dictionary<string, Dictionary<string, bool>>> keep_dataset_sampleName_gene_dict = new Dictionary<string, Dictionary<string, Dictionary<string, bool>>>();
            Dictionary<string, Dictionary<string, bool>> dataset_sampleName_with_atLeastOneNonInfinityValue1stOr2nd_dict = new Dictionary<string, Dictionary<string, bool>>();
            List<KPMP_standardized_dataset_line_class> keep = new List<KPMP_standardized_dataset_line_class>();
            List<KPMP_standardized_dataset_line_class> always_keep = new List<KPMP_standardized_dataset_line_class>();
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                kpmp_standardized_data_line = this.KPMP_data[indexData];
                if (  (kpmp_standardized_data_line.Value_type_1st.Equals(value_type_1st))
                    &&(kpmp_standardized_data_line.Value_type_2nd.Equals(value_type_2nd)))
                {
                    if (kpmp_standardized_data_line.Fractional_rank_based_on_value_types <= max_fractional_rank)
                    {
                        keep.Add(kpmp_standardized_data_line);
                        if ((!double.IsInfinity(kpmp_standardized_data_line.Value_1st))
                            || (!double.IsInfinity(kpmp_standardized_data_line.Value_2nd)))
                        {
                            if (!dataset_sampleName_with_atLeastOneNonInfinityValue1stOr2nd_dict.ContainsKey(kpmp_standardized_data_line.Dataset))
                            {
                                dataset_sampleName_with_atLeastOneNonInfinityValue1stOr2nd_dict.Add(kpmp_standardized_data_line.Dataset, new Dictionary<string, bool>());
                            }
                            if (!dataset_sampleName_with_atLeastOneNonInfinityValue1stOr2nd_dict[kpmp_standardized_data_line.Dataset].ContainsKey(kpmp_standardized_data_line.Cell_segment))
                            {
                                dataset_sampleName_with_atLeastOneNonInfinityValue1stOr2nd_dict[kpmp_standardized_data_line.Dataset].Add(kpmp_standardized_data_line.Cell_segment, true);
                            }
                        }
                    }
                }
                else
                {
                    always_keep.Add(kpmp_standardized_data_line);
                }
            }
            this.KPMP_data = keep.ToArray();
            data_length = this.KPMP_data.Length;
            keep.Clear();
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                kpmp_standardized_data_line = this.KPMP_data[indexData];
                if (dataset_sampleName_with_atLeastOneNonInfinityValue1stOr2nd_dict[kpmp_standardized_data_line.Dataset].ContainsKey(kpmp_standardized_data_line.Cell_segment))
                {
                    keep.Add(kpmp_standardized_data_line);
                }
                else
                {
                    throw new Exception();
                }
            }
            keep.AddRange(always_keep);
            this.KPMP_data = keep.ToArray();
            this.Data_filtered_by_cutoff = true;
        }

        public void Keep_only_indicated_integration_types_and_datasetPatients_with_all_terms(string[] keep_integrations)
        {
            keep_integrations = keep_integrations.OrderBy(l => l).ToArray();
            string keep_integration;
            int keep_integrations_length = keep_integrations.Length;
            int indexKeep = 0;
            this.KPMP_data = this.KPMP_data.OrderBy(l => l.Dataset).ThenBy(l => l.PatientId).ThenBy(l => l.KPMP_data_integration_term).ToArray();
            int kpmp_data_length = this.KPMP_data.Length;
            KPMP_standardized_dataset_line_class kpmp_data_line;
            List<KPMP_standardized_dataset_line_class> final_keep = new List<KPMP_standardized_dataset_line_class>();
            List<KPMP_standardized_dataset_line_class> current_keep = new List<KPMP_standardized_dataset_line_class>();
            int stringCompare = -2;
            bool[] integration_type_exists = new bool[keep_integrations_length];
            bool all_integration_types_exist;
            for (int indexKPMP = 0; indexKPMP < kpmp_data_length; indexKPMP++)
            {
                kpmp_data_line = this.KPMP_data[indexKPMP];
                if ((indexKPMP == 0)
                    || (!kpmp_data_line.Dataset.Equals(this.KPMP_data[indexKPMP - 1].Dataset))
                    || (!kpmp_data_line.PatientId.Equals(this.KPMP_data[indexKPMP - 1].PatientId)))
                {
                    indexKeep = 0;
                    current_keep.Clear();
                    for (int indexInt = 0; indexInt < keep_integrations_length; indexInt++)
                    {
                        integration_type_exists[indexInt] = false;
                    }
                }
                stringCompare = -2;
                while ((indexKeep < keep_integrations_length) && (stringCompare < 0))
                {
                    keep_integration = keep_integrations[indexKeep];
                    stringCompare = keep_integration.CompareTo(kpmp_data_line.KPMP_data_integration_term);
                    if (stringCompare < 0)
                    {
                        indexKeep++;
                    }
                }
                if (stringCompare == 0)
                {
                    integration_type_exists[indexKeep] = true;
                    current_keep.Add(kpmp_data_line);
                }
                if ((indexKPMP == kpmp_data_length - 1)
                    || (!kpmp_data_line.Dataset.Equals(this.KPMP_data[indexKPMP + 1].Dataset))
                    || (!kpmp_data_line.PatientId.Equals(this.KPMP_data[indexKPMP + 1].PatientId)))
                {
                    all_integration_types_exist = true;
                    for (int indexCheck = 0; indexCheck < keep_integrations_length; indexCheck++)
                    {
                        if (!integration_type_exists[indexCheck])
                        {
                            all_integration_types_exist = false;
                        }
                    }
                    if (all_integration_types_exist)
                    {
                        final_keep.AddRange(current_keep);
                    }
                }
            }
            this.KPMP_data = final_keep.ToArray();
        }

        public void Keep_only_indicated_gene_symbols(string[] keep_gene_symbols)
        {
            Dictionary<string, bool> keep_gene_symbol_dict = new Dictionary<string, bool>();
            foreach (string keep_gene_symbol in keep_gene_symbols)
            {
                keep_gene_symbol_dict.Add(keep_gene_symbol.ToUpper(),true);
            }

            int kpmp_data_length = this.KPMP_data.Length;
            KPMP_standardized_dataset_line_class kpmp_data_line;
            List<KPMP_standardized_dataset_line_class> final_keep = new List<KPMP_standardized_dataset_line_class>();
            for (int indexKPMP = 0; indexKPMP < kpmp_data_length; indexKPMP++)
            {
                kpmp_data_line = this.KPMP_data[indexKPMP];
                if (keep_gene_symbol_dict.ContainsKey(kpmp_data_line.Gene_symbol.ToUpper()))
                {
                    final_keep.Add(kpmp_data_line);
                }
            }
            this.KPMP_data = final_keep.ToArray();
        }

        public void Keep_only_indicated_datasets(params string[] keep_datasets)
        {
            int keep_integrations_length = keep_datasets.Length;
            int kpmp_data_length = this.KPMP_data.Length;
            KPMP_standardized_dataset_line_class kpmp_data_line;
            List<KPMP_standardized_dataset_line_class> final_keep = new List<KPMP_standardized_dataset_line_class>();
            for (int indexKPMP = 0; indexKPMP < kpmp_data_length; indexKPMP++)
            {
                kpmp_data_line = this.KPMP_data[indexKPMP];
                if (keep_datasets.Contains(kpmp_data_line.Dataset))
                {
                    final_keep.Add(kpmp_data_line);
                }
            }
            this.KPMP_data = final_keep.ToArray();
        }

        public void Keep_only_indicated_valueTypes_1st(params KPMP_value_type_enum[] value_types_1st)
        {
            int kpmp_data_length = this.KPMP_data.Length;
            KPMP_standardized_dataset_line_class kpmp_data_line;
            List<KPMP_standardized_dataset_line_class> final_keep = new List<KPMP_standardized_dataset_line_class>();
            for (int indexKPMP = 0; indexKPMP < kpmp_data_length; indexKPMP++)
            {
                kpmp_data_line = this.KPMP_data[indexKPMP];
                if (value_types_1st.Contains(kpmp_data_line.Value_type_1st))
                {
                    final_keep.Add(kpmp_data_line);
                }
            }
            this.KPMP_data = final_keep.ToArray();
        }

        public void Keep_only_lines_with_1stValue_aboveEqual_cutoff_and_check_if_1st_value_is_among_input_values(double cutoff, params KPMP_value_type_enum[] mandatory_first_values)
        {
            List<KPMP_standardized_dataset_line_class> keep = new List<KPMP_standardized_dataset_line_class>();
            foreach (KPMP_standardized_dataset_line_class dataset_line in this.KPMP_data)
            {
                if (!mandatory_first_values.Contains(dataset_line.Value_type_1st)) { throw new Exception(); }
                if (dataset_line.Value_1st>=cutoff)
                {
                    keep.Add(dataset_line);
                }
            }
            this.KPMP_data = keep.ToArray();
            this.Data_filtered_by_cutoff = true;
        }
        #endregion

        #region Get
        public string[] Get_all_genes_that_have_a_non_zero_value_of_indicated_value_type_1st_in_all_datasets(KPMP_value_type_enum value_type_1st)
        {
            string[] all_datasets = Get_all_unique_ordered_datasets();
            int all_datasets_length = all_datasets.Length;
            string current_dataset;
            int indexDataSet = 0;
            bool gene_in_all_datasets = false;
            this.KPMP_data = this.KPMP_data.OrderBy(l => l.Value_type_1st).ThenBy(l => l.Gene_symbol).ThenBy(l => l.Dataset).ToArray();
            int data_length = this.KPMP_data.Length;
            KPMP_standardized_dataset_line_class standardized_line;
            List<string> genes = new List<string>();
            int stringCompare = -2;
            List<string> genes_not_in_dataset = new List<string>();
            bool[] gene_in_dataset = new bool[all_datasets_length];
            for (int indexD = 0; indexD < data_length; indexD++)
            {
                standardized_line = this.KPMP_data[indexD];
                if (standardized_line.Value_type_1st.Equals(value_type_1st))
                {
                    if ((indexD == 0)
                        || (!standardized_line.Value_type_1st.Equals(this.KPMP_data[indexD - 1].Value_type_1st))
                        || (!standardized_line.Gene_symbol.Equals(this.KPMP_data[indexD - 1].Gene_symbol)))
                    {
                        indexDataSet = 0;
                        for (int indexDatasetLoop=0; indexDatasetLoop < all_datasets_length; indexDatasetLoop++)
                        {
                            gene_in_dataset[indexDatasetLoop] = false;
                        }
                    }
                    current_dataset = all_datasets[indexDataSet];
                    stringCompare = current_dataset.CompareTo(standardized_line.Dataset);
                    if (stringCompare < 0) { indexDataSet++; }
                    else if (stringCompare==0)
                    {
                        if (standardized_line.Value_1st != 0)
                        {
                            gene_in_dataset[indexDataSet] = true;
                        }
                    }
                    if ((indexD == data_length - 1)
                        || (!standardized_line.Value_type_1st.Equals(this.KPMP_data[indexD + 1].Value_type_1st))
                        || (!standardized_line.Gene_symbol.Equals(this.KPMP_data[indexD + 1].Gene_symbol)))
                    {
                        gene_in_all_datasets = true;
                        for (int indexDatasetLoop = 0; indexDatasetLoop < all_datasets_length; indexDatasetLoop++)
                        {
                            if (!gene_in_dataset[indexDatasetLoop])
                            {
                                gene_in_all_datasets = false;
                                break;
                            }
                        }
                        if (gene_in_all_datasets)
                        {
                            genes.Add(standardized_line.Gene_symbol);
                        }
                        else
                        {
                            genes_not_in_dataset.Add(standardized_line.Gene_symbol);
                        }
                    }
                }
            }
            string[] overlap = Overlap_class.Get_intersection(genes.ToArray(), genes_not_in_dataset.ToArray());
            if (overlap.Length!=0) { throw new Exception(); }
            return genes.ToArray();
        }

        public string[] Get_all_unique_ordered_datasets()
        {
            this.KPMP_data = this.KPMP_data.OrderBy(l => l.Dataset).ToArray();
            List<string> unique_datasets = new List<string>();
            int kpmp_length = this.KPMP_data.Length;
            KPMP_standardized_dataset_line_class dataset_line;
            for (int indexKPMP = 0; indexKPMP < kpmp_length; indexKPMP++)
            {
                dataset_line = this.KPMP_data[indexKPMP];
                if ((indexKPMP == 0)
                    || (!dataset_line.Dataset.Equals(this.KPMP_data[indexKPMP - 1].Dataset)))
                {
                    unique_datasets.Add(dataset_line.Dataset);
                }
            }
            return unique_datasets.ToArray();
        }

        public string[] Get_all_unique_ordered_genes()
        {
            this.KPMP_data = this.KPMP_data.OrderBy(l => l.Gene_symbol).ToArray();
            List<string> unique_genes = new List<string>();
            int kpmp_length = this.KPMP_data.Length;
            KPMP_standardized_dataset_line_class dataset_line;
            for (int indexKPMP=0;indexKPMP<kpmp_length;indexKPMP++)
            {
                dataset_line = this.KPMP_data[indexKPMP];
                if ((indexKPMP==0)
                    ||(!dataset_line.Gene_symbol.Equals(this.KPMP_data[indexKPMP-1].Gene_symbol)))
                {
                    unique_genes.Add(dataset_line.Gene_symbol);
                }
            }
            return unique_genes.ToArray();
        }
        #endregion

        public KPMP_standardized_dataset_class Deep_copy()
        {
            KPMP_standardized_dataset_class copy = (KPMP_standardized_dataset_class)this.MemberwiseClone();
            int data_length = this.KPMP_data.Length;
            copy.KPMP_data = new KPMP_standardized_dataset_line_class[data_length];
            for (int indexKPMP=0; indexKPMP<data_length;indexKPMP++)
            {
                copy.KPMP_data[indexKPMP] = this.KPMP_data[indexKPMP].Deep_copy();
            }
            string[] datasets = this.Dataset_bgGenesProteins_dict.Keys.ToArray();
            string dataset;
            int datasets_length = datasets.Length;
            copy.Dataset_bgGenesProteins_dict = new Dictionary<string, string[]>();
            for (int indexD=0;indexD<datasets_length;indexD++)
            {
                dataset = datasets[indexD];
                copy.Dataset_bgGenesProteins_dict.Add(dataset, Array_class.Deep_copy_string_array(this.Dataset_bgGenesProteins_dict[dataset]));
            }
            return copy;
        }

        private void Read_and_fill_bgGenesProteins_dictionary(string directory, string fileName)
        {
            string bgSet_fileName = Path.GetFileNameWithoutExtension(fileName) + "_bgGenes.txt";
            KPMP_standardized_dataset_bgGenesProteins_readWriteOptions_class bg_readWriteOptions = new KPMP_standardized_dataset_bgGenesProteins_readWriteOptions_class(directory, bgSet_fileName);
            KPMP_standardized_bgGeneProtein_line_class[] bgGenesProteins = ReadWriteClass.ReadRawData_and_FillArray<KPMP_standardized_bgGeneProtein_line_class>(bg_readWriteOptions);
            Dataset_bgGenesProteins_dict = new Dictionary<string, string[]>();
            bgGenesProteins = bgGenesProteins.OrderBy(l => l.Dataset).ThenBy(l => l.BgGeneProtein).ToArray();
            int bg_length = bgGenesProteins.Length;
            KPMP_standardized_bgGeneProtein_line_class bgGenesProtein_line;
            List<string> currentDataset_genes = new List<string>();
            for (int indexBg = 0; indexBg < bg_length; indexBg++)
            {
                bgGenesProtein_line = bgGenesProteins[indexBg];
                if ((indexBg == 0) || (!bgGenesProtein_line.Dataset.Equals(bgGenesProteins[indexBg - 1].Dataset)))
                {
                    currentDataset_genes.Clear();
                }
                currentDataset_genes.Add(bgGenesProtein_line.BgGeneProtein);
                if ((indexBg == bg_length - 1) || (!bgGenesProtein_line.Dataset.Equals(bgGenesProteins[indexBg + 1].Dataset)))
                {
                    Dataset_bgGenesProteins_dict.Add(bgGenesProtein_line.Dataset, currentDataset_genes.ToArray());
                }
            }
        }

        private void Write_bgGenesProteins_of_selected_dataset(string directory, string dataset, string fileName)
        {
            string[] genes = Dataset_bgGenesProteins_dict[dataset];
            string bgSet_fileName = Path.GetFileNameWithoutExtension(fileName) + "_bgGenes" + Path.GetExtension(fileName);
            ReadWriteClass.WriteArray<string>(genes, directory + bgSet_fileName);
        }

        private void Write_bgGenesProteins_as_one_file(string subdirectory, string fileName)
        {
            string[] datasets = Dataset_bgGenesProteins_dict.Keys.ToArray();
            string dataset;
            int datasets_length = datasets.Length;
            KPMP_standardized_bgGeneProtein_line_class bgGene_line;
            List<KPMP_standardized_bgGeneProtein_line_class> bgGene_lines = new List<KPMP_standardized_bgGeneProtein_line_class>();
            string[] bgGenes;
            for (int indexD = 0; indexD < datasets_length; indexD++)
            {
                dataset = datasets[indexD];
                bgGenes = Dataset_bgGenesProteins_dict[dataset];
                foreach (string bgGene in bgGenes)
                {
                    bgGene_line = new KPMP_standardized_bgGeneProtein_line_class();
                    bgGene_line.BgGeneProtein = (string)bgGene.Clone();
                    bgGene_line.Dataset = (string)dataset.Clone();
                    bgGene_lines.Add(bgGene_line);
                }
            }
            KPMP_standardized_dataset_bgGenesProteins_readWriteOptions_class readWriteOptions = new KPMP_standardized_dataset_bgGenesProteins_readWriteOptions_class(subdirectory, fileName);
            ReadWriteClass.WriteData(bgGene_lines, readWriteOptions);
        }

        public void Read(string subdirectory, string fileName, bool readBg_proteins)
        {
            KPMP_standardized_dataset_readWriteOptions_class readWriteOptions = new KPMP_standardized_dataset_readWriteOptions_class(subdirectory, fileName);
            this.KPMP_data = ReadWriteClass.ReadRawData_and_FillArray<KPMP_standardized_dataset_line_class>(readWriteOptions);
            if (readBg_proteins) { Read_and_fill_bgGenesProteins_dictionary(subdirectory, fileName); }
        }

        public void Write(string subdirectory, string fileName, bool write_bgProteins)
        {
            KPMP_standardized_dataset_readWriteOptions_class readWriteOptions = new KPMP_standardized_dataset_readWriteOptions_class(subdirectory, fileName);
            ReadWriteClass.WriteData(this.KPMP_data, readWriteOptions);
            if (write_bgProteins) { Write_bgGenesProteins_as_one_file(subdirectory, Path.GetFileNameWithoutExtension(fileName) + "_bgGenes.txt"); }
        }

        public void Write_one_file_foreach_dataset(string directory, string fileName_end, bool write_bgProteins)
        {
            int kpmp_data_length = this.KPMP_data.Length;
            KPMP_standardized_dataset_line_class data_line;
            List<KPMP_standardized_dataset_line_class> dataset_kpmp_data = new List<KPMP_standardized_dataset_line_class>();
            this.KPMP_data = this.KPMP_data.OrderBy(l => l.Dataset).ToArray();
            string current_dataset_fileName;
            string current_fileName_end;
            KPMP_value_type_enum value_type = KPMP_value_type_enum.E_m_p_t_y;
            for (int indexKPMP = 0; indexKPMP < kpmp_data_length; indexKPMP++)
            {
                data_line = this.KPMP_data[indexKPMP];
                if ((indexKPMP == 0) || (!data_line.Dataset.Equals(this.KPMP_data[indexKPMP - 1].Dataset)))
                {
                    dataset_kpmp_data.Clear();
                    value_type = data_line.Value_type_1st;
                }
                if (!data_line.Value_type_1st.Equals(value_type)) { throw new Exception(); }
                dataset_kpmp_data.Add(data_line);
                if ((indexKPMP == kpmp_data_length - 1) || (!data_line.Dataset.Equals(this.KPMP_data[indexKPMP + 1].Dataset)))
                {
                    current_fileName_end = "";
                    switch (value_type)
                    {
                        case KPMP_value_type_enum.Minus_log10_pvalue:
                            current_fileName_end = "pvalue" + fileName_end;
                            break;
                        case KPMP_value_type_enum.Minus_log10_pvalue_adjusted:
                            current_fileName_end = "adjPvalue" + fileName_end;
                            break;
                        default:
                            throw new Exception();
                    }
                    current_dataset_fileName = data_line.Dataset + "_" + current_fileName_end;
                    KPMP_standardized_dataset_readWriteOptions_class readWriteOptions = new KPMP_standardized_dataset_readWriteOptions_class(directory, current_dataset_fileName);
                    ReadWriteClass.WriteData(dataset_kpmp_data.ToArray(), readWriteOptions);
                    if (write_bgProteins)
                    { Write_bgGenesProteins_of_selected_dataset(directory, data_line.Dataset, current_dataset_fileName); }
                }
            }
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class KPMP_standardized_dataset_collapseOnMethod_line_class
    {
        public string Gene_symbol { get; set; }
        public string Method { get; set; }
        public string KPMP_data_integration_term { get; set; }
        public double ScoreOfInterest_median { get; set; }
        public double ScoreOfInterest_max { get; set; }
        public double ScoreOfInterest_min { get; set; }
        public double ScoreOfInterest_average { get; set; }
        public double ScoreOfInterest_populationSD { get; set; }
        public double ScoreOfInterest_populationCV { get; set; }
        public KPMP_value_type_enum Value_type { get; set; }
        public Statistical_score_of_interest_type_enum Score_of_interest_type { get; set; }

        public KPMP_standardized_dataset_collapseOnMethod_line_class Deep_copy()
        {
            KPMP_standardized_dataset_collapseOnMethod_line_class copy = (KPMP_standardized_dataset_collapseOnMethod_line_class)this.MemberwiseClone();
            copy.Gene_symbol = (string)this.Gene_symbol.Clone();
            copy.Method = (string)this.Method.Clone();
            copy.KPMP_data_integration_term = (string)this.KPMP_data_integration_term.Clone();
            return copy;
        }
    }

    class KPMP_standardized_dataset_collapsOnMethod_readWriteOptions : ReadWriteOptions_base
    {
        public KPMP_standardized_dataset_collapsOnMethod_readWriteOptions(string directory, string fileName)
        {
            ReadWriteClass.Create_directory_if_it_does_not_exist(directory);
            this.File = directory + fileName;
            this.Key_propertyNames = new string[] { "Gene_symbol", "KPMP_data_integration_term", "Method", "ScoreOfInterest_median", "ScoreOfInterest_max", "ScoreOfInterest_min", "ScoreOfInterest_average", "ScoreOfInterest_populationSD", "ScoreOfInterest_populationCV", "Score_of_interest_type" };
            this.Key_columnNames = this.Key_propertyNames;
            this.HeadlineDelimiters = new char[] { Global_class.Tab };
            this.LineDelimiters = new char[] { Global_class.Tab };
            this.File_has_headline = true;
            this.Report = ReadWrite_report_enum.Report_main;
        }
    }

    class KPMP_standardized_dataset_collapseOnMethod_class
    {
        public KPMP_standardized_dataset_collapseOnMethod_line_class[] Collapse { get; set; }

        public KPMP_standardized_dataset_collapseOnMethod_class()
        {
            this.Collapse = new KPMP_standardized_dataset_collapseOnMethod_line_class[0];
        }

        private void Add_to_array(KPMP_standardized_dataset_collapseOnMethod_line_class[] add_collapse)
        {
            int this_length = this.Collapse.Length;
            int add_length = add_collapse.Length;
            int new_length = this_length + add_length;
            int indexNew = -1;
            KPMP_standardized_dataset_collapseOnMethod_line_class[] new_collapse = new KPMP_standardized_dataset_collapseOnMethod_line_class[new_length];
            for (int indexThis = 0; indexThis < this_length; indexThis++)
            {
                indexNew++;
                new_collapse[indexNew] = this.Collapse[indexThis];
            }
            for (int indexAdd = 0; indexAdd < add_length; indexAdd++)
            {
                indexNew++;
                new_collapse[indexNew] = add_collapse[indexAdd];
            }
            this.Collapse = new_collapse;
        }

        public void Generate_from_kpmp_standardized_dataset(KPMP_standardized_dataset_class standardized_dataset)
        {
            standardized_dataset.KPMP_data = standardized_dataset.KPMP_data.OrderBy(l => l.Value_type_1st).ThenBy(l=>l.Score_of_interest_type).ThenBy(l => l.KPMP_data_integration_term).ThenBy(l => l.Gene_symbol).ThenBy(l => l.Dataset).ThenBy(l => l.PatientId).ThenBy(l=>l.Cell_segment).ToArray();
            int standardized_kpmp_data_length = standardized_dataset.KPMP_data.Length;
            KPMP_standardized_dataset_line_class standardized_line;
            List<double> currentITGeneDatasetPatient_allCellSegmentScoreOfInterest = new List<double>();
            List<double> currentITGeneDataset_allPatientScoreOfInterest_basedOnAverage = new List<double>();
            List<double> currentITGene_allProteomics_allDatasetScoreOfInterest_basedOnAverage = new List<double>();
            List<double> currentITGene_allRNASeq_allDatasetScoreOfInterest_basedOnAverage = new List<double>();
            List<double> currentITGeneDataset_allPatientScoreOfInterest_basedOnMedian = new List<double>();
            List<double> currentITGene_allProteomics_allDatasetScoreOfInterest_basedOnMedian = new List<double>();
            List<double> currentITGene_allRNASeq_allDatasetScoreOfInterest_basedOnMedian = new List<double>();
            double currentITGeneDataset_averageAllPatientScoreOfInterest;
            double currentITGeneDatasetPatient_averageAllCellSegmentScoreOfInterest;
            double currentITGeneDataset_medianAllPatientScoreOfInterest;
            double currentITGeneDatasetPatient_medianAllCellSegmentScoreOfInterest;
            double currentITGeneDataset_populationSDAllPatientScoreOfInterest;
            double currentITGeneDataset_maxAllPatientScoreOfInterest;
            double currentITGeneDataset_minAllPatientScoreOfInterest;
            double averageScoreOfInterest;
            double popSDScoreOfInterest;
            double medianScoreOfInterest;
            double maxScoreOfInterest;
            double minScoreOfInterest;
            Dictionary<string,Dictionary<string,Dictionary<string, int>>> dataset_patientID_integrationTerm_count_dict = new Dictionary<string, Dictionary<string, Dictionary<string, int>>>();
            Dictionary<string,int> dataset_patientIDCount_dict = new Dictionary<string, int>();

            KPMP_standardized_dataset_collapseOnMethod_line_class standardized_dataset_collapseMethod_line;
            List<KPMP_standardized_dataset_collapseOnMethod_line_class> collapse_list = new List<KPMP_standardized_dataset_collapseOnMethod_line_class>();

            for (int indexData=0; indexData<standardized_kpmp_data_length;indexData++)
            {
                standardized_line = standardized_dataset.KPMP_data[indexData];
                if (!standardized_line.Value_type_2nd.Equals(KPMP_value_type_enum.No_selection)) { throw new Exception(); }
                if ((indexData == 0)
                    || (!standardized_line.Score_of_interest_type.Equals(standardized_dataset.KPMP_data[indexData - 1].Score_of_interest_type))
                    || (!standardized_line.Value_type_1st.Equals(standardized_dataset.KPMP_data[indexData - 1].Value_type_1st))
                    || (!standardized_line.KPMP_data_integration_term.Equals(standardized_dataset.KPMP_data[indexData - 1].KPMP_data_integration_term))
                    || (!standardized_line.Gene_symbol.Equals(standardized_dataset.KPMP_data[indexData - 1].Gene_symbol)))
                {
                    currentITGene_allProteomics_allDatasetScoreOfInterest_basedOnAverage.Clear();
                    currentITGene_allRNASeq_allDatasetScoreOfInterest_basedOnAverage.Clear();
                    currentITGene_allProteomics_allDatasetScoreOfInterest_basedOnMedian.Clear();
                    currentITGene_allRNASeq_allDatasetScoreOfInterest_basedOnMedian.Clear();
                }

                if ((indexData == 0)
                    || (!standardized_line.Score_of_interest_type.Equals(standardized_dataset.KPMP_data[indexData - 1].Score_of_interest_type))
                    || (!standardized_line.Value_type_1st.Equals(standardized_dataset.KPMP_data[indexData - 1].Value_type_1st))
                    || (!standardized_line.KPMP_data_integration_term.Equals(standardized_dataset.KPMP_data[indexData - 1].KPMP_data_integration_term))
                    || (!standardized_line.Gene_symbol.Equals(standardized_dataset.KPMP_data[indexData - 1].Gene_symbol))
                    || (!standardized_line.Dataset.Equals(standardized_dataset.KPMP_data[indexData - 1].Dataset)))
                {
                    currentITGeneDataset_allPatientScoreOfInterest_basedOnAverage.Clear();
                    currentITGeneDataset_allPatientScoreOfInterest_basedOnMedian.Clear();
                }

                if ((indexData == 0)
                    || (!standardized_line.Score_of_interest_type.Equals(standardized_dataset.KPMP_data[indexData - 1].Score_of_interest_type))
                    || (!standardized_line.Value_type_1st.Equals(standardized_dataset.KPMP_data[indexData - 1].Value_type_1st))
                    || (!standardized_line.KPMP_data_integration_term.Equals(standardized_dataset.KPMP_data[indexData - 1].KPMP_data_integration_term))
                    || (!standardized_line.Gene_symbol.Equals(standardized_dataset.KPMP_data[indexData - 1].Gene_symbol))
                    || (!standardized_line.Dataset.Equals(standardized_dataset.KPMP_data[indexData - 1].Dataset))
                    || (!standardized_line.PatientId.Equals(standardized_dataset.KPMP_data[indexData - 1].PatientId)))
                {
                    currentITGeneDatasetPatient_allCellSegmentScoreOfInterest.Clear();
                }

                currentITGeneDatasetPatient_allCellSegmentScoreOfInterest.Add(standardized_line.Score_of_interest);

                if ((indexData == standardized_kpmp_data_length - 1)
                    || (!standardized_line.Score_of_interest_type.Equals(standardized_dataset.KPMP_data[indexData + 1].Score_of_interest_type))
                    || (!standardized_line.Value_type_1st.Equals(standardized_dataset.KPMP_data[indexData + 1].Value_type_1st))
                    || (!standardized_line.KPMP_data_integration_term.Equals(standardized_dataset.KPMP_data[indexData + 1].KPMP_data_integration_term))
                    || (!standardized_line.Gene_symbol.Equals(standardized_dataset.KPMP_data[indexData + 1].Gene_symbol))
                    || (!standardized_line.Dataset.Equals(standardized_dataset.KPMP_data[indexData + 1].Dataset))
                    || (!standardized_line.PatientId.Equals(standardized_dataset.KPMP_data[indexData + 1].PatientId)))
                {
                    currentITGeneDatasetPatient_averageAllCellSegmentScoreOfInterest = Math_class.Get_average(currentITGeneDatasetPatient_allCellSegmentScoreOfInterest.ToArray());
                    currentITGeneDataset_allPatientScoreOfInterest_basedOnAverage.Add(currentITGeneDatasetPatient_averageAllCellSegmentScoreOfInterest);
                    currentITGeneDatasetPatient_medianAllCellSegmentScoreOfInterest = Math_class.Get_median(currentITGeneDatasetPatient_allCellSegmentScoreOfInterest.ToArray());
                    currentITGeneDataset_allPatientScoreOfInterest_basedOnMedian.Add(currentITGeneDatasetPatient_medianAllCellSegmentScoreOfInterest);
                    if (!dataset_patientID_integrationTerm_count_dict.ContainsKey(standardized_line.Dataset))
                    {
                        dataset_patientID_integrationTerm_count_dict.Add(standardized_line.Dataset, new Dictionary<string, Dictionary<string, int>>());
                    }
                    if (!dataset_patientID_integrationTerm_count_dict[standardized_line.Dataset].ContainsKey(standardized_line.PatientId))
                    {
                        dataset_patientID_integrationTerm_count_dict[standardized_line.Dataset].Add(standardized_line.PatientId, new Dictionary<string, int>());
                    }
                    if (!dataset_patientID_integrationTerm_count_dict[standardized_line.Dataset][standardized_line.PatientId].ContainsKey(standardized_line.KPMP_data_integration_term))
                    {
                        dataset_patientID_integrationTerm_count_dict[standardized_line.Dataset][standardized_line.PatientId].Add(standardized_line.KPMP_data_integration_term, currentITGeneDatasetPatient_allCellSegmentScoreOfInterest.Count);
                    }
                    else if (dataset_patientID_integrationTerm_count_dict[standardized_line.Dataset][standardized_line.PatientId][standardized_line.KPMP_data_integration_term] != currentITGeneDatasetPatient_allCellSegmentScoreOfInterest.Count)
                    {
                        throw new Exception();
                    }
                }

                if ((indexData == standardized_kpmp_data_length-1)
                    || (!standardized_line.Score_of_interest_type.Equals(standardized_dataset.KPMP_data[indexData + 1].Score_of_interest_type))
                    || (!standardized_line.Value_type_1st.Equals(standardized_dataset.KPMP_data[indexData + 1].Value_type_1st))
                    || (!standardized_line.KPMP_data_integration_term.Equals(standardized_dataset.KPMP_data[indexData + 1].KPMP_data_integration_term))
                    || (!standardized_line.Gene_symbol.Equals(standardized_dataset.KPMP_data[indexData + 1].Gene_symbol))
                    || (!standardized_line.Dataset.Equals(standardized_dataset.KPMP_data[indexData + 1].Dataset)))
                {
                    Math_class.Get_mean_and_population_sd(currentITGeneDataset_allPatientScoreOfInterest_basedOnAverage.ToArray(), out currentITGeneDataset_averageAllPatientScoreOfInterest, out currentITGeneDataset_populationSDAllPatientScoreOfInterest);
                    Math_class.Get_max_min_of_array(currentITGeneDataset_allPatientScoreOfInterest_basedOnAverage.ToArray(), out currentITGeneDataset_maxAllPatientScoreOfInterest, out currentITGeneDataset_minAllPatientScoreOfInterest);
                    currentITGeneDataset_medianAllPatientScoreOfInterest = Math_class.Get_median(currentITGeneDataset_allPatientScoreOfInterest_basedOnMedian.ToArray());
                    switch (KPMP_dataset_name_class.Get_datasetClass(standardized_line.Dataset))
                    {
                        case KPMP_dataset_name_class.Proteomics_datasetClass:
                            currentITGene_allProteomics_allDatasetScoreOfInterest_basedOnAverage.Add(currentITGeneDataset_averageAllPatientScoreOfInterest);
                            currentITGene_allProteomics_allDatasetScoreOfInterest_basedOnMedian.Add(currentITGeneDataset_medianAllPatientScoreOfInterest);
                            break;
                        case KPMP_dataset_name_class.RNASeq_datasetClass:
                            currentITGene_allRNASeq_allDatasetScoreOfInterest_basedOnAverage.Add(currentITGeneDataset_averageAllPatientScoreOfInterest);
                            currentITGene_allRNASeq_allDatasetScoreOfInterest_basedOnMedian.Add(currentITGeneDataset_medianAllPatientScoreOfInterest);
                            break;
                        default:
                            throw new Exception();
                    }
                    if (!dataset_patientIDCount_dict.ContainsKey(standardized_line.Dataset))
                    {
                        dataset_patientIDCount_dict.Add(standardized_line.Dataset, currentITGeneDataset_allPatientScoreOfInterest_basedOnAverage.Count);
                    }
                    else if (dataset_patientIDCount_dict[standardized_line.Dataset]!= currentITGeneDataset_allPatientScoreOfInterest_basedOnAverage.Count)
                    {
                        throw new Exception();
                    }
                    standardized_dataset_collapseMethod_line = new KPMP_standardized_dataset_collapseOnMethod_line_class();
                    standardized_dataset_collapseMethod_line.Gene_symbol = (string)standardized_line.Gene_symbol.Clone();
                    standardized_dataset_collapseMethod_line.Score_of_interest_type = standardized_line.Score_of_interest_type;
                    standardized_dataset_collapseMethod_line.KPMP_data_integration_term = (string)standardized_line.KPMP_data_integration_term.Clone();
                    standardized_dataset_collapseMethod_line.Value_type = standardized_line.Value_type_1st;
                    standardized_dataset_collapseMethod_line.Method = (string)standardized_line.Dataset.Clone();

                    standardized_dataset_collapseMethod_line.ScoreOfInterest_average = currentITGeneDataset_averageAllPatientScoreOfInterest;
                    standardized_dataset_collapseMethod_line.ScoreOfInterest_populationSD = currentITGeneDataset_populationSDAllPatientScoreOfInterest;
                    standardized_dataset_collapseMethod_line.ScoreOfInterest_populationCV = standardized_dataset_collapseMethod_line.ScoreOfInterest_populationSD / standardized_dataset_collapseMethod_line.ScoreOfInterest_average;
                    standardized_dataset_collapseMethod_line.ScoreOfInterest_median = currentITGeneDataset_medianAllPatientScoreOfInterest;
                    standardized_dataset_collapseMethod_line.ScoreOfInterest_min = currentITGeneDataset_minAllPatientScoreOfInterest;
                    standardized_dataset_collapseMethod_line.ScoreOfInterest_max = currentITGeneDataset_maxAllPatientScoreOfInterest;
                    collapse_list.Add(standardized_dataset_collapseMethod_line);
                }

                if ((indexData == standardized_kpmp_data_length - 1)
                    || (!standardized_line.Score_of_interest_type.Equals(standardized_dataset.KPMP_data[indexData + 1].Score_of_interest_type))
                    || (!standardized_line.Value_type_1st.Equals(standardized_dataset.KPMP_data[indexData + 1].Value_type_1st))
                    || (!standardized_line.KPMP_data_integration_term.Equals(standardized_dataset.KPMP_data[indexData + 1].KPMP_data_integration_term))
                    || (!standardized_line.Gene_symbol.Equals(standardized_dataset.KPMP_data[indexData + 1].Gene_symbol)))
                {
                    Math_class.Get_mean_and_population_sd(currentITGene_allProteomics_allDatasetScoreOfInterest_basedOnAverage.ToArray(), out averageScoreOfInterest, out popSDScoreOfInterest);
                    medianScoreOfInterest = Math_class.Get_median(currentITGene_allProteomics_allDatasetScoreOfInterest_basedOnMedian.ToArray());
                    Math_class.Get_max_min_of_array(currentITGeneDatasetPatient_allCellSegmentScoreOfInterest.ToArray(), out maxScoreOfInterest, out minScoreOfInterest);
                    standardized_dataset_collapseMethod_line = new KPMP_standardized_dataset_collapseOnMethod_line_class();
                    standardized_dataset_collapseMethod_line.Gene_symbol = (string)standardized_line.Gene_symbol.Clone();
                    standardized_dataset_collapseMethod_line.Score_of_interest_type = standardized_line.Score_of_interest_type;
                    standardized_dataset_collapseMethod_line.KPMP_data_integration_term = (string)standardized_line.KPMP_data_integration_term.Clone();
                    standardized_dataset_collapseMethod_line.Value_type = standardized_line.Value_type_1st;

                    standardized_dataset_collapseMethod_line.Method = (string)KPMP_dataset_name_class.Proteomics_datasetClass.Clone();
                    standardized_dataset_collapseMethod_line.ScoreOfInterest_average = averageScoreOfInterest;
                    standardized_dataset_collapseMethod_line.ScoreOfInterest_populationSD = popSDScoreOfInterest;
                    standardized_dataset_collapseMethod_line.ScoreOfInterest_populationCV = standardized_dataset_collapseMethod_line.ScoreOfInterest_populationSD / standardized_dataset_collapseMethod_line.ScoreOfInterest_average;
                    standardized_dataset_collapseMethod_line.ScoreOfInterest_median = medianScoreOfInterest;
                    standardized_dataset_collapseMethod_line.ScoreOfInterest_min = minScoreOfInterest;
                    standardized_dataset_collapseMethod_line.ScoreOfInterest_max = maxScoreOfInterest;
                    collapse_list.Add(standardized_dataset_collapseMethod_line);

                    Math_class.Get_mean_and_population_sd(currentITGene_allRNASeq_allDatasetScoreOfInterest_basedOnAverage.ToArray(), out averageScoreOfInterest, out popSDScoreOfInterest);
                    medianScoreOfInterest = Math_class.Get_median(currentITGene_allRNASeq_allDatasetScoreOfInterest_basedOnMedian.ToArray());
                    Math_class.Get_max_min_of_array(currentITGene_allRNASeq_allDatasetScoreOfInterest_basedOnAverage.ToArray(), out maxScoreOfInterest, out minScoreOfInterest);
                    standardized_dataset_collapseMethod_line = new KPMP_standardized_dataset_collapseOnMethod_line_class();
                    standardized_dataset_collapseMethod_line.Gene_symbol = (string)standardized_line.Gene_symbol.Clone();
                    standardized_dataset_collapseMethod_line.Score_of_interest_type = standardized_line.Score_of_interest_type;
                    standardized_dataset_collapseMethod_line.KPMP_data_integration_term = (string)standardized_line.KPMP_data_integration_term.Clone();
                    standardized_dataset_collapseMethod_line.Value_type = standardized_line.Value_type_1st;

                    standardized_dataset_collapseMethod_line.Method = (string)KPMP_dataset_name_class.RNASeq_datasetClass.Clone();
                    standardized_dataset_collapseMethod_line.ScoreOfInterest_average = averageScoreOfInterest;
                    standardized_dataset_collapseMethod_line.ScoreOfInterest_populationSD = popSDScoreOfInterest;
                    standardized_dataset_collapseMethod_line.ScoreOfInterest_populationCV = standardized_dataset_collapseMethod_line.ScoreOfInterest_populationSD / standardized_dataset_collapseMethod_line.ScoreOfInterest_average;
                    standardized_dataset_collapseMethod_line.ScoreOfInterest_median = medianScoreOfInterest;
                    standardized_dataset_collapseMethod_line.ScoreOfInterest_min = minScoreOfInterest;
                    standardized_dataset_collapseMethod_line.ScoreOfInterest_max = maxScoreOfInterest;
                    collapse_list.Add(standardized_dataset_collapseMethod_line);
                }
            }
            this.Collapse = collapse_list.ToArray();
        }

        public void Write(string directory, string fileName)
        {
            KPMP_standardized_dataset_collapsOnMethod_readWriteOptions readWriteOptions = new KPMP_standardized_dataset_collapsOnMethod_readWriteOptions(directory, fileName);
            ReadWriteClass.WriteData(this.Collapse, readWriteOptions);
        }

    }
}
