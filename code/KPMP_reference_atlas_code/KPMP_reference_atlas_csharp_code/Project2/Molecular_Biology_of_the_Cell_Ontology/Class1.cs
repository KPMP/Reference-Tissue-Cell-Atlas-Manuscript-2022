using System;
using System.Collections.Generic;
using System.Linq;
using System.Drawing;
using Enumerations;
using ReadWrite;
using KPMP;

namespace MBCO_input_for_windows_form_application
{
    enum Significance_value_type_enum { E_m_p_t_y, Minus_log10_pvalue, Minus_log10_adjusted_pvalue, Signed_minus_log10_pvalue, Log_fold_change, Log2_fold_change }

    class Color_conversion_class
    {
        public static string Get_color_string(System.Drawing.Color color)
        {
            string color_string = color.ToString().Replace("Color ", "").Replace("[", "").Replace("]", "");
            return color_string;
        }

        public static System.Drawing.Color Set_color_from_string(string color_string)
        {
            System.Drawing.Color return_color = System.Drawing.Color.FromName(color_string);
            return return_color;
        }
    }

    class MBCO_windows_form_input_line_class
    {
        public string Center { get; set; }
        public string Dataset_name { get; set; }
        public string Integration_group { get; set; }
        public string NCBI_official_gene_symbol { get; set; }
        public double Minus_log10_pval_or_adj_pval { get; set; }
        public Significance_value_type_enum Significance_value_type { get; set; }
        public double Log2_fold_change { get; set; }
        public Significance_value_type_enum Fold_change_value_type { get; set; }
        public float Timepoint { get; set; }
        public Color Dataset_color_struct { get; set; }
        public string Dataset_color
        {
            get { return Color_conversion_class.Get_color_string(Dataset_color_struct); }
            set { Dataset_color_struct = Color_conversion_class.Set_color_from_string(value); }
        }

        public MBCO_windows_form_input_line_class Deep_copy()
        {
            MBCO_windows_form_input_line_class copy = (MBCO_windows_form_input_line_class)this.MemberwiseClone();
            copy.Dataset_name = (string)this.Dataset_name.Clone();
            copy.Integration_group = (string)this.Integration_group.Clone();
            return copy;
        }
    }

    class MBCO_windows_form_input_readWriteOptions_class : ReadWriteOptions_base
    {
        public MBCO_windows_form_input_readWriteOptions_class(string directory, string fileName)
        {
            this.File = directory + fileName;
            List<string> key_properties_list = new List<string>();
            key_properties_list.Add("Integration_group");
            key_properties_list.Add("Center");
            key_properties_list.Add("Dataset_name");
            key_properties_list.Add("Dataset_color");
            key_properties_list.Add("NCBI_official_gene_symbol");
            key_properties_list.Add("Minus_log10_pval_or_adj_pval");
            key_properties_list.Add("Significance_value_type");
            key_properties_list.Add("Log2_fold_change");
            key_properties_list.Add("Fold_change_value_type");
            this.Key_propertyNames = key_properties_list.ToArray();
            this.Key_columnNames = this.Key_propertyNames;
            this.HeadlineDelimiters = new char[] { Global_class.Tab };
            this.LineDelimiters = new char[] { Global_class.Tab };
            this.File_has_headline = true;
            this.Report = ReadWrite_report_enum.Report_main;
        }
    }

    class MBCO_windows_form_input_class
    {
        public MBCO_windows_form_input_line_class[] MBCO_input_lines { get; set; }

        public Dictionary<string, string[]> BgGenes_dict { get; set; }
        public MBCO_windows_form_input_class()
        {
            this.MBCO_input_lines = new MBCO_windows_form_input_line_class[0];
            this.BgGenes_dict = new Dictionary<string, string[]>();
        }

        private void Add_to_array(MBCO_windows_form_input_line_class[] add_mbco_input_lines)
        {
            int this_length = this.MBCO_input_lines.Length;
            int add_length = add_mbco_input_lines.Length;
            int new_length = this_length + add_length;
            MBCO_windows_form_input_line_class[] new_mbco_input_lines = new MBCO_windows_form_input_line_class[new_length];
            int indexNew = -1;
            for (int indexThis = 0; indexThis < this_length; indexThis++)
            {
                indexNew++;
                new_mbco_input_lines[indexNew] = this.MBCO_input_lines[indexNew];
            }
            for (int indexAdd = 0; indexAdd < add_length; indexAdd++)
            {
                indexNew++;
                new_mbco_input_lines[indexNew] = add_mbco_input_lines[indexAdd];
            }
            this.MBCO_input_lines = new_mbco_input_lines;
        }

        private void Add_to_bgGenes_dict(string bgGenes_name, string[] bgGenes)
        {
            //   if (BgGenes_dict.Keys.ToArray().Length == 1) { throw new Exception(); }
            if (!BgGenes_dict.ContainsKey(bgGenes_name))
            {
                BgGenes_dict.Add((string)bgGenes_name.Clone(), Array_class.Deep_copy_string_array(bgGenes));
            }
            else
            {
                foreach (string bgGene in bgGenes)
                {
                    if (!Array_class.Equal_arrays(BgGenes_dict[bgGenes_name], bgGenes)) { throw new Exception(); }
                }
            }
        }

        public void Clear()
        {
            this.MBCO_input_lines = new MBCO_windows_form_input_line_class[0];
            this.BgGenes_dict.Clear();
        }

        public void Set_all_dataset_names(string dataset_name)
        {
            foreach (MBCO_windows_form_input_line_class mbco_windows_form_line in this.MBCO_input_lines)
            {
                mbco_windows_form_line.Dataset_name = (string)dataset_name.Clone();
            }
        }

        public void Set_all_colors(Color new_color)
        {
            foreach (MBCO_windows_form_input_line_class mbco_windows_form_line in this.MBCO_input_lines)
            {
                mbco_windows_form_line.Dataset_color_struct = new_color;
            }
        }

        public void Set_color_for_indicated_constellation(string datasetName, float timepoint, int direction_of_change, Color new_color)
        {
            bool constellation_found = false;
            foreach (MBCO_windows_form_input_line_class mbco_windows_form_line in this.MBCO_input_lines)
            {
                if ((mbco_windows_form_line.Dataset_name.Equals(datasetName))
                    && (mbco_windows_form_line.Timepoint.Equals(timepoint)))
                {
                    if ((direction_of_change > 0) && (mbco_windows_form_line.Log2_fold_change > 0))
                    {
                        mbco_windows_form_line.Dataset_color_struct = new_color;
                        constellation_found = true;
                    }
                    else if ((direction_of_change < 0) && (mbco_windows_form_line.Log2_fold_change < 0))
                    {
                        mbco_windows_form_line.Dataset_color_struct = new_color;
                        constellation_found = true;
                    }
                    else if (mbco_windows_form_line.Log2_fold_change == 0) { throw new Exception(); }
                }
            }
            if (!constellation_found) { throw new Exception(); }
        }

        public void Set_integrationGroup_for_indicated_constellation(string datasetName, float timepoint, int direction_of_change, string integrationGroup)
        {
            bool constellation_found = false;
            foreach (MBCO_windows_form_input_line_class mbco_windows_form_line in this.MBCO_input_lines)
            {
                if ((mbco_windows_form_line.Dataset_name.Equals(datasetName))
                    && (mbco_windows_form_line.Timepoint.Equals(timepoint)))
                {
                    if ((direction_of_change > 0) && (mbco_windows_form_line.Log2_fold_change > 0))
                    {
                        mbco_windows_form_line.Integration_group = (string)integrationGroup.Clone();
                        constellation_found = true;
                    }
                    else if ((direction_of_change < 0) && (mbco_windows_form_line.Log2_fold_change < 0))
                    {
                        mbco_windows_form_line.Integration_group = (string)integrationGroup.Clone();
                        constellation_found = true;
                    }
                    else if (mbco_windows_form_line.Log2_fold_change == 0) { throw new Exception(); }
                }
            }
            if (!constellation_found) { throw new Exception(); }
        }


        public void Generate_from_KPMP_standardized_dataset_instance_and_add_to_array(KPMP_standardized_dataset_class kpmp_standad)
        {
            int kpmp_length = kpmp_standad.KPMP_data.Length;
            KPMP_standardized_dataset_line_class kpmp_data_line;
            MBCO_windows_form_input_line_class new_mbco_windows_line;
            List<MBCO_windows_form_input_line_class> new_mbco_windows_lines = new List<MBCO_windows_form_input_line_class>();
            bool pvalue_set;
            bool log2_fc_set;
            for (int indexKPMP = 0; indexKPMP < kpmp_length; indexKPMP++)
            {
                kpmp_data_line = kpmp_standad.KPMP_data[indexKPMP];
                new_mbco_windows_line = new MBCO_windows_form_input_line_class();
                pvalue_set = false;
                log2_fc_set = false;
                new_mbco_windows_line.Dataset_name = (string)kpmp_data_line.Dataset.Clone() + " - " + (string)kpmp_data_line.Cell_segment.Clone();
                new_mbco_windows_line.Center = (string)kpmp_data_line.Dataset.Clone();
                new_mbco_windows_line.Integration_group = (string)kpmp_data_line.KPMP_data_integration_term.Replace("-allPatients", "");
                new_mbco_windows_line.NCBI_official_gene_symbol = (string)kpmp_data_line.Gene_symbol.Clone();
                switch (kpmp_data_line.Value_type_1st)
                {
                    case Highthroughput_data.KPMP_value_type_enum.Minus_log10_pvalue:
                        if (pvalue_set) { throw new Exception(); }
                        new_mbco_windows_line.Significance_value_type = Significance_value_type_enum.Minus_log10_pvalue;
                        new_mbco_windows_line.Minus_log10_pval_or_adj_pval = kpmp_data_line.Value_1st;
                        pvalue_set = true;
                        break;
                    case Highthroughput_data.KPMP_value_type_enum.Minus_log10_pvalue_adjusted:
                        if (pvalue_set) { throw new Exception(); }
                        new_mbco_windows_line.Significance_value_type = Significance_value_type_enum.Minus_log10_pvalue;
                        new_mbco_windows_line.Minus_log10_pval_or_adj_pval = kpmp_data_line.Value_1st;
                        pvalue_set = true;
                        break;
                    case Highthroughput_data.KPMP_value_type_enum.Log2_ratioavg:
                        if (log2_fc_set) { throw new Exception(); }
                        new_mbco_windows_line.Log2_fold_change = kpmp_data_line.Value_1st;
                        new_mbco_windows_line.Fold_change_value_type = Significance_value_type_enum.Log2_fold_change;
                        log2_fc_set = true;
                        break;
                    case Highthroughput_data.KPMP_value_type_enum.Log_ratioavg:
                        if (log2_fc_set) { throw new Exception(); }
                        new_mbco_windows_line.Log2_fold_change = kpmp_data_line.Value_1st;
                        new_mbco_windows_line.Fold_change_value_type = Significance_value_type_enum.Log_fold_change;
                        log2_fc_set = true;
                        break;
                    default:
                        throw new Exception();
                }
                switch (kpmp_data_line.Value_type_2nd)
                {
                    case Highthroughput_data.KPMP_value_type_enum.Minus_log10_pvalue:
                        if (pvalue_set) { throw new Exception(); }
                        new_mbco_windows_line.Significance_value_type = Significance_value_type_enum.Minus_log10_pvalue;
                        new_mbco_windows_line.Minus_log10_pval_or_adj_pval = kpmp_data_line.Value_2nd;
                        pvalue_set = true;
                        break;
                    case Highthroughput_data.KPMP_value_type_enum.Minus_log10_pvalue_adjusted:
                        if (pvalue_set) { throw new Exception(); }
                        new_mbco_windows_line.Significance_value_type = Significance_value_type_enum.Minus_log10_pvalue;
                        new_mbco_windows_line.Minus_log10_pval_or_adj_pval = kpmp_data_line.Value_2nd;
                        pvalue_set = true;
                        break;
                    case Highthroughput_data.KPMP_value_type_enum.Log2_ratioavg:
                        if (log2_fc_set) { throw new Exception(); }
                        new_mbco_windows_line.Log2_fold_change = kpmp_data_line.Value_2nd;
                        new_mbco_windows_line.Fold_change_value_type = Significance_value_type_enum.Log2_fold_change;
                        log2_fc_set = true;
                        break;
                    case Highthroughput_data.KPMP_value_type_enum.Log_ratioavg:
                        if (log2_fc_set) { throw new Exception(); }
                        new_mbco_windows_line.Log2_fold_change = kpmp_data_line.Value_2nd;
                        new_mbco_windows_line.Fold_change_value_type = Significance_value_type_enum.Log_fold_change;
                        log2_fc_set = true;
                        break;
                    default:
                        throw new Exception();
                }
                if ((!pvalue_set) || (!log2_fc_set)) { throw new Exception(); }
                new_mbco_windows_lines.Add(new_mbco_windows_line);
            }
            string[] bgGeneList_names = kpmp_standad.Dataset_bgGenesProteins_dict.Keys.ToArray();
            foreach (string bgGeneList_name in bgGeneList_names)
            {
                Add_to_bgGenes_dict(bgGeneList_name, kpmp_standad.Dataset_bgGenesProteins_dict[bgGeneList_name]);
            }
            Add_to_array(new_mbco_windows_lines.ToArray());
        }

        public void Add_colors(Dictionary<string,Color> datasetname_color_dict)
        {
            foreach (MBCO_windows_form_input_line_class mbco_input_line in this.MBCO_input_lines)
            {
                mbco_input_line.Dataset_color = Color_conversion_class.Get_color_string(datasetname_color_dict[mbco_input_line.Center]);
            }
        }

        public void Add_other(MBCO_windows_form_input_class other)
        {
            this.Add_to_array(other.MBCO_input_lines);
            string[] bgGeneNames = other.BgGenes_dict.Keys.ToArray();
            foreach (string bgGeneName in bgGeneNames)
            {
                Add_to_bgGenes_dict(bgGeneName, other.BgGenes_dict[bgGeneName]);
            }
        }

        public void Write_one_file_for_each_center(string directory)
        {
            this.MBCO_input_lines = this.MBCO_input_lines.OrderBy(l => l.Center).ThenBy(l=>l.Dataset_name).ThenByDescending(l => l.Minus_log10_pval_or_adj_pval).ThenBy(l => Math.Abs(l.Log2_fold_change)).ToArray();
            MBCO_windows_form_input_line_class mbco_input_line;
            int mbco_input_length = this.MBCO_input_lines.Length;
            List<MBCO_windows_form_input_line_class> sameDataset_lines = new List<MBCO_windows_form_input_line_class>();
            for (int indexMBCO=0; indexMBCO<mbco_input_length;indexMBCO++)
            {
                mbco_input_line = this.MBCO_input_lines[indexMBCO];
                if ((indexMBCO==0)||(!mbco_input_line.Center.Equals(this.MBCO_input_lines[indexMBCO-1].Center)))
                {
                    sameDataset_lines.Clear();
                }
                sameDataset_lines.Add(mbco_input_line);
                if ((indexMBCO == mbco_input_length-1) || (!mbco_input_line.Center.Equals(this.MBCO_input_lines[indexMBCO + 1].Center)))
                {
                    MBCO_windows_form_input_readWriteOptions_class readWriteOptions = new MBCO_windows_form_input_readWriteOptions_class(directory, mbco_input_line.Center + ".txt");
                    ReadWriteClass.WriteData(sameDataset_lines.ToArray(), readWriteOptions);
                }
            }

            string[] bgGeneList_names = this.BgGenes_dict.Keys.ToArray();
            foreach (string bgGeneList_name in bgGeneList_names)
            {
                ReadWriteClass.WriteArray(this.BgGenes_dict[bgGeneList_name], directory + bgGeneList_name + "_bgGenes.txt");
            }
        }
    }
}
