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
using System.IO;
using ReadWrite;
using Statistic;
using Enumerations;
using KPMP;
using Enrichment_2018;
using Gene_databases;

namespace Highthroughput_data
{
    enum KPMP_value_type_enum { E_m_p_t_y, Single_value, Scaled_single_value, Average, Median, Ratioavg, Log_ratioavg, Log2_ratioavg, Log2_single_value, Log2_ratio_singlepatient, Ratio_singlepatient, Scaled_log2_single_value, Scaled_log2_ratio_singlepatient, Scaled_ratio_singlepatient, Minus_log10_pvalue, Pvalue, Minus_log10_pvalue_adjusted, Pvalue_adjusted, Fdr, Minus_log10_fdr, No_selection, Correlation }
    enum KPMP_DEG_value_type_enum { E_m_p_t_y, Minus_log10_pvalue, Log_fold_change, Abs_expression }
    enum KPMP_metabolite_id_type_enum {  E_m_p_t_y, Hmdb, Kegg_id }
    enum KPMP_analysis_set_enum {  E_m_p_t_y, Pilot_data, Disease_data_2021 }
    enum SingleCell_valueType_enum { E_m_p_t_y, Signed_minus_log10_pvalue, Log_fold_change, Signed_pct1, Signed_pct2, Signed_fractionalRank_pvalue, Signed_fractionalRank_absAvgLogFc, }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class Global_kpmp_class
    {
        public const KPMP_analysis_set_enum Analysis_set = KPMP_analysis_set_enum.Pilot_data;

        public static string Get_bgGenesInUpperCase_completeFileName(string subdirectory, string dataset_name)
        {
            string complete_fileName = Global_directory_class.Experimental_data_directory + subdirectory + "BgGenesInUpperCase_" + dataset_name + ".txt";
            return complete_fileName;
        }

    }


    class KPMP_singleCell_annotation_line_class
    {
        public string Barcode { get; set; }
        public string Lmd_segment { get; set; }
        public int Cluster_no { get; set; }
        public string Cell_subtype { get; set; }

        public KPMP_singleCell_annotation_line_class Deep_copy()
        {
            KPMP_singleCell_annotation_line_class copy = (KPMP_singleCell_annotation_line_class)this.MemberwiseClone();
            copy.Barcode = (string)this.Barcode.Clone();
            copy.Lmd_segment = (string)this.Lmd_segment.Clone();
            copy.Cell_subtype = (string)this.Cell_subtype.Clone();
            return copy;
        }
    }

    class KPMP_singleCell_annotation_readWriteOptions_class : ReadWriteOptions_base
    {
        public KPMP_singleCell_annotation_readWriteOptions_class(string fileName)
        {
            this.File = Global_directory_class.Seurat_singleCell_data_directory + fileName;
            this.Key_propertyNames = new string[] { "Barcode", "Lmd_segment", "Cell_subtype" };
            this.Key_columnNames = new string[] { "Barcode", "Lmd_segment", "Cluster" };
            this.HeadlineDelimiters = new char[] { Global_class.Tab };
            this.LineDelimiters = new char[] { Global_class.Tab };
            this.File_has_headline = true;
            this.Report = ReadWrite_report_enum.Report_main;
        }
    }

    class KPMP_singleCell_annotation_class
    {
        public KPMP_singleCell_annotation_line_class[] Annotations { get; set; }

        public KPMP_singleCell_annotation_class()
        { }

        private void Adjust_annotations_premiere()
        {
            int annotations_length = Annotations.Length;
            KPMP_singleCell_annotation_line_class singleCell_annotation_line;
            for (int indexA=0; indexA<annotations_length;indexA++)
            {
                singleCell_annotation_line = this.Annotations[indexA];
                switch (singleCell_annotation_line.Cell_subtype)
                {
                    case "B cell":
                        singleCell_annotation_line.Cell_subtype = "B Cell";
                        break;
                    case "PEC":
                        singleCell_annotation_line.Cell_subtype = "PEC/LOH";
                        break;
                    case "T-MEM":
                        singleCell_annotation_line.Cell_subtype = "T-CYT-MEM";
                        break;
                    case "T-ACT":
                        singleCell_annotation_line.Cell_subtype = "T Cell";
                        break;
                    default:
                        break;
                }
            }
        }

        private void Adjust_annotations_ucsd()
        {
            int annotations_length = Annotations.Length;
            KPMP_singleCell_annotation_line_class singleCell_annotation_line;
            for (int indexA = 0; indexA < annotations_length; indexA++)
            {
                singleCell_annotation_line = this.Annotations[indexA];
                switch (singleCell_annotation_line.Cell_subtype)
                {
                    case "DL":
                        singleCell_annotation_line.Cell_subtype = "DTL";
                        break;
                    default:
                        break;
                }
            }
        }

        public void Generate_by_reading(string fileName)
        {
            Read_file(fileName);
            if (fileName.IndexOf("premiere") != -1)
            {
                Adjust_annotations_premiere();
            }
            if (fileName.IndexOf("ucsd") != -1)
            {
                Adjust_annotations_ucsd();
            }
        }

        private void Read_file(string fileName)
        {
            KPMP_singleCell_annotation_readWriteOptions_class readWriteOptions = new KPMP_singleCell_annotation_readWriteOptions_class(fileName);
            this.Annotations = ReadWriteClass.ReadRawData_and_FillArray<KPMP_singleCell_annotation_line_class>(readWriteOptions);
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class KPMP_singleCell_subtype_count_line_class
    {
        public string Center { get; set; }
        public string Cell_subtype { get; set; }
        public string Cell_type { get; set; }
        public string KPMP_integration_term { get; set; }
        public int Cell_counts { get; set; }

        public KPMP_singleCell_subtype_count_line_class Deep_copy()
        {
            KPMP_singleCell_subtype_count_line_class copy = (KPMP_singleCell_subtype_count_line_class)this.MemberwiseClone();
            copy.Center = (string)this.Center.Clone();
            copy.Cell_subtype = (string)this.Cell_subtype.Clone();
            copy.Cell_type = (string)this.Cell_type.Clone();
            copy.KPMP_integration_term = (string)this.KPMP_integration_term.Clone();
            return copy;
        }
    }

    class KPMP_singleCell_subtype_count_readWriteOptions_class : ReadWriteOptions_base
    {
        public KPMP_singleCell_subtype_count_readWriteOptions_class(string directory, string fileName)
        {
            this.File = directory + fileName;
            this.Key_propertyNames = new string[] { "Center", "Cell_subtype", "Cell_type", "KPMP_integration_term", "Cell_counts" };
            this.Key_columnNames = this.Key_propertyNames;
            this.HeadlineDelimiters = new char[] { Global_class.Tab };
            this.LineDelimiters = new char[] { Global_class.Tab };
            this.File_has_headline = true;
            this.Report = ReadWrite_report_enum.Report_main;
        }
    }

    class KPMP_singleCell_subtype_count_class
    {
        public KPMP_singleCell_subtype_count_line_class[] Subtype_counts { get; set; }

        public KPMP_singleCell_subtype_count_class()
        {
            this.Subtype_counts = new KPMP_singleCell_subtype_count_line_class[0];
        }

        private void Add_to_array(KPMP_singleCell_subtype_count_line_class[] add_subtype_counts)
        {
            int this_length = this.Subtype_counts.Length;
            int add_length = add_subtype_counts.Length;
            int new_length = this_length + add_length;
            int indexNew = -1;
            KPMP_singleCell_subtype_count_line_class[] new_subtype_counts = new KPMP_singleCell_subtype_count_line_class[new_length];
            for (int indexThis = 0; indexThis < this_length; indexThis++)
            {
                indexNew++;
                new_subtype_counts[indexNew] = this.Subtype_counts[indexThis];
            }
            for (int indexAdd = 0; indexAdd < add_length; indexAdd++)
            {
                indexNew++;
                new_subtype_counts[indexNew] = add_subtype_counts[indexAdd];
            }
            this.Subtype_counts = new_subtype_counts;
        }

        #region Generate
        private KPMP_singleCell_subtype_count_line_class[] Generate_from_annotation_and_return_initial_array(KPMP_singleCell_annotation_class annotation, string center)
        {
            List<KPMP_singleCell_subtype_count_line_class> add_subtype_counts = new List<KPMP_singleCell_subtype_count_line_class>();
            KPMP_singleCell_subtype_count_line_class new_singleCell_subtype_count_line;
            annotation.Annotations = annotation.Annotations.OrderBy(l => l.Cell_subtype).ThenBy(l=>l.Barcode).ToArray();
            int annotations_length = annotation.Annotations.Length;
            KPMP_singleCell_annotation_line_class annotation_line;
            int barcodes_count = 0;
            for (int indexAnn=0; indexAnn<annotations_length;indexAnn++)
            {
                annotation_line = annotation.Annotations[indexAnn];
                if ((indexAnn==0)||(!annotation_line.Cell_subtype.Equals(annotation.Annotations[indexAnn-1].Cell_subtype)))
                {
                    barcodes_count = 0;
                }
                if (  (indexAnn != 0) 
                    &&(annotation_line.Barcode.Equals(annotation.Annotations[indexAnn - 1].Barcode)))
                { throw new Exception(); }
                barcodes_count++;
                if ((indexAnn == annotations_length-1) || (!annotation_line.Cell_subtype.Equals(annotation.Annotations[indexAnn + 1].Cell_subtype)))
                {
                    new_singleCell_subtype_count_line = new KPMP_singleCell_subtype_count_line_class();
                    new_singleCell_subtype_count_line.Cell_subtype = (string)annotation_line.Cell_subtype.Clone();
                    new_singleCell_subtype_count_line.Cell_type = "";
                    new_singleCell_subtype_count_line.Center = (string)center.Clone();
                    new_singleCell_subtype_count_line.KPMP_integration_term = "";
                    new_singleCell_subtype_count_line.Cell_counts = barcodes_count;
                    add_subtype_counts.Add(new_singleCell_subtype_count_line);
                }
            }
            return add_subtype_counts.ToArray();
        }

        private KPMP_singleCell_subtype_count_line_class[] Add_kpmp_data_integration_term(KPMP_singleCell_subtype_count_line_class[] single_cell_subtype_counts, string add_toIntegrationTerm)
        {
            int single_cell_subtype_counts_length = single_cell_subtype_counts.Length;
            KPMP_singleCell_subtype_count_line_class subtype_line;
            for (int indexST =0; indexST<single_cell_subtype_counts_length;indexST++)
            {
                subtype_line = single_cell_subtype_counts[indexST];
                subtype_line.KPMP_integration_term = KPMP_data_integration_class.Get_kpmp_integration_term(subtype_line.Cell_subtype) + add_toIntegrationTerm;
            }
            return single_cell_subtype_counts;
        }

        public void Generate_from_annotations_add_to_array(KPMP_singleCell_annotation_class annotation, string center, string addToKPMPIntegrationTerm)
        {
            KPMP_singleCell_subtype_count_line_class[] subtype_count_lines = Generate_from_annotation_and_return_initial_array(annotation, center);
            subtype_count_lines = Add_kpmp_data_integration_term(subtype_count_lines, addToKPMPIntegrationTerm);
            Add_to_array(subtype_count_lines);
        }
        #endregion

        public void Reset_kpmpIntegration_term_as_indicated_in_dictionary(Dictionary<string, Dictionary<string, string>> center_cellSubType_setIntegrationTerm_dict)
        {
            int subtype_counts_length = this.Subtype_counts.Length;
            KPMP_singleCell_subtype_count_line_class subtype_count_line;
            for (int indexSC=0; indexSC<subtype_counts_length; indexSC++)
            {
                subtype_count_line = this.Subtype_counts[indexSC];
                if (center_cellSubType_setIntegrationTerm_dict.ContainsKey(subtype_count_line.Center))
                {
                    if (center_cellSubType_setIntegrationTerm_dict[subtype_count_line.Center].ContainsKey(subtype_count_line.Cell_subtype))
                    {
                        subtype_count_line.KPMP_integration_term = (string)center_cellSubType_setIntegrationTerm_dict[subtype_count_line.Center][subtype_count_line.Cell_subtype].Clone();
                    }
                }
            }
        }

        public Dictionary<string, Dictionary<string, int>> Get_center_cell_subtype_counts_dict()
        {
            this.Subtype_counts = this.Subtype_counts.OrderBy(l => l.Center).ThenBy(l => l.Cell_subtype).ToArray();
            int subtype_counts_length = this.Subtype_counts.Length;
            KPMP_singleCell_subtype_count_line_class subtype_line;
            Dictionary<string, Dictionary<string, int>> center_cellSubtype_counts_dict = new Dictionary<string, Dictionary<string, int>>();
            for (int indexSC = 0; indexSC < subtype_counts_length; indexSC++)
            {
                subtype_line = this.Subtype_counts[indexSC];
                if (!center_cellSubtype_counts_dict.ContainsKey(subtype_line.Center)) { center_cellSubtype_counts_dict.Add(subtype_line.Center, new Dictionary<string, int>()); }
                center_cellSubtype_counts_dict[subtype_line.Center].Add(subtype_line.Cell_subtype, subtype_line.Cell_counts);
            }
            return center_cellSubtype_counts_dict;
        }

        public Dictionary<string, Dictionary<string, int>> Get_center_kpmpIntegrationTerm_counts_dict()
        {
            this.Subtype_counts = this.Subtype_counts.OrderBy(l => l.Center).ThenBy(l => l.KPMP_integration_term).ToArray();
            int subtype_counts_length = this.Subtype_counts.Length;
            KPMP_singleCell_subtype_count_line_class subtype_line;
            Dictionary<string, Dictionary<string, int>> center_kpmpIntegrationTerm_counts_dict = new Dictionary<string, Dictionary<string, int>>();
            int current_cell_counts = 0;
            for (int indexSC = 0; indexSC < subtype_counts_length; indexSC++)
            {
                subtype_line = this.Subtype_counts[indexSC];
                if (subtype_line.Cell_counts==0) { throw new Exception(); }
                if (  (indexSC==0)
                    || (!subtype_line.Center.Equals(this.Subtype_counts[indexSC-1].Center))
                    || (!subtype_line.KPMP_integration_term.Equals(this.Subtype_counts[indexSC-1].KPMP_integration_term)))
                {
                    current_cell_counts = 0;
                }
                current_cell_counts += subtype_line.Cell_counts;
                if ((indexSC == subtype_counts_length - 1)
                    || (!subtype_line.Center.Equals(this.Subtype_counts[indexSC + 1].Center))
                    || (!subtype_line.KPMP_integration_term.Equals(this.Subtype_counts[indexSC + 1].KPMP_integration_term)))
                {
                    if (!center_kpmpIntegrationTerm_counts_dict.ContainsKey(subtype_line.Center)) { center_kpmpIntegrationTerm_counts_dict.Add(subtype_line.Center, new Dictionary<string, int>()); }
                    center_kpmpIntegrationTerm_counts_dict[subtype_line.Center].Add(subtype_line.KPMP_integration_term, current_cell_counts);
                }
            }
            return center_kpmpIntegrationTerm_counts_dict;
        }

        public void Write(string directory, string fileName)
        {
            KPMP_singleCell_subtype_count_readWriteOptions_class readWriteOptions = new KPMP_singleCell_subtype_count_readWriteOptions_class(directory, fileName);
            ReadWriteClass.WriteData(this.Subtype_counts, readWriteOptions);
        }

        public KPMP_singleCell_subtype_count_class Deep_copy()
        {
            KPMP_singleCell_subtype_count_class copy = (KPMP_singleCell_subtype_count_class)this.MemberwiseClone();
            int subtype_counts_length = this.Subtype_counts.Length;
            copy.Subtype_counts = new KPMP_singleCell_subtype_count_line_class[subtype_counts_length];
            for (int indexST =0; indexST<subtype_counts_length;indexST++)
            {
                copy.Subtype_counts[indexST] = this.Subtype_counts[indexST].Deep_copy();
            }
            return copy;
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class KPMP_singleCell_averageMinusLog10Pvalue_line_class
    {
        public string Center { get; set; }
        public string Cell_type { get; set; }
        public string KPMP_data_integration_term { get; set; }
        public Ontology_type_enum Ontology { get; set; }
        public string Scp { get; set; }
        public double Averaged_minusLog10Pvalue { get; set; }
        public double Percentage_averaged_minusLog10Pvalue { get; set; }
        public int Averaged_entities_count { get; set; }
        public string[] Averaged_entities { get; set; }

        public string ReadWrite_averaged_entities
        {
            get { return ReadWriteClass.Get_writeLine_from_array(this.Averaged_entities, KPMP_singleCell_averageMinusLog10Pvalue_readWriteOptions_class.Array_delimiter); }
            set { this.Averaged_entities = ReadWriteClass.Get_array_from_readLine<string>(value, KPMP_singleCell_averageMinusLog10Pvalue_readWriteOptions_class.Array_delimiter); }
        }

        public KPMP_singleCell_averageMinusLog10Pvalue_line_class()
        {
            this.Averaged_entities = new string[0];
        }

        public KPMP_singleCell_averageMinusLog10Pvalue_line_class Deep_copy()
        {
            KPMP_singleCell_averageMinusLog10Pvalue_line_class copy = (KPMP_singleCell_averageMinusLog10Pvalue_line_class)this.MemberwiseClone();
            copy.Cell_type = (string)this.Cell_type.Clone();
            copy.Center = (string)this.Center.Clone();
            copy.KPMP_data_integration_term = (string)this.KPMP_data_integration_term.Clone();
            copy.Scp = (string)this.Scp.Clone();
            copy.Averaged_entities = Array_class.Deep_copy_string_array(this.Averaged_entities);
            return copy;
        }
    }

    class KPMP_singleCell_averageMinusLog10Pvalue_readWriteOptions_class: ReadWriteOptions_base
    {
        public static char Array_delimiter { get { return ';'; } }
        public KPMP_singleCell_averageMinusLog10Pvalue_readWriteOptions_class(string directory, string fileName)
        {
            this.File = directory + fileName;
            this.Key_propertyNames = new string[] { "KPMP_data_integration_term", "Center","Cell_type", "Ontology","Scp", "Averaged_minusLog10Pvalue", "Percentage_averaged_minusLog10Pvalue" ,"Averaged_entities_count","ReadWrite_averaged_entities"};
            this.Key_columnNames = this.Key_propertyNames;
            this.HeadlineDelimiters = new char[] { Global_class.Tab };
            this.LineDelimiters = new char[] { Global_class.Tab };
            this.File_has_headline = true;
            this.Report = ReadWrite_report_enum.Report_main;
        }
    }

    class KPMP_singleCell_averageMinusLog10Pvalue_class
    {
        public KPMP_singleCell_averageMinusLog10Pvalue_line_class[] Averaged_minuslog10p { get; set; }

        public KPMP_singleCell_averageMinusLog10Pvalue_class()
        {
            this.Averaged_minuslog10p = new KPMP_singleCell_averageMinusLog10Pvalue_line_class[0];
        }

        public void Add_to_array(KPMP_singleCell_averageMinusLog10Pvalue_line_class[] add_averaged_minuslog10p)
        {
            int this_length = this.Averaged_minuslog10p.Length;
            int add_length = add_averaged_minuslog10p.Length;
            int new_length = this_length + add_length;
            KPMP_singleCell_averageMinusLog10Pvalue_line_class[] new_averaged_minusLog10p = new KPMP_singleCell_averageMinusLog10Pvalue_line_class[new_length];
            int indexNew = -1;
            for (int indexThis = 0; indexThis < this_length; indexThis++)
            {
                indexNew++;
                new_averaged_minusLog10p[indexNew] = this.Averaged_minuslog10p[indexThis];
            }
            for (int indexAdd = 0; indexAdd < add_length; indexAdd++)
            {
                indexNew++;
                new_averaged_minusLog10p[indexNew] = add_averaged_minuslog10p[indexAdd];
            }
            this.Averaged_minuslog10p = new_averaged_minusLog10p;
        }

        private KPMP_singleCell_averageMinusLog10Pvalue_line_class[] Generate_from_dictionaries(Dictionary<string, Dictionary<string, double>> center_kpmpIntegrationTerm_minusLog10PvalueSum_dict, Dictionary<string, Dictionary<string, int>> center_kpmpIntegrationTerm_counts_dict, Ontology_type_enum ontology, string scp)
        {
            Dictionary<string, double> kpmpIntegrationTerm_minusLog10PvalueSum_dict;
            string[] centers = center_kpmpIntegrationTerm_minusLog10PvalueSum_dict.Keys.ToArray();
            string center;
            int centers_length = centers.Length;
            string[] kpmpIntegrationTerms;
            string kpmpIntegrationTerm;
            int kpmpIntegrationTerms_length;
            double current_minusLog10Pvalue_sum;
            int current_total_cell_counts =-1;
            List<KPMP_singleCell_averageMinusLog10Pvalue_line_class> singleCell_averageMinusLog10Pvalues = new List<KPMP_singleCell_averageMinusLog10Pvalue_line_class>();
            KPMP_singleCell_averageMinusLog10Pvalue_line_class new_averageMinusLog10Pvalue_line;
            for (int indexCenter = 0; indexCenter < centers_length; indexCenter++)
            {
                center = centers[indexCenter];
                kpmpIntegrationTerm_minusLog10PvalueSum_dict = center_kpmpIntegrationTerm_minusLog10PvalueSum_dict[center];
                kpmpIntegrationTerms = kpmpIntegrationTerm_minusLog10PvalueSum_dict.Keys.ToArray();
                kpmpIntegrationTerms_length = kpmpIntegrationTerms.Length;
                for (int indexInt = 0; indexInt < kpmpIntegrationTerms_length; indexInt++)
                {
                    kpmpIntegrationTerm = kpmpIntegrationTerms[indexInt];
                    current_minusLog10Pvalue_sum = kpmpIntegrationTerm_minusLog10PvalueSum_dict[kpmpIntegrationTerm];
                    current_total_cell_counts = -1;
                    if (center_kpmpIntegrationTerm_counts_dict.ContainsKey(center))
                    {
                        current_total_cell_counts = center_kpmpIntegrationTerm_counts_dict[center][kpmpIntegrationTerm];
                    }
                    else
                    {
                        current_total_cell_counts = 1;
                    }

                    new_averageMinusLog10Pvalue_line = new KPMP_singleCell_averageMinusLog10Pvalue_line_class();
                    new_averageMinusLog10Pvalue_line.Center = (string)center.Clone();
                    new_averageMinusLog10Pvalue_line.KPMP_data_integration_term = (string)kpmpIntegrationTerm.Clone();
                    new_averageMinusLog10Pvalue_line.Scp = (string)scp.Clone();
                    new_averageMinusLog10Pvalue_line.Ontology = ontology;
                    new_averageMinusLog10Pvalue_line.Averaged_minusLog10Pvalue = current_minusLog10Pvalue_sum / (double)current_total_cell_counts;
                    new_averageMinusLog10Pvalue_line.Averaged_entities_count = current_total_cell_counts;
                    if (Double.IsNaN(new_averageMinusLog10Pvalue_line.Averaged_minusLog10Pvalue)) { throw new Exception(); }
                    singleCell_averageMinusLog10Pvalues.Add(new_averageMinusLog10Pvalue_line);
                }
            }
            return singleCell_averageMinusLog10Pvalues.ToArray();
        }

        private KPMP_singleCell_averageMinusLog10Pvalue_line_class[] Generate_from_enrichment_results_and_cell_subtype_counts(Enrichment2018_results_class enrichment, KPMP_singleCell_subtype_count_class singleCell_counts)
        {
            Dictionary<string, Dictionary<string, int>> center_cellSubtype_counts_dict = singleCell_counts.Get_center_cell_subtype_counts_dict();
            Dictionary<string, Dictionary<string, int>> center_kpmpIntegrationTerm_counts_dict = singleCell_counts.Get_center_kpmpIntegrationTerm_counts_dict();

            Dictionary<string, Dictionary<string, double>> center_kpmpIntegrationTerm_minusLog10PvalueSum_dict = new Dictionary<string, Dictionary<string, double>>();
            Enrichment2018_results_line_class enrichment_results_line;
            enrichment.Enrichment_results = enrichment.Enrichment_results.OrderBy(l => l.Scp).ToArray();
            int enrichment_length = enrichment.Enrichment_results.Length;
            string[] splitStrings;
            string current_center;
            string current_kpmpIntegrationTerm;
            string current_cell_subtype="error";
            int current_cell_subtype_counts;
            bool Continue;
            int indexSplitString;
            List<string> missing_cellTypes = new List<string>();
            List<KPMP_singleCell_averageMinusLog10Pvalue_line_class> singleCell_averageMinusLog10Pvalues = new List<KPMP_singleCell_averageMinusLog10Pvalue_line_class>();
            KPMP_singleCell_averageMinusLog10Pvalue_line_class[] new_averageMinusLog10Pvalue_lines;
            for (int indexE=0; indexE<enrichment_length;indexE++)
            {
                enrichment_results_line = enrichment.Enrichment_results[indexE];
                if ((indexE == 0) || (!enrichment_results_line.Scp.Equals(enrichment.Enrichment_results[indexE-1].Scp)))
                {
                    center_kpmpIntegrationTerm_minusLog10PvalueSum_dict.Clear();
                }

                #region Identify center, kpmpIntegrationTerm and cell subtype
                splitStrings = enrichment_results_line.Sample_name.Split('-');
                indexSplitString = 0;
                current_center = splitStrings[indexSplitString];
                Continue = true;
                current_cell_subtype = "";
                while (Continue)
                {
                    indexSplitString++;
                    if (current_cell_subtype.Length > 0) { current_cell_subtype += "-"; }
                    current_cell_subtype += splitStrings[indexSplitString];
                    if (  (splitStrings[indexSplitString + 1].Equals("allPatients"))
                        || (splitStrings[indexSplitString + 1].Equals("aP")))
                    {
                        Continue = false;
                    }
                }
                indexSplitString += 2;
                Continue = true;
                current_kpmpIntegrationTerm = "";
                while (Continue)
                {
                    indexSplitString++;
                    if (current_kpmpIntegrationTerm.Length > 0) { current_kpmpIntegrationTerm += "-"; }
                    current_kpmpIntegrationTerm += splitStrings[indexSplitString];
                    if (indexSplitString == splitStrings.Length-1)
                    {
                        Continue = false;
                    }
                }
                #endregion
                current_cell_subtype_counts = -9999;
                if (center_cellSubtype_counts_dict.ContainsKey(current_center))
                {
                    if (center_cellSubtype_counts_dict[current_center].ContainsKey(current_cell_subtype))
                    {
                        current_cell_subtype_counts = center_cellSubtype_counts_dict[current_center][current_cell_subtype];
                    }
                    else
                    {
                        current_cell_subtype_counts = -9999;
                        missing_cellTypes.Add(current_center + ": " + current_cell_subtype);
                    }
                }
                else { current_cell_subtype_counts = 1; }
                if (!center_kpmpIntegrationTerm_minusLog10PvalueSum_dict.ContainsKey(current_center))
                {
                    center_kpmpIntegrationTerm_minusLog10PvalueSum_dict.Add(current_center, new Dictionary<string, double>());
                }
                if (!center_kpmpIntegrationTerm_minusLog10PvalueSum_dict[current_center].ContainsKey(current_kpmpIntegrationTerm))
                {
                    center_kpmpIntegrationTerm_minusLog10PvalueSum_dict[current_center].Add(current_kpmpIntegrationTerm, 0);
                }
                center_kpmpIntegrationTerm_minusLog10PvalueSum_dict[current_center][current_kpmpIntegrationTerm] += enrichment_results_line.Minus_log10_pvalue * current_cell_subtype_counts;
                if ((indexE == enrichment_length-1) || (!enrichment_results_line.Scp.Equals(enrichment.Enrichment_results[indexE + 1].Scp)))
                {
                    new_averageMinusLog10Pvalue_lines = Generate_from_dictionaries(center_kpmpIntegrationTerm_minusLog10PvalueSum_dict, center_kpmpIntegrationTerm_counts_dict,enrichment_results_line.Ontology, enrichment_results_line.Scp);
                    singleCell_averageMinusLog10Pvalues.AddRange(new_averageMinusLog10Pvalue_lines);
                }
            }
            missing_cellTypes = missing_cellTypes.Distinct().ToList();
            if (missing_cellTypes.Count>0) { throw new Exception(); }
            return singleCell_averageMinusLog10Pvalues.ToArray();
        }

        private KPMP_singleCell_averageMinusLog10Pvalue_line_class[] Calculate_percentage_of_minusLog10Pvalues_over_integration_term_and_center(KPMP_singleCell_averageMinusLog10Pvalue_line_class[] new_averageMinusLog10Pvalue_lines)
        {
            int averaged_minusLog10p_length = new_averageMinusLog10Pvalue_lines.Length;
            KPMP_singleCell_averageMinusLog10Pvalue_line_class average_line;
            KPMP_singleCell_averageMinusLog10Pvalue_line_class inner_average_line;
            new_averageMinusLog10Pvalue_lines = new_averageMinusLog10Pvalue_lines.OrderBy(l => l.Center).ThenBy(l => l.KPMP_data_integration_term).ToArray();
            int firstIndexSameGroup = 0;
            double sum_of_minusLog10Pvalues = 0;
            for (int indexSC=0; indexSC< averaged_minusLog10p_length; indexSC++)
            {
                average_line = new_averageMinusLog10Pvalue_lines[indexSC];
                if (  (indexSC==0)
                    || (!average_line.KPMP_data_integration_term.Equals(new_averageMinusLog10Pvalue_lines[indexSC - 1].KPMP_data_integration_term))
                    || (!average_line.Center.Equals(new_averageMinusLog10Pvalue_lines[indexSC-1].Center)))
                {
                    firstIndexSameGroup = indexSC;
                    sum_of_minusLog10Pvalues = 0;
                }
                sum_of_minusLog10Pvalues += average_line.Averaged_minusLog10Pvalue;
                if ((indexSC == averaged_minusLog10p_length-1)
                    || (!average_line.KPMP_data_integration_term.Equals(new_averageMinusLog10Pvalue_lines[indexSC+1].KPMP_data_integration_term))
                    || (!average_line.Center.Equals(new_averageMinusLog10Pvalue_lines[indexSC+1].Center)))
                {
                    for (int indexInner=firstIndexSameGroup;indexInner<=indexSC;indexInner++)
                    {
                        inner_average_line = new_averageMinusLog10Pvalue_lines[indexInner];
                        if (sum_of_minusLog10Pvalues != 0)
                        {
                            inner_average_line.Percentage_averaged_minusLog10Pvalue = (double)100 * inner_average_line.Averaged_minusLog10Pvalue / sum_of_minusLog10Pvalues;
                        }
                        else
                        {
                            inner_average_line.Percentage_averaged_minusLog10Pvalue = 0;
                        }
                    }
                }
            }
            return new_averageMinusLog10Pvalue_lines;
        }

        public void Generate_from_enrichment_results_and_cell_subtype_counts_and_add_to_array(Enrichment2018_results_class enrichment, KPMP_singleCell_subtype_count_class singleCell_counts)
        {
            KPMP_singleCell_averageMinusLog10Pvalue_line_class[] new_averageMinusLog10Pvalue_lines = Generate_from_enrichment_results_and_cell_subtype_counts(enrichment, singleCell_counts);
            new_averageMinusLog10Pvalue_lines = Calculate_percentage_of_minusLog10Pvalues_over_integration_term_and_center(new_averageMinusLog10Pvalue_lines);
            Add_to_array(new_averageMinusLog10Pvalue_lines);
        }

        public void Keep_only_indicated_centers(params string[] keep_centers)
        {
            keep_centers = keep_centers.Distinct().ToArray();
            Dictionary<string, bool> keep_center_dict = new Dictionary<string, bool>();
            foreach (string keep_center in keep_centers)
            {
                keep_center_dict.Add(keep_center, true);
            }
            List<KPMP_singleCell_averageMinusLog10Pvalue_line_class> keep = new List<KPMP_singleCell_averageMinusLog10Pvalue_line_class>();
            foreach (KPMP_singleCell_averageMinusLog10Pvalue_line_class average_line in this.Averaged_minuslog10p)
            {
                if (keep_center_dict.ContainsKey(average_line.Center))
                {
                    keep.Add(average_line);
                }
            }
            this.Averaged_minuslog10p = keep.ToArray();
        }

        private Dictionary<string,int> Generate_integrationTerm_centerCounts_dict()
        {
            Dictionary<string, int> integrationTerm_centerCounts_dict = new Dictionary<string, int>();
            this.Averaged_minuslog10p = this.Averaged_minuslog10p.OrderBy(l => l.KPMP_data_integration_term).ThenBy(l => l.Center).ToArray();
            int average_length = this.Averaged_minuslog10p.Length;
            KPMP_singleCell_averageMinusLog10Pvalue_line_class average_line;
            int center_counts = 0;
            for (int indexA=0; indexA<average_length;indexA++)
            {
                average_line = this.Averaged_minuslog10p[indexA];
                if ((indexA == 0)
                    || (!average_line.KPMP_data_integration_term.Equals(this.Averaged_minuslog10p[indexA - 1].KPMP_data_integration_term)))
                {
                    center_counts = 0;
                }
                if ((indexA == 0)
                    || (!average_line.Center.Equals(this.Averaged_minuslog10p[indexA - 1].Center))
                    || (!average_line.KPMP_data_integration_term.Equals(this.Averaged_minuslog10p[indexA - 1].KPMP_data_integration_term)))
                {
                    center_counts++;
                }
                if ((indexA == average_length-1)
                    || (!average_line.KPMP_data_integration_term.Equals(this.Averaged_minuslog10p[indexA + 1].KPMP_data_integration_term)))
                {
                    integrationTerm_centerCounts_dict.Add(average_line.KPMP_data_integration_term, center_counts);
                }
            }
            return integrationTerm_centerCounts_dict;
        }

        public void Generate_new_enrichment_results_lines_by_averaging_average_values_over_all_assays_of_same_integration_term_and_add_to_array(string new_center_name)
        {
            Dictionary<string, int> integrationTerm_centerCounts_dict = Generate_integrationTerm_centerCounts_dict();
            int average_length = this.Averaged_minuslog10p.Length;
            KPMP_singleCell_averageMinusLog10Pvalue_line_class average_line;
            KPMP_singleCell_averageMinusLog10Pvalue_line_class new_average_line;
            List<KPMP_singleCell_averageMinusLog10Pvalue_line_class> new_average_lines_list = new List<KPMP_singleCell_averageMinusLog10Pvalue_line_class>();
            this.Averaged_minuslog10p = this.Averaged_minuslog10p.OrderBy(l=>l.Ontology).ThenBy(l => l.KPMP_data_integration_term).ThenBy(l=>l.Scp).ToArray();
            List<double> current_averageMinusLog10p = new List<double>();
            List<double> current_percent_averageMinusLog10p = new List<double>();
            for (int indexAvg=0; indexAvg<average_length;indexAvg++)
            {
                average_line = this.Averaged_minuslog10p[indexAvg];
                if (  (indexAvg==0)
                    || (!average_line.Ontology.Equals(this.Averaged_minuslog10p[indexAvg - 1].Ontology))
                    || (!average_line.Scp.Equals(this.Averaged_minuslog10p[indexAvg - 1].Scp))
                    || (!average_line.KPMP_data_integration_term.Equals(this.Averaged_minuslog10p[indexAvg - 1].KPMP_data_integration_term)))
                {
                    current_averageMinusLog10p.Clear();
                    current_percent_averageMinusLog10p.Clear();
                }
                current_percent_averageMinusLog10p.Add(average_line.Percentage_averaged_minusLog10Pvalue);
                current_averageMinusLog10p.Add(average_line.Averaged_minusLog10Pvalue);
                if ((indexAvg == average_length-1)
                    || (!average_line.Ontology.Equals(this.Averaged_minuslog10p[indexAvg + 1].Ontology))
                    || (!average_line.Scp.Equals(this.Averaged_minuslog10p[indexAvg + 1].Scp))
                    || (!average_line.KPMP_data_integration_term.Equals(this.Averaged_minuslog10p[indexAvg + 1].KPMP_data_integration_term)))
                {
                    if (current_percent_averageMinusLog10p.Count > integrationTerm_centerCounts_dict[average_line.KPMP_data_integration_term]) { throw new Exception(); }
                    if (current_averageMinusLog10p.Count > integrationTerm_centerCounts_dict[average_line.KPMP_data_integration_term]) { throw new Exception(); }
                    while (current_averageMinusLog10p.Count< integrationTerm_centerCounts_dict[average_line.KPMP_data_integration_term])
                    {
                        current_percent_averageMinusLog10p.Add(0);
                        current_averageMinusLog10p.Add(0);
                    }
                    new_average_line = new KPMP_singleCell_averageMinusLog10Pvalue_line_class();
                    new_average_line.KPMP_data_integration_term = (string)average_line.KPMP_data_integration_term.Clone();
                    new_average_line.Center = (string)new_center_name.Clone();
                    new_average_line.Ontology = average_line.Ontology;
                    new_average_line.Scp = (string)average_line.Scp.Clone();
                    new_average_line.Averaged_minusLog10Pvalue = Math_class.Get_average(current_averageMinusLog10p.ToArray());
                    new_average_line.Averaged_entities_count = current_averageMinusLog10p.Count;
                    new_average_lines_list.Add(new_average_line);
                }
            }
            KPMP_singleCell_averageMinusLog10Pvalue_line_class[] new_average_lines = Calculate_percentage_of_minusLog10Pvalues_over_integration_term_and_center(new_average_lines_list.ToArray());
            Add_to_array(new_average_lines);
        }

        public void Write(string directory,string fileName)
        {
            KPMP_singleCell_averageMinusLog10Pvalue_readWriteOptions_class readWriteOptions = new KPMP_singleCell_averageMinusLog10Pvalue_readWriteOptions_class(directory, fileName);
            ReadWriteClass.WriteData(this.Averaged_minuslog10p, readWriteOptions);
        }

    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class KPMP_proteomic_exchangeGeneSymbol_documentation_line_class
    {
        public string Accession { get; set; }
        public string Dataset { get; set; }
        public string Comparison { get; set; }
        public string Full_gene_name_proteomicCore { get; set; }
        public string Full_gene_name_mssm { get; set; }
        public string Gene_symbol_proteomicCore { get; set; }
        public string Gene_symbol_mssm { get; set; }
        public string Final_gene_symbol { get; set; }
        public string Final_gene_symbol_origin { get; set; }
        public string Comment { get; set; }

        public KPMP_proteomic_exchangeGeneSymbol_documentation_line_class()
        {
            Accession = "";
            Dataset = "";
            Comparison = "";
            Full_gene_name_proteomicCore = "";
            Full_gene_name_mssm = "";
            Gene_symbol_proteomicCore = "";
            Gene_symbol_mssm = "";
            Final_gene_symbol = "";
            Final_gene_symbol_origin = "";
            Comment = "";
        }

        public KPMP_proteomic_exchangeGeneSymbol_documentation_line_class Deep_copy()
        {
            KPMP_proteomic_exchangeGeneSymbol_documentation_line_class copy = (KPMP_proteomic_exchangeGeneSymbol_documentation_line_class)this.MemberwiseClone();
            copy.Dataset = (string)this.Dataset.Clone();
            copy.Comparison = (string)this.Comparison.Clone();
            copy.Full_gene_name_proteomicCore = (string)this.Full_gene_name_proteomicCore.Clone();
            copy.Gene_symbol_proteomicCore = (string)this.Gene_symbol_proteomicCore.Clone();
            copy.Gene_symbol_mssm = (string)this.Gene_symbol_mssm.Clone();
            copy.Final_gene_symbol = (string)this.Final_gene_symbol.Clone();
            return copy;
        }
    }

    class KPMP_proteomic_exchangeGeneSymbol_documentation_readWriteOptions_class : ReadWriteOptions_base
    {
        public KPMP_proteomic_exchangeGeneSymbol_documentation_readWriteOptions_class(string fileName)
        {
            this.File = Global_directory_class.Results_directory + fileName;
            this.Key_propertyNames = new string[] { "Dataset", "Accession", "Full_gene_name_proteomicCore", "Gene_symbol_proteomicCore", "Full_gene_name_mssm", "Gene_symbol_mssm", "Final_gene_symbol", "Final_gene_symbol_origin","Comparison","Comment" };
            this.Key_columnNames = this.Key_propertyNames;
            this.File_has_headline = true;
            this.HeadlineDelimiters = new char[] { Global_class.Tab };
            this.LineDelimiters = new char[] { Global_class.Tab };
            this.Report = ReadWrite_report_enum.Report_main;
        }
    }

    class KPMP_proteomic_exchangeGeneSymbol_documentation_class
    {
        public KPMP_proteomic_exchangeGeneSymbol_documentation_line_class[] Proteomic_exchanges { get; set; }

        public KPMP_proteomic_exchangeGeneSymbol_documentation_class()
        {
            this.Proteomic_exchanges = new KPMP_proteomic_exchangeGeneSymbol_documentation_line_class[0];
        }

        private void Remove_duplicates()
        {
            List<KPMP_proteomic_exchangeGeneSymbol_documentation_line_class> keep = new List<KPMP_proteomic_exchangeGeneSymbol_documentation_line_class>();
            this.Proteomic_exchanges = this.Proteomic_exchanges.OrderBy(l => l.Dataset).ThenBy(l => l.Comparison).ThenBy(l => l.Final_gene_symbol).ThenBy(l => l.Full_gene_name_mssm).ThenBy(l => l.Full_gene_name_proteomicCore).ThenBy(l => l.Gene_symbol_mssm).ThenBy(l => l.Gene_symbol_proteomicCore).ToArray();
            int proteomic_exchanges_length = this.Proteomic_exchanges.Length;
            KPMP_proteomic_exchangeGeneSymbol_documentation_line_class documentation_line;
            for (int indexD=0; indexD<proteomic_exchanges_length;indexD++)
            {
                documentation_line = this.Proteomic_exchanges[indexD];
                if (  (indexD==0)
                    || (!documentation_line.Dataset.Equals(this.Proteomic_exchanges[indexD - 1].Dataset))
                    || (!documentation_line.Comparison.Equals(this.Proteomic_exchanges[indexD - 1].Comparison))
                    || (!documentation_line.Final_gene_symbol.Equals(this.Proteomic_exchanges[indexD - 1].Final_gene_symbol))
                    || (!documentation_line.Full_gene_name_mssm.Equals(this.Proteomic_exchanges[indexD - 1].Full_gene_name_mssm))
                    || (!documentation_line.Full_gene_name_proteomicCore.Equals(this.Proteomic_exchanges[indexD - 1].Full_gene_name_proteomicCore))
                    || (!documentation_line.Gene_symbol_mssm.Equals(this.Proteomic_exchanges[indexD - 1].Gene_symbol_mssm))
                    || (!documentation_line.Gene_symbol_proteomicCore.Equals(this.Proteomic_exchanges[indexD - 1].Gene_symbol_proteomicCore)))
                {
                    keep.Add(documentation_line);
                }
            }
            this.Proteomic_exchanges = keep.ToArray();
        }

        private void Add_deep_copy_to_array(KPMP_proteomic_exchangeGeneSymbol_documentation_line_class[] other_proteomic_exchanges)
        {
            int other_length = other_proteomic_exchanges.Length;
            int this_length = this.Proteomic_exchanges.Length;
            int new_length = this_length + other_length;
            KPMP_proteomic_exchangeGeneSymbol_documentation_line_class[] new_proteomic_exchanges = new KPMP_proteomic_exchangeGeneSymbol_documentation_line_class[new_length];
            int indexNew = -1;
            for (int indexThis = 0; indexThis < this_length; indexThis++)
            {
                indexNew++;
                new_proteomic_exchanges[indexNew] = this.Proteomic_exchanges[indexThis];
            }
            for (int indexOther = 0; indexOther < other_length; indexOther++)
            {
                indexNew++;
                new_proteomic_exchanges[indexNew] = other_proteomic_exchanges[indexOther];
            }
            this.Proteomic_exchanges = new_proteomic_exchanges.ToArray();
        }

        public void Generate_by_adding_to_array(KPMP_proteomic_exchangeGeneSymbol_documentation_line_class[] proteomic_exchanges)
        {
            Add_deep_copy_to_array(proteomic_exchanges);
            Remove_duplicates();
        }

        public void Write(string fileName)
        {
            KPMP_proteomic_exchangeGeneSymbol_documentation_readWriteOptions_class readWriteOptions = new KPMP_proteomic_exchangeGeneSymbol_documentation_readWriteOptions_class(fileName);
            ReadWriteClass.WriteData(this.Proteomic_exchanges, readWriteOptions);
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class KPMP_dataset_name_class
    {
        public const string Disease_singleCellRNAseq = "Disease SC RNAseq";
        public const string SingleCell_premiere = "SC RNASeq";//"SC RNASeq PREMIERE";
        public const string SingleCell_ucsf = "SC RNASeq II";//"SC RNASeq UCSF";
        public const string SingleNucleus_ucsd = "SN RNASeq";//"SN RNASeq UCSDWU";
        public const string Lmd_rnaseq_iuosu = "LMD RNASeq";//"LMD RNASeq OSUIU";
        public const string Lmd_proteomics_iuosu = "LMD Proteomics";//"LMD Proteomics OSUIU";
        public const string Codex_imaging_ucsf = "Codex imaging";//"Codex imaging UCSF";
        public const string NearSingleCell_proteomics_ucsf = "NSC Proteomics";//"NSC Proteomics UCSF";
        public const string Spatial_metabolomics = "Spatial metabolomics";//"Spatial metabolomics UTHSA-PNNL-EMBL";
        public const string Proteomics_datasetClass = "Proteomics";
        public const string RNASeq_datasetClass = "RNASeq";

        public static string[] Get_all_rnaSeq_datasets()
        {
            return new string[] { SingleCell_premiere, SingleCell_ucsf, SingleNucleus_ucsd, Lmd_rnaseq_iuosu };
        }

        public static string[] Get_all_proteomics_datasets()
        {
            return new string[] { Lmd_proteomics_iuosu, NearSingleCell_proteomics_ucsf };
        }

        public static string Get_datasetClass(string dataset)
        {
            if (Get_all_proteomics_datasets().Contains(dataset))
            {
                return Proteomics_datasetClass;
            }
            else if (Get_all_rnaSeq_datasets().Contains(dataset))
            {
                return RNASeq_datasetClass;
            }
            else
            {
                throw new Exception();
            }
        }

        public static string Get_bgGenes_fileName_for_dataset(string dataset)
        {
            return "BgGenesInUpperCase_" + dataset + ".txt";
        }
    }

    class KPMP_data_integration_class
    {
        public static string Get_podocyte_glomerulus_integration_term()
        {
            return "Podocyte_glomerulus";
        }

        public static string Get_proximal_tubule_integration_term()
        {
            return "Proximal_tubule";
        }

        public static string Get_not_assigned_label()
        {
            return "Not assigned";
        }

        public static string Rearrange_cellType_order_if_multiple_cellTypes_discovered(string multipleCellType_string)
        {
            int indexBackslash = multipleCellType_string.IndexOf('/');
            int indexForwardSlash = multipleCellType_string.IndexOf('\\');
            if (indexForwardSlash != -1) { throw new Exception(); }
            string returnCellType = (string)multipleCellType_string.Clone();
            if (indexBackslash != -1)
            {
                switch (multipleCellType_string)
                {
                    case "MC/vSMC/P":
                    case "vSMC/P/MC":
                        returnCellType = "MC/vSMC/P";
                        break;
                    case "DTL/vSMC/P":
                    case "vSMC/P/DTL":
                        returnCellType = "DTL/vSMC/P";
                        break;
                    case "Tcell/Bcell":
                    case "Bcell/Tcell":
                        returnCellType = "Bcell/Tcell";
                        break;
                    case "PC/IC":
                    case "IC/PC":
                        returnCellType = "PC/IC";
                        break;
                    case "INT/MC":
                    case "MC/INT":
                        returnCellType = "MC/INT";
                        break;
                    case "PT/DTL":
                    case "DTL/PT":
                        returnCellType = "PT/DTL";
                        break;
                    case "PT/ATL":
                    case "ATL/PT":
                        returnCellType = "PT/ATL";
                        break;
                    case "DCT/TAL":
                    case "TAL/DCT":
                        returnCellType = "TAL/DCT";
                        break;
                    case "DCT/ATL":
                    case "ATL/DCT":
                        returnCellType = "ATL/DCT";
                        break;
                    case "INT/DTL":
                    case "DTL/INT":
                        returnCellType = "DTL/INT";
                        break;
                    case "MON/MAC":
                    case "MAC/MON":
                        returnCellType = "MAC/MON";
                        break;
                    case "EC-AVR/EC-AEA-DVR":
                    case "EC-AEA-DVR/EC-AVR":
                        returnCellType = "EC-AVR/EC-AEA-DVR";
                        break;
                    case "DTL/EC-AEA-DVR":
                    case "EC-AEA-DVR/DTL":
                        returnCellType = "DTL/EC-AEA-DVR";
                        break;
                    case "vSMC/P/EC-AEA-DVR":
                    case "EC-AEA-DVR/vSMC/P":
                        returnCellType = "EC-AEA-DVR/vSMC/P";
                        break;
                    case "ATL/PC":
                    case "PC/ATL":
                        returnCellType = "ATL/PC";
                        break;
                    case "DCT/PC":
                    case "PC/DCT":
                        returnCellType = "DCT/PC";
                        break;
                    case "ATL/TAL":
                    case "TAL/ATL":
                        returnCellType = "ATL/TAL";
                        break;
                    case "MC/DTL":
                    case "DTL/MC":
                        returnCellType = "MC/DTL";
                        break;
                    case "TAL/IC":
                    case "IC/TAL":
                        returnCellType = "TAL/IC";
                        break;
                    case "DCT/IC":
                    case "IC/DCT":
                        returnCellType = "DCT/IC";
                        break;
                    case "PC/TAL":
                    case "TAL/PC":
                        returnCellType = "TAL/PC";
                        break;
                    case "PC/CNT":
                    case "CNT/PC":
                        returnCellType = "CNT/PC";
                        break;
                    case "DCT/CNT":
                    case "CNT/DCT":
                        returnCellType = "DCT/CNT";
                        break;
                    case "Bcell/MON":
                    case "MON/Bcell":
                        returnCellType = "Bcell/MON";
                        break;
                    case "Tcell/EC-AEA-DVR":
                    case "EC-AEA-DVR/Tcell":
                        returnCellType = "EC-AEA-DVR/Tcell";
                        break;
                    case "INT/EC-AEA-DVR":
                    case "EC-AEA-DVR/INT":
                        returnCellType = "EC-AEA-DVR/INT";
                        break;
                    case "POD/INT":
                    case "INT/POD":
                        returnCellType = "POD/INT";
                        break;
                    case "vSMC/P/INT":
                    case "INT/vSMC/P":
                        returnCellType = "INT/vSMC/P";////new
                        break;
                    case "PT/TAL":
                        returnCellType = "PT/TAL";////new
                        break;
                    case "DTL/TAL":
                        returnCellType = "DTL/TAL";////new
                        break;
                    case "EC-AEA-DVR/EC-GLO":
                        returnCellType = "EC-GLO/EC-AEA-DVR";////new
                        break;
                    case "DCT/PT":
                        returnCellType = "PT/DCT";////new
                        break;
                    case "EC-GLO/MC":
                        returnCellType = "MC/EC-GLO";////new
                        break;
                    case "EC-AVR/EC-GLO":
                        returnCellType = "EC-GLO/EC-AVR";////new
                        break;
                    case "CNT/DTL":
                        returnCellType = "DTL/CNT";////new
                        break;
                    case "TAL/vSMC/P":
                        returnCellType = "TAL/vSMC/P";////new
                        break;
                    case "CNT/PT":
                        returnCellType = "PT/CNT";////new
                        break;
                    case "CNT/TAL":
                        returnCellType = "TAL/CNT";////new
                        break;
                    case "DCT/POD":
                        returnCellType = "POD/DCT";////new
                        break;
                    case "POD/PT":
                        returnCellType = "POD/PT";////new
                        break;
                    case "PT/vSMC/P":
                        returnCellType = "PT/vSMC/P";////new
                        break;
                    case "CNT/IC":
                        returnCellType = "CNT/IC";////new
                        break;
                    case "DTC/IC":
                        returnCellType = "DTC/IC";////new
                        break;
                    case "DTL/PC":
                        returnCellType = "DTL/PC";////new
                        break;
                    case "ATL/CNT":
                        returnCellType = "ATL/CNT";////new
                        break;
                    case "PC/vSMC/P":
                        returnCellType = "PC/vSMC/P";////new
                        break;
                    case "DCT/DTL":
                        returnCellType = "DTL/DCT";////new
                        break;
                    case "EC-GLO/PT":
                        returnCellType = "PT/EC-GLO";////new
                        break;
                    case "MAC/TAL":
                        returnCellType = "TAL/MAC";////new
                        break;
                    case "INT/PT":
                        returnCellType = "PT/INT";////new
                        break;
                    case "CNT/POD":
                        returnCellType = "POD/CNT";////new
                        break;
                    case "EC-GLO/INT":
                        returnCellType = "INT/EC-GLO";////new
                        break;
                    case "PC/PT":
                        returnCellType = "PT/PC";////new
                        break;
                    case "CNT/INT":
                        returnCellType = "CNT/INT";////new
                        break;
                    case "CNT/EC-GLO":
                        returnCellType = "CNT/EC-GLO";////new
                        break;
                    case "ATL/vSMC/P":
                        returnCellType = "ATL/vSMC/P";////new
                        break;
                    case "CNT/vSMC/P":
                        returnCellType = "CNT/vSMC/P";////new
                        break;
                    case "DCT/EC-GLO":
                        returnCellType = "DCT/EC-GLO";////new
                        break;
                    case "EC-AVR/PT":
                        returnCellType = "PT/EC-AVR";////new
                        break;
                    case "DTL/IC":
                        returnCellType = "DTL/IC";////new
                        break;
                    case "vSMC/P":
                        returnCellType = "vSMC/P";
                        break;
                    default:
                        throw new Exception();
                }
            }
            if (multipleCellType_string.Equals(returnCellType)) { }
            else if (multipleCellType_string.IndexOf("vSMC/P") != -1)
            {
                string[] splitStrings = multipleCellType_string.Split('/');
                if (splitStrings[1].Equals("P"))
                {
                    if (!("vSMC/P/" + splitStrings[2]).Equals(returnCellType)) { throw new Exception(); }
                }
                else if (splitStrings[1].Equals("vSMC"))
                {
                    if (!(splitStrings[0] + "vSMC/P/").Equals(returnCellType)) { throw new Exception(); }
                }
                else { throw new Exception(); }
            }
            else
            {
                string[] splitStrings = multipleCellType_string.Split('/');
                if (!returnCellType.Equals(splitStrings[1] + "/" + splitStrings[0])) { throw new Exception(); }
            }
            return returnCellType;
        }

        public static string Get_kpmp_integration_term_for_pilot_data(string cell)
        {
            string kpmp_integration_term = "error";
            switch (cell)
            {
                case "Bulk":
                case "BULK":
                    kpmp_integration_term = "bulk";
                    break;
                case "1":
                case "Ratio GvsTI - Glom":  //Subsegmental proteomics (own labeling)
                case "RATIO GvsTI - Glom":  //Subsegmental proteomics (own labeling)
                case "Glom - POD":          //Subsegmental RNASeq (own labeling)
                case "Gloms":
                case "Glom":
                case "G":
                case "Podocyte":            //UMICH single cell  //UMICH 2
                case "podocytes":           //UCSF mdscRNASeq
                case "POD":                 //UCSD singleNucleusRNASeq: Podocytes
                case "POD_PTPRO":           //UCSD singleNucleusRNASeq: Podocytes
                    kpmp_integration_term = Get_podocyte_glomerulus_integration_term();
                    break;
                case "Parietal epithelial": //UMICH 2
                    kpmp_integration_term = "Parietal epithelial";
                    break;
                case "0":
                case "Ratio TIvsG - PT":    //Subsegmental proteomics
                case "RATIO TIvsG - PT":    //Subsegmental proteomics
                case "RATIO TIVSG - PT":    //Subsegmental proteomics
                case "PT":
                case "ProxTub":
                case "ProxTub-1":
                case "ProxTub-2":
                case "ProxTub-3":
                case "ProxTub-4":
                case "ProxTub-5":
                case "ProxTub-6":
                case "ProxTub-7":
                case "ProxTub-8":
                case "ProxTub-9":
                case "Proximal":            //UMICH single cell  //UMICH 2
                case "Proximal S2":            //UMICH single cell  //UMICH 2
                case "Proximal S1":            //UMICH single cell  //UMICH 2
                case "TI - PT":
                case "PT1":                  
                case "PT2":                  
                case "PT3":                  
                case "PT4":                  
                case "PT5":                  
                case "PT6":                 
                case "PT7":                 
                case "PT8":                  
                case "PT9_HighRibogenes ":                 
                case "PT-DTL":                  
                case "PT-1":                //UCSD singleNucleusRNASeq: Proximal Tubule Epithelial Cells (S1)
                case "PT_S1_SLC5A12":       //UCSD singleNucleusRNASeq: Proximal Tubule Epithelial Cells (S1)
                case "PT-2":                //UCSD singleNucleusRNASeq: Proximal Tubule Epithelial Cells (S2)
                case "PT_S2_ACSM3":         //UCSD singleNucleusRNASeq: Proximal Tubule Epithelial Cells (S2)
                case "PT-3":                //UCSD singleNucleusRNASeq: Proximal Tubule Epithelial Cells (Stress/Inflam)
                case "PT_S1_IL36B":         //UCSD singleNucleusRNASeq: Proximal Tubule Epithelial Cells (Stress/Inflam)
                case "PT-4":                //UCSD singleNucleusRNASeq: Proximal Tubule Epithelial Cells (Fibrinogen+ (S3))
                case "PT_S3_PDZK1IP1":      //UCSD singleNucleusRNASeq: Proximal Tubule Epithelial Cells (Fibrinogen+ (S3))
                case "PT-5":                //UCSD singleNucleusRNASeq: Proximal Tubule Epithelial Cells (S3)
                case "PT_S3_MT1G":          //UCSD singleNucleusRNASeq: Proximal Tubule Epithelial Cells (S3)
                case "Unk":                 //UCSD singleNucleusRNASeq: Unknown - Novel PT CFH+ Subpopulation (S2)
                case "PT_S2_LINC01435":     //UCSD singleNucleusRNASeq: Unknown - Novel PT CFH+ Subpopulation (S2) 
                case "PT-6":
                case "PT-7":
                case "PT-8":
                case "PT-9":
                case "PT-10":
                    kpmp_integration_term = Get_proximal_tubule_integration_term();
                    break;
                case "DL":                  //UCSD singleNucleusRNASeq: Descending limb
                case "DL_CPE":              //UCSD singleNucleusRNASeq: Descending limb
                case "DTL":                 //UMICH single cell 2 - PREMIERE
                case "DTL1":                 //AKI UMich
                case "DTL2":                 //AKI UMich
                case "DTL3":                 //AKI UMich
                case "DTL4_injured":                 //AKI UMich
                    kpmp_integration_term = "Descending limb";
                    break;
                case "TAL1":
                case "TAL2":
                case "TAL3":
                case "TAL4":
                case "TAL5":
                case "TAL6":
                case "TAL7_HighATPgenes":
                case "TAL-1":               //UCSD singleNucleusRNASeq: thick ascending limb
                case "TAL_KNG1":            //UCSD singleNucleusRNASeq: thick ascending limb
                case "TAL-2":               //UCSD singleNucleusRNASeq: thick ascending limb
                case "TAL-3":               //UCSD singleNucleusRNASeq: thick ascending limb
                case "TAL_CASR":            //UCSD singleNucleusRNASeq: thick ascending limb
                case "TAL":                 //UMICH single cell
                case "tAL":                 //UMICH single cell
                case "Thick ascending LOH": //UMICH single cell 2
                case "TI - TAL":
                    kpmp_integration_term = "Thick_ascending_limb";
                    break;
                case "ATL":
                case "ATL1":
                case "ATL2":
                case "ATL3_HighATPgenes":
                case "ATL4_Injured":
                case "ATL-1":               //UCSD singleNucleusRNASeq: Thin ascending limb
                case "ATL_AKR1B1":          //UCSD singleNucleusRNASeq: Thin ascending limb
                case "ATL-2":               //UCSD singleNucleusRNASeq: Thin ascending limb
                case "ATL_SERPINA1":        //UCSD singleNucleusRNASeq: Thin ascending limb
                case "ATL-3":               //UCSD singleNucleusRNASeq: Thin ascending limb
                case "ATL_CLDN10-AS1":      //UCSD singleNucleusRNASeq: Thin ascending limb
                case "dLOH/aLOH (thin aLOH)": //UMICH 2
                    kpmp_integration_term = "Thin_ascending_limb";
                    break;
                case "DCT":                 //UCSD singleNucleusRNASeq: Distal convoluted tubule
                case "DCT_SLC12A3":         //UCSD singleNucleusRNASeq: Distal convoluted tubule
                case "DCT1":
                case "DCT2":
                case "DCT3_HighMito ":
                    kpmp_integration_term = "Distal Convoluted Tubule";
                    break;
                case "PC1":
                case "PC2":
                case "PC3":
                case "PC4_HighMito":
                case "PC5_HighATPgenes":
                case "PC":                 //UCSD singleNucleusRNASeq
                case "PC-1":               //UCSD singleNucleusRNASeq: Collecting duct - Principal Cells (cortex)
                case "PC_PWRN1":           //UCSD singleNucleusRNASeq: Collecting duct - Principal Cells (cortex)
                case "PC-2":               //UCSD singleNucleusRNASeq: Collecting duct - Principal Cells - Stressed Dissoc Subset
                case "PC_FOSB":            //UCSD singleNucleusRNASeq: Collecting duct - Principal Cells - Stressed Dissoc Subset
                case "PC-3":               //UCSD singleNucleusRNASeq: Collecting duct - Principal Cells (medulla)
                case "PC_AQP3":            //UCSD singleNucleusRNASeq: Collecting duct - Principal Cells (medulla)
                case "PC-CNT":             //UMICH single cell 2 - PREMIERE
                case "CD":                 //Subsegmental RNASeq
                case "CD - PC":             //UCSF mdscRNASeq //Subsegmental RNASeq (own labeling)
                case "Principal cells":     //UMICH single cell  //UMICH single cell 2
                case "TI - PC":
                    kpmp_integration_term = "Principal cells";
                    break;
                case "ICA1":
                case "ICA2":
                case "ICA/ICB":
                case "IC-A1":               //UCSD singleNucleusRNASeq: Collecting duct - Intercalated cells Type A (cortex)
                case "ICA_CLNK":            //UCSD singleNucleusRNASeq: Collecting duct - Intercalated cells Type A (cortex)
                case "IC-A2":               //UCSD singleNucleusRNASeq: Collecting duct - Intercalated cells Type A (medulla)
                case "ICA_CALCA":           //UCSD singleNucleusRNASeq: Collecting duct - Intercalated cells Type A (medulla)
                case "ICB_SLC26A4":         //UCSD singleNucleusRNASeq
                case "IC-A":                //UCSD singleNucleusRNASeq: Collecting duct - Intercalated cells Type B
                case "IC-B":                //UCSD singleNucleusRNASeq: Collecting duct - Intercalated cells Type B
                case "ICB-SLC26A4":         //UCSD singleNucleusRNASeq: Collecting duct - Intercalated cells Type B
                case "Intercalated":        //UMICH single cell  //UMICH single cell 2
                case "IC":                  //UCSD singleNucleusRNASeq
                case "IC-1":
                case "IC-2":
                case "tPC-IC":              //UMICH single cell 2 - PREMIERE
                case "tPC_IC":              //UMICH single cell 2 - PREMIERE
                case "CD - IC":             //UCSF mdscRNASeq //Subsegmental RNASeq (own labeling)
                case "TI - IC":
                    kpmp_integration_term = "Intercalated cells";
                    break;
                case "TI":
                case "Ti":
                    kpmp_integration_term = "Tubulointerstitial_compartment";
                    break;
                case "EC_GC_EMCN":          //UCSD singleNucleusRNASeq: Endothelial Cells - glomerular capillaries
                case "EC-GLO":              //UMICH single cell
                case "GC-EC":              //UMICH single cell
                case "Glom - EC":           //Subsegmental RNASeq (own labeling)
                case "Glomerular endothelial":  //UMICH 2
                case "endothelial":         //UCSF mdscRNASeq
                case "endothelial  ":       //UCSF mdscRNASeq
                case "Ratio GvsTI - End - Glom":
                case "RATIO GvsTI - End - Glom":
                    kpmp_integration_term = "Endothelial glomerula";
                    break;
                case "EC1":
                case "EC2":
                case "EC3":
                case "EC4":
                case "EC-1":
                case "EC-2":                //UCSD singleNucleusRNASeq: Endothelial Cells - AVR
                case "EC-3":                //UCSD singleNucleusRNASeq: Endothelial Cells - AVR
                case "EC-4":                //UCSD singleNucleusRNASeq: Endothelial Cells - AVR
                case "EC_AVR_ADGRL4":       //UCSD singleNucleusRNASeq: Endothelial Cells - AVR
                case "EC_AEA_DVR_A2M":      //UCSD singleNucleusRNASeq: Endothelial Cells - AEA & DVR
                case "EC_MMRN1":            //UCSD singleNucleusRNASeq: Endothelial Cells (unassigned)
                case "EC-DVR":              //UMICH single cell
                case "EC-AVR":              //UMICH single cell
                case "EC-AEA-DVR":          //UCSD singleNucleusRNASeq
                case "EC-AEA":
                case "EC-PT":
                case "Peritubular endothelial": //UMICH 2
                case "Arteriolar endothelial": //UMICH 2
                    kpmp_integration_term = "Endothelial";
                    break;
                case "dLOH":                //UMICH single cell
                case "LOH - DL":            //UCSF mdscRNASeq
                case "short descending loop of Henle (sdLOH)":  //UMICH 2
                    kpmp_integration_term = "Descending Loop of Henle";
                    break;
                case "aLOH":                //UMICH single cell
                case "LOH - AL":            //UCSF mdscRNASeq
                    kpmp_integration_term = "Ascending Loop of Henle";
                    break;
                case "Interstitial":        //UMICH single cell
                case "Int":
                case "INT":                         //UCSD singleNucleusRNASeq: Interstitium
                case "INT_C7":                      //UCSD singleNucleusRNASeq: Interstitium
                case "Ratio TIvsG - Interstitium":  //Subsegmental proteomics
                case "RATIO TIvsG - Interstitium":  //Subsegmental proteomics
                case "RATIO TIVSG - Interstitium":  //Subsegmental proteomics
                case "RATIO TIVSG - INTERSTITIUM":  //Subsegmental proteomics
                case "Fibroblast":  //UMICH 2
                case "TI - Interstitium":
                case "FIB":                 //UMICH single cell 2 - PREMIERE
                case "INT - FIB":
                    kpmp_integration_term = "Interstitium";
                    break;
                case "Ratio GvsTI - Mesangial":    //Subsegmental proteomics
                case "RATIO GvsTI - Mesangial":    //Subsegmental proteomics
                case "Mesangial-VSMC":      //UMICH single cell
                case "Mesangial":           //UCSF mdscRNASeq
                case "MC":                  //UCSD singleNucleusRNASeq Mesangial cells
                case "MC_EBF1":             //UCSD singleNucleusRNASeq Mesangial cells
                case "Glom - Mesangial":    //Subsegmental RNASeq (own labeling)
                case "Vascular smooth muscle cells (VSMCs) /Mesangial": //UMICH single cells 2
                case "vSMC/MC":
                    kpmp_integration_term = "Mesangial";
                    break;
                case "vSMC1":
                case "vSMC2":
                case "vSMC/P":
                    kpmp_integration_term = "vSMC";
                    break;
                case "CNT1":
                case "CNT":                 //UMICH single cell: Connecting tubule
                case "CNT-DCT":             //UMICH single cell: Connecting tubule
                case "CNT-PC":
                case "tCNT_PC":
                case "Connecting /Distal tubule":  //UMICH single cell 2
                case "CNT_LINC01099":       //UCSD singleNucleusRNASeq
                    kpmp_integration_term = "Connecting tubule";
                    break;
                case "Macrophage":          //UMICH single cell
                case "MAC":                 //UCSD singleNucleusRNASeq 
                case "TI - MAC":
                case "IMM":                 //UCSD singleNucleusRNASeq: Immune Cells - Macrophages
                case "IMM_Mac_CD163":       //UCSD singleNucleusRNASeq: Immune Cells - Macrophages
                case "MON":                 //Monocyte //UMICH single cell 2 - PREMIERE
                case "Monocyte":            //UMICH single cell //UMICH single cell 2
                case "INT - MAC":
                    kpmp_integration_term = "Macrophage";
                    break;
                case "vSMC-P":              //UCSD singleNucleusRNASeq: Vascular Smooth Muscle Cells and pericytes
                case "vSMC-P_MYH11":        //UCSD singleNucleusRNASeq: Vascular Smooth Muscle Cells and pericytes
                case "VSMC/P":
                    kpmp_integration_term = "VSMC and pericytes";
                    break;
                case "AEA":
                    kpmp_integration_term = "AEA";
                    break;
                case "T cells":             //UMICH single cell //UCSF mdscRNASeq
                case "Tcells1":             //UMICH single cell //UCSF mdscRNASeq
                case "Tcells2":             //UMICH single cell //UCSF mdscRNASeq
                case "Tcells3":             //UMICH single cell //UCSF mdscRNASeq
                case "Tcells4":             //UMICH single cell //UCSF mdscRNASeq
                case "Tcells5_HighRibo":             //UMICH single cell //UCSF mdscRNASeq
                case "T cells (CD8 positive)": //UMICH single cell 2
                case "T cells (activated)": //UMICH 2
                case "T-CYT":               //UMICH single cell 2 - PREMIERE
                case "T Cell":              //UMICH single cell 2 - PREMIERE
                case "T-CYT-MEM":           //UMICH single cell 2 - PREMIERE
                case "Tcells":
                case "T-ACT":
                case "T-MEM":
                    kpmp_integration_term = "T cells";
                    break;
                case "Bcells1":             //UCSF mdscRNASeq //UMICH 2
                case "Bcells2":             //UCSF mdscRNASeq //UMICH 2
                case "Bcells3":             //UCSF mdscRNASeq //UMICH 2
                case "Bcells4":             //UCSF mdscRNASeq //UMICH 2
                case "B cells":             //UCSF mdscRNASeq //UMICH 2
                case "B cell":             //UCSF mdscRNASeq //UMICH 2
                case "B Cell":              //UMICH single cell 2 - PREMIERE
                case "Bcells":
                    kpmp_integration_term = "B cells";
                    break;
                case "Natural killer cells"://UMICH single cell
                case "Natural killer cells (NKC)":  //UMICH single cell 2
                case "Natural killer T cells": //UMICH 2
                case "NKC": //UMICH 2
                    kpmp_integration_term = "Natural killer cells";
                    break;
                case "EPC":                 //UCSD singleNucleusRNASeq: Epithelial cells (unassigned)
                case "PEC/LOH":                 //UCSD singleNucleusRNASeq: Epithelial cells (unassigned)
                case "EPC_SPATA22":         //UCSD singleNucleusRNASeq: Epithelial cells (unassigned)
                case "PEC":                 //AKI UMich
                case "Glom - EPC":
                    kpmp_integration_term = "Epithelial Cell";
                    break;
                case "Myeloid1":
                case "Myeloid2":
                case "Myeloid":
                    kpmp_integration_term = "Myeloid";
                    break;
                //case "endothelial and immune"://UCSF mdscRNASeq
                //case "immune":              //UCSF mdscRNASeq
                //case "maybe - glomerular (VCAM1) and tubular (GATM) maybe": //UCSF mdscRNASeq
                //case "RBCs":                //UCSF mdscRNASeq
                //case "MT":                  //UCSF mdscRNASeq
                //case "Ratio GvsTI":
                //case "Ratio TIvsG":
                //case "RATIO GvsTI":
                //case "RATIO TIvsG":
                //case "T-HP":
                case "PAX8positivecells":
                case "LOH/DCT/IC":
                  //  kpmp_integration_term = "no integration term";
                    break;
                default:
                    string cell_error = cell;
                    throw new Exception();
            }
            return kpmp_integration_term;
        }

        public static string Get_kpmp_integration_term_for_disease_data(string cell_input)
        {
            string[] char_in_front_of_cell_types_array = new string[] { "_", " " };
            string[] char_after_of_cell_types_array = new string[] { ".", "." };

            string char_in_front_of_cell_types;
            string char_after_cell_types;
            int delimiters_length = char_in_front_of_cell_types_array.Length;
            //if (indexVS!=-1) { cell = cell.Substring(0, indexVS); }
            string kpmp_integration_term = "error";
            bool integrationGroup_found = false;
            string cell = "_" + cell_input + ".";
            for (int indexD = 0; indexD < delimiters_length; indexD++)
            {
                char_in_front_of_cell_types = char_in_front_of_cell_types_array[indexD];
                char_after_cell_types = char_after_of_cell_types_array[indexD];
                if ((cell.IndexOf(char_in_front_of_cell_types + "POD" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "dPOD" + char_after_cell_types) != -1))
                {
                    if (integrationGroup_found) { throw new Exception(); }
                    integrationGroup_found = true;
                    kpmp_integration_term = "Podocyte";
                }
                if ((cell.IndexOf(char_in_front_of_cell_types + "MC" + char_after_cell_types) != -1))
                {
                    if (integrationGroup_found) { throw new Exception(); }
                    integrationGroup_found = true;
                    kpmp_integration_term = "Mesangial_cell";
                }
                if ((cell.IndexOf(char_in_front_of_cell_types + "PEC" + char_after_cell_types) != -1))
                {
                    if (integrationGroup_found) { throw new Exception(); }
                    integrationGroup_found = true;
                    kpmp_integration_term = "Parietal_epithelial_cell";
                }
                if ((cell.IndexOf(char_in_front_of_cell_types + "aPT" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "cycPT" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "PT-S3" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "PT-S1-2" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "dPT" + char_after_cell_types) != -1))
                {
                    if (integrationGroup_found) { throw new Exception(); }
                    integrationGroup_found = true;
                    kpmp_integration_term = "Proximal_tubule";
                }
                if ((cell.IndexOf(char_in_front_of_cell_types + "DTL1" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "DTL2" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "DTL3" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "dDTL" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "dDTL3" + char_after_cell_types) != -1))
                {
                    if (integrationGroup_found) { throw new Exception(); }
                    integrationGroup_found = true;
                    kpmp_integration_term = "Descending_limb";
                }
                if ((cell.IndexOf(char_in_front_of_cell_types + "ATL" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "dATL" + char_after_cell_types) != -1))
                {
                    if (integrationGroup_found) { throw new Exception(); }
                    integrationGroup_found = true;
                    kpmp_integration_term = "Thin_ascending_limb";
                }
                if ((cell.IndexOf(char_in_front_of_cell_types + "aTAL1" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "aTAL2" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "C-TAL" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "dC-TAL" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "dM-TAL" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "M-TAL" + char_after_cell_types) != -1))
                {
                    if (integrationGroup_found) { throw new Exception(); }
                    integrationGroup_found = true;
                    kpmp_integration_term = "Thick_ascending_limb";
                }
                if ((cell.IndexOf(char_in_front_of_cell_types + "DCT1" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "DCT2" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "dDCT" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "cycDCT" + char_after_cell_types) != -1))
                {
                    if (integrationGroup_found) { throw new Exception(); }
                    integrationGroup_found = true;
                    kpmp_integration_term = "Distal_convoluted_tubule";
                }
                if ((cell.IndexOf(char_in_front_of_cell_types + "CNT" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "dCNT" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "cycCNT" + char_after_cell_types) != -1))
                {
                    if (integrationGroup_found) { throw new Exception(); }
                    integrationGroup_found = true;
                    kpmp_integration_term = "Connecting_tubule";
                }
                if ((cell.IndexOf(char_in_front_of_cell_types + "CNT-PC" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "CCD-PC" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "dOMCD-PC" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "OMCD-PC" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "IMCD-PC" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "dIMCD-PC" + char_after_cell_types) != -1))
                {
                    if (integrationGroup_found) { throw new Exception(); }
                    integrationGroup_found = true;
                    kpmp_integration_term = "Principal_cell";
                }
                if ((cell.IndexOf(char_in_front_of_cell_types + "dC-IC-A" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "IC-B" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "OMCD-IC-A" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "CCD-IC-A" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "CCD-IC" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "CNT-IC" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "tPC-IC" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "CNT-IC-A" + char_after_cell_types) != -1))
                {
                    if (integrationGroup_found) { throw new Exception(); }
                    integrationGroup_found = true;
                    kpmp_integration_term = "Intercalated_cell";
                }
                if (  (cell.IndexOf(char_in_front_of_cell_types + "dIMCD" + char_after_cell_types) != -1)
                    ||(cell.IndexOf(char_in_front_of_cell_types + "IMCD" + char_after_cell_types) != -1))
                {
                    if (integrationGroup_found) { throw new Exception(); }
                    integrationGroup_found = true;
                    kpmp_integration_term = "Collecting_duct";
                }
                if ((cell.IndexOf(char_in_front_of_cell_types + "MD" + char_after_cell_types) != -1))
                {
                    if (integrationGroup_found) { throw new Exception(); }
                    integrationGroup_found = true;
                    kpmp_integration_term = "Macula_densa";
                }
                if ((cell.IndexOf(char_in_front_of_cell_types + "REN" + char_after_cell_types) != -1))
                {
                    if (integrationGroup_found) { throw new Exception(); }
                    integrationGroup_found = true;
                    kpmp_integration_term = "Juxtaglomerular";
                }
                if ((cell.IndexOf(char_in_front_of_cell_types + "FIB" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "aFIB" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "dFIB" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "M-FIB" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "MYOF" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "cycMYOF" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "dM-FIB" + char_after_cell_types) != -1))
                {
                    if (integrationGroup_found) { throw new Exception(); }
                    integrationGroup_found = true;
                    kpmp_integration_term = "Fibroblast";
                }
                if ((cell.IndexOf(char_in_front_of_cell_types + "EC-PTC" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "EC-GC" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "dEC" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "cycEC" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "EC-LYM" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "dEC-PTC" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "EC-AVR" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "EC-AEA-DVR" + char_after_cell_types) != -1))
                {
                    if (integrationGroup_found) { throw new Exception(); }
                    integrationGroup_found = true;
                    kpmp_integration_term = "Endothelial";
                }
                if ((cell.IndexOf(char_in_front_of_cell_types + "MAST" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "cDC" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "MAC-M2" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "PL" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "T" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "N" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "B" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "MDC" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "ncMON" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "cycMNP" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "NKC-T" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "pDC" + char_after_cell_types) != -1))
                {
                    if (integrationGroup_found) { throw new Exception(); }
                    integrationGroup_found = true;
                    kpmp_integration_term = "Immune_cell";
                }
                if ((cell.IndexOf(char_in_front_of_cell_types + "PapE" + char_after_cell_types) != -1))
                {
                    if (integrationGroup_found) { throw new Exception(); }
                    integrationGroup_found = true;
                    kpmp_integration_term = "Papilla";
                }
                if ((cell.IndexOf(char_in_front_of_cell_types + "VSMC" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "dVSMC" + char_after_cell_types) != -1)
                    || (cell.IndexOf(char_in_front_of_cell_types + "VSMC-P" + char_after_cell_types) != -1))
                {
                    if (integrationGroup_found) { throw new Exception(); }
                    integrationGroup_found = true;
                    kpmp_integration_term = "Vascular_smooth_muscle";
                }
                if (  (cell.IndexOf(char_in_front_of_cell_types + "SC/NEU" + char_after_cell_types) != -1)
                    ||(cell.IndexOf(char_in_front_of_cell_types + "SC-NEU" + char_after_cell_types) != -1))
                {
                    if (integrationGroup_found) { throw new Exception(); }
                    integrationGroup_found = true;
                    kpmp_integration_term = "Neuronal";
                }
            }
            if (!integrationGroup_found) { throw new Exception(); }
            return kpmp_integration_term;
        }

        public static string Get_kpmp_integration_term(string cell)
        {
            KPMP_analysis_set_enum analysis_set = Global_kpmp_class.Analysis_set;
            switch (analysis_set)
            {
                case KPMP_analysis_set_enum.Pilot_data:
                    return Get_kpmp_integration_term_for_pilot_data(cell);
                case KPMP_analysis_set_enum.Disease_data_2021:
                    return Get_kpmp_integration_term_for_disease_data(cell);
                default:
                    throw new Exception();
            }
        }

        public static string Get_kpmp_integration_term_plus_patientID(string cell, string patientID)
        {
            return Get_kpmp_integration_term(cell) + "-" + patientID.Replace('-', '_');
        }

        public static string Get_combined_patients_label()
        {
            return "allPatients";
        }

        public static string[] Get_names_for_de_in_correct_order(string cellSegment, string patientID, string dataset, KPMP_value_type_enum value_type, string data_integration_term)
        {
            string[] names_for_de = new string[] { (string)dataset.Clone() + "$" + (string)cellSegment.Clone() + "$" + (string)patientID.Clone() + "$" + value_type.ToString() + "$" + (string)data_integration_term.Clone() };
            return names_for_de;
        }

        public static string[] Get_names_for_de_in_correct_order_separated_by_minus(string cellSegment, string patientID, string dataset, KPMP_value_type_enum value_type, string data_integration_term)
        {
            string[] names_for_de = new string[] { (string)dataset.Clone() + "-" + (string)cellSegment.Clone() + "-" + (string)patientID.Clone() + "-" + value_type.ToString() + "-" + (string)data_integration_term.Clone() };
            return names_for_de;
        }

        public static string[] Get_names_for_de_in_correct_order_separated_by_x(string cellSegment, string patientID, string dataset, KPMP_value_type_enum value_type, string data_integration_term)
        {
            string[] names_for_de = new string[] { (string)dataset.Clone() + "x" + (string)cellSegment.Clone() + "x" + (string)patientID.Clone() + "x" + value_type.ToString() + "x" + (string)data_integration_term.Clone() };
            return names_for_de;
        }

    }

    class KPMP_subsegmental_rnaSeq_rttest_line_class
    {
        public string Dataset { get; set; }
        public string Segment { get; set; }
        public string Gene { get; set; }
        public double Pvalue { get; set; }
        public double MinusLog10Pvalue { get; set; }
        public double FoldChange { get; set; }
        public double Log2_fc { get; set; }

        public KPMP_subsegmental_rnaSeq_rttest_line_class()
        { }

        public KPMP_subsegmental_rnaSeq_rttest_line_class Deep_copy()
        {
            KPMP_subsegmental_rnaSeq_rttest_line_class copy = (KPMP_subsegmental_rnaSeq_rttest_line_class)this.MemberwiseClone();
            copy.Segment = (string)this.Segment.Clone();
            copy.Gene = (string)this.Gene.Clone();
            return copy;
        }
    }

    class KPMP_subsegmental_rnaSeq_rttest_readWriteOptions_class : ReadWriteOptions_base
    {
        public KPMP_subsegmental_rnaSeq_rttest_readWriteOptions_class(string subdirecty, string fileName)
        {
            this.File = Global_directory_class.Experimental_data_directory + subdirecty + fileName;
            this.Key_propertyNames = new string[] { "Segment", "Gene", "Pvalue", "MinusLog10Pvalue", "FoldChange", "Log2_fc" };
            this.Key_columnNames = this.Key_propertyNames;
            this.HeadlineDelimiters = new char[] { Global_class.Tab };
            this.LineDelimiters = new char[] { Global_class.Tab };
            this.File_has_headline = true;
            this.Report = ReadWrite_report_enum.Report_main;
        }
    }

    class KPMP_subsegmental_rnaSeq_rttest_class
    {
        public KPMP_subsegmental_rnaSeq_rttest_line_class[] Ttest_data { get; set; }

        private void Add_dataset(string dataset)
        {
            int ttest_data_length = this.Ttest_data.Length;
            KPMP_subsegmental_rnaSeq_rttest_line_class rnaSeq_line;
            for (int indexTTest=0; indexTTest<ttest_data_length;indexTTest++)
            {
                rnaSeq_line = this.Ttest_data[indexTTest];
                rnaSeq_line.Dataset = (string)dataset.Clone();
            }
        }

        public void Generate_by_reading(string subdirectory, string fileName, string dataset)
        {
            Read_routput(subdirectory, fileName);
            Add_dataset(dataset);
        }

        private void Read_routput(string subdirectory, string fileName)
        {
            KPMP_subsegmental_rnaSeq_rttest_readWriteOptions_class readWriteOptions = new KPMP_subsegmental_rnaSeq_rttest_readWriteOptions_class(subdirectory,fileName);
            this.Ttest_data = ReadWriteClass.ReadRawData_and_FillArray<KPMP_subsegmental_rnaSeq_rttest_line_class>(readWriteOptions);
        }
    }

    class KPMP_subsegmental_rnaSeq_line_class
    {
        public string Gene_symbol { get; set; }
        public string Sample_name { get; set; }
        public double Expression_value { get; set; }
        public string PatientId { get; set; }
        public string Subsegment { get; set; }
        public string Kpmp_data_integration_term { get; set; }
        public string Dataset { get; set; }
        public KPMP_value_type_enum Value_type { get; set; }

        #region Equal set
        public static KPMP_subsegmental_rnaSeq_line_class[] Order_by_set(KPMP_subsegmental_rnaSeq_line_class[] array)
        {
            array = array.OrderBy(l => l.Sample_name).ThenBy(l => l.PatientId).ThenBy(l => l.Subsegment).ThenBy(l => l.Kpmp_data_integration_term).ThenBy(l => l.Dataset).ThenBy(l => l.Value_type).ToArray();
            return array;
        }

        public static KPMP_subsegmental_rnaSeq_line_class[] Order_by_set_and_descending_expression_value(KPMP_subsegmental_rnaSeq_line_class[] array)
        {
            array = array.OrderBy(l => l.Sample_name).ThenBy(l => l.PatientId).ThenBy(l => l.Subsegment).ThenBy(l => l.Kpmp_data_integration_term).ThenBy(l => l.Dataset).ThenBy(l => l.Value_type).ThenByDescending(l => l.Expression_value).ToArray();
            return array;
        }

        public static KPMP_subsegmental_rnaSeq_line_class[] Order_by_set_and_expression_value(KPMP_subsegmental_rnaSeq_line_class[] array)
        {
            array = array.OrderBy(l => l.Sample_name).ThenBy(l => l.PatientId).ThenBy(l => l.Subsegment).ThenBy(l => l.Kpmp_data_integration_term).ThenBy(l => l.Dataset).ThenBy(l => l.Value_type).ThenBy(l => l.Expression_value).ToArray();
            return array;
        }

        public bool Equal_set(KPMP_subsegmental_rnaSeq_line_class other)
        {
            bool equal = (this.Sample_name.Equals(other.Sample_name))
                          && (this.PatientId.Equals(other.PatientId))
                          && (this.Subsegment.Equals(other.Subsegment))
                          && (this.Kpmp_data_integration_term.Equals(other.Kpmp_data_integration_term))
                          && (this.Dataset.Equals(other.Dataset))
                          && (this.Value_type.Equals(other.Value_type));
            return equal;
        }
        #endregion

        public static KPMP_subsegmental_rnaSeq_line_class[] Order_by_dataset_patientId_subsegment_geneSymbol_valueType(KPMP_subsegmental_rnaSeq_line_class[] data)
        {
            Dictionary<string, Dictionary<string, Dictionary<string, Dictionary<string, Dictionary<KPMP_value_type_enum, List<KPMP_subsegmental_rnaSeq_line_class>>>>>> dataset_patientID_subsegment_geneSymbol_valueType_dict = new Dictionary<string, Dictionary<string, Dictionary<string, Dictionary<string, Dictionary<KPMP_value_type_enum, List<KPMP_subsegmental_rnaSeq_line_class>>>>>>();
            Dictionary<string, Dictionary<string, Dictionary<string, Dictionary<KPMP_value_type_enum, List<KPMP_subsegmental_rnaSeq_line_class>>>>> patientID_subsegment_geneSymbol_valueType_dict = new Dictionary<string, Dictionary<string, Dictionary<string, Dictionary<KPMP_value_type_enum, List<KPMP_subsegmental_rnaSeq_line_class>>>>>();
            Dictionary<string, Dictionary<string, Dictionary<KPMP_value_type_enum, List<KPMP_subsegmental_rnaSeq_line_class>>>> subsegment_geneSymbol_valueType_dict = new Dictionary<string, Dictionary<string, Dictionary<KPMP_value_type_enum, List<KPMP_subsegmental_rnaSeq_line_class>>>>();
            Dictionary<string, Dictionary<KPMP_value_type_enum, List<KPMP_subsegmental_rnaSeq_line_class>>> geneSymbol_valueType_dict = new Dictionary<string, Dictionary<KPMP_value_type_enum, List<KPMP_subsegmental_rnaSeq_line_class>>>();
            Dictionary<KPMP_value_type_enum, List<KPMP_subsegmental_rnaSeq_line_class>> valueType_dict = new Dictionary<KPMP_value_type_enum, List<KPMP_subsegmental_rnaSeq_line_class>>();
            KPMP_subsegmental_rnaSeq_line_class subsegmental_rnaSeq_line;
            int data_length = data.Length;
            for (int indexData = 0; indexData<data_length;indexData++)
            {
                subsegmental_rnaSeq_line = data[indexData];
                if (!dataset_patientID_subsegment_geneSymbol_valueType_dict.ContainsKey(subsegmental_rnaSeq_line.Dataset))
                {
                    dataset_patientID_subsegment_geneSymbol_valueType_dict.Add(subsegmental_rnaSeq_line.Dataset, new Dictionary<string, Dictionary<string, Dictionary<string, Dictionary<KPMP_value_type_enum, List<KPMP_subsegmental_rnaSeq_line_class>>>>>());
                }
                if (!dataset_patientID_subsegment_geneSymbol_valueType_dict[subsegmental_rnaSeq_line.Dataset].ContainsKey(subsegmental_rnaSeq_line.PatientId))
                {
                    dataset_patientID_subsegment_geneSymbol_valueType_dict[subsegmental_rnaSeq_line.Dataset].Add(subsegmental_rnaSeq_line.PatientId, new Dictionary<string, Dictionary<string, Dictionary<KPMP_value_type_enum, List<KPMP_subsegmental_rnaSeq_line_class>>>>());
                }
                if (!dataset_patientID_subsegment_geneSymbol_valueType_dict[subsegmental_rnaSeq_line.Dataset][subsegmental_rnaSeq_line.PatientId].ContainsKey(subsegmental_rnaSeq_line.Subsegment))
                {
                    dataset_patientID_subsegment_geneSymbol_valueType_dict[subsegmental_rnaSeq_line.Dataset][subsegmental_rnaSeq_line.PatientId].Add(subsegmental_rnaSeq_line.Subsegment, new Dictionary<string, Dictionary<KPMP_value_type_enum, List<KPMP_subsegmental_rnaSeq_line_class>>>());
                }
                if (!dataset_patientID_subsegment_geneSymbol_valueType_dict[subsegmental_rnaSeq_line.Dataset][subsegmental_rnaSeq_line.PatientId][subsegmental_rnaSeq_line.Subsegment].ContainsKey(subsegmental_rnaSeq_line.Gene_symbol))
                {
                    dataset_patientID_subsegment_geneSymbol_valueType_dict[subsegmental_rnaSeq_line.Dataset][subsegmental_rnaSeq_line.PatientId][subsegmental_rnaSeq_line.Subsegment].Add(subsegmental_rnaSeq_line.Gene_symbol, new Dictionary<KPMP_value_type_enum, List<KPMP_subsegmental_rnaSeq_line_class>>());
                }
                if (!dataset_patientID_subsegment_geneSymbol_valueType_dict[subsegmental_rnaSeq_line.Dataset][subsegmental_rnaSeq_line.PatientId][subsegmental_rnaSeq_line.Subsegment][subsegmental_rnaSeq_line.Gene_symbol].ContainsKey(subsegmental_rnaSeq_line.Value_type))
                {
                    dataset_patientID_subsegment_geneSymbol_valueType_dict[subsegmental_rnaSeq_line.Dataset][subsegmental_rnaSeq_line.PatientId][subsegmental_rnaSeq_line.Subsegment][subsegmental_rnaSeq_line.Gene_symbol].Add(subsegmental_rnaSeq_line.Value_type, new List<KPMP_subsegmental_rnaSeq_line_class>());
                }
                dataset_patientID_subsegment_geneSymbol_valueType_dict[subsegmental_rnaSeq_line.Dataset][subsegmental_rnaSeq_line.PatientId][subsegmental_rnaSeq_line.Subsegment][subsegmental_rnaSeq_line.Gene_symbol][subsegmental_rnaSeq_line.Value_type].Add(subsegmental_rnaSeq_line);
            }

            string[] datasets = dataset_patientID_subsegment_geneSymbol_valueType_dict.Keys.ToArray();
            string dataset;
            int datasets_length = datasets.Length;
            string[] patientIDs;
            string patientID;
            int patientIDs_length;
            string[] subsegments;
            string subsegment;
            int subsegments_length;
            string[] geneSymbols;
            string geneSymbol;
            int geneSymbols_length;
            KPMP_value_type_enum[] valueTypes;
            KPMP_value_type_enum valueType;
            int valueTypes_length;
            List<KPMP_subsegmental_rnaSeq_line_class> ordered_data = new List<KPMP_subsegmental_rnaSeq_line_class>();
            datasets = datasets.OrderBy(l => l).ToArray();
            for (int indexData=0; indexData<datasets_length;indexData++)
            {
                dataset = datasets[indexData];
                patientID_subsegment_geneSymbol_valueType_dict = dataset_patientID_subsegment_geneSymbol_valueType_dict[dataset];
                patientIDs = patientID_subsegment_geneSymbol_valueType_dict.Keys.ToArray();
                patientIDs_length = patientIDs.Length;
                patientIDs = patientIDs.OrderBy(l => l).ToArray();
                for (int indexPatientID=0; indexPatientID<patientIDs_length;indexPatientID++)
                {
                    patientID = patientIDs[indexPatientID];
                    subsegment_geneSymbol_valueType_dict = patientID_subsegment_geneSymbol_valueType_dict[patientID];
                    subsegments = subsegment_geneSymbol_valueType_dict.Keys.ToArray();
                    subsegments_length = subsegments.Length;
                    subsegments = subsegments.OrderBy(l => l).ToArray();
                    for (int indexSubsegment=0; indexSubsegment<subsegments_length; indexSubsegment++)
                    {
                        subsegment = subsegments[indexSubsegment];
                        geneSymbol_valueType_dict = subsegment_geneSymbol_valueType_dict[subsegment];
                        geneSymbols = geneSymbol_valueType_dict.Keys.ToArray();
                        geneSymbols_length = geneSymbols.Length;
                        geneSymbols = geneSymbols.OrderBy(l => l).ToArray();
                        for (int indexGS =0; indexGS<geneSymbols_length;indexGS++)
                        {
                            geneSymbol = geneSymbols[indexGS];
                            valueType_dict = geneSymbol_valueType_dict[geneSymbol];
                            valueTypes = valueType_dict.Keys.ToArray();
                            valueTypes_length = valueTypes.Length;
                            valueTypes = valueTypes.OrderBy(l => l).ToArray();
                            for (int indexVT=0; indexVT<valueTypes_length;indexVT++)
                            {
                                valueType = valueTypes[indexVT];
                                ordered_data.AddRange(valueType_dict[valueType]);
                            }
                        }
                    }
                }
            }

            if (Global_class.Check_ordering)
            {
                int ordered_data_length = ordered_data.Count;
                if (ordered_data_length != data_length) { throw new Exception(); }
                KPMP_subsegmental_rnaSeq_line_class this_line;
                KPMP_subsegmental_rnaSeq_line_class previous_line;
                for (int indexO = 1; indexO < ordered_data_length; indexO++)
                {
                    previous_line = ordered_data[indexO - 1];
                    this_line = ordered_data[indexO];
                    if ((this_line.Dataset.CompareTo(previous_line.Dataset) < 0)) { throw new Exception(); }
                    else if ((this_line.Dataset.Equals(previous_line.Dataset))
                             && (this_line.PatientId.CompareTo(previous_line.PatientId) < 0)) { throw new Exception(); }
                    else if ((this_line.Dataset.Equals(previous_line.Dataset))
                             && (this_line.PatientId.Equals(previous_line.PatientId))
                             && (this_line.Subsegment.CompareTo(previous_line.Subsegment) < 0)) { throw new Exception(); }
                    else if ((this_line.Dataset.Equals(previous_line.Dataset))
                             && (this_line.PatientId.Equals(previous_line.PatientId))
                             && (this_line.Subsegment.Equals(previous_line.Subsegment))
                             && (this_line.Gene_symbol.CompareTo(previous_line.Gene_symbol) < 0)) { throw new Exception(); }
                    else if ((this_line.Dataset.Equals(previous_line.Dataset))
                             && (this_line.PatientId.Equals(previous_line.PatientId))
                             && (this_line.Subsegment.Equals(previous_line.Subsegment))
                             && (this_line.Gene_symbol.Equals(previous_line.Gene_symbol))
                             && (this_line.Value_type.CompareTo(previous_line.Value_type) < 0)) { throw new Exception(); }
                }
            }
            return ordered_data.ToArray();
        }


        public KPMP_subsegmental_rnaSeq_line_class()
        {
            this.Gene_symbol = "";
            this.Sample_name = "";
            this.Subsegment = "";
            this.PatientId = "";
            this.Kpmp_data_integration_term = "";
            this.Dataset = "";
        }

        public KPMP_subsegmental_rnaSeq_line_class Deep_copy()
        {
            KPMP_subsegmental_rnaSeq_line_class copy = (KPMP_subsegmental_rnaSeq_line_class)this.MemberwiseClone();
            copy.Gene_symbol = (string)this.Gene_symbol.Clone();
            copy.Sample_name = (string)this.Sample_name.Clone();
            copy.PatientId = (string)this.PatientId.Clone();
            copy.Kpmp_data_integration_term = (string)this.Kpmp_data_integration_term.Clone();
            return copy;
        }
    }

    class KPMP_subsegmental_rnaSeq_class
    {
        public KPMP_subsegmental_rnaSeq_line_class[] Data { get; set; }
        KPMP_subsegmental_rnaSeq_line_class[] Average_data { get; set; }
        public KPMP_integration_paper_metadata_class Dataset_patient { get; set; }
        string[] Bg_genes_in_upperCase { get; set; }

        public KPMP_subsegmental_rnaSeq_class()
        {
            this.Data = new KPMP_subsegmental_rnaSeq_line_class[0];
            this.Average_data = new KPMP_subsegmental_rnaSeq_line_class[0];
            this.Bg_genes_in_upperCase = new string[0];
            this.Dataset_patient = new KPMP_integration_paper_metadata_class();
        }

        //private void Calculate_average_data_over_sample_types()
        //{
        //    Data = Data.OrderBy(l => l.Sample_type).ThenBy(l=>l.Gene_symbol).ToArray();
        //    int data_length = Data.Length;
        //    KPMP_subsegmental_rnaSeq_line_class rnaSeq_line;
        //    KPMP_subsegmental_rnaSeq_line_class average_rnaSeq_line;
        //    List<KPMP_subsegmental_rnaSeq_line_class> average_list = new List<KPMP_subsegmental_rnaSeq_line_class>();
        //    List<float> current_symbol_sample_type_expression_values = new List<float>();
        //    for (int indexData=0; indexData<data_length; indexData++)
        //    {
        //        rnaSeq_line = Data[indexData];
        //        if ((indexData == 0)
        //            || (!rnaSeq_line.Gene_symbol.Equals(Data[indexData - 1].Gene_symbol))
        //            || (!rnaSeq_line.Sample_type.Equals(Data[indexData - 1].Sample_type)))
        //        {
        //            current_symbol_sample_type_expression_values.Clear();
        //        }
        //        current_symbol_sample_type_expression_values.Add(rnaSeq_line.Expression_value);
        //        if ((indexData == data_length-1)
        //            || (!rnaSeq_line.Gene_symbol.Equals(Data[indexData + 1].Gene_symbol))
        //            || (!rnaSeq_line.Sample_type.Equals(Data[indexData + 1].Sample_type)))
        //        {
        //            average_rnaSeq_line = new KPMP_tissueSubsection_rnaSeq_line_class();
        //            average_rnaSeq_line.Gene_symbol = (string)rnaSeq_line.Gene_symbol.Clone();
        //            average_rnaSeq_line.Expression_value = Math_class.Get_average(current_symbol_sample_type_expression_values.ToArray());
        //            average_rnaSeq_line.Sample_name = (string)rnaSeq_line.Sample_type.Clone();
        //            average_list.Add(average_rnaSeq_line);
        //        }
        //    }
        //    this.Average_data = average_list.ToArray();
        //}

        private void Add_to_array(KPMP_subsegmental_rnaSeq_line_class[] add_data)
        {
            int this_length = this.Data.Length;
            int add_length = add_data.Length;
            int new_length = this_length + add_length;
            int indexNew = -1;
            KPMP_subsegmental_rnaSeq_line_class[] new_data = new KPMP_subsegmental_rnaSeq_line_class[new_length];
            for (int indexThis = 0; indexThis < this_length; indexThis++)
            {
                indexNew++;
                new_data[indexNew] = this.Data[indexThis];
            }
            for (int indexAdd = 0; indexAdd < add_length; indexAdd++)
            {
                indexNew++;
                new_data[indexNew] = add_data[indexAdd];
            }
            this.Data = new_data;
        }

        private void Set_bg_genes_in_upperCase()
        {
            this.Data = this.Data.OrderBy(l => l.Gene_symbol).ToArray();
            int data_length = this.Data.Length;
            KPMP_subsegmental_rnaSeq_line_class subsegmental_line;
            List<string> bg_genes = new List<string>();
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                subsegmental_line = this.Data[indexData];
                if ((indexData == 0)
                    || (!subsegmental_line.Gene_symbol.Equals(this.Data[indexData - 1].Gene_symbol)))
                {
                    bg_genes.Add(subsegmental_line.Gene_symbol.ToUpper());
                }
            }
            this.Bg_genes_in_upperCase = bg_genes.Distinct().OrderBy(l => l).ToArray();
        }

        public string[] Get_deep_copy_of_bg_genes_in_upperCase()
        {
            int bg_genes_length = this.Bg_genes_in_upperCase.Length;
            string[] bg_genes = new string[bg_genes_length];
            for (int indexBG = 0; indexBG < bg_genes_length; indexBG++)
            {
                bg_genes[indexBG] = (string)this.Bg_genes_in_upperCase[indexBG].Clone();
            }
            return bg_genes;
        }

        private void Set_value_type_patientId_and_subsegments()
        {
            int data_length = this.Data.Length;
            KPMP_subsegmental_rnaSeq_line_class subsegmental_line;
            List<KPMP_subsegmental_rnaSeq_line_class> new_subsegmental_lines = new List<KPMP_subsegmental_rnaSeq_line_class>();
            string[] splitStrings;
            //string[] secondSplitStrings;
            List<string> patient_samples = new List<string>();
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                subsegmental_line = this.Data[indexData];
                splitStrings = subsegmental_line.Sample_name.Split('_');
                subsegmental_line.Subsegment = (string)splitStrings[0].Clone();
                subsegmental_line.PatientId = (string)splitStrings[1].Clone();
                subsegmental_line.Value_type = KPMP_value_type_enum.Single_value; 
                patient_samples.Add(subsegmental_line.Sample_name);
            }
            patient_samples = patient_samples.Distinct().OrderBy(l => l).ToList();
        }

        private void Set_kpmp_data_integration()
        {
            int data_length = this.Data.Length;
            KPMP_subsegmental_rnaSeq_line_class subsegmental_rnaSeq_line;
            for (int indexP = 0; indexP < data_length; indexP++)
            {
                subsegmental_rnaSeq_line = this.Data[indexP];
                subsegmental_rnaSeq_line.Kpmp_data_integration_term = KPMP_data_integration_class.Get_kpmp_integration_term_plus_patientID(subsegmental_rnaSeq_line.Subsegment, subsegmental_rnaSeq_line.PatientId);
            }
        }

        private void Duplicate_lines_add_add_cell_type_specificity()
        {
            int data_length = this.Data.Length;
            KPMP_subsegmental_rnaSeq_line_class data_line;
            KPMP_subsegmental_rnaSeq_line_class new_data_line;
            List<KPMP_subsegmental_rnaSeq_line_class> new_data_list = new List<KPMP_subsegmental_rnaSeq_line_class>();
            for (int indexD = 0; indexD < data_length; indexD++)
            {
                data_line = this.Data[indexD];
                if (data_line.Subsegment.Equals("CD"))
                {
                    new_data_line = data_line.Deep_copy();
                    new_data_line.Sample_name = new_data_line.Sample_name + " - IC";
                    new_data_line.Subsegment = new_data_line.Subsegment + " - IC";
                    new_data_list.Add(new_data_line);
                    data_line.Sample_name = data_line.Subsegment + " - PC";
                    data_line.Subsegment = data_line.Subsegment + " - PC";
                }
                else if (data_line.Subsegment.Equals("Glom"))
                {
                    new_data_line = data_line.Deep_copy();
                    new_data_line.Sample_name = new_data_line.Sample_name + " - POD";
                    new_data_line.Subsegment = new_data_line.Subsegment + " - POD";
                    new_data_list.Add(new_data_line);
                    new_data_line = data_line.Deep_copy();
                    new_data_line.Sample_name = new_data_line.Sample_name + " - EC";
                    new_data_line.Subsegment = new_data_line.Subsegment + " - EC";
                    new_data_list.Add(new_data_line);
                    new_data_line = data_line.Deep_copy();
                    new_data_line.Sample_name = new_data_line.Sample_name + " - EPC";
                    new_data_line.Subsegment = new_data_line.Subsegment + " - EPC";
                    new_data_list.Add(new_data_line);
                    data_line.Sample_name = data_line.Sample_name + " - Mesangial";
                    data_line.Subsegment = data_line.Subsegment + " - Mesangial";
                }
                else if (data_line.Subsegment.Equals("INT"))
                {
                    new_data_line = data_line.Deep_copy();
                    new_data_line.Sample_name = new_data_line.Sample_name + " - MAC";
                    new_data_line.Subsegment = new_data_line.Subsegment + " - MAC";
                    new_data_list.Add(new_data_line);
                    data_line.Sample_name = data_line.Sample_name + " - FIB";
                    data_line.Subsegment = data_line.Subsegment + " - FIB";
                }
            }
            Add_to_array(new_data_list.ToArray());
        }

        private void Remove_duplicated_patientId_valueType_subsegments_gene_symbol_throw_an_exception_if_anything_is_removed()
        {
            int data_length = this.Data.Length;
            KPMP_subsegmental_rnaSeq_line_class subsegmental_line;
            KPMP_subsegmental_rnaSeq_line_class inner_subsegmental_line;
            KPMP_subsegmental_rnaSeq_line_class inner_subsegmental_line2;
            List<KPMP_subsegmental_rnaSeq_line_class> keep = new List<KPMP_subsegmental_rnaSeq_line_class>();
            List<KPMP_subsegmental_rnaSeq_line_class> remove = new List<KPMP_subsegmental_rnaSeq_line_class>();
            this.Data = this.Data.OrderBy(l => l.Dataset).ThenBy(l => l.PatientId).ThenBy(l => l.Value_type).ThenBy(l => l.Subsegment).ThenBy(l => l.Gene_symbol).ThenBy(l => l.Sample_name).ToArray();
            int firstIndex_same_patientID_subsegment_geneSymbol = -1;
            for (int indexD = 0; indexD < data_length; indexD++)
            {
                subsegmental_line = this.Data[indexD];
                if ((indexD == 0)
                    || (!subsegmental_line.Dataset.Equals(this.Data[indexD - 1].Dataset))
                    || (!subsegmental_line.PatientId.Equals(this.Data[indexD - 1].PatientId))
                    || (!subsegmental_line.Value_type.Equals(this.Data[indexD - 1].Value_type))
                    || (!subsegmental_line.Subsegment.Equals(this.Data[indexD - 1].Subsegment))
                    || (!subsegmental_line.Gene_symbol.Equals(this.Data[indexD - 1].Gene_symbol)))
                {
                    firstIndex_same_patientID_subsegment_geneSymbol = indexD;
                }
                if ((indexD == data_length - 1)
                    || (!subsegmental_line.Dataset.Equals(this.Data[indexD + 1].Dataset))
                    || (!subsegmental_line.PatientId.Equals(this.Data[indexD + 1].PatientId))
                    || (!subsegmental_line.Value_type.Equals(this.Data[indexD + 1].Value_type))
                    || (!subsegmental_line.Subsegment.Equals(this.Data[indexD + 1].Subsegment))
                    || (!subsegmental_line.Gene_symbol.Equals(this.Data[indexD + 1].Gene_symbol)))
                {
                    if ((firstIndex_same_patientID_subsegment_geneSymbol != indexD)
                        && (subsegmental_line.PatientId.Equals(KPMP_data_integration_class.Get_combined_patients_label())))
                    {
                        throw new Exception();
                    }
                    for (int indexInner = firstIndex_same_patientID_subsegment_geneSymbol; indexInner <= indexD; indexInner++)
                    {
                        inner_subsegmental_line = this.Data[indexInner];
                        if ((indexInner != firstIndex_same_patientID_subsegment_geneSymbol)
                            && ((!inner_subsegmental_line.Sample_name.Equals(this.Data[indexInner - 1].Sample_name))
                               || (!inner_subsegmental_line.Expression_value.Equals(this.Data[indexInner - 1].Expression_value))))
                        {
                            inner_subsegmental_line2 = Data[indexInner - 1];
                            if (inner_subsegmental_line.PatientId.Equals(KPMP_data_integration_class.Get_combined_patients_label()))
                            {
                                throw new Exception();
                            }
                            //Console.WriteLine();
                            //throw new Exception();
                        }
                        if ((indexInner != firstIndex_same_patientID_subsegment_geneSymbol)
                            && ((inner_subsegmental_line.Sample_name.Equals(this.Data[indexInner - 1].Sample_name))
                                && (!inner_subsegmental_line.Expression_value.Equals(this.Data[indexInner - 1].Expression_value))))
                        {
                            inner_subsegmental_line2 = Data[indexInner - 1];
                            if (inner_subsegmental_line.PatientId.Equals(KPMP_data_integration_class.Get_combined_patients_label()))
                            {
                                throw new Exception();
                            }
                            //Console.WriteLine();
                            //throw new Exception();
                        }
                    }
                    keep.Add(subsegmental_line);
                }
                else
                {
                    throw new Exception();
                }
            }
            this.Data = keep.ToArray();
        }

        private void Remove_bulk_and_ti_segments()
        {
            List<KPMP_subsegmental_rnaSeq_line_class> keep = new List<KPMP_subsegmental_rnaSeq_line_class>();
            int data_length = this.Data.Length;
            KPMP_subsegmental_rnaSeq_line_class data_line;
            for (int indexD=0; indexD<data_length;indexD++)
            {
                data_line = this.Data[indexD];
                switch (data_line.Subsegment)
                {
                    case "Bulk":
                    case "TI":
                    case "Ti":
                        break;
                    case "Glom":
                    case "ProxTub":
                    case "TAL":
                    case "DCT":
                    case "CD":
                    case "INT":
                        keep.Add(data_line);
                        break;
                    default:
                        throw new Exception();
                }
            }
            this.Data = keep.ToArray();
        }

        private void Check_for_duplicates_for_de_instance()
        {
            int data_length = this.Data.Length;
            KPMP_subsegmental_rnaSeq_line_class subsegmental_line;
            this.Data = this.Data.OrderBy(l => l.Dataset).ThenBy(l => l.Subsegment).ThenBy(l => l.Value_type).ThenBy(l => l.PatientId).ThenBy(l => l.Kpmp_data_integration_term).ThenBy(l => l.Gene_symbol).ToArray();
            for (int indexD = 0; indexD < data_length; indexD++)
            {
                subsegmental_line = this.Data[indexD];
                if ((indexD != 0)
                    && (subsegmental_line.Dataset.Equals(this.Data[indexD - 1].Dataset))
                    && (subsegmental_line.Subsegment.Equals(this.Data[indexD - 1].Subsegment))
                    && (subsegmental_line.Value_type.Equals(this.Data[indexD - 1].Value_type))
                    && (subsegmental_line.PatientId.Equals(this.Data[indexD - 1].PatientId))
                    && (subsegmental_line.Kpmp_data_integration_term.Equals(this.Data[indexD - 1].Kpmp_data_integration_term))
                    && (subsegmental_line.Gene_symbol.Equals(this.Data[indexD - 1].Gene_symbol)))
                {
                    throw new Exception();
                }
            }
        }

        private void Set_all_geneSymbols_to_upperCase()
        {
            foreach (KPMP_subsegmental_rnaSeq_line_class subsegmental_line in this.Data)
            {
                subsegmental_line.Gene_symbol = subsegmental_line.Gene_symbol.ToUpper();
            }
        }

        private void Remove_all_genes_that_contain_only_NaN_minusLog10pvalues_after_setting_of_bg_genes_and_throw_an_exception_if_you_do_so()
        {
            Dictionary<string, bool> keep_genes = new Dictionary<string, bool>();
            int data_length = this.Data.Length;
            KPMP_subsegmental_rnaSeq_line_class lmd_rnaSeq_line;
            for (int indexData=0; indexData<data_length;indexData++)
            {
                lmd_rnaSeq_line = this.Data[indexData];
                if (lmd_rnaSeq_line.Value_type.Equals(KPMP_value_type_enum.Minus_log10_pvalue))
                {
                    if (!Double.IsNaN(lmd_rnaSeq_line.Expression_value))
                    {
                        if (!keep_genes.ContainsKey(lmd_rnaSeq_line.Gene_symbol))
                        {
                            keep_genes.Add(lmd_rnaSeq_line.Gene_symbol, true);
                        }
                    }
                }
            }
            List<KPMP_subsegmental_rnaSeq_line_class> keep = new List<KPMP_subsegmental_rnaSeq_line_class>();
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                lmd_rnaSeq_line = this.Data[indexData];
                if (keep_genes.ContainsKey(lmd_rnaSeq_line.Gene_symbol))
                {
                    keep.Add(lmd_rnaSeq_line);
                }
            }
            if (keep.Count != data_length) { throw new Exception(); }
            this.Data = keep.ToArray();
        }

        private void Remove_lines_with_negative_log2foldchanges()
        {
            int data_length = this.Data.Length;
            KPMP_subsegmental_rnaSeq_line_class subsegmental_rnaSeq_line;
            KPMP_subsegmental_rnaSeq_line_class inner_subsegmental_rnaSeq_line;
            List<KPMP_subsegmental_rnaSeq_line_class> keep = new List<KPMP_subsegmental_rnaSeq_line_class>();
            this.Data = this.Data.OrderBy(l => l.Dataset).ThenBy(l => l.PatientId).ThenBy(l => l.Subsegment).ThenBy(l => l.Gene_symbol).ThenBy(l => l.Expression_value).ToArray();
            int notset0_keep1_remove2 = 0;
            int indexFirstSameGroup = -1;
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                subsegmental_rnaSeq_line = this.Data[indexData];
                if ((indexData == 0)
                    || (!subsegmental_rnaSeq_line.Dataset.Equals(this.Data[indexData - 1].Dataset))
                    || (!subsegmental_rnaSeq_line.Gene_symbol.Equals(this.Data[indexData - 1].Gene_symbol))
                    || (!subsegmental_rnaSeq_line.Subsegment.Equals(this.Data[indexData - 1].Subsegment))
                    || (!subsegmental_rnaSeq_line.PatientId.Equals(this.Data[indexData - 1].PatientId)))
                {
                    indexFirstSameGroup = indexData;
                    notset0_keep1_remove2 = 0;
                }
                switch (subsegmental_rnaSeq_line.Value_type)
                {
                    case KPMP_value_type_enum.Log2_ratioavg:
                        if (subsegmental_rnaSeq_line.Expression_value > 0) { notset0_keep1_remove2 = 1; }
                        else { notset0_keep1_remove2 = 2; }
                        break;
                    case KPMP_value_type_enum.Minus_log10_pvalue:
                    case KPMP_value_type_enum.Minus_log10_pvalue_adjusted:
                        break;
                    default:
                        throw new Exception();
                }
                if ((indexData == data_length - 1)
                    || (!subsegmental_rnaSeq_line.Dataset.Equals(this.Data[indexData + 1].Dataset))
                    || (!subsegmental_rnaSeq_line.Gene_symbol.Equals(this.Data[indexData + 1].Gene_symbol))
                    || (!subsegmental_rnaSeq_line.Subsegment.Equals(this.Data[indexData + 1].Subsegment))
                    || (!subsegmental_rnaSeq_line.PatientId.Equals(this.Data[indexData + 1].PatientId)))
                {
                    switch (notset0_keep1_remove2)
                    {
                        case 1:
                            for (int indexInner = indexFirstSameGroup; indexInner <= indexData; indexInner++)
                            {
                                inner_subsegmental_rnaSeq_line = this.Data[indexInner];
                                keep.Add(inner_subsegmental_rnaSeq_line);
                            }
                            break;
                        case 2:
                            break;
                        case 0:
                        default:
                            throw new Exception();
                    }
                }
            }
            this.Data = keep.ToArray();
        }

        private void Add_2_to_all_singleValues()
        {
            int data_length = this.Data.Length;
            KPMP_subsegmental_rnaSeq_line_class data_line;
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                data_line = this.Data[indexData];
                switch (data_line.Value_type)
                {
                    case KPMP_value_type_enum.Single_value:
                        data_line.Expression_value += 2;
                        break;
                    default:
                        throw new Exception();
                }
            }
        }

        private void Generate_by_reading_r_output_private(string subdirectory, string fileName, string dataset)
        {
            KPMP_subsegmental_rnaSeq_rttest_class rttest = new KPMP_subsegmental_rnaSeq_rttest_class();
            rttest.Generate_by_reading(subdirectory, fileName, dataset);

            KPMP_subsegmental_rnaSeq_rttest_line_class rnaSeq_rttest_line;

            int ttest_length = rttest.Ttest_data.Length;
            KPMP_subsegmental_rnaSeq_line_class rnaSeq_line;
            List<KPMP_subsegmental_rnaSeq_line_class> add = new List<KPMP_subsegmental_rnaSeq_line_class>();
            for (int indexTtest=0; indexTtest<ttest_length;indexTtest++)
            {
                rnaSeq_rttest_line = rttest.Ttest_data[indexTtest];
                rnaSeq_line = new KPMP_subsegmental_rnaSeq_line_class();
                rnaSeq_line.Dataset = (string)rnaSeq_rttest_line.Dataset.Clone();
                rnaSeq_line.Expression_value = rnaSeq_rttest_line.Log2_fc;
                rnaSeq_line.Gene_symbol = (string)rnaSeq_rttest_line.Gene.Clone();
                rnaSeq_line.PatientId = KPMP_data_integration_class.Get_combined_patients_label();
                rnaSeq_line.Sample_name = KPMP_data_integration_class.Get_combined_patients_label();
                rnaSeq_line.Subsegment = (string)rnaSeq_rttest_line.Segment.Clone();
                rnaSeq_line.Value_type = KPMP_value_type_enum.Log2_ratioavg;
                add.Add(rnaSeq_line);

                rnaSeq_line = new KPMP_subsegmental_rnaSeq_line_class();
                rnaSeq_line.Dataset = (string)rnaSeq_rttest_line.Dataset.Clone();
                rnaSeq_line.Expression_value = rnaSeq_rttest_line.MinusLog10Pvalue;
                rnaSeq_line.Gene_symbol = (string)rnaSeq_rttest_line.Gene.Clone();
                rnaSeq_line.PatientId = KPMP_data_integration_class.Get_combined_patients_label();
                rnaSeq_line.Sample_name = KPMP_data_integration_class.Get_combined_patients_label();
                rnaSeq_line.Subsegment = (string)rnaSeq_rttest_line.Segment.Clone();
                rnaSeq_line.Value_type = KPMP_value_type_enum.Minus_log10_pvalue;
                add.Add(rnaSeq_line);
            }
            this.Data = add.ToArray();
        }

        public void Generate_by_reading_r_output(string subdirectory, string fileName, string dataset)
        {
            Generate_by_reading_r_output_private(subdirectory, fileName, dataset);
            Set_bg_genes_in_upperCase();
            Remove_all_genes_that_contain_only_NaN_minusLog10pvalues_after_setting_of_bg_genes_and_throw_an_exception_if_you_do_so();
            Keep_only_lines_with_indicated_patientIDs(KPMP_data_integration_class.Get_combined_patients_label());
            Remove_lines_with_negative_log2foldchanges();
            Set_all_geneSymbols_to_upperCase();
            Remove_duplicated_patientId_valueType_subsegments_gene_symbol_throw_an_exception_if_anything_is_removed();
            Duplicate_lines_add_add_cell_type_specificity();
            Set_kpmp_data_integration();
            Check_for_duplicates_for_de_instance();
        }

        public KPMP_standardized_dataset_class Generate_standardized_dataset_instance(KPMP_value_type_enum value_type_1st, KPMP_value_type_enum value_type_2nd)
        {
            KPMP_standardized_dataset_line_class standardized_data_line;
            List<KPMP_standardized_dataset_line_class> standardized_data_list = new List<KPMP_standardized_dataset_line_class>();
            int data_length = this.Data.Length;
            KPMP_subsegmental_rnaSeq_line_class subsegmental_line;
            this.Data = this.Data.OrderBy(l => l.Dataset).ThenBy(l => l.Subsegment).ThenBy(l => l.PatientId).ThenBy(l => l.Kpmp_data_integration_term).ThenBy(l => l.Gene_symbol).ThenBy(l => l.Value_type).ToArray();
            string dataset = (string)this.Data[0].Dataset.Clone();
            double current_value_1st = -1;
            double current_value_2nd = -1;
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                subsegmental_line = this.Data[indexData];
                if (!subsegmental_line.Dataset.Equals(dataset)) { throw new Exception(); }
                if ((indexData != 0)
                    && (subsegmental_line.Dataset.Equals(this.Data[indexData - 1].Dataset))
                    && (subsegmental_line.Subsegment.Equals(this.Data[indexData - 1].Subsegment))
                    && (subsegmental_line.PatientId.Equals(this.Data[indexData - 1].PatientId))
                    && (subsegmental_line.Kpmp_data_integration_term.Equals(this.Data[indexData - 1].Kpmp_data_integration_term))
                    && (subsegmental_line.Gene_symbol.Equals(this.Data[indexData - 1].Gene_symbol))
                    && (subsegmental_line.Value_type.Equals(this.Data[indexData - 1].Value_type)))
                {
                    throw new Exception();
                }
                if ((indexData == 0)
                    || (!subsegmental_line.Dataset.Equals(this.Data[indexData - 1].Dataset))
                    || (!subsegmental_line.Subsegment.Equals(this.Data[indexData - 1].Subsegment))
                    || (!subsegmental_line.PatientId.Equals(this.Data[indexData - 1].PatientId))
                    || (!subsegmental_line.Kpmp_data_integration_term.Equals(this.Data[indexData - 1].Kpmp_data_integration_term))
                    || (!subsegmental_line.Gene_symbol.Equals(this.Data[indexData - 1].Gene_symbol)))
                {
                    current_value_1st = -1;
                    current_value_2nd = -1;
                }
                if (subsegmental_line.Value_type.Equals(value_type_1st))
                {
                    if (current_value_1st != -1) { throw new Exception(); }
                    current_value_1st = subsegmental_line.Expression_value;
                }
                if (subsegmental_line.Value_type.Equals(value_type_2nd))
                {
                    if (current_value_2nd != -1) { throw new Exception(); }
                    current_value_2nd = subsegmental_line.Expression_value;
                }
                if ((indexData == data_length - 1)
                    || (!subsegmental_line.Dataset.Equals(this.Data[indexData + 1].Dataset))
                    || (!subsegmental_line.Subsegment.Equals(this.Data[indexData + 1].Subsegment))
                    || (!subsegmental_line.PatientId.Equals(this.Data[indexData + 1].PatientId))
                    || (!subsegmental_line.Kpmp_data_integration_term.Equals(this.Data[indexData + 1].Kpmp_data_integration_term))
                    || (!subsegmental_line.Gene_symbol.Equals(this.Data[indexData + 1].Gene_symbol)))
                {
                    if (current_value_1st == -1) { throw new Exception(); }
                    if (current_value_2nd == -1) { throw new Exception(); }
                    standardized_data_line = new KPMP_standardized_dataset_line_class();
                    standardized_data_line.Cell_segment = (string)subsegmental_line.Subsegment.Clone();
                    standardized_data_line.Dataset = (string)subsegmental_line.Dataset.Clone();
                    standardized_data_line.Value_1st = current_value_1st;
                    standardized_data_line.Value_2nd = current_value_2nd;
                    standardized_data_line.KPMP_data_integration_term = (string)subsegmental_line.Kpmp_data_integration_term.Clone();
                    standardized_data_line.Value_type_1st = value_type_1st;
                    standardized_data_line.Value_type_2nd = value_type_2nd;
                    standardized_data_line.PatientId = (string)subsegmental_line.PatientId.Clone();
                    standardized_data_line.Gene_symbol = (string)subsegmental_line.Gene_symbol.Clone();
                    standardized_data_list.Add(standardized_data_line);
                }
            }
            Dictionary<string, string[]> dataset_bgProteins_dict = new Dictionary<string, string[]>();
            dataset_bgProteins_dict.Add(dataset, Get_deep_copy_of_bg_genes_in_upperCase());
            KPMP_standardized_dataset_class standardized_data = new KPMP_standardized_dataset_class();
            standardized_data.Add_to_existing_instances(standardized_data_list.ToArray(), dataset_bgProteins_dict);
            return standardized_data;

            #region Old commented code
            //Report_class.Write("{0}: Generate standardized dataset instance:", typeof(KPMP_subsegmental_rnaSeq_class).Name);
            //int data_length = this.Data.Length;
            //KPMP_subsegmental_rnaSeq_line_class subsegmental_line;
            //KPMP_standardized_dataset_line_class standardized_data_line;
            //List<KPMP_standardized_dataset_line_class> standardized_data_list = new List<KPMP_standardized_dataset_line_class>();
            //this.Data = this.Data.OrderBy(l => l.Value_type).ThenBy(l => l.Kpmp_data_integration_term).ThenBy(l => l.PatientId).ThenBy(l => l.Subsegment).ThenBy(l => l.Dataset).ThenBy(l => l.Gene_symbol).ToArray();
            //List<KPMP_value_type_enum> value_types = new List<KPMP_value_type_enum>();
            //string dataset = (string)this.Data[0].Dataset.Clone();
            //for (int indexD = 0; indexD < data_length; indexD++)
            //{
            //    subsegmental_line = this.Data[indexD];
            //    if (!subsegmental_line.Dataset.Equals(dataset)) { throw new Exception(); }
            //    if ((indexD != 0)
            //        && (subsegmental_line.Dataset.Equals(this.Data[indexD - 1].Dataset))
            //        && (subsegmental_line.PatientId.Equals(this.Data[indexD - 1].PatientId))
            //        && (subsegmental_line.Subsegment.Equals(this.Data[indexD - 1].Subsegment))
            //        && (subsegmental_line.Kpmp_data_integration_term.Equals(this.Data[indexD - 1].Kpmp_data_integration_term))
            //        && (subsegmental_line.Gene_symbol.Equals(this.Data[indexD - 1].Gene_symbol))
            //        && (subsegmental_line.Value_type.Equals(this.Data[indexD - 1].Value_type)))
            //    {
            //        throw new Exception();
            //    }
            //    if ((indexD == 0) || (!subsegmental_line.Value_type.Equals(this.Data[indexD - 1].Value_type)))
            //    {
            //        value_types.Add(subsegmental_line.Value_type);
            //    }
            //    standardized_data_line = new KPMP_standardized_dataset_line_class();
            //    standardized_data_line.Dataset = (string)subsegmental_line.Dataset.Clone();
            //    standardized_data_line.PatientId = (string)subsegmental_line.PatientId.Clone();
            //    standardized_data_line.Cell_segment = (string)subsegmental_line.Subsegment.Clone();
            //    standardized_data_line.KPMP_data_integration_term = (string)subsegmental_line.Kpmp_data_integration_term.Clone();
            //    standardized_data_line.Gene_symbol = (string)subsegmental_line.Gene_symbol.Clone();
            //    standardized_data_line.Value_type = subsegmental_line.Value_type;

            //    standardized_data_line.Value = subsegmental_line.Expression_value;
            //    standardized_data_list.Add(standardized_data_line);
            //}
            //Dictionary<string, string[]> dataset_bgGenesProteins_dict = new Dictionary<string, string[]>();
            //dataset_bgGenesProteins_dict.Add(dataset, Get_deep_copy_of_bg_genes_in_upperCase());

            //KPMP_standardized_dataset_class standardized_data = new KPMP_standardized_dataset_class();
            //standardized_data.Add_to_existing_instances(standardized_data_list.ToArray(), dataset_bgGenesProteins_dict);

            //foreach (KPMP_value_type_enum value_type in value_types)
            //{
            //    Report_class.Write(" {0}", value_type);
            //}
            //Report_class.WriteLine();
            //return standardized_data;
            #endregion
        }

        public void Keep_top_x_genes_for_each_average_sample(int top_x)
        {
            this.Average_data = this.Average_data.OrderBy(l => l.Sample_name).ThenByDescending(l => l.Expression_value).ToArray();
            int data_length = Average_data.Length;
            KPMP_subsegmental_rnaSeq_line_class kpmp_line;
            List<KPMP_subsegmental_rnaSeq_line_class> keep_list = new List<KPMP_subsegmental_rnaSeq_line_class>();
            int current_kept_count = 0;
            for (int indexKPMP = 0; indexKPMP < data_length; indexKPMP++)
            {
                kpmp_line = this.Average_data[indexKPMP];
                if ((indexKPMP == 0) || (!kpmp_line.Sample_name.Equals(Average_data[indexKPMP - 1].Sample_name)))
                {
                    current_kept_count = 0;
                }
                if (current_kept_count < top_x)
                {
                    keep_list.Add(kpmp_line);
                    current_kept_count++;
                }
            }
            this.Average_data = keep_list.ToArray();
        }

        public void Keep_only_indicated_samples_for_average_data(params string[] samples)
        {
            samples = samples.Distinct().OrderBy(l => l).ToArray();
            string sample;
            int samples_length = samples.Length;
            int indexSample = 0;

            int average_length = Average_data.Length;
            int stringCompare = -2;
            KPMP_subsegmental_rnaSeq_line_class average_data_line;
            List<KPMP_subsegmental_rnaSeq_line_class> keep = new List<KPMP_subsegmental_rnaSeq_line_class>();
            this.Average_data = this.Average_data.OrderBy(l => l.Sample_name).ToArray();
            for (int indexA = 0; indexA < average_length; indexA++)
            {
                average_data_line = this.Average_data[indexA];
                stringCompare = -2;
                while ((indexSample < samples_length) && (stringCompare < 0))
                {
                    sample = samples[indexSample];
                    stringCompare = sample.CompareTo(average_data_line.Sample_name);
                    if (stringCompare < 0)
                    {
                        indexSample++;
                    }
                    else if (stringCompare == 0)
                    {
                        keep.Add(average_data_line);
                    }
                }
            }
            this.Average_data = keep.ToArray();
        }

        public void Keep_only_lines_with_indicated_patientIDs(params string[] patientIDs)
        {
            patientIDs = patientIDs.Distinct().OrderBy(l => l).ToArray();
            int indexPatientID = 0;
            int patientIDs_length = patientIDs.Length;
            string patientID;
            int stringCompare = -2;

            int data_length = this.Data.Length;
            List<KPMP_subsegmental_rnaSeq_line_class> keep = new List<KPMP_subsegmental_rnaSeq_line_class>();
            KPMP_subsegmental_rnaSeq_line_class data_line;
            this.Data = this.Data.OrderBy(l => l.PatientId).ToArray();
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                data_line = this.Data[indexData];
                stringCompare = -2;
                while ((indexPatientID < patientIDs_length) && (stringCompare < 0))
                {
                    patientID = patientIDs[indexPatientID];
                    stringCompare = patientID.CompareTo(data_line.PatientId);
                    if (stringCompare < 0)
                    {
                        indexPatientID++;
                    }
                    else if (stringCompare == 0)
                    {
                        keep.Add(data_line);
                    }
                }
            }
            this.Data = keep.ToArray();
        }

        public void Keep_only_indicated_valueTypes(params KPMP_value_type_enum[] keep_value_types)
        {
            int data_length = Data.Length;
            KPMP_subsegmental_rnaSeq_line_class rnaSet_line;
            List<KPMP_subsegmental_rnaSeq_line_class> keep = new List<KPMP_subsegmental_rnaSeq_line_class>();
            for (int indexP = 0; indexP < data_length; indexP++)
            {
                rnaSet_line = this.Data[indexP];
                if (keep_value_types.Contains(rnaSet_line.Value_type))
                {
                    keep.Add(rnaSet_line);
                }
            }
            this.Data = keep.ToArray();
        }

        public void Keep_only_lines_with_minimum_selectedValueType_and_all_other_values_of_the_same_lines(float minimum_value, KPMP_value_type_enum selected_value_type)
        {
            this.Data = this.Data.OrderBy(l => l.Gene_symbol).ThenBy(l => l.PatientId).ThenBy(l => l.Subsegment).ToArray();
            int data_length = this.Data.Length;
            KPMP_subsegmental_rnaSeq_line_class rnaSet_line;
            List<KPMP_subsegmental_rnaSeq_line_class> keep = new List<KPMP_subsegmental_rnaSeq_line_class>();
            List<KPMP_subsegmental_rnaSeq_line_class> remove = new List<KPMP_subsegmental_rnaSeq_line_class>();
            bool keep_line_block = false;
            int firstIndexData = -1;
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                rnaSet_line = this.Data[indexData];
                if ((indexData == 0)
                    || (!rnaSet_line.Gene_symbol.Equals(this.Data[indexData - 1].Gene_symbol))
                    || (!rnaSet_line.PatientId.Equals(this.Data[indexData - 1].PatientId))
                    || (!rnaSet_line.Subsegment.Equals(this.Data[indexData - 1].Subsegment)))
                {
                    keep_line_block = false;
                    firstIndexData = indexData;
                }
                if (rnaSet_line.Value_type.Equals(selected_value_type))
                {
                    if (rnaSet_line.Expression_value >= minimum_value)
                    {
                        keep_line_block = true;
                    }
                }
                if ((indexData == data_length - 1)
                    || (!rnaSet_line.Gene_symbol.Equals(this.Data[indexData + 1].Gene_symbol))
                    || (!rnaSet_line.PatientId.Equals(this.Data[indexData + 1].PatientId))
                    || (!rnaSet_line.Subsegment.Equals(this.Data[indexData + 1].Subsegment)))
                {
                    if (keep_line_block)
                    {
                        for (int indexInner = firstIndexData; indexInner <= indexData; indexInner++)
                        {
                            keep.Add(this.Data[indexInner]);
                        }
                    }
                    else
                    {
                        for (int indexInner = firstIndexData; indexInner <= indexData; indexInner++)
                        {
                            remove.Add(this.Data[indexInner]);
                        }
                    }
                }
            }
            this.Data = keep.ToArray();
        }


        public void Keep_only_lines_with_expression_value_above_minuslog10_alpha(float alpha)
        {
            float minusLog10_alpha = -(float)Math.Log10(alpha);
            int data_length = Data.Length;
            KPMP_subsegmental_rnaSeq_line_class rnaSet_line;
            List<KPMP_subsegmental_rnaSeq_line_class> keep = new List<KPMP_subsegmental_rnaSeq_line_class>();
            for (int indexP = 0; indexP < data_length; indexP++)
            {
                rnaSet_line = this.Data[indexP];
                if (rnaSet_line.Expression_value >= minusLog10_alpha)
                {
                    keep.Add(rnaSet_line);
                }
            }
            this.Data = keep.ToArray();
        }

        public void Keep_only_lines_with_expression_value_above_input_value(float input_value)
        {
            int data_length = Data.Length;
            KPMP_subsegmental_rnaSeq_line_class rnaSet_line;
            List<KPMP_subsegmental_rnaSeq_line_class> keep = new List<KPMP_subsegmental_rnaSeq_line_class>();
            for (int indexP = 0; indexP < data_length; indexP++)
            {
                rnaSet_line = this.Data[indexP];
                if (rnaSet_line.Expression_value >= input_value)
                {
                    keep.Add(rnaSet_line);
                }
            }
            this.Data = keep.ToArray();
        }

        public KPMP_integration_paper_metadata_class Get_deep_copy_of_dataset_patient_instance()
        {
            return this.Dataset_patient.Deep_copy();
        }

        private void Generate_dataset_patient_instance()
        {
            Dictionary<string, Dictionary<string, bool>> dataset_patient_considered_dict = new Dictionary<string, Dictionary<string, bool>>();
            int data_length = this.Data.Length;
            KPMP_integration_paper_metadata_line_class new_dataset_patient_line;
            List<KPMP_integration_paper_metadata_line_class> new_dataset_patient_lines = new List<KPMP_integration_paper_metadata_line_class>();
            KPMP_subsegmental_rnaSeq_line_class singleClister_line;
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                singleClister_line = this.Data[indexData];
                if (!dataset_patient_considered_dict.ContainsKey(singleClister_line.Dataset))
                {
                    dataset_patient_considered_dict.Add(singleClister_line.Dataset, new Dictionary<string, bool>());
                }
                if (!dataset_patient_considered_dict[singleClister_line.Dataset].ContainsKey(singleClister_line.PatientId))
                {
                    dataset_patient_considered_dict[singleClister_line.Dataset].Add(singleClister_line.PatientId, true);
                    new_dataset_patient_line = new KPMP_integration_paper_metadata_line_class();
                    new_dataset_patient_line.Dataset = (string)singleClister_line.Dataset.Clone();
                    new_dataset_patient_line.Libraries = new string[] { (string)singleClister_line.PatientId.Clone() };
                    new_dataset_patient_lines.Add(new_dataset_patient_line);
                }
            }
            if (Dataset_patient.Documentations.Length!=0) { throw new Exception(); }
            Dataset_patient = new KPMP_integration_paper_metadata_class();
            Dataset_patient.Add_to_array(new_dataset_patient_lines.ToArray());
        }

        public DE_class Generate_de_instance_based_on_average_data()
        {
            Fill_de_line_class fill_de_line;
            List<Fill_de_line_class> fill_de_list = new List<Fill_de_line_class>();
            foreach (KPMP_subsegmental_rnaSeq_line_class rnaSeq_line in Average_data)
            {
                fill_de_line = new Fill_de_line_class();
                fill_de_line.Names_for_de = new string[] { (string)rnaSeq_line.Sample_name.Clone() };
                fill_de_line.Symbols_for_de = new string[] { (string)rnaSeq_line.Gene_symbol.Clone() };
                fill_de_line.Value_for_de = rnaSeq_line.Expression_value;
                fill_de_list.Add(fill_de_line);
            }
            DE_class de = new DE_class();
            de.Fill_with_data_alternatively(fill_de_list.ToArray());
            return de;
        }

        private void Read_and_fill_array_old(string subdirectory, string dataset)
        {
            string complete_directory = Global_directory_class.Experimental_data_directory + subdirectory;
            string[] complete_fileNames = Directory.GetFiles(complete_directory);
            string complete_file_name;
            int complete_fileNames_length = complete_fileNames.Length;
            KPMP_subsegmental_rnaSeq_line_class new_rnaSeq_line;
            List<KPMP_subsegmental_rnaSeq_line_class> rnaSeq_list = new List<KPMP_subsegmental_rnaSeq_line_class>();
            for (int indexC = 0; indexC < complete_fileNames_length; indexC++)
            {
                complete_file_name = complete_fileNames[indexC];

                char delimiter = Global_class.Comma;
                StreamReader reader = new StreamReader(complete_file_name);
                string inputLine;
                string[] splitStrings;
                string[] columnNames;
                string gene_symbol;
                string expression_value_string;

                inputLine = reader.ReadLine();
                columnNames = inputLine.Split(delimiter);
                int columns_length = columnNames.Length;

                while ((inputLine = reader.ReadLine()) != null)
                {
                    splitStrings = inputLine.Split(delimiter);
                    if (splitStrings.Length != columns_length) { throw new Exception(); }
                    gene_symbol = splitStrings[0];

                    for (int indexCol = 1; indexCol < columns_length; indexCol++)
                    {
                        new_rnaSeq_line = new KPMP_subsegmental_rnaSeq_line_class();
                        new_rnaSeq_line.Gene_symbol = (string)gene_symbol.Clone();
                        new_rnaSeq_line.Sample_name = (string)columnNames[indexCol].Clone();
                        new_rnaSeq_line.Dataset = dataset;
                        expression_value_string = splitStrings[indexCol];
                        if (expression_value_string.Equals("#DIV/0!"))
                        {
                            new_rnaSeq_line.Expression_value = double.NaN;
                        }
                        else
                        {
                            new_rnaSeq_line.Expression_value = double.Parse(splitStrings[indexCol]);
                        }
                        rnaSeq_list.Add(new_rnaSeq_line);
                    }
                }
            }
            this.Data = rnaSeq_list.ToArray();
        }

        private void Read_r_output_and_fill_array(string subdirectory)
        {
            string complete_directory = Global_directory_class.Experimental_data_directory + subdirectory;
            string complete_fileName = complete_directory + "Sub-segmental_RNASeq_analysis_ttest.txt";
            StreamReader reader = new StreamReader(complete_fileName);



        }

        private void Read_and_fill_array(string subdirectory, string dataset)
        {
            string complete_directory = Global_directory_class.Experimental_data_directory + subdirectory;
            string[] complete_fileNames = Directory.GetFiles(complete_directory);
            string complete_file_name;
            int complete_fileNames_length = complete_fileNames.Length;
            KPMP_subsegmental_rnaSeq_line_class new_rnaSeq_line;
            List<KPMP_subsegmental_rnaSeq_line_class> rnaSeq_list = new List<KPMP_subsegmental_rnaSeq_line_class>();
            for (int indexC = 0; indexC < complete_fileNames_length; indexC++)
            {
                complete_file_name = complete_fileNames[indexC];

                char delimiter = Global_class.Tab;
                StreamReader reader = new StreamReader(complete_file_name);
                string inputLine;
                string[] splitStrings;
                string[] columnNames;
                string gene_symbol;
                string expression_value_string;

                inputLine = reader.ReadLine();
                columnNames = inputLine.Split(delimiter);
                int columns_length = columnNames.Length;

                while ((inputLine = reader.ReadLine()) != null)
                {
                    splitStrings = inputLine.Split(delimiter);
                    if (splitStrings.Length != columns_length) { throw new Exception(); }
                    gene_symbol = splitStrings[0];

                    for (int indexCol = 1; indexCol < columns_length; indexCol++)
                    {
                        new_rnaSeq_line = new KPMP_subsegmental_rnaSeq_line_class();
                        new_rnaSeq_line.Gene_symbol = (string)gene_symbol.Clone();
                        new_rnaSeq_line.Sample_name = (string)columnNames[indexCol].Clone();
                        new_rnaSeq_line.Dataset = dataset;
                        expression_value_string = splitStrings[indexCol];
                        new_rnaSeq_line.Expression_value = double.Parse(splitStrings[indexCol]);
                        rnaSeq_list.Add(new_rnaSeq_line);
                    }
                }
            }
            this.Data = rnaSeq_list.ToArray();
        }


        public KPMP_subsegmental_rnaSeq_class Deep_copy()
        {
            int data_length = this.Data.Length;
            KPMP_subsegmental_rnaSeq_class copy = (KPMP_subsegmental_rnaSeq_class)this.MemberwiseClone();
            copy.Data = new KPMP_subsegmental_rnaSeq_line_class[data_length];
            for (int indexD = 0; indexD < data_length; indexD++)
            {
                copy.Data[indexD] = this.Data[indexD].Deep_copy();
            }
            int average_data_length = this.Average_data.Length;
            copy.Average_data = new KPMP_subsegmental_rnaSeq_line_class[average_data_length];
            for (int indexD = 0; indexD < average_data_length; indexD++)
            {
                copy.Data[indexD] = this.Average_data[indexD].Deep_copy();
            }
            int bg_genes_length = this.Bg_genes_in_upperCase.Length;
            for (int indexBg = 0; indexBg < bg_genes_length; indexBg++)
            {
                copy.Bg_genes_in_upperCase[indexBg] = (string)this.Bg_genes_in_upperCase[indexBg].Clone();
            }

            return copy;
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class KPMP_compare_datasets_line_class
    {
        public string Dataset0 { get; set; }
        public string Dataset1 { get; set; }
        public string Segment { get; set; }
        public string Patient { get; set; }
        public string Gene { get; set; }

        public KPMP_value_type_enum Value_type { get; set; }
        public double Expression_value0 { get; set; }
        public double Expression_value1 { get; set; }
        public double Abs_diff_expression_value { get; set; }

        public KPMP_compare_datasets_line_class Deep_copy()
        {
            KPMP_compare_datasets_line_class copy = (KPMP_compare_datasets_line_class)this.MemberwiseClone();
            copy.Dataset0 = (string)this.Dataset0.Clone();
            copy.Dataset1 = (string)this.Dataset1.Clone();
            copy.Segment = (string)this.Segment.Clone();
            copy.Patient = (string)this.Patient.Clone();
            copy.Gene = (string)this.Gene.Clone();
            return copy;
        }
    }

    class KPMP_compare_dataset_readWriteOptions_class : ReadWriteOptions_base
    {
        public KPMP_compare_dataset_readWriteOptions_class(string fileName)
        {
            this.File = Global_directory_class.Results_directory + fileName;
            this.Key_propertyNames = new string[] { "Dataset0", "Dataset1", "Segment", "Patient", "Gene", "Value_type","Expression_value0", "Expression_value1", "Abs_diff_expression_value" };
            this.Key_columnNames = this.Key_propertyNames;
            this.HeadlineDelimiters = new char[] { Global_class.Tab };
            this.LineDelimiters = new char[] { Global_class.Tab };
            this.File_has_headline = true;
            this.Report = ReadWrite_report_enum.Report_main;
        }
    }

    class KPMP_compare_dataset_class
    {
        public KPMP_compare_datasets_line_class[] KPMP_compare { get; set; }

        public void Write_file(string fileName)
        {
            KPMP_compare_dataset_readWriteOptions_class readWriteOptions = new KPMP_compare_dataset_readWriteOptions_class(fileName);
            ReadWriteClass.WriteData(this.KPMP_compare, readWriteOptions);
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class KPMP_subsegmental_proteomics_line_class
    {
        public string Description { get; set; }
      //  public string Full_gene_name_from_description { get; set; }
        public string Gene_symbol { get; set; }
        public string SampleName { get; set; }
        public string Accession { get; set; }
        public string PatientId { get; set; }
        public string Subsegment { get; set; }
        public string Dataset { get; set; }
        public string Kpmp_data_integration_term { get; set; }
        public double Value { get; set; }
        public KPMP_value_type_enum Value_type { get; set; }

        #region Equal set
        public static KPMP_subsegmental_proteomics_line_class[] Order_by_set(KPMP_subsegmental_proteomics_line_class[] array)
        {
            array = array.OrderBy(l => l.Dataset).ThenBy(l => l.PatientId).ThenBy(l => l.Subsegment).ThenBy(l => l.SampleName).ThenBy(l => l.Kpmp_data_integration_term).ThenBy(l => l.Value_type).ToArray();
            return array;
        }

        public static KPMP_subsegmental_proteomics_line_class[] Order_by_set_and_descending_value(KPMP_subsegmental_proteomics_line_class[] array)
        {
            array = array.OrderBy(l => l.Dataset).ThenBy(l => l.PatientId).ThenBy(l => l.Subsegment).ThenBy(l => l.SampleName).ThenBy(l => l.Kpmp_data_integration_term).ThenBy(l => l.Value_type).ThenByDescending(l => l.Value).ToArray();
            return array;
        }

        public static KPMP_subsegmental_proteomics_line_class[] Order_by_set_gene_symbol_and_descending_value(KPMP_subsegmental_proteomics_line_class[] array)
        {
            array = array.OrderBy(l => l.Dataset).ThenBy(l => l.PatientId).ThenBy(l => l.Subsegment).ThenBy(l => l.SampleName).ThenBy(l => l.Kpmp_data_integration_term).ThenBy(l => l.Value_type).ThenBy(l => l.Gene_symbol).ThenByDescending(l => l.Value).ToArray();
            return array;
        }

        public bool Equal_set(KPMP_subsegmental_proteomics_line_class other)
        {
            bool equal = (this.Dataset.Equals(other.Dataset))
                         && (this.PatientId.Equals(other.PatientId))
                         && (this.Subsegment.Equals(other.Subsegment))
                         && (this.SampleName.Equals(other.SampleName))
                         && (this.Kpmp_data_integration_term.Equals(other.Kpmp_data_integration_term))
                         && (this.Value_type.Equals(other.Value_type));
            return equal;
        }

        public bool Equal_set_and_gene_symbol(KPMP_subsegmental_proteomics_line_class other)
        {
            bool equal = (this.Dataset.Equals(other.Dataset))
                         && (this.PatientId.Equals(other.PatientId))
                         && (this.Subsegment.Equals(other.Subsegment))
                         && (this.SampleName.Equals(other.SampleName))
                         && (this.Kpmp_data_integration_term.Equals(other.Kpmp_data_integration_term))
                         && (this.Value_type.Equals(other.Value_type))
                         && (this.Gene_symbol.Equals(other.Gene_symbol));
            return equal;
        }
        #endregion

        public KPMP_subsegmental_proteomics_line_class()
        {
            Description = "";
            Gene_symbol = "";
            SampleName = "";
            PatientId = "";
            Subsegment = "";
            Dataset = "";
            Kpmp_data_integration_term = "";
        }

        public KPMP_subsegmental_proteomics_line_class Deep_copy()
        {
            KPMP_subsegmental_proteomics_line_class copy = (KPMP_subsegmental_proteomics_line_class)this.MemberwiseClone();
            copy.Description = (string)this.Description.Clone();
            copy.Gene_symbol = (string)this.Gene_symbol.Clone();
            copy.Dataset = (string)this.Dataset.Clone();
            copy.SampleName = (string)this.SampleName.Clone();
            copy.Subsegment = (string)this.Subsegment.Clone();
            copy.PatientId = (string)this.PatientId.Clone();
            copy.Kpmp_data_integration_term = (string)this.Kpmp_data_integration_term.Clone();
            return copy;
        }
    }

    class KPMP_subsegmental_proteomics_class
    {
        public KPMP_subsegmental_proteomics_line_class[] Data { get; set; }
        public KPMP_integration_paper_metadata_class Dataset_patient { get; set; }
        public string[] Bg_proteins_in_upper_case { get; set; }

        public KPMP_subsegmental_proteomics_class()
        {
            this.Data = new KPMP_subsegmental_proteomics_line_class[0];
            this.Bg_proteins_in_upper_case = new string[0];
            this.Dataset_patient = new KPMP_integration_paper_metadata_class();
        }

        private void Add_to_array(KPMP_subsegmental_proteomics_line_class[] add_data)
        {
            int this_length = this.Data.Length;
            int add_length = add_data.Length;
            int new_length = this_length + add_length;
            KPMP_subsegmental_proteomics_line_class[] new_data = new KPMP_subsegmental_proteomics_line_class[new_length];
            int indexNew = -1;
            for (int indexThis = 0; indexThis < this_length; indexThis++)
            {
                indexNew++;
                new_data[indexNew] = this.Data[indexThis];
            }
            for (int indexAdd = 0; indexAdd < add_length; indexAdd++)
            {
                indexNew++;
                new_data[indexNew] = add_data[indexAdd];
            }
            this.Data = new_data;
        }

        private void Generate_bg_proteins_in_upper_case()
        {
            Data = Data.OrderBy(l => l.Gene_symbol).ToArray();
            int data_length = Data.Length;
            KPMP_subsegmental_proteomics_line_class proteomics_line;
            List<string> bg_proteins = new List<string>();
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                proteomics_line = Data[indexData];
                if ((indexData == 0)
                    || (!proteomics_line.Gene_symbol.Equals(Data[indexData - 1].Gene_symbol)))
                {
                    bg_proteins.Add(proteomics_line.Gene_symbol.ToUpper());
                }
            }
            this.Bg_proteins_in_upper_case = bg_proteins.ToArray();
        }

        private void Add_subsegment_valueType_patientID_harmonize_statisticsSampleNames_replace_pvalue_by_minusLog10pvalue_and_replace_fdr_by_minusLog10FDR()
        {
            int data_length = Data.Length;
            KPMP_subsegmental_proteomics_line_class proteomics_line;
            string[] splitStrings;
            int splitStrings_length;
            StringBuilder sb = new StringBuilder();
            //string gene_symbol_splitString;
            //string originally_assigned_gene_symbol_string;
            //string[] gene_symbol_splitString_splitStrings;
            //int indexGS = 0;
            //int indexFirstEqualSign = 0;
            //int indexEndOfFullGeneName = 0;
            List<double> sp_pvalues = new List<double>();
            //bool move_back = true;
            //bool letter_label_passed = false;
            for (int indexD = 0; indexD < data_length; indexD++)
            {
                proteomics_line = this.Data[indexD];
                 if (((proteomics_line.SampleName.IndexOf("-") != -1) && (proteomics_line.SampleName.IndexOf("P") != 0))
                    || (proteomics_line.SampleName.IndexOf("OSU") != -1)
                    || (proteomics_line.SampleName.IndexOf("IU") != -1))
                {
                    splitStrings = proteomics_line.SampleName.Split('-');
                    splitStrings_length = splitStrings.Length;
                    if (splitStrings_length != 2) { throw new Exception(); }
                    proteomics_line.PatientId = (string)splitStrings[1].Clone();// splitStrings[0].Replace("M", "");
                    proteomics_line.Subsegment = (string)splitStrings[0].Clone();
                    proteomics_line.Value_type = KPMP_value_type_enum.Single_value;
                }
                else if (proteomics_line.SampleName.ToUpper().IndexOf("AVG") == 0)
                {
                    proteomics_line.PatientId = KPMP_data_integration_class.Get_combined_patients_label(); ;
                    proteomics_line.Subsegment = (string)proteomics_line.SampleName.Split(' ')[1];
                    proteomics_line.Value_type = KPMP_value_type_enum.Average;
                    proteomics_line.SampleName = "Statistics";
                }
                else if (proteomics_line.SampleName.ToUpper().IndexOf("MEDIAN") == 0)
                {
                    proteomics_line.PatientId = KPMP_data_integration_class.Get_combined_patients_label(); ;
                    proteomics_line.Subsegment = (string)proteomics_line.SampleName.Split(' ')[1];
                    proteomics_line.Value_type = KPMP_value_type_enum.Median;
                    proteomics_line.SampleName = "Statistics";
                }
                else if (proteomics_line.SampleName.ToUpper().IndexOf("FOLD CHANGE G/T") == 0)
                {
                    proteomics_line.PatientId = KPMP_data_integration_class.Get_combined_patients_label(); ;
                    proteomics_line.Subsegment = "Glom";
                    proteomics_line.Value_type = KPMP_value_type_enum.Ratioavg;
                    proteomics_line.SampleName = "Statistics";
                }
                else if (proteomics_line.SampleName.ToUpper().IndexOf("P_VALUE") == 0)
                {
                    proteomics_line.PatientId = KPMP_data_integration_class.Get_combined_patients_label(); ;
                    proteomics_line.Subsegment = (string)proteomics_line.SampleName.Clone();
                    proteomics_line.Value = -(float)Math.Log10(proteomics_line.Value);
                    proteomics_line.Value_type = KPMP_value_type_enum.Minus_log10_pvalue;
                    proteomics_line.SampleName = "Statistics";
                }
                else if (proteomics_line.SampleName.ToUpper().IndexOf("FDR")==0)
                {
                    proteomics_line.PatientId = KPMP_data_integration_class.Get_combined_patients_label(); ;
                    proteomics_line.Subsegment = (string)proteomics_line.SampleName.Clone();
                    proteomics_line.Value = -(float)Math.Log10(proteomics_line.Value);
                    proteomics_line.Value_type = KPMP_value_type_enum.Minus_log10_fdr;
                    proteomics_line.SampleName = "Statistics";
                }
                else
                {
                    splitStrings = proteomics_line.SampleName.Split('-');
                    proteomics_line.PatientId = splitStrings[0].Substring(1, splitStrings[0].Length - 1);
                    proteomics_line.Subsegment = (string)splitStrings[1].Clone();
                    throw new Exception();
                }
            }
        }

        private void Remove_all_singleValues()
        {
            int data_length = Data.Length;
            KPMP_subsegmental_proteomics_line_class proteomics_line;
            List<KPMP_subsegmental_proteomics_line_class> keep = new List<KPMP_subsegmental_proteomics_line_class>();
            for (int indexD = 0; indexD < data_length; indexD++)
            {
                proteomics_line = this.Data[indexD];
                if (!proteomics_line.Value_type.Equals(KPMP_value_type_enum.Single_value))
                {
                    keep.Add(proteomics_line);
                }
            }
            this.Data = keep.ToArray();
        }

        public void Keep_top_x_lines_per_sampleName_based_on_descending_indicated_valueType_and_keep_all_other_valueTypes_of_kept_genes(int keep_top_x, KPMP_value_type_enum selective_value_type)
        {
            this.Data = this.Data.OrderBy(l => l.Dataset).ThenBy(l => l.PatientId).ThenBy(l => l.Subsegment).ThenBy(l => l.Kpmp_data_integration_term).ThenBy(l => l.Value_type).ThenByDescending(l => l.Value).ToArray();
            int data_length = this.Data.Length;
            KPMP_subsegmental_proteomics_line_class proteomics_line;
            KPMP_subsegmental_proteomics_line_class inner_proteomics_line;
            List<KPMP_subsegmental_proteomics_line_class> keep = new List<KPMP_subsegmental_proteomics_line_class>();
            int kept_lines_count = 0;
            int firstIndex_sameSet = -1;
            List<string> sampleNames = new List<string>();
            List<string> subsegments = new List<string>();
            Dictionary<string, bool> current_keep_proteins_dict = new Dictionary<string, bool>();
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                proteomics_line = this.Data[indexData];
                if (  (indexData == 0)
                    ||(!proteomics_line.Dataset.Equals(Data[indexData - 1].Dataset))
                    ||(!proteomics_line.PatientId.Equals(Data[indexData - 1].PatientId))
                    ||(!proteomics_line.Subsegment.Equals(Data[indexData - 1].Subsegment))
                    ||(!proteomics_line.Kpmp_data_integration_term.Equals(Data[indexData - 1].Kpmp_data_integration_term)))
                {
                    kept_lines_count = 0;
                    firstIndex_sameSet = indexData;
                    current_keep_proteins_dict.Clear();
                }
                if ((proteomics_line.Value_type.Equals(selective_value_type))&&(kept_lines_count < keep_top_x))
                {
                    current_keep_proteins_dict.Add(proteomics_line.Gene_symbol, true);
                    kept_lines_count++;
                }
                if ((indexData == data_length-1)
                    || (!proteomics_line.Dataset.Equals(Data[indexData + 1].Dataset))
                    || (!proteomics_line.PatientId.Equals(Data[indexData + 1].PatientId))
                    || (!proteomics_line.Subsegment.Equals(Data[indexData + 1].Subsegment))
                    || (!proteomics_line.Kpmp_data_integration_term.Equals(Data[indexData + 1].Kpmp_data_integration_term)))
                {
                    for (int indexInner = firstIndex_sameSet; indexInner <= indexData; indexInner++)
                    {
                        inner_proteomics_line = this.Data[indexInner];
                        if (current_keep_proteins_dict.ContainsKey(inner_proteomics_line.Gene_symbol))
                        {
                            keep.Add(inner_proteomics_line);
                        }
                    }
                }
            }
            this.Data = keep.ToArray();
        }

        public void Keep_only_lines_with_minimum_fold_change(float minimum_fold_change)
        {
            this.Data = this.Data.OrderBy(l => l.SampleName).ThenByDescending(l => l.Value).ToArray();
            int data_length = this.Data.Length;
            KPMP_subsegmental_proteomics_line_class proteomics_line;
            List<KPMP_subsegmental_proteomics_line_class> keep = new List<KPMP_subsegmental_proteomics_line_class>();
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                proteomics_line = this.Data[indexData];
                if (proteomics_line.Value >= minimum_fold_change)
                {
                    keep.Add(proteomics_line);
                }
            }
            this.Data = keep.ToArray();
        }

        public void Keep_only_lines_with_minimum_selectedValueType_and_all_other_values_of_the_same_lines(float minimum_value, KPMP_value_type_enum selected_value_type)
        {
            this.Data = this.Data.OrderBy(l => l.Accession).ThenBy(l => l.Description).ThenBy(l => l.PatientId).ThenBy(l => l.Subsegment).ToArray();
            int data_length = this.Data.Length;
            KPMP_subsegmental_proteomics_line_class proteomics_line;
            List<KPMP_subsegmental_proteomics_line_class> keep = new List<KPMP_subsegmental_proteomics_line_class>();
            List<KPMP_subsegmental_proteomics_line_class> remove = new List<KPMP_subsegmental_proteomics_line_class>();
            bool keep_line_block = false;
            int firstIndexData = -1;
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                proteomics_line = this.Data[indexData];
                if ((indexData == 0)
                    || (!proteomics_line.Accession.Equals(this.Data[indexData - 1].Accession))
                    || (!proteomics_line.Description.Equals(this.Data[indexData - 1].Description))
                    || (!proteomics_line.PatientId.Equals(this.Data[indexData - 1].PatientId))
                    || (!proteomics_line.Subsegment.Equals(this.Data[indexData - 1].Subsegment)))
                {
                    keep_line_block = false;
                    firstIndexData = indexData;
                }
                if (proteomics_line.Value_type.Equals(selected_value_type))
                {
                    if (proteomics_line.Value >= minimum_value)
                    {
                        keep_line_block = true;
                    }
                }
                if ((indexData == data_length - 1)
                    || (!proteomics_line.Accession.Equals(this.Data[indexData + 1].Accession))
                    || (!proteomics_line.Description.Equals(this.Data[indexData + 1].Description))
                    || (!proteomics_line.PatientId.Equals(this.Data[indexData + 1].PatientId))
                    || (!proteomics_line.Subsegment.Equals(this.Data[indexData + 1].Subsegment)))
                {
                    if (keep_line_block)
                    {
                        for (int indexInner = firstIndexData; indexInner <= indexData; indexInner++)
                        {
                            keep.Add(this.Data[indexInner]);
                        }
                    }
                    else
                    {
                        for (int indexInner = firstIndexData; indexInner <= indexData; indexInner++)
                        {
                            remove.Add(this.Data[indexInner]);
                        }
                    }
                }
            }
            this.Data = keep.ToArray();
        }

        public void Keep_only_lines_with_minimum_log2foldChange_and_all_other_values_of_the_same_lines(float minimum_log2foldchange)
        {
            this.Data = this.Data.OrderBy(l => l.Accession).ThenBy(l => l.Description).ThenBy(l => l.PatientId).ThenBy(l => l.Subsegment).ToArray();
            int data_length = this.Data.Length;
            KPMP_subsegmental_proteomics_line_class proteomics_line;
            List<KPMP_subsegmental_proteomics_line_class> keep = new List<KPMP_subsegmental_proteomics_line_class>();
            List<KPMP_subsegmental_proteomics_line_class> remove = new List<KPMP_subsegmental_proteomics_line_class>();
            bool keep_line_block = false;
            int firstIndexData = -1;
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                proteomics_line = this.Data[indexData];
                if ((indexData == 0)
                    || (!proteomics_line.Accession.Equals(this.Data[indexData - 1].Accession))
                    || (!proteomics_line.Description.Equals(this.Data[indexData - 1].Description))
                    || (!proteomics_line.PatientId.Equals(this.Data[indexData - 1].PatientId))
                    || (!proteomics_line.Subsegment.Equals(this.Data[indexData - 1].Subsegment)))
                {
                    keep_line_block = false;
                    firstIndexData = indexData;
                }
                if ((proteomics_line.Value_type.Equals(KPMP_value_type_enum.Log2_ratioavg))
                    || (proteomics_line.Value_type.Equals(KPMP_value_type_enum.Log2_ratio_singlepatient))
                    || (proteomics_line.Value_type.Equals(KPMP_value_type_enum.Log2_single_value)))
                {
                    if (proteomics_line.Value >= minimum_log2foldchange)
                    {
                        keep_line_block = true;
                    }
                }
                if ((indexData == data_length - 1)
                    || (!proteomics_line.Accession.Equals(this.Data[indexData + 1].Accession))
                    || (!proteomics_line.Description.Equals(this.Data[indexData + 1].Description))
                    || (!proteomics_line.PatientId.Equals(this.Data[indexData + 1].PatientId))
                    || (!proteomics_line.Subsegment.Equals(this.Data[indexData + 1].Subsegment)))
                {
                    if (keep_line_block)
                    {
                        for (int indexInner = firstIndexData; indexInner <= indexData; indexInner++)
                        {
                            keep.Add(this.Data[indexInner]);
                        }
                    }
                    else
                    {
                        for (int indexInner = firstIndexData; indexInner <= indexData; indexInner++)
                        {
                            remove.Add(this.Data[indexInner]);
                        }
                    }
                }
            }
            if (remove.Count > 0) { throw new Exception(); }
            this.Data = keep.ToArray();
        }

        private void Set_kpmp_integration_term()
        {
            int data_length = this.Data.Length;
            KPMP_subsegmental_proteomics_line_class proteomics_line;
            for (int indexD = 0; indexD < data_length; indexD++)
            {
                proteomics_line = this.Data[indexD];
                proteomics_line.Kpmp_data_integration_term = KPMP_data_integration_class.Get_kpmp_integration_term_plus_patientID(proteomics_line.Subsegment, proteomics_line.PatientId);
            }
        }

        private void Remove_duplicated_lines_based_on_dataset_geneSymbol_sampleName_by_keeping_higher_selected_value_type(KPMP_value_type_enum selected_value_type)
        {
            //Report_class.Write_code_imperfect_line("{0}: Still Remove_duplicated_lines_based_on_geneSymbol_by_keeping_higher_value", typeof(KPMP_subsegmental_proteomics_class).Name);
            this.Data = this.Data.OrderBy(l => l.Dataset).ThenBy(l => l.Gene_symbol).ThenBy(l=>l.SampleName).ThenBy(l => l.Value_type).ThenByDescending(l => l.Value).ToArray();
            KPMP_subsegmental_proteomics_line_class proteomics_line;
            Dictionary<string, Dictionary<string, Dictionary<string,Dictionary<string, bool>>>> keep_dataset_sampleName_description_accession_dict = new Dictionary<string, Dictionary<string, Dictionary<string, Dictionary<string, bool>>>>();
            int data_length = this.Data.Length;
            bool kept_line_for_gene_identified = false;
            for (int indexD = 0; indexD < data_length; indexD++)
            {
                proteomics_line = this.Data[indexD];
                if ((indexD == 0)
                    || (!proteomics_line.Dataset.Equals(this.Data[indexD - 1].Dataset))
                    || (!proteomics_line.Gene_symbol.Equals(this.Data[indexD - 1].Gene_symbol))
                    || (!proteomics_line.SampleName.Equals(this.Data[indexD-1].SampleName)))
                {
                    kept_line_for_gene_identified = false;
                }
                if ((proteomics_line.Value_type.Equals(selected_value_type))
                    && (!kept_line_for_gene_identified))
                {
                    if (!keep_dataset_sampleName_description_accession_dict.ContainsKey(proteomics_line.Dataset))
                    {
                        keep_dataset_sampleName_description_accession_dict.Add(proteomics_line.Dataset, new Dictionary<string, Dictionary<string, Dictionary<string, bool>>>());
                    }
                    if (!keep_dataset_sampleName_description_accession_dict[proteomics_line.Dataset].ContainsKey(proteomics_line.SampleName))
                    {
                        keep_dataset_sampleName_description_accession_dict[proteomics_line.Dataset].Add(proteomics_line.SampleName, new Dictionary<string, Dictionary<string, bool>>());
                    }
                    if (!keep_dataset_sampleName_description_accession_dict[proteomics_line.Dataset][proteomics_line.SampleName].ContainsKey(proteomics_line.Description))
                    {
                        keep_dataset_sampleName_description_accession_dict[proteomics_line.Dataset][proteomics_line.SampleName].Add(proteomics_line.Description, new Dictionary<string, bool>());
                    }
                    keep_dataset_sampleName_description_accession_dict[proteomics_line.Dataset][proteomics_line.SampleName][proteomics_line.Description].Add(proteomics_line.Accession, true);
                    kept_line_for_gene_identified = true;
                }
                if ((indexD == data_length-1)
                    || (!proteomics_line.Dataset.Equals(this.Data[indexD + 1].Dataset))
                    || (!proteomics_line.Gene_symbol.Equals(this.Data[indexD + 1].Gene_symbol))
                    || (!proteomics_line.SampleName.Equals(this.Data[indexD + 1].SampleName)))
                {
                    if (!kept_line_for_gene_identified) { throw new Exception(); }
                }
            }

            List<KPMP_subsegmental_proteomics_line_class> keep = new List<KPMP_subsegmental_proteomics_line_class>();
            List<KPMP_subsegmental_proteomics_line_class> remove = new List<KPMP_subsegmental_proteomics_line_class>();
            for (int indexD = 0; indexD < data_length; indexD++)
            {
                proteomics_line = this.Data[indexD];
                if ((keep_dataset_sampleName_description_accession_dict.ContainsKey(proteomics_line.Dataset))
                    && (keep_dataset_sampleName_description_accession_dict[proteomics_line.Dataset].ContainsKey(proteomics_line.SampleName))
                    && (keep_dataset_sampleName_description_accession_dict[proteomics_line.Dataset][proteomics_line.SampleName].ContainsKey(proteomics_line.Description))
                    && (keep_dataset_sampleName_description_accession_dict[proteomics_line.Dataset][proteomics_line.SampleName][proteomics_line.Description].ContainsKey(proteomics_line.Accession))
                    && (keep_dataset_sampleName_description_accession_dict[proteomics_line.Dataset][proteomics_line.SampleName][proteomics_line.Description][proteomics_line.Accession].Equals(true)))
                {
                    keep.Add(proteomics_line);
                }
                else
                {
                    remove.Add(proteomics_line);
                }
            }
            this.Data = keep.ToArray();
        }

        private void Generate_ratios_and_log2ratios_for_single_patients_for_each_geneSymbol()
        {
            int data_length = this.Data.Length;
            KPMP_subsegmental_proteomics_line_class proteomics_line;
            KPMP_subsegmental_proteomics_line_class new_proteomics_line;
            List<KPMP_subsegmental_proteomics_line_class> new_lines = new List<KPMP_subsegmental_proteomics_line_class>();
            this.Data = this.Data.OrderBy(l => l.PatientId).ThenBy(l => l.Gene_symbol).ToArray();
            int firstIndex_sameGenePatient = -1;
            double glom_value = -1;
            double ti_value = -1;
            for (int indexD = 0; indexD < data_length; indexD++)
            {
                proteomics_line = Data[indexD];
                if ((proteomics_line.SampleName.ToUpper().IndexOf("AVG") == -1)
                    && (proteomics_line.SampleName.ToUpper().IndexOf("RATIO") == -1)
                    && (proteomics_line.SampleName.ToUpper().IndexOf("P-VALUE") == -1)
                    && (proteomics_line.SampleName.ToUpper().IndexOf("MEDIAN") == -1)
                    && (proteomics_line.SampleName.ToUpper().IndexOf("FOLD CHANGE") == -1))
                {
                    if ((indexD == 0)
                            || (!proteomics_line.PatientId.Equals(Data[indexD - 1].PatientId))
                            || (!proteomics_line.Gene_symbol.Equals(Data[indexD - 1].Gene_symbol)))
                    {
                        firstIndex_sameGenePatient = indexD;
                        glom_value = -1;
                        ti_value = -1;
                    }
                    if (proteomics_line.Subsegment.Equals("Glom"))
                    {
                        if (glom_value != -1) { throw new Exception(); }
                        glom_value = proteomics_line.Value;
                    }
                    else if (proteomics_line.Subsegment.Equals("TI"))
                    {
                        if (ti_value != -1) { throw new Exception(); }
                        ti_value = proteomics_line.Value;
                    }
                    else { throw new Exception(); }
                    if ((indexD == data_length - 1)
                        || (!proteomics_line.PatientId.Equals(Data[indexD + 1].PatientId))
                        || (!proteomics_line.Gene_symbol.Equals(Data[indexD + 1].Gene_symbol)))
                    {
                        if ((ti_value == -1) || (glom_value == -1)) { throw new Exception(); }
                        ti_value = ti_value + 1;
                        glom_value = glom_value + 1;
                        new_proteomics_line = proteomics_line.Deep_copy();
                        new_proteomics_line.SampleName = "Ratio Glom vs ProxTub-" + new_proteomics_line.PatientId;
                        new_proteomics_line.Subsegment = "Ratio Glom vs ProxTub";
                        new_proteomics_line.Value = glom_value / ti_value;
                        new_proteomics_line.Value_type = KPMP_value_type_enum.Ratio_singlepatient;
                        new_lines.Add(new_proteomics_line);
                        new_proteomics_line = proteomics_line.Deep_copy();
                        new_proteomics_line.SampleName = "Ratio ProxTub vs Glom-" + new_proteomics_line.PatientId;
                        new_proteomics_line.Subsegment = "Ratio ProxTub vs Glom";
                        new_proteomics_line.Value = ti_value / glom_value;
                        new_proteomics_line.Value_type = KPMP_value_type_enum.Ratio_singlepatient;
                        new_lines.Add(new_proteomics_line);

                        new_proteomics_line = proteomics_line.Deep_copy();
                        new_proteomics_line.SampleName = "Log2Ratio Glom vs ProxTub-" + new_proteomics_line.PatientId;
                        new_proteomics_line.Subsegment = "Log2Ratio Glom vs ProxTub";
                        new_proteomics_line.Value = (float)Math.Log(glom_value / ti_value, 2);
                        new_proteomics_line.Value_type = KPMP_value_type_enum.Log2_ratio_singlepatient;
                        new_lines.Add(new_proteomics_line);
                        new_proteomics_line = proteomics_line.Deep_copy();
                        new_proteomics_line.SampleName = "Log2Ratio ProxTub vs Glom-" + new_proteomics_line.PatientId;
                        new_proteomics_line.Subsegment = "Log2Ratio ProxTub vs Glom";
                        new_proteomics_line.Value = (float)Math.Log(ti_value / glom_value, 2);
                        new_proteomics_line.Value_type = KPMP_value_type_enum.Log2_ratio_singlepatient;
                        new_lines.Add(new_proteomics_line);
                    }
                }
            }
            new_lines.AddRange(this.Data);
            this.Data = new_lines.ToArray();
        }

        private void Generate_new_lines_by_inverting_ratio_for_average_ratios()
        {
            int data_length = this.Data.Length;
            KPMP_subsegmental_proteomics_line_class proteomics_line;
            KPMP_subsegmental_proteomics_line_class new_proteomics_line;
            List<KPMP_subsegmental_proteomics_line_class> new_lines = new List<KPMP_subsegmental_proteomics_line_class>();
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                proteomics_line = Data[indexData];
                if (proteomics_line.SampleName.ToUpper().Equals("Ratio GvsTI".ToUpper()))
                {
                    new_proteomics_line = proteomics_line.Deep_copy();
                    new_proteomics_line.Value = 1F / proteomics_line.Value;
                    new_proteomics_line.SampleName = "Ratio TIvsG".ToUpper();
                    new_proteomics_line.Subsegment = "TI".ToUpper();
                    new_lines.Add(new_proteomics_line);
                }
            }
            new_lines.AddRange(this.Data);
            this.Data = new_lines.ToArray();
        }

        private void Generate_new_lines_by_inverting_ratios()
        {
            int data_length = this.Data.Length;
            KPMP_subsegmental_proteomics_line_class proteomics_line;
            KPMP_subsegmental_proteomics_line_class new_proteomics_line;
            List<KPMP_subsegmental_proteomics_line_class> new_lines = new List<KPMP_subsegmental_proteomics_line_class>();
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                proteomics_line = Data[indexData];
                if ((proteomics_line.Value_type.Equals(KPMP_value_type_enum.Ratioavg))
                    || (proteomics_line.Value_type.Equals(KPMP_value_type_enum.Ratio_singlepatient)))
                {
                    new_proteomics_line = proteomics_line.Deep_copy();
                    new_proteomics_line.Value = 1F / proteomics_line.Value;
                    new_proteomics_line.SampleName = "Ratio TIvsG".ToUpper();
                    new_proteomics_line.Subsegment = "TI".ToUpper();
                    new_lines.Add(new_proteomics_line);
                }
            }
            new_lines.AddRange(this.Data);
            this.Data = new_lines.ToArray();
        }

        private void Duplicate_lines_and_add_cell_specificity()
        {
            int data_length = Data.Length;
            KPMP_subsegmental_proteomics_line_class data_line;
            KPMP_subsegmental_proteomics_line_class new_data_line;
            List<KPMP_subsegmental_proteomics_line_class> add = new List<KPMP_subsegmental_proteomics_line_class>();
            for (int indexD = 0; indexD < data_length; indexD++)
            {
                data_line = this.Data[indexD];
                if (data_line.Subsegment.ToUpper().Equals("GLOM"))
                {
                    new_data_line = data_line.Deep_copy();
                    new_data_line.Subsegment = new_data_line.Subsegment + " - POD";
                    add.Add(new_data_line);
                    new_data_line = data_line.Deep_copy();
                    new_data_line.Subsegment = new_data_line.Subsegment + " - Mesangial";
                    add.Add(new_data_line);
                    new_data_line = data_line.Deep_copy();
                    new_data_line.Subsegment = new_data_line.Subsegment + " - EPC";
                    add.Add(new_data_line);
                    data_line.Subsegment = data_line.Subsegment + " - EC";
                }
                else if (data_line.Subsegment.ToUpper().Equals("TI"))
                {
                    //new_data_line = data_line.Deep_copy();
                    //new_data_line.Subsegment = new_data_line.Subsegment + " - TAL";
                    //add.Add(new_data_line);
                    //new_data_line = data_line.Deep_copy();
                    //new_data_line.Subsegment = new_data_line.Subsegment + " - PC";
                    //add.Add(new_data_line);
                    //new_data_line = data_line.Deep_copy();
                    //new_data_line.Subsegment = new_data_line.Subsegment + " - IC";
                    //add.Add(new_data_line);
                    //new_data_line = data_line.Deep_copy();
                    //new_data_line.Subsegment = new_data_line.Subsegment + " - MAC";
                    //add.Add(new_data_line);
                    data_line.Subsegment = data_line.Subsegment + " - PT";
                }
                else
                {
                    throw new Exception();
                }
            }
            Add_to_array(add.ToArray());
        }

        private void Check_if_patientID_is_there()
        {
            bool is_there = false;
            string patientID = "18-342";
            foreach (KPMP_subsegmental_proteomics_line_class proteomics_line in this.Data)
            {
                if (proteomics_line.PatientId.Equals(patientID))
                {
                    is_there = true;
                }
            }
            if (!is_there)
            {
                throw new Exception();
            }
        }

        private void Specify_subsegment_based_on_fold_change()
        {
            int data_length = Data.Length;
            KPMP_subsegmental_proteomics_line_class data_line;
            KPMP_subsegmental_proteomics_line_class inner_data_line;
            this.Data = this.Data.OrderBy(l => l.Dataset).ThenBy(l => l.Accession).ThenBy(l => l.Description).ThenBy(l => l.Subsegment).ToArray();
            int firstIndex_same_description = -1;
            double fold_change_glom_over_ti = -1;
            double reverse_fold_change_only_to_check = -1;
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                data_line = this.Data[indexData];
                if ((indexData == 0)
                    || (!data_line.Dataset.Equals(this.Data[indexData - 1].Dataset))
                    || (!data_line.Accession.Equals(this.Data[indexData - 1].Accession))
                    || (!data_line.Description.Equals(this.Data[indexData - 1].Description)))
                {
                    firstIndex_same_description = indexData;
                    fold_change_glom_over_ti = -1;
                    reverse_fold_change_only_to_check = -1;
                }
                if ((data_line.Value_type.Equals(KPMP_value_type_enum.Ratioavg))
                    && (data_line.Subsegment.ToUpper().Equals("GLOM")))
                {
                    if (fold_change_glom_over_ti != -1) { throw new Exception(); }
                    fold_change_glom_over_ti = data_line.Value;
                }
                if (data_line.Subsegment.ToUpper().Equals("TI"))
                {
                    if (reverse_fold_change_only_to_check != -1) { throw new Exception(); }
                    reverse_fold_change_only_to_check = data_line.Value;  // this does not exist here
                }
                if ((indexData == data_length - 1)
                    || (!data_line.Dataset.Equals(this.Data[indexData + 1].Dataset))
                    || (!data_line.Accession.Equals(this.Data[indexData + 1].Accession))
                    || (!data_line.Description.Equals(this.Data[indexData + 1].Description)))
                {
                    if (fold_change_glom_over_ti == -1) { throw new Exception(); }
                    if (reverse_fold_change_only_to_check != -1) { throw new Exception(); } 
                    for (int indexInner = firstIndex_same_description; indexInner <= indexData; indexInner++)
                    {
                        inner_data_line = this.Data[indexInner];
                        switch (inner_data_line.Value_type)
                        {
                            case KPMP_value_type_enum.Minus_log10_pvalue:
                            case KPMP_value_type_enum.Minus_log10_fdr:
                            case KPMP_value_type_enum.Log2_ratioavg:
                            case KPMP_value_type_enum.Ratioavg:
                                if (fold_change_glom_over_ti > 1)
                                {
                                    inner_data_line.Subsegment = "Glom";
                                }
                                else
                                {
                                    inner_data_line.Subsegment = "TI";
                                }
                                break;
                            case KPMP_value_type_enum.Average:
                            case KPMP_value_type_enum.Ratio_singlepatient:
                            case KPMP_value_type_enum.Single_value:
                            case KPMP_value_type_enum.Median:
                                break;
                            default:
                                throw new Exception();
                        }
                    }
                }
            }
        }

        private void Reverse_ratioAVA_and_log2RatioAVG_for_ti()
        {
            int data_length = Data.Length;
            KPMP_subsegmental_proteomics_line_class data_line;
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                data_line = this.Data[indexData];
                if (data_line.Subsegment.Equals("TI"))
                {
                    switch (data_line.Value_type)
                    {
                        case KPMP_value_type_enum.Log2_ratioavg:
                            data_line.Value = -data_line.Value;
                            break;
                        case KPMP_value_type_enum.Ratioavg:
                            data_line.Value = (double)1 / data_line.Value;
                            break;
                        case KPMP_value_type_enum.Average:
                        case KPMP_value_type_enum.Median:
                        case KPMP_value_type_enum.Minus_log10_fdr:
                        case KPMP_value_type_enum.Minus_log10_pvalue:
                            break;
                        default:
                            throw new Exception();
                    }
                }
                else if (data_line.Subsegment.Equals("Glom")) { }
                else { throw new Exception(); }
            }
        }

        private void Keep_only_lines_that_indicate_ratioAVGs()
        {
            int data_length = Data.Length;
            KPMP_subsegmental_proteomics_line_class data_line;
            List<KPMP_subsegmental_proteomics_line_class> keep = new List<KPMP_subsegmental_proteomics_line_class>();
            for (int indexD = 0; indexD < data_length; indexD++)
            {
                data_line = this.Data[indexD];
                if (data_line.SampleName.ToUpper().IndexOf("RATIO") == 0)
                {
                    keep.Add(data_line);
                }
            }
            this.Data = keep.ToArray();
        }

        private void Keep_only_lines_that_indicate_minusLog10Pvalues()
        {
            int data_length = Data.Length;
            KPMP_subsegmental_proteomics_line_class data_line;
            List<KPMP_subsegmental_proteomics_line_class> keep = new List<KPMP_subsegmental_proteomics_line_class>();
            for (int indexD = 0; indexD < data_length; indexD++)
            {
                data_line = this.Data[indexD];
                if (data_line.Value_type.Equals(KPMP_value_type_enum.Minus_log10_pvalue))
                {
                    keep.Add(data_line);
                }
            }
            this.Data = keep.ToArray();
        }

        private void Keep_only_lines_that_indicate_selected_valueTypes(KPMP_value_type_enum[] selected_value_types)
        {
            int data_length = Data.Length;
            KPMP_subsegmental_proteomics_line_class data_line;
            List<KPMP_subsegmental_proteomics_line_class> keep = new List<KPMP_subsegmental_proteomics_line_class>();
            for (int indexD = 0; indexD < data_length; indexD++)
            {
                data_line = this.Data[indexD];
                if (selected_value_types.Contains(data_line.Value_type))
                {
                    keep.Add(data_line);
                }
            }
            this.Data = keep.ToArray();
        }

        private void Check_fo_duplicates()
        {
            int data_length = this.Data.Length;
            KPMP_subsegmental_proteomics_line_class proteomics_line;
            KPMP_subsegmental_proteomics_line_class previous_proteomics_line;
            this.Data = KPMP_subsegmental_proteomics_line_class.Order_by_set(Data).OrderBy(l => l.Gene_symbol).ToArray();
            for (int indexD = 0; indexD < data_length; indexD++)
            {
                proteomics_line = this.Data[indexD];
                if (indexD != 0)
                {
                    previous_proteomics_line = this.Data[indexD - 1];
                    if ((proteomics_line.Equal_set(previous_proteomics_line))
                       && (proteomics_line.Gene_symbol.Equals(previous_proteomics_line.Gene_symbol)))
                    {
                        throw new Exception();
                    }

                }
            }
        }

        private void Check_for_duplicates_based_on_names_for_de_instance()
        {
            int data_length = this.Data.Length;
            KPMP_subsegmental_proteomics_line_class proteomics_line;
            this.Data = this.Data.OrderBy(l => l.Value_type).ThenBy(l => l.Dataset).ThenBy(l => l.PatientId).ThenBy(l => l.Subsegment).ThenBy(l => l.Kpmp_data_integration_term).ThenBy(l => l.Gene_symbol).ToArray();
            for (int indexD = 0; indexD < data_length; indexD++)
            {
                proteomics_line = this.Data[indexD];
                if ((indexD != 0)
                    && (proteomics_line.Value_type.Equals(this.Data[indexD - 1].Value_type))
                    && (proteomics_line.Dataset.Equals(this.Data[indexD - 1].Dataset))
                    && (proteomics_line.PatientId.Equals(this.Data[indexD - 1].PatientId))
                    && (proteomics_line.Subsegment.Equals(this.Data[indexD - 1].Subsegment))
                    && (proteomics_line.Kpmp_data_integration_term.Equals(this.Data[indexD - 1].Kpmp_data_integration_term))
                    && (proteomics_line.Gene_symbol.Equals(this.Data[indexD - 1].Gene_symbol)))
                {
                    throw new Exception();
                }
            }
        }

        private void Check_if_each_gene_is_only_glom_or_TI_and_that_log2Ratioavg_minusLog10Pvalue_and_minusLog10FDR_are_largerEqual_zero()
        {
            int data_length = this.Data.Length;
            KPMP_subsegmental_proteomics_line_class proteomics_line;
            this.Data = Data.OrderBy(l => l.Dataset).ThenBy(l => l.Description).ThenBy(l => l.Accession).ThenBy(l=>l.PatientId).ToArray();
            bool current_glom_logratioavg_set = false;
            bool current_ti_logratioavg_set = false;
            bool current_glom_minusLog10Pvalue_set = false;
            bool current_ti_minusLog10Pvalue_set = false;
            bool current_glom_minusLog10FDR_set = false;
            bool current_ti_minusLog10FDR_set = false;
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                proteomics_line = this.Data[indexData];
                if ((indexData == 0)
                    || (!proteomics_line.Dataset.Equals(this.Data[indexData - 1].Dataset))
                    || (!proteomics_line.Description.Equals(this.Data[indexData - 1].Description))
                    || (!proteomics_line.Accession.Equals(this.Data[indexData - 1].Accession))
                    || (!proteomics_line.PatientId.Equals(this.Data[indexData - 1].PatientId)))
                {
                    current_glom_logratioavg_set = false;
                    current_ti_logratioavg_set = false;
                    current_glom_minusLog10Pvalue_set = false;
                    current_ti_minusLog10Pvalue_set = false;
                    current_glom_minusLog10FDR_set = false;
                    current_ti_minusLog10FDR_set = false;
                }
                switch (proteomics_line.Value_type)
                {
                    case KPMP_value_type_enum.Log2_ratioavg:
                    case KPMP_value_type_enum.Minus_log10_pvalue:
                    case KPMP_value_type_enum.Minus_log10_fdr:
                        if (proteomics_line.Value<0) { throw new Exception(); }
                        break;
                    default:
                        throw new Exception();

                }

                switch (proteomics_line.Subsegment.ToUpper())
                {
                    case "GLOM":
                        switch (proteomics_line.Value_type)
                        {
                            case KPMP_value_type_enum.Log2_ratioavg:
                                if (current_glom_logratioavg_set) { throw new Exception(); }
                                current_glom_logratioavg_set = true;
                                break;
                            case KPMP_value_type_enum.Minus_log10_pvalue:
                                if (current_glom_minusLog10Pvalue_set) { throw new Exception(); }
                                current_glom_minusLog10Pvalue_set = true;
                                break;
                            case KPMP_value_type_enum.Minus_log10_fdr:
                                if (current_glom_minusLog10FDR_set) { throw new Exception(); }
                                current_glom_minusLog10FDR_set = true;
                                break;
                            default:
                                break;
                        }
                        break;
                    case "TI":
                        switch (proteomics_line.Value_type)
                        {
                            case KPMP_value_type_enum.Log2_ratioavg:
                                if (current_ti_logratioavg_set) { throw new Exception(); }
                                current_ti_logratioavg_set = true;
                                break;
                            case KPMP_value_type_enum.Minus_log10_pvalue:
                                if (current_ti_minusLog10Pvalue_set) { throw new Exception(); }
                                current_ti_minusLog10Pvalue_set = true;
                                break;
                            case KPMP_value_type_enum.Minus_log10_fdr:
                                if (current_ti_minusLog10FDR_set) { throw new Exception(); }
                                current_ti_minusLog10FDR_set = true;
                                break;
                            default:
                                break;
                        }
                        break;
                    default:
                        throw new Exception();
                }
                if ((indexData == data_length-1)
                    || (!proteomics_line.Dataset.Equals(this.Data[indexData + 1].Dataset))
                    || (!proteomics_line.Description.Equals(this.Data[indexData + 1].Description))
                    || (!proteomics_line.Accession.Equals(this.Data[indexData + 1].Accession))
                    || (!proteomics_line.PatientId.Equals(this.Data[indexData + 1].PatientId)))
                {
                    if ((current_glom_minusLog10Pvalue_set) && (current_ti_minusLog10Pvalue_set)) { throw new Exception(); }
                    //if ((current_glom_minusLog10FDR_set) && (current_ti_minusLog10FDR_set)) { throw new Exception(); }
                    if ((current_glom_logratioavg_set) && (current_ti_logratioavg_set)) { throw new Exception(); }
                    if ((!current_glom_minusLog10Pvalue_set) && (!current_ti_minusLog10Pvalue_set)) { throw new Exception(); }
                    //if ((!current_glom_minusLog10FDR_set) && (!current_ti_minusLog10FDR_set)) { throw new Exception(); }
                    if ((!current_glom_logratioavg_set) && (!current_ti_logratioavg_set)) { throw new Exception(); }
                }
            }
        }

        private void Remove_lines_with_empty_gene_symbols()
        {
            List<KPMP_subsegmental_proteomics_line_class> keep = new List<KPMP_subsegmental_proteomics_line_class>();
            foreach (KPMP_subsegmental_proteomics_line_class line in this.Data)
            {
                if (!String.IsNullOrEmpty(line.Gene_symbol))
                {
                    keep.Add(line);
                }
            }
            this.Data = keep.ToArray();
        }

        private void Remove_lines_with_empty_descriptions()
        {
            List<KPMP_subsegmental_proteomics_line_class> keep = new List<KPMP_subsegmental_proteomics_line_class>();
            foreach (KPMP_subsegmental_proteomics_line_class line in this.Data)
            {
                if (!String.IsNullOrEmpty(line.Description))
                {
                    keep.Add(line);
                }
                else
                {
                    //string ist = "Ist";
                }
            }
            this.Data = keep.ToArray();
        }

        public void Reverse_foldChanges_and_log2foldChanges_for_ti()
        {
            foreach (KPMP_subsegmental_proteomics_line_class subseg_line in this.Data)
            {
                switch (subseg_line.Subsegment)
                {
                    case "Glom":
                        switch (subseg_line.Value_type)
                        {
                            case KPMP_value_type_enum.Log2_ratioavg:
                            case KPMP_value_type_enum.Log2_ratio_singlepatient:
                            case KPMP_value_type_enum.Ratioavg:
                            case KPMP_value_type_enum.Ratio_singlepatient:
                            case KPMP_value_type_enum.Minus_log10_pvalue:
                            case KPMP_value_type_enum.Minus_log10_fdr:
                                break;
                            default:
                                throw new Exception();
                        }
                        break;
                    case "TI":
                        switch (subseg_line.Value_type)
                        {
                            case KPMP_value_type_enum.Log2_ratioavg:
                            case KPMP_value_type_enum.Log2_ratio_singlepatient:
                                subseg_line.Value = -subseg_line.Value;
                                break;
                            case KPMP_value_type_enum.Ratioavg:
                            case KPMP_value_type_enum.Ratio_singlepatient:
                                subseg_line.Value = 1F / subseg_line.Value;
                                break;
                            case KPMP_value_type_enum.Minus_log10_pvalue:
                            case KPMP_value_type_enum.Minus_log10_fdr:
                                break;
                            default:
                                throw new Exception();
                        }
                        break;
                    default:
                        throw new Exception();
                }
            }
        }

        private void Set_all_symbols_to_upperCase()
        {
            foreach (KPMP_subsegmental_proteomics_line_class subseg_prot_line in this.Data)
            {
                subseg_prot_line.Gene_symbol = subseg_prot_line.Gene_symbol.ToUpper();
            }
        }

        public void Generate_second_part_process_data_for_leavePatientsOut_including_reversal_of_foldChanges_for_ti()
        {
            Remove_all_singleValues();
            Specify_subsegment_based_on_fold_change();
            Keep_only_indicated_valueTypes(new KPMP_value_type_enum[] { KPMP_value_type_enum.Minus_log10_pvalue, KPMP_value_type_enum.Log2_ratioavg, KPMP_value_type_enum.Minus_log10_fdr });
            Remove_duplicated_lines_based_on_dataset_geneSymbol_sampleName_by_keeping_higher_selected_value_type(KPMP_value_type_enum.Minus_log10_pvalue);
            Check_fo_duplicates();
            Set_kpmp_integration_term();
            Check_for_duplicates_based_on_names_for_de_instance();
            Reverse_foldChanges_and_log2foldChanges_for_ti();
            Check_if_each_gene_is_only_glom_or_TI_and_that_log2Ratioavg_minusLog10Pvalue_and_minusLog10FDR_are_largerEqual_zero();
            Keep_only_lines_with_minimum_selectedValueType_and_all_other_values_of_the_same_lines(-(float)Math.Log10(0.05), KPMP_value_type_enum.Minus_log10_fdr);
           // Keep_only_lines_with_minimum_log2foldChange_and_all_other_values_of_the_same_lines(0);
        }

        private void Remove_spike_accessions()
        {
            int subsegmental_prot_length = this.Data.Length;
            KPMP_subsegmental_proteomics_line_class subsegmental_proteomics_line;
            List<KPMP_subsegmental_proteomics_line_class> keep = new List<KPMP_subsegmental_proteomics_line_class>();
            for (int indexS=0; indexS < subsegmental_prot_length; indexS++)
            {
                subsegmental_proteomics_line = this.Data[indexS];
                if (!subsegmental_proteomics_line.Accession.Equals("sp"))
                {
                    keep.Add(subsegmental_proteomics_line);
                }
            }
            this.Data = keep.ToArray();
        }

        private void Remove_space_from_geneSymbols_and_check_official_ncbi_gene_symbols_by_comparing_with_uniprot_gene_names1()
        {
            int subsegmental_prot_length = this.Data.Length;
            KPMP_subsegmental_proteomics_line_class subsegmental_proteomics_line;
            List<string> accessions_with_mismatching_gene_symbols = new List<string>();
            List<string> geneSymbols_with_mismatching_gene_symbols = new List<string>();
            List<string> gn_geneSymbols_with_mismatching_gene_symbols = new List<string>();
            List<string> accessions_without_gn_gene_symbols = new List<string>();
            List<string> geneSymbols_without_gn_gene_symbols = new List<string>();
            List<string> geneSymbol_with_space = new List<string>();
            int indexGN = -1;
            for (int indexSP=0; indexSP<subsegmental_prot_length;indexSP++)
            {
                subsegmental_proteomics_line = this.Data[indexSP];
                int length = subsegmental_proteomics_line.Gene_symbol.Length;
                subsegmental_proteomics_line.Gene_symbol = Text2_class.Remove_space_comma_semicolon_colon_underline_from_end_and_beginning_of_text(subsegmental_proteomics_line.Gene_symbol);
                if (subsegmental_proteomics_line.Gene_symbol.Length!=length)
                {
                    geneSymbol_with_space.Add(subsegmental_proteomics_line.Gene_symbol);
                }
                indexGN = subsegmental_proteomics_line.Description.IndexOf("GN=");
                if (indexGN==-1)
                {
                    accessions_without_gn_gene_symbols.Add(subsegmental_proteomics_line.Accession);
                    geneSymbols_without_gn_gene_symbols.Add(subsegmental_proteomics_line.Gene_symbol);
                }
                else
                {
                    string substring = subsegmental_proteomics_line.Description.Substring(indexGN+3, subsegmental_proteomics_line.Description.Length - indexGN-3);
                    string geneSymbol_from_gn = substring.Split(' ')[0];
                    if (!subsegmental_proteomics_line.Gene_symbol.Equals(geneSymbol_from_gn))
                    {
                        accessions_with_mismatching_gene_symbols.Add(subsegmental_proteomics_line.Accession);
                        geneSymbols_with_mismatching_gene_symbols.Add(subsegmental_proteomics_line.Gene_symbol);
                        gn_geneSymbols_with_mismatching_gene_symbols.Add(geneSymbol_from_gn);
                    }
                }
            }
            geneSymbols_without_gn_gene_symbols = geneSymbols_without_gn_gene_symbols.Distinct().OrderBy(l => l).ToList();
            accessions_without_gn_gene_symbols = accessions_without_gn_gene_symbols.Distinct().OrderBy(l => l).ToList();
            accessions_with_mismatching_gene_symbols = accessions_with_mismatching_gene_symbols.Distinct().OrderBy(l => l).ToList();
            geneSymbols_with_mismatching_gene_symbols = geneSymbols_with_mismatching_gene_symbols.Distinct().ToList();
            gn_geneSymbols_with_mismatching_gene_symbols = gn_geneSymbols_with_mismatching_gene_symbols.Distinct().ToList();
            if (accessions_with_mismatching_gene_symbols.Count>0) { throw new Exception(); }
        }

        private void Add_geneSymbols_from_description_gn_and_remove_lines_without_gn_symbol()
        {
            int subsegmental_prot_length = this.Data.Length;
            KPMP_subsegmental_proteomics_line_class subsegmental_proteomics_line;
            List<KPMP_subsegmental_proteomics_line_class> keep = new List<KPMP_subsegmental_proteomics_line_class>();
            int indexGN = -1;
            string substring;
            string geneSymbol_from_gn;
            for (int indexSP = 0; indexSP < subsegmental_prot_length; indexSP++)
            {
                subsegmental_proteomics_line = this.Data[indexSP];
                indexGN = subsegmental_proteomics_line.Description.IndexOf("GN=");
                if (indexGN != -1)
                {
                    substring = subsegmental_proteomics_line.Description.Substring(indexGN + 3, subsegmental_proteomics_line.Description.Length - indexGN - 3);
                    geneSymbol_from_gn = substring.Split(' ')[0];
                    subsegmental_proteomics_line.Gene_symbol = (string)geneSymbol_from_gn.Clone();
                    keep.Add(subsegmental_proteomics_line);
                }
            }
            this.Data = keep.ToArray();
        }


        //private void Update_official_ncbi_gene_symbols_by_using_my_my_uniprot_version_and_document()
        //{
        //    KPMP_proteomic_exchangeGeneSymbol_documentation_line_class documentation_line;
        //    List<KPMP_proteomic_exchangeGeneSymbol_documentation_line_class> documentations = new List<KPMP_proteomic_exchangeGeneSymbol_documentation_line_class>();
        //    string comparison = "GN from Description in LMD proteomics vs MSSM GN from Uniprot";

        //    UniProt_protein_entries_class uniprot = new UniProt_protein_entries_class();
        //    //uniprot.Generate_de_novo_and_write(Organism_enum.Homo_sapiens);
        //    uniprot.Read();

        //    int uniprot_length = uniprot.UniProts.Length;
        //    UniProt_protein_entries_line_class uniprot_line = new UniProt_protein_entries_line_class();
        //    int indexUniprot = 0;
        //    uniprot.UniProts = uniprot.UniProts.OrderBy(l => l.Accession_number).ThenBy(l => l.Protein_name).ToArray();

        //    KPMP_subsegmental_proteomics_line_class subsegmental_proteomics_line;
        //    int subsegmental_prot_length = this.Data.Length;
        //    this.Data = this.Data.OrderBy(l => l.Accession).ThenBy(l => l.Description).ThenBy(l => l.Gene_symbol).ToArray();
        //    int stringCompare = -2;

        //    List<KPMP_subsegmental_proteomics_line_class> missing_uniprot_ids = new List<KPMP_subsegmental_proteomics_line_class>();

        //    for (int indexSub=0; indexSub<subsegmental_prot_length;indexSub++)
        //    {
        //        subsegmental_proteomics_line = this.Data[indexSub];
        //        stringCompare = -2;
        //        while (stringCompare<0)
        //        {
        //            uniprot_line = uniprot.UniProts[indexUniprot];
        //            if ((indexUniprot!=0)&&(uniprot_line.Accession_number.Equals(uniprot.UniProts[indexUniprot-1].Accession_number))) { throw new Exception(); }
        //            if (!uniprot_line.Organism.Equals(Organism_enum.Homo_sapiens)) { throw new Exception(); }
        //            stringCompare = uniprot_line.Accession_number.CompareTo(subsegmental_proteomics_line.Accession);
        //            if (stringCompare<0)
        //            {
        //                indexUniprot++;
        //            }
        //        }
        //        if (stringCompare != 0)
        //        {
        //            documentation_line = new KPMP_proteomic_exchangeGeneSymbol_documentation_line_class();
        //            documentation_line.Accession = (string)subsegmental_proteomics_line.Accession.Clone();
        //            documentation_line.Comparison = (string)comparison.Clone();
        //            documentation_line.Dataset = (string)subsegmental_proteomics_line.Dataset.Clone();
        //            documentation_line.Full_gene_name_mssm = "";
        //            documentation_line.Gene_symbol_mssm = "not in MSSM Uniprot database";
        //            documentation_line.Full_gene_name_proteomicCore = (string)subsegmental_proteomics_line.Full_gene_name_from_description.Clone();
        //            documentation_line.Gene_symbol_proteomicCore = (string)subsegmental_proteomics_line.Gene_symbol.Clone();

        //            documentation_line.Final_gene_symbol_origin = "Uniprot LMD proteomics";
        //            documentation_line.Final_gene_symbol = (string)subsegmental_proteomics_line.Gene_symbol.Clone();
        //            documentations.Add(documentation_line);
        //        }
        //        else if (uniprot_line.Gene_name.Equals(Global_class.Empty_entry))
        //        {
        //            documentation_line = new KPMP_proteomic_exchangeGeneSymbol_documentation_line_class();
        //            documentation_line.Accession = (string)subsegmental_proteomics_line.Accession.Clone();
        //            documentation_line.Comparison = (string)comparison.Clone();
        //            documentation_line.Dataset = (string)subsegmental_proteomics_line.Dataset.Clone();
        //            documentation_line.Full_gene_name_mssm = "";
        //            documentation_line.Gene_symbol_mssm = "no GN in MSSM Uniprot database";
        //            documentation_line.Full_gene_name_proteomicCore = (string)subsegmental_proteomics_line.Full_gene_name_from_description.Clone();
        //            documentation_line.Gene_symbol_proteomicCore = (string)subsegmental_proteomics_line.Gene_symbol.Clone();

        //            documentation_line.Final_gene_symbol_origin = "Uniprot LMD proteomics";
        //            documentation_line.Final_gene_symbol = (string)subsegmental_proteomics_line.Gene_symbol.Clone();
        //            documentations.Add(documentation_line);
        //        }
        //        else if (!subsegmental_proteomics_line.Gene_symbol.ToUpper().Equals(uniprot_line.Gene_name.ToUpper()))
        //        {
        //            documentation_line = new KPMP_proteomic_exchangeGeneSymbol_documentation_line_class();
        //            documentation_line.Accession = (string)subsegmental_proteomics_line.Accession.Clone();
        //            documentation_line.Comparison = (string)comparison.Clone();
        //            documentation_line.Dataset = (string)subsegmental_proteomics_line.Dataset.Clone();
        //            documentation_line.Gene_symbol_mssm = (string)uniprot_line.Gene_name.Clone();
        //            documentation_line.Full_gene_name_mssm = (string)uniprot_line.Protein_name.Clone();
        //            documentation_line.Gene_symbol_proteomicCore = (string)subsegmental_proteomics_line.Gene_symbol.Clone();
        //            documentation_line.Full_gene_name_proteomicCore = (string)subsegmental_proteomics_line.Full_gene_name_from_description.Clone();
        //            if (subsegmental_proteomics_line.Gene_symbol.IndexOf("SEP") == -1)
        //            {
        //                documentation_line.Final_gene_symbol_origin = "Uniprot LMD proteomics";
        //                documentation_line.Final_gene_symbol = (string)subsegmental_proteomics_line.Gene_symbol.Clone();
        //                documentations.Add(documentation_line);
        //            }
        //            else
        //            {
        //                documentation_line.Final_gene_symbol_origin = "Uniprot LMD proteomics";
        //                documentation_line.Final_gene_symbol = (string)subsegmental_proteomics_line.Gene_symbol.Clone();
        //                documentation_line.Comment = "Sept was recently updated to Septin - older name was kept to ensure matching with other databases";
        //                documentations.Add(documentation_line);
        //            }
        //        }
        //        else if (subsegmental_proteomics_line.Gene_symbol.ToUpper().Equals(uniprot_line.Gene_name.ToUpper()))
        //        { }
        //        else
        //        { 
        //            throw new Exception();
        //        }
        //        //else if (subsegmental_proteomics_line.Full_gene_name_from_description.IndexOf("Uncharacterized protein")!=-1) { }
        //        //else if ((subsegmental_proteomics_line.Gene_symbol.Equals("NONE")) && (uniprot_line.Gene_name.Equals(Global_class.Empty_entry))) { }
        //        //else if ((subsegmental_proteomics_line.Gene_symbol.Equals("none")) && (uniprot_line.Gene_name.Equals(Global_class.Empty_entry))) { }
        //        //else if ((subsegmental_proteomics_line.Gene_symbol.Equals("NA")) && (uniprot_line.Gene_name.Equals(Global_class.Empty_entry))) { }
        //        //else if ((subsegmental_proteomics_line.Gene_symbol.Equals("na")) && (uniprot_line.Gene_name.Equals(Global_class.Empty_entry))) { }
        //        //else if (subsegmental_proteomics_line.Gene_symbol.Equals("Tetrspanin".ToUpper())) { }
        //        //else if (subsegmental_proteomics_line.Description.Equals("CD59 glycoprotein OS=Homo sapiens OX=9606 PE=1 SV=1")) { }
        //        //else if (subsegmental_proteomics_line.Description.Equals("2-oxoisovalerate dehydrogenase subunit alpha OS=Homo sapiens OX=9606 PE=3 SV=1")) { }
        //        //else if (subsegmental_proteomics_line.Description.Equals("Protein O-GlcNAcase OS=Homo sapiens OX=9606 GN=MGEA5 PE=1 SV=2")) { }
        //        //else if (subsegmental_proteomics_line.Description.Equals("ATP synthase subunit d, mitochondrial OS=Homo sapiens OX=9606 GN=ATP5H PE=1 SV=3")) { }
        //        //else if (subsegmental_proteomics_line.Description.Equals("ATP synthase subunit g, mitochondrial OS=Homo sapiens OX=9606 GN=ATP5L PE=1 SV=3")) { }
        //        //else if (subsegmental_proteomics_line.Description.Equals("ATP synthase-coupling factor 6, mitochondrial OS=Homo sapiens OX=9606 GN=ATP5J PE=1 SV=1")) { }
        //        //else if (subsegmental_proteomics_line.Description.Equals("ATP synthase F(0) complex subunit B1, mitochondrial OS=Homo sapiens OX=9606 GN=ATP5F1 PE=1 SV=2")) { }
        //        //else if (subsegmental_proteomics_line.Description.Equals("ATP synthase subunit O, mitochondrial OS=Homo sapiens OX=9606 GN=ATP5O PE=1 SV=1")) { }
        //        //else if (subsegmental_proteomics_line.Description.Equals("Up-regulated during skeletal muscle growth protein 5 OS=Homo sapiens OX=9606 GN=USMG5 PE=1 SV=1")) { }
        //        //else if (subsegmental_proteomics_line.Description.Equals("ATP synthase subunit O, mitochondrial OS=Homo sapiens OX=9606 GN=ATP5O PE=1 SV=1")) { }
        //        //else if (subsegmental_proteomics_line.Description.Equals("ATP synthase subunit O, mitochondrial OS=Homo sapiens OX=9606 GN=ATP5O PE=1 SV=1")) { }
        //        //else if (subsegmental_proteomics_line.Description.Equals("ATP synthase subunit O, mitochondrial OS=Homo sapiens OX=9606 GN=ATP5O PE=1 SV=1")) { }
        //        //else if (!subsegmental_proteomics_line.Gene_symbol.ToUpper().Equals(uniprot_line.Gene_name.ToUpper())) { throw new Exception(); }
        //    }
        //    KPMP_proteomic_exchangeGeneSymbol_documentation_class proteomic_exchangeGeneSymbol = new KPMP_proteomic_exchangeGeneSymbol_documentation_class();
        //    proteomic_exchangeGeneSymbol.Generate_by_adding_to_array(documentations.ToArray());
        //    proteomic_exchangeGeneSymbol.Write("Update_of_LMD_gene_symbols.txt");
        //}

        private void Generate_first_part_after_reading_data_for_both_analyses(string subdirectory, string dataset)
        {
            Read_and_fill_array(subdirectory, dataset);
            //Calculate_FDR();
            Add_subsegment_valueType_patientID_harmonize_statisticsSampleNames_replace_pvalue_by_minusLog10pvalue_and_replace_fdr_by_minusLog10FDR();
            Remove_spike_accessions();
            Add_geneSymbols_from_description_gn_and_remove_lines_without_gn_symbol();
            Generate_dataset_patient_instance();
            Set_all_symbols_to_upperCase();
            Remove_lines_with_empty_gene_symbols();
            Generate_bg_proteins_in_upper_case();
        }

        public void Generate_second_part_for_leaveOneOutAnalysis_including_reversal_of_fold_changes_for_ti()
        {
            Remove_all_singleValues();
            Remove_duplicated_lines_based_on_dataset_geneSymbol_sampleName_by_keeping_higher_selected_value_type(KPMP_value_type_enum.Minus_log10_pvalue);
            Specify_subsegment_based_on_fold_change();
            Reverse_ratioAVA_and_log2RatioAVG_for_ti();
            Duplicate_lines_and_add_cell_specificity();
            Check_fo_duplicates();
            Set_kpmp_integration_term();
            Check_for_duplicates_based_on_names_for_de_instance();
        }

        public void Generate_firstPart_generation_for_leaveOneOut_analysis_and_keep_only_singleValues(string subdirectory, string dataset)
        {
            Generate_first_part_after_reading_data_for_both_analyses(subdirectory, dataset);
            Keep_only_indicated_valueTypes(new KPMP_value_type_enum[] { KPMP_value_type_enum.Single_value });
        }

        private void Calculate_log2ratioavg_from_ratioavg()
        {
            int data_length = Data.Length;
            KPMP_subsegmental_proteomics_line_class data_line;
            KPMP_subsegmental_proteomics_line_class new_data_line;
            List<KPMP_subsegmental_proteomics_line_class> new_lines = new List<KPMP_subsegmental_proteomics_line_class>();
            this.Data = this.Data.OrderBy(l => l.Dataset).ThenBy(l => l.Accession).ThenBy(l => l.Description).ThenBy(l => l.Subsegment).ToArray();
            bool log2ratioavg_calculated_and_added = false;
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                data_line = this.Data[indexData];
                if ((indexData == 0)
                    || (!data_line.Dataset.Equals(this.Data[indexData - 1].Dataset))
                    || (!data_line.Accession.Equals(this.Data[indexData - 1].Accession))
                    || (!data_line.Description.Equals(this.Data[indexData - 1].Description)))
                {
                    log2ratioavg_calculated_and_added = false;
                }
                if (data_line.Value_type.Equals(KPMP_value_type_enum.Ratioavg))
                {
                    if (log2ratioavg_calculated_and_added) { throw new Exception(); }
                    new_data_line = data_line.Deep_copy();
                    new_data_line.Value = Math.Log(new_data_line.Value, 2);
                    new_data_line.Value_type = KPMP_value_type_enum.Log2_ratioavg;
                    new_lines.Add(new_data_line);
                    log2ratioavg_calculated_and_added = true;
                }
                if ((indexData == data_length - 1)
                    || (!data_line.Dataset.Equals(this.Data[indexData + 1].Dataset))
                    || (!data_line.Accession.Equals(this.Data[indexData + 1].Accession))
                    || (!data_line.Description.Equals(this.Data[indexData + 1].Description)))
                {
                    if (!log2ratioavg_calculated_and_added) { throw new Exception(); }
                }
            }
            Add_to_array(new_lines.ToArray());
        }

        private void Set_segment_to_eihter_glom_or_ti_based_on_avgRatio_and_logAvgRatio()
        {
            int data_length = Data.Length;
            KPMP_subsegmental_proteomics_line_class data_line;
            KPMP_subsegmental_proteomics_line_class inner_data_line;
            int firstIndexSame_specification = -1;
            this.Data = this.Data.OrderBy(l => l.Dataset).ThenBy(l => l.Accession).ThenBy(l => l.Description).ThenBy(l => l.Subsegment).ToArray();
            string subsegment = "not_sepcified";
            string current_subsegment = "";

            for (int indexData = 0; indexData < data_length; indexData++)
            {
                data_line = this.Data[indexData];
                if ((indexData == 0)
                    || (!data_line.Dataset.Equals(this.Data[indexData - 1].Dataset))
                    || (!data_line.Accession.Equals(this.Data[indexData - 1].Accession))
                    || (!data_line.Description.Equals(this.Data[indexData - 1].Description)))
                {
                    subsegment = "not_yet_specified";
                    firstIndexSame_specification = indexData;
                }
                switch (data_line.Value_type)
                {
                    case KPMP_value_type_enum.Log2_ratioavg:
                        if (data_line.Value < 0) { current_subsegment = "TI"; }
                        else { current_subsegment = "Glom"; }
                        break;
                    case KPMP_value_type_enum.Ratioavg:
                        if (data_line.Value < 1) { data_line.Subsegment = "TI"; }
                        else { data_line.Subsegment = "Glom"; }
                        break;
                    case KPMP_value_type_enum.Minus_log10_pvalue:
                    case KPMP_value_type_enum.Minus_log10_fdr:
                        break;
                    default:
                        throw new Exception();
                }
                if (subsegment.Equals("not_yet_specified"))
                {
                    subsegment = current_subsegment;
                }
                else if (!subsegment.Equals(current_subsegment)) { throw new Exception(); }
                if ((indexData == data_length-1)
                    || (!data_line.Dataset.Equals(this.Data[indexData + 1].Dataset))
                    || (!data_line.Accession.Equals(this.Data[indexData + 1].Accession))
                    || (!data_line.Description.Equals(this.Data[indexData + 1].Description)))
                {
                    if (subsegment.Equals("")) { throw new Exception(); }
                    for (int indexInner=firstIndexSame_specification;indexInner<=indexData;indexInner++)
                    {
                        inner_data_line = this.Data[indexInner];
                        inner_data_line.Subsegment = (string)subsegment.Clone();
                    }
                }
            }
        }

        private void Calculate_FDR()
        {
            int data_length = this.Data.Length;
            KPMP_subsegmental_proteomics_line_class proteomics_line;
            Dictionary<string, bool> accession_dict = new Dictionary<string, bool>();
            int tested_proteins_count = 0;
            for (int indexD=0; indexD<data_length;indexD++)
            {
                proteomics_line = this.Data[indexD];
                if (!accession_dict.ContainsKey(proteomics_line.Description))
                {
                    accession_dict.Add(proteomics_line.Description, true);
                }
                if (proteomics_line.SampleName.ToUpper().IndexOf("P_VALUE") == 0)
                {
                    tested_proteins_count++;
                }
            }


            KPMP_subsegmental_proteomics_line_class fdr_proteomics_line;

            this.Data = this.Data.OrderBy(l => l.Dataset).ThenBy(l => l.SampleName).ThenBy(l=>l.PatientId).ThenBy(l=>l.Subsegment).ThenBy(l => l.Value_type).ThenBy(l => l.Value).ToArray();
            List<KPMP_subsegmental_proteomics_line_class> fdr_values = new List<KPMP_subsegmental_proteomics_line_class>();
            int current_pvalue_rank = -1;
            for (int indexD=0; indexD<data_length;indexD++)
            {
                proteomics_line = this.Data[indexD];
                if ( (indexD==0)
                    || (!proteomics_line.Dataset.Equals(this.Data[indexD - 1].Dataset))
                    || (!proteomics_line.Subsegment.Equals(this.Data[indexD - 1].Subsegment))
                    || (!proteomics_line.PatientId.Equals(this.Data[indexD - 1].PatientId))
                    || (!proteomics_line.SampleName.Equals(this.Data[indexD - 1].SampleName)))
                {
                current_pvalue_rank = 0;
                }
                if ((proteomics_line.SampleName.ToUpper().IndexOf("P_VALUE") == 0))
                {
                    current_pvalue_rank++;
                    fdr_proteomics_line = proteomics_line.Deep_copy();
                    fdr_proteomics_line.SampleName = fdr_proteomics_line.SampleName.Replace("P_Value","FDR");
                    fdr_proteomics_line.Value_type = KPMP_value_type_enum.Fdr;
                    fdr_proteomics_line.Value = fdr_proteomics_line.Value * (tested_proteins_count / current_pvalue_rank);
                    fdr_values.Add(fdr_proteomics_line);
                }
            }
            int fdr_values_length = fdr_values.Count;
            fdr_values = fdr_values.OrderBy(l => l.Dataset).ThenBy(l => l.SampleName).ThenBy(l => l.PatientId).ThenBy(l => l.Subsegment).ThenBy(l => l.Value_type).ThenByDescending(l => l.Value).ToList();
            double current_lowest_fdr = 1;
            int fdr_of_one_set_count = 0;
            for (int indexFDR=0; indexFDR<fdr_values_length;indexFDR++)
            {
                fdr_proteomics_line = fdr_values[indexFDR];
                if (  (indexFDR==0)
                    || (!fdr_proteomics_line.Subsegment.Equals(this.Data[indexFDR - 1].Subsegment))
                    || (!fdr_proteomics_line.PatientId.Equals(this.Data[indexFDR - 1].PatientId))
                    || (!fdr_proteomics_line.Dataset.Equals(fdr_values[indexFDR-1].Dataset))
                    || (!fdr_proteomics_line.SampleName.Equals(fdr_values[indexFDR-1].SampleName)))
                {
                    current_lowest_fdr = 1;
                }
                fdr_of_one_set_count++;
                if (fdr_proteomics_line.Value<current_lowest_fdr) { current_lowest_fdr = fdr_proteomics_line.Value; }
                if (!fdr_proteomics_line.Value_type.Equals(KPMP_value_type_enum.Fdr)) { throw new Exception(); }
                fdr_proteomics_line.Value = current_lowest_fdr;
                fdr_values[indexFDR] = fdr_proteomics_line;
                if ((indexFDR == fdr_values_length-1)
                    || (!fdr_proteomics_line.Subsegment.Equals(this.Data[indexFDR + 1].Subsegment))
                    || (!fdr_proteomics_line.PatientId.Equals(this.Data[indexFDR + 1].PatientId))
                    || (!fdr_proteomics_line.Dataset.Equals(fdr_values[indexFDR + 1].Dataset))
                    || (!fdr_proteomics_line.SampleName.Equals(fdr_values[indexFDR + 1].SampleName)))
                {
                    if (fdr_of_one_set_count!=tested_proteins_count) { throw new Exception(); }
                }
            }
            Add_to_array(fdr_values.ToArray());
        }

        public void Generate_for_regular_analysis(string subdirectory, string dataset)
        {
            Generate_first_part_after_reading_data_for_both_analyses(subdirectory, dataset);
            Remove_all_singleValues();
            Remove_duplicated_lines_based_on_dataset_geneSymbol_sampleName_by_keeping_higher_selected_value_type(KPMP_value_type_enum.Minus_log10_pvalue);

            Calculate_log2ratioavg_from_ratioavg();
            Keep_only_indicated_valueTypes(new KPMP_value_type_enum[] { KPMP_value_type_enum.Minus_log10_pvalue, KPMP_value_type_enum.Log2_ratioavg, KPMP_value_type_enum.Minus_log10_fdr });
            Set_segment_to_eihter_glom_or_ti_based_on_avgRatio_and_logAvgRatio();
            Reverse_foldChanges_and_log2foldChanges_for_ti();
            Check_if_each_gene_is_only_glom_or_TI_and_that_log2Ratioavg_minusLog10Pvalue_and_minusLog10FDR_are_largerEqual_zero();
            Duplicate_lines_and_add_cell_specificity();
            Check_fo_duplicates();
            Set_kpmp_integration_term();
            Check_for_duplicates_based_on_names_for_de_instance();
        }

        public void Replace_all_segment0s_by_segment1s(string segment0, string segment1)
        {
            foreach (KPMP_subsegmental_proteomics_line_class subsegmental_proteomic_line in this.Data)
            {
                if (subsegmental_proteomic_line.Subsegment.Equals(segment1)) { throw new Exception(); }
                if (subsegmental_proteomic_line.Subsegment.Equals(segment0))
                {
                    subsegmental_proteomic_line.Subsegment = (string)segment1.Clone();
                }
            }
        }

        private void Add_lines_with_zscore_normalized_expression_of_each_gene()
        {
            List<double> current_expression_values = new List<double>();
            double population_sd = -1;
            double mean = -1;
            int firstIndex_current_datasetCelltypeValueTypeSymbol = -1;
            int data_length = this.Data.Length;
            KPMP_subsegmental_proteomics_line_class new_data_line;
            List<KPMP_subsegmental_proteomics_line_class> new_data_list = new List<KPMP_subsegmental_proteomics_line_class>();
            KPMP_subsegmental_proteomics_line_class data_line;
            KPMP_subsegmental_proteomics_line_class inner_data_line;
            this.Data = this.Data.OrderBy(l => l.Dataset).ThenBy(l => l.Value_type).ThenBy(l => l.Gene_symbol).ToArray();
            Dictionary<string, bool> considered_cellTypes = new Dictionary<string, bool>();
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                data_line = this.Data[indexData];
                if ((indexData == 0)
                    || (!data_line.Dataset.Equals(this.Data[indexData - 1].Dataset))
                    || (!data_line.Value_type.Equals(this.Data[indexData - 1].Value_type))
                    || (!data_line.Gene_symbol.Equals(this.Data[indexData - 1].Gene_symbol)))
                {
                    firstIndex_current_datasetCelltypeValueTypeSymbol = indexData;
                    current_expression_values.Clear();
                }
                current_expression_values.Add(data_line.Value);
                if ((indexData == data_length - 1)
                    || (!data_line.Dataset.Equals(this.Data[indexData + 1].Dataset))
                    || (!data_line.Value_type.Equals(this.Data[indexData + 1].Value_type))
                    || (!data_line.Gene_symbol.Equals(this.Data[indexData + 1].Gene_symbol)))
                {
                    Math_class.Get_mean_and_population_sd(current_expression_values.ToArray(), out mean, out population_sd);
                    for (int indexInner = firstIndex_current_datasetCelltypeValueTypeSymbol; indexInner <= indexData; indexInner++)
                    {
                        inner_data_line = this.Data[indexInner];
                        new_data_line = inner_data_line.Deep_copy();
                        new_data_line.Value = (new_data_line.Value - mean) / population_sd;
                        switch (new_data_line.Value_type)
                        {
                            case KPMP_value_type_enum.Log2_ratio_singlepatient:
                                new_data_line.Value_type = KPMP_value_type_enum.Scaled_log2_ratio_singlepatient;
                                new_data_list.Add(new_data_line);
                                break;
                            case KPMP_value_type_enum.Ratio_singlepatient:
                                new_data_line.Value_type = KPMP_value_type_enum.Scaled_ratio_singlepatient;
                                new_data_list.Add(new_data_line);
                                break;
                            case KPMP_value_type_enum.Single_value:
                                new_data_line.Value_type = KPMP_value_type_enum.Scaled_single_value;
                                new_data_list.Add(new_data_line);
                                break;
                            default:
                                break;

                        }
                    }
                }
            }
            Add_to_array(new_data_list.ToArray());
        }

        public KPMP_integration_paper_metadata_class Get_deep_copy_of_dataset_patient_instance()
        {
            return Dataset_patient.Deep_copy();
        }

        private void Generate_dataset_patient_instance()
        {
            Dictionary<string, Dictionary<string, bool>> dataset_patient_considered_dict = new Dictionary<string, Dictionary<string, bool>>();
            int data_length = this.Data.Length;
            KPMP_integration_paper_metadata_line_class new_dataset_patient_line;
            List<KPMP_integration_paper_metadata_line_class> new_dataset_patient_lines = new List<KPMP_integration_paper_metadata_line_class>();
            KPMP_subsegmental_proteomics_line_class singleClister_line;
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                singleClister_line = this.Data[indexData];
                if (!dataset_patient_considered_dict.ContainsKey(singleClister_line.Dataset))
                {
                    dataset_patient_considered_dict.Add(singleClister_line.Dataset, new Dictionary<string, bool>());
                }
                if (!dataset_patient_considered_dict[singleClister_line.Dataset].ContainsKey(singleClister_line.PatientId))
                {
                    dataset_patient_considered_dict[singleClister_line.Dataset].Add(singleClister_line.PatientId, true);
                    new_dataset_patient_line = new KPMP_integration_paper_metadata_line_class();
                    new_dataset_patient_line.Dataset = (string)singleClister_line.Dataset.Clone();
                    new_dataset_patient_line.Libraries = new string[] { (string)singleClister_line.PatientId.Clone() };
                    new_dataset_patient_lines.Add(new_dataset_patient_line);
                }
            }
            Dataset_patient = new KPMP_integration_paper_metadata_class();
            Dataset_patient.Add_to_array(new_dataset_patient_lines.ToArray());
        }

        public void Generate_for_single_patient_analysis(string subdirectory, string dataset)
        {
            Read_and_fill_array(subdirectory, dataset);
            //Remove_lines_with_empty_descriptions();
            Add_subsegment_valueType_patientID_harmonize_statisticsSampleNames_replace_pvalue_by_minusLog10pvalue_and_replace_fdr_by_minusLog10FDR();
            Remove_spike_accessions();
            Add_geneSymbols_from_description_gn_and_remove_lines_without_gn_symbol();
            Generate_dataset_patient_instance();
            Set_all_symbols_to_upperCase();
            Remove_lines_with_empty_gene_symbols();
            Generate_bg_proteins_in_upper_case();
            //Specify_subsegment_based_on_fold_change();
            //Generate_new_lines_by_inverting_ratios();
            Keep_only_indicated_valueTypes(KPMP_value_type_enum.Single_value);
            Remove_duplicated_lines_based_on_dataset_geneSymbol_sampleName_by_keeping_higher_selected_value_type(KPMP_value_type_enum.Single_value);
            //Generate_ratios_and_log2ratios_for_single_patients_for_each_geneSymbol();
            //Duplicate_lines_and_add_cell_specificity();
            //Set_kpmp_integration_term();
            Check_for_duplicates_based_on_names_for_de_instance();
            //Keep_only_lines_with_indicated_subsegments(new string[] { "Glom", "TI", "Ratio ProxTub vs Glom", "Ratio Glom vs ProxTub", "Log2Ratio Glom vs ProxTub", "Log2Ratio ProxTub vs Glom" });
            //Add_lines_with_zscore_normalized_expression_of_each_gene();
        }

        public void Keep_only_lines_with_indicated_patientIDs(params string[] patientIDs)
        {
            patientIDs = patientIDs.Distinct().OrderBy(l => l).ToArray();
            int indexPatientID = 0;
            int patientIDs_length = patientIDs.Length;
            string patientID;
            int stringCompare = -2;
            bool patientID_found = false;
            int data_length = this.Data.Length;
            List<KPMP_subsegmental_proteomics_line_class> keep = new List<KPMP_subsegmental_proteomics_line_class>();
            KPMP_subsegmental_proteomics_line_class data_line;
            this.Data = this.Data.OrderBy(l => l.PatientId).ToArray();
            List<string> removed_patients = new List<string>();
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                data_line = this.Data[indexData];
                stringCompare = -2;
                while ((indexPatientID < patientIDs_length) && (stringCompare < 0))
                {
                    patientID = patientIDs[indexPatientID];
                    stringCompare = patientID.CompareTo(data_line.PatientId);
                    if (stringCompare < 0)
                    {
                        if (!patientID_found) { throw new Exception(); }
                        indexPatientID++;
                        patientID_found = false;
                    }
                    else if (stringCompare == 0)
                    {
                        patientID_found = true;
                        keep.Add(data_line);
                    }
                }
                if (stringCompare != 0)
                {
                    removed_patients.Add(data_line.PatientId);
                }
            }
            removed_patients = removed_patients.Distinct().OrderBy(l => l).ToList();
            this.Data = keep.ToArray();
        }

        public void Keep_only_lines_with_indicated_subsegments(params string[] subsegments)
        {
            subsegments = subsegments.Distinct().OrderBy(l => l).ToArray();
            int indexSubsegment = 0;
            int subsegments_length = subsegments.Length;
            string subsegment;
            int stringCompare = -2;

            int data_length = this.Data.Length;
            List<KPMP_subsegmental_proteomics_line_class> keep = new List<KPMP_subsegmental_proteomics_line_class>();
            KPMP_subsegmental_proteomics_line_class data_line;
            this.Data = this.Data.OrderBy(l => l.Subsegment).ToArray();
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                data_line = this.Data[indexData];
                stringCompare = -2;
                while ((indexSubsegment < subsegments_length) && (stringCompare < 0))
                {
                    subsegment = subsegments[indexSubsegment];
                    stringCompare = subsegment.CompareTo(data_line.Subsegment);
                    if (stringCompare < 0)
                    {
                        indexSubsegment++;
                    }
                    else if (stringCompare == 0)
                    {
                        keep.Add(data_line);
                    }
                }
            }
            this.Data = keep.ToArray();
        }

        public void Keep_only_lines_with_indicated_symbols(params string[] symbols)
        {
            symbols = symbols.Distinct().OrderBy(l => l).ToArray();
            int indexSymbol = 0;
            int symbols_length = symbols.Length;
            string symbol;
            int stringCompare = -2;

            int data_length = this.Data.Length;
            List<KPMP_subsegmental_proteomics_line_class> keep = new List<KPMP_subsegmental_proteomics_line_class>();
            KPMP_subsegmental_proteomics_line_class data_line;
            this.Data = this.Data.OrderBy(l => l.Gene_symbol).ToArray();
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                data_line = this.Data[indexData];
                stringCompare = -2;
                while ((indexSymbol < symbols_length) && (stringCompare < 0))
                {
                    symbol = symbols[indexSymbol];
                    stringCompare = symbol.CompareTo(data_line.Gene_symbol);
                    if (stringCompare < 0)
                    {
                        indexSymbol++;
                    }
                    else if (stringCompare == 0)
                    {
                        keep.Add(data_line);
                    }
                }
            }
            this.Data = keep.ToArray();
        }

        public void Keep_only_indicated_valueTypes(params KPMP_value_type_enum[] keep_value_types)
        {
            int data_length = Data.Length;
            KPMP_subsegmental_proteomics_line_class proteomics_line;
            List<KPMP_subsegmental_proteomics_line_class> keep = new List<KPMP_subsegmental_proteomics_line_class>();
            for (int indexP = 0; indexP < data_length; indexP++)
            {
                proteomics_line = this.Data[indexP];
                if (keep_value_types.Contains(proteomics_line.Value_type))
                {
                    keep.Add(proteomics_line);
                }
            }
            this.Data = keep.ToArray();
        }

        public void Keep_only_glom_subsegments_or_minusLog10Pvalues()
        {
            int data_length = Data.Length;
            KPMP_subsegmental_proteomics_line_class proteomics_line;
            List<KPMP_subsegmental_proteomics_line_class> keep = new List<KPMP_subsegmental_proteomics_line_class>();
            for (int indexP = 0; indexP < data_length; indexP++)
            {
                proteomics_line = this.Data[indexP];
                if (  (proteomics_line.Subsegment.Equals("Glom"))
                    ||(proteomics_line.Value_type.Equals(KPMP_value_type_enum.Minus_log10_pvalue)))
                {
                    keep.Add(proteomics_line);
                }
            }
            this.Data = keep.ToArray();
        }

        public string[] Get_all_proteins()
        {
            this.Data = this.Data.OrderBy(l => l.Gene_symbol).ToArray();
            int proteomic_length = this.Data.Length;
            KPMP_subsegmental_proteomics_line_class proteomic_line;
            List<string> all_proteins = new List<string>();
            for (int indexP = 0; indexP < proteomic_length; indexP++)
            {
                proteomic_line = this.Data[indexP];
                if ((indexP == 0) || (!proteomic_line.Gene_symbol.Equals(this.Data[indexP - 1].Gene_symbol)))
                {
                    all_proteins.Add(proteomic_line.Gene_symbol);
                }
            }
            return all_proteins.ToArray();
        }

        public string[] Get_all_subsegments()
        {
            this.Data = this.Data.OrderBy(l => l.Subsegment).ToArray();
            int proteomic_length = this.Data.Length;
            KPMP_subsegmental_proteomics_line_class proteomic_line;
            List<string> all_subsegments = new List<string>();
            for (int indexP = 0; indexP < proteomic_length; indexP++)
            {
                proteomic_line = this.Data[indexP];
                if ((indexP == 0) || (!proteomic_line.Subsegment.Equals(this.Data[indexP - 1].Subsegment)))
                {
                    all_subsegments.Add(proteomic_line.Subsegment);
                }
            }
            return all_subsegments.ToArray();
        }

        public string[] Get_all_patientIDs()
        {
            this.Data = this.Data.OrderBy(l => l.PatientId).ToArray();
            int proteomic_length = this.Data.Length;
            KPMP_subsegmental_proteomics_line_class proteomic_line;
            List<string> all_patientIDs = new List<string>();
            for (int indexP = 0; indexP < proteomic_length; indexP++)
            {
                proteomic_line = this.Data[indexP];
                if ((indexP == 0) || (!proteomic_line.PatientId.Equals(this.Data[indexP - 1].PatientId)))
                {
                    all_patientIDs.Add(proteomic_line.PatientId);
                }
            }
            return all_patientIDs.ToArray();
        }

        public string[] Get_all_singleValue_patientIDs()
        {
            this.Data = this.Data.OrderBy(l => l.PatientId).ToArray();
            int proteomic_length = this.Data.Length;
            KPMP_subsegmental_proteomics_line_class proteomic_line;
            List<string> all_patientIDs = new List<string>();
            for (int indexP = 0; indexP < proteomic_length; indexP++)
            {
                proteomic_line = this.Data[indexP];
                if ((indexP == 0) || (!proteomic_line.PatientId.Equals(this.Data[indexP - 1].PatientId)))
                {
                    if (proteomic_line.Value_type.Equals(KPMP_value_type_enum.Single_value))
                    {
                        all_patientIDs.Add(proteomic_line.PatientId);
                    }
                }
            }
            return all_patientIDs.ToArray();
        }

        public string[] Get_deepCopy_of_bg_proteins_in_upper_case()
        {
            int bg_proteins_length = this.Bg_proteins_in_upper_case.Length;
            string[] bg_proteins = new string[bg_proteins_length];
            for (int indexP = 0; indexP < bg_proteins_length; indexP++)
            {
                bg_proteins[indexP] = (string)this.Bg_proteins_in_upper_case[indexP].Clone();
            }
            return bg_proteins;
        }

        public KPMP_standardized_dataset_class Generate_kpmp_standardized_dataset_instance(KPMP_value_type_enum value_type_1st, KPMP_value_type_enum value_type_2nd)
        {
            KPMP_standardized_dataset_line_class standardized_data_line;
            List<KPMP_standardized_dataset_line_class> standardized_data_list = new List<KPMP_standardized_dataset_line_class>();
            int data_length = this.Data.Length;
            KPMP_subsegmental_proteomics_line_class proteomics_line;
            this.Data = this.Data.OrderBy(l => l.Dataset).ThenBy(l => l.Subsegment).ThenBy(l => l.PatientId).ThenBy(l => l.Kpmp_data_integration_term).ThenBy(l => l.Gene_symbol).ThenBy(l=>l.Value_type).ToArray();
            string dataset = (string)this.Data[0].Dataset.Clone();
            double current_value_1st = -1;
            double current_value_2nd = -1;
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                proteomics_line = this.Data[indexData];
                if (!proteomics_line.Dataset.Equals(dataset)) { throw new Exception(); }
                if ((indexData != 0)
                    && (proteomics_line.Dataset.Equals(this.Data[indexData - 1].Dataset))
                    && (proteomics_line.Subsegment.Equals(this.Data[indexData - 1].Subsegment))
                    && (proteomics_line.PatientId.Equals(this.Data[indexData - 1].PatientId))
                    && (proteomics_line.Kpmp_data_integration_term.Equals(this.Data[indexData - 1].Kpmp_data_integration_term))
                    && (proteomics_line.Gene_symbol.Equals(this.Data[indexData - 1].Gene_symbol))
                    && (proteomics_line.Value_type.Equals(this.Data[indexData - 1].Value_type)))
                {
                    throw new Exception();
                }
                if ((indexData == 0)
                    || (!proteomics_line.Dataset.Equals(this.Data[indexData - 1].Dataset))
                    || (!proteomics_line.Subsegment.Equals(this.Data[indexData - 1].Subsegment))
                    || (!proteomics_line.PatientId.Equals(this.Data[indexData - 1].PatientId))
                    || (!proteomics_line.Kpmp_data_integration_term.Equals(this.Data[indexData - 1].Kpmp_data_integration_term))
                    || (!proteomics_line.Gene_symbol.Equals(this.Data[indexData - 1].Gene_symbol)))
                {
                    current_value_1st = -1;
                    switch (value_type_2nd)
                    {
                        case KPMP_value_type_enum.No_selection:
                            current_value_2nd = 0;
                            break;
                        default:
                            current_value_2nd = -1;
                            break;
                    }
                }
                if (proteomics_line.Value_type.Equals(value_type_1st))
                {
                    if (current_value_1st != -1) { throw new Exception(); }
                    current_value_1st = proteomics_line.Value;
                }

                if (proteomics_line.Value_type.Equals(value_type_2nd))
                {
                    if (current_value_2nd != -1) { throw new Exception(); }
                    current_value_2nd = proteomics_line.Value;
                }
                if ((indexData == data_length-1)
                    || (!proteomics_line.Dataset.Equals(this.Data[indexData + 1].Dataset))
                    || (!proteomics_line.Subsegment.Equals(this.Data[indexData + 1].Subsegment))
                    || (!proteomics_line.PatientId.Equals(this.Data[indexData + 1].PatientId))
                    || (!proteomics_line.Kpmp_data_integration_term.Equals(this.Data[indexData + 1].Kpmp_data_integration_term))
                    || (!proteomics_line.Gene_symbol.Equals(this.Data[indexData + 1].Gene_symbol)))
                {
                    if (current_value_1st == -1) { throw new Exception(); }
                    if (current_value_2nd == -1) { throw new Exception(); }
                    standardized_data_line = new KPMP_standardized_dataset_line_class();
                    standardized_data_line.Cell_segment = (string)proteomics_line.Subsegment.Clone();
                    standardized_data_line.Dataset = (string)proteomics_line.Dataset.Clone();
                    standardized_data_line.Value_1st = current_value_1st;
                    standardized_data_line.Value_2nd = current_value_2nd;
                    standardized_data_line.KPMP_data_integration_term = (string)proteomics_line.Kpmp_data_integration_term.Clone();
                    standardized_data_line.Value_type_1st = value_type_1st;
                    standardized_data_line.Value_type_2nd = value_type_2nd;
                    standardized_data_line.PatientId = (string)proteomics_line.PatientId.Clone();
                    standardized_data_line.Gene_symbol = (string)proteomics_line.Gene_symbol.Clone();
                    standardized_data_list.Add(standardized_data_line);
                }
            }
            Dictionary<string, string[]> dataset_bgProteins_dict = new Dictionary<string, string[]>();
            dataset_bgProteins_dict.Add(dataset, Get_deepCopy_of_bg_proteins_in_upper_case());
            KPMP_standardized_dataset_class standardized_data = new KPMP_standardized_dataset_class();
            standardized_data.Add_to_existing_instances(standardized_data_list.ToArray(), dataset_bgProteins_dict);
            return standardized_data;
        }

        public KPMP_standardized_dataset_class Generate_standardized_dataset_instance_with_all_values_filled_in_value_1st()
        {
            KPMP_standardized_dataset_line_class new_standardized_data_line;
            List<KPMP_standardized_dataset_line_class> standardized_data_list = new List<KPMP_standardized_dataset_line_class>();
            int data_length = this.Data.Length;
            KPMP_subsegmental_proteomics_line_class proteomics_line;
            this.Data = this.Data.OrderBy(l => l.Dataset).ThenBy(l => l.Subsegment).ThenBy(l => l.PatientId).ThenBy(l => l.Kpmp_data_integration_term).ThenBy(l => l.Value_type).ThenBy(l => l.Gene_symbol).ToArray();
            for (int indexR = 0; indexR < data_length; indexR++)
            {
                proteomics_line = this.Data[indexR];
                if ((indexR != 0)
                    && (proteomics_line.Dataset.Equals(this.Data[indexR - 1].Dataset))
                    && (proteomics_line.Subsegment.Equals(this.Data[indexR - 1].Subsegment))
                    && (proteomics_line.PatientId.Equals(this.Data[indexR - 1].PatientId))
                    && (proteomics_line.Kpmp_data_integration_term.Equals(this.Data[indexR - 1].Kpmp_data_integration_term))
                    && (proteomics_line.Gene_symbol.Equals(this.Data[indexR - 1].Gene_symbol))
                    && (proteomics_line.Value_type.Equals(this.Data[indexR - 1].Value_type)))
                {
                    throw new Exception();
                }
                new_standardized_data_line = new KPMP_standardized_dataset_line_class();
                new_standardized_data_line.Cell_segment = (string)proteomics_line.Subsegment.Clone();
                new_standardized_data_line.Dataset = (string)proteomics_line.Dataset.Clone();
                new_standardized_data_line.Gene_symbol = (string)proteomics_line.Gene_symbol.Clone();
                new_standardized_data_line.PatientId = (string)proteomics_line.PatientId.Clone();
                new_standardized_data_line.Value_1st = proteomics_line.Value;
                new_standardized_data_line.Value_2nd = 0;
                new_standardized_data_line.Value_type_1st = proteomics_line.Value_type;
                new_standardized_data_line.Value_type_2nd = KPMP_value_type_enum.No_selection;
                standardized_data_list.Add(new_standardized_data_line);
            }
            string dataset = KPMP_dataset_name_class.Lmd_proteomics_iuosu;
            Dictionary<string, string[]> dataset_bgSymbol_dict = new Dictionary<string, string[]>();
            dataset_bgSymbol_dict.Add(dataset, Get_deepCopy_of_bg_proteins_in_upper_case());
            KPMP_standardized_dataset_class standard = new KPMP_standardized_dataset_class();
            standard.Add_to_existing_instances(standardized_data_list.ToArray(), dataset_bgSymbol_dict);
            return standard;
        }

        private void Read_and_fill_array(string subdirectory, string dataset)
        {
            string complete_directory = Global_directory_class.Experimental_data_directory + subdirectory;
            string[] complete_fileNames = Directory.GetFiles(complete_directory);
            string complete_fileName;
            int complete_fileNames_length = complete_fileNames.Length;
            KPMP_subsegmental_proteomics_line_class new_proteomic_line;
            List<KPMP_subsegmental_proteomics_line_class> data_list = new List<KPMP_subsegmental_proteomics_line_class>();
            for (int indexC = 0; indexC < complete_fileNames_length; indexC++)
            {
                complete_fileName = complete_fileNames[indexC];

                char delimiter = Global_class.Tab;
                StreamReader reader = new StreamReader(complete_fileName);
                string inputLine;
                string[] splitStrings;
                string[] columnNames;
                string description;
                string accession;
                string geneID;
                string expression_value_string;

                inputLine = reader.ReadLine();
                columnNames = inputLine.Split(delimiter);
                int columns_length = columnNames.Length;
                int readLines_count = 0;

                while ((inputLine = reader.ReadLine()) != null)
                {
                    readLines_count++;
                    splitStrings = inputLine.Split(delimiter);
                    if (splitStrings.Length != columns_length) { throw new Exception(); }
                    accession = splitStrings[0];
                    geneID = splitStrings[1];
                    description = splitStrings[2];
                    for (int indexCol = 3; indexCol < columns_length; indexCol++)
                    {
                        new_proteomic_line = new KPMP_subsegmental_proteomics_line_class();
                        new_proteomic_line.Description = (string)description.Clone();
                        new_proteomic_line.Accession = (string)accession.Clone();
                        new_proteomic_line.Gene_symbol = (string)geneID.Clone();
                        new_proteomic_line.SampleName = (string)columnNames[indexCol].Clone();
                        new_proteomic_line.Dataset = (string)dataset.Clone(); 
                        expression_value_string = splitStrings[indexCol];
                        if (expression_value_string.Equals("#DIV/0!"))
                        {
                            new_proteomic_line.Value = double.MaxValue;
                            throw new Exception();
                        }
                        else
                        {
                            new_proteomic_line.Value = double.Parse(splitStrings[indexCol]);
                            data_list.Add(new_proteomic_line);
                        }
                    }
                }
            }
            this.Data = data_list.ToArray();
        }

        public KPMP_subsegmental_proteomics_class Deep_copy()
        {
            int data_length = this.Data.Length;
            KPMP_subsegmental_proteomics_class copy = (KPMP_subsegmental_proteomics_class)this.MemberwiseClone();
            copy.Data = new KPMP_subsegmental_proteomics_line_class[data_length];
            for (int indexD = 0; indexD < data_length; indexD++)
            {
                copy.Data[indexD] = this.Data[indexD].Deep_copy();
            }
            int bg_genes_length = this.Bg_proteins_in_upper_case.Length;
            for (int indexBg = 0; indexBg < bg_genes_length; indexBg++)
            {
                copy.Bg_proteins_in_upper_case[indexBg] = (string)this.Bg_proteins_in_upper_case[indexBg].Clone();
            }

            return copy;
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class KPMP_nearSingleCell_proteomics_UCSF_line_class
    {
        public string Gene_symbol { get; set; }
        public string Subsegment { get; set; }
        public string SampleName { get; set; }
        public string PatientID { get; set; }
        public string Dataset { get; set; }
        public string Swissprot_id { get; set; }
        public string Protein_name { get; set; }
        public string Uniprot_accession { get; set; }
        public string Kpmp_data_integration_term { get; set; }
        public double Value { get; set; }
        public KPMP_value_type_enum Value_type { get; set; }

        #region equal
        public static KPMP_nearSingleCell_proteomics_UCSF_line_class[] Order_by_set(KPMP_nearSingleCell_proteomics_UCSF_line_class[] array)
        {
            array = array.OrderBy(l => l.Dataset).ThenBy(l => l.Subsegment).ThenBy(l => l.SampleName).ThenBy(l => l.PatientID).ThenBy(l => l.Kpmp_data_integration_term).ThenBy(l => l.Value_type).ToArray();
            return array;
        }

        public static KPMP_nearSingleCell_proteomics_UCSF_line_class[] Order_by_set_gene_descending_value(KPMP_nearSingleCell_proteomics_UCSF_line_class[] array)
        {
            array = array.OrderBy(l => l.Dataset).ThenBy(l => l.Subsegment).ThenBy(l => l.SampleName).ThenBy(l => l.PatientID).ThenBy(l => l.Kpmp_data_integration_term).ThenBy(l => l.Value_type).ThenBy(l => l.Gene_symbol).ThenByDescending(l => l.Value).ToArray();
            return array;
        }

        public static KPMP_nearSingleCell_proteomics_UCSF_line_class[] Order_by_set_descending_value(KPMP_nearSingleCell_proteomics_UCSF_line_class[] array)
        {
            array = array.OrderBy(l => l.Dataset).ThenBy(l => l.Subsegment).ThenBy(l => l.SampleName).ThenBy(l => l.PatientID).ThenBy(l => l.Kpmp_data_integration_term).ThenBy(l => l.Value_type).ThenByDescending(l => l.Value).ToArray();
            return array;
        }

        public bool Equal_set(KPMP_nearSingleCell_proteomics_UCSF_line_class other)
        {
            bool equal = (this.Dataset.Equals(other.Dataset))
                         && (this.SampleName.Equals(other.SampleName))
                         && (this.Subsegment.Equals(other.Subsegment))
                         && (this.PatientID.Equals(other.PatientID))
                         && (this.Kpmp_data_integration_term.Equals(other.Kpmp_data_integration_term))
                         && (this.Value_type.Equals(other.Value_type));
            return equal;
        }
        #endregion

        public KPMP_nearSingleCell_proteomics_UCSF_line_class()
        {
            Gene_symbol = "";
            Subsegment = "";
            SampleName = "";
            PatientID = "";
            Dataset = "";
            Kpmp_data_integration_term = "";
        }

        public KPMP_nearSingleCell_proteomics_UCSF_line_class Deep_copy()
        {
            KPMP_nearSingleCell_proteomics_UCSF_line_class copy = (KPMP_nearSingleCell_proteomics_UCSF_line_class)this.MemberwiseClone();
            copy.Gene_symbol = (string)this.Gene_symbol.Clone();
            copy.Subsegment = (string)this.Subsegment.Clone();
            copy.SampleName = (string)this.SampleName.Clone();
            copy.PatientID = (string)this.PatientID.Clone();
            copy.Dataset = (string)this.Dataset.Clone();
            copy.Kpmp_data_integration_term = (string)this.Kpmp_data_integration_term.Clone();
            return copy;
        }
    }

    class KPMP_nearSingleCell_proteomics_UCSF_class
    {
        public KPMP_nearSingleCell_proteomics_UCSF_line_class[] Data { get; set; }
        public KPMP_integration_paper_metadata_class Dataset_patient { get; set; }
        public string[] Bg_genes_in_upperCase { get; set; }

        public KPMP_nearSingleCell_proteomics_UCSF_class()
        {
            Data = new KPMP_nearSingleCell_proteomics_UCSF_line_class[0];
            Dataset_patient = new KPMP_integration_paper_metadata_class();
            Bg_genes_in_upperCase = new string[0];
        }

        private void Check_for_duplicates_for_de_instance()
        {
            this.Data = this.Data.OrderBy(l => l.Dataset).ThenBy(l => l.PatientID).ThenBy(l => l.Subsegment).ThenBy(l => l.Kpmp_data_integration_term).ThenBy(l => l.Gene_symbol).ThenBy(l=>l.Value_type).ToArray();
            int data_length = this.Data.Length;
            KPMP_nearSingleCell_proteomics_UCSF_line_class proteomics_line;
            for (int indexD = 0; indexD < data_length; indexD++)
            {
                proteomics_line = this.Data[indexD];
                if ((indexD != 0)
                    && (proteomics_line.Dataset.Equals(this.Data[indexD - 1].Dataset))
                    && (proteomics_line.PatientID.Equals(this.Data[indexD - 1].PatientID))
                    && (proteomics_line.Subsegment.Equals(this.Data[indexD - 1].Subsegment))
                    && (proteomics_line.Kpmp_data_integration_term.Equals(this.Data[indexD - 1].Kpmp_data_integration_term))
                    && (proteomics_line.Gene_symbol.Equals(this.Data[indexD - 1].Gene_symbol))
                    && (proteomics_line.Value_type.Equals(this.Data[indexD - 1].Value_type)))
                {
                    throw new Exception();
                }
            }
        }

        private void Add_to_array(KPMP_nearSingleCell_proteomics_UCSF_line_class[] add_data)
        {
            int this_length = this.Data.Length;
            int add_length = add_data.Length;
            int new_length = this_length + add_length;
            KPMP_nearSingleCell_proteomics_UCSF_line_class[] new_data = new KPMP_nearSingleCell_proteomics_UCSF_line_class[new_length];
            int indexNew = -1;
            for (int indexThis = 0; indexThis < this_length; indexThis++)
            {
                indexNew++;
                new_data[indexNew] = this.Data[indexThis];
            }
            for (int indexAdd = 0; indexAdd < add_length; indexAdd++)
            {
                indexNew++;
                new_data[indexNew] = add_data[indexAdd];
            }
            this.Data = new_data;
        }

        private void Set_segment_patient_id()
        {
            int data_length = Data.Length;
            KPMP_nearSingleCell_proteomics_UCSF_line_class proteomics_line;
            int indexDot;
            int indexFirstMinus;
            int indexLastMinus;
            List<string> non_single_value_sampleNames = new List<string>();
            List<string> single_value_sampleNames = new List<string>();
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                proteomics_line = this.Data[indexData];
                indexDot = 0;
                indexFirstMinus = proteomics_line.SampleName.IndexOf('-');
                indexLastMinus = proteomics_line.SampleName.Length;
                if ((indexDot != -1) && (indexFirstMinus != -1))
                {
                    if (!proteomics_line.Value_type.Equals(KPMP_value_type_enum.Single_value)) { throw new Exception(); }
                    single_value_sampleNames.Add(proteomics_line.SampleName);
                    proteomics_line.PatientID = proteomics_line.SampleName.Substring(indexFirstMinus + 1, indexLastMinus - indexFirstMinus - 1);
                    proteomics_line.Subsegment = proteomics_line.SampleName.Substring(indexDot, indexFirstMinus);
                }
                else if ((indexDot == 0) && (indexFirstMinus == -1))
                {
                    if (proteomics_line.Value_type.Equals(KPMP_value_type_enum.Single_value)) { throw new Exception(); }
                    non_single_value_sampleNames.Add(proteomics_line.SampleName);
                    proteomics_line.PatientID = KPMP_data_integration_class.Get_combined_patients_label();
                    proteomics_line.Subsegment = (string)proteomics_line.SampleName.Clone();
                }
                else
                {
                    throw new Exception();
                }
            }
            non_single_value_sampleNames = non_single_value_sampleNames.Distinct().OrderBy(l => l).ToList();
            single_value_sampleNames = single_value_sampleNames.Distinct().OrderBy(l => l).ToList();
        }

        private void Set_segment_patient_id_old()
        {
            int data_length = Data.Length;
            KPMP_nearSingleCell_proteomics_UCSF_line_class proteomics_line;
            int indexDot;
            int indexFirstMinus;
            int indexLastMinus;
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                proteomics_line = this.Data[indexData];
                indexDot = proteomics_line.SampleName.IndexOf('.');
                indexFirstMinus = proteomics_line.SampleName.IndexOf('-');
                indexLastMinus = proteomics_line.SampleName.LastIndexOf('-');
                if ((indexDot != -1) && (indexFirstMinus != -1))
                {
                    proteomics_line.PatientID = proteomics_line.SampleName.Substring(indexFirstMinus + 1, indexLastMinus - indexFirstMinus - 1);
                    proteomics_line.Subsegment = proteomics_line.SampleName.Substring(indexDot + 1, 1);
                }
                else if ((indexDot == -1) && (indexFirstMinus == -1))
                {
                    proteomics_line.PatientID = KPMP_data_integration_class.Get_combined_patients_label();
                    proteomics_line.Subsegment = (string)proteomics_line.SampleName.Clone();
                }
                else
                {
                    throw new Exception();
                }
            }
        }

        private void Add_kpmp_integration_term()
        {
            int data_length = Data.Length;
            KPMP_nearSingleCell_proteomics_UCSF_line_class proteomics_line;
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                proteomics_line = this.Data[indexData];
                proteomics_line.Kpmp_data_integration_term = KPMP_data_integration_class.Get_kpmp_integration_term_plus_patientID(proteomics_line.Subsegment, proteomics_line.PatientID);
            }
        }

        private void Remove_all_lines_with_empty_symbols()
        {
            List<KPMP_nearSingleCell_proteomics_UCSF_line_class> keep = new List<KPMP_nearSingleCell_proteomics_UCSF_line_class>();
            int data_length = this.Data.Length;
            KPMP_nearSingleCell_proteomics_UCSF_line_class line;
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                line = this.Data[indexData];
                if (!String.IsNullOrEmpty(line.Gene_symbol))
                {
                    keep.Add(line);
                }
            }
            this.Data = keep.ToArray();
        }

        private void Check_if_duplicated_gene_Symbols()
        {
            this.Data = KPMP_nearSingleCell_proteomics_UCSF_line_class.Order_by_set_gene_descending_value(this.Data);
            int data_length = this.Data.Length;
            KPMP_nearSingleCell_proteomics_UCSF_line_class line;
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                line = this.Data[indexData];
                if ((indexData != 0)
                    && (line.Equal_set(Data[indexData - 1]))
                    && (line.Gene_symbol.Equals(Data[indexData - 1].Gene_symbol)))
                {
                    throw new Exception();
                }
            }
        }

        private void Merge_duplicated_gene_Symbols_and_throw_exception_if_duplicated_gene_symbols()
        {
            KPMP_nearSingleCell_proteomics_UCSF_line_class subsegmental_proteomics_line;
            int data_length;

            #region Merge duplicated symbols
            this.Data = this.Data.OrderBy(l => l.Value_type).ThenBy(l => l.Dataset).ThenBy(l => l.Gene_symbol).ThenBy(l => l.SampleName).ToArray();
            List<KPMP_nearSingleCell_proteomics_UCSF_line_class> keep = new List<KPMP_nearSingleCell_proteomics_UCSF_line_class>();
            data_length = this.Data.Length;
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                subsegmental_proteomics_line = this.Data[indexData];
                if (subsegmental_proteomics_line.Value_type.Equals(KPMP_value_type_enum.Single_value))
                {
                    if ((indexData == data_length - 1)
                        || (!subsegmental_proteomics_line.Dataset.Equals(Data[indexData + 1].Dataset))
                        || (!subsegmental_proteomics_line.Gene_symbol.Equals(Data[indexData + 1].Gene_symbol))
                        || (!subsegmental_proteomics_line.SampleName.Equals(Data[indexData + 1].SampleName)))
                    {
                        keep.Add(subsegmental_proteomics_line);
                    }
                    else
                    {
                        this.Data[indexData + 1].Value += subsegmental_proteomics_line.Value;
                        throw new Exception();
                    }
                }
                else
                {
                    keep.Add(subsegmental_proteomics_line);
                }
            }
            this.Data = keep.ToArray();
            if (this.Data.Length != data_length)
            {
                throw new Exception();
            }
            #endregion
        }

        private void Add_individual_patient_ratios()
        {
            double glomerular_value = -1;
            double tubular_value = -1;
            int data_length = this.Data.Length;
            KPMP_nearSingleCell_proteomics_UCSF_line_class proteomics_line;
            KPMP_nearSingleCell_proteomics_UCSF_line_class new_proteomics_line;
            List<KPMP_nearSingleCell_proteomics_UCSF_line_class> new_proteomics_lines = new List<KPMP_nearSingleCell_proteomics_UCSF_line_class>();
            this.Data = this.Data.OrderBy(l => l.Dataset).ThenBy(l => l.PatientID).ThenBy(l => l.Gene_symbol).ToArray();
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                proteomics_line = this.Data[indexData];
                if (proteomics_line.Value_type.Equals(KPMP_value_type_enum.Single_value))
                {
                    if ((indexData == 0)
                        || (!proteomics_line.PatientID.Equals(this.Data[indexData - 1].PatientID))
                        || (!proteomics_line.Gene_symbol.Equals(this.Data[indexData - 1].Gene_symbol)))
                    {
                        glomerular_value = -1;
                        tubular_value = -1;
                    }
                    if (proteomics_line.Subsegment.Equals("G"))
                    {
                        if (glomerular_value != -1) { throw new Exception(); }
                        glomerular_value = proteomics_line.Value;
                    }
                    else if (proteomics_line.Subsegment.Equals("T"))
                    {
                        if (tubular_value != -1) { throw new Exception(); }
                        tubular_value = proteomics_line.Value;
                    }
                    if ((indexData == data_length - 1)
                        || (!proteomics_line.PatientID.Equals(this.Data[indexData + 1].PatientID))
                        || (!proteomics_line.Gene_symbol.Equals(this.Data[indexData + 1].Gene_symbol)))
                    {
                        if ((tubular_value == -1) || (glomerular_value == -1)) { throw new Exception(); }
                        new_proteomics_line = proteomics_line.Deep_copy();
                        new_proteomics_line.Subsegment = "Ratio GvsTI";
                        new_proteomics_line.Value_type = KPMP_value_type_enum.Ratioavg;
                        new_proteomics_line.Value = (glomerular_value + 1) / (tubular_value + 1);
                        new_proteomics_lines.Add(new_proteomics_line);

                        new_proteomics_line = proteomics_line.Deep_copy();
                        new_proteomics_line.Subsegment = "Ratio TIvsG";
                        new_proteomics_line.Value_type = KPMP_value_type_enum.Ratioavg;
                        new_proteomics_line.Value = (tubular_value + 1) / (glomerular_value + 1);
                        new_proteomics_lines.Add(new_proteomics_line);
                    }
                }
            }
            this.Add_to_array(new_proteomics_lines.ToArray());
        }

        private void Add_reverted_ratios_and_adjust_subsegment_names()
        {
            int data_length = this.Data.Length;
            KPMP_nearSingleCell_proteomics_UCSF_line_class subsegmental_line;
            KPMP_nearSingleCell_proteomics_UCSF_line_class new_subsegmental_line;
            List<KPMP_nearSingleCell_proteomics_UCSF_line_class> new_lines = new List<KPMP_nearSingleCell_proteomics_UCSF_line_class>();
            for (int indexD = 0; indexD < data_length; indexD++)
            {
                subsegmental_line = this.Data[indexD];
                if (subsegmental_line.SampleName.Equals("G/T"))
                {
                    new_subsegmental_line = subsegmental_line.Deep_copy();
                    new_subsegmental_line.SampleName = "T/G";
                    new_subsegmental_line.Subsegment = "Ratio TIvsG";
                    new_subsegmental_line.Value = 1F / subsegmental_line.Value;
                    new_lines.Add(new_subsegmental_line);

                    subsegmental_line.Subsegment = "Ratio GvsTI";
                }
            }
            Add_to_array(new_lines.ToArray());
        }

        private void Specify_subsegment_for_pvalues_based_on_fold_change()
        {
            int data_length = Data.Length;
            KPMP_nearSingleCell_proteomics_UCSF_line_class data_line;
            KPMP_nearSingleCell_proteomics_UCSF_line_class inner_data_line;
            this.Data = this.Data.OrderBy(l => l.Dataset).ThenBy(l => l.Gene_symbol).ThenBy(l => l.Subsegment).ToArray();
            int firstIndex_same_description = -1;
            double fold_change_glom_over_ti = -1;
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                data_line = this.Data[indexData];
                if ((indexData == 0)
                    || (!data_line.Dataset.Equals(this.Data[indexData - 1].Dataset))
                    || (!data_line.Gene_symbol.Equals(this.Data[indexData - 1].Gene_symbol)))
                {
                    firstIndex_same_description = indexData;
                    fold_change_glom_over_ti = -1;
                }
                if ((data_line.Value_type.Equals(KPMP_value_type_enum.Ratioavg))
                    && (data_line.Subsegment.Equals("Fold Glom/PT")))
                {
                    if (fold_change_glom_over_ti != -1) { throw new Exception(); }
                    fold_change_glom_over_ti = data_line.Value;
                }
                if ((indexData == data_length - 1)
                    || (!data_line.Dataset.Equals(this.Data[indexData + 1].Dataset))
                    || (!data_line.Gene_symbol.Equals(this.Data[indexData + 1].Gene_symbol)))
                {
                    if (fold_change_glom_over_ti == -1) { throw new Exception(); }
                    for (int indexInner = firstIndex_same_description; indexInner <= indexData; indexInner++)
                    {
                        inner_data_line = this.Data[indexInner];
                        if (  (inner_data_line.Value_type.Equals(KPMP_value_type_enum.Minus_log10_pvalue))
                            || (inner_data_line.Value_type.Equals(KPMP_value_type_enum.Log2_ratioavg))
                            || (inner_data_line.Value_type.Equals(KPMP_value_type_enum.Ratioavg)))
                        {
                            if (fold_change_glom_over_ti > 1)
                            {
                                inner_data_line.Subsegment = "Glom";
                            }
                            else
                            {
                                inner_data_line.Subsegment = "PT";
                            }
                        }
                    }
                }
            }
        }

        private void Reverse_foldChanges_for_pt()
        {
            int data_length = Data.Length;
            KPMP_nearSingleCell_proteomics_UCSF_line_class data_line;
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                data_line = this.Data[indexData];
                if (data_line.Subsegment.Equals("PT"))
                {
                    switch (data_line.Value_type)
                    {
                        case KPMP_value_type_enum.Log2_ratioavg:
                            data_line.Value = -data_line.Value;
                            break;
                        case KPMP_value_type_enum.Ratioavg:
                            data_line.Value = (double)1 / data_line.Value;
                            break;
                        case KPMP_value_type_enum.Average:
                        case KPMP_value_type_enum.Median:
                        case KPMP_value_type_enum.Minus_log10_fdr:
                        case KPMP_value_type_enum.Minus_log10_pvalue:
                        case KPMP_value_type_enum.Single_value:
                            break;
                        default:
                            throw new Exception();
                    }
                }
                else if (!data_line.Subsegment.Equals("Glom"))
                {
                    throw new Exception();
                }
            }
        }


        private void Duplicate_lines_and_add_cell_specificity()
        {
            int data_length = Data.Length;
            KPMP_nearSingleCell_proteomics_UCSF_line_class data_line;
            KPMP_nearSingleCell_proteomics_UCSF_line_class new_data_line;
            List<KPMP_nearSingleCell_proteomics_UCSF_line_class> add = new List<KPMP_nearSingleCell_proteomics_UCSF_line_class>();
            for (int indexD = 0; indexD < data_length; indexD++)
            {
                data_line = this.Data[indexD];
                if (data_line.Subsegment.Equals("Glom"))
                {
                    new_data_line = data_line.Deep_copy();
                    new_data_line.Subsegment = new_data_line.Subsegment + " - POD";
                    add.Add(new_data_line);
                    new_data_line = data_line.Deep_copy();
                    new_data_line.Subsegment = new_data_line.Subsegment + " - Mesangial";
                    add.Add(new_data_line);
                    new_data_line = data_line.Deep_copy();
                    new_data_line.Subsegment = new_data_line.Subsegment + " - EPC";
                    add.Add(new_data_line);
                    data_line.Subsegment = data_line.Subsegment + " - EC";
                }
                else if (!data_line.Subsegment.Equals("PT"))
                {
                    throw new Exception();
                }
            }
            Add_to_array(add.ToArray());
        }

        private void Set_bg_proteins_in_upperCase()
        {
            this.Data = this.Data.OrderBy(l => l.Gene_symbol).ToArray();
            List<string> bg_genes_list = new List<string>();
            int data_length = this.Data.Length;
            KPMP_nearSingleCell_proteomics_UCSF_line_class proteomics_line;
            for (int indexD = 0; indexD < data_length; indexD++)
            {
                proteomics_line = this.Data[indexD];
                if ((indexD == 0) || (!proteomics_line.Gene_symbol.Equals(this.Data[indexD - 1].Gene_symbol)))
                {
                    bg_genes_list.Add(proteomics_line.Gene_symbol.ToUpper());
                }
            }
            this.Bg_genes_in_upperCase = bg_genes_list.ToArray();
        }

        private void Set_all_geneSymbols_to_upperCase()
        {
            foreach (KPMP_nearSingleCell_proteomics_UCSF_line_class nearSingleCell_proteomics_line in this.Data)
            {
                nearSingleCell_proteomics_line.Gene_symbol = nearSingleCell_proteomics_line.Gene_symbol.ToUpper();
            }
        }

        public string[] Get_deep_copy_of_bg_proteins_in_upper_case()
        {
            int bg_genes_length = this.Bg_genes_in_upperCase.Length;
            string[] copy = new string[bg_genes_length];
            for (int indexBg = 0; indexBg < bg_genes_length; indexBg++)
            {
                copy[indexBg] = (string)this.Bg_genes_in_upperCase[indexBg].Clone();
            }
            return copy;
        }

        private void Check_if_each_gene_is_only_glom_or_PT_and_that_log2Ratioavg_and_minusLog10Pvalue_are_largerEqual_zero()
        {
            int data_length = this.Data.Length;
            KPMP_nearSingleCell_proteomics_UCSF_line_class proteomics_line;
            this.Data = Data.OrderBy(l => l.Dataset).ThenBy(l => l.Uniprot_accession).ThenBy(l => l.Gene_symbol).ThenBy(l => l.PatientID).ToArray();
            bool current_glom_logratioavg_set = false;
            bool current_ti_logratioavg_set = false;
            bool current_glom_minusLog10Pvalue_set = false;
            bool current_ti_minusLog10Pvalue_set = false;
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                proteomics_line = this.Data[indexData];
                if ((indexData == 0)
                    || (!proteomics_line.Dataset.Equals(this.Data[indexData - 1].Dataset))
                    || (!proteomics_line.Uniprot_accession.Equals(this.Data[indexData - 1].Uniprot_accession))
                    || (!proteomics_line.Gene_symbol.Equals(this.Data[indexData - 1].Gene_symbol))
                    || (!proteomics_line.PatientID.Equals(this.Data[indexData - 1].PatientID)))
                {
                    current_glom_logratioavg_set = false;
                    current_ti_logratioavg_set = false;
                    current_glom_minusLog10Pvalue_set = false;
                    current_ti_minusLog10Pvalue_set = false;
                }
                switch (proteomics_line.Value_type)
                {
                    case KPMP_value_type_enum.Log2_ratioavg:
                    case KPMP_value_type_enum.Minus_log10_pvalue:
                        if (proteomics_line.Value < 0) { throw new Exception(); }
                        break;
                    default:
                        throw new Exception();

                }

                switch (proteomics_line.Subsegment.ToUpper())
                {
                    case "GLOM":
                        switch (proteomics_line.Value_type)
                        {
                            case KPMP_value_type_enum.Log2_ratioavg:
                                if (current_glom_logratioavg_set) { throw new Exception(); }
                                current_glom_logratioavg_set = true;
                                break;
                            case KPMP_value_type_enum.Minus_log10_pvalue:
                                if (current_glom_minusLog10Pvalue_set) { throw new Exception(); }
                                current_glom_minusLog10Pvalue_set = true;
                                break;
                            default:
                                break;
                        }
                        break;
                    case "PT":
                        switch (proteomics_line.Value_type)
                        {
                            case KPMP_value_type_enum.Log2_ratioavg:
                                if (current_ti_logratioavg_set) { throw new Exception(); }
                                current_ti_logratioavg_set = true;
                                break;
                            case KPMP_value_type_enum.Minus_log10_pvalue:
                                if (current_ti_minusLog10Pvalue_set) { throw new Exception(); }
                                current_ti_minusLog10Pvalue_set = true;
                                break;
                            default:
                                break;
                        }
                        break;
                    default:
                        throw new Exception();
                }
                if ((indexData == data_length - 1)
                    || (!proteomics_line.Dataset.Equals(this.Data[indexData + 1].Dataset))
                    || (!proteomics_line.Uniprot_accession.Equals(this.Data[indexData + 1].Uniprot_accession))
                    || (!proteomics_line.Gene_symbol.Equals(this.Data[indexData + 1].Gene_symbol))
                    || (!proteomics_line.PatientID.Equals(this.Data[indexData + 1].PatientID)))
                {
                    if ((current_glom_minusLog10Pvalue_set) && (current_ti_minusLog10Pvalue_set)) { throw new Exception(); }
                    if ((current_glom_logratioavg_set) && (current_ti_logratioavg_set)) { throw new Exception(); }
                    if ((!current_glom_minusLog10Pvalue_set) && (!current_ti_minusLog10Pvalue_set)) { throw new Exception(); }
                    if ((!current_glom_logratioavg_set) && (!current_ti_logratioavg_set)) { throw new Exception(); }
                }
            }
        }


        public void Generate_including_reversal_of_fold_changes_for_pt(string subdirectory, string dataset)
        {
            Read_and_add_minusLog10Pvalue_and_add_log2ratioavg(subdirectory, dataset);
            Set_all_geneSymbols_to_upperCase();
            Merge_duplicated_gene_Symbols_and_throw_exception_if_duplicated_gene_symbols();
            Set_segment_patient_id();
            Generate_dataset_patient_instance();
            Remove_all_lines_with_empty_symbols();
            Set_bg_proteins_in_upperCase();
            Check_if_duplicated_gene_Symbols();
            Specify_subsegment_for_pvalues_based_on_fold_change();
            Reverse_foldChanges_for_pt();
            Keep_only_indicated_value_types(new KPMP_value_type_enum[] { KPMP_value_type_enum.Minus_log10_pvalue, KPMP_value_type_enum.Log2_ratioavg });
            Check_if_each_gene_is_only_glom_or_PT_and_that_log2Ratioavg_and_minusLog10Pvalue_are_largerEqual_zero();
            //Add_individual_patient_ratios();
            //Add_reverted_ratios_and_adjust_subsegment_names();
            Duplicate_lines_and_add_cell_specificity();
            Add_kpmp_integration_term();
            Check_for_duplicates_for_de_instance();
        }

        public void Replace_all_segment0s_by_segment1s(string segment0, string segment1)
        {
            foreach (KPMP_nearSingleCell_proteomics_UCSF_line_class subsegmental_proteomic_line in this.Data)
            {
                if (subsegmental_proteomic_line.Subsegment.Equals(segment1)) { throw new Exception(); }
                if (subsegmental_proteomic_line.Subsegment.Equals(segment0))
                {
                    subsegmental_proteomic_line.Subsegment = (string)segment1.Clone();
                }
            }
        }

        public void Generate_for_single_patient_analysis(string subdirectory, string dataset)
        {
            Read_and_add_minusLog10Pvalue_and_add_log2ratioavg(subdirectory, dataset);
            Set_all_geneSymbols_to_upperCase();
            Merge_duplicated_gene_Symbols_and_throw_exception_if_duplicated_gene_symbols();
            Set_segment_patient_id();
            Generate_dataset_patient_instance();
            Remove_all_lines_with_empty_symbols();
            Set_bg_proteins_in_upperCase();
            Check_if_duplicated_gene_Symbols();
            Keep_only_indicated_value_types(KPMP_value_type_enum.Single_value);
            Check_for_duplicates_for_de_instance();
        }

        public KPMP_standardized_dataset_class Generate_kpmp_standardized_dataset_instance(KPMP_value_type_enum value_type_1st, KPMP_value_type_enum value_type_2nd)
        {
            KPMP_standardized_dataset_line_class standardized_data_line;
            List<KPMP_standardized_dataset_line_class> standardized_data_list = new List<KPMP_standardized_dataset_line_class>();
            int data_length = this.Data.Length;
            KPMP_nearSingleCell_proteomics_UCSF_line_class proteomics_line;
            this.Data = this.Data.OrderBy(l => l.Dataset).ThenBy(l => l.Subsegment).ThenBy(l => l.PatientID).ThenBy(l => l.Kpmp_data_integration_term).ThenBy(l => l.Gene_symbol).ThenBy(l => l.Value_type).ToArray();
            string dataset = (string)this.Data[0].Dataset.Clone();
            double current_value_1st = -1;
            double current_value_2nd = -1;
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                proteomics_line = this.Data[indexData];
                if (!proteomics_line.Dataset.Equals(dataset)) { throw new Exception(); }
                if ((indexData != 0)
                    && (proteomics_line.Dataset.Equals(this.Data[indexData - 1].Dataset))
                    && (proteomics_line.Subsegment.Equals(this.Data[indexData - 1].Subsegment))
                    && (proteomics_line.PatientID.Equals(this.Data[indexData - 1].PatientID))
                    && (proteomics_line.Kpmp_data_integration_term.Equals(this.Data[indexData - 1].Kpmp_data_integration_term))
                    && (proteomics_line.Gene_symbol.Equals(this.Data[indexData - 1].Gene_symbol))
                    && (proteomics_line.Value_type.Equals(this.Data[indexData - 1].Value_type)))
                {
                    throw new Exception();
                }
                if ((indexData == 0)
                    || (!proteomics_line.Dataset.Equals(this.Data[indexData - 1].Dataset))
                    || (!proteomics_line.Subsegment.Equals(this.Data[indexData - 1].Subsegment))
                    || (!proteomics_line.PatientID.Equals(this.Data[indexData - 1].PatientID))
                    || (!proteomics_line.Kpmp_data_integration_term.Equals(this.Data[indexData - 1].Kpmp_data_integration_term))
                    || (!proteomics_line.Gene_symbol.Equals(this.Data[indexData - 1].Gene_symbol)))
                {
                    current_value_1st = -1;
                    switch (value_type_2nd)
                    {
                        case KPMP_value_type_enum.No_selection:
                            current_value_2nd = 0;
                            break;
                        default:
                            current_value_2nd = -1;
                            break;
                    }
                }
                if (proteomics_line.Value_type.Equals(value_type_1st))
                {
                    if (current_value_1st != -1) { throw new Exception(); }
                    current_value_1st = proteomics_line.Value;
                }

                if (proteomics_line.Value_type.Equals(value_type_2nd))
                {
                    if (current_value_2nd != -1) { throw new Exception(); }
                    current_value_2nd = proteomics_line.Value;
                }
                if ((indexData == data_length - 1)
                    || (!proteomics_line.Dataset.Equals(this.Data[indexData + 1].Dataset))
                    || (!proteomics_line.Subsegment.Equals(this.Data[indexData + 1].Subsegment))
                    || (!proteomics_line.PatientID.Equals(this.Data[indexData + 1].PatientID))
                    || (!proteomics_line.Kpmp_data_integration_term.Equals(this.Data[indexData + 1].Kpmp_data_integration_term))
                    || (!proteomics_line.Gene_symbol.Equals(this.Data[indexData + 1].Gene_symbol)))
                {
                    if (current_value_1st == -1) { throw new Exception(); }
                    if (current_value_2nd == -1) { throw new Exception(); }
                    standardized_data_line = new KPMP_standardized_dataset_line_class();
                    standardized_data_line.Cell_segment = (string)proteomics_line.Subsegment.Clone();
                    standardized_data_line.Dataset = (string)proteomics_line.Dataset.Clone();
                    standardized_data_line.Value_1st = current_value_1st;
                    standardized_data_line.Value_2nd = current_value_2nd;
                    standardized_data_line.KPMP_data_integration_term = (string)proteomics_line.Kpmp_data_integration_term.Clone();
                    standardized_data_line.Value_type_1st = value_type_1st;
                    standardized_data_line.Value_type_2nd = value_type_2nd;
                    standardized_data_line.PatientId = (string)proteomics_line.PatientID.Clone();
                    standardized_data_line.Gene_symbol = (string)proteomics_line.Gene_symbol.Clone();
                    standardized_data_list.Add(standardized_data_line);
                }
            }
            Dictionary<string, string[]> dataset_bgProteins_dict = new Dictionary<string, string[]>();
            dataset_bgProteins_dict.Add(dataset, Get_deep_copy_of_bg_proteins_in_upper_case());
            KPMP_standardized_dataset_class standardized_data = new KPMP_standardized_dataset_class();
            standardized_data.Add_to_existing_instances(standardized_data_list.ToArray(), dataset_bgProteins_dict);
            return standardized_data;

            #region old commented script
            //Report_class.Write("{0}: Generate_kpmp_standardized_dataset_instance", typeof(KPMP_nearSingleCell_proteomics_class).Name);
            //KPMP_standardized_dataset_line_class standardized_data_line;
            //List<KPMP_standardized_dataset_line_class> standardized_data_list = new List<KPMP_standardized_dataset_line_class>();
            //int data_length = this.Data.Length;
            //KPMP_nearSingleCell_proteomics_UCSF_line_class proteomics_line;
            //List<KPMP_value_type_enum> value_types = new List<KPMP_value_type_enum>();
            //this.Data = this.Data.OrderBy(l => l.Value_type).ThenBy(l => l.Dataset).ThenBy(l => l.Subsegment).ThenBy(l => l.PatientID).ThenBy(l => l.Kpmp_data_integration_term).ThenBy(l => l.Gene_symbol).ToArray();
            //string dataset = (string)this.Data[0].Dataset.Clone();
            //for (int indexData = 0; indexData < data_length; indexData++)
            //{
            //    proteomics_line = this.Data[indexData];
            //    if ((indexData == 0) || (!proteomics_line.Value_type.Equals(this.Data[indexData - 1].Value_type)))
            //    {
            //        value_types.Add(proteomics_line.Value_type);
            //    }
            //    if (!proteomics_line.Dataset.Equals(dataset)) { throw new Exception(); }
            //    if ((indexData != 0)
            //        && (!proteomics_line.Dataset.Equals(this.Data[indexData - 1].Dataset))
            //        && (!proteomics_line.Subsegment.Equals(this.Data[indexData - 1].Subsegment))
            //        && (!proteomics_line.PatientID.Equals(this.Data[indexData - 1].PatientID))
            //        && (!proteomics_line.Kpmp_data_integration_term.Equals(this.Data[indexData - 1].Kpmp_data_integration_term))
            //        && (!proteomics_line.Gene_symbol.Equals(this.Data[indexData - 1].Gene_symbol))
            //        && (!proteomics_line.Value_type.Equals(this.Data[indexData - 1].Value_type)))
            //    {
            //        throw new Exception();
            //    }
            //    standardized_data_line = new KPMP_standardized_dataset_line_class();
            //    standardized_data_line.Cell_segment = (string)proteomics_line.Subsegment.Clone();
            //    standardized_data_line.Dataset = (string)proteomics_line.Dataset.Clone();
            //    standardized_data_line.Value = proteomics_line.Value;
            //    standardized_data_line.KPMP_data_integration_term = (string)proteomics_line.Kpmp_data_integration_term.Clone();
            //    standardized_data_line.Value_type = proteomics_line.Value_type;
            //    standardized_data_line.PatientId = (string)proteomics_line.PatientID.Clone();
            //    standardized_data_line.Gene_symbol = (string)proteomics_line.Gene_symbol.Clone();
            //    standardized_data_list.Add(standardized_data_line);
            //}
            //Dictionary<string, string[]> dataset_bgGenesProteins_dict = new Dictionary<string, string[]>();
            //dataset_bgGenesProteins_dict.Add(dataset, Get_deep_copy_of_bg_proteins_in_upper_case());

            //KPMP_standardized_dataset_class standardized_data = new KPMP_standardized_dataset_class();
            //standardized_data.Add_to_existing_instances(standardized_data_list.ToArray(), dataset_bgGenesProteins_dict);
            //foreach (KPMP_value_type_enum value_type in value_types)
            //{
            //    Report_class.Write(" {0}", value_type);
            //}
            //Report_class.WriteLine();
            //return standardized_data;
            #endregion
        }

        public KPMP_standardized_dataset_class Generate_standardized_dataset_instance_with_all_values_filled_in_value_1st()
        {
            KPMP_standardized_dataset_line_class new_standardized_data_line;
            List<KPMP_standardized_dataset_line_class> standardized_data_list = new List<KPMP_standardized_dataset_line_class>();
            int data_length = this.Data.Length;
            KPMP_nearSingleCell_proteomics_UCSF_line_class proteomics_line;
            this.Data = this.Data.OrderBy(l => l.Dataset).ThenBy(l => l.Subsegment).ThenBy(l => l.PatientID).ThenBy(l => l.Kpmp_data_integration_term).ThenBy(l => l.Value_type).ThenBy(l => l.Gene_symbol).ToArray();
            for (int indexR = 0; indexR < data_length; indexR++)
            {
                proteomics_line = this.Data[indexR];
                if ((indexR != 0)
                    && (proteomics_line.Dataset.Equals(this.Data[indexR - 1].Dataset))
                    && (proteomics_line.Subsegment.Equals(this.Data[indexR - 1].Subsegment))
                    && (proteomics_line.PatientID.Equals(this.Data[indexR - 1].PatientID))
                    && (proteomics_line.Kpmp_data_integration_term.Equals(this.Data[indexR - 1].Kpmp_data_integration_term))
                    && (proteomics_line.Value_type.Equals(this.Data[indexR - 1].Value_type))
                    && (proteomics_line.Gene_symbol.Equals(this.Data[indexR - 1].Gene_symbol)))
                {
                    throw new Exception();
                }
                new_standardized_data_line = new KPMP_standardized_dataset_line_class();
                new_standardized_data_line.Cell_segment = (string)proteomics_line.Subsegment.Clone();
                new_standardized_data_line.Dataset = (string)proteomics_line.Dataset.Clone();
                new_standardized_data_line.Gene_symbol = (string)proteomics_line.Gene_symbol.Clone();
                new_standardized_data_line.PatientId = (string)proteomics_line.PatientID.Clone();
                new_standardized_data_line.Value_1st = proteomics_line.Value;
                new_standardized_data_line.Value_2nd = 0;
                new_standardized_data_line.Value_type_1st = proteomics_line.Value_type;
                new_standardized_data_line.Value_type_2nd = KPMP_value_type_enum.No_selection;
                standardized_data_list.Add(new_standardized_data_line);
            }
            string dataset = KPMP_dataset_name_class.NearSingleCell_proteomics_ucsf;
            Dictionary<string, string[]> dataset_bgSymbol_dict = new Dictionary<string, string[]>();
            dataset_bgSymbol_dict.Add(dataset, Get_deep_copy_of_bg_proteins_in_upper_case());
            KPMP_standardized_dataset_class standard = new KPMP_standardized_dataset_class();
            standard.Add_to_existing_instances(standardized_data_list.ToArray(), dataset_bgSymbol_dict);
            return standard;
        }


        public DE_class Generate_de_instance()
        {
            Fill_de_line_class fill_de_line;
            List<Fill_de_line_class> fill_de_list = new List<Fill_de_line_class>();
            int data_length = this.Data.Length;
            KPMP_nearSingleCell_proteomics_UCSF_line_class subsegmental_proteomics_line;
            for (int indexD = 0; indexD < data_length; indexD++)
            {
                subsegmental_proteomics_line = this.Data[indexD];
                fill_de_line = new Fill_de_line_class();
                fill_de_line.Names_for_de = KPMP_data_integration_class.Get_names_for_de_in_correct_order(subsegmental_proteomics_line.Subsegment, subsegmental_proteomics_line.PatientID, subsegmental_proteomics_line.Dataset, subsegmental_proteomics_line.Value_type, subsegmental_proteomics_line.Kpmp_data_integration_term);
                fill_de_line.Symbols_for_de = new string[] { (string)subsegmental_proteomics_line.Gene_symbol.Clone() };
                fill_de_line.Value_for_de = subsegmental_proteomics_line.Value;
                fill_de_list.Add(fill_de_line);
            }
            DE_class de = new DE_class();
            de.Fill_with_data_alternatively(fill_de_list.ToArray());
            return de;
        }

        public KPMP_integration_paper_metadata_class Get_deep_copy_of_dataset_patient_instance()
        {
            return Dataset_patient.Deep_copy();
        }

        public void Generate_dataset_patient_instance()
        {
            Dictionary<string, Dictionary<string, bool>> dataset_patient_considered_dict = new Dictionary<string, Dictionary<string, bool>>();
            int data_length = this.Data.Length;
            KPMP_integration_paper_metadata_line_class new_dataset_patient_line;
            List<KPMP_integration_paper_metadata_line_class> new_dataset_patient_lines = new List<KPMP_integration_paper_metadata_line_class>();
            KPMP_nearSingleCell_proteomics_UCSF_line_class singleClister_line;
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                singleClister_line = this.Data[indexData];
                if (!dataset_patient_considered_dict.ContainsKey(singleClister_line.Dataset))
                {
                    dataset_patient_considered_dict.Add(singleClister_line.Dataset, new Dictionary<string, bool>());
                }
                if (!dataset_patient_considered_dict[singleClister_line.Dataset].ContainsKey(singleClister_line.PatientID))
                {
                    dataset_patient_considered_dict[singleClister_line.Dataset].Add(singleClister_line.PatientID, true);
                    new_dataset_patient_line = new KPMP_integration_paper_metadata_line_class();
                    new_dataset_patient_line.Dataset = (string)singleClister_line.Dataset.Clone();
                    new_dataset_patient_line.Libraries = new string[] { (string)singleClister_line.PatientID.Clone() };
                    new_dataset_patient_lines.Add(new_dataset_patient_line);
                }
            }
            if (Dataset_patient.Documentations.Length!=0) { throw new Exception(); }
            Dataset_patient = new KPMP_integration_paper_metadata_class();
            Dataset_patient.Add_to_array(new_dataset_patient_lines.ToArray());
        }


        public void Keep_only_indicated_value_types(params KPMP_value_type_enum[] keep_value_types)
        {
            keep_value_types = keep_value_types.Distinct().OrderBy(l => l).ToArray();
            KPMP_value_type_enum keep_value_type;
            int keep_value_types_length = keep_value_types.Length;
            int indexKeep = 0;
            int enumCompare = -2;
            int data_length = this.Data.Length;
            this.Data = this.Data.OrderBy(l => l.Value_type).ToArray();
            KPMP_nearSingleCell_proteomics_UCSF_line_class proteomics_line;
            List<KPMP_nearSingleCell_proteomics_UCSF_line_class> keep = new List<KPMP_nearSingleCell_proteomics_UCSF_line_class>();
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                proteomics_line = this.Data[indexData];
                enumCompare = -2;
                while ((indexKeep < keep_value_types_length) && (enumCompare < 0))
                {
                    keep_value_type = keep_value_types[indexKeep];
                    enumCompare = keep_value_type.CompareTo(proteomics_line.Value_type);
                    if (enumCompare < 0)
                    {
                        indexKeep++;
                    }
                    else if (enumCompare == 0)
                    {
                        keep.Add(proteomics_line);
                    }
                }
            }
            this.Data = keep.ToArray();
        }

        public void Keep_only_indicated_patientIDs(params string[] keep_patientIDs)
        {
            keep_patientIDs = keep_patientIDs.Distinct().OrderBy(l => l).ToArray();
            string keep_patientID;
            int keep_patientIDs_length = keep_patientIDs.Length;
            int indexKeep = 0;
            int enumCompare = -2;
            int data_length = this.Data.Length;
            this.Data = this.Data.OrderBy(l => l.PatientID).ToArray();
            KPMP_nearSingleCell_proteomics_UCSF_line_class proteomics_line;
            List<KPMP_nearSingleCell_proteomics_UCSF_line_class> keep = new List<KPMP_nearSingleCell_proteomics_UCSF_line_class>();
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                proteomics_line = this.Data[indexData];
                enumCompare = -2;
                while ((indexKeep < keep_patientIDs_length) && (enumCompare < 0))
                {
                    keep_patientID = keep_patientIDs[indexKeep];
                    enumCompare = keep_patientID.CompareTo(proteomics_line.PatientID);
                    if (enumCompare < 0)
                    {
                        indexKeep++;
                    }
                    else if (enumCompare == 0)
                    {
                        keep.Add(proteomics_line);
                    }
                }
            }
            this.Data = keep.ToArray();
        }

        public void Keep_only_lines_with_minimum_value(float minimum_value)
        {
            int data_length = this.Data.Length;
            this.Data = this.Data.OrderBy(l => l.PatientID).ToArray();
            KPMP_nearSingleCell_proteomics_UCSF_line_class proteomics_line;
            List<KPMP_nearSingleCell_proteomics_UCSF_line_class> keep = new List<KPMP_nearSingleCell_proteomics_UCSF_line_class>();
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                proteomics_line = this.Data[indexData];
                if (proteomics_line.Value >= minimum_value)
                {
                    keep.Add(proteomics_line);
                }
            }
            this.Data = keep.ToArray();
        }
        public void Keep_only_lines_with_minimum_selectedValueType_and_all_other_values_of_the_same_lines(float minimum_value, KPMP_value_type_enum selected_value_type)
        {
            this.Data = this.Data.OrderBy(l => l.Uniprot_accession).ThenBy(l => l.Protein_name).ThenBy(l => l.PatientID).ThenBy(l => l.Subsegment).ToArray();
            int data_length = this.Data.Length;
            KPMP_nearSingleCell_proteomics_UCSF_line_class proteomics_line;
            List<KPMP_nearSingleCell_proteomics_UCSF_line_class> keep = new List<KPMP_nearSingleCell_proteomics_UCSF_line_class>();
            List<KPMP_nearSingleCell_proteomics_UCSF_line_class> remove = new List<KPMP_nearSingleCell_proteomics_UCSF_line_class>();
            bool keep_line_block = false;
            int firstIndexData = -1;
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                proteomics_line = this.Data[indexData];
                if ((indexData == 0)
                    || (!proteomics_line.Uniprot_accession.Equals(this.Data[indexData - 1].Uniprot_accession))
                    || (!proteomics_line.Protein_name.Equals(this.Data[indexData - 1].Protein_name))
                    || (!proteomics_line.PatientID.Equals(this.Data[indexData - 1].PatientID))
                    || (!proteomics_line.Subsegment.Equals(this.Data[indexData - 1].Subsegment)))
                {
                    keep_line_block = false;
                    firstIndexData = indexData;
                }
                if (proteomics_line.Value_type.Equals(selected_value_type))
                {
                    if (proteomics_line.Value >= minimum_value)
                    {
                        keep_line_block = true;
                    }
                }
                if ((indexData == data_length - 1)
                    || (!proteomics_line.Uniprot_accession.Equals(this.Data[indexData + 1].Uniprot_accession))
                    || (!proteomics_line.Protein_name.Equals(this.Data[indexData + 1].Protein_name))
                    || (!proteomics_line.PatientID.Equals(this.Data[indexData + 1].PatientID))
                    || (!proteomics_line.Subsegment.Equals(this.Data[indexData + 1].Subsegment)))
                {
                    if (keep_line_block)
                    {
                        for (int indexInner = firstIndexData; indexInner <= indexData; indexInner++)
                        {
                            keep.Add(this.Data[indexInner]);
                        }
                    }
                    else
                    {
                        for (int indexInner = firstIndexData; indexInner <= indexData; indexInner++)
                        {
                            remove.Add(this.Data[indexInner]);
                        }
                    }
                }
            }
            this.Data = keep.ToArray();
        }
        public void Keep_top_x_lines_of_each_set_based_on_indicated_descending_valueType_and_keep_all_other_valueType_of_same_gene(int keep_top_x, KPMP_value_type_enum selected_valueType)
        {
            int data_length = this.Data.Length;
            KPMP_nearSingleCell_proteomics_UCSF_line_class proteomics_line;
            KPMP_nearSingleCell_proteomics_UCSF_line_class inner_proteomics_line;
            List<KPMP_nearSingleCell_proteomics_UCSF_line_class> keep = new List<KPMP_nearSingleCell_proteomics_UCSF_line_class>();
            this.Data = this.Data.OrderBy(l => l.Dataset).ThenBy(l => l.Subsegment).ThenBy(l => l.PatientID).ThenBy(l => l.Kpmp_data_integration_term).ThenBy(l => l.Value_type).ThenByDescending(l => l.Value).ToArray();
            int kept_top_x = 0;
            int firstIndex_sameSet = -1;
            Dictionary<string, bool> currentSet_keepGenes_dict = new Dictionary<string, bool>();
            for (int indexD = 0; indexD < data_length; indexD++)
            {
                proteomics_line = this.Data[indexD];
                if ((indexD == 0)
                    || (!proteomics_line.Dataset.Equals(this.Data[indexD - 1].Dataset))
                    || (!proteomics_line.Subsegment.Equals(this.Data[indexD - 1].Subsegment))
                    || (!proteomics_line.PatientID.Equals(this.Data[indexD - 1].PatientID))
                    || (!proteomics_line.Kpmp_data_integration_term.Equals(this.Data[indexD - 1].Kpmp_data_integration_term)))
                {
                    firstIndex_sameSet = indexD;
                    kept_top_x = 0;
                    currentSet_keepGenes_dict.Clear();
                }
                if ((proteomics_line.Value_type.Equals(selected_valueType)) && (kept_top_x < keep_top_x))
                {
                    kept_top_x++;
                    currentSet_keepGenes_dict.Add(proteomics_line.Gene_symbol, true);
                }
                if ((indexD == data_length - 1)
                    || (!proteomics_line.Dataset.Equals(this.Data[indexD + 1].Dataset))
                    || (!proteomics_line.Subsegment.Equals(this.Data[indexD + 1].Subsegment))
                    || (!proteomics_line.PatientID.Equals(this.Data[indexD + 1].PatientID))
                    || (!proteomics_line.Kpmp_data_integration_term.Equals(this.Data[indexD + 1].Kpmp_data_integration_term)))
                {
                    for (int indexInner = firstIndex_sameSet; indexInner <= indexD; indexInner++)
                    {
                        inner_proteomics_line = this.Data[indexInner];
                        if (currentSet_keepGenes_dict.ContainsKey(inner_proteomics_line.Gene_symbol))
                        {
                            keep.Add(inner_proteomics_line);
                        }
                    }
                }
            }
            this.Data = keep.ToArray();
        }
        private void Read_and_add_minusLog10Pvalue_and_add_log2ratioavg(string subdirectory, string dataset)
        {
            string directory = Global_directory_class.Experimental_data_directory + subdirectory;
            string[] complete_fileNames = Directory.GetFiles(directory);
            string complete_fileName;
            int complete_fileNames_length = complete_fileNames.Length;
            KPMP_nearSingleCell_proteomics_UCSF_line_class new_proteomics_line;
            List<KPMP_nearSingleCell_proteomics_UCSF_line_class> proteomics_list = new List<KPMP_nearSingleCell_proteomics_UCSF_line_class>();
            for (int indexC = 0; indexC < complete_fileNames_length; indexC++)
            {
                complete_fileName = complete_fileNames[indexC];
                StreamReader reader = new StreamReader(complete_fileName);
                char delimiter = Global_class.Tab;

                string headline = reader.ReadLine();
                string[] columnNames = headline.Split(delimiter);
                string columnName;
                int columnNames_length = columnNames.Length;

                int indexGeneSymbol = -1;
                int indexRatio_glomPT = -1;
                int indexAverageGlom = -1;
                int indexAveragePT = -1;
                int indexPvalue = -1;
                int indexSwissprot = -1;
                int indexUniprot_accession = -1;
                int indexProtein_name = -1;
                List<int> indexesIntensities_list = new List<int>();

                for (int indexCol = 0; indexCol < columnNames_length; indexCol++)
                {
                    columnName = columnNames[indexCol];
                    if (columnName.Equals("Average Glom"))
                    {
                        if (indexAverageGlom != -1) { throw new Exception(); }
                        indexAverageGlom = indexCol;
                    }
                    else if (columnName.Equals("Average PT"))
                    {
                        if (indexAveragePT != -1) { throw new Exception(); }
                        indexAveragePT = indexCol;
                    }
                    else if (columnName.Equals("Fold Glom/PT"))
                    {
                        if (indexRatio_glomPT != -1) { throw new Exception(); }
                        indexRatio_glomPT = indexCol;
                    }
                    else if (columnName.IndexOf("Intensity") != -1)
                    {
                        indexesIntensities_list.Add(indexCol);
                    }
                    else if (columnName.IndexOf("P_value") != -1)
                    {
                        if (indexPvalue != -1) { throw new Exception(); }
                        indexPvalue = indexCol;
                    }
                    else if (columnName.Equals("Gene Symbol"))
                    {
                        indexGeneSymbol = indexCol;
                    }
                    else if (columnName.IndexOf("Glom") == 0)
                    {
                        indexesIntensities_list.Add(indexCol);
                    }
                    else if (columnName.IndexOf("PT") == 0)
                    {
                        indexesIntensities_list.Add(indexCol);
                    }
                    else if (columnName.IndexOf("SWISSPROT") == 0)
                    {
                        if (indexSwissprot != -1) { throw new Exception(); }
                        indexSwissprot = indexCol;
                    }
                    else if (columnName.IndexOf("UNIPROT ACCESSION") == 0)
                    {
                        if (indexUniprot_accession != -1) { throw new Exception(); }
                        indexUniprot_accession = indexCol;
                    }
                    else if (columnName.IndexOf("Protein Name") == 0)
                    {
                        indexProtein_name = indexCol;
                    }
                    else
                    {
                        throw new Exception();
                    }
                }

                int[] indexesIntensity = indexesIntensities_list.ToArray();
                int indexIntensity;
                int indexesIntensity_length = indexesIntensity.Length;
                string inputLine;
                string[] columnEntries;
                string columnEntry;
                int columnEntries_length;
                string uniprot_accession;
                string swissprot_id;
                string gene_symbol;
                string protein_name;

                while ((inputLine = reader.ReadLine()) != null)
                {
                    columnEntries = inputLine.Split(delimiter);
                    columnEntries_length = columnEntries.Length;
                    if (columnEntries_length!=columnNames_length) { throw new Exception(); }
                    gene_symbol = columnEntries[indexGeneSymbol];
                    uniprot_accession = columnEntries[indexUniprot_accession];
                    swissprot_id = columnEntries[indexSwissprot];
                    protein_name = columnEntries[indexProtein_name];
                    for (int indexIndex = 0; indexIndex < indexesIntensity_length; indexIndex++)
                    {
                        indexIntensity = indexesIntensity[indexIndex];
                        columnName = columnNames[indexIntensity];
                        columnEntry = columnEntries[indexIntensity];

                        if (String.IsNullOrEmpty(columnEntry)) { columnEntry = "0"; }
                        new_proteomics_line = new KPMP_nearSingleCell_proteomics_UCSF_line_class();
                        new_proteomics_line.Dataset = (string)dataset.Clone();
                        new_proteomics_line.Gene_symbol = (string)gene_symbol.Clone();
                        new_proteomics_line.Protein_name = (string)protein_name.Clone();
                        new_proteomics_line.SampleName = (string)columnName.Clone();
                        new_proteomics_line.Swissprot_id = (string)swissprot_id.Clone();
                        new_proteomics_line.Uniprot_accession = (string)uniprot_accession.Clone();
                        new_proteomics_line.Value = double.Parse(columnEntry);
                        new_proteomics_line.Value_type = KPMP_value_type_enum.Single_value;
                        proteomics_list.Add(new_proteomics_line);
                    }
                    //new_proteomics_line = new KPMP_subsegmental_proteomics_UCSF_line_class();
                    //new_proteomics_line.Dataset = (string)dataset.Clone();
                    //new_proteomics_line.Gene_symbol = (string)gene_symbol.Clone();
                    //new_proteomics_line.SampleName = (string)columnNames[indexAverage].Clone();
                    //columnEntry = columnEntries[indexAverage];
                    ////if (!String.IsNullOrEmpty(columnEntry)) { columnEntry = "0"; }
                    //new_proteomics_line.Value = float.Parse(columnEntry);
                    //new_proteomics_line.PatientID = KPMP_data_integration_class.Get_combined_patients_label(); ;
                    //new_proteomics_line.Value_type = KPMP_value_type_enum.Average;
                    //proteomics_list.Add(new_proteomics_line);

                    new_proteomics_line = new KPMP_nearSingleCell_proteomics_UCSF_line_class();
                    new_proteomics_line.Dataset = (string)dataset.Clone();
                    new_proteomics_line.Gene_symbol = (string)gene_symbol.Clone();
                    new_proteomics_line.Protein_name = (string)protein_name.Clone();
                    new_proteomics_line.SampleName = (string)columnNames[indexPvalue].Clone();
                    new_proteomics_line.Swissprot_id = (string)swissprot_id.Clone();
                    new_proteomics_line.Uniprot_accession = (string)uniprot_accession.Clone();
                    columnEntry = columnEntries[indexPvalue];
                    //if (!String.IsNullOrEmpty(columnEntry)) { columnEntry = "0"; }
                    new_proteomics_line.Value = -(double)Math.Log10(double.Parse(columnEntry));
                    new_proteomics_line.PatientID = KPMP_data_integration_class.Get_combined_patients_label(); ;
                    new_proteomics_line.Value_type = KPMP_value_type_enum.Minus_log10_pvalue;
                    proteomics_list.Add(new_proteomics_line);

                    new_proteomics_line = new KPMP_nearSingleCell_proteomics_UCSF_line_class();
                    new_proteomics_line.Dataset = (string)dataset.Clone();
                    new_proteomics_line.Gene_symbol = (string)gene_symbol.Clone();
                    new_proteomics_line.Protein_name = (string)protein_name.Clone();
                    new_proteomics_line.SampleName = (string)columnNames[indexPvalue].Clone();
                    new_proteomics_line.Swissprot_id = (string)swissprot_id.Clone();
                    new_proteomics_line.Uniprot_accession = (string)uniprot_accession.Clone();
                    columnEntry = columnEntries[indexPvalue];
                    //if (!String.IsNullOrEmpty(columnEntry)) { columnEntry = "0"; }
                    new_proteomics_line.Value = double.Parse(columnEntry);
                    new_proteomics_line.PatientID = KPMP_data_integration_class.Get_combined_patients_label(); ;
                    new_proteomics_line.Value_type = KPMP_value_type_enum.Pvalue;
                    //proteomics_list.Add(new_proteomics_line);

                    new_proteomics_line = new KPMP_nearSingleCell_proteomics_UCSF_line_class();
                    new_proteomics_line.Dataset = (string)dataset.Clone();
                    new_proteomics_line.Gene_symbol = (string)gene_symbol.Clone();
                    new_proteomics_line.Protein_name = (string)protein_name.Clone();
                    new_proteomics_line.SampleName = (string)columnNames[indexRatio_glomPT].Clone();
                    new_proteomics_line.Swissprot_id = (string)swissprot_id.Clone();
                    new_proteomics_line.Uniprot_accession = (string)uniprot_accession.Clone();
                    columnEntry = columnEntries[indexRatio_glomPT];
                    //if (!String.IsNullOrEmpty(columnEntry)) { columnEntry = "0"; }
                    new_proteomics_line.Value = double.Parse(columnEntry);
                    new_proteomics_line.PatientID = KPMP_data_integration_class.Get_combined_patients_label(); ;
                    new_proteomics_line.Value_type = KPMP_value_type_enum.Ratioavg;
                    proteomics_list.Add(new_proteomics_line);

                    new_proteomics_line = new KPMP_nearSingleCell_proteomics_UCSF_line_class();
                    new_proteomics_line.Dataset = (string)dataset.Clone();
                    new_proteomics_line.Gene_symbol = (string)gene_symbol.Clone();
                    new_proteomics_line.Protein_name = (string)protein_name.Clone();
                    new_proteomics_line.SampleName = "Log2 ratio avg own";
                    new_proteomics_line.Swissprot_id = (string)swissprot_id.Clone();
                    new_proteomics_line.Uniprot_accession = (string)uniprot_accession.Clone();
                    columnEntry = columnEntries[indexRatio_glomPT];
                    //if (!String.IsNullOrEmpty(columnEntry)) { columnEntry = "0"; }
                    new_proteomics_line.Value = Math.Log(double.Parse(columnEntry),2);
                    new_proteomics_line.PatientID = KPMP_data_integration_class.Get_combined_patients_label(); ;
                    new_proteomics_line.Value_type = KPMP_value_type_enum.Log2_ratioavg;
                    proteomics_list.Add(new_proteomics_line);
                }
                reader.Close();
            }
            this.Data = proteomics_list.ToArray();
        }
        public KPMP_nearSingleCell_proteomics_UCSF_class Deep_copy()
        {
            KPMP_nearSingleCell_proteomics_UCSF_class copy = (KPMP_nearSingleCell_proteomics_UCSF_class)this.MemberwiseClone();
            int data_length = this.Data.Length;
            copy.Data = new KPMP_nearSingleCell_proteomics_UCSF_line_class[data_length];
            for (int indexD = 0; indexD < data_length; indexD++)
            {
                copy.Data[indexD] = this.Data[indexD].Deep_copy();
            }
            return copy;
        }
    }

    class KPMP_nearSingleCell_proteomics_line_class
    {
        public string Gene_symbol { get; set; }
        public string Subsegment { get; set; }
        public string SampleName { get; set; }
        public string PatientId { get; set; }
        public string Dataset { get; set; }
        public string Kpmp_data_integration_term { get; set; }
        public float Value { get; set; }
        public KPMP_value_type_enum Value_type { get; set; }

        public KPMP_nearSingleCell_proteomics_line_class()
        {
            Gene_symbol = "";
            Subsegment = "";
            SampleName = "";
            PatientId = "";
            Dataset = "";
            Kpmp_data_integration_term = "";
        }

        public KPMP_nearSingleCell_proteomics_line_class Deep_copy()
        {
            KPMP_nearSingleCell_proteomics_line_class copy = (KPMP_nearSingleCell_proteomics_line_class)this.MemberwiseClone();
            copy.Gene_symbol = (string)this.Gene_symbol.Clone();
            copy.Dataset = (string)this.Dataset.Clone();
            copy.Subsegment = (string)this.Subsegment.Clone();
            copy.SampleName = (string)this.SampleName.Clone();
            copy.PatientId = (string)this.PatientId.Clone();
            copy.Kpmp_data_integration_term = (string)this.Kpmp_data_integration_term.Clone();
            return copy;
        }
    }

    class KPMP_nearSingleCell_proteomics_class
    {
        public KPMP_nearSingleCell_proteomics_line_class[] Data { get; set; }
        public string[] Bg_proteins_in_upper_case { get; set; }

        private void Add_to_array(KPMP_nearSingleCell_proteomics_line_class[] add_data)
        {
            int this_length = this.Data.Length;
            int add_length = add_data.Length;
            int new_length = this_length + add_length;
            KPMP_nearSingleCell_proteomics_line_class[] new_data = new KPMP_nearSingleCell_proteomics_line_class[new_length];
            int indexNew = -1;
            for (int indexThis = 0; indexThis < this_length; indexThis++)
            {
                indexNew++;
                new_data[indexNew] = this.Data[indexThis];
            }
            for (int indexAdd = 0; indexAdd < add_length; indexAdd++)
            {
                indexNew++;
                new_data[indexNew] = add_data[indexAdd];
            }
            this.Data = new_data;
        }

        private void Generate_bg_proteins_in_upper_case()
        {
            Data = Data.OrderBy(l => l.Gene_symbol).ToArray();
            int data_length = Data.Length;
            KPMP_nearSingleCell_proteomics_line_class proteomics_line;
            List<string> bg_proteins = new List<string>();
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                proteomics_line = Data[indexData];
                if ((indexData == 0)
                    || (!proteomics_line.Gene_symbol.Equals(Data[indexData - 1].Gene_symbol)))
                {
                    bg_proteins.Add(proteomics_line.Gene_symbol.ToUpper());
                }
            }
            this.Bg_proteins_in_upper_case = bg_proteins.ToArray();
        }

        public void Keep_top_x_lines_per_sampleName(int keep_top_x)
        {
            this.Data = this.Data.OrderBy(l => l.SampleName).ThenByDescending(l => l.Value).ToArray();
            int data_length = this.Data.Length;
            KPMP_nearSingleCell_proteomics_line_class proteomics_line;
            List<KPMP_nearSingleCell_proteomics_line_class> keep = new List<KPMP_nearSingleCell_proteomics_line_class>();
            int kept_lines_count = 0;
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                proteomics_line = this.Data[indexData];
                if ((indexData == 0) || (!proteomics_line.SampleName.Equals(Data[indexData - 1].SampleName)))
                {
                    kept_lines_count = 0;
                }
                if (kept_lines_count < keep_top_x)
                {
                    keep.Add(proteomics_line);
                    kept_lines_count++;
                }
            }
            this.Data = keep.ToArray();
        }

        public void Keep_only_lines_with_minimum_fold_change(float minimum_fold_change)
        {
            this.Data = this.Data.OrderBy(l => l.SampleName).ThenByDescending(l => l.Value).ToArray();
            int data_length = this.Data.Length;
            KPMP_nearSingleCell_proteomics_line_class proteomics_line;
            List<KPMP_nearSingleCell_proteomics_line_class> keep = new List<KPMP_nearSingleCell_proteomics_line_class>();
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                proteomics_line = this.Data[indexData];
                if (proteomics_line.Value >= minimum_fold_change)
                {
                    keep.Add(proteomics_line);
                }
            }
            this.Data = keep.ToArray();
        }

        private void Set_kpmp_integration_term()
        {
            int data_length = this.Data.Length;
            KPMP_nearSingleCell_proteomics_line_class proteomics_line;
            for (int indexD = 0; indexD < data_length; indexD++)
            {
                proteomics_line = this.Data[indexD];
                proteomics_line.Kpmp_data_integration_term = KPMP_data_integration_class.Get_kpmp_integration_term_plus_patientID(proteomics_line.Subsegment, proteomics_line.PatientId);
            }
        }

        private void Remove_duplicated_lines_by_keeping_higher_value()
        {
            this.Data = this.Data.OrderBy(l => l.Gene_symbol).ThenBy(l => l.Subsegment).ThenBy(l => l.SampleName).ThenByDescending(l => l.Value).ToArray();
            KPMP_nearSingleCell_proteomics_line_class proteomics_line;
            List<KPMP_nearSingleCell_proteomics_line_class> keep = new List<KPMP_nearSingleCell_proteomics_line_class>();
            int data_length = this.Data.Length;
            for (int indexD = 0; indexD < data_length; indexD++)
            {
                proteomics_line = this.Data[indexD];
                if ((indexD == 0)
                    || (!proteomics_line.Gene_symbol.Equals(Data[indexD - 1].Gene_symbol))
                    || (!proteomics_line.SampleName.Equals(Data[indexD - 1].SampleName))
                    || (!proteomics_line.Subsegment.Equals(Data[indexD - 1].Subsegment)))
                {
                    keep.Add(proteomics_line);
                }
            }
            this.Data = keep.ToArray();
        }

        public void Generate()
        {
            //Read_and_fill_array();
            Set_kpmp_integration_term();
            Remove_duplicated_lines_by_keeping_higher_value();
            Generate_bg_proteins_in_upper_case();
        }

        public string[] Get_bg_proteins_in_upper_case()
        {
            int bg_proteins_length = this.Bg_proteins_in_upper_case.Length;
            string[] bg_proteins = new string[bg_proteins_length];
            for (int indexP = 0; indexP < bg_proteins_length; indexP++)
            {
                bg_proteins[indexP] = (string)this.Bg_proteins_in_upper_case[indexP].Clone();
            }
            return bg_proteins;
        }


        public KPMP_nearSingleCell_proteomics_class Deep_copy()
        {
            int data_length = this.Data.Length;
            KPMP_nearSingleCell_proteomics_class copy = (KPMP_nearSingleCell_proteomics_class)this.MemberwiseClone();
            copy.Data = new KPMP_nearSingleCell_proteomics_line_class[data_length];
            for (int indexD = 0; indexD < data_length; indexD++)
            {
                copy.Data[indexD] = this.Data[indexD].Deep_copy();
            }
            int bg_genes_length = this.Bg_proteins_in_upper_case.Length;
            for (int indexBg = 0; indexBg < bg_genes_length; indexBg++)
            {
                copy.Bg_proteins_in_upper_case[indexBg] = (string)this.Bg_proteins_in_upper_case[indexBg].Clone();
            }

            return copy;
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class KPMP_singleRNASeqCluster_line_class
    {
        public string Gene_symbol { get; set; }
        public double Pvalue { get; set; }
        public double Avg_logFC { get; set; }
        public double Pvalue_adjusted { get; set; }

        public float Pct1 { get; set; }
        public float Pct2 { get; set; }
        public float Fractional_rank_pvalue { get; set; }
        public float Fractional_rank_absAvgLogFc { get; set; }
        public int Cluster_no { get; set; }
        public string Cluster_cell_type { get; set; }
        public string PatientId { get; set; }
        public string KPMP_data_integration_term { get; set; }
        public string Dataset { get; set; }
        public KPMP_value_type_enum Value_type { get; set; }

        #region Equal set
        public static KPMP_singleRNASeqCluster_line_class[] Order_by_set(KPMP_singleRNASeqCluster_line_class[] array)
        {
            array = array.OrderBy(l => l.Cluster_no).ThenBy(l => l.Cluster_cell_type).ThenBy(l => l.PatientId).ThenBy(l => l.KPMP_data_integration_term).ThenBy(l => l.Dataset).ThenBy(l => l.Value_type).ToArray();
            return array;
        }

        public static KPMP_singleRNASeqCluster_line_class[] Order_by_set_and_pvalue(KPMP_singleRNASeqCluster_line_class[] array)
        {
            array = array.OrderBy(l => l.Cluster_no).ThenBy(l => l.Cluster_cell_type).ThenBy(l => l.PatientId).ThenBy(l => l.KPMP_data_integration_term).ThenBy(l => l.Dataset).ThenBy(l => l.Value_type).ThenBy(l => l.Pvalue).ToArray();
            return array;
        }

        public static KPMP_singleRNASeqCluster_line_class[] Order_by_set_and_gene_symbol(KPMP_singleRNASeqCluster_line_class[] array)
        {
            array = array.OrderBy(l => l.Cluster_no).ThenBy(l => l.Cluster_cell_type).ThenBy(l => l.PatientId).ThenBy(l => l.KPMP_data_integration_term).ThenBy(l => l.Dataset).ThenBy(l => l.Value_type).ThenBy(l => l.Gene_symbol).ToArray();
            return array;
        }

        public bool Equal_set(KPMP_singleRNASeqCluster_line_class other)
        {
            bool equal = (this.Cluster_no.Equals(other.Cluster_no))
                         && (this.Cluster_cell_type.Equals(other.Cluster_cell_type))
                         && (this.PatientId.Equals(other.PatientId))
                         && (this.KPMP_data_integration_term.Equals(other.KPMP_data_integration_term))
                         && (this.Dataset.Equals(other.Dataset))
                         && (this.Value_type.Equals(other.Value_type));
            return equal;
        }
        #endregion

        public KPMP_singleRNASeqCluster_line_class()
        {
            this.Dataset = "";
            this.PatientId = "";
            this.KPMP_data_integration_term = "";
        }

        public KPMP_singleRNASeqCluster_line_class Deep_copy()
        {
            KPMP_singleRNASeqCluster_line_class copy = (KPMP_singleRNASeqCluster_line_class)this.MemberwiseClone();
            copy.Gene_symbol = (string)this.Gene_symbol.Clone();
            copy.Cluster_cell_type = (string)this.Cluster_cell_type.Clone();
            copy.PatientId = (string)this.PatientId.Clone();
            copy.Dataset = (string)this.Dataset.Clone();
            return copy;
        }
    }

    class KPMP_singleRNASeqCluster_readWriteOptions : ReadWriteOptions_base
    {
        public static string Get_bgGenesInUpperCase_complete_fileName(string center)
        {
            string bgGenes_complete_fileName = Global_directory_class.Experimental_data_directory + "SingleCellNucleus_bgGenes\\" + KPMP_dataset_name_class.Get_bgGenes_fileName_for_dataset(center);
            return bgGenes_complete_fileName;
        }

        public KPMP_singleRNASeqCluster_readWriteOptions(string completeFileName, bool read_cluster_number)
        {
            this.File = completeFileName;
            if (read_cluster_number)
            {
                this.Key_propertyNames = new string[] { "Gene_symbol", "Pvalue", "Avg_logFC", "Pvalue_adjusted", "Pct1", "Pct2", "Cluster_no", "Cluster_cell_type" };
                this.Key_columnNames = new string[] { "gene", "p_val", "avg_logFC", "p_val_adj", "pct.1", "pct.2", "cluster", "Cluster_cell type" };
                
                this.Key_propertyNames = new string[] { "Gene_symbol", "Pvalue", "Avg_logFC", "Pvalue_adjusted", "Cluster_no", "Cluster_cell_type" };
                this.Key_columnNames = new string[] { "Gene_symbol", "p_val", "avg_log2FC", "p_val_adj", "ClusterNo", "Cluster_cell type" };

                this.Key_propertyNames = new string[] { "Gene_symbol", "Pvalue", "Avg_logFC", "Pvalue_adjusted", "Cluster_no", "Cluster_cell_type" };
                this.Key_columnNames = new string[] { "gene", "p_val", "avg_logFC", "p_val_adj", "cluster", "Cluster_cell type" };
            }
            else
            {
                this.Key_propertyNames = new string[] { "Gene_symbol", "Pvalue", "Avg_logFC", "Pvalue_adjusted", "Pct1", "Pct2", "Cluster_cell_type" };
                this.Key_columnNames = new string[] { "gene", "p_val", "avg_logFC", "p_val_adj", "pct.1", "pct.2", "Cluster_cell type" };

                KPMP_analysis_set_enum analysis_set = Global_kpmp_class.Analysis_set;
                switch (analysis_set)
                {
                    case KPMP_analysis_set_enum.Disease_data_2021:
                        this.Key_propertyNames = new string[] { "Gene_symbol", "Pvalue", "Avg_logFC", "Pvalue_adjusted", "Cluster_cell_type" };
                        this.Key_columnNames = new string[] { "Gene_symbol", "p_val", "avg_log2FC", "p_val_adj", "Cluster_cell_type" };
                        break;
                    case KPMP_analysis_set_enum.Pilot_data:
                        this.Key_propertyNames = new string[] { "Gene_symbol", "Pvalue", "Avg_logFC", "Pvalue_adjusted", "Cluster_cell_type" };
                        this.Key_columnNames = new string[] { "gene", "p_val", "avg_logFC", "p_val_adj", "Cluster_cell type" };
                        break;
                    default:
                        break;
                }
            }

            this.HeadlineDelimiters = new char[] { Global_class.Tab, Global_class.Comma };
            this.LineDelimiters = new char[] { Global_class.Tab, Global_class.Comma };
            this.File_has_headline = true;
            this.Report = ReadWrite_report_enum.Report_main;
        }
    }

    class KPMP_own_singleRNASeqCluster_readWriteOptions : ReadWriteOptions_base
    {
        public static string Get_bgGenesInUpperCase_complete_fileName(string center)
        {
            string bgGenes_complete_fileName = Global_directory_class.Experimental_data_directory + "SingleCellNucleus_bgGenes\\" + KPMP_dataset_name_class.Get_bgGenes_fileName_for_dataset(center);
            return bgGenes_complete_fileName;
        }

        public KPMP_own_singleRNASeqCluster_readWriteOptions(string completeFileName, bool read_cluster_number)
        {
            this.File = completeFileName;
            if (read_cluster_number)
            {
                this.Key_propertyNames = new string[] { "Gene_symbol", "Pvalue", "Avg_logFC", "Pvalue_adjusted", "Pct1","Pct2", "Cluster_no", "Cluster_cell_type" };
                this.Key_columnNames = new string[] { "Gene_symbol", "p_val", "avg_logFC", "p_val_adj", "pct.1","pct.2","cluster", "Cluster_name" };
            }
            else
            {
                this.Key_propertyNames = new string[] { "Gene_symbol", "Pvalue", "Avg_logFC", "Pvalue_adjusted", "Pct1", "Pct2", "Cluster_cell_type" };
                this.Key_columnNames = new string[] { "Gene_symbol", "p_val", "avg_log2FC", "p_val_adj", "pct.1","pct.2","Cluster" };
            }

            this.HeadlineDelimiters = new char[] { Global_class.Tab, Global_class.Comma };
            this.LineDelimiters = new char[] { Global_class.Tab, Global_class.Comma };
            this.File_has_headline = true;
            this.Report = ReadWrite_report_enum.Report_main;
        }
    }

    class KPMP_singleRNASeqCluster_class
    {
        public KPMP_singleRNASeqCluster_line_class[] Data { get; set; }
        public string[] Bg_genes_in_upperCase { get; set; }
        public KPMP_integration_paper_metadata_class Dataset_patient { get; set; }

        public KPMP_singleRNASeqCluster_class()
        {
            this.Bg_genes_in_upperCase = new string[0];
            this.Data = new KPMP_singleRNASeqCluster_line_class[0];
            this.Dataset_patient = new KPMP_integration_paper_metadata_class();
        }

        public void Check_for_duplicated_symbols_within_the_same_set()
        {
            int data_length = Data.Length;
            KPMP_singleRNASeqCluster_line_class singleRNASeqCluster_line;
            KPMP_singleRNASeqCluster_line_class previous_singleRNASeqCluster_line;
            this.Data = KPMP_singleRNASeqCluster_line_class.Order_by_set_and_gene_symbol(this.Data);
            for (int indexD = 1; indexD < data_length; indexD++)
            {
                previous_singleRNASeqCluster_line = this.Data[indexD - 1];
                singleRNASeqCluster_line = this.Data[indexD];
                if ((singleRNASeqCluster_line.Equal_set(previous_singleRNASeqCluster_line))
                    && (singleRNASeqCluster_line.Gene_symbol.Equals(previous_singleRNASeqCluster_line.Gene_symbol)))
                { throw new Exception(); }
            }
        }

        public bool Remove_complete_duplicated_symbols_within_the_same_set_and_check_if_values_are_the_same()
        {
            int data_length = Data.Length;
            KPMP_singleRNASeqCluster_line_class singleRNASeqCluster_line;
            KPMP_singleRNASeqCluster_line_class previous_singleRNASeqCluster_line;
            this.Data = KPMP_singleRNASeqCluster_line_class.Order_by_set_and_gene_symbol(this.Data);
            List<KPMP_singleRNASeqCluster_line_class> keep = new List<KPMP_singleRNASeqCluster_line_class>();
            bool ok = true;
            for (int indexD = 1; indexD < data_length; indexD++)
            {
                previous_singleRNASeqCluster_line = this.Data[indexD - 1];
                singleRNASeqCluster_line = this.Data[indexD];
                if (  (indexD==0)
                    ||(!singleRNASeqCluster_line.Equal_set(previous_singleRNASeqCluster_line))
                    ||(!singleRNASeqCluster_line.Gene_symbol.Equals(previous_singleRNASeqCluster_line.Gene_symbol)))
                {
                    keep.Add(singleRNASeqCluster_line);
                }
                else
                {
                    if (   (!singleRNASeqCluster_line.Pvalue.Equals(previous_singleRNASeqCluster_line.Pvalue))
                        || (!singleRNASeqCluster_line.Pvalue_adjusted.Equals(previous_singleRNASeqCluster_line.Pvalue_adjusted))
                        || (!singleRNASeqCluster_line.Pct1.Equals(previous_singleRNASeqCluster_line.Pct1))
                        || (!singleRNASeqCluster_line.Pct2.Equals(previous_singleRNASeqCluster_line.Pct2))
                        || (!singleRNASeqCluster_line.Avg_logFC.Equals(previous_singleRNASeqCluster_line.Avg_logFC)))
                    {
                        ok = false;
                        throw new Exception();
                    }
                }
            }
            this.Data = keep.ToArray();
            return ok;
        }

        public void Calculate_fractional_ranks_for_pvalue_and_avg_logfc()
        {
            int data_length = this.Data.Length;
            KPMP_singleRNASeqCluster_line_class singleRNASeqCluster_line;
            KPMP_singleRNASeqCluster_line_class inner_singleRNASeqCluster_line;

            this.Data = this.Data.OrderBy(l => l.Dataset).ThenBy(l => l.Cluster_cell_type).ThenBy(l => l.Pvalue).ThenByDescending(l => Math.Abs(l.Avg_logFC)).ToArray();
            List<float> current_ranks = new List<float>();
            int current_rank = 0;
            float fractional_rank;
            int indexFirstIndexSameValues = 0;

            for (int indexD = 0; indexD < data_length; indexD++)
            {
                singleRNASeqCluster_line = this.Data[indexD];
                if ((indexD == 0)
                    || (!singleRNASeqCluster_line.Dataset.Equals(this.Data[indexD - 1].Dataset))
                    || (!singleRNASeqCluster_line.Cluster_cell_type.Equals(this.Data[indexD - 1].Cluster_cell_type)))
                {
                    current_rank = 0;
                }
                if ((indexD == 0)
                    || (!singleRNASeqCluster_line.Dataset.Equals(this.Data[indexD - 1].Dataset))
                    || (!singleRNASeqCluster_line.Cluster_cell_type.Equals(this.Data[indexD - 1].Cluster_cell_type))
                    || (!singleRNASeqCluster_line.Pvalue.Equals(this.Data[indexD - 1].Pvalue))
                    || (!Math.Abs(singleRNASeqCluster_line.Avg_logFC).Equals(Math.Abs(this.Data[indexD - 1].Avg_logFC))))
                {
                    current_ranks.Clear();
                    indexFirstIndexSameValues = indexD;
                }
                current_rank++;
                current_ranks.Add(current_rank);
                if ((indexD == data_length - 1)
                    || (!singleRNASeqCluster_line.Dataset.Equals(this.Data[indexD + 1].Dataset))
                    || (!singleRNASeqCluster_line.Cluster_cell_type.Equals(this.Data[indexD + 1].Cluster_cell_type))
                    || (!singleRNASeqCluster_line.Pvalue.Equals(this.Data[indexD + 1].Pvalue))
                    || (!Math.Abs(singleRNASeqCluster_line.Avg_logFC).Equals(Math.Abs(this.Data[indexD + 1].Avg_logFC))))
                {
                    fractional_rank = Math_class.Get_average(current_ranks.ToArray());
                    for (int indexInner = indexFirstIndexSameValues; indexInner <= indexD; indexInner++)
                    {
                        inner_singleRNASeqCluster_line = this.Data[indexInner];
                        inner_singleRNASeqCluster_line.Fractional_rank_pvalue = fractional_rank;
                    }
                }
            }

            this.Data = this.Data.OrderBy(l => l.Dataset).ThenBy(l => l.Cluster_cell_type).ThenByDescending(l => Math.Abs(l.Avg_logFC)).ThenBy(l => l.Pvalue).ToArray();

            for (int indexD = 0; indexD < data_length; indexD++)
            {
                singleRNASeqCluster_line = this.Data[indexD];
                if ((indexD == 0)
                    || (!singleRNASeqCluster_line.Dataset.Equals(this.Data[indexD - 1].Dataset))
                    || (!singleRNASeqCluster_line.Cluster_cell_type.Equals(this.Data[indexD - 1].Cluster_cell_type)))
                {
                    current_rank = 0;
                }
                if ((indexD == 0)
                    || (!singleRNASeqCluster_line.Dataset.Equals(this.Data[indexD - 1].Dataset))
                    || (!singleRNASeqCluster_line.Cluster_cell_type.Equals(this.Data[indexD - 1].Cluster_cell_type))
                    || (!singleRNASeqCluster_line.Pvalue.Equals(this.Data[indexD - 1].Pvalue))
                    || (!Math.Abs(singleRNASeqCluster_line.Avg_logFC).Equals(Math.Abs(this.Data[indexD - 1].Avg_logFC))))
                {
                    current_ranks.Clear();
                    indexFirstIndexSameValues = indexD;
                }
                current_rank++;
                current_ranks.Add(current_rank);
                if ((indexD == data_length - 1)
                    || (!singleRNASeqCluster_line.Dataset.Equals(this.Data[indexD + 1].Dataset))
                    || (!singleRNASeqCluster_line.Cluster_cell_type.Equals(this.Data[indexD + 1].Cluster_cell_type))
                    || (!singleRNASeqCluster_line.Pvalue.Equals(this.Data[indexD + 1].Pvalue))
                    || (!Math.Abs(singleRNASeqCluster_line.Avg_logFC).Equals(Math.Abs(this.Data[indexD + 1].Avg_logFC))))
                {
                    fractional_rank = Math_class.Get_average(current_ranks.ToArray());
                    for (int indexInner = indexFirstIndexSameValues; indexInner <= indexD; indexInner++)
                    {
                        inner_singleRNASeqCluster_line = this.Data[indexInner];
                        inner_singleRNASeqCluster_line.Fractional_rank_absAvgLogFc = fractional_rank;
                    }
                }
            }
        }

        private void Read_bg_genes_in_upperCase(string center)
        {
            string complete_bgGenes_directory = KPMP_singleRNASeqCluster_readWriteOptions.Get_bgGenesInUpperCase_complete_fileName(center);
            this.Bg_genes_in_upperCase = ReadWriteClass.Read_string_array(complete_bgGenes_directory);
            this.Bg_genes_in_upperCase = this.Bg_genes_in_upperCase.Distinct().ToArray();
        }

        private void Set_kpmp_data_integration_term()
        {
            int data_length = this.Data.Length;
            KPMP_singleRNASeqCluster_line_class data_line;
            for (int indexD = 0; indexD < data_length; indexD++)
            {
                data_line = this.Data[indexD];
                data_line.KPMP_data_integration_term = KPMP_data_integration_class.Get_kpmp_integration_term_plus_patientID(data_line.Cluster_cell_type, data_line.PatientId);
            }
        }

        private void Set_all_geneSymbols_to_upperCase()
        {
            foreach (KPMP_singleRNASeqCluster_line_class singleCell_line in this.Data)
            {
                singleCell_line.Gene_symbol = singleCell_line.Gene_symbol.ToUpper();
            }

        }

        public void Generate_and_override_current_array(string directory, string university)
        {
            Read_all_files_in_directory(directory, university);
            Generate_dataset_patient_instance();
            Set_all_geneSymbols_to_upperCase();
            Set_kpmp_data_integration_term();
            Read_bg_genes_in_upperCase(university);
        }

        public void Remove_selected_cluster(string[] remove_cluster_cell_types)
        {
            remove_cluster_cell_types = remove_cluster_cell_types.Distinct().ToArray();
            Dictionary<string, bool> remove_dict = new Dictionary<string, bool>();
            foreach (string remove_cluster in remove_cluster_cell_types)
            {
                remove_dict.Add(remove_cluster, true);
            }
            int data_length = this.Data.Length;
            KPMP_singleRNASeqCluster_line_class single_cell_line;
            List<KPMP_singleRNASeqCluster_line_class> keep = new List<KPMP_singleRNASeqCluster_line_class>();
            for (int indexData=0; indexData<data_length;indexData++)
            {
                single_cell_line = this.Data[indexData];
                if (!remove_dict.ContainsKey(single_cell_line.Cluster_cell_type))
                {
                    keep.Add(single_cell_line);
                }
            }
            this.Data = keep.ToArray();
        }

        public void Collapse_by_calculating_geometric_means_for_pvalue_and_fold_change_set_minimum_pvalue_as_10Eminus308(params int[] forced_dataset_length)
        {
            string[] datasets = Get_all_unique_datasets();
            int datasets_length = datasets.Length;
            if (forced_dataset_length.Length>1) { throw new Exception(); }
            else if (forced_dataset_length.Length==1) { datasets_length = forced_dataset_length[0]; }
            int data_length = this.Data.Length;
            KPMP_singleRNASeqCluster_line_class cluster_line;
            KPMP_singleRNASeqCluster_line_class collapsed_cluster_line;
            List<KPMP_singleRNASeqCluster_line_class> collapsed_data = new List<KPMP_singleRNASeqCluster_line_class>();
            this.Data = this.Data.OrderBy(l => l.Gene_symbol).ThenBy(l=>l.Cluster_cell_type).ToArray();
            List<double> minLog10Pvalues = new List<double>();
            List<double> nonZero_minLog10Pvalues = new List<double>();
            List<double> minLog10AdjPvalues = new List<double>();
            List<double> logfoldchanges = new List<double>();
            List<float> pct1s = new List<float>();
            List<float> pct2s = new List<float>();
            double add_minLog10Pvalue;
            double add_minLog10AdjPvalue;
            double average_minLog10Pvalue;
            double average_minLog10AdjPvalue;
            double average_logfoldchange;
            float average_pct1;
            float average_pct2;
            int at_least_two_nonZero_pvalues_are_identical_count = 0;
            int genes_count = 0;
            for (int indexData=0; indexData<data_length;indexData++)
            {
                cluster_line = this.Data[indexData];
                if (  (indexData==0)
                    || (!cluster_line.Gene_symbol.Equals(this.Data[indexData - 1].Gene_symbol))
                    || (!cluster_line.Cluster_cell_type.Equals(this.Data[indexData - 1].Cluster_cell_type)))
                {
                    minLog10Pvalues.Clear();
                    minLog10AdjPvalues.Clear();
                    logfoldchanges.Clear();
                    pct1s.Clear();
                    pct2s.Clear();
                    genes_count++;
                }
                if (cluster_line.Pvalue != 0) { add_minLog10Pvalue = -Math.Log10(cluster_line.Pvalue); }
                else { add_minLog10Pvalue = 308; }
                if (add_minLog10Pvalue > 308) { add_minLog10Pvalue = 308; }
                minLog10Pvalues.Add(add_minLog10Pvalue);
                if (cluster_line.Pvalue_adjusted != 0) { add_minLog10AdjPvalue  = - Math.Log10(cluster_line.Pvalue_adjusted); }
                else { add_minLog10AdjPvalue = 308; }
                if (add_minLog10AdjPvalue>308) { add_minLog10AdjPvalue = 308; }
                minLog10AdjPvalues.Add(add_minLog10AdjPvalue);

                logfoldchanges.Add(cluster_line.Avg_logFC);
                pct1s.Add(cluster_line.Pct1);
                pct2s.Add(cluster_line.Pct2);
                if (   (indexData == data_length-1)
                    || (!cluster_line.Gene_symbol.Equals(this.Data[indexData+1].Gene_symbol))
                    || (!cluster_line.Cluster_cell_type.Equals(this.Data[indexData + 1].Cluster_cell_type)))
                {
                    nonZero_minLog10Pvalues.Clear();
                    foreach (double minLog10Pvalue in minLog10AdjPvalues)
                    {
                        if (minLog10Pvalue!=0) { nonZero_minLog10Pvalues.Add(minLog10Pvalue); }
                    }
                    if (nonZero_minLog10Pvalues.Distinct().ToArray().Length!=nonZero_minLog10Pvalues.Count) { at_least_two_nonZero_pvalues_are_identical_count++; }
                    while (minLog10AdjPvalues.Count<datasets_length)
                    {
                        minLog10Pvalues.Add(0);
                        minLog10AdjPvalues.Add(0);
                        logfoldchanges.Add(0);
                        //pct1s.Add(cluster_line.Pct1);
                        //pct2s.Add(cluster_line.Pct2);
                    }
                    average_minLog10Pvalue = Math_class.Get_average(minLog10Pvalues.ToArray());
                    average_minLog10AdjPvalue = Math_class.Get_average(minLog10AdjPvalues.ToArray());
                    average_logfoldchange = Math_class.Get_average(logfoldchanges.ToArray());
                    average_pct1 = Math_class.Get_average(pct1s.ToArray());
                    average_pct2 = Math_class.Get_average(pct2s.ToArray());
                    collapsed_cluster_line = cluster_line.Deep_copy();
                    collapsed_cluster_line.Fractional_rank_absAvgLogFc = -1;
                    collapsed_cluster_line.Fractional_rank_pvalue = -1;
                    collapsed_cluster_line.Pvalue = Math.Pow(10, -average_minLog10Pvalue);
                    collapsed_cluster_line.Pvalue_adjusted = Math.Pow(10, -average_minLog10AdjPvalue);
                    collapsed_cluster_line.Avg_logFC = average_logfoldchange;
                    collapsed_cluster_line.Pct1 = average_pct1;
                    collapsed_cluster_line.Pct2 = average_pct2;
                    collapsed_data.Add(collapsed_cluster_line);
                }
            }
            this.Data = collapsed_data.ToArray();
            double fraction_of_at_least_two_nonZero_pvalues_are_identical_count = (double)at_least_two_nonZero_pvalues_are_identical_count / (double)genes_count;
            if (fraction_of_at_least_two_nonZero_pvalues_are_identical_count > 0.2) { throw new Exception(); }
        }

        public string[] Get_all_unique_datasets()
        {
            List<string> datasets = new List<string>();
            foreach (KPMP_singleRNASeqCluster_line_class data_line in this.Data)
            {
                datasets.Add(data_line.Dataset);
            }
            return datasets.Distinct().OrderBy(l => l).ToArray();
        }

        public string[] Get_deep_copy_of_bg_genes_in_upperCase()
        {
            int bg_genes_length = this.Bg_genes_in_upperCase.Length;
            string[] bg_genes = new string[bg_genes_length];
            for (int indexBG = 0; indexBG < bg_genes_length; indexBG++)
            {
                bg_genes[indexBG] = (string)this.Bg_genes_in_upperCase[indexBG].Clone();
            }
            return bg_genes;
        }

        private void Read_all_files_in_directory(string directory, string dataset)
        {
            KPMP_singleRNASeqCluster_line_class[] add_data;
            List<KPMP_singleRNASeqCluster_line_class> data_list = new List<KPMP_singleRNASeqCluster_line_class>();
            string[] completeFileNames = Directory.GetFiles(directory);
            string completeFileName;
            int completeFileNames_length = completeFileNames.Length;
            string patientId;
            string[] splitStrings;
            bool readCluster_no = false;
            for (int indexC = 0; indexC < completeFileNames_length; indexC++)
            {
                completeFileName = completeFileNames[indexC];
                //if (completeFileName.ToUpper().IndexOf("AllPatients".ToUpper()) != -1) { readCluster_no = true; }
                //else { readCluster_no = false; }
                KPMP_singleRNASeqCluster_readWriteOptions readWriteOptions = new KPMP_singleRNASeqCluster_readWriteOptions(completeFileName, readCluster_no);
                add_data = ReadWriteClass.ReadRawData_and_FillArray<KPMP_singleRNASeqCluster_line_class>(readWriteOptions);
                splitStrings = Path.GetFileNameWithoutExtension(completeFileName).Split('_');
                patientId = splitStrings[splitStrings.Length - 1];
                foreach (KPMP_singleRNASeqCluster_line_class add_data_line in add_data)
                {
                    add_data_line.PatientId = (string)patientId.Clone();
                    add_data_line.Dataset = (string)dataset.Clone();
                    add_data_line.Value_type = KPMP_value_type_enum.Single_value;
                }
                data_list.AddRange(add_data);
            }
            this.Data = data_list.ToArray();
        }

        public void Read_without_clusterCellType(string subdirectory, string fileName)
        {
            string completeFileName = Global_directory_class.Experimental_data_directory + subdirectory + fileName;
            KPMP_singleRNASeqCluster_readWriteOptions readWriteOptions = new KPMP_singleRNASeqCluster_readWriteOptions(completeFileName, true);
            List<string> columnNames = new List<string>();
            List<string> propertyNames = new List<string>();
            string removePropertyName = "Cluster_cell_type";
            int propertyNames_length = readWriteOptions.Key_propertyNames.Length;
            string propertyName;
            string columnName;
            for (int indexProp = 0; indexProp < propertyNames_length; indexProp++)
            {
                propertyName = readWriteOptions.Key_propertyNames[indexProp];
                columnName = readWriteOptions.Key_columnNames[indexProp];
                if (!propertyName.Equals(removePropertyName))
                {
                    propertyNames.Add(propertyName);
                    columnNames.Add(columnName);
                }
            }
            readWriteOptions.Key_propertyNames = propertyNames.ToArray();
            readWriteOptions.Key_columnNames = columnNames.ToArray();
            this.Data = ReadWriteClass.ReadRawData_and_FillArray<KPMP_singleRNASeqCluster_line_class>(readWriteOptions);
        }

        public void Read(string subdirectory, string fileName)
        {
            string completeFileName = Global_directory_class.Experimental_data_directory + subdirectory + fileName;
            KPMP_singleRNASeqCluster_readWriteOptions readWriteOptions = new KPMP_singleRNASeqCluster_readWriteOptions(completeFileName, true);
            this.Data = ReadWriteClass.ReadRawData_and_FillArray<KPMP_singleRNASeqCluster_line_class>(readWriteOptions);
        }

        public void Read_own(string subdirectory, params string[] fileNames)
        {
            int fileNames_length = fileNames.Length;
            string dataset;
            KPMP_singleRNASeqCluster_line_class[] add_to_array;
            List<KPMP_singleRNASeqCluster_line_class> new_data = new List<KPMP_singleRNASeqCluster_line_class>();
            for (int indexFileName = 0; indexFileName < fileNames_length; indexFileName++)
            {
                string completeFileName = Global_directory_class.Experimental_data_directory + subdirectory + fileNames[indexFileName];
                dataset = Path.GetFileNameWithoutExtension(completeFileName);
                KPMP_own_singleRNASeqCluster_readWriteOptions readWriteOptions = new KPMP_own_singleRNASeqCluster_readWriteOptions(completeFileName, false);
                add_to_array = ReadWriteClass.ReadRawData_and_FillArray<KPMP_singleRNASeqCluster_line_class>(readWriteOptions);
                foreach (KPMP_singleRNASeqCluster_line_class add_line in add_to_array)
                {
                    add_line.Dataset = (string)dataset.Clone();
                }
                new_data.AddRange(add_to_array);
            }
            this.Data = new_data.ToArray();
        }

        public void Write(string subdirectory, string fileName)
        {
            string completeFileName = Global_directory_class.Experimental_data_directory + subdirectory + fileName;
            KPMP_singleRNASeqCluster_readWriteOptions readWriteOptions = new KPMP_singleRNASeqCluster_readWriteOptions(completeFileName, true);
            ReadWriteClass.WriteData(this.Data, readWriteOptions);
        }

        public void Write_into_results_directory(string subdirectory, string fileName)
        {
            string completeFileName = Global_directory_class.Results_directory + subdirectory + fileName;
            KPMP_singleRNASeqCluster_readWriteOptions readWriteOptions = new KPMP_singleRNASeqCluster_readWriteOptions(completeFileName, true);
            ReadWriteClass.WriteData(this.Data, readWriteOptions);
        }

        public void Keep_only_top_x_lines_per_sample_based_on_pvalue_and_descending_abs_logfc(int top_x)
        {
            int data_length = this.Data.Length;
            List<KPMP_singleRNASeqCluster_line_class> keep = new List<KPMP_singleRNASeqCluster_line_class>();
            KPMP_singleRNASeqCluster_line_class data_line;
            this.Data = this.Data.OrderBy(l => l.Cluster_no).ThenBy(l => l.Cluster_cell_type).ThenBy(l => l.PatientId).ThenBy(l => l.KPMP_data_integration_term).ThenBy(l => l.Dataset).ThenBy(l => l.Value_type).ThenBy(l => l.Pvalue).ThenByDescending(l=>Math.Abs(l.Avg_logFC)).ToArray();
            int kept_lines_count = 0;
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                data_line = this.Data[indexData];
                if (  (indexData == 0)
                    ||(!data_line.Cluster_no.Equals(this.Data[indexData-1].Cluster_no))
                    ||(!data_line.Cluster_cell_type.Equals(this.Data[indexData - 1].Cluster_cell_type))
                    ||(!data_line.PatientId.Equals(this.Data[indexData - 1].PatientId))
                    ||(!data_line.KPMP_data_integration_term.Equals(this.Data[indexData - 1].KPMP_data_integration_term))
                    ||(!data_line.Dataset.Equals(this.Data[indexData - 1].Dataset))
                    ||(!data_line.Value_type.Equals(this.Data[indexData - 1].Value_type)))
                {
                    kept_lines_count = 0;
                }
                if (kept_lines_count < top_x)
                {
                    keep.Add(data_line);
                    kept_lines_count++;
                }
                else if (kept_lines_count==top_x)
                {
                    //if (data_line.Pvalue.Equals(this.Data[indexData-1].Pvalue)) { throw new Exception(); }
                    if (double.IsInfinity(data_line.Pvalue)) { throw new Exception(); }
                }
            }
            this.Data = keep.ToArray();
        }

        public void Keep_only_top_x_lines_per_sample_based_on_descending_selected_valueType(int keep_x, KPMP_value_type_enum selected_valueType)
        {
            switch (selected_valueType)
            {
                case KPMP_value_type_enum.Minus_log10_pvalue:
                    Keep_only_top_x_lines_per_sample_based_on_pvalue_and_descending_abs_logfc(keep_x);
                    break;
                default:
                    throw new Exception();
            }
        }

        public void Keep_only_lines_with_pvalue_below_alpha(float alpha)
        {
            int data_length = this.Data.Length;
            List<KPMP_singleRNASeqCluster_line_class> keep = new List<KPMP_singleRNASeqCluster_line_class>();
            KPMP_singleRNASeqCluster_line_class data_line;
            this.Data = this.Data.OrderBy(l => l.PatientId).ThenBy(l => l.Cluster_no).ThenBy(l => l.Pvalue).ToArray();
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                data_line = this.Data[indexData];
                if (data_line.Pvalue <= alpha)
                {
                    keep.Add(data_line);
                }
            }
            this.Data = keep.ToArray();
        }

        public void Keep_only_upregulated_genes()
        {
            int data_length = this.Data.Length;
            List<KPMP_singleRNASeqCluster_line_class> keep = new List<KPMP_singleRNASeqCluster_line_class>();
            KPMP_singleRNASeqCluster_line_class data_line;
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                data_line = this.Data[indexData];
                if (data_line.Avg_logFC>0)
                {
                    keep.Add(data_line);
                }
            }
            this.Data = keep.ToArray();
        }

        public void Keep_only_lines_with_adjusted_pvalue_below_alpha(double alpha)
        {
            int data_length = this.Data.Length;
            List<KPMP_singleRNASeqCluster_line_class> keep = new List<KPMP_singleRNASeqCluster_line_class>();
            KPMP_singleRNASeqCluster_line_class data_line;
            this.Data = this.Data.OrderBy(l => l.PatientId).ThenBy(l => l.Cluster_no).ThenBy(l => l.Pvalue).ToArray();
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                data_line = this.Data[indexData];
                if (data_line.Pvalue_adjusted <= alpha)
                {
                    keep.Add(data_line);
                }
            }
            this.Data = keep.ToArray();
        }

        public void Keep_only_lines_with_indicated_patientIDs(params string[] patientIDs)
        {
            patientIDs = patientIDs.Distinct().OrderBy(l => l).ToArray();
            int indexPatientID = 0;
            int patientIDs_length = patientIDs.Length;
            string patientID;
            int stringCompare = -2;

            int data_length = this.Data.Length;
            List<KPMP_singleRNASeqCluster_line_class> keep = new List<KPMP_singleRNASeqCluster_line_class>();
            KPMP_singleRNASeqCluster_line_class data_line;
            this.Data = this.Data.OrderBy(l => l.PatientId).ThenBy(l => l.Cluster_no).ThenBy(l => l.Pvalue).ToArray();
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                data_line = this.Data[indexData];
                stringCompare = -2;
                while ((indexPatientID < patientIDs_length) && (stringCompare < 0))
                {
                    patientID = patientIDs[indexPatientID];
                    stringCompare = patientID.CompareTo(data_line.PatientId);
                    if (stringCompare < 0)
                    {
                        indexPatientID++;
                    }
                    else if (stringCompare == 0)
                    {
                        keep.Add(data_line);
                    }
                }
            }
            this.Data = keep.ToArray();
        }

        public void Keep_only_lines_with_cluster_cell_types_containing_substrings(params string[] cluster_cell_type_substrings)
        {
            cluster_cell_type_substrings = cluster_cell_type_substrings.Distinct().OrderBy(l => l).ToArray();
            int cluster_cell_type_substrings_length = cluster_cell_type_substrings.Length;
            string cluster_cell_type_substring;
            bool keep;
            int data_length = this.Data.Length;
            List<KPMP_singleRNASeqCluster_line_class> keep_lines = new List<KPMP_singleRNASeqCluster_line_class>();
            KPMP_singleRNASeqCluster_line_class data_line;
            this.Data = this.Data.OrderBy(l => l.PatientId).ThenBy(l => l.Cluster_no).ThenBy(l => l.Pvalue).ToArray();
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                data_line = this.Data[indexData];
                keep = false;
                for (int indexSubstring=0; indexSubstring < cluster_cell_type_substrings_length; indexSubstring++)
                {
                    cluster_cell_type_substring = cluster_cell_type_substrings[indexSubstring];
                    if (data_line.Cluster_cell_type.IndexOf(cluster_cell_type_substring)!=-1)
                    {
                        keep = true;
                        break;
                    }
                }
                if (keep)
                {
                    keep_lines.Add(data_line);
                }
            }
            this.Data = keep_lines.ToArray();
        }

        public KPMP_integration_paper_metadata_class Get_deep_copy_of_dataset_patient_instance()
        {
            return Dataset_patient.Deep_copy();
        }

        private void Generate_dataset_patient_instance()
        {
            Dictionary<string, Dictionary<string, bool>> dataset_patient_considered_dict = new Dictionary<string, Dictionary<string, bool>>();
            int data_length = this.Data.Length;
            KPMP_integration_paper_metadata_line_class new_dataset_patient_line;
            List<KPMP_integration_paper_metadata_line_class> new_dataset_patient_lines = new List<KPMP_integration_paper_metadata_line_class>();
            KPMP_singleRNASeqCluster_line_class singleClister_line;
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                singleClister_line = this.Data[indexData];
                if (!dataset_patient_considered_dict.ContainsKey(singleClister_line.Dataset))
                {
                    dataset_patient_considered_dict.Add(singleClister_line.Dataset, new Dictionary<string, bool>());
                }
                if (!dataset_patient_considered_dict[singleClister_line.Dataset].ContainsKey(singleClister_line.PatientId))
                {
                    dataset_patient_considered_dict[singleClister_line.Dataset].Add(singleClister_line.PatientId, true);
                    new_dataset_patient_line = new KPMP_integration_paper_metadata_line_class();
                    new_dataset_patient_line.Dataset = (string)singleClister_line.Dataset.Clone();
                    new_dataset_patient_line.Libraries =  new string[] { (string)singleClister_line.PatientId.Clone() };
                    new_dataset_patient_lines.Add(new_dataset_patient_line);
                }
            }
            if (Dataset_patient.Documentations.Length!=0) { throw new Exception(); }
            Dataset_patient = new KPMP_integration_paper_metadata_class();
            Dataset_patient.Add_to_array(new_dataset_patient_lines.ToArray());
        }

        public KPMP_standardized_dataset_class Generate_standardized_dataset_on_expression_values(KPMP_value_type_enum value_type_1st, KPMP_value_type_enum value_type_2nd)
        {
            int data_length = this.Data.Length;
            List<KPMP_value_type_enum> value_types = new List<KPMP_value_type_enum>();
            this.Data = this.Data.OrderBy(l => l.Value_type).ThenBy(l => l.Dataset).ToArray();
            KPMP_singleRNASeqCluster_line_class data_line;
            KPMP_standardized_dataset_line_class standardized_data_line;
            List<KPMP_standardized_dataset_line_class> standardized_data_list = new List<KPMP_standardized_dataset_line_class>();
            List<string> current_dataset_symbols = new List<string>();
            string dataset = (string)this.Data[0].Dataset.Clone();
            for (int indexD = 0; indexD < data_length; indexD++)
            {
                data_line = this.Data[indexD];
                if ((indexD == 0) || (!data_line.Value_type.Equals(this.Data[indexD - 1].Value_type)))
                {
                    value_types.Add(data_line.Value_type);
                }
                if (!data_line.Dataset.Equals(dataset)) { throw new Exception(); }
                standardized_data_line = new KPMP_standardized_dataset_line_class();
                standardized_data_line.Cell_segment = (string)data_line.Cluster_cell_type.Clone();
                standardized_data_line.Dataset = (string)data_line.Dataset.Clone();
                standardized_data_line.KPMP_data_integration_term = (string)data_line.KPMP_data_integration_term.Clone();
                standardized_data_line.PatientId = (string)data_line.PatientId.Clone();
                standardized_data_line.Gene_symbol = (string)data_line.Gene_symbol.Clone();

                standardized_data_line.Value_type_1st = value_type_1st;
                switch (standardized_data_line.Value_type_1st)
                {
                    case KPMP_value_type_enum.Log_ratioavg:
                        standardized_data_line.Value_1st = data_line.Avg_logFC;
                        break;
                    case KPMP_value_type_enum.Minus_log10_pvalue:
                        standardized_data_line.Value_1st = -Math.Log10(data_line.Pvalue); ;
                        break;
                    case KPMP_value_type_enum.Minus_log10_pvalue_adjusted:
                        standardized_data_line.Value_1st = -Math.Log10(data_line.Pvalue_adjusted); 
                        break;
                    default:
                        throw new Exception();
                }

                standardized_data_line.Value_type_2nd = value_type_2nd;
                switch (standardized_data_line.Value_type_2nd)
                {
                    case KPMP_value_type_enum.Log_ratioavg:
                        standardized_data_line.Value_2nd = data_line.Avg_logFC;
                        break;
                    case KPMP_value_type_enum.Minus_log10_pvalue:
                        standardized_data_line.Value_2nd = -Math.Log10(data_line.Pvalue); ;
                        break;
                    case KPMP_value_type_enum.Minus_log10_pvalue_adjusted:
                        standardized_data_line.Value_2nd = -Math.Log10(data_line.Pvalue_adjusted);
                        break;
                    default:
                        throw new Exception();
                }
                standardized_data_list.Add(standardized_data_line);
            }
            Dictionary<string, string[]> dataset_bgGenesProteins_dict = new Dictionary<string, string[]>();
            dataset_bgGenesProteins_dict.Add(dataset, this.Bg_genes_in_upperCase);

            KPMP_standardized_dataset_class standardized_data = new KPMP_standardized_dataset_class();
            standardized_data.Add_to_existing_instances(standardized_data_list.ToArray(), dataset_bgGenesProteins_dict);
            return standardized_data;
        }

        public DE_class Generate_de_instance_and_fill_with_minusLog10Pvalue()
        {
            int data_length = this.Data.Length;
            Fill_de_line_class fill_de_line;
            List<Fill_de_line_class> fill_de_list = new List<Fill_de_line_class>();
            KPMP_singleRNASeqCluster_line_class data_line;
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                data_line = this.Data[indexData];
                fill_de_line = new Fill_de_line_class();
                fill_de_line.Names_for_de = KPMP_data_integration_class.Get_names_for_de_in_correct_order(data_line.Cluster_cell_type + "_no" + data_line.Cluster_no, data_line.PatientId, data_line.Dataset, KPMP_value_type_enum.Minus_log10_pvalue, data_line.KPMP_data_integration_term);
                fill_de_line.Value_for_de = (float)-Math.Log10(data_line.Pvalue);
                fill_de_line.Symbols_for_de = new string[] { (string)data_line.Gene_symbol.Clone() };
                fill_de_list.Add(fill_de_line);
            }
            DE_class de = new DE_class();
            de.Fill_with_data_alternatively(fill_de_list.ToArray());
            de.Set_symbols_to_upper_case_letters();
            return de;
        }

        public KPMP_singleRNASeqCluster_class Deep_copy()
        {
            KPMP_singleRNASeqCluster_class copy = (KPMP_singleRNASeqCluster_class)this.MemberwiseClone();
            int single_length = this.Data.Length;
            copy.Data = new KPMP_singleRNASeqCluster_line_class[single_length];
            for (int indexD = 0; indexD < single_length; indexD++)
            {
                copy.Data[indexD] = this.Data[indexD].Deep_copy();
            }
            int bg_length = this.Bg_genes_in_upperCase.Length;
            copy.Bg_genes_in_upperCase = new string[bg_length];
            for (int indexBg = 0; indexBg < bg_length; indexBg++)
            {
                copy.Bg_genes_in_upperCase[indexBg] = (string)this.Bg_genes_in_upperCase[indexBg].Clone();
            }
            return copy;
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    #region Old
    class KPMP_single_rnaSeq_line_class
    {
        public string Gene { get; set; }
        public float p_val { get; set; }
        public float p_val_adj { get; set; }
        public float avg_logFC { get; set; }
        public int cluster { get; set; }
        public string Cell_name { get; set; }
        public string Kpmp_data_integration_term { get; set; }

        public KPMP_single_rnaSeq_line_class Deep_copy()
        {
            KPMP_single_rnaSeq_line_class copy = (KPMP_single_rnaSeq_line_class)this.MemberwiseClone();
            copy.Gene = (string)this.Gene.Clone();
            copy.Cell_name = (string)this.Cell_name.Clone();
            copy.Kpmp_data_integration_term = (string)this.Kpmp_data_integration_term.Clone();
            return copy;
        }
    }

    class KPMP_single_rnaSeq_readOptions_class : ReadWriteOptions_base
    {
        public KPMP_single_rnaSeq_readOptions_class(string file_name)
        {
            this.File = Global_directory_class.Experimental_data_directory + "RNASeq\\" + file_name;
            this.Key_propertyNames = new string[] { "Gene", "p_val", "p_val_adj", "avg_logFC", "cluster" };
            this.Key_columnNames = this.Key_propertyNames;
            HeadlineDelimiters = new char[] { Global_class.Tab };
            LineDelimiters = new char[] { Global_class.Tab };
            File_has_headline = true;
            Report = ReadWrite_report_enum.Report_main;
        }
    }

    class KPMP_single_rnaSeq_class
    {
        public KPMP_single_rnaSeq_line_class[] Single_data { get; set; }

        public void Generate_by_reading(string file_name)
        {
            Read(file_name);
        }

        public void Keep_only_indicated_clusters(params int[] keep_cluters)
        {
            List<KPMP_single_rnaSeq_line_class> keep = new List<KPMP_single_rnaSeq_line_class>();
            foreach (KPMP_single_rnaSeq_line_class rnaSeq_line in Single_data)
            {
                if (keep_cluters.Contains(rnaSeq_line.cluster))
                {
                    keep.Add(rnaSeq_line);
                }
            }
            this.Single_data = keep.ToArray();
        }

        public void Keep_top_s_genes_based_on_pvalue_per_cluster(int top_x)
        {
            int single_data_length = Single_data.Length;
            KPMP_single_rnaSeq_line_class single_line;
            List<KPMP_single_rnaSeq_line_class> keep = new List<KPMP_single_rnaSeq_line_class>();
            Single_data = Single_data.OrderBy(l => l.cluster).ThenBy(l => l.p_val).ToArray();
            int kept_lines_count = 0;
            for (int indexSingle = 0; indexSingle < single_data_length; indexSingle++)
            {
                single_line = this.Single_data[indexSingle];
                if ((indexSingle == 0)
                    || (!single_line.cluster.Equals(this.Single_data[indexSingle - 1].cluster)))
                {
                    kept_lines_count = 0;
                }
                if (kept_lines_count < top_x)
                {
                    kept_lines_count++;
                    keep.Add(single_line);
                }
            }
            this.Single_data = keep.ToArray();
        }

        private void Read(string file_name)
        {
            KPMP_single_rnaSeq_readOptions_class readOptions = new KPMP_single_rnaSeq_readOptions_class(file_name);
            this.Single_data = ReadWriteClass.ReadRawData_and_FillArray<KPMP_single_rnaSeq_line_class>(readOptions);
        }
    }
    #endregion

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class KPMP_metabolite_line_class
    {
        public string[] Associated_proteins { get; set; }
        public string[] Potential_metabolites { get; set; }
        public string Molecular_formula { get; set; }
        public float Glom_correlation { get; set; }
        public string[] Hmdb_ids { get; set; }
        public string[] Kegg_ids { get; set; }
        public int Is_glomerula_origin { get; set; }
        public string Molecular_class { get; set; }
        public string KPMP_data_integration_term { get; set; }
        public string PatientID { get; set; }

        public string ReadWrite_hmdb_ids
        {
            get { return ReadWriteClass.Get_writeLine_from_array(this.Hmdb_ids, KPMP_metabolite_results_readWriteOptions.Delimiter); }
            set
            {
                string input_line = value;
                string[] splitStrings = input_line.Split(';');
                List<string> hmdb_ids_list = new List<string>();
                foreach (string splitString in splitStrings)
                {
                    if (splitString.IndexOf("HMDB") != -1)
                    {
                        hmdb_ids_list.Add(splitString);
                    }
                }
                this.Hmdb_ids = hmdb_ids_list.ToArray();
            }
        }
        public string ReadWrite_associated_proteins
        {
            get { return ReadWriteClass.Get_writeLine_from_array(this.Associated_proteins, KPMP_metabolite_results_readWriteOptions.Delimiter); }
            set { this.Associated_proteins = ReadWriteClass.Get_array_from_readLine<string>(value, KPMP_metabolite_results_readWriteOptions.Delimiter); }
        }
        public string ReadWrite_potential_metabolites
        {
            get { return ReadWriteClass.Get_writeLine_from_array(this.Potential_metabolites, KPMP_metabolite_results_readWriteOptions.Delimiter); }
            set { this.Potential_metabolites = ReadWriteClass.Get_array_from_readLine<string>(value, KPMP_metabolite_results_readWriteOptions.Delimiter); }
        }


        #region Equal set
        public static KPMP_metabolite_line_class[] Order_by_set(KPMP_metabolite_line_class[] array)
        {
            array = array.OrderBy(l => l.PatientID).ThenBy(l => l.KPMP_data_integration_term).ThenBy(l => l.Is_glomerula_origin).ToArray();
            return array;
        }
        #endregion

        public KPMP_metabolite_line_class()
        {
            this.Associated_proteins = new string[0];
            this.Potential_metabolites = new string[0];
            this.Hmdb_ids = new string[0];
            Molecular_class = "";
            Molecular_formula = "";
            KPMP_data_integration_term = "";
        }

        public KPMP_metabolite_line_class Deep_copy()
        {
            KPMP_metabolite_line_class copy = (KPMP_metabolite_line_class)this.MemberwiseClone();
            copy.Molecular_class = (string)this.Molecular_class.Clone();
            copy.Molecular_formula = (string)this.Molecular_formula.Clone();
            copy.Hmdb_ids = Array_class.Deep_copy_string_array(this.Hmdb_ids);
            copy.Associated_proteins = Array_class.Deep_copy_string_array(this.Associated_proteins);
            copy.Potential_metabolites = Array_class.Deep_copy_string_array(this.Potential_metabolites);
            copy.KPMP_data_integration_term = (string)this.KPMP_data_integration_term.Clone();
            return copy;
        }
    }

    class KPMP_metabolite_inputReadOptions : ReadWriteOptions_base
    {
        public KPMP_metabolite_inputReadOptions(string fileName)
        {
            this.File = Global_directory_class.Experimental_data_directory + "Metabolomics_2018December15\\" + fileName;
            this.Key_propertyNames = new string[] { "Molecular_formula", "Glom_correlation", "ReadWrite_hmdb_ids", "Is_glomerula_origin" };
            this.Key_columnNames = new string[] { "Molecular formula", "Glom correlation", "HMDB ID(s)", "Glomerular (1: yes 0: no)" };
            this.HeadlineDelimiters = new char[] { Global_class.Comma };
            this.LineDelimiters = new char[] { Global_class.Comma };
            this.File_has_headline = true;
            this.Report = ReadWrite_report_enum.Report_main;
        }
    }

    class KPMP_metabolite_results_readWriteOptions : ReadWriteOptions_base
    {
        public static char Delimiter { get { return ';'; } }

        public KPMP_metabolite_results_readWriteOptions(string fileName)
        {
            this.File = Global_directory_class.Results_directory + fileName;
            this.Key_propertyNames = new string[] { "Molecular_formula", "Glom_correlation", "ReadWrite_hmdb_ids", "Is_glomerula_origin", "ReadWrite_associated_proteins", "ReadWrite_potential_metabolites" };
            this.Key_columnNames = this.Key_propertyNames;
            this.HeadlineDelimiters = new char[] { Global_class.Tab };
            this.LineDelimiters = new char[] { Global_class.Tab };
            this.File_has_headline = true;
            this.Report = ReadWrite_report_enum.Report_main;
        }
    }

    class KPMP_metabolite_class
    {
        KPMP_metabolite_line_class[] Metabolites { get; set; }
        string[] Bg_genes_in_upperCase { get; set; }

        public KPMP_metabolite_class()
        {
            this.Metabolites = new KPMP_metabolite_line_class[0];
            Bg_genes_in_upperCase = new string[0];
        }

        private void Add_KPMP_integration_term()
        {
            int met_length = this.Metabolites.Length;
            KPMP_metabolite_line_class met_line;
            for (int indexMet = 0; indexMet < met_length; indexMet++)
            {
                met_line = this.Metabolites[indexMet];
                met_line.KPMP_data_integration_term = KPMP_data_integration_class.Get_kpmp_integration_term_plus_patientID(met_line.Is_glomerula_origin.ToString(), met_line.PatientID);
            }
        }

        private void Add_associated_proteins_and_potential_metabolites_and_bg_genes()
        {
            HMDB_metabolite_class metabolites = new HMDB_metabolite_class();
            metabolites.Generate_by_reading_safed_file();

            this.Bg_genes_in_upperCase = metabolites.Get_all_ordered_distinct_accociated_genes_in_uppserCase();

            HMDB_metabolite_line_class reference_line;

            metabolites.Metabolites = metabolites.Metabolites.OrderBy(l => l.Chemical_formula).ToArray();

            Metabolites = Metabolites.OrderBy(l => l.Molecular_formula).ToArray();
            int this_length = this.Metabolites.Length;

            int reference_length = metabolites.Metabolites.Length;
            KPMP_metabolite_line_class metabolite_line;
            int indexRef = 0;
            int stringCompare = -2;

            List<string> associated_proteins = new List<string>();
            List<string> potential_metabolites = new List<string>();

            for (int indexThis = 0; indexThis < this_length; indexThis++)
            {
                metabolite_line = this.Metabolites[indexThis];
                associated_proteins.Clear();
                potential_metabolites.Clear();
                stringCompare = -2;
                while ((indexRef < reference_length) && (stringCompare <= 0))
                {
                    reference_line = metabolites.Metabolites[indexRef];
                    stringCompare = reference_line.Chemical_formula.CompareTo(metabolite_line.Molecular_formula);
                    if (stringCompare < 0)
                    {
                        indexRef++;
                    }
                    else if (stringCompare == 0)
                    {
                        associated_proteins.AddRange(reference_line.Protein_associations);
                        potential_metabolites.Add(reference_line.Primary_hmdb_accession);// + " (" + reference_line.Metabolite_name + ")");
                        indexRef++;
                    }
                }
                metabolite_line.Associated_proteins = associated_proteins.Distinct().OrderBy(l => l).ToArray();
                metabolite_line.Potential_metabolites = potential_metabolites.ToArray();
            }
        }

        private void Add_to_array(KPMP_metabolite_line_class[] add_metabolites)
        {
            int add_length = add_metabolites.Length;
            int this_length = this.Metabolites.Length;
            int new_length = add_length + this_length;
            KPMP_metabolite_line_class[] new_metabolites = new KPMP_metabolite_line_class[new_length];
            int indexNew = -1;
            for (int indexThis = 0; indexThis < this_length; indexThis++)
            {
                indexNew++;
                new_metabolites[indexNew] = this.Metabolites[indexThis];
            }
            for (int indexAdd = 0; indexAdd < add_length; indexAdd++)
            {
                indexNew++;
                new_metabolites[indexNew] = add_metabolites[indexAdd];
            }
            this.Metabolites = new_metabolites;
        }

        private string[] Get_all_participants()
        {
            Dictionary<string, bool> participant_exists_dict = new Dictionary<string, bool>();
            foreach (KPMP_metabolite_line_class metabolite_line in this.Metabolites)
            {
                if (!participant_exists_dict.ContainsKey(metabolite_line.PatientID))
                {
                    participant_exists_dict.Add(metabolite_line.PatientID, true);
                }
            }
            return participant_exists_dict.Keys.OrderBy(l => l).ToArray();
        }

        private void Add_avarage_correlation()
        {
            this.Metabolites = this.Metabolites.OrderBy(l => l.Molecular_formula).ToArray();
            int patients_count = Get_all_participants().Length;
            List<float> correlations = new List<float>();
            int metabolites_length = this.Metabolites.Length;
            KPMP_metabolite_line_class metabolite_line;
            KPMP_metabolite_line_class new_metabolite_line;
            List<KPMP_metabolite_line_class> add = new List<KPMP_metabolite_line_class>();
            for (int indexM = 0; indexM < metabolites_length; indexM++)
            {
                metabolite_line = this.Metabolites[indexM];
                if ((indexM == 0)
                    || (!metabolite_line.Molecular_formula.Equals(this.Metabolites[indexM - 1].Molecular_formula)))
                {
                    correlations.Clear();
                }
                correlations.Add(metabolite_line.Is_glomerula_origin);
                if ((indexM == metabolites_length - 1)
                    || (!metabolite_line.Molecular_formula.Equals(this.Metabolites[indexM + 1].Molecular_formula)))
                {
                    new_metabolite_line = metabolite_line.Deep_copy();
                    new_metabolite_line.PatientID = "allPatients";
                    while (correlations.Count < patients_count)
                    {
                        correlations.Add(0);
                    }
                    new_metabolite_line.Is_glomerula_origin = (int)Math.Ceiling(Math_class.Get_average(correlations.ToArray()));
                    add.Add(new_metabolite_line);
                    correlations.Clear();
                }
            }
            Add_to_array(add.ToArray());
        }

        public void Keep_top_x_correlations_per_patientID_glomerular_origin(int keep_top_x)
        {
            List<KPMP_metabolite_line_class> keep_metabolites = new List<KPMP_metabolite_line_class>();
            int metabolites_length = this.Metabolites.Length;
            KPMP_metabolite_line_class metabolite_line;
            this.Metabolites = this.Metabolites.OrderBy(l => l.PatientID).ThenBy(l => l.Is_glomerula_origin).ThenByDescending(l => l.Glom_correlation).ToArray();
            int kept_metabolites_of_current_condition = 0;
            for (int indexM = 0; indexM < metabolites_length; indexM++)
            {
                metabolite_line = this.Metabolites[indexM];
                if ((indexM == 0)
                    || (!metabolite_line.PatientID.Equals(this.Metabolites[indexM - 1].PatientID))
                    || (!metabolite_line.Is_glomerula_origin.Equals(this.Metabolites[indexM - 1].Is_glomerula_origin)))
                {
                    kept_metabolites_of_current_condition = 0;
                }
                if (kept_metabolites_of_current_condition < keep_top_x)
                {
                    kept_metabolites_of_current_condition++;
                    keep_metabolites.Add(metabolite_line);
                }
            }
            this.Metabolites = keep_metabolites.ToArray();
        }

        public void Keep_only_indicated_patientIDs(params string[] keep_patientIDs)
        {
            keep_patientIDs = keep_patientIDs.Distinct().OrderBy(l => l).ToArray();
            string keep_patientID;
            int keep_patientIDs_length = keep_patientIDs.Length;
            int indexPatientID = 0;
            int stringCompare = -2;
            int metabolites_length = this.Metabolites.Length;
            KPMP_metabolite_line_class metabolite_line;
            List<KPMP_metabolite_line_class> metabolite_list = new List<KPMP_metabolite_line_class>();
            for (int indexM = 0; indexM < metabolites_length; indexM++)
            {
                metabolite_line = this.Metabolites[indexM];
                stringCompare = -2;
                while ((indexPatientID < keep_patientIDs_length) && (stringCompare < 0))
                {
                    keep_patientID = keep_patientIDs[indexPatientID];
                    stringCompare = keep_patientID.CompareTo(metabolite_line.PatientID);
                    if (stringCompare < 0)
                    {
                        indexPatientID++;
                    }
                    else if (stringCompare == 0)
                    {
                        metabolite_list.Add(metabolite_line);
                    }
                }
            }
            this.Metabolites = metabolite_list.ToArray();
        }

        private void Remove_brackets_and_quotation_marks()
        {
            int metabolites_length = this.Metabolites.Length;
            KPMP_metabolite_line_class metabolite_line;
            int hmdbs_length;
            for (int indexM = 0; indexM < metabolites_length; indexM++)
            {
                metabolite_line = this.Metabolites[indexM];
                hmdbs_length = metabolite_line.Hmdb_ids.Length;
                for (int indexHmbdb = 0; indexHmbdb < hmdbs_length; indexHmbdb++)
                {
                    metabolite_line.Hmdb_ids[indexHmbdb] = Text2_class.Remove_indicated_characters_from_end_and_beginning_of_text(metabolite_line.Hmdb_ids[indexHmbdb], new char[] { '[', ']', '\'' });
                }
            }
        }

        public void Generate(params string[] patientIDs)
        {
            Read_and_add_to_array(patientIDs);
            Add_avarage_correlation();
            Remove_brackets_and_quotation_marks();
            //  Add_associated_proteins_and_potential_metabolites_and_bg_genes();
        }

        public DE_class Generate_de_instance_with_associated_proteins_based_total_gene_counts(string firstName_for_de)
        {
            Fill_de_line_class fill_de_line;
            List<Fill_de_line_class> fill_de_list = new List<Fill_de_line_class>();
            foreach (KPMP_metabolite_line_class line in this.Metabolites)
            {
                if (line.Is_glomerula_origin == 1)
                {
                    Console.WriteLine();
                }
                if (line.Is_glomerula_origin == 0)
                {
                    Console.WriteLine();
                }
                foreach (string protein in line.Associated_proteins)
                {
                    if (!String.IsNullOrEmpty(protein))
                    {
                        fill_de_line = new Fill_de_line_class();
                        fill_de_line.Names_for_de = new string[] { line.Is_glomerula_origin.ToString() + "$" + KPMP_data_integration_class.Get_kpmp_integration_term_plus_patientID("Gloms", line.PatientID) };
                        fill_de_line.Symbols_for_de = new string[] { (string)protein.Clone() };
                        fill_de_line.Value_for_de = 1;
                        fill_de_list.Add(fill_de_line);
                    }
                }
            }
            fill_de_list = fill_de_list.OrderBy(l => l.Names_for_de[0]).ThenBy(l => l.Symbols_for_de[0]).ThenByDescending(l => l.Value_for_de).ToList();
            List<Fill_de_line_class> final_fill_de_list = new List<Fill_de_line_class>();
            int fill_count = fill_de_list.Count;
            for (int indexFill = 0; indexFill < fill_count; indexFill++)
            {
                fill_de_line = fill_de_list[indexFill];
                if ((indexFill == fill_count - 1)
                    || (!fill_de_line.Names_for_de[0].Equals(fill_de_list[indexFill + 1].Names_for_de[0]))
                    || (!fill_de_line.Symbols_for_de[0].Equals(fill_de_list[indexFill + 1].Symbols_for_de[0])))
                {
                    final_fill_de_list.Add(fill_de_line);
                }
                else
                {
                    fill_de_list[indexFill + 1].Value_for_de += fill_de_line.Value_for_de;
                }
            }

            final_fill_de_list = final_fill_de_list.OrderBy(l => l.Symbols_for_de[0]).ThenByDescending(l => l.Names_for_de[0]).ToList();
            int final_fill_de_list_count = final_fill_de_list.Count;
            Fill_de_line_class new_fill_de_line;
            List<Fill_de_line_class> new_fill_de_list = new List<Fill_de_line_class>();
            double non_glomerular_value = -1;
            double glomerular_value = -1;
            for (int indexF = 0; indexF < final_fill_de_list_count; indexF++)
            {
                fill_de_line = final_fill_de_list[indexF];
                if ((indexF == 0) || (!fill_de_line.Symbols_for_de[0].Equals(final_fill_de_list[indexF - 1].Symbols_for_de[0])))
                {
                    non_glomerular_value = -1;
                    glomerular_value = -1;
                }
                if (fill_de_line.Names_for_de[0].IndexOf("0") == 0)
                {
                    if (non_glomerular_value != -1) { throw new Exception(); }
                    non_glomerular_value = fill_de_line.Value_for_de;
                }
                if (fill_de_line.Names_for_de[0].IndexOf("1") == 0)
                {
                    if (glomerular_value != -1) { throw new Exception(); }
                    glomerular_value = fill_de_line.Value_for_de;
                }
                if ((indexF == final_fill_de_list_count - 1) || (!fill_de_line.Symbols_for_de[0].Equals(final_fill_de_list[indexF + 1].Symbols_for_de[0])))
                {
                    new_fill_de_line = new Fill_de_line_class();
                    new_fill_de_line.Symbols_for_de = new string[] { (string)fill_de_line.Symbols_for_de[0].Clone() };
                    new_fill_de_line.Names_for_de = new string[] { firstName_for_de + "$" + "Glom - POD" + "$" + "Diff protein association" + "$" + "Podocyte_glomerulus_allPatients" };
                    new_fill_de_line.Value_for_de = (float)(glomerular_value + 1) / (float)(non_glomerular_value + 1);
                    new_fill_de_list.Add(new_fill_de_line);

                    new_fill_de_line = new_fill_de_line.Deep_copy();
                    new_fill_de_line.Names_for_de = new string[] { firstName_for_de + "$" + "Glom - Mesangial" + "$" + "Diff protein association" + "$" + "Mesangial_allPatients" };
                    new_fill_de_list.Add(new_fill_de_line);

                    new_fill_de_line = new_fill_de_line.Deep_copy();
                    new_fill_de_line.Names_for_de = new string[] { firstName_for_de + "$" + "Glom - EC" + "$" + "Diff protein association" + "$" + "Endothelial glomerula_allPatients" };
                    new_fill_de_list.Add(new_fill_de_line);

                    new_fill_de_line = new Fill_de_line_class();
                    new_fill_de_line.Symbols_for_de = new string[] { (string)fill_de_line.Symbols_for_de[0].Clone() };
                    new_fill_de_line.Names_for_de = new string[] { firstName_for_de + "$" + "TI - PT" + "$" + "Diff protein association" + "$" + "Proximal_tubule_allPatients" };
                    new_fill_de_line.Value_for_de = (float)(non_glomerular_value + 1) / (float)(glomerular_value + 1);
                    new_fill_de_list.Add(new_fill_de_line);

                    new_fill_de_line = new_fill_de_line.Deep_copy();
                    new_fill_de_line.Names_for_de = new string[] { firstName_for_de + "$" + "TI - Interstitium" + "$" + "Diff protein association" + "$" + "Interstitium_allPatients" };
                    new_fill_de_list.Add(new_fill_de_line);
                }
            }
            final_fill_de_list.AddRange(new_fill_de_list);
            DE_class de = new DE_class();
            de.Fill_with_data_alternatively(new_fill_de_list.ToArray());
            return de;
        }

        public DE_class Generate_de_instance_with_metabolite_hmdbIDs_filled_with_highest_rank_for_associated_molecular_formula(string firstName_for_de)
        {
            Fill_de_line_class fill_de_line;
            List<Fill_de_line_class> fill_de_list = new List<Fill_de_line_class>();
            int metabolite_length = this.Metabolites.Length;
            KPMP_metabolite_line_class metabolite_line;
            this.Metabolites = this.Metabolites.OrderBy(l => l.Is_glomerula_origin).ThenBy(l => l.PatientID).ThenByDescending(l => l.Glom_correlation).ToArray();
            int current_molecular_formula_rank = 0;
            for (int indexM = 0; indexM < metabolite_length; indexM++)
            {
                metabolite_line = this.Metabolites[indexM];
                if ((indexM == 0)
                    || (!metabolite_line.Is_glomerula_origin.Equals(this.Metabolites[indexM - 1].Is_glomerula_origin))
                    || (!metabolite_line.PatientID.Equals(this.Metabolites[indexM - 1].PatientID)))
                {
                    current_molecular_formula_rank = 0;
                }
                current_molecular_formula_rank++;
                foreach (string hmdbID in metabolite_line.Hmdb_ids)
                {
                    if (!String.IsNullOrEmpty(hmdbID))
                    {
                        fill_de_line = new Fill_de_line_class();
                        fill_de_line.Names_for_de = new string[] { metabolite_line.Is_glomerula_origin.ToString() + "$" + metabolite_line.PatientID };
                        fill_de_line.Symbols_for_de = new string[] { (string)hmdbID.Clone() };
                        fill_de_line.Value_for_de = current_molecular_formula_rank;
                        fill_de_list.Add(fill_de_line);
                    }
                }
            }
            fill_de_list = fill_de_list.OrderBy(l => l.Names_for_de[0]).ThenBy(l => l.Symbols_for_de[0]).ThenBy(l => l.Value_for_de).ToList();
            List<Fill_de_line_class> final_fill_de_list = new List<Fill_de_line_class>();
            int fill_count = fill_de_list.Count;
            for (int indexFill = 0; indexFill < fill_count; indexFill++)
            {
                fill_de_line = fill_de_list[indexFill];
                if ((indexFill == fill_count - 1)
                    || (!fill_de_line.Names_for_de[0].Equals(fill_de_list[indexFill + 1].Names_for_de[0]))
                    || (!fill_de_line.Symbols_for_de[0].Equals(fill_de_list[indexFill + 1].Symbols_for_de[0])))
                {
                    final_fill_de_list.Add(fill_de_line);
                }
                else
                {
                }
            }

            DE_class de = new DE_class();
            de.Fill_with_data_alternatively(final_fill_de_list.ToArray());
            return de;
        }

 
        public KPMP_standardized_dataset_class Generate_kpmp_standardized_dataset_instance_for_allPatients(KPMP_metabolite_id_type_enum metabolite_name_type)
        {
            int metabolites_length = this.Metabolites.Length;
            KPMP_metabolite_line_class metabolite_line;
            KPMP_standardized_dataset_line_class standardized_dataset_line;
            List<KPMP_standardized_dataset_line_class> standardized_dataset_list = new List<KPMP_standardized_dataset_line_class>();
            string[] metabolite_ids = new string[0];
            for (int indexM=0; indexM<metabolites_length; indexM++)
            {
                metabolite_line = this.Metabolites[indexM];
                standardized_dataset_line = new KPMP_standardized_dataset_line_class();
                switch (metabolite_name_type)
                {
                    case KPMP_metabolite_id_type_enum.Hmdb:
                        metabolite_ids = metabolite_line.Hmdb_ids;
                        break;
                    case KPMP_metabolite_id_type_enum.Kegg_id:
                        metabolite_ids = metabolite_line.Kegg_ids;
                        break;
                    default:
                        throw new Exception();
                }

                if (metabolite_line.PatientID.Equals("allPatients"))
                {
                    foreach (string metabolite_id in metabolite_ids)
                    {
                        if (!metabolite_id.Equals("NA"))
                        {
                            standardized_dataset_line = new KPMP_standardized_dataset_line_class();
                            if (metabolite_line.Is_glomerula_origin == 1)
                            {
                                standardized_dataset_line.Cell_segment = "Glom";
                            }
                            else if (metabolite_line.Is_glomerula_origin == 0)
                            {
                                standardized_dataset_line.Cell_segment = "TI";
                            }
                            else { throw new Exception(); }
                            standardized_dataset_line.Dataset = (string)KPMP_dataset_name_class.Spatial_metabolomics.Clone();
                            standardized_dataset_line.Gene_symbol = (string)metabolite_id.Clone();
                            standardized_dataset_line.KPMP_data_integration_term = "";
                            standardized_dataset_line.PatientId = (string)metabolite_line.PatientID.Clone();
                            standardized_dataset_line.Value_1st = metabolite_line.Glom_correlation;
                            standardized_dataset_line.Value_type_1st = KPMP_value_type_enum.Correlation;
                            standardized_dataset_line.Value_type_2nd = KPMP_value_type_enum.Log2_ratioavg;
                            standardized_dataset_line.Value_2nd = 1;
                            standardized_dataset_list.Add(standardized_dataset_line);
                        }
                    }
                }
            }
            KPMP_standardized_dataset_class standardized_kpmp = new KPMP_standardized_dataset_class();
            standardized_kpmp.Add_deep_copy_to_array(standardized_dataset_list.ToArray());
            return standardized_kpmp;
        }

        public string[] Get_deep_copy_of_bg_gene_in_upper_case()
        {
            return Array_class.Deep_copy_string_array(this.Bg_genes_in_upperCase);
        }

        public void Read_and_add_to_array(string[] patientIDs)
        {
            patientIDs = patientIDs.Distinct().ToArray();
            int patientIDs_length = patientIDs.Length;
            string patientID;
            for (int indexP = 0; indexP < patientIDs_length; indexP++)
            {
                patientID = patientIDs[indexP];
                string fileName = "Metabolomics_2018December15_" + patientID + ".csv";
                KPMP_metabolite_inputReadOptions inputReadOptions = new KPMP_metabolite_inputReadOptions(fileName);
                KPMP_metabolite_line_class[] add_metabolites = ReadWriteClass.ReadRawData_and_FillArray<KPMP_metabolite_line_class>(inputReadOptions);
                foreach (KPMP_metabolite_line_class add_metabolite in add_metabolites)
                {
                    add_metabolite.PatientID = (string)patientID.Clone();
                }
                Add_to_array(add_metabolites);
            }
        }

        public void Write_into_results_directory(string fileName)
        {
            KPMP_metabolite_results_readWriteOptions readWriteOptions = new KPMP_metabolite_results_readWriteOptions(fileName);
            ReadWriteClass.WriteData(this.Metabolites, readWriteOptions);
        }
    }

}
