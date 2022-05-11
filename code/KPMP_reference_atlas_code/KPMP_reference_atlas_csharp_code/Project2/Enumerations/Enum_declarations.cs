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

namespace Enumerations
{
    class Global_class
    {
        #region Const
        private const Input_data_enum input_data = Input_data_enum.Kidney; //Check, if right organism and if gene ontology is switched to organism
        private const Organism_enum organism = Organism_enum.Homo_sapiens;
        private const string empty_entry = "E_m_p_t_y";  //check enums, Empty has to be the same!!
        private const char tab = '\t';
        private const char comma = ',';
        private const string space_text = "S_p_a_c_e";
        private const int inf_number_abs = 99999999;
        public const bool Check_ordering = true;
        public const int Csharp_script_number = 0;
        public static string Csharp_script_label = "_script" + Csharp_script_number;
        #endregion

        public static Input_data_enum Input_data
        {
            get
            {
                Input_data_enum Input_data = input_data;
                return Input_data;
            }
        }

        public static Organism_enum Organism
        {
            get { return organism; }
        }

        public static int Inf_number_abs
        { get { return inf_number_abs; } }

        public static int Inf_number_pos
        { get { return Inf_number_abs; } }

        public static int Inf_number_neg
        { get { return -Inf_number_abs; } }


        public static string Empty_entry
        { get { return empty_entry; } }

        public static char Comma
        { get { return comma; } }

        public static char Tab
        { get { return tab; } }

        public static string Space_text
        { get { return space_text; } }
    }

    class Global_directory_class
    {
        private const string hard_drive = "D:\\";
        private const string experimental_data_subdirectory = "Experimental_data\\";
        private const string additional_sc_sn_subdirectory = "Additional_singleCellNucleus_datasets\\";
        private const string mbco_dataset_subdirectory = "MBCO_datasets\\";
        private const string go_dataset_subdirectory = "GeneOntology_datasets\\";
        private const string results_subdirectory = "Results\\";
        private const string mbco_enrichment_results_subdirectory = "MBCO_enrichment_results\\";
        private const string mbco_scp_networks_subdirectory = "MBCO_scp_networks\\";
        private const string marker_genes_proteins_subdirectory = "Marker_genes_and_proteins\\";
        private const string marker_genes_proteins_for_MBCO_subdirectory = "MBCO_enrichment_input\\";
        private const string downloaded_datasets_subdirectory = "Downloaded_datasets\\";
        private const string metabolic_enrichment_results_directory = "Metabolic_activities\\";
        private const string sample_metadata_additional_datasets_subdirectory = "Sample_metadata_additional_datasets\\";

        public static string Major_directory
        { 
            get 
            { 
                string current_directory = System.Environment.CurrentDirectory + "\\";
                int indexKPMP_directory = current_directory.ToUpper().IndexOf("KPMP_reference_atlas_csharp_code".ToUpper());
                if (indexKPMP_directory == -1) { throw new Exception(); }//Directory that contains c# script needs to match name above
                return current_directory.Substring(0, indexKPMP_directory);
            }
            
        }
        public static string Get_report_fileName()
        {
            return Global_directory_class.Results_directory + "Report.txt";
        }

        public static string Mbco_dataset_directory
        { get { return Major_directory + mbco_dataset_subdirectory; } }
        public static string Go_dataset_directory
        { get { return Major_directory + go_dataset_subdirectory; } }

        public static string Gene_library_directory
        { get { return Major_directory + mbco_dataset_subdirectory; } }

        public static string GeneNameDatabases_directory
        { get { return Major_directory + downloaded_datasets_subdirectory; } }

        public static string Complete_human_mbco_association_v11_fileName
        { get { return Mbco_dataset_directory + "MBCO_v1.1_gene-SCP_associations_human.txt"; } }

        public static string Complete_mbco_inferred_scp_relationships_v11_fileName
        {
            get { return Mbco_dataset_directory + "MBCO_v1.1_inferred_SCP_relationships.txt"; }
        }

        public static string Experimental_data_directory
        { get { return Major_directory + experimental_data_subdirectory; } }

        public static string Seurat_singleCell_cluster_directory
        { get { return Experimental_data_directory + "SingleCellCluster\\"; } }
        public static string PostHoc_power_seurat_results_directory
        { get { return Results_directory + "PostHocPower_seurat\\"; } }

        public static string Seurat_integratedCluster_avgExpression_directory
        { get { return Experimental_data_directory + "IntegratedCluster_AvgExpression_RNAcounts\\"; } }
        public static string Seurat_singleCell_data_directory
        {
            get
            {
                return Experimental_data_directory + "SingleCell_datasets\\";
            }
        }
        public static string Sample_metadata_additional_datasets
        { get { return Experimental_data_directory + sample_metadata_additional_datasets_subdirectory; } }
        public static string SingleCellNucleus_additional_datasets
        { get { return Experimental_data_directory + additional_sc_sn_subdirectory; } }
        public static string Results_directory
        { get { return Major_directory + results_subdirectory; } }
        public static string MBCO_enrichment_results_directory
        { get { return Results_directory + mbco_enrichment_results_subdirectory; } }
        public static string MBCO_scp_networks_results_directory
        { get { return MBCO_enrichment_results_directory + mbco_scp_networks_subdirectory; } }
        public static string Metabolic_activities_results_directory
        { get { return Results_directory + metabolic_enrichment_results_directory; } }
        public static string Marker_genes_proteins_directory
        { get { return Results_directory + marker_genes_proteins_subdirectory; } }
        public static string Marker_genes_proteins_for_MBCO_directory
        { get { return Marker_genes_proteins_directory + marker_genes_proteins_for_MBCO_subdirectory;} }
        public static string MBCO_results_directory
        { get { return Results_directory + mbco_enrichment_results_subdirectory; } }
    }
    public enum Input_data_enum { E_m_p_t_y, Kidney }
    public enum Organism_enum { E_m_p_t_y = 0, Homo_sapiens = 9606}

    public enum ReadWrite_report_enum { E_m_p_t_y = 0, Report_main = 1, Report_everything = 2, Report_nothing }
    public enum Ontology_type_enum
    {
        E_m_p_t_y, Metabolic_energy_generation_pathways, Mbco_level3
    }
    public enum Data_value_signs_of_interest_enum { E_m_p_t_y, Upregulated, Downregulated, Combined }

    class Deep_copy_class
    {
        public static string[] Deep_copy_string_array(string[] input)
        {
            int length = input.Length;
            string[] copy = new string[length];
            for (int indexS = 0; indexS < length; indexS++)
            {
                copy[indexS] = (string)input[indexS].Clone();
            }
            return copy;
        }
    }
    class Text2_class
    {
        public static string Remove_space_comma_semicolon_colon_underline_from_end_and_beginning_of_text(string text)
        {
            int text_length = text.Length;
            bool space_comma_semicolon_colon_at_beginning = true;
            bool space_comma_semicolon_colon_at_end = true;
            while (((space_comma_semicolon_colon_at_beginning) || (space_comma_semicolon_colon_at_end)) && ((!String.IsNullOrEmpty(text)) && (text_length >= 2)))
            {
                text_length = text.Length;
                space_comma_semicolon_colon_at_beginning = text[0].Equals(' ') || (text[0].Equals(',')) || (text[0].Equals('_')) || (text[0].Equals(';')) || (text[0].Equals(':'));
                space_comma_semicolon_colon_at_end = (text[text_length - 1].Equals(' ')) || (text[text_length - 1].Equals(',')) || (text[text_length - 1].Equals('_')) || (text[text_length - 1].Equals(';')) || (text[text_length - 1].Equals(':'));
                if (space_comma_semicolon_colon_at_beginning && space_comma_semicolon_colon_at_end)
                {
                    text = text.Substring(1, text_length - 2);
                }
                else if (space_comma_semicolon_colon_at_beginning)
                {
                    text = text.Substring(1, text_length - 1);
                }
                else if (space_comma_semicolon_colon_at_end)
                {
                    text = text.Substring(0, text_length - 1);
                }
            }
            return text;
        }

        public static string Remove_indicated_characters_from_end_and_beginning_of_text(string text, char[] remove_chars)
        {
            int text_length = text.Length;
            bool removeChar_at_beginning = true;
            bool removeChar_at_end = true;
            int removeChars_length = remove_chars.Length;
            while (((removeChar_at_beginning) || (removeChar_at_end)) && ((!String.IsNullOrEmpty(text)) && (text_length >= 2)))
            {
                text_length = text.Length;
                removeChar_at_beginning = text[0].Equals(' ') || (text[0].Equals(',')) || (text[0].Equals('_')) || (text[0].Equals(';')) || (text[0].Equals(':'));
                removeChar_at_end = false;
                for (int indexRemove = 0; indexRemove < removeChars_length; indexRemove++)
                {
                    if (text[text_length - 1].Equals(remove_chars[indexRemove])) { removeChar_at_end = true; }
                }
                for (int indexRemove = 0; indexRemove < removeChars_length; indexRemove++)
                {
                    if (text[0].Equals(remove_chars[indexRemove])) { removeChar_at_beginning = true; }
                }
                if (removeChar_at_beginning && removeChar_at_end)
                {
                    text = text.Substring(1, text_length - 2);
                }
                else if (removeChar_at_beginning)
                {
                    text = text.Substring(1, text_length - 1);
                }
                else if (removeChar_at_end)
                {
                    text = text.Substring(0, text_length - 1);
                }
            }
            return text;
        }
    }

    class Array_class
    {
        public static bool Array_order_dependent_equal<T>(T[] array1, T[] array2)
        {
            bool equal = true;
            int array1_length = array1.Length;
            int array2_length = array2.Length;
            if (array1_length != array2_length)
            {
                equal = false;
            }
            else
            {
                for (int indexA = 0; indexA < array1_length; indexA++)
                {
                    if (!array1[indexA].Equals(array2[indexA]))
                    {
                        equal = false;
                        break;
                    }
                }
            }
            return equal;
        }

        public static string[] Deep_copy_string_array(string[] array)
        {
            int array_length = array.Length;
            string[] copy = new string[0];
            if (array_length > 0)
            {
                copy = new string[array_length];
                for (int indexA = 0; indexA < array_length; indexA++)
                {
                    copy[indexA] = (string)array[indexA].Clone();
                }
            }
            return copy;
        }

        public static bool Equal_arrays<T>(T[] array1, T[] array2) where T : IComparable<T>
        {
            int array1_length = array1.Length;
            int array2_length = array2.Length;
            bool equal = true;
            if (array1_length != array2_length)
            {
                equal = false;
            }
            else
            {
                for (int indexA = 0; indexA < array1_length; indexA++)
                {
                    if (!array1[indexA].Equals(array2[indexA]))
                    {
                        equal = false;
                        break;
                    }
                }
            }
            return equal;
        }

        public static T[] Deep_copy_array<T>(T[] array)
        {
            T[] copy = new T[0];
            if (typeof(T) == typeof(string))
            {
                throw new Exception();
            }
            else if (array != null)
            {
                int array_length = array.Length;
                copy = new T[array_length];
                for (int indexA = 0; indexA < array_length; indexA++)
                {
                    copy[indexA] = (T)array[indexA];
                }
            }
            else
            {
                copy = array;
            }
            return copy;
        }
    }
}
