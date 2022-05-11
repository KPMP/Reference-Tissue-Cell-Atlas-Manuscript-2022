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
using System.Linq;
using Enumerations;
using ReadWrite;


namespace MBCO_data
{
    class Custom_data_line_class : IAdd_to_data
    {
        public string Dataset { get; set; }
        public string KPMP_data_integration_term { get; set; }
        public string Cell_segment { get; set; }
        public string Gene_symbol { get; set; }
        public double Value { get; set; }

        public string NCBI_official_symbol_for_data { get { return Gene_symbol; } }
        public string SampleName_for_data { get { return KPMP_data_integration_term + "$" + Dataset + "$" + Cell_segment; } }
        public double Value_for_data { get { return Value; } }

        public Custom_data_line_class Deep_copy()
        {
            Custom_data_line_class copy = (Custom_data_line_class)this.MemberwiseClone();
            copy.Dataset = (string)this.Dataset.Clone();
            copy.KPMP_data_integration_term = (string)this.KPMP_data_integration_term.Clone();
            copy.Cell_segment = (string)this.Cell_segment.Clone();
            copy.Gene_symbol = (string)this.Gene_symbol.Clone();
            return copy;
        }
    }

    class Custom_data_readWriteOptions_class : ReadWriteOptions_base
    {
        public Custom_data_readWriteOptions_class(string file_name)
        {
            string directory = Global_directory_class.Marker_genes_proteins_for_MBCO_directory;
            this.File = directory + file_name;
            Key_propertyNames = new string[] { "Dataset", "KPMP_data_integration_term", "Cell_segment","Gene_symbol","Value" };
            Key_columnNames = new string[] { "Dataset", "KPMP_data_integration_term", "Cell_segment", "Gene_symbol", "Value_2nd" };
            HeadlineDelimiters = new char[] { Global_class.Tab };
            LineDelimiters = new char[] { Global_class.Tab };
            File_has_headline = true;
            Report = ReadWrite_report_enum.Report_nothing;
        }
    }

    class Custom_data_class
    {
        public Custom_data_line_class[] Custom_data { get; set; }

        public Custom_data_class()
        {
            Custom_data = new Custom_data_line_class[0];
        }



        private void Check_for_duplicates()
        {
            int custom_data_length = Custom_data.Length;
            Custom_data_line_class custom_data_line;
            Custom_data_line_class previous_custom_data_line;
            this.Custom_data = this.Custom_data.OrderBy(l => l.Dataset).ThenBy(l=>l.Cell_segment).ThenBy(l => l.Gene_symbol).ToArray();
            for (int indexC = 1; indexC < custom_data_length; indexC++)
            {
                custom_data_line = this.Custom_data[indexC];
                previous_custom_data_line = this.Custom_data[indexC - 1];
                if (   (custom_data_line.Dataset.Equals(previous_custom_data_line.Dataset))
                    && (custom_data_line.Cell_segment.Equals(previous_custom_data_line.Cell_segment))
                    && (custom_data_line.Gene_symbol.Equals(previous_custom_data_line.Gene_symbol)))
                {
                   throw new Exception();
                }
            }
        }

        public void Add_to_array(Custom_data_line_class[] add_custom_data)
        {
            int this_length = this.Custom_data.Length;
            int add_length = add_custom_data.Length;
            int new_length = this_length + add_length;
            Custom_data_line_class[] new_custom_data = new Custom_data_line_class[new_length];
            int indexNew = -1;
            for (int indexThis = 0; indexThis < this_length; indexThis++)
            {
                indexNew++;
                new_custom_data[indexNew] = this.Custom_data[indexThis];
            }
            for (int indexAdd = 0; indexAdd < add_length; indexAdd++)
            {
                indexNew++;
                new_custom_data[indexNew] = add_custom_data[indexAdd];
            }
            this.Custom_data = new_custom_data;
        }

        public void Clear_custom_data()
        {
            this.Custom_data = new Custom_data_line_class[0];
        }

        public void Generate_custom_data_instance(Custom_data_readWriteOptions_class readWriteOptions)
        {
            Read(readWriteOptions);
            Check_for_duplicates();
        }

        public Data_class Generate_new_data_instance()
        {
            Check_for_duplicates();
            Data_class data = new Data_class();
            data.Add_to_data_instance(this.Custom_data);
            data.Set_all_ncbi_official_gene_symbols_to_upper_case();
            return data;
        }

        private void Read(Custom_data_readWriteOptions_class readWriteOptions)
        {
            this.Custom_data = ReadWriteClass.ReadRawData_and_FillArray<Custom_data_line_class>(readWriteOptions);
        }
    }
}
