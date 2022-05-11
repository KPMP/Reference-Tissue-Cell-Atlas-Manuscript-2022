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

using ReadWrite;
using Enumerations;

namespace KPMP
{
    class KPMP_integration_paper_metadata_line_class
    {
        public string Figure { get; set; }
        public string[] Libraries { get; set; }
        public string Tis { get; set; }
        public string Dataset { get; set; }
        public string Analysis { get; set; }
        public string Patient_id { get; set; }
        public string[] Tissue_types { get; set; }
        public string Tissue_collection { get; set; }
        public string Tissue_interrogation_site { get; set; }
        public string Pubmed_id { get; set; }

        public string Analysis_dataset { get { return Analysis + " - " + Dataset; } }

        public string ReadWrite_tissue_types
        {
            get { return ReadWriteClass.Get_writeLine_from_array(this.Tissue_types, KPMP_integration_paper_metadata_input_readWriteOptions_class.Delimiter); }
            set { this.Tissue_types = ReadWriteClass.Get_array_from_readLine<string>(value, KPMP_integration_paper_metadata_input_readWriteOptions_class.Delimiter); }
        }

        public string ReadWrite_libraries
        {
            get { return ReadWriteClass.Get_writeLine_from_array(this.Libraries, KPMP_integration_paper_metadata_input_readWriteOptions_class.Delimiter); }
            set { this.Libraries = ReadWriteClass.Get_array_from_readLine<string>(value, KPMP_integration_paper_metadata_input_readWriteOptions_class.Delimiter); }
        }

        public KPMP_integration_paper_metadata_line_class()
        {
            Figure = "";
            Libraries = new string[0];
            Dataset = "";
            Analysis = "";
            Patient_id = "";
            Tissue_types = new string[] { "" };
            Tissue_collection = "";
            Tissue_interrogation_site = "";
            Pubmed_id = "";
        }

        public KPMP_integration_paper_metadata_line_class Deep_copy()
        {
            KPMP_integration_paper_metadata_line_class copy = (KPMP_integration_paper_metadata_line_class)this.MemberwiseClone();
            copy.Libraries = Array_class.Deep_copy_string_array(this.Libraries);
            copy.Dataset = (string)this.Dataset.Clone();
            copy.Analysis = (string)this.Analysis.Clone();
            copy.Patient_id = (string)this.Patient_id.Clone();
            copy.Tissue_interrogation_site = (string)this.Tissue_interrogation_site.Clone();
            copy.Tissue_collection = (string)this.Tissue_collection.Clone();
            copy.Tissue_types = Array_class.Deep_copy_string_array(this.Tissue_types);
            copy.Pubmed_id = (string)this.Pubmed_id.Clone();
            copy.Figure = (string)this.Figure.Clone();
            return copy;
        }
    }

    class KPMP_integration_paper_metadata_input_readWriteOptions_class : ReadWriteOptions_base
    {
        public static char Delimiter { get { return ';'; } }

        public KPMP_integration_paper_metadata_input_readWriteOptions_class(string fileName)
        {
            string complete_directory = Global_directory_class.Experimental_data_directory + "Sample_metadata//";
            this.File = complete_directory + fileName;
            this.Key_propertyNames = new string[] { "Figure", "ReadWrite_libraries", "Dataset", "Analysis", "Patient_id", "ReadWrite_tissue_types", "Tissue_collection", "Tissue_interrogation_site", "Pubmed_id" };
            this.Key_columnNames = this.Key_propertyNames;
            this.HeadlineDelimiters = new char[] { Global_class.Tab };
            this.LineDelimiters = new char[] { Global_class.Tab };
            this.File_has_headline = true;
            this.Report = ReadWrite_report_enum.Report_main;
        }
    }

    class KPMP_integration_paper_metadata_results_readWriteOptions_class : ReadWriteOptions_base
    {
        public static char Delimiter { get { return ';'; } }

        public KPMP_integration_paper_metadata_results_readWriteOptions_class(string subdirectory, string fileName)
        {
            string complete_directory = Global_directory_class.Results_directory + subdirectory;
            ReadWriteClass.Create_directory_if_it_does_not_exist(complete_directory);
            this.File = complete_directory + fileName;
            this.Key_propertyNames = new string[] { "Figure", "ReadWrite_libraries", "Dataset", "Analysis", "Patient_id", "Tis","ReadWrite_tissue_types", "Tissue_collection", "Tissue_interrogation_site", "Pubmed_id" };
            this.Key_columnNames = this.Key_propertyNames;
            this.HeadlineDelimiters = new char[] { Global_class.Tab };
            this.LineDelimiters = new char[] { Global_class.Tab };
            this.File_has_headline = true;
            this.Report = ReadWrite_report_enum.Report_main;
        }
    }

    class KPMP_integration_paper_metadata_class
    {
        public KPMP_integration_paper_metadata_line_class[] Documentations { get; set; }

        public KPMP_integration_paper_metadata_class()
        {
            this.Documentations = new KPMP_integration_paper_metadata_line_class[0];
        }

        public void Add_to_array(KPMP_integration_paper_metadata_line_class[] add_documentations)
        {
            int this_documentations_length = this.Documentations.Length;
            int add_documentations_length = add_documentations.Length;
            int new_documentations_length = this_documentations_length + add_documentations_length;
            KPMP_integration_paper_metadata_line_class[] new_documentations = new KPMP_integration_paper_metadata_line_class[new_documentations_length];
            int indexNew = -1;
            for (int indexThis = 0; indexThis < this_documentations_length; indexThis++)
            {
                indexNew++;
                new_documentations[indexNew] = this.Documentations[indexThis];
            }
            for (int indexAdd = 0; indexAdd < add_documentations_length; indexAdd++)
            {
                indexNew++;
                new_documentations[indexNew] = add_documentations[indexAdd];
            }
            this.Documentations = new_documentations;
        }

        public void Add_analysis_to_all_lines(string analysis)
        {
            foreach (KPMP_integration_paper_metadata_line_class patient_line in this.Documentations)
            {
                patient_line.Analysis = (string)analysis.Clone();
            }
        }
 
        public void Add_other(KPMP_integration_paper_metadata_class other)
        {
            Add_to_array(other.Documentations);
        }

        public void Write_into_results_directory(string subdirectory, string fileName)
        {
            KPMP_integration_paper_metadata_results_readWriteOptions_class readWriteOptions = new KPMP_integration_paper_metadata_results_readWriteOptions_class(subdirectory, fileName);
            ReadWriteClass.WriteData(this.Documentations, readWriteOptions);
        }

        public KPMP_integration_paper_metadata_class Deep_copy()
        {
            KPMP_integration_paper_metadata_class copy = (KPMP_integration_paper_metadata_class)this.MemberwiseClone();
            int data_length = this.Documentations.Length;
            copy.Documentations = new KPMP_integration_paper_metadata_line_class[data_length];
            for (int indexD = 0; indexD < data_length; indexD++)
            {
                copy.Documentations[indexD] = this.Documentations[indexD].Deep_copy();
            }
            return copy;
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

}
