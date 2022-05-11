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

using System.Collections.Generic;
using System.Linq;
using ReadWrite;
using Enumerations;

namespace Gene_databases
{
    class HMDB_metabolite_self_readWriteOptions_class : ReadWriteOptions_base
    {
        public static char Delimiter { get { return '@'; } }

        public HMDB_metabolite_self_readWriteOptions_class()
        {
            File = Global_directory_class.GeneNameDatabases_directory + "Self\\Slightly modified Hmdb metabolites.jns";
            Key_propertyNames = new string[] { "Metabolite_name", "Primary_hmdb_accession", "Chemical_formula", "Organism", "ReadWrite_synonyms", "ReadWrite_secondary_hmdb_accesions", "ReadWrite_protein_associations" };
            Key_columnNames = Key_propertyNames;
            LineDelimiters = new char[] { Global_class.Tab };
            HeadlineDelimiters = new char[] { Global_class.Tab };
            File_has_headline = true;
            Report = ReadWrite_report_enum.Report_main;
        }
    }

    class HMDB_metabolite_line_class
    {
        public string Primary_hmdb_accession { get; set; }
        public string Metabolite_name { get; set; }
        public string Chemical_formula { get; set; }
        public Organism_enum Organism { get; set; }
        public string[] Secondary_hmdb_accessions { get; set; }
        public string[] Synonyms { get; set; }
        public string[] Protein_associations { get; set; }

        public string ReadWrite_synonyms
        {
            get { return ReadWriteClass.Get_writeLine_from_array(Synonyms, HMDB_metabolite_self_readWriteOptions_class.Delimiter); }
            set { Synonyms = ReadWriteClass.Get_array_from_readLine<string>(value, HMDB_metabolite_self_readWriteOptions_class.Delimiter); }
        }

        public string ReadWrite_protein_associations
        {
            get { return ReadWriteClass.Get_writeLine_from_array(Protein_associations, HMDB_metabolite_self_readWriteOptions_class.Delimiter); }
            set { Protein_associations = ReadWriteClass.Get_array_from_readLine<string>(value, HMDB_metabolite_self_readWriteOptions_class.Delimiter); }
        }

        public string ReadWrite_secondary_hmdb_accesions
        {
            get { return ReadWriteClass.Get_writeLine_from_array(Secondary_hmdb_accessions, HMDB_metabolite_self_readWriteOptions_class.Delimiter); }
            set { Secondary_hmdb_accessions = ReadWriteClass.Get_array_from_readLine<string>(value, HMDB_metabolite_self_readWriteOptions_class.Delimiter); }
        }

        public HMDB_metabolite_line_class()
        {
            this.Protein_associations = new string[0];
        }

        public HMDB_metabolite_line_class Deep_copy()
        {
            HMDB_metabolite_line_class copy = (HMDB_metabolite_line_class)this.MemberwiseClone();
            copy.Metabolite_name = (string)this.Metabolite_name.Clone();
            copy.Synonyms = Array_class.Deep_copy_string_array(this.Synonyms);
            copy.Primary_hmdb_accession = (string)this.Primary_hmdb_accession.Clone();
            copy.Secondary_hmdb_accessions = Array_class.Deep_copy_string_array(this.Secondary_hmdb_accessions);
            copy.Protein_associations = Array_class.Deep_copy_string_array(this.Protein_associations);
            return copy;
        }
    }

    class HMDB_metabolite_class
    {
        public HMDB_metabolite_line_class[] Metabolites { get; set; }

        public HMDB_metabolite_class()
        {
        }

        public string[] Get_all_ordered_distinct_accociated_genes_in_uppserCase()
        {
            HMDB_metabolite_line_class metabolite_line;
            int metabolites_length = this.Metabolites.Length;
            List<string> genes = new List<string>();
            for (int indexMet=0; indexMet<metabolites_length; indexMet++)
            {
                metabolite_line = Metabolites[indexMet];
                foreach (string protein_association in metabolite_line.Protein_associations)
                {
                    genes.Add(protein_association.ToUpper());
                }
            }
            return genes.Distinct().OrderBy(l => l).ToArray();
        }

        public void Generate_by_reading_safed_file()
        {
            Read();
        }

        #region Write
        private void Write()
        {
            HMDB_metabolite_self_readWriteOptions_class readWriteOptions = new HMDB_metabolite_self_readWriteOptions_class();
            ReadWriteClass.WriteData(Metabolites, readWriteOptions);
        }

        private void Read()
        {
            HMDB_metabolite_self_readWriteOptions_class readWriteOptions = new HMDB_metabolite_self_readWriteOptions_class();
            Metabolites = ReadWriteClass.ReadRawData_and_FillArray<HMDB_metabolite_line_class>(readWriteOptions);
        }

        #endregion

    }
}
