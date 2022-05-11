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
using MBCO_global;

namespace MBCO_gene_scp_associations
{
    class MBCO_association_line_class
    {
     
        #region Fields
        public int ProcessLevel { get; set; }
        public string ProcessID { get; set; }
        public string ProcessName { get; set; }
        public string Parent_processName { get; set; }
        public string Symbol { get; set; }
        public string Description { get; set; }
        public string[] References { get; set; }
        public Manual_validation_enum Manual_validation { get; set; }

        public string ReadWrite_references
        {
            get { return ReadWriteClass.Get_writeLine_from_array(this.References, MBCO_association_readOptions_class.Array_delimiter); }
            set { this.References = ReadWriteClass.Get_array_from_readLine<string>(value, MBCO_association_readOptions_class.Array_delimiter); }
        }
        #endregion

        public MBCO_association_line_class()
        {
            ProcessLevel = -1;
            ProcessID = Global_class.Empty_entry;
            ProcessName = Global_class.Empty_entry;
            Parent_processName = Global_class.Empty_entry;
            Description = "";
            References = new string[0];
        }

        public MBCO_association_line_class Deep_copy()
        {
            MBCO_association_line_class copy = (MBCO_association_line_class)this.MemberwiseClone();
            copy.ProcessID = (string)this.ProcessID.Clone();
            copy.ProcessName = (string)this.ProcessName.Clone();
            copy.Parent_processName = (string)this.Parent_processName.Clone();
            copy.Symbol = (string)this.Symbol.Clone();
            copy.Description = (string)this.Description.Clone();
            copy.References = Array_class.Deep_copy_string_array(this.References);
            return copy;
        }
    }

    class MBCO_association_readOptions_class : ReadWriteOptions_base
    {
        public static char Array_delimiter { get { return ';'; } }

        public MBCO_association_readOptions_class(Ontology_type_enum ontology)
        {
            switch (ontology)
            {
                case Ontology_type_enum.Mbco_level3:
                    this.File = Global_directory_class.Complete_human_mbco_association_v11_fileName;
                    break;
                default:
                    throw new Exception();
            }
            Key_propertyNames = new string[] { "ProcessLevel", "Parent_processName", "ProcessID", "ProcessName", "Symbol" };
            Key_columnNames = Key_propertyNames;
            this.HeadlineDelimiters = new char[] { Global_class.Tab };
            this.LineDelimiters = new char[] { Global_class.Tab };
            this.File_has_headline = true;
            this.Report = ReadWrite_report_enum.Report_main;
        }
    }

    class MBCO_association_writeOptions_class : ReadWriteOptions_base
    {
        public MBCO_association_writeOptions_class(Ontology_type_enum ontology)
        {
            switch (ontology)
            {
                case Ontology_type_enum.Mbco_level3:
                    this.File = Global_directory_class.Mbco_dataset_directory;
                    break;
                default:
                    throw new Exception();
            }
            Key_propertyNames = new string[] { "ProcessLevel", "ProcessName", "ProcessID", "Symbol", "Manual_validation", "ReadWrite_references", "Parent_processName" };
            Key_columnNames = Key_propertyNames;
            this.HeadlineDelimiters = new char[] { Global_class.Tab };
            this.LineDelimiters = new char[] { Global_class.Tab };
            this.File_has_headline = true;
            this.Report = ReadWrite_report_enum.Report_main;
        }
    }

    class MBCO_association_class
    {
        public MBCO_association_line_class[] MBCO_associations { get; set; }
        public Ontology_type_enum Ontology { get; private set; }

        public MBCO_association_class()
        {
            MBCO_associations = new MBCO_association_line_class[0];
        }

        #region Generate
        public void Generate_by_reading_safed_file(Ontology_type_enum ontology)
        {
            this.Ontology = ontology;
            Read_mbco_associations();
        }
        #endregion

        #region Order
        public void Order_by_symbol_processName()
        {
            MBCO_associations = MBCO_associations.OrderBy(l => l.Symbol).ThenBy(l => l.ProcessName).ToArray();
        }
        public void Order_by_processName_symbol()
        {
            MBCO_associations = MBCO_associations.OrderBy(l => l.ProcessName).ThenBy(l => l.Symbol).ToArray();
        }
        #endregion

        #region Get
        public string[] Get_all_distinct_ordered_symbols()
        {
            this.Order_by_symbol_processName();
            int onto_length = MBCO_associations.Length;
            MBCO_association_line_class onto_line;
            List<string> all_distinct_ordered_symbols = new List<string>();
            for (int indexOnto = 0; indexOnto < onto_length; indexOnto++)
            {
                onto_line = MBCO_associations[indexOnto];
                if ((indexOnto == 0)
                    || (!onto_line.Symbol.Equals(MBCO_associations[indexOnto - 1].Symbol)))
                {
                    all_distinct_ordered_symbols.Add(onto_line.Symbol);
                }
            }
            return all_distinct_ordered_symbols.OrderBy(l => l).ToArray();
        }

        public string[] Get_all_symbols_of_process_names(params string[] process_names)
        {
            process_names = process_names.Distinct().OrderBy(l => l).ToArray();
            int process_names_length = process_names.Length;
            string process_name;

            int onto_length = MBCO_associations.Length;
            this.MBCO_associations = this.MBCO_associations.OrderBy(l => l.ProcessName).ToArray();
            MBCO_association_line_class onto_association_line;
            List<string> process_symbols_list = new List<string>();
            int indexOnto = 0;
            int stringCompare = -2;

            bool process_name_exists = false;
            for (int indexProcessName = 0; indexProcessName < process_names_length; indexProcessName++)
            {
                process_name = process_names[indexProcessName];
                process_name_exists = false;
                stringCompare = -2;
                while ((indexOnto < onto_length) && (stringCompare <= 0))
                {
                    onto_association_line = MBCO_associations[indexOnto];
                    stringCompare = onto_association_line.ProcessName.CompareTo(process_name);
                    if (stringCompare < 0)
                    {
                        indexOnto++;
                    }
                    else if (stringCompare == 0)
                    {
                        process_symbols_list.Add(onto_association_line.Symbol);
                        indexOnto++;
                        process_name_exists = true;
                    }
                }
                if (!process_name_exists) { throw new Exception("process name does not exist"); }
            }
            return process_symbols_list.ToArray();
        }
        #endregion

        #region Keep, Remove
        public void Keep_only_bg_symbols(string[] bg_symbols)
        {
            bg_symbols = bg_symbols.Distinct().OrderBy(l => l).ToArray();
            string bg_symbol;
            int bg_symbols_length = bg_symbols.Length;
            int indexSymbol = 0;

            this.Order_by_symbol_processName();
            int mbco_associations_length = this.MBCO_associations.Length;
            MBCO_association_line_class mbco_association_line;
            int stringCompare = -2;
            List<MBCO_association_line_class> keep = new List<MBCO_association_line_class>();
            for (int indexMBCO = 0; indexMBCO < mbco_associations_length; indexMBCO++)
            {
                mbco_association_line = this.MBCO_associations[indexMBCO];
                stringCompare = -2;
                while ((indexSymbol < bg_symbols_length) && (stringCompare < 0))
                {
                    bg_symbol = bg_symbols[indexSymbol];
                    stringCompare = bg_symbol.CompareTo(mbco_association_line.Symbol);
                    if (stringCompare < 0)
                    {
                        indexSymbol++;
                    }
                    else if (stringCompare == 0)
                    {
                        keep.Add(mbco_association_line);
                    }
                }
            }
            this.MBCO_associations = keep.ToArray();
        }

        public void Remove_background_genes_scp()
        {
            List<MBCO_association_line_class> keep = new List<MBCO_association_line_class>();
            foreach (MBCO_association_line_class mbco_association_line in this.MBCO_associations)
            {
                if (!mbco_association_line.ProcessName.Equals(MBCO_global_class.Background_genes_scpName))
                {
                    keep.Add(mbco_association_line);
                }
            }
            this.MBCO_associations = keep.ToArray();
        }

        #endregion

        #region Read write copy
        private void Read_mbco_associations()
        {
            MBCO_association_readOptions_class readOptions = new MBCO_association_readOptions_class(this.Ontology);
            this.MBCO_associations = ReadWriteClass.ReadRawData_and_FillArray<MBCO_association_line_class>(readOptions);
        }
        public MBCO_association_class Deep_copy()
        {
            MBCO_association_class copy = (MBCO_association_class)this.MemberwiseClone();
            int associations_length = MBCO_associations.Length;
            copy.MBCO_associations = new MBCO_association_line_class[associations_length];
            for (int indexA = 0; indexA < associations_length; indexA++)
            {
                copy.MBCO_associations[indexA] = this.MBCO_associations[indexA].Deep_copy();
            }
            return copy;
        }
        #endregion

    }
}

