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
using Enumerations;
using ReadWrite;

namespace Enrichment_2018
{
 
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////

    class Ontology_libary_global_options_class
    {
        public static string EnrichR_download_date { get { return "2018July17"; } }
        public static string X2K_download_date { get { return "x2k_v161207_2018July19"; } } //2013October25
        public static string X2Kweb_download_date { get { return "2018August08"; } } //2013October25
    }

    class Ontology_library_line_class
    {
        public int Level { get; set; }
        public string Scp { get; set; }
        public string Target_gene_symbol { get; set; }
        public string Organism_string { get; set; }

        public string Additional_information { get; set; }
        public float Target_gene_score { get; set; }
        public Organism_enum Organism { get; set; }

        public Ontology_type_enum Ontology { get; set; }

        public static bool Check_if_ordered_correctly { get { return Global_class.Check_ordering; } }

        public Ontology_library_line_class()
        {
            this.Level = -1;
            this.Scp = "";
            this.Target_gene_symbol = "";
            this.Organism_string = "";
            this.Additional_information = "";
            this.Target_gene_score = -1;
        }

        #region Order
        public static Ontology_library_line_class[] Order_by_scp_targetGeneSymbol(Ontology_library_line_class[] lines)
        {
            Dictionary<string, Dictionary<string, List<Ontology_library_line_class>>> scp_targetGeneSymbol_dict = new Dictionary<string, Dictionary<string, List<Ontology_library_line_class>>>();
            Dictionary<string, List<Ontology_library_line_class>> targetGeneSymbol_dict = new Dictionary<string, List<Ontology_library_line_class>>();
            int lines_length = lines.Length;
            Ontology_library_line_class library_line;
            for (int indexL = 0; indexL < lines_length; indexL++)
            {
                library_line = lines[indexL];
                if (!scp_targetGeneSymbol_dict.ContainsKey(library_line.Scp))
                {
                    scp_targetGeneSymbol_dict.Add(library_line.Scp, new Dictionary<string, List<Ontology_library_line_class>>());
                }
                if (!scp_targetGeneSymbol_dict[library_line.Scp].ContainsKey(library_line.Target_gene_symbol))
                {
                    scp_targetGeneSymbol_dict[library_line.Scp].Add(library_line.Target_gene_symbol, new List<Ontology_library_line_class>());
                }
                scp_targetGeneSymbol_dict[library_line.Scp][library_line.Target_gene_symbol].Add(library_line);
            }

            string[] scps = scp_targetGeneSymbol_dict.Keys.ToArray();
            string scp;
            int scps_length = scps.Length;
            string[] geneSymbols;
            string geneSymbol;
            int geneSymbols_length;
            scps = scps.OrderBy(l => l).ToArray();
            List<Ontology_library_line_class> ordered_lines = new List<Ontology_library_line_class>();
            for (int indexScp = 0; indexScp < scps_length; indexScp++)
            {
                scp = scps[indexScp];
                targetGeneSymbol_dict = scp_targetGeneSymbol_dict[scp];
                geneSymbols = targetGeneSymbol_dict.Keys.ToArray();
                geneSymbols_length = geneSymbols.Length;
                geneSymbols = geneSymbols.OrderBy(l => l).ToArray();
                for (int indexGS = 0; indexGS < geneSymbols_length; indexGS++)
                {
                    geneSymbol = geneSymbols[indexGS];
                    ordered_lines.AddRange(targetGeneSymbol_dict[geneSymbol]);
                }
            }

            if (Check_if_ordered_correctly)
            {
                #region Check if ordered correctly
                int ordered_length = ordered_lines.Count;
                Ontology_library_line_class previous_line;
                Ontology_library_line_class current_line;
                for (int indexO = 1; indexO < ordered_length; indexO++)
                {
                    previous_line = ordered_lines[indexO - 1];
                    current_line = ordered_lines[indexO];
                    if ((current_line.Scp.CompareTo(previous_line.Scp) < 0)) { throw new Exception(); }
                    if ((current_line.Scp.Equals(previous_line.Scp))
                        && (current_line.Target_gene_symbol.CompareTo(previous_line.Target_gene_symbol) < 0)) { throw new Exception(); }
                }
                #endregion
            }
            return ordered_lines.ToArray();
        }

        public static Ontology_library_line_class[] Order_by_targetGeneSymbol_scp(Ontology_library_line_class[] lines)
        {
            Dictionary<string, Dictionary<string, List<Ontology_library_line_class>>> targetGeneSymbol_scp_dict = new Dictionary<string, Dictionary<string, List<Ontology_library_line_class>>>();
            Dictionary<string, List<Ontology_library_line_class>> scp_dict = new Dictionary<string, List<Ontology_library_line_class>>();
            int lines_length = lines.Length;
            Ontology_library_line_class library_line;
            for (int indexL = 0; indexL < lines_length; indexL++)
            {
                library_line = lines[indexL];
                if (!targetGeneSymbol_scp_dict.ContainsKey(library_line.Target_gene_symbol))
                {
                    targetGeneSymbol_scp_dict.Add(library_line.Target_gene_symbol, new Dictionary<string, List<Ontology_library_line_class>>());
                }
                if (!targetGeneSymbol_scp_dict[library_line.Target_gene_symbol].ContainsKey(library_line.Scp))
                {
                    targetGeneSymbol_scp_dict[library_line.Target_gene_symbol].Add(library_line.Scp, new List<Ontology_library_line_class>());
                }
                targetGeneSymbol_scp_dict[library_line.Target_gene_symbol][library_line.Scp].Add(library_line);
            }

            string[] geneSymbols = targetGeneSymbol_scp_dict.Keys.ToArray();
            string geneSymbol;
            int geneSymbols_length = geneSymbols.Length; ;
            string[] scps;
            string scp;
            int scps_length;
            geneSymbols = geneSymbols.OrderBy(l => l).ToArray();
            List<Ontology_library_line_class> ordered_lines = new List<Ontology_library_line_class>();
            for (int indexGS = 0; indexGS < geneSymbols_length; indexGS++)
            {
                geneSymbol = geneSymbols[indexGS];
                scp_dict = targetGeneSymbol_scp_dict[geneSymbol];
                scps = scp_dict.Keys.ToArray();
                scps_length = scps.Length;
                scps = scps.OrderBy(l => l).ToArray();
                for (int indexScp = 0; indexScp < scps_length; indexScp++)
                {
                    scp = scps[indexScp];
                    ordered_lines.AddRange(scp_dict[scp]);
                }
            }

            if (Check_if_ordered_correctly)
            {
                #region Check if ordered correctly
                int ordered_length = ordered_lines.Count;
                Ontology_library_line_class previous_line;
                Ontology_library_line_class current_line;
                for (int indexO = 1; indexO < ordered_length; indexO++)
                {
                    previous_line = ordered_lines[indexO - 1];
                    current_line = ordered_lines[indexO];
                    if ((current_line.Target_gene_symbol.CompareTo(previous_line.Target_gene_symbol) < 0)) { throw new Exception(); }
                    if ((current_line.Target_gene_symbol.Equals(previous_line.Target_gene_symbol))
                        && (current_line.Scp.CompareTo(previous_line.Scp) < 0)) { throw new Exception(); }
                }
                #endregion
            }
            return ordered_lines.ToArray();
        }


        #endregion

        public Ontology_library_line_class Deep_copy()
        {
            Ontology_library_line_class copy = (Ontology_library_line_class)this.MemberwiseClone();
            copy.Target_gene_symbol = (string)this.Target_gene_symbol.Clone();
            copy.Organism_string = (string)this.Organism_string.Clone();
            copy.Scp = (string)this.Scp.Clone();
            copy.Additional_information = (string)this.Additional_information.Clone();
            return copy;
        }
    }

    class Ontology_libary_readWriteOptions_class : ReadWriteOptions_base
    {
        public Ontology_libary_readWriteOptions_class(Ontology_type_enum ontology, Organism_enum organism)
        {
            string organism_string = "";
            switch (organism)
            {
                case Organism_enum.Homo_sapiens:
                    organism_string = "human";
                    break;
                default:
                    throw new Exception();
            }
            this.File = Global_directory_class.Gene_library_directory + ontology + "_" + organism_string + ".txt";
            this.Key_propertyNames = new string[] { "Ontology", "Scp", "Target_gene_symbol", "Organism" };
            this.Key_columnNames = this.Key_propertyNames;
            this.LineDelimiters = new char[] { Global_class.Tab };
            this.HeadlineDelimiters = new char[] { Global_class.Tab };
            this.File_has_headline = true;
            this.Report = ReadWrite_report_enum.Report_main;
        }
    }

    class Ontology_library_class
    {
        public Ontology_library_line_class[] Library { get; set; }

        public Ontology_library_class()
        {
            this.Library = new Ontology_library_line_class[0];
        }

        #region Order
        public void Order_by_scp_target_gene_symbol()
        {
            this.Library = Ontology_library_line_class.Order_by_scp_targetGeneSymbol(this.Library);
            //this.Library = this.Library.OrderBy(l => l.Scp).ThenBy(l => l.Target_gene_symbol).ToArray();
        }
        #endregion

        public void Add_to_array(Ontology_library_line_class[] add_library)
        {
            int this_length = this.Library.Length;
            int add_length = add_library.Length;
            int new_length = this_length + add_length;
            int indexNew = -1;
            Ontology_library_line_class[] new_library = new Ontology_library_line_class[new_length];
            for (int indexThis = 0; indexThis < this_length; indexThis++)
            {
                indexNew++;
                new_library[indexNew] = this.Library[indexThis];
            }
            for (int indexAdd = 0; indexAdd < add_length; indexAdd++)
            {
                indexNew++;
                new_library[indexNew] = add_library[indexAdd];
            }
            this.Library = new_library;
        }

        #region Generate by reading
        public void Remove_duplicates_based_on_scp_gene_target_symbols()
        {
            this.Library = this.Library.OrderBy(l => l.Scp).ThenBy(l => l.Target_gene_symbol).ToArray();
            int library_length = this.Library.Length;
            Ontology_library_line_class library_line;
            List<Ontology_library_line_class> library_list = new List<Ontology_library_line_class>();
            for (int indexL = 0; indexL < library_length; indexL++)
            {
                library_line = this.Library[indexL];
                if ((indexL == 0)
                    || (!library_line.Scp.Equals(this.Library[indexL - 1].Scp))
                    || (!library_line.Target_gene_symbol.Equals(this.Library[indexL - 1].Target_gene_symbol)))
                {
                    library_line.Target_gene_score = -1;
                    library_line.Additional_information = "only unique lines are kept";
                    library_list.Add(library_line);
                }
            }
            this.Library = library_list.ToArray();
        }

        public void Generate_by_reading(Ontology_type_enum ontology, Organism_enum organism, params string[] bg_genes_in_upperCase)
        {
            Read_and_override_library(ontology, organism);
            Remove_duplicates_based_on_scp_gene_target_symbols();
            if (bg_genes_in_upperCase.Length > 0) { Keep_only_indicated_genes(bg_genes_in_upperCase); }
        }

        private void Read_and_override_library(Ontology_type_enum ontology, Organism_enum organism)
        {
            Ontology_libary_readWriteOptions_class readWriteOptions = new Ontology_libary_readWriteOptions_class(ontology, organism);
            this.Library = ReadWriteClass.ReadRawData_and_FillArray<Ontology_library_line_class>(readWriteOptions);
        }
        #endregion
 
        public void Keep_only_indicated_genes(string[] keep_genes_in_upperCase)
        {
            keep_genes_in_upperCase = keep_genes_in_upperCase.Distinct().OrderBy(l => l).ToArray();
            int keep_genes_length = keep_genes_in_upperCase.Length;
            int indexKeep = 0;
            string keep_gene;

            int stringCompare = -2;
            Ontology_library_line_class library_line;
            int library_length = this.Library.Length;
            List<Ontology_library_line_class> keep = new List<Ontology_library_line_class>();
            Library = Library.OrderBy(l => l.Target_gene_symbol).ToArray();

            for (int indexL = 0; indexL < library_length; indexL++)
            {
                library_line = Library[indexL];
                stringCompare = -2;
                while ((indexKeep < keep_genes_length) && (stringCompare < 0))
                {
                    keep_gene = keep_genes_in_upperCase[indexKeep];
                    stringCompare = keep_gene.CompareTo(library_line.Target_gene_symbol);
                    if (stringCompare < 0)
                    {
                        indexKeep++;
                    }
                    else if (stringCompare == 0)
                    {
                        keep.Add(library_line);
                    }
                }
            }
            this.Library = keep.ToArray();
        }

        #region Get scps and target gene symbols
        public Ontology_type_enum Get_current_ontology_and_check_if_only_one()
        {
            Ontology_type_enum ontology = this.Library[0].Ontology;
            foreach (Ontology_library_line_class ontology_library_line in this.Library)
            {
                if (!ontology_library_line.Ontology.Equals(ontology))
                {
                    throw new Exception();
                }
            }
            return ontology;
        }
        public string[] Get_all_ordered_unique_gene_symbols()
        {
            int library_length = this.Library.Length;
            Ontology_library_line_class library_line;
            List<string> all_symbols = new List<string>();
            for (int indexL = 0; indexL < library_length; indexL++)
            {
                library_line = this.Library[indexL];
                all_symbols.Add(library_line.Target_gene_symbol);
            }
            return all_symbols.Distinct().OrderBy(l=>l).ToArray();
        }
        #endregion
    }
}
