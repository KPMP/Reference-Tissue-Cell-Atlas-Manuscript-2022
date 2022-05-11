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

namespace MBCO_data
{
    interface IAdd_to_data
    {
        string NCBI_official_symbol_for_data { get; }
        string SampleName_for_data { get; }
        double Value_for_data { get; }
    }

    ///////////////////////////////////////////////////////////////

    class Colchar_column_line_class
    {
        #region Fields
        public string SampleName { get; set; }
        #endregion

        public Colchar_column_line_class()
        {
        }

        public Colchar_column_line_class(string sampleName)
        {
            SampleName = (string)sampleName.Clone();
        }

        #region Standard
        public static Colchar_column_line_class[] Order_in_standard_way(Colchar_column_line_class[] columns)
        {
            return columns.OrderBy(l => l.SampleName).ToArray();
        }

        public bool Equal_in_standard_way(IAdd_to_data other)
        {
            return this.SampleName.Equals(other.SampleName_for_data);
        }

        public bool Equal_in_standard_way(Colchar_column_line_class other)
        {
            return this.SampleName.Equals(other.SampleName);
        }

        public static bool Equal_in_standard_way(Colchar_column_line_class line1, Colchar_column_line_class line2)
        {
            return line1.SampleName.Equals(line2.SampleName);
        }

        public static IAdd_to_data[] Order_in_standard_way(IAdd_to_data[] add_lines)
        {
            return add_lines.OrderBy(l => l.SampleName_for_data).ToArray();
        }
        #endregion

        public Colchar_column_line_class Deep_copy()
        {
            Colchar_column_line_class copy = (Colchar_column_line_class)this.MemberwiseClone();
            copy.SampleName = (string)this.SampleName.Clone();
            return copy;
        }
    }

    class Colchar_class
    {
        #region Fields
        bool column_rearrangements_adopted;

        public Colchar_column_line_class[] Columns { get; set; }
        public int Columns_length { get { return Columns.Length; } }
        public bool Column_rearrangements_adopted
        {
            get { return column_rearrangements_adopted; }
            set
            {
                if (value.Equals(false))
                {
                    SampleName_index_dict.Clear();
                }
                column_rearrangements_adopted = value;
            }
        }
        private Dictionary<string, int> SampleName_index_dict { get; set; }
        #endregion

        public Colchar_class()
        {
            Column_rearrangements_adopted = true;
            Columns = new Colchar_column_line_class[0];
            SampleName_index_dict = new Dictionary<string, int>();
        }

        #region Check
        public void Correctness_check()
        {
            int col_length = Columns_length;
            Colchar_column_line_class[] columns_copy = Deep_copy_columns();
            columns_copy = Colchar_column_line_class.Order_in_standard_way(columns_copy);
            for (int indexC = 1; indexC < col_length; indexC++)
            {
                if (columns_copy[indexC].Equal_in_standard_way(columns_copy[indexC - 1]))
                {
                    throw new Exception("duplicated column charactarization");
                }
            }
            if (!Column_rearrangements_adopted)
            {
                throw new Exception("column rearrangements are not adopted");
            }
        }
        #endregion

        #region Get
        public int Get_index_of_iadd_data_line_from_dictionary(IAdd_to_data add_data_line)
        {
            return SampleName_index_dict[add_data_line.SampleName_for_data];
        }
        #endregion

        public void Generate_sampleName_index_dict()
        {
            SampleName_index_dict.Clear();
            int col_length = this.Columns_length;
            Colchar_column_line_class column_line;
            for (int indexCol = 0; indexCol < col_length; indexCol++)
            {
                column_line = this.Columns[indexCol];
                if (!SampleName_index_dict.ContainsKey(column_line.SampleName))
                {
                    SampleName_index_dict.Add(column_line.SampleName, indexCol);
                }
                else
                {
                    throw new Exception();
                }
            }
        }


        #region Keep
        public void Keep_only_input_columns(params int[] inputColumns)
        {
            Correctness_check();
            int inputColumns_length = inputColumns.Length;
            Colchar_column_line_class[] new_columns = new Colchar_column_line_class[inputColumns_length];
            for (int indexInput = 0; indexInput < inputColumns_length; indexInput++)
            {
                new_columns[indexInput] = this.Columns[inputColumns[indexInput]].Deep_copy();
            }
            Columns = new_columns;
            Column_rearrangements_adopted = false;
        }
        #endregion

        #region Add
        public void Identify_new_columns_and_add_at_right_site(Colchar_column_line_class[] add_data)
        {
            Correctness_check();

            int old_columns_length = Columns_length;
            Colchar_column_line_class[] add_columns;

            #region Identify add columns
            List<Colchar_column_line_class> add_columns_list = new List<Colchar_column_line_class>();
            add_data = Colchar_column_line_class.Order_in_standard_way(add_data);
            int add_data_length = add_data.Length;
            Colchar_column_line_class add_line;
            bool old_columns_contain_new_combination = false;
            Colchar_column_line_class old_column_line;
            Colchar_column_line_class add_column_line;
            for (int indexAdd = 0; indexAdd < add_data_length; indexAdd++)
            {
                add_line = add_data[indexAdd];
                if ((indexAdd == 0) || (!Colchar_column_line_class.Equal_in_standard_way(add_line, add_data[indexAdd - 1])))
                {
                    old_columns_contain_new_combination = false;
                    for (int indexOld = 0; indexOld < old_columns_length; indexOld++)
                    {
                        old_column_line = Columns[indexOld];
                        if (old_column_line.Equal_in_standard_way(add_line))
                        {
                            old_columns_contain_new_combination = true;
                            break;
                        }
                    }
                    if (!old_columns_contain_new_combination)
                    {
                        add_column_line = new Colchar_column_line_class(add_line.SampleName);
                        add_columns_list.Add(add_column_line);
                    }
                }
            }
            add_columns = add_columns_list.ToArray();
            #endregion

            #region Add add columns at right site
            int add_columns_length = add_columns.Length;
            int new_columns_length = add_columns_length + old_columns_length;
            Colchar_column_line_class[] new_columns = new Colchar_column_line_class[new_columns_length];
            int indexNew = -1;
            for (int indexOld = 0; indexOld < old_columns_length; indexOld++)
            {
                indexNew++;
                new_columns[indexNew] = this.Columns[indexOld].Deep_copy();
            }
            for (int indexAdd = 0; indexAdd < add_columns_length; indexAdd++)
            {
                indexNew++;
                new_columns[indexNew] = add_columns[indexAdd].Deep_copy();
            }
            Columns = new_columns;
            #endregion

            Correctness_check();
        }
        #endregion

        #region Copy
        private Colchar_column_line_class[] Deep_copy_columns()
        {
            int columns_length = Columns.Length;
            Colchar_column_line_class[] copy_columns = new Colchar_column_line_class[columns_length];
            for (int indexC = 0; indexC < columns_length; indexC++)
            {
                copy_columns[indexC] = this.Columns[indexC].Deep_copy();
            }
            return copy_columns;
        }

        public Colchar_class Deep_copy()
        {
            Colchar_class copy = (Colchar_class)this.MemberwiseClone();
            int columns_length = Columns.Length;
            copy.Columns = Deep_copy_columns();
            return copy;
        }
        #endregion
    }

    ///////////////////////////////////////////////////////////////

    class Data_line_class
    {
        public static double Empty_entry { get { return 0; } }
        public double[] Columns { get; set; }
        public int Columns_length { get { return Columns.Length; } }
        public string NCBI_official_symbol { get; set; }

        public string NCBI_description { get; set; }

        public Data_line_class(string ncbi_symbol, int columns_length)
        {
            NCBI_official_symbol = (string)ncbi_symbol.Clone();
            NCBI_description = "";
            Columns = new double[columns_length];
            for (int indexCol = 0; indexCol < columns_length; indexCol++)
            {
                Columns[indexCol] = Empty_entry;
            }
        }

        public void Add_to_this_line_after_checking_if_this_line_has_empty_entry(IAdd_to_data add_line, int indexCol)
        {
            if (!NCBI_official_symbol.Equals(add_line.NCBI_official_symbol_for_data))
            {
                throw new Exception("rowNames do not match");
            }
            if (Columns[indexCol].Equals(Empty_entry))
            {
                Columns[indexCol] = add_line.Value_for_data;
            }
            else
            {
                throw new Exception("Position already filled");
            }
        }

        public void Add_nonemtpyt_values_of_other_line_to_this_line_if_this_line_has_empty_value(Data_line_class other)
        {
            int columns_length = Columns_length;
            for (int indexC = 0; indexC < columns_length; indexC++)
            {
                if ((this.Columns[indexC] != Empty_entry)
                    && (other.Columns[indexC] != Empty_entry))
                {
                    throw new Exception();
                }
                else if (other.Columns[indexC] != Empty_entry)
                {
                    this.Columns[indexC] = other.Columns[indexC];
                }
            }
        }

        public void Keep_columns(params int[] kept_columnIndexes)
        {
            List<double> keep_columns = new List<double>();
            int columns_length = this.Columns_length;
            for (int indexColumn = 0; indexColumn < columns_length; indexColumn++)
            {
                if (kept_columnIndexes.Contains(indexColumn))
                {
                    keep_columns.Add(this.Columns[indexColumn]);
                }
            }
            this.Columns = keep_columns.ToArray();
        }

        #region Copy
        private double[] Deep_copy_columns()
        {
            int columns_length = Columns.Length;
            double[] copy = new double[columns_length];
            for (int indexC = 0; indexC < columns_length; indexC++)
            {
                copy[indexC] = Columns[indexC];
            }
            return copy;
        }

        public Data_line_class Deep_copy()
        {
            Data_line_class copy = (Data_line_class)this.MemberwiseClone();
            copy.NCBI_official_symbol = (string)this.NCBI_official_symbol.Clone();
            copy.NCBI_description = (string)this.NCBI_description.Clone();
            copy.Columns = Deep_copy_columns();
            return copy;
        }
        #endregion
    }

    class Data_class
    {
        #region Fields
        public Data_line_class[] Data { get; set; }
        public Data_line_class Column_entries_count_line { get; set; }
        public int Data_length { get { return Data.Length; } }
        public Colchar_class ColChar { get; set; }
        #endregion

        public Data_class()
        {
            ColChar = new Colchar_class();
            Data = new Data_line_class[0];
            Column_entries_count_line = new Data_line_class("Column entries count", 0);
        }

        #region Check
        public void Correctness_check()
        {
            ColChar.Correctness_check();
            int col_length = ColChar.Columns.Length;
            int data_length = Data.Length;
            Data_line_class data_line;
            for (int indexD = 0; indexD < data_length; indexD++)
            {
                data_line = Data[indexD];
                if (data_line.Columns.Length != col_length)
                {
                    throw new Exception("Column lengths do not match");
                }
            }
        }
        #endregion

        #region Order
        public void Order_by_ncbiOfficialSymbol()
        {
            Data = Data.OrderBy(l => l.NCBI_official_symbol).ToArray();
        }
        #endregion

        #region Fill de instance
        public void Add_to_data_instance(IAdd_to_data[] add_data)
        {
            int old_columns_length = ColChar.Columns_length;
            int add_data_length = add_data.Length;
            IAdd_to_data add_data_interface_line;

            #region Analyze, if new columns needed and add to colChar
            add_data = add_data.OrderBy(l => l.SampleName_for_data).ToArray();
            List<Colchar_column_line_class> colChar_column_list = new List<Colchar_column_line_class>();
            Colchar_column_line_class colChar_column_line;
            for (int indexAdd = 0; indexAdd < add_data_length; indexAdd++)
            {
                add_data_interface_line = add_data[indexAdd];
                if ((indexAdd == 0)
                    || (!add_data_interface_line.SampleName_for_data.Equals(add_data[indexAdd - 1].SampleName_for_data)))
                {
                    colChar_column_line = new Colchar_column_line_class();
                    colChar_column_line.SampleName = (string)add_data_interface_line.SampleName_for_data.Clone();
                    colChar_column_list.Add(colChar_column_line);
                }
            }
            ColChar.Identify_new_columns_and_add_at_right_site(colChar_column_list.ToArray());
            ColChar.Generate_sampleName_index_dict();
            int new_columns_length = ColChar.Columns.Length;
            #endregion

            int data_length = this.Data_length;
            Data_line_class data_line;

            #region Extend data lines by adding zero values at right site
            double[] new_columns;
            for (int indexD = 0; indexD < data_length; indexD++)
            {
                data_line = this.Data[indexD];
                new_columns = new double[new_columns_length];
                for (int indexOld = 0; indexOld < old_columns_length; indexOld++)
                {
                    new_columns[indexOld] = data_line.Columns[indexOld];
                }
                data_line.Columns = new_columns.ToArray();
            }
            #endregion

            #region Add to data and generate ordered data add lines to data
            add_data = add_data.OrderBy(l => l.NCBI_official_symbol_for_data).ToArray();

            int indexData = 0;
            int stringCompare;
            Data_line_class new_data_line;
            List<Data_line_class> add_lines_to_data_list = new List<Data_line_class>();
            int indexCol = -1;

            for (int indexAdd = 0; indexAdd < add_data_length; indexAdd++)
            {
                add_data_interface_line = add_data[indexAdd];
                indexCol = ColChar.Get_index_of_iadd_data_line_from_dictionary(add_data_interface_line);
                stringCompare = -2;
                while ((indexData < data_length) && (stringCompare < 0))
                {
                    data_line = Data[indexData];
                    stringCompare = data_line.NCBI_official_symbol.CompareTo(add_data_interface_line.NCBI_official_symbol_for_data);
                    if (stringCompare < 0)
                    {
                        indexData++;
                    }
                    else if (stringCompare == 0)
                    {
                        data_line.Add_to_this_line_after_checking_if_this_line_has_empty_entry(add_data_interface_line, indexCol);
                    }
                }
                if (stringCompare != 0)
                {
                    new_data_line = new Data_line_class(add_data_interface_line.NCBI_official_symbol_for_data, new_columns_length);
                    new_data_line.Columns[indexCol] = add_data_interface_line.Value_for_data;
                    add_lines_to_data_list.Add(new_data_line);
                }
            }
            Data_line_class[] add_lines_to_data = add_lines_to_data_list.OrderBy(l => l.NCBI_official_symbol).ToArray();
            #endregion

            #region Combine add_lines_to_data to generate final add_data_lines
            int add_lines_length = add_lines_to_data.Length;
            int firstIndex_same_rowName = -1;
            Data_line_class add_data_line;
            Data_line_class inner_add_data_line;
            List<Data_line_class> final_add_data_list = new List<Data_line_class>();
            for (int indexAL = 0; indexAL < add_lines_length; indexAL++)
            {
                add_data_line = add_lines_to_data[indexAL];
                if ((indexAL == 0) || (!add_data_line.NCBI_official_symbol.Equals(add_lines_to_data[indexAL - 1].NCBI_official_symbol)))
                {
                    firstIndex_same_rowName = indexAL;
                }
                if ((indexAL == add_lines_length - 1) || (!add_data_line.NCBI_official_symbol.Equals(add_lines_to_data[indexAL + 1].NCBI_official_symbol)))
                {
                    for (int indexInner = firstIndex_same_rowName; indexInner < indexAL; indexInner++)
                    {
                        inner_add_data_line = add_lines_to_data[indexInner];
                        add_data_line.Add_nonemtpyt_values_of_other_line_to_this_line_if_this_line_has_empty_value(inner_add_data_line);
                    }
                    final_add_data_list.Add(add_data_line);
                }
            }
            #endregion

            #region Add final data line to data
            int final_data_length = final_add_data_list.Count;
            int new_data_length = final_data_length + data_length;
            Data_line_class[] new_data = new Data_line_class[new_data_length];
            int indexNew = -1;
            for (int indexOld = 0; indexOld < data_length; indexOld++)
            {
                indexNew++;
                new_data[indexNew] = Data[indexOld];
            }
            for (int indexFinal = 0; indexFinal < final_data_length; indexFinal++)
            {
                indexNew++;
                new_data[indexNew] = final_add_data_list[indexFinal];
            }
            Data = new_data;
            #endregion

            Keep_only_lines_that_contain_at_least_one_non_zero_value();
        }
        #endregion

        #region Keep
        public void Keep_only_input_rowNames(string[] input_rowNames)
        {
            List<Data_line_class> keep_data = new List<Data_line_class>();
            List<Data_line_class> remove_keep_data = new List<Data_line_class>();
            input_rowNames = input_rowNames.Distinct().OrderBy(l => l).ToArray();
            string input_rowName;
            int indexInput = 0;
            int input_rowNames_length = input_rowNames.Length;
            this.Order_by_ncbiOfficialSymbol();
            int data_length = Data_length;
            int stringCompare = -2;
            Data_line_class data_line;
            bool kept = false;
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                data_line = Data[indexData];
                kept = false;
                stringCompare = -2;
                while ((indexInput < input_rowNames_length) && (stringCompare < 0))
                {
                    input_rowName = input_rowNames[indexInput];
                    stringCompare = input_rowName.CompareTo(data_line.NCBI_official_symbol);
                    if (stringCompare < 0)
                    {
                        indexInput++;
                    }
                    else if (stringCompare == 0)
                    {
                        keep_data.Add(data_line);
                        kept = true;
                    }
                }
                if (!kept)
                {
                    remove_keep_data.Add(data_line);
                }
            }
            Data = keep_data.ToArray();
        }

        public void Remove_empty_rows_and_columns()
        {
            int column_length = this.ColChar.Columns_length;
            List<int> keepColumns = new List<int>();
            bool keep_row;
            List<Data_line_class> newDE = new List<Data_line_class>();
            foreach (Data_line_class line in Data)
            {
                keep_row = false;
                for (int indexCol = 0; indexCol < column_length; indexCol++)
                {
                    if (line.Columns[indexCol] != 0)
                    {
                        keepColumns.Add(indexCol);
                        keep_row = true;
                    }
                }
                if (keep_row)
                {
                    newDE.Add(line);
                }
            }

            if (keepColumns.Count() < column_length)
            {
                ColChar.Keep_only_input_columns(keepColumns.ToArray());
                foreach (Data_line_class line in Data)
                {
                    line.Keep_columns(keepColumns.ToArray());
                }
                ColChar.Column_rearrangements_adopted = true;
            }
            Data = newDE.ToArray();
        }

        public void Keep_only_input_columns_and_remove_all_rows_that_are_left_over_with_only_zero_values(params int[] inputColumns)
        {
            inputColumns = inputColumns.Distinct().OrderBy(l => l).ToArray();
            ColChar.Keep_only_input_columns(inputColumns);
            inputColumns = inputColumns.Distinct().OrderBy(l => l).ToArray();
            int inputColumns_length = inputColumns.Length;
            int data_length = Data.Length;
            Data_line_class data_line;
            double[] new_columns;
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                data_line = Data[indexData];
                new_columns = new double[inputColumns_length];
                for (int indexI = 0; indexI < inputColumns_length; indexI++)
                {
                    new_columns[indexI] = data_line.Columns[inputColumns[indexI]];
                }
                data_line.Columns = new_columns;
            }
            Remove_empty_rows_and_columns();
            ColChar.Column_rearrangements_adopted = true;
        }
        #endregion

        #region Set to upper case
        public void Set_all_ncbi_official_gene_symbols_to_upper_case()
        {
            foreach (Data_line_class data_line in Data)
            {
                data_line.NCBI_official_symbol = data_line.NCBI_official_symbol.ToUpper();
            }
        }
        #endregion

        #region Keep only above or below cutoff
        private void Keep_only_lines_that_contain_at_least_one_non_zero_value()
        {
            int data_length = this.Data_length;
            int column_length = this.ColChar.Columns_length;
            Data_line_class data_line;
            List<Data_line_class> kept_data_list = new List<Data_line_class>();
            bool keep_line;
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                data_line = this.Data[indexData];
                keep_line = false;
                for (int indexCol = 0; indexCol < column_length; indexCol++)
                {
                    if (data_line.Columns[indexCol] != 0) { keep_line = true; }
                }
                if (keep_line) { kept_data_list.Add(data_line); }
            }
            Data = kept_data_list.ToArray();
        }
        #endregion

        #region Write read copy
        public Data_class Deep_copy()
        {
            Data_class copy = (Data_class)this.MemberwiseClone();
            int data_length = Data.Length;
            copy.Data = new Data_line_class[data_length];
            copy.Column_entries_count_line = this.Column_entries_count_line.Deep_copy();
            for (int indexD = 0; indexD < data_length; indexD++)
            {
                copy.Data[indexD] = this.Data[indexD].Deep_copy();
            }
            copy.ColChar = this.ColChar.Deep_copy();
            return copy;
        }
        #endregion
    }
}
