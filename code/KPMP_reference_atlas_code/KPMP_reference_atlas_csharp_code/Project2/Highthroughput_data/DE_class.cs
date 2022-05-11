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
using ReadWrite;
using Enumerations;

namespace Highthroughput_data
{
    interface IFill_de
    {
        string[] Symbols_for_de { get; }
        string[] Names_for_de { get; }
        double Value_for_de { get; }
    }

    class Fill_de_line_class : IFill_de
    {
        public string[] Symbols_for_de { get; set; }
        public string[] Names_for_de { get; set; }
        public double Value_for_de { get; set; }

        public static bool Check_if_ordered { get { return true; } }

        public int Column { get; set; }

        public string Combined_names { get; set; }

        public string Combined_symbols { get; set; }

        public string ReadWrite_symbols_for_de
        {
            get { return ReadWriteClass.Get_writeLine_from_array(Symbols_for_de, Fill_de_readWriteOptions.Delimiter); }
            set { this.Symbols_for_de = ReadWriteClass.Get_array_from_readLine<string>(value, Fill_de_readWriteOptions.Delimiter); }
        }

        public string ReadWrite_names_for_de
        {
            get { return ReadWriteClass.Get_writeLine_from_array(Names_for_de, Fill_de_readWriteOptions.Delimiter); }
            set { this.Names_for_de = ReadWriteClass.Get_array_from_readLine<string>(value, Fill_de_readWriteOptions.Delimiter); }
        }

        public Fill_de_line_class()
        {
            Symbols_for_de = new string[] { "" };
            Names_for_de = new string[] { "" };
            Combined_names = "";
            Combined_symbols = "";
        }

        #region Order
        public static Fill_de_line_class[] Order_by_combinedNames(Fill_de_line_class[] fill_de_lines)
        {
            Dictionary<string, List<Fill_de_line_class>> combinedNames_dict = new Dictionary<string, List<Fill_de_line_class>>();
            int fill_de_lines_length = fill_de_lines.Length;
            Fill_de_line_class fill_de_line;
            for (int indexFill = 0; indexFill < fill_de_lines_length; indexFill++)
            {
                fill_de_line = fill_de_lines[indexFill];
                if (!combinedNames_dict.ContainsKey(fill_de_line.Combined_names))
                {
                    combinedNames_dict.Add(fill_de_line.Combined_names, new List<Fill_de_line_class>());
                }
                combinedNames_dict[fill_de_line.Combined_names].Add(fill_de_line);
            }
            string[] combinedNames = combinedNames_dict.Keys.ToArray();
            string combinedName;
            int combinedNames_length = combinedNames.Length;
            List<Fill_de_line_class> ordered_fill_de_list = new List<Fill_de_line_class>();
            combinedNames = combinedNames.OrderBy(l => l).ToArray();
            for (int indexC = 0; indexC < combinedNames_length; indexC++)
            {
                combinedName = combinedNames[indexC];
                ordered_fill_de_list.AddRange(combinedNames_dict[combinedName]);
            }

            if (Check_if_ordered)
            {
                #region Check if ordered
                int ordered_length = ordered_fill_de_list.Count;
                Fill_de_line_class previous_line;
                Fill_de_line_class current_line;
                for (int indexO = 1; indexO < ordered_length; indexO++)
                {
                    previous_line = ordered_fill_de_list[indexO - 1];
                    current_line = ordered_fill_de_list[indexO];
                    if (current_line.Combined_names.CompareTo(previous_line.Combined_names) < 0) { throw new Exception(); }
                }
            }
            #endregion
            return ordered_fill_de_list.ToArray();
        }
        #endregion


        private bool Equal_symbols(string[] other_symbols_for_de)
        {
            bool equals = true;
            int this_symbols_length = this.Symbols_for_de.Length;
            int other_symbols_length = other_symbols_for_de.Length;
            if (this_symbols_length != other_symbols_length)
            {
                equals = false;
            }
            else
            {
                for (int indexS = 0; indexS < this_symbols_length; indexS++)
                {
                    if (this.Symbols_for_de[indexS].Equals(other_symbols_for_de[indexS]))
                    {
                        equals = false;
                        break;
                    }
                }
            }
            return equals;
        }

        public void Set_combined_symbols_and_names()
        {
            int symbols_length = this.Symbols_for_de.Length;
            int names_length = this.Names_for_de.Length;
            StringBuilder sb = new StringBuilder();
            sb.Clear();
            for (int indexS = 0; indexS < symbols_length; indexS++)
            {
                if (indexS != 0) { sb.AppendFormat("@"); }
                sb.AppendFormat("{0}", Symbols_for_de[indexS]);
            }
            this.Combined_symbols = sb.ToString();
            sb.Clear();
            for (int indexN = 0; indexN < names_length; indexN++)
            {
                if (indexN != 0) { sb.AppendFormat("@"); }
                sb.AppendFormat("{0}", Names_for_de[indexN]);
            }
            this.Combined_names = sb.ToString();
        }

        public bool Equals(Fill_de_line_class other)
        {
            return Equal_symbols(other.Symbols_for_de);
        }

        public Fill_de_line_class Deep_copy()
        {
            Fill_de_line_class copy = (Fill_de_line_class)this.MemberwiseClone();
            copy.Names_for_de = Array_class.Deep_copy_string_array(this.Names_for_de);
            copy.Symbols_for_de = Array_class.Deep_copy_string_array(this.Symbols_for_de);
            copy.Combined_names = (string)this.Combined_names.Clone();
            copy.Combined_symbols = (string)this.Combined_symbols.Clone();
            return copy;
        }
    }

    class Fill_de_readWriteOptions : ReadWriteOptions_base
    {
        public static char Delimiter { get { return ';'; } }

        public Fill_de_readWriteOptions(string subdirectory, string fileName)
        {
            string directory = Global_directory_class.Results_directory + subdirectory;
            ReadWriteClass.Create_directory_if_it_does_not_exist(directory);
            this.File = directory + fileName;
            this.Key_propertyNames = new string[] { "ReadWrite_names_for_de", "Entry_type_for_de", "Timepoint_for_de", "ReadWrite_symbols_for_de", "Value_for_de" };
            this.Key_columnNames = this.Key_propertyNames;
            this.HeadlineDelimiters = new char[] { Global_class.Tab };
            this.LineDelimiters = new char[] { Global_class.Tab };
            this.Report = ReadWrite_report_enum.Report_main;
            this.File_has_headline = true;
        }
    }

    class DE_columns_line_class
    {
        public double Value { get; set; }

        public DE_columns_line_class Deep_copy()
        {
            DE_columns_line_class copy = (DE_columns_line_class)this.MemberwiseClone();
            return copy;
        }
    }

    class DE_column_characterization_line_class
    {
        #region Fields
        public int Index { get; set; }
        public int IndexOld { get; set; }
        public string[] Names { get; set; }

        public string Combined_names
        {
            get
            {
                StringBuilder sb = new StringBuilder();
                foreach (string name in Names)
                {
                    sb.AppendFormat(name);
                }
                return sb.ToString();
            }
        }

        public string ReadWrite_names
        {
            get { return ReadWriteClass.Get_writeLine_from_array(this.Names, ';'); }
            set { this.Names = ReadWriteClass.Get_array_from_readLine<string>(value, ';'); }
        }
        #endregion

        public DE_column_characterization_line_class()
        {
            Names = new string[] { "" };
        }

        public string[] Get_names()
        {
            int names_length = Names.Length;
            string[] copy = new string[names_length];
            for (int indexN = 0; indexN < names_length; indexN++)
            {
                copy[indexN] = (string)Names[indexN].Clone();
            }
            return copy;
        }

        public string Get_full_column_name()
        {
            int names_length = Names.Length;
            StringBuilder name = new StringBuilder();
            for (int indexN = 0; indexN < names_length; indexN++)
            {
                if (indexN == 0)
                {
                    name.AppendFormat("{0}", Names[indexN]);
                }
                else
                {
                    name.AppendFormat("-{0}", Names[indexN]);
                }
            }
            return name.ToString();
        }

        public string Get_full_column_label()
        {
            StringBuilder name = new StringBuilder();
            bool first_add = true;
            int names_length = Names.Length;
            if (first_add)
            {
                name.AppendFormat("{0}", Get_full_column_name());
            }
            else
            {
                name.AppendFormat("-{0}", Get_full_column_name());
            }
            return name.ToString();
        }

        public void Add_new_names(params string[] add_names)
        {
            int this_length = Names.Length;
            int add_length = add_names.Length;
            int new_length = this_length + add_length;
            string[] newNames = new string[new_length];
            int indexNew = -1;
            for (int indexThis = 0; indexThis < this_length; indexThis++)
            {
                indexNew++;
                newNames[indexNew] = (string)Names[indexThis].Clone();
            }
            for (int indexAdd = 0; indexAdd < add_length; indexAdd++)
            {
                indexNew++;
                newNames[indexNew] = (string)add_names[indexAdd].Clone();
            }
            this.Names = newNames;
        }

        public static DE_column_characterization_line_class[] Order_in_standard_way_for_equal_comparison(DE_column_characterization_line_class[] column_char_lines)
        {
            return column_char_lines.OrderBy(l => l.Combined_names).ToArray();
        }

        public bool Equals_other(DE_column_characterization_line_class other)
        {
            bool equal = true;
            int this_names_length = this.Names.Length;
            int other_names_length = other.Names.Length;
            if (this_names_length != other_names_length)
            {
                equal = false;
            }
            else
            {
                for (int indexN = 0; indexN < this_names_length; indexN++)
                {
                    if (!this.Names[indexN].Equals(other.Names[indexN]))
                    {
                        equal = false;
                        break;
                    }
                }
            }
            return equal;
        }

        public DE_column_characterization_line_class Deep_copy()
        {
            DE_column_characterization_line_class copy = (DE_column_characterization_line_class)this.MemberwiseClone();
            copy.Names = Deep_copy_class.Deep_copy_string_array(this.Names);
            return copy;
        }
    }

    class DE_column_characterization_class
    {
        #region Fields
        public List<DE_column_characterization_line_class> Columns { get; set; }
        public bool Index_changes_adopted { get; set; }
        #endregion

        #region Constructors
        public DE_column_characterization_class()
        {
            Columns = new List<DE_column_characterization_line_class>();
            Index_changes_adopted = true;
        }
        #endregion

        #region Check
        public bool Small_correctness_check()
        {
            bool everything_correct = true;
            if (!Index_changes_adopted)
            {
                throw new Exception();
            }
            return everything_correct;
        }

        public bool Check_if_inputColumns_do_not_already_exist(DE_column_characterization_line_class[] otherColChar)
        {
            bool ok = true;
            int other_length = otherColChar.Length;
            int this_count = this.Columns.Count;
            DE_column_characterization_line_class this_colChar_line;
            DE_column_characterization_line_class other_colChar_line;
            for (int indexThis = 0; indexThis < this_count; indexThis++)
            {
                this_colChar_line = Columns[indexThis];
                for (int indexOther = 0; indexOther < other_length; indexOther++)
                {
                    other_colChar_line = otherColChar[indexOther];
                    if (other_colChar_line.Equals_other(this_colChar_line))
                    {
                        throw new Exception();
                        ok = false;
                    }
                }
            }
            return ok;
        }
        #endregion

        #region bool
        public bool Does_column_already_exist(DE_column_characterization_line_class other_colChar_line)
        {
            bool other_column_already_exists = false;
            int this_count = this.Columns.Count;
            DE_column_characterization_line_class this_colChar_line;
            for (int indexThis = 0; indexThis < this_count; indexThis++)
            {
                this_colChar_line = Columns[indexThis];
                if (other_colChar_line.Equals_other(this_colChar_line))
                {
                    other_column_already_exists = true;
                    break;
                }
            }
            return other_column_already_exists;
        }

        public bool Equals_other(DE_column_characterization_class other)
        {
            int this_col_count = this.Columns.Count;
            int other_col_count = other.Columns.Count;
            if (this_col_count != other_col_count)
            {
                return false;
            }
            else
            {
                bool equal = true;
                DE_column_characterization_line_class this_line;
                DE_column_characterization_line_class other_line;
                int this_names_length;
                int other_names_length;
                for (int indexCol = 0; indexCol < this_col_count; indexCol++)
                {
                    this_line = this.Columns[indexCol];
                    other_line = other.Columns[indexCol];
                    this_names_length = this_line.Names.Length;
                    other_names_length = other_line.Names.Length;
                    if (this_names_length != other_names_length)
                    {
                        equal = false;
                        break;
                    }
                    for (int indexN = 0; indexN < this_names_length; indexN++)
                    {
                        if (!this_line.Names[indexN].Equals(other_line.Names[indexN]))
                        {
                            equal = false;
                            break;
                        }
                    }
                }
                return equal;
            }
        }
        #endregion

        public void Clear()
        {
            this.Columns = new List<DE_column_characterization_line_class>();
        }

        #region Get colIndexes

        public int[] Get_colIndexes_of_column_characterization_lines(params DE_column_characterization_line_class[] column_characterization_lines)
        {
            column_characterization_lines = DE_column_characterization_line_class.Order_in_standard_way_for_equal_comparison(column_characterization_lines);
            this.Columns = DE_column_characterization_line_class.Order_in_standard_way_for_equal_comparison(this.Columns.ToArray()).ToList();
            int input_length = column_characterization_lines.Length;
            DE_column_characterization_line_class input_line;
            DE_column_characterization_line_class this_line;
            int this_length = this.Columns.Count;
            List<int> colIndexes_list = new List<int>();
            for (int indexInput = 0; indexInput < input_length; indexInput++)
            {
                input_line = column_characterization_lines[indexInput];
                for (int indexThis = 0; indexThis < this_length; indexThis++)
                {
                    this_line = this.Columns[indexThis];
                    if (this_line.Equals_other(input_line))
                    {
                        colIndexes_list.Add(indexThis);
                    }
                }
            }
            return colIndexes_list.ToArray();
        }

        public int[] Get_colIndexes_of_exact_names(params string[][] names_of_interest)
        {
            List<int> indexes = new List<int>();
            int columns_count = Columns.Count;
            int names_group_length = names_of_interest.Length;
            for (int indexCol = 0; indexCol < columns_count; indexCol++)
            {
                for (int indexG = 0; indexG < names_group_length; indexG++)
                {
                    if (Array_class.Array_order_dependent_equal<string>(Columns[indexCol].Names, names_of_interest[indexG]))
                    {
                        indexes.Add(Columns[indexCol].Index);
                    }
                }
            }
            return indexes.ToArray();
        }

        public int[] Get_colIndexes_that_countain_at_least_one_name(params string[] names)
        {
            List<int> indexes = new List<int>();
            int columns_count = Columns.Count;
            int names_length = names.Length;
            for (int indexCol = 0; indexCol < columns_count; indexCol++)
            {
                for (int indexN = 0; indexN < names_length; indexN++)
                {
                    if (Columns[indexCol].Names.Contains(names[indexN]))
                    {
                        indexes.Add(indexCol);
                        break;
                    }
                }
            }
            return indexes.ToArray();
        }

        public int[] Get_colIndexes_that_countain_all_names(params string[] names)
        {
            List<int> indexes = new List<int>();
            int columns_count = Columns.Count;
            int names_length = names.Length;
            bool contains_all_names;
            for (int indexCol = 0; indexCol < columns_count; indexCol++)
            {
                contains_all_names = true;
                for (int indexN = 0; indexN < names_length; indexN++)
                {
                    if (!Columns[indexCol].Names.Contains(names[indexN]))
                    {
                        contains_all_names = false;
                        break;
                    }
                }
                if (contains_all_names) { indexes.Add(indexCol); }
            }
            return indexes.ToArray();
        }

        public int[] Get_colIndexes_of_first_names(params string[] first_names_ofInterest)
        {
            List<int> indexes = new List<int>();
            int columns_count = Columns.Count;
            for (int i = 0; i < columns_count; i++)
            {
                if (first_names_ofInterest.Contains(Columns[i].Names[0]))
                {
                    indexes.Add(Columns[i].Index);
                }
            }
            return indexes.ToArray();
        }

        public int[] Get_colIndexes_with_indexed_names_containing_at_least_one_inputString(int nameOfInterest_index, params string[] names_inputStrings)
        {
            string names_inputString;
            int names_inputStrings_length = names_inputStrings.Length;
            List<int> indexes = new List<int>();
            int columns_count = Columns.Count;
            for (int i = 0; i < columns_count; i++)
            {
                for (int indexInputString = 0; indexInputString < names_inputStrings_length; indexInputString++)
                {
                    names_inputString = names_inputStrings[indexInputString];
                    if (Columns[i].Names[nameOfInterest_index].IndexOf(names_inputString) != -1)
                    {
                        indexes.Add(Columns[i].Index);
                        break;
                    }
                }
            }
            return indexes.ToArray();
        }
        public int Get_max_index()
        {
            int maxIndex = -1;
            foreach (DE_column_characterization_line_class line in Columns)
            {
                if (line.Index > maxIndex) { maxIndex = line.Index; }
            }
            if (maxIndex != Columns.Count - 1)
            {
                throw new Exception();
            }
            return maxIndex;
        }
        #endregion

        public int Get_columns_count()
        {
            return Columns.Count;
        }

        #region Keep remove columns
        public void Keep_columns(int[] keepColIndexes)
        {
            Small_correctness_check();
            Index_changes_adopted = false;
            Array.Sort(keepColIndexes);
            int col_count = Columns.Count;
            List<DE_column_characterization_line_class> newColumns = new List<DE_column_characterization_line_class>();
            int newIndex = -1;
            foreach (DE_column_characterization_line_class line in Columns)
            {
                if (keepColIndexes.Contains(line.Index))
                {
                    newIndex++;
                    line.Index = newIndex;
                    newColumns.Add(line);
                }
            }
            Columns = newColumns;
        }
        #endregion

        #region Get column labels
        public string[] Get_all_names()
        {
            List<string> names = new List<string>();
            foreach (DE_column_characterization_line_class line in this.Columns)
            {
                if (!string.IsNullOrEmpty(line.Names[0]))
                {
                    names.Add(line.Names[0]);
                }
            }
            return names.ToArray();
        }
        public string Get_complete_column_label(int indexColumn)
        {
            StringBuilder sb = new StringBuilder();
            DE_column_characterization_line_class colChar_line = Columns[indexColumn];
            int names_length = colChar_line.Names.Length;
            for (int indexN = 0; indexN < names_length; indexN++)
            {
                if (sb.Length > 0) { sb.AppendFormat("-"); }
                sb.AppendFormat("{0}", colChar_line.Names[indexN]);
            }
            return sb.ToString();
        }
        #endregion

        #region Order
        private void Reset_indexes()
        {
            int col_count = Columns.Count;
            for (int indexCol = 0; indexCol < col_count; indexCol++)
            {
                Columns[indexCol].IndexOld = Columns[indexCol].Index;
                Columns[indexCol].Index = indexCol;
            }
        }

        public void Order_columns_by_names()
        {
            foreach (DE_column_characterization_line_class col_line in Columns)
            {
                col_line.IndexOld = col_line.Index;
            }
            Index_changes_adopted = false;
            int names_length = Columns[0].Names.Length;
            for (int indexN = names_length - 1; indexN >= 0; indexN--)
            {
                Columns = Columns.OrderBy(l => l.Names[indexN]).ToList();
            }
            Reset_indexes();
        }

        public void Order_columns_by_given_order_of_columnNames(string[] combinedColumnNames)
        {
            if (combinedColumnNames.Distinct().ToArray().Length != combinedColumnNames.Length) { throw new Exception(); }
            int columns_count = Columns.Count;
            if (columns_count != combinedColumnNames.Length) { throw new Exception(); }
            bool[] column_considered = new bool[columns_count];
            Index_changes_adopted = false;
            DE_column_characterization_line_class col_line;
            for (int indexCol = 0; indexCol < columns_count; indexCol++)
            {
                col_line = Columns[indexCol];
                col_line.IndexOld = col_line.Index;
                col_line.Index = -1;
                for (int indexNewCol = 0; indexNewCol < columns_count; indexNewCol++)
                {
                    if (col_line.Combined_names.Equals(combinedColumnNames[indexNewCol]))
                    {
                        col_line.Index = indexNewCol;
                        break;
                    }
                }
                if (col_line.Index == -1) { throw new Exception(); }
            }
        }
        #endregion

        #region Extend
        public int Add_columns_and_count_added_columns(params DE_column_characterization_line_class[] addColchar_lines)
        {
            int col_count = Columns.Count;
            bool line_already_exists;
            int actIndex = Get_max_index();
            int added_cols = 0;
            foreach (DE_column_characterization_line_class line in addColchar_lines)
            {
                line_already_exists = false;
                for (int colIndex = 0; colIndex < col_count; colIndex++)
                {
                    if (line.Equals_other(Columns[colIndex]))
                    {
                        line_already_exists = true;
                        break;
                    }
                }
                if (line_already_exists)
                {
                    //Report_class.Write_error_line("{0}: Timepoint \"{1}\", EntryType \"{2}\" and Name \"{3}\" already exist", typeof(DE_column_characterization_class).Name, line.Timepoint, line.EntryType, line.Name);
                    //Report_class.Write_error_line("{0}: line will not be added", typeof(DE_column_characterization_class).Name, line.Timepoint, line.EntryType);
                }
                else
                {
                    actIndex++;
                    added_cols++;
                    line.Index = actIndex;
                    Columns.Add(line.Deep_copy());
                }
            }
            return added_cols;
        }
        #endregion

        #region Copy
        public DE_column_characterization_class Deep_copy_without_columns()
        {
            DE_column_characterization_class copy = (DE_column_characterization_class)this.MemberwiseClone();
            copy.Columns = new List<DE_column_characterization_line_class>();
            return copy;
        }

        public DE_column_characterization_class Deep_copy()
        {
            DE_column_characterization_class copy = Deep_copy_without_columns();
            foreach (DE_column_characterization_line_class line in this.Columns)
            {
                copy.Columns.Add(line.Deep_copy());
            }
            return copy;
        }
        #endregion
    }

    ///////////////////////////////////////////////////////////////////////////

    class DE_readWriteOptions_class : ReadWriteOptions_base
    {
        #region Fields
        public const char array_delimiter = ';';
        public const char de_array_delimiter = '\t';
        public const char de_symbols_delimiter = ';';

        public static char Array_delimiter { get { return array_delimiter; } }
        public static char DE_array_delimiter { get { return de_array_delimiter; } }
        public static char DE_symbols_delimiter { get { return de_symbols_delimiter; } }
        #endregion

        public DE_readWriteOptions_class(string directory, string file_name, DE_column_characterization_class colChar)
        {
            ReadWriteClass.Create_directory_if_it_does_not_exist(directory);
            File = directory + file_name;
            File_has_headline = true;
            LineDelimiters = new char[] { Global_class.Tab };
            HeadlineDelimiters = new char[] { Global_class.Tab };

            StringBuilder de_headline = new StringBuilder();
            de_headline.AppendFormat("{0}", colChar.Columns[0].Get_full_column_label());
            for (int i = 1; i < colChar.Columns.Count; i++)
            {
                de_headline.AppendFormat("{0}{1}", DE_array_delimiter, colChar.Columns[i].Get_full_column_label());
            }

            Key_propertyNames = new string[] { "ReadWrite_symbols", "ReadWrite_columns"  };
            Key_columnNames = new string[] { "Symbol", de_headline.ToString() };

            Report = ReadWrite_report_enum.Report_main;
        }

    }


    class DE_line_class
    {
        #region Fields
        public string[] Symbols { get; set; }
        public string Description { get; set; }
        public int GeneID { get; set; }
        public DE_columns_line_class[] Columns { get; set; }
        public string Gene_symbol { get { return Symbols[0]; } set { Symbols[0] = value; } }

        #region readWrite def
        public string ReadWrite_symbols
        {
            get { return ReadWriteClass.Get_writeLine_from_array<string>(Symbols, DE_readWriteOptions_class.DE_symbols_delimiter); }
            set { Symbols = ReadWriteClass.Get_array_from_readLine<string>(value, DE_readWriteOptions_class.DE_symbols_delimiter); }
        }

        public string ReadWrite_symbols_with_quote
        {
            get
            {
                int symbols_length = Symbols.Length;
                StringBuilder sb = new StringBuilder();
                for (int indexS = 0; indexS < symbols_length; indexS++)
                {
                    if (indexS != 0) { sb.AppendFormat("{0}", DE_readWriteOptions_class.DE_symbols_delimiter); }
                    sb.AppendFormat("'{0}'", Symbols[indexS]);
                }
                return sb.ToString();
            }
        }

        public string ReadWrite_columns
        {
            get
            {
                StringBuilder write_columns = new StringBuilder();
                int columns_count = Columns.Length;
                write_columns.AppendFormat("{0}", Columns[0].Value);
                for (int indexCol = 1; indexCol < columns_count; indexCol++)
                {
                    write_columns.AppendFormat("{0}{1}", DE_readWriteOptions_class.DE_array_delimiter, Columns[indexCol].Value);
                }
                return write_columns.ToString();
            }
            set
            {
                string[] entries = value.Split(DE_readWriteOptions_class.DE_array_delimiter);
                int entries_length = entries.Length;
                Columns = new DE_columns_line_class[entries_length];
                for (int indexE = 0; indexE < entries_length; indexE++)
                {
                    Columns[indexE] = new DE_columns_line_class();
                    Columns[indexE].Value = Convert.ToDouble(entries[indexE]);
                }
            }
        }
        #endregion
        #endregion

        #region Constructor
        private DE_line_class()
        {
            Symbols = new string[1];
            Columns = new DE_columns_line_class[0];
            Gene_symbol = Global_class.Empty_entry;
            Description = Global_class.Empty_entry;
        }

        public DE_line_class(int col_count)
            : this()
        {
            Extend_columns(col_count);
        }
        #endregion

        #region operations on columns
        public void Extend_columns(int additional_column_count)
        {
            DE_columns_line_class[] add = new DE_columns_line_class[additional_column_count];
            int this_length = this.Columns.Length;
            int add_length = additional_column_count;
            int new_length = this_length + add_length;
            DE_columns_line_class[] new_columns = new DE_columns_line_class[new_length];
            DE_columns_line_class new_column;
            int indexNew = -1;
            for (int indexThis = 0; indexThis < this_length; indexThis++)
            {
                indexNew++;
                new_columns[indexNew] = this.Columns[indexThis];
            }
            for (int indexAdd = 0; indexAdd < add_length; indexAdd++)
            {
                indexNew++;
                new_column = new DE_columns_line_class();
                new_column.Value = 0;
                new_columns[indexNew] = new_column;
            }
            this.Columns = new_columns;
        }

        public void Keep_columns(Dictionary<int, bool> keepCols_dict, int kept_columns_length)
        {
            DE_columns_line_class[] newColumns = new DE_columns_line_class[kept_columns_length];
            int indexNew = -1;
            int col_count = Columns.Length;
            for (int indexCol = 0; indexCol < col_count; indexCol++)
            {
                if (keepCols_dict.ContainsKey(indexCol))
                {
                    indexNew++;
                    newColumns[indexNew] = Columns[indexCol];
                }
            }
            if (indexNew != kept_columns_length - 1) { throw new Exception(); }
            Columns = newColumns;
        }

        public void Keep_columns(int[] keep_columns)
        {
            Dictionary<int, bool> keepCols_dict = new Dictionary<int, bool>();
            foreach (int keep_column in keep_columns)
            {
                keepCols_dict.Add(keep_column, true);
            }
            Keep_columns(keepCols_dict, keep_columns.Length);
        }
        #endregion

        #region Compare equals
        public int Compare_this_de_line_symbols_with_other_de_line_symbols(DE_line_class other)
        {
            int de_line_symbols_length = this.Symbols.Length;
            int data_line_symbols_length = other.Symbols.Length;
            if (de_line_symbols_length != data_line_symbols_length)
            {
                throw new Exception();
            }
            int stringCompare = 1;
            for (int indexS = 0; indexS < de_line_symbols_length; indexS++)
            {
                stringCompare = this.Symbols[indexS].CompareTo(other.Symbols[indexS]);
                if (stringCompare != 0)
                {
                    break;
                }
            }
            return stringCompare;

        }
        #endregion
 
        #region Copy
        public string[] Deep_copy_symbols()
        {
            int symbols_length = this.Symbols.Length;
            string[] copy = new string[symbols_length];
            for (int indexS = 0; indexS < symbols_length; indexS++)
            {
                copy[indexS] = (string)this.Symbols[indexS].Clone();
            }
            return copy;
        }

        public DE_line_class Deep_copy_without_columns()
        {
            DE_line_class newLine = (DE_line_class)this.MemberwiseClone();
            int symbols_length = this.Symbols.Length;
            newLine.Symbols = new string[symbols_length];
            for (int indexS = 0; indexS < symbols_length; indexS++)
            {
                newLine.Symbols[indexS] = (string)this.Symbols[indexS].Clone();
            }
            if (this.Description != null) { newLine.Description = (string)this.Description.Clone(); }
            return newLine;
        }

        public DE_line_class Deep_copy()
        {
            DE_line_class newLine = this.Deep_copy_without_columns();
            int columns_length = this.Columns.Length;
            newLine.Columns = new DE_columns_line_class[columns_length];
            for (int indexC = 0; indexC < columns_length; indexC++)
            {
                newLine.Columns[indexC] = this.Columns[indexC].Deep_copy();
            }
            return newLine;
        }
        #endregion
    }

    class DE_class
    {
        #region Fields
        public List<DE_line_class> DE { get; set; }
        public int Symbols_length { get; set; }
        public DE_line_class ColSums { get; set; }
        public DE_line_class ColSums_abs { get; set; }
        public DE_line_class ColMeans { get; set; }
        public DE_line_class SampleColSDs { get; set; }
        public DE_line_class PopulationColSDs { get; set; }
        public DE_line_class Non_zero_entry_count_in_cols { get; set; }
        public DE_column_characterization_class ColChar { get; set; }
        #endregion

        #region Constructor
        public DE_class()
        {
            DE = new List<DE_line_class>();
            ColChar = new DE_column_characterization_class();
            ColMeans = new DE_line_class(0);
            ColSums = new DE_line_class(0);
            ColSums_abs = new DE_line_class(0);
            SampleColSDs = new DE_line_class(0);
            PopulationColSDs = new DE_line_class(0);
        }
        #endregion

        #region Check
        public void Check_for_duplicates()
        {
            this.Order_by_symbol();
            int de_length = this.DE.Count;
            DE_line_class de_line;
            DE_line_class previous_de_line;
            for (int indexDE = 1; indexDE < de_length; indexDE++)
            {
                previous_de_line = DE[indexDE - 1];
                de_line = DE[indexDE];
                if (previous_de_line.Symbols[0].Equals(de_line.Symbols[0]))
                {
                    throw new Exception();
                }
            }

        }
        #endregion

        #region Order
        public void Order_by_symbol()
        {
            if (DE.Count > 0)
            {
                int symbols_length = DE[0].Symbols.Length;
                for (int indexS = symbols_length - 1; indexS >= 0; indexS--)
                {
                    DE = DE.OrderBy(l => l.Symbols[indexS]).ToList();
                }
            }
        }

        public void Order_by_descending_column_values(int indexCol)
        {
            DE = DE.OrderByDescending(l => l.Columns[indexCol].Value).ToList();
        }

        public void Order_by_column_values(int indexCol)
        {
            DE = DE.OrderBy(l => l.Columns[indexCol].Value).ToList();
        }
        #endregion


        #region Fill with data
        private Fill_de_line_class[] Set_combined_symbols_and_combined_names(Fill_de_line_class[] inputData)
        {
            int inputData_length = inputData.Length;
            Fill_de_line_class fill_de_line;
            StringBuilder sb = new StringBuilder();
            for (int indexInput = 0; indexInput < inputData_length; indexInput++)
            {
                fill_de_line = inputData[indexInput];
                fill_de_line.Set_combined_symbols_and_names();
            }
            return inputData;
        }

        public int Fill_with_data_alternatively(Fill_de_line_class[] inputData)
        {
            System.Diagnostics.Stopwatch stopwatch = new System.Diagnostics.Stopwatch();
            stopwatch.Start();

            inputData = Set_combined_symbols_and_combined_names(inputData);

            List<DE_column_characterization_line_class> colChars = new List<DE_column_characterization_line_class>();
            DE_column_characterization_line_class new_colChar_line;
            Fill_de_line_class fill_de_line;
            int inputData_length = inputData.Length;

            int indexCol = -1;
            int currentCol = -1;

            inputData = Fill_de_line_class.Order_by_combinedNames(inputData);
            for (int indexInput = 0; indexInput < inputData_length; indexInput++)
            {
                fill_de_line = inputData[indexInput];
                if ((indexInput == 0)
                    || (!fill_de_line.Combined_names.Equals(inputData[indexInput - 1].Combined_names)))
                {
                    new_colChar_line = new DE_column_characterization_line_class();
                    new_colChar_line.Names = Array_class.Deep_copy_string_array(fill_de_line.Names_for_de);
                    indexCol++;
                    new_colChar_line.Index = indexCol;
                    colChars.Add(new_colChar_line);
                    currentCol++;
                }
                fill_de_line.Column = currentCol;
            }
            this.ColChar.Columns = colChars;
            int col_count = this.ColChar.Columns.Count;
            indexCol = 0;

            Dictionary<string, int> symbol_rowIndex_dict = new Dictionary<string, int>();

            DE_line_class new_de_line = new DE_line_class(col_count);
            List<DE_line_class> new_de_list = new List<DE_line_class>();
            int current_lastIndex_of_new_de_list = -1;
            int row_of_symbol;
            for (int indexInput = 0; indexInput < inputData_length; indexInput++)
            {
                fill_de_line = inputData[indexInput];
                if (!symbol_rowIndex_dict.ContainsKey(fill_de_line.Combined_symbols))
                {
                    new_de_line = new DE_line_class(col_count);
                    new_de_line.Symbols = Array_class.Deep_copy_string_array(fill_de_line.Symbols_for_de);
                    current_lastIndex_of_new_de_list++;
                    symbol_rowIndex_dict.Add(fill_de_line.Combined_symbols, current_lastIndex_of_new_de_list);
                    new_de_list.Add(new_de_line);
                }
                row_of_symbol = symbol_rowIndex_dict[fill_de_line.Combined_symbols];
                new_de_list[row_of_symbol].Columns[fill_de_line.Column].Value = fill_de_line.Value_for_de;
            }
            this.DE = new_de_list;
            Remove_empty_rows_and_columns();
            stopwatch.Stop();
            return col_count;
        }
        #endregion

        #region Modify symbols
        public void Set_symbols_to_upper_case_letters()
        {
            foreach (DE_line_class line in DE)
            {
                line.Gene_symbol = line.Gene_symbol.ToUpper();
            }

        }
        #endregion

        public int[] Get_non_zero_counts_of_each_column_in_indexOrder()
        {
            int col_count = ColChar.Columns.Count;
            DE_line_class de_line;
            DE_columns_line_class column_line;
            int de_count = DE.Count;
            int[] non_zero_counts = new int[col_count];
            for (int indexDE = 0; indexDE < de_count; indexDE++)
            {
                de_line = DE[indexDE];
                col_count = de_line.Columns.Length;
                for (int indexCol = 0; indexCol < col_count; indexCol++)
                {
                    column_line = de_line.Columns[indexCol];
                    if (column_line.Value != 0)
                    {
                        non_zero_counts[indexCol]++;
                    }
                }
            }
            return non_zero_counts;
        }

        #region Keep and remove
        public void Remove_empty_rows_and_columns()
        {
            int column_count = ColChar.Columns.Count;
            List<int> keepColumns = new List<int>();
            bool keep_row;
            List<DE_line_class> newDE = new List<DE_line_class>();
            int de_length = DE.Count;
            DE_line_class line;
            for (int indexDE = 0; indexDE < de_length; indexDE++)
            {
                line = DE[indexDE];
                keep_row = false;
                for (int indexCol = 0; indexCol < column_count; indexCol++)
                {
                    if (line.Columns[indexCol].Value != 0)
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
            keepColumns = keepColumns.Distinct().OrderBy(l => l).ToList();
            int keepColumns_length = keepColumns.Count;
            Dictionary<int, bool> keepCols_dict = new Dictionary<int, bool>();
            foreach (int keepColumn in keepColumns)
            {
                keepCols_dict.Add(keepColumn, true);
            }

            if (keepColumns.Count() < ColChar.Columns.Count)
            {
                ColChar.Keep_columns(keepColumns.ToArray());
                de_length = newDE.Count;
                for (int indexDE = 0; indexDE < de_length; indexDE++)
                {
                    line = newDE[indexDE];
                    line.Keep_columns(keepCols_dict, keepColumns_length);
                    newDE[indexDE] = line;
                }
            }
            DE = newDE;
            this.ColChar.Index_changes_adopted = true;
        }
        public void Keep_only_stated_symbols(string[] keepSymbols)
        {
            keepSymbols = keepSymbols.Distinct().OrderBy(l => l).ToArray();
            int keep_length = keepSymbols.Length;
            Order_by_symbol();
            int de_symbol_count = DE.Count;
            int indexDE = 0;
            int stringCompare = 0;
            DE_line_class de_line;
            List<DE_line_class> newDE = new List<DE_line_class>();
            for (int indexKeep = 0; indexKeep < keep_length; indexKeep++)
            {
                stringCompare = -2;
                while ((indexDE < de_symbol_count) && (stringCompare < 0))
                {
                    de_line = DE[indexDE];
                    stringCompare = de_line.Gene_symbol.CompareTo(keepSymbols[indexKeep]);
                    if (stringCompare < 0)
                    {
                        indexDE++;
                    }
                    else if (stringCompare == 0)
                    {
                        newDE.Add(de_line);
                    }
                }
            }
            DE = newDE;
            Remove_empty_rows_and_columns();
        }
        #endregion

        #region Write Deep copy
        public void Write_file(string directory, string file_name)
        {
            DE_readWriteOptions_class options = new DE_readWriteOptions_class(directory, file_name, ColChar);
            ReadWriteClass.WriteData(DE, options);
        }
        public DE_class Deep_copy_without_DE_lines()
        {
            DE_class copy = (DE_class)this.MemberwiseClone();
            copy.ColChar = this.ColChar.Deep_copy();
            copy.DE = new List<DE_line_class>();
            return copy;
        }
        public DE_class Deep_copy()
        {
            DE_class copy = Deep_copy_without_DE_lines();
            copy.ColSums = this.ColSums.Deep_copy();
            copy.ColSums_abs = this.ColSums_abs.Deep_copy();
            copy.ColMeans = this.ColMeans.Deep_copy();
            copy.SampleColSDs = this.SampleColSDs.Deep_copy();
            copy.PopulationColSDs = this.PopulationColSDs.Deep_copy();
            foreach (DE_line_class line in DE)
            {
                copy.DE.Add(line.Deep_copy());
            }
            return copy;
        }
        #endregion
    }

}
