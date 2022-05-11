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

namespace Statistic
{
    class Fisher_exact_test_class
    {
        private double[] log_factorials;
        private int max_size;
 
        //         a     b
        //         c     d


        public Fisher_exact_test_class(int input_max_size, bool report)
        {
            max_size = input_max_size;
            log_factorials = new double[max_size+1];
            log_factorials[0] = 0;
            for (int i = 1; i < max_size+1; i++)
            {
                log_factorials[i] = log_factorials[i - 1] + Math.Log(i);
            }
        }

        private bool Check_if_n_not_larger_than_max_size(int a, int b, int c, int d)
        {
            bool smaller = true;
            int n = a+b+c+d;
            if (n > max_size + 1)
            {
                throw new Exception();
            }
            return smaller;
        }

        private double Get_specific_log_p_value(int a, int b, int c, int d)
        {
            int n = a + b + c + d;
            double log_p = log_factorials[a + b] + log_factorials[c + d] + log_factorials[a + c] + log_factorials[b + d] - log_factorials[n] - log_factorials[a] - log_factorials[b] - log_factorials[c] - log_factorials[d];
            return log_p;
        }

        private double Get_specific_p_value(int a, int b, int c, int d)
        {
            double log_p = Get_specific_log_p_value(a, b, c, d);
            return Math.Exp(log_p);
        }

        public double Get_rightTailed_p_value(int a, int b, int c, int d)
        {
            double p;
            if (Check_if_n_not_larger_than_max_size(a, b, c, d))
            {
                p = Get_specific_p_value(a, b, c, d);
                int min = (c < b) ? c : b;
                for (int i = 0; i < min; i++)
                {
                    p += Get_specific_p_value(++a, --b, --c, ++d);
                }
            }
            else { p = -1; };
            if (p > 1) { p = 1; }
            return p;
        }

        public double Get_rightTailed_log10_p_value_based_on_log(int a, int b, int c, int d)
        {
            double p = 0;
            double specific_log_p_value = 1; 
            if (Check_if_n_not_larger_than_max_size(a, b, c, d))
            {
                int min = (c < b) ? c : b;
                double[] specific_log_p_values = new double[min+1];
                specific_log_p_values[0] = Get_specific_log_p_value(a, b, c, d);
                for (int i = 1; i < min+1; i++)
                {
                    specific_log_p_values[i] = Get_specific_log_p_value(++a, --b, --c, ++d);
                }
                specific_log_p_values = specific_log_p_values.OrderByDescending(l => l).ToArray();
                int specific_log_p_values_length = specific_log_p_values.Length;
                int indexMiddleValue = (int)Math.Round((double)specific_log_p_values_length / (double)2);

                double min_specific_log_p_value = specific_log_p_values[0];
                double p_add;
                for (int i = 0; i < min+1; i++)
                {
                    specific_log_p_value = specific_log_p_values[i];
                    specific_log_p_value -= min_specific_log_p_value;
                    p_add = Math.Exp(specific_log_p_value);
                    if (!Double.IsInfinity(p_add))
                    {
                        p += p_add;
                    }
                    else
                    {
                      //  throw new Exception();
                    }
                }
                double log_10 = Math.Log(10);
                specific_log_p_value = (1.0/log_10) * (Math.Log(p) + min_specific_log_p_value);
            }
            else { p = -1; };
            return (double)specific_log_p_value;
        }

        public double Get_leftTailed_log10_p_value_based_on_log(int a, int b, int c, int d)
        {
            double p = 0;
            double specific_log_p_value = 1;
            if (Check_if_n_not_larger_than_max_size(a, b, c, d))
            {
                int min = (a < d) ? a : d;
                double[] specific_log_p_values = new double[min + 1];
                specific_log_p_values[0] = Get_specific_log_p_value(a, b, c, d);
                for (int i = 1; i < min+1; i++)
                {
                    specific_log_p_values[i] = Get_specific_log_p_value(--a, ++b, ++c, --d);
                }
                specific_log_p_values = specific_log_p_values.OrderByDescending(l => l).ToArray();
                double min_specific_log_p_value = specific_log_p_values[0];
                double p_add;
                for (int i = 0; i < min + 1; i++)
                {
                    specific_log_p_value = specific_log_p_values[i];
                    specific_log_p_value -= min_specific_log_p_value;
                    p_add = Math.Exp(specific_log_p_value);
                    if (!Double.IsInfinity(p_add))
                    {
                        p += p_add;
                    }
                    else
                    {
                        //  throw new Exception();
                    }
                }
                double log_10 = Math.Log(10);
                specific_log_p_value = (1.0 / log_10) * (Math.Log(p) + min_specific_log_p_value);
            }
            else { p = -1; };
            return (double)specific_log_p_value;
        }

        public double Get_leftTailed_p_value(int a, int b, int c, int d)
        {
            double p;
            if (Check_if_n_not_larger_than_max_size(a, b, c, d))
            {
                p = Get_specific_p_value(a, b, c, d);
                int min = (a < d) ? a : d;
                for (int i = 0; i < min; i++)
                {
                    p += Get_specific_p_value(--a, ++b, ++c, --d);
                }
            }
            else { p = -1; };
            if (p > 1) { p = 1; }
            return p;
        }
    }

    class Overlap_class
    {
        public static bool Check_if_identical_string_arrays(string[] array1, string[] array2)
        {
            int array_length = array1.Length;
            if (array_length!=array2.Length) { throw new Exception(); }
            array1 = array1.OrderBy(l => l).ToArray();
            array2 = array2.OrderBy(l => l).ToArray();
            for (int indexA=0; indexA<array_length;indexA++)
            {
                if (!array1[indexA].Equals(array2[indexA]))
                {
                    throw new Exception();
                }
            }
            return true;
        }

        public static string[] Get_intersection(string[] list1, string[] list2)
        {
            list1 = list1.Distinct().OrderBy(l => l).ToArray();
            list2 = list2.Distinct().OrderBy(l => l).ToArray();
            int list1_length = list1.Length;
            int list2_length = list2.Length;
            int index1 = 0;
            int index2 = 0;
            int stringCompare;
            List<string> intersection = new List<string>();
            while ((index1 < list1_length) && (index2 < list2_length))
            {
                stringCompare = list2[index2].CompareTo(list1[index1]);
                if (stringCompare < 0) { index2++; }
                else if (stringCompare > 0) { index1++; }
                else
                {
                    intersection.Add(list1[index1]);
                    index1++;
                    index2++;
                }
            }
            return intersection.ToArray();
        }

        public static string[] Get_intersection(string[][] lists)
        {
            int lists_length = lists.Length;
            if (lists_length==1)
            {
                throw new Exception();
            }
            string[] intersection = Get_intersection(lists[0], lists[1]);
            for (int indexList=2; indexList<lists_length;indexList++)
            {
                intersection = Get_intersection(intersection, lists[indexList]);
            }
            return intersection;
        }

        public static string[] Get_union(string[] list1, params string[] list2)
        {
            string[] union = new string[0];
            if ((list1 != null) && (list2 != null))
            {
                list1 = list1.Distinct().OrderBy(l => l).ToArray();
                list2 = list2.Distinct().OrderBy(l => l).ToArray();
                int list1_length = list1.Length;
                int list2_length = list2.Length;
                int index1 = 0;
                int index2 = 0;
                int stringCompare;
                List<string> union_list = new List<string>();
                while ((index1 < list1_length) || (index2 < list2_length))
                {
                    if ((index1 < list1_length) && (index2 < list2_length))
                    {
                        stringCompare = list2[index2].CompareTo(list1[index1]);
                        if (stringCompare < 0)
                        {
                            union_list.Add(list2[index2]);
                            index2++;
                        }
                        else if (stringCompare > 0)
                        {
                            union_list.Add(list1[index1]);
                            index1++;
                        }
                        else
                        {
                            union_list.Add(list1[index1]);
                            index1++;
                            index2++;
                        }
                    }
                    else if (index1 < list1_length)
                    {
                        union_list.Add(list1[index1]);
                        index1++;
                    }
                    else if (index2 < list2_length)
                    {
                        union_list.Add(list2[index2]);
                        index2++;
                    }
                }
                union = union_list.ToArray();
            }
            else if (list1 != null)
            {
                union = list1;
            }
            else if (list2 != null)
            {
                union = list2;
            }
            return union;
        }

        public static string[] Get_union_of_distinct_ordered_symbols(params string[][] lists)
        {
            List<string> union = new List<string>();
            foreach (string[] list in lists)
            {
                union.AddRange(list);
            }
            return union.Distinct().OrderBy(l => l).ToArray();
        }

        public static string[] Add_to_deep_copy_of_list2_to_end_of_deep_copy_of_list1(string[] list1, params string[] list2)
        {
            int list1_count = list1.Length;
            int list2_count = list2.Length;
            int new_list_count = list1_count + list2_count;
            string[] new_list = new string[new_list_count];
            int indexNew = -1;
            for (int indexList1 = 0; indexList1 < list1_count; indexList1++)
            {
                indexNew++;
                new_list[indexNew] = (string)list1[indexList1].Clone();
            }
            for (int indexList2 = 0; indexList2 < list2_count; indexList2++)
            {
                indexNew++;
                new_list[indexNew] = (string)list2[indexList2].Clone();
            }
            return new_list;
        }

        public static string[] Get_part_of_list1_but_not_of_list2(string[] list1, string[] list2)
        {
            list1 = list1.Distinct().OrderBy(l => l).ToArray();
            list2 = list2.Distinct().OrderBy(l => l).ToArray();
            List<string> not = new List<string>();
            int list1_length = list1.Length;
            int list2_length = list2.Length;
            int index2=0;
            int stringCompare;
            for (int index1 = 0; index1 < list1_length; index1++)
            {
                stringCompare = -2;
                while ((index2 < list2_length) && (stringCompare < 0))
                {
                    stringCompare = list2[index2].CompareTo(list1[index1]);
                    if (stringCompare < 0) { index2++; }
                }
                if ((stringCompare > 0) || (index2 == list2_length))
                {
                    not.Add(list1[index1]);
                }
            }
            return not.ToArray();
        }

        public static List<string> Get_part_of_list1_but_not_of_list2(List<string> list1, List<string> list2)
        {
            list1 = list1.Distinct().OrderBy(l => l).ToList();
            list2 = list2.Distinct().OrderBy(l => l).ToList();
            List<string> not = new List<string>();
            int list1_length = list1.Count;
            int list2_length = list2.Count;
            int index2 = 0;
            int stringCompare;
            for (int index1 = 0; index1 < list1_length; index1++)
            {
                stringCompare = -2;
                while ((index2 < list2_length) && (stringCompare < 0))
                {
                    stringCompare = list2[index2].CompareTo(list1[index1]);
                    if (stringCompare < 0) { index2++; }
                }
                if ((stringCompare > 0) || (index2 == list2_length))
                {
                    not.Add(list1[index1]);
                }
            }
            return not;
        }


        public static int[] Get_intersection(int[] list1, int[] list2) 
        {
            list1 = list1.Distinct().OrderBy(l => l).ToArray();
            list2 = list2.Distinct().OrderBy(l => l).ToArray();
            int list1_length = list1.Length;
            int list2_length = list2.Length;
            int index1 = 0;
            int index2 = 0;
            int stringCompare;
            List<int> intersection = new List<int>();
            while ((index1 < list1_length) && (index2 < list2_length))
            {
                stringCompare = list2[index2].CompareTo(list1[index1]);
                if (stringCompare < 0) { index2++; }
                else if (stringCompare > 0) { index1++; }
                else
                {
                    intersection.Add(list1[index1]);
                    index1++;
                    index2++;
                }
            }
            return intersection.ToArray();
        }

        public static int[] Get_union(int[] list1, int[] list2)
        {
            int[] union = new int[0];
            if ((list1 != null) && (list2 != null))
            {
                list1 = list1.Distinct().OrderBy(l => l).ToArray();
                list2 = list2.Distinct().OrderBy(l => l).ToArray();
                int list1_length = list1.Length;
                int list2_length = list2.Length;
                int index1 = 0;
                int index2 = 0;
                int intCompare;
                List<int> union_list = new List<int>();
                while ((index1 < list1_length) || (index2 < list2_length))
                {
                    if ((index1 < list1_length) && (index2 < list2_length))
                    {
                        intCompare = list2[index2] - list1[index1];
                        if (intCompare < 0)
                        {
                            union_list.Add(list2[index2]);
                            index2++;
                        }
                        else if (intCompare > 0)
                        {
                            union_list.Add(list1[index1]);
                            index1++;
                        }
                        else
                        {
                            union_list.Add(list1[index1]);
                            index1++;
                            index2++;
                        }
                    }
                    else if (index1 < list1_length)
                    {
                        union_list.Add(list1[index1]);
                        index1++;
                    }
                    else if (index2 < list2_length)
                    {
                        union_list.Add(list2[index2]);
                        index2++;
                    }
                }
                union = union_list.ToArray();
            }
            else if (list1 != null)
            {
                union = list1;
            }
            else if (list2 != null)
            {
                union = list2;
            }
            return union;
        }

        public static T[] Get_union<T>(T[] list1, params T[] list2) where T:IComparable
        {
            T[] union = new T[0];
            if ((list1 != null) && (list2 != null))
            {
                list1 = list1.Distinct().OrderBy(l => l).ToArray();
                list2 = list2.Distinct().OrderBy(l => l).ToArray();
                int list1_length = list1.Length;
                int list2_length = list2.Length;
                int index1 = 0;
                int index2 = 0;
                int intCompare;
                List<T> union_list = new List<T>();
                while ((index1 < list1_length) || (index2 < list2_length))
                {
                    if ((index1 < list1_length) && (index2 < list2_length))
                    {
                        intCompare = list2[index2].CompareTo(list1[index1]);
                        if (intCompare < 0)
                        {
                            union_list.Add(list2[index2]);
                            index2++;
                        }
                        else if (intCompare > 0)
                        {
                            union_list.Add(list1[index1]);
                            index1++;
                        }
                        else
                        {
                            union_list.Add(list1[index1]);
                            index1++;
                            index2++;
                        }
                    }
                    else if (index1 < list1_length)
                    {
                        union_list.Add(list1[index1]);
                        index1++;
                    }
                    else if (index2 < list2_length)
                    {
                        union_list.Add(list2[index2]);
                        index2++;
                    }
                }
                union = union_list.ToArray();
            }
            else if (list1 != null)
            {
                union = list1;
            }
            else if (list2 != null)
            {
                union = list2;
            }
            return union;
        }

        public static T[] Get_intersection<T>(T[] list1, T[] list2) where T:IComparable
        {
            list1 = list1.Distinct().OrderBy(l => l).ToArray();
            list2 = list2.Distinct().OrderBy(l => l).ToArray();
            int list1_length = list1.Length;
            int list2_length = list2.Length;
            int index1 = 0;
            int index2 = 0;
            int stringCompare;
            List<T> intersection = new List<T>();
            while ((index1 < list1_length) && (index2 < list2_length))
            {
                stringCompare = list2[index2].CompareTo(list1[index1]);
                if (stringCompare < 0) { index2++; }
                else if (stringCompare > 0) { index1++; }
                else
                {
                    intersection.Add(list1[index1]);
                    index1++;
                    index2++;
                }
            }
            return intersection.ToArray();
        }

        public static T[] Get_part_of_list1_but_not_of_list2<T>(T[] list1, T[] list2) where T:IComparable
        {
            list1 = list1.Distinct().OrderBy(l => l).ToArray();
            list2 = list2.Distinct().OrderBy(l => l).ToArray();
            List<T> not = new List<T>();
            int list1_length = list1.Length;
            int list2_length = list2.Length;
            int index2 = 0;
            int stringCompare;
            for (int index1 = 0; index1 < list1_length; index1++)
            {
                stringCompare = -2;
                while ((index2 < list2_length) && (stringCompare < 0))
                {
                    stringCompare = list2[index2].CompareTo(list1[index1]);
                    if (stringCompare < 0) { index2++; }
                }
                if ((stringCompare > 0) || (index2 == list2_length))
                {
                    not.Add(list1[index1]);
                }
            }
            return not.ToArray();
        }
    }

    class Math_class
    {
        private static string Get_hexadecimal_sign(int number)
        {
            string sign = "no value";
            switch (number)
            {
                case 0:
                case 1:
                case 2:
                case 3:
                case 4:
                case 5:
                case 6:
                case 7:
                case 8:
                case 9:
                    sign = number.ToString();
                    break;
                case 10:
                    sign = "A";
                    break;
                case 11:
                    sign = "B";
                    break;
                case 12:
                    sign = "C";
                    break;
                case 13:
                    sign = "D";
                    break;
                case 14:
                    sign = "E";
                    break;
                case 15:
                    sign = "F";
                    break;
                default:
                    throw new Exception();
            }
            return sign;
        }

        public static string Convert_into_two_digit_hexadecimal(int number)
        {
            if ((number > 255) || (number < 0))
            {
                throw new Exception();
            }
            else
            {
                int multiples_of_16 = (int)Math.Floor((double)number / (double)16);
                int modulus = number % 16;
                return Get_hexadecimal_sign(multiples_of_16) + Get_hexadecimal_sign(modulus);
            }
        }

        public static string Get_hexadecimal_code(int red, int green, int blue)
        {
            return "#" + Convert_into_two_digit_hexadecimal(red) + Convert_into_two_digit_hexadecimal(green) + Convert_into_two_digit_hexadecimal(blue);
        }

        public static string Get_hexadecimal_light_blue()
        {
            return Get_hexadecimal_code(51,153,255);
        }

        public static string Get_hexadecimal_light_red()
        {
            return Get_hexadecimal_code(230, 0, 0);
        }

        public static string Get_hexadecimal_white()
        {
            return Get_hexadecimal_code(255, 255, 255);
        }

        public static string Get_hexadecimal_light_green()
        {
            return Get_hexadecimal_code(134, 196, 64);
        }

        public static string Get_hexadecimal_bright_green()
        {
            return Get_hexadecimal_code(0, 255, 0);
        }

        public static string Get_hexadecimal_dark_green()
        {
            return Get_hexadecimal_code(0, 100, 0);
        }

        public static string Get_hexadecimal_black()
        {
            return Get_hexadecimal_code(0, 0, 0);
        }

        public static string Get_hexadecimal_light_gray()
        {
            return Get_hexadecimal_code(211, 211, 211);
        }

        public static string Get_hexadecimal_orange()
        {
            return Get_hexadecimal_code(255, 128, 0);
        }

        public static double Get_median(double[] values)
        {
            values = values.OrderBy(l => l).ToArray();
            int values_length = values.Length;
            if (values_length == 0)
            {
                throw new InvalidOperationException("Empty collection");
            }
            else if (values_length % 2 == 0)
            {
                // count is even, average two middle elements
                double a = values[(values_length / 2) - 1];
                double b = values[(values_length / 2)];
                return (a + b) / 2;
            }
            else if (values_length == 1)
            {
                return values[0];
            }
            else
            {
                // count is odd, return the middle element
                double return_value = values[(int)Math.Floor((double)values_length / (double)2)];
                return return_value;
            }
        }

        public static float Get_median(float[] values)
        {
            values = values.OrderBy(l => l).ToArray();
            int values_length = values.Length;
            if (values_length == 0)
            {
                throw new InvalidOperationException("Empty collection");
            }
            else if (values_length % 2 == 0)
            {
                // count is even, average two middle elements
                float a = values[(values_length / 2) - 1];
                float b = values[(values_length / 2)];
                return (a + b) / 2;
            }
            else if (values_length == 1)
            {
                return values[0];
            }
            else
            {
                // count is odd, return the middle element
                float return_value = values[(int)Math.Floor((float)values_length / (float)2)];
                return return_value;
            }
        }

        public static float Get_average(float[] values)
        {
            int values_length = values.Length;
            float sum = 0;
            for (int indexV = 0; indexV < values_length; indexV++)
            {
                sum += values[indexV];
            }
            return sum / (float)values_length;
        }

        public static double Get_average(double[] values)
        {
            int values_length = values.Length;
            double sum = 0;
            for (int indexV = 0; indexV < values_length; indexV++)
            {
                sum += values[indexV];
            }
            return sum / (double)values_length;
        }

        public static float Get_geometric_mean(float[] values)
        {
            int values_length = values.Length;
            float product = 1;
            for (int indexV = 0; indexV < values_length; indexV++)
            {
                product *= values[indexV];
            }
            float geometric_mean = -1;
            checked { geometric_mean = (float)Math.Pow(product, 1.0 / (double)values_length); }
            return geometric_mean;
        }

        public static double Get_geometric_mean(double[] values)
        {
            int values_length = values.Length;
            double product = 1;
            for (int indexV = 0; indexV < values_length; indexV++)
            {
                product *= values[indexV];
            }
            double geometric_mean = Math.Pow(product, 1.0 / (double)values_length);
            return geometric_mean;
        }

        public static double Get_harmonic_mean(double[] values)
        {
            int values_length = values.Length;
            double sum_of_inverse_values = 0;
            for (int indexV = 0; indexV < values_length; indexV++)
            {
                sum_of_inverse_values += (double)1/values[indexV];
            }
            double harmonic_mean = (double)values_length / sum_of_inverse_values;
            return harmonic_mean;
        }

        public static void Get_mean_and_population_sd(float[] values, out float mean, out float sd)
        {
            int values_length = values.Length;
            double sum = 0;
            double sum_of_squares = 0;
            for (int indexV = 0; indexV < values_length; indexV++)
            {
                sum += values[indexV];
                sum_of_squares += Math.Pow(values[indexV], 2);
            }
            double mean_double = sum / (double)values_length;
            sd = (float)Math.Sqrt(sum_of_squares / (float)values_length - Math.Pow(mean_double, 2));
            mean = (float)mean_double;
            // if (  (float.IsNaN(sd))
            //     ||(float.IsInfinity(sd)))
            // {
            //     throw new Exception();
            // }
        }

        public static void Get_mean_and_population_sd(double[] values, out double mean, out double sd)
        {
            int values_length = values.Length;
            double sum = 0;
            double sum_of_squares = 0;
            for (int indexV = 0; indexV < values_length; indexV++)
            {
                sum += values[indexV];
                sum_of_squares += Math.Pow(values[indexV], 2);
            }
            double mean_double = sum / (double)values_length;
            sd = (double)Math.Sqrt(sum_of_squares / (double)values_length - Math.Pow(mean_double, 2));
            mean = (double)mean_double;
            // if (  (float.IsNaN(sd))
            //     ||(float.IsInfinity(sd)))
            // {
            //     throw new Exception();
            // }
        }

        public static void Get_mean_and_sample_sd(float[] values, out float mean, out float sd)
        {
            int values_length = values.Length;
            Get_mean_and_population_sd(values, out mean, out sd);
            if (values_length ==1) { sd = 0; }
            else
            {
                sd = sd * (float)Math.Sqrt((float)(values_length) / (float)(values_length - 1));
            }
        }

        public static void Get_mean_and_sample_sd(double[] values, out double mean, out double sd)
        {
            double population_sd;
            int values_length = values.Length;
            Get_mean_and_population_sd(values, out mean, out population_sd);
            if (values_length == 1) { sd = 0; }
            else
            {
                sd = (double)population_sd * Math.Sqrt((double)(values_length) / (double)(values_length - 1));
            }
        }

        public static void Get_mean_and_sample_sd_using_long(double[] values, out double mean, out double sd)
        {
            int values_length = values.Length;
            mean = 0;
            for (int indexV = 0; indexV < values_length; indexV++)
            {
                mean += values[indexV] / (double)values_length;
                if (double.IsNaN(mean))
                {
                    Console.WriteLine();
                }
            }
            double sum_of_deviations = 0;
            for (int indexV = 0; indexV < values_length; indexV++)
            {
                sum_of_deviations += Math.Pow((long)values[indexV]-mean,2);
            }
            double sd_double = Math.Sqrt((double)sum_of_deviations / (values_length-1));
            if (values_length == 1) { sd = 0; }
            else
            {
                sd = (double)(sd_double * Math.Sqrt((double)(values_length) / (double)values_length - 1));
            }
            throw new Exception();
        }


        public static void Get_max_min_of_array(float[] array, out float max, out float min)
        {
            int array_length = array.Length;
            max = -1;
            min = -1;
            float array_entry;
            for (int indexA = 0; indexA < array_length; indexA++)
            {
                array_entry = array[indexA];
                if ((max == -1) || (array_entry > max))
                {
                    max = array_entry;
                }
                if ((min == -1) || (array_entry < min))
                {
                    min = array_entry;
                }
            }
        }

        public static void Get_max_min_of_array(double[] array, out double max, out double min)
        {
            int array_length = array.Length;
            max = -1;
            min = -1;
            double array_entry;
            for (int indexA = 0; indexA < array_length; indexA++)
            {
                array_entry = array[indexA];
                if ((max == -1) || (array_entry > max))
                {
                    max = array_entry;
                }
                if ((min == -1) || (array_entry < min))
                {
                    min = array_entry;
                }
            }
        }

        public static void Get_max_min_of_array(int[] array, out int max, out int min)
        {
            int array_length = array.Length;
            max = -1;
            min = -1;
            int array_entry;
            for (int indexA = 0; indexA < array_length; indexA++)
            {
                array_entry = array[indexA];
                if ((max == -1) || (array_entry > max))
                {
                    max = array_entry;
                }
                if ((min == -1) || (array_entry < min))
                {
                    min = array_entry;
                }
            }
        }

    }
}
