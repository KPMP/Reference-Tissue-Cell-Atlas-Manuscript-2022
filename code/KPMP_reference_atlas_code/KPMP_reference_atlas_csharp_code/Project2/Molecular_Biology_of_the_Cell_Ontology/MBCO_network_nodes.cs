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

Please see http://www.mbc-ontology.org/ or https://github.com/SBCNY/Molecular-Biology-of-the-Cell for a windows application offering the same and additional functionalities as the script below.
Please acknowledge the MBC Ontology in your publications by citing the following reference:
Jens Hansen, David Meretzky, Simeneh Woldesenbet, Gustavo Stolovitzky, Ravi Iyengar: 
A flexible ontology for inference of emergent whole cell function from relationships between subcellular processes
Sci Rep. 2017 Dec18th
*/

using System;
using System.Collections.Generic;
using System.Linq;
using Enumerations;
using System.Drawing;

namespace MBCO_network_and_nodes
{
    public enum Obo_interaction_type_enum { E_m_p_t_y, Is_a, Part_of }
    public enum Ontology_namespace_enum { E_m_p_t_y, Disease_ontology }
    public enum Ontology_direction_enum { E_m_p_t_y, Parent_child, Child_parent }

    class Color_specification_line_class
    {
        public Color Fill_color { get; set; }
        public float Size { get; set; }

        public Color_specification_line_class()
        {
            Size = 1;
        }

        public Color_specification_line_class Deep_copy()
        {
            Color_specification_line_class copy = (Color_specification_line_class)this.MemberwiseClone();
            return copy;
        }
    }

    class Add_node_line_class
    {
        public string Name { get; set; }
        public Ontology_type_enum Ontology_type { get; set; }
        public Ontology_namespace_enum Ontology_namespace { get; set; }

        public static Add_node_line_class[] Order_by_standard_way(Add_node_line_class[] add_nodes)
        {
            add_nodes = add_nodes.OrderBy(l => l.Name).ThenBy(l => l.Ontology_namespace).ThenBy(l => l.Ontology_type).ToArray();
            return add_nodes;
        }

        public bool Equal_in_standard_way(Add_node_line_class other)
        {
            bool equal = ( (this.Name.Equals(other.Name))
                          && (this.Ontology_type.Equals(other.Ontology_type))
                          && (this.Ontology_namespace.Equals(other.Ontology_namespace)));
            return equal;
        }
    }

    class NetworkNodes_set_line_class
    {
        public string Set_name { get; set; }
        public float Minus_log10_pvalue { get; set; }

        public NetworkNodes_set_line_class Deep_copy()
        {
            NetworkNodes_set_line_class copy = (NetworkNodes_set_line_class)this.MemberwiseClone();
            copy.Set_name = (string)this.Set_name.Clone();
            return copy;
        }
    }

    class NetworkNode_line_class
    {
        #region Fields
        const string empty_entry = "E_m_p_t_y";

        public string Name { get; set; }
        public int NW_index { get; set; }
        public int NW_index_old { get; set; }
        public int Level { get; set; }
        public NetworkNodes_set_line_class[] Sets_that_contain_nodes { get; set; }

        public static string Empty_entry { get { return empty_entry; } }
        #endregion

        #region Constructor
        public NetworkNode_line_class()
        {
            Name = (string)Empty_entry.Clone();
            NW_index_old = -1;
            NW_index = -1;
            Level = -1;
            Sets_that_contain_nodes = new NetworkNodes_set_line_class[0];
        }
        #endregion

        public bool Equal_in_standard_way(NetworkNode_line_class other)
        {
            bool equal = (this.Name.Equals(other.Name));
            return equal;
        }

        public static NetworkNode_line_class[] Order_in_standard_way(NetworkNode_line_class[] nodes)
        {
            return nodes.OrderBy(l => l.Name).ToArray();
        }

        #region Copy
        public NetworkNode_line_class Deep_copy()
        {
            NetworkNode_line_class node_line = (NetworkNode_line_class)this.MemberwiseClone();
            node_line.Name = (string)this.Name.Clone();
            int sets_length = this.Sets_that_contain_nodes.Length;
            node_line.Sets_that_contain_nodes = new NetworkNodes_set_line_class[sets_length];
            for (int indexSet = 0; indexSet < sets_length; indexSet++)
            {
                node_line.Sets_that_contain_nodes[indexSet] = this.Sets_that_contain_nodes[indexSet].Deep_copy();
            }
            return node_line;
        }
        #endregion
    }

    class NetworkNode_class
    {
        #region Fields
        public NetworkNode_line_class[] Nodes { get; set; }
        public bool Index_change_adopted { get; set; }
        public int Nodes_length { get { return Nodes.Length; } }
        public Ontology_direction_enum Direction { get; set; }
        #endregion

        public NetworkNode_class()
        {
            Nodes = new NetworkNode_line_class[0];
            Index_change_adopted = true;
        }

        #region Check
        public bool Correctness_check()
        {
            bool ok = true;
            if (!Index_change_adopted)
            {
                string text = "Index changes are not adopted";
                throw new Exception(text);
            }
            return ok;
        }
        #endregion

        #region Order
        public void Order_by_name()
        {
            Nodes = Nodes.OrderBy(l => l.Name).ToArray();
        }

        public void Order_by_nw_index_old()
        {
            Nodes = Nodes.OrderBy(l => l.NW_index_old).ToArray();
        }

        public void Order_by_nw_index()
        {
            Nodes = Nodes.OrderBy(l => l.NW_index).ToArray();
        }

        public void Order_by_level()
        {
            Nodes = Nodes.OrderBy(l => l.Level).ToArray();
        }
        #endregion

        #region Get
        public NetworkNode_line_class Get_indexed_node_line_if_index_is_correct(int indexNode)
        {
            NetworkNode_line_class node_line = Nodes[indexNode];
            if (node_line.NW_index != indexNode)
            {
                node_line = new NetworkNode_line_class();
                throw new Exception();
            }
            return node_line;
        }

        public NetworkNode_line_class Get_indexed_old_node_line_if_index_old_is_correct(int indexNode_old)
        {
            NetworkNode_line_class node_line = Nodes[indexNode_old];
            if (node_line.NW_index_old != indexNode_old)
            {
                node_line = new NetworkNode_line_class();
                string text = typeof(NetworkNode_class).Name + ": Get indexed old node line, Indexes_old do not match (" + indexNode_old + " <-> " + node_line.NW_index_old + ")";
                throw new Exception(text);
            }
            return node_line;
        }

        public int Get_max_nw_index()
        {
            int max_index = -1;
            foreach (NetworkNode_line_class node_line in Nodes)
            {
                max_index = Math.Max(max_index, node_line.NW_index);
            }
            return max_index;
        }

        public string[] Get_all_nodeNames_of_indicated_levels(params int[] levels)
        {
            levels = levels.Distinct().OrderBy(l => l).ToArray();
            int nodes_length = Nodes_length;
            int levels_length = levels.Length;
            int indexLevel = 0;
            int levelCompare;
            List<string> nodeNames = new List<string>();
            NetworkNode_line_class node_line;
            this.Order_by_level();
            for (int indexN = 0; indexN < nodes_length; indexN++)
            {
                node_line = Nodes[indexN];
                if (node_line.Level == -1)
                {
                    throw new Exception();
                }
                levelCompare = -2;
                while ((indexLevel < levels_length) && (levelCompare < 0))
                {
                    levelCompare = levels[indexLevel].CompareTo(node_line.Level);
                    if (levelCompare < 0)
                    {
                        indexLevel++;
                    }
                    else if (levelCompare == 0)
                    {
                        nodeNames.Add(node_line.Name);
                    }
                }
            }
            return nodeNames.ToArray();
        }

        public string[] Get_all_nodeNames()
        {
            NetworkNode_line_class node_line;
            this.Order_by_level();
            int nodes_length = this.Nodes_length;
            List<string> nodeNames = new List<string>();
            for (int indexN = 0; indexN < nodes_length; indexN++)
            {
                node_line = Nodes[indexN];
                nodeNames.Add(node_line.Name);
            }
            return nodeNames.ToArray();

        }
        #endregion

        #region Get dictionaries
        public Dictionary<string, int> Get_name_index_dictionary()
        {
            Dictionary<string, int> name_index_dict = new Dictionary<string, int>();
            if (Correctness_check())
            {
                foreach (NetworkNode_line_class node_line in Nodes)
                {
                    name_index_dict.Add(node_line.Name, node_line.NW_index);
                }
            }
            return name_index_dict;
        }
        #endregion

        #region Set scp levels
        public void Set_processLevel_for_all_nodes_based_on_dictionary(Dictionary<string,int> scp_level_dict)
        {
            int nodes_length = this.Nodes.Length;
            NetworkNode_line_class node_line;
            for (int indexNodes = 0; indexNodes < nodes_length; indexNodes++)
            {
                node_line = this.Nodes[indexNodes];
                if (scp_level_dict.ContainsKey(node_line.Name))
                {
                    node_line.Level = scp_level_dict[node_line.Name];
                }
            }
        }

        public void Set_processLevel_for_all_nodes_to_level3()
        {
            int nodes_length = this.Nodes.Length;
            NetworkNode_line_class node_line;
            for (int indexNodes = 0; indexNodes < nodes_length; indexNodes++)
            {
                node_line = this.Nodes[indexNodes];
                node_line.Level = 3;
            }
        }

        public void Set_level_for_all_nodes(int level)
        {
            int nodes_length = this.Nodes_length;
            NetworkNode_line_class node_line;
            for (int indexNode = 0; indexNode < nodes_length; indexNode++)
            {
                node_line = this.Nodes[indexNode];
                node_line.Level = level;
            }
        }
        #endregion

        #region Generate node colors
        public yed_node_color_line_class[] Get_yED_node_colors_based_on_sets_if_not_indicated_different_in_dictionary(Dictionary<string, Color[]> nodeLabel_colors_dict)
        {
            NetworkNode_line_class node_line;
            int nodes_length = this.Nodes_length;

            yed_node_color_line_class[] yED_node_color_lines = new yed_node_color_line_class[nodes_length];
            yed_node_color_line_class new_yED_node_color_line;
            int current_colors_length;
            Color_specification_line_class color_speci_line;
            for (int indexN = 0; indexN < nodes_length; indexN++)
            {
                node_line = this.Nodes[indexN];
                new_yED_node_color_line = new yed_node_color_line_class();
                new_yED_node_color_line.NodeName = (string)node_line.Name.Clone();
                if (!nodeLabel_colors_dict.ContainsKey(new_yED_node_color_line.NodeName))
                {
                    throw new Exception();
                }
                current_colors_length = nodeLabel_colors_dict[new_yED_node_color_line.NodeName].Length;
                new_yED_node_color_line.Color_specifications = new Color_specification_line_class[current_colors_length];
                for (int indexC = 0; indexC < current_colors_length; indexC++)
                {
                    color_speci_line = new Color_specification_line_class();
                    color_speci_line.Fill_color = nodeLabel_colors_dict[new_yED_node_color_line.NodeName][indexC];
                    new_yED_node_color_line.Color_specifications[indexC] = color_speci_line;
                }
                yED_node_color_lines[indexN] = new_yED_node_color_line;
            }
            return yED_node_color_lines;
        }
        #endregion

        #region Keep node lines
        public void Keep_only_input_nodeNames_and_reindex(string[] input_node_Names)
        {
            input_node_Names = input_node_Names.Distinct().OrderBy(l => l).ToArray();
            if (input_node_Names.Length == 0)
            {
                throw new Exception("no nodes");
            }
            string input_nodeID;
            int input_nodes_length = input_node_Names.Length;
            int this_length = Nodes_length;
            NetworkNode_line_class node_line;
            int stringCompare;
            int indexInput = 0;
            List<NetworkNode_line_class> keep_nodes = new List<NetworkNode_line_class>();
            this.Order_by_name();
            for (int indexThis = 0; indexThis < this_length; indexThis++)
            {
                node_line = Nodes[indexThis];
                stringCompare = -2;
                while ((indexInput < input_nodes_length) && (stringCompare < 0))
                {
                    input_nodeID = input_node_Names[indexInput];
                    stringCompare = input_nodeID.CompareTo(node_line.Name);
                    if (stringCompare < 0)
                    {
                        indexInput++;
                    }
                    else if (stringCompare == 0)
                    {
                        keep_nodes.Add(node_line);
                    }
                }
            }
            Nodes = keep_nodes.ToArray();
            Reindex_nodes_and_set_index_old();
        }
        #endregion

        #region Direction
        public void Reverse_direction()
        {
            switch (Direction)
            {
                case Ontology_direction_enum.Child_parent:
                    Direction = Ontology_direction_enum.Parent_child;
                    break;
                case Ontology_direction_enum.Parent_child:
                    Direction = Ontology_direction_enum.Child_parent;
                    break;
                default:
                    break;
            }
        }
        #endregion

        #region Generate
        private NetworkNode_line_class[] Get_new_nodes_with_indexes_above_max_index(string[] add_nodeNames)
        {
            int add_index = Get_max_nw_index();
            add_nodeNames = add_nodeNames.Distinct().OrderBy(l => l).ToArray();
            string add_nodeName;
            int add_nodeNames_length = add_nodeNames.Length;
            List<NetworkNode_line_class> new_nodes = new List<NetworkNode_line_class>();
            NetworkNode_line_class new_node_line;
            for (int indexN = 0; indexN < add_nodeNames_length; indexN++)
            {
                add_nodeName = add_nodeNames[indexN];
                add_index++;
                new_node_line = new NetworkNode_line_class();
                new_node_line.Name = (string)add_nodeName.Clone();
                new_node_line.NW_index = add_index;
                new_nodes.Add(new_node_line);
            }
            return new_nodes.ToArray();
        }

        private void Add_nodes(NetworkNode_line_class[] add_nodes)
        {
            int add_length = add_nodes.Length;
            int this_length = this.Nodes_length;
            int new_length = add_length + this_length;
            NetworkNode_line_class[] new_nodes = new NetworkNode_line_class[new_length];
            int indexNew = -1;
            for (int indexThis = 0; indexThis < this_length; indexThis++)
            {
                indexNew++;
                new_nodes[indexNew] = this.Nodes[indexThis];
            }
            for (int indexAdd = 0; indexAdd < add_length; indexAdd++)
            {
                indexNew++;
                new_nodes[indexNew] = add_nodes[indexAdd].Deep_copy();
            }
            Nodes = new_nodes;
        }

        private void Remove_duplicates()
        {
            Nodes = NetworkNode_line_class.Order_in_standard_way(Nodes);
            int nodes_length = Nodes_length;
            NetworkNode_line_class kept_node_line;
            List<NetworkNode_line_class> kept_nodes_list = new List<NetworkNode_line_class>();
            for (int indexN = 0; indexN < nodes_length; indexN++)
            {
                kept_node_line = this.Nodes[indexN];
                if ((indexN == 0) || (!kept_node_line.Equal_in_standard_way(this.Nodes[indexN - 1])))
                {
                    kept_nodes_list.Add(kept_node_line);
                }
                else
                {
                    //throw new Exception();
                }
            }
            Nodes = kept_nodes_list.ToArray();
        }

        public void Reindex_nodes_and_set_index_old()
        {
            Nodes = NetworkNode_line_class.Order_in_standard_way(Nodes);
            NetworkNode_line_class node_line;
            int nodes_length = Nodes_length;
            for (int indexN = 0; indexN < nodes_length; indexN++)
            {
                node_line = Nodes[indexN];
                node_line.NW_index_old = node_line.NW_index;
                node_line.NW_index = indexN;
            }
            Index_change_adopted = false;
        }

        public NetworkNode_class Merge_this_nodes_with_other_nodes_and_get_new_indexes_of_other_UN(NetworkNode_class inputOther)
        {
            NetworkNode_class other = inputOther.Deep_copy();
            this.Order_by_name();
            other.Order_by_name();
            int this_length = this.Nodes.Length;
            int indexThis = 0;
            int other_length = other.Nodes.Length;
            int indexOther = 0;
            int stringCompare = 0;
            List<NetworkNode_line_class> mergedUN = new List<NetworkNode_line_class>();
            NetworkNode_line_class this_line = new NetworkNode_line_class();
            NetworkNode_line_class other_line = new NetworkNode_line_class();
            int overlapping_nodes_count = 0;
            int indexOtherNew = this_length - 1;
            List<NetworkNodes_set_line_class> new_set_list = new List<NetworkNodes_set_line_class>();
            while ((indexThis < this_length) || (indexOther < other_length))
            {
                if ((indexThis < this_length) && (indexOther < other_length))
                {
                    this_line = this.Nodes[indexThis];
                    other_line = other.Nodes[indexOther];
                    this_line.NW_index_old = this_line.NW_index;
                    stringCompare = other_line.Name.CompareTo(this_line.Name);
                }
                else if (indexThis < this_length)
                {
                    this_line = this.Nodes[indexThis];
                    this_line.NW_index_old = this_line.NW_index;
                    stringCompare = 2;
                }
                else // if (indexOther < other_count)
                {
                    other_line = other.Nodes[indexOther];
                    stringCompare = -2;
                }
                if (stringCompare < 0) //other_line is not in this.un
                {
                    indexOtherNew++;
                    other_line.NW_index_old = other_line.NW_index;
                    other_line.NW_index = indexOtherNew;
                    NetworkNode_line_class newOtherLine = other_line.Deep_copy();
                    newOtherLine.NW_index = indexOtherNew;
                    newOtherLine.NW_index_old = indexOtherNew;
                    mergedUN.Add(newOtherLine);
                    indexOther++;
                }
                else if (stringCompare == 0)
                {
                    new_set_list.Clear();
                    new_set_list.AddRange(this_line.Sets_that_contain_nodes);
                    new_set_list.AddRange(other_line.Sets_that_contain_nodes);
                    this_line.Sets_that_contain_nodes = new_set_list.ToArray();
                    mergedUN.Add(this_line);
                    other_line.NW_index_old = other_line.NW_index;
                    other_line.NW_index = this_line.NW_index;
                    overlapping_nodes_count++;
                    indexOther++;
                    indexThis++;
                }
                else // (stringCompare > 0)
                {
                    mergedUN.Add(this_line);
                    indexThis++;
                }
            }
            Nodes = mergedUN.ToArray();
            other.Index_change_adopted = false;

            return other;
        }

        private void Add_new_obo_node_lines_remvoe_duplicates_and_reindex(NetworkNode_line_class[] add_node_lines)
        {
            Add_nodes(add_node_lines);
            Remove_duplicates();
            Reindex_nodes_and_set_index_old();
        }

        public void Add_new_nodes_remove_duplicates_and_reindex(string[] newNodes_names)
        {
            if (Correctness_check())
            {
                NetworkNode_line_class[] add_node_lines = Get_new_nodes_with_indexes_above_max_index(newNodes_names);
                Add_new_obo_node_lines_remvoe_duplicates_and_reindex(add_node_lines);
            }
        }
        #endregion

        #region Read write copy
        public NetworkNode_class Deep_copy()
        {
            NetworkNode_class copy = (NetworkNode_class)this.MemberwiseClone();
            int nodes_length = Nodes_length;
            copy.Nodes = new NetworkNode_line_class[nodes_length];
            for (int indexN = 0; indexN < nodes_length; indexN++)
            {
                copy.Nodes[indexN] = this.Nodes[indexN].Deep_copy();
            }
            return copy;
        }
        #endregion
    }
}
