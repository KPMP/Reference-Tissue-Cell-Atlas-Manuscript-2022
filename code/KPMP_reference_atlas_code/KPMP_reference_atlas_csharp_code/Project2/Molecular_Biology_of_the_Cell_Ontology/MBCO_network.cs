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
    class Network_target_line_class
    {
        #region Fields
        public int NW_index { get; set; }
        public Obo_interaction_type_enum Interaction_type { get; set; }
        public NWedge_type_enum Edge_type { get; set; }
        public float Width { get; set; }
        public string Label { get; set; }
        #endregion

        public Network_target_line_class()
        {
            Width = 1;
            Label = "";
            Edge_type = NWedge_type_enum.Arrow;
        }

        public Network_target_line_class(int nw_index, Obo_interaction_type_enum interaction_type) : this()
        {
            this.NW_index = nw_index;
            this.Interaction_type = interaction_type;
        }

        #region Equals, compare, order in standard way
        public bool Equals_in_standard_way(Network_target_line_class other)
        {
            bool equal = ((this.NW_index.Equals(other.NW_index))
                          && (this.Interaction_type.Equals(other.Interaction_type)));
            return equal;
        }

        public int Compare_in_standard_way(Network_target_line_class other)
        {
            int targetCompare = this.NW_index - other.NW_index;
            if (targetCompare == 0)
            {
                targetCompare = this.Interaction_type.CompareTo(other.Interaction_type);
            }
            return targetCompare;
        }

        public static Network_target_line_class[] Order_in_standard_way(Network_target_line_class[] target_lines)
        {
            return target_lines.OrderBy(l => l.NW_index).ThenBy(l => l.Interaction_type).ToArray();
        }
        #endregion

        #region Copy
        public Network_target_line_class Deep_copy()
        {
            Network_target_line_class copy = (Network_target_line_class)this.MemberwiseClone();
            copy.Label = (string)this.Label.Clone();
            return copy;
        }
        #endregion
    }

    class Network_line_class
    {
        public Network_target_line_class[] Targets { get; set; }
        public int Targets_length { get { return Targets.Length; } }

        public Network_line_class()
        {
            Targets = new Network_target_line_class[0];
        }

        public void Add_not_existing_targets_and_order_in_standard_way(params Network_target_line_class[] new_target_lines)
        {
            new_target_lines = Network_target_line_class.Order_in_standard_way(new_target_lines);
            Network_target_line_class new_target_line;
            Network_target_line_class this_target_line;
            int new_target_length = new_target_lines.Length;
            int this_length = Targets_length;
            int indexThis = 0;
            int targetCompare;

            List<Network_target_line_class> new_targets = new List<Network_target_line_class>();

            #region Get new targets that do not exist among old targets
            for (int indexNew = 0; indexNew < new_target_length; indexNew++)
            {
                new_target_line = new_target_lines[indexNew];
                if ((indexNew == 0)
                    || (!new_target_line.Equals_in_standard_way(new_target_lines[indexNew - 1])))
                {
                    targetCompare = -2;
                    while ((indexThis < this_length) && (targetCompare < 0))
                    {
                        this_target_line = this.Targets[indexThis];
                        targetCompare = this_target_line.Compare_in_standard_way(new_target_line);
                        if (targetCompare < 0)
                        {
                            indexThis++;
                        }
                    }
                    if (targetCompare != 0)
                    {
                        new_targets.Add(new_target_line);
                    }
                }
            }
            #endregion

            #region Add new targets and order
            new_targets.AddRange(this.Targets);
            this.Targets = Network_target_line_class.Order_in_standard_way(new_targets.ToArray());
            #endregion
        }

        public Network_line_class Deep_copy()
        {
            Network_line_class copy = (Network_line_class)this.MemberwiseClone();
            int targets_length = Targets.Length;
            copy.Targets = new Network_target_line_class[targets_length];
            for (int indexT = 0; indexT < targets_length; indexT++)
            {
                copy.Targets[indexT] = (Network_target_line_class)this.Targets[indexT].Deep_copy();
            }
            return copy;
        }
    }

    class Network_class
    {
        #region Fields
        public Network_line_class[] NW { get; set; }
        public NetworkNode_class Nodes { get; set; }
        public int NW_length { get { return NW.Length; } }
        #endregion

        public Network_class()
        {
            NW = new Network_line_class[0];
            Nodes = new NetworkNode_class();
        }

        public void Switch_all_edge_types_to_input_type(NWedge_type_enum edge)
        {
            foreach (Network_line_class nw_line in NW)
            {
                foreach (Network_target_line_class target in nw_line.Targets)
                {
                    target.Edge_type = edge;
                }
            }
        }

        protected void Reindex_network_based_on_nodes_and_add_new_nw_lines_if_neccessary()
        {
            int nw_length = NW_length;

            #region Get indexOld - indexNew array
            int[] indexOld_indexNew = new int[nw_length];
            for (int indexArray = 0; indexArray < nw_length; indexArray++)
            {
                indexOld_indexNew[indexArray] = -1;
            }

            Nodes.Order_by_nw_index_old();
            int nodes_length = Nodes.Nodes_length;
            NetworkNode_line_class node_line;
            for (int indexN = 0; indexN < nodes_length; indexN++)
            {
                node_line = Nodes.Nodes[indexN];
                if (node_line.NW_index_old < nw_length)
                {
                    if (indexOld_indexNew[node_line.NW_index_old] != -1)
                    {
                        throw new Exception("index new has already been assigned");
                    }
                    else
                    {
                        indexOld_indexNew[node_line.NW_index_old] = node_line.NW_index;
                    }
                }
            }
            #endregion

            int new_nw_length = Nodes.Get_max_nw_index() + 1;
            Network_line_class[] new_nw = new Network_line_class[new_nw_length];
            Network_line_class source_nw_line;
            Network_target_line_class old_target_line;
            int targets_length;
            Nodes.Order_by_nw_index();
            NetworkNode_line_class source_node_line;
            NetworkNode_line_class target_node_line;
            int newIndex;
            List<Network_target_line_class> kept_target_lines = new List<Network_target_line_class>();
            for (int indexNodeNew = 0; indexNodeNew < new_nw_length; indexNodeNew++)
            {
                source_node_line = Nodes.Get_indexed_node_line_if_index_is_correct(indexNodeNew);
                if (source_node_line.NW_index_old < nw_length)
                {
                    source_nw_line = NW[source_node_line.NW_index_old];
                    targets_length = source_nw_line.Targets.Length;
                    kept_target_lines.Clear();
                    for (int indexT = 0; indexT < targets_length; indexT++)
                    {
                        old_target_line = source_nw_line.Targets[indexT];
                        newIndex = indexOld_indexNew[old_target_line.NW_index];
                        if (newIndex != -1)
                        {
                            target_node_line = Nodes.Get_indexed_node_line_if_index_is_correct(newIndex);
                            if (target_node_line.NW_index_old != old_target_line.NW_index) { throw new Exception("index old does not match"); }
                            old_target_line.NW_index = target_node_line.NW_index;
                            kept_target_lines.Add(old_target_line);
                        }
                    }
                    source_nw_line.Targets = kept_target_lines.OrderBy(l => l.NW_index).ToArray();
                    new_nw[source_node_line.NW_index] = source_nw_line;
                }
                else
                {
                    new_nw[source_node_line.NW_index] = new Network_line_class();
                }
            }
            NW = new_nw;
            Nodes.Index_change_adopted = true;
        }

        private Dictionary<string, Visualization_nw_edge_characterisation_line_class[]> Get_visualization_sourceName_nw_edge_characterisation_dictionary()
        {
            Dictionary<string, Visualization_nw_edge_characterisation_line_class[]> source_name_target_name = new Dictionary<string, Visualization_nw_edge_characterisation_line_class[]>();
            Visualization_nw_edge_characterisation_line_class new_edge_line;
            int nw_length = NW_length;
            Network_line_class nw_line;
            Network_target_line_class target_line;
            NetworkNode_line_class source_node_line;
            NetworkNode_line_class target_node_line;
            int targets_length;
            Nodes.Order_by_nw_index();
            List<Visualization_nw_edge_characterisation_line_class> nw_edge_characterisation_list = new List<Visualization_nw_edge_characterisation_line_class>();
            Visualization_nw_node_characterisation_line_class new_nw_node_characterisation_line;
            string source_name;
            for (int indexNW = 0; indexNW < nw_length; indexNW++)
            {
                nw_line = NW[indexNW];
                source_node_line = Nodes.Get_indexed_node_line_if_index_is_correct(indexNW);
                nw_edge_characterisation_list.Clear();
                targets_length = nw_line.Targets_length;
                source_name = source_node_line.Name;
                for (int indexT = 0; indexT < targets_length; indexT++)
                {
                    target_line = nw_line.Targets[indexT];
                    target_node_line = Nodes.Get_indexed_node_line_if_index_is_correct(target_line.NW_index);
                    new_edge_line = new Visualization_nw_edge_characterisation_line_class();
                    new_edge_line.Edge_width = target_line.Width;
                    new_edge_line.Edge_label = (string)target_line.Label.Clone();
                    new_edge_line.Target = (string)target_node_line.Name.Clone();
                    switch (target_line.Edge_type)
                    {
                        case NWedge_type_enum.Arrow:
                            new_edge_line.EdgeArrow_type = EdgeArrow_type_enum.Arrow;
                            break;
                        case NWedge_type_enum.Dotted_line:
                            new_edge_line.EdgeArrow_type = EdgeArrow_type_enum.Dotted_line;
                            break;
                        case NWedge_type_enum.Thick_dotted_line:
                            new_edge_line.EdgeArrow_type = EdgeArrow_type_enum.Thick_dotted_line;
                            break;
                        case NWedge_type_enum.Dashed_line:
                            new_edge_line.EdgeArrow_type = EdgeArrow_type_enum.Dashed_line;
                            break;
                        default:
                            throw new Exception();
                    }
                    nw_edge_characterisation_list.Add(new_edge_line);
                }
                new_nw_node_characterisation_line = new Visualization_nw_node_characterisation_line_class();
                new_nw_node_characterisation_line.NodeName = (string)source_node_line.Name.Clone();
                new_nw_node_characterisation_line.Level = source_node_line.Level;
                source_name_target_name.Add(source_name, nw_edge_characterisation_list.ToArray());
            }
            return source_name_target_name;
        }

        #region Generate
        private void Reindex_and_add_nw_connections_form_networkTable_lines(NetworkTable_line_class[] networkTable_lines)
        {
            Reindex_network_based_on_nodes_and_add_new_nw_lines_if_neccessary();
            int sigNW_length = networkTable_lines.Length;
            NetworkTable_line_class networkTable_line;
            Dictionary<string, int> name_nw_index_dict = Nodes.Get_name_index_dictionary();
            string source_name;
            string target_name;
            int source_index;
            int target_index;
            Network_target_line_class new_target_line;
            for (int indexSigNW = 0; indexSigNW < sigNW_length; indexSigNW++)
            {
                networkTable_line = networkTable_lines[indexSigNW];
                target_name = networkTable_line.Target;
                source_name = networkTable_line.Source;
                if ((!source_name.Equals(Global_class.Empty_entry)) && (!target_name.Equals(Global_class.Empty_entry)))
                {
                    source_index = name_nw_index_dict[source_name];
                    target_index = name_nw_index_dict[target_name];
                    new_target_line = new Network_target_line_class(target_index, Obo_interaction_type_enum.E_m_p_t_y);
                    new_target_line.Width = networkTable_line.Width;
                    new_target_line.Label = (string)networkTable_line.Edge_label.Clone();
                    new_target_line.Edge_type = networkTable_line.Edge_type;
                    NW[source_index].Add_not_existing_targets_and_order_in_standard_way(new_target_line);
                }
            }
        }

        public void Add_from_networkTable_lines(NetworkTable_line_class[] networkTable_lines)
        {
            #region Get all nodes
            List<string> all_nodes = new List<string>();
            foreach (NetworkTable_line_class networkTable_line in networkTable_lines)
            {
                all_nodes.Add(networkTable_line.Source);
                all_nodes.Add(networkTable_line.Target);
            }
            all_nodes = all_nodes.Distinct().ToList();
            #endregion

            Nodes.Add_new_nodes_remove_duplicates_and_reindex(all_nodes.ToArray());
            Reindex_and_add_nw_connections_form_networkTable_lines(networkTable_lines);
        }

        public void Add_single_nodes(params string[] nodeNames)
        {
            Nodes.Add_new_nodes_remove_duplicates_and_reindex(nodeNames.ToArray());
            Reindex_network_based_on_nodes_and_add_new_nw_lines_if_neccessary();
        }
        #endregion

        #region Keep
        public void Keep_only_input_nodeNames(string[] input_node_names)
        {
            Nodes.Keep_only_input_nodeNames_and_reindex(input_node_names);
            Reindex_network_based_on_nodes_and_add_new_nw_lines_if_neccessary();
        }
        #endregion

        #region Direction
        public void Reverse_direction()
        {
            int nw_length = NW.Length;
            Network_line_class nw_line;
            Network_line_class[] reversed_nw = new Network_line_class[nw_length];
            for (int indexRevNW = 0; indexRevNW < nw_length; indexRevNW++)
            {
                reversed_nw[indexRevNW] = new Network_line_class();
            }
            Network_target_line_class target_line;
            Network_target_line_class reversed_target_line;
            int target_length;
            for (int indexNW = 0; indexNW < nw_length; indexNW++)
            {
                nw_line = NW[indexNW];
                target_length = nw_line.Targets_length;
                for (int indexT = 0; indexT < target_length; indexT++)
                {
                    target_line = nw_line.Targets[indexT];
                    reversed_target_line = new Network_target_line_class();
                    reversed_target_line.NW_index = indexNW;
                    reversed_target_line.Interaction_type = target_line.Interaction_type;
                    reversed_target_line.Edge_type = target_line.Edge_type;
                    reversed_nw[target_line.NW_index].Add_not_existing_targets_and_order_in_standard_way(reversed_target_line);
                }
            }
            Nodes.Reverse_direction();
            NW = reversed_nw;
        }

        public void Transform_into_undirected_double_network()
        {
            Network_class copy = this.Deep_copy();
            copy.Reverse_direction();
            Network_line_class this_line;
            Network_line_class copy_line;
            int nw_length = this.NW_length;
            for (int indexNW = 0; indexNW < nw_length; indexNW++)
            {
                this_line = this.NW[indexNW];
                copy_line = copy.NW[indexNW];
                this_line.Add_not_existing_targets_and_order_in_standard_way(copy_line.Targets);
            }
        }

        public void Transform_into_undirected_single_network_and_set_all_widths_to_one()
        {
            Transform_into_undirected_double_network();
            int nw_length = this.NW_length;
            Network_target_line_class current_target;
            int targets_length;

            List<Network_target_line_class> current_kept_network_target_lines = new List<Network_target_line_class>();
            Network_line_class this_line;
            for (int indexNW = 0; indexNW < nw_length; indexNW++)
            {
                this_line = this.NW[indexNW];
                targets_length = this_line.Targets_length;
                current_kept_network_target_lines.Clear();
                for (int indexT = 0; indexT < targets_length; indexT++)
                {
                    current_target = this_line.Targets[indexT];
                    if (current_target.NW_index <= indexNW)
                    {
                        current_target.Width = 1;
                        current_kept_network_target_lines.Add(current_target);
                    }
                }
                this_line.Targets = current_kept_network_target_lines.ToArray();
            }
        }
        #endregion

        public void Merge_this_network_with_other_network(Network_class input_otherNW)
        {
            int this_nodes_length_old = this.Nodes.Nodes.Length;
            Network_class otherNW = input_otherNW.Deep_copy();
            otherNW.Nodes = this.Nodes.Merge_this_nodes_with_other_nodes_and_get_new_indexes_of_other_UN(otherNW.Nodes);

            List<Network_target_line_class> target_nodes_list = new List<Network_target_line_class>();
            NetworkNode_line_class source_node_line;
            NetworkNode_line_class target_node_line;
            Network_target_line_class target_line;
            int targets_length;

            otherNW.Nodes.Order_by_nw_index_old();
            Dictionary<int, Network_target_line_class[]> other_nw_nodes_dict = new Dictionary<int, Network_target_line_class[]>();
            int other_nw_length = otherNW.NW_length;
            Network_line_class other_nw_line;
            for (int indexOtherOld = 0; indexOtherOld < other_nw_length; indexOtherOld++)
            {
                other_nw_line = otherNW.NW[indexOtherOld];
                source_node_line = otherNW.Nodes.Get_indexed_old_node_line_if_index_old_is_correct(indexOtherOld);
                targets_length = other_nw_line.Targets_length;
                target_nodes_list.Clear();
                for (int indexT = 0; indexT < targets_length; indexT++)
                {
                    target_line = other_nw_line.Targets[indexT];
                    target_node_line = otherNW.Nodes.Get_indexed_old_node_line_if_index_old_is_correct(target_line.NW_index);
                    target_line.NW_index = target_node_line.NW_index;
                    target_nodes_list.Add(target_line);
                }
                other_nw_nodes_dict.Add(source_node_line.NW_index, target_nodes_list.ToArray());
            }

            this.Nodes.Order_by_nw_index_old();
            Dictionary<int, Network_target_line_class[]> this_nw_nodes_dict = new Dictionary<int, Network_target_line_class[]>();
            int this_nw_length = this.NW_length;
            Network_line_class this_nw_line;
            for (int indexThisOld = 0; indexThisOld < this_nw_length; indexThisOld++)
            {
                this_nw_line = this.NW[indexThisOld];
                source_node_line = this.Nodes.Get_indexed_old_node_line_if_index_old_is_correct(indexThisOld);
                targets_length = this_nw_line.Targets_length;
                target_nodes_list.Clear();
                for (int indexT = 0; indexT < targets_length; indexT++)
                {
                    target_line = this_nw_line.Targets[indexT];
                    target_node_line = this.Nodes.Get_indexed_old_node_line_if_index_old_is_correct(target_line.NW_index);
                    target_line.NW_index = target_node_line.NW_index;
                    target_nodes_list.Add(target_line);
                }
                this_nw_nodes_dict.Add(source_node_line.NW_index, target_nodes_list.ToArray());
            }

            this.Nodes.Order_by_nw_index();
            int new_nw_length = this.Nodes.Nodes_length;
            Network_line_class[] new_nw = new Network_line_class[new_nw_length];
            for (int indexNew = 0; indexNew < new_nw_length; indexNew++)
            {
                target_nodes_list.Clear();
                if (this_nw_nodes_dict.ContainsKey(indexNew))
                {
                    target_nodes_list.AddRange(this_nw_nodes_dict[indexNew]);
                }
                if (other_nw_nodes_dict.ContainsKey(indexNew))
                {
                    target_nodes_list.AddRange(other_nw_nodes_dict[indexNew]);
                }
                new_nw[indexNew] = new Network_line_class();
                new_nw[indexNew].Targets = target_nodes_list.ToArray();
            }
            this.NW = new_nw;
            Nodes.Index_change_adopted = true;
        }

        public string[] Get_all_scps()
        {
            return Nodes.Get_all_nodeNames();
        }

        #region Breadth first search
        public void Get_direct_neighbors_with_breadth_first_search(out int[] direct_neighbors, out int[] distances, int max_distance, params int[] seedNodeIndexes)
        {
            if (seedNodeIndexes.Length == 0)
            {
                throw new Exception("no seed nodes");
            }

            int nw_length = NW.Length;
            int[] all_nodes = new int[nw_length];
            int[] all_distances = new int[nw_length];

            //Generate and fill queue and distance arrays
            bool[] already_visited = new bool[nw_length];
            int seedNode_length = seedNodeIndexes.Length;
            for (int i = 0; i < seedNode_length; i++)
            {
                all_nodes[i] = seedNodeIndexes[i];
                all_distances[i] = 0;
                already_visited[seedNodeIndexes[i]] = true;
            }

            //Do breadth first search
            int writePointer = seedNode_length;
            int readPointer = 0;
            int actNode;
            Network_target_line_class[] target_lines;

            while ((readPointer != writePointer) && (all_distances[readPointer] + 1 <= max_distance))
            {
                actNode = all_nodes[readPointer];
                target_lines = NW[actNode].Targets;
                foreach (Network_target_line_class target_line in target_lines)
                {
                    if (!already_visited[target_line.NW_index])
                    {
                        already_visited[target_line.NW_index] = true;
                        all_nodes[writePointer] = target_line.NW_index;
                        all_distances[writePointer] = all_distances[readPointer] + 1;
                        writePointer++;
                    }
                }
                readPointer++;
            }

            int results_length = writePointer;
            direct_neighbors = new int[results_length];
            distances = new int[results_length];
            for (int indexR = 0; indexR < results_length; indexR++)
            {
                direct_neighbors[indexR] = all_nodes[indexR];
                distances[indexR] = all_distances[indexR];
            }
        }

        public void Keep_only_nodes_with_indicated_levels(params int[] levels)
        {
            string[] keep_nodeNames = this.Nodes.Get_all_nodeNames_of_indicated_levels(levels);
            Keep_only_input_nodeNames(keep_nodeNames);
        }
        #endregion

        #region Write yed network
        public void Write_yED_nw_in_results_directory_with_nodes_colored_by_set_and_sized_by_number_of_different_colors_and_sameLevel_processes_grouped(string complete_nw_name_without_extension, Shape_enum standard_node_shape, Dictionary<string, Shape_enum> nodeLabel_nodeShape_dict, Dictionary<string, Color[]> nodeLabel_colors_dict)
        {
            Dictionary<string, Visualization_nw_edge_characterisation_line_class[]> source_nw_edge_characterisation_dict = Get_visualization_sourceName_nw_edge_characterisation_dictionary();

            string complete_file_name_without_extension = complete_nw_name_without_extension;
            yED_class yed = new yED_class();
            yed.Options.Node_size_determinant = Yed_network_node_size_determinant_enum.No_of_different_colors;
            yed.Options.Group_same_level_processes = true;
            yed.Options.NodeLabel_nodeShape_dict = nodeLabel_nodeShape_dict;
            yed.Options.Node_shape = standard_node_shape;
            yed_node_color_line_class[] node_colors = this.Nodes.Get_yED_node_colors_based_on_sets_if_not_indicated_different_in_dictionary(nodeLabel_colors_dict);
            yed.Write_yED_file(this.Nodes.Nodes, source_nw_edge_characterisation_dict, complete_file_name_without_extension, node_colors, new yed_node_color_line_class[0]);
        }
        #endregion

        #region copy
        public Network_class Deep_copy()
        {
            Network_class copy = (Network_class)this.MemberwiseClone();
            int nw_length = NW_length;
            copy.NW = new Network_line_class[nw_length];
            for (int indexNW = 0; indexNW < nw_length; indexNW++)
            {
                copy.NW[indexNW] = this.NW[indexNW].Deep_copy();
            }
            copy.Nodes = this.Nodes.Deep_copy();
            return copy;
        }
        #endregion
    }
}
