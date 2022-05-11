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


namespace MBCO_network_and_nodes
{

    class Leave_out_line_class
    {
        public int ProcessLevel { get; set; }
        public string ProcessName { get; set; }
        public string Left_out_processName { get; set; }
        public float New_symbols_rank_increase_divided_by_original_process_size { get; set; }
        public float Edge_weight { get; set; }

        public Leave_out_line_class Deep_copy()
        {
            Leave_out_line_class copy = (Leave_out_line_class)this.MemberwiseClone();
            copy.ProcessName = (string)this.ProcessName.Clone();
            copy.Left_out_processName = (string)this.Left_out_processName.Clone();
            return copy;
        }

    }

    class Leave_out_input_readWriteOptions_class : ReadWriteOptions_base
    {
        public Leave_out_input_readWriteOptions_class(Ontology_type_enum ontology)
        {
            switch (ontology)
            {
                case Ontology_type_enum.Mbco_level3:
                    File = Global_directory_class.Complete_mbco_inferred_scp_relationships_v11_fileName;
                    break;
                default:
                    throw new Exception();
            }
            Key_propertyNames = new string[] { "ProcessLevel", "Left_out_processName", "ProcessName", "New_symbols_rank_increase_divided_by_original_process_size" };
            Key_columnNames = Key_propertyNames;
            File_has_headline = true;
            HeadlineDelimiters = new char[] { Global_class.Tab };
            LineDelimiters = new char[] { Global_class.Tab };
            Report = ReadWrite_report_enum.Report_main;
        }
    }

    class Leave_out_readWriteOptions_class : ReadWriteOptions_base
    {
        public Leave_out_readWriteOptions_class(string fileName)
        {
            string directory = Global_directory_class.Results_directory;
            ReadWriteClass.Create_directory_if_it_does_not_exist(directory);
            File = directory + fileName;
            Key_propertyNames = new string[] { "ProcessLevel", "Left_out_processName", "ProcessID", "ProcessName", "Original_process_size", "Process_size_after_left_out", "New_symbols_count", "New_symbols_count_divided_by_original_process_size", "New_symbols_rank_increase", "New_symbols_rank_increase_divided_by_original_process_size", "Edge_weight" };
            Key_columnNames = Key_propertyNames;
            File_has_headline = true;
            HeadlineDelimiters = new char[] { Global_class.Tab };
            LineDelimiters = new char[] { Global_class.Tab };
            Report = ReadWrite_report_enum.Report_main;
        }
    }

    class LeaveOut_options_class
    {
        public Ontology_type_enum Ontology { get; private set; }

        public LeaveOut_options_class(Ontology_type_enum ontology)
        {
            this.Ontology = ontology;
        }
    }

    class Leave_out_class
    {
        #region Fields
        public Leave_out_line_class[] Leave_out_lines { get; set; }
        public LeaveOut_options_class Options { get; set; }
        #endregion

        public Leave_out_class(Ontology_type_enum ontology)
        {
            Options = new LeaveOut_options_class(ontology);
        }

        #region Order
        public void Order_by_processLevel()
        {
            this.Leave_out_lines = this.Leave_out_lines.OrderBy(l => l.ProcessLevel).ToArray();
        }

        public void Order_by_processLevel_descending_newSymbolsRankIncreaseDividedByOriginalProcessSize()
        {
            this.Leave_out_lines = this.Leave_out_lines.OrderBy(l => l.ProcessLevel).ThenByDescending(l => l.New_symbols_rank_increase_divided_by_original_process_size).ToArray();
        }

        public void Order_by_left_out_processName_processName()
        {
            this.Leave_out_lines = this.Leave_out_lines.OrderBy(l => l.Left_out_processName).ThenBy(l => l.ProcessName).ToArray();
        }

        public void Order_by_processName()
        {
            this.Leave_out_lines = this.Leave_out_lines.OrderBy(l => l.ProcessName).ToArray();
        }
        #endregion

        #region Check
        private void Check_for_duplicated_scp_scp_connections()
        {
            int leave_out_length = this.Leave_out_lines.Length;
            Leave_out_line_class leave_out_line;
            this.Leave_out_lines = this.Leave_out_lines.OrderBy(l => l.Left_out_processName).ThenBy(l => l.ProcessName).ToArray();
            for (int indexL = 0; indexL < leave_out_length; indexL++)
            {
                leave_out_line = this.Leave_out_lines[indexL];
                if ((indexL == leave_out_length - 1)
                    || (!leave_out_line.Left_out_processName.Equals(this.Leave_out_lines[indexL + 1].Left_out_processName))
                    || (!leave_out_line.ProcessName.Equals(this.Leave_out_lines[indexL + 1].ProcessName)))
                {
                }
                else
                {
                    throw new Exception();
                }
            }
        }
        #endregion

        private void Calculate_edge_weights_as_highest_fractional_new_symbols_rank_increase()
        {
            int leave_out_length = this.Leave_out_lines.Length;
            Leave_out_line_class leave_out_line;
            for (int indexLeaveOut = 0; indexLeaveOut < leave_out_length; indexLeaveOut++)
            {
                leave_out_line = this.Leave_out_lines[indexLeaveOut];
                leave_out_line.Edge_weight = leave_out_line.New_symbols_rank_increase_divided_by_original_process_size;
            }
        }

        public void Generate_by_reading_safed_file()
        {
            Read();
            Check_for_duplicated_scp_scp_connections();
            Calculate_edge_weights_as_highest_fractional_new_symbols_rank_increase();
        }

        #region Read copy
        private void Read()
        {
            Leave_out_input_readWriteOptions_class readWriteOptions = new Leave_out_input_readWriteOptions_class(this.Options.Ontology);
            this.Leave_out_lines = ReadWriteClass.ReadRawData_and_FillArray<Leave_out_line_class>(readWriteOptions);
        }

        public Leave_out_class Deep_copy()
        {
            Leave_out_class copy = (Leave_out_class)this.MemberwiseClone();
            int leave_out_length = this.Leave_out_lines.Length;
            copy.Leave_out_lines = new Leave_out_line_class[leave_out_length];
            for (int indexL = 0; indexL < leave_out_length; indexL++)
            {
                copy.Leave_out_lines[indexL] = this.Leave_out_lines[indexL].Deep_copy();
            }
            return copy;
        }
        #endregion
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class Leave_out_scp_scp_network_options_class
    {
        public Ontology_type_enum Ontology { get; set; }
        public float[] Top_quantile_of_considered_SCP_interactions_per_level { get; set; }

        public Leave_out_scp_scp_network_options_class(Ontology_type_enum ontology)
        {
            Top_quantile_of_considered_SCP_interactions_per_level = new float[] { 1F, 1F, 1F, 1F, 1F };
            Ontology = ontology;
        }

        public Leave_out_scp_scp_network_options_class Deep_copy()
        {
            Leave_out_scp_scp_network_options_class copy = (Leave_out_scp_scp_network_options_class)this.MemberwiseClone();
            copy.Top_quantile_of_considered_SCP_interactions_per_level = Array_class.Deep_copy_array(this.Top_quantile_of_considered_SCP_interactions_per_level);
            return copy;
        }
    }

    class Leave_out_scp_scp_network_class
    {
        public Leave_out_scp_scp_network_options_class Options { get; set; }
        public Network_class Scp_nw { get; set; }

        public Leave_out_scp_scp_network_class(Ontology_type_enum ontology)
        {
            this.Options = new Leave_out_scp_scp_network_options_class(ontology);
            this.Scp_nw = new Network_class();
        }

        #region Generate from leave out instance
        private NetworkTable_line_class[] Generate_networkTable_lines_for_one_level(Leave_out_line_class[] leave_out_lines, Dictionary<string, bool> scps_that_shall_not_be_connected_with_eachOther_dict)
        {
            int level = leave_out_lines[0].ProcessLevel;
            List<NetworkTable_line_class> networkTable_list = new List<NetworkTable_line_class>();
            NetworkTable_line_class networkTable_line;
            leave_out_lines = leave_out_lines.OrderByDescending(l => l.Edge_weight).ToArray();
            Leave_out_line_class leave_out_line;
            int leave_out_length = leave_out_lines.Length;
            float top_quantile_for_cutoff = this.Options.Top_quantile_of_considered_SCP_interactions_per_level[level];
            int top_quantile_lines_count = (int)Math.Round((float)leave_out_length * (top_quantile_for_cutoff));
            if (top_quantile_lines_count == leave_out_length) { top_quantile_lines_count = leave_out_length - 1; }
            for (int indexL = 0; indexL <= top_quantile_lines_count; indexL++)
            {
                leave_out_line = leave_out_lines[indexL];
                if (leave_out_line.ProcessLevel != level) { throw new Exception("Leave out contains more than one level"); }
                if ((!scps_that_shall_not_be_connected_with_eachOther_dict.ContainsKey(leave_out_line.ProcessName))
                    || (!scps_that_shall_not_be_connected_with_eachOther_dict.ContainsKey(leave_out_line.Left_out_processName)))
                {
                    networkTable_line = new NetworkTable_line_class();
                    networkTable_line.Source = (string)leave_out_line.Left_out_processName.Clone();
                    networkTable_line.Target = (string)leave_out_line.ProcessName.Clone();
                    networkTable_line.Width = leave_out_line.Edge_weight * 2;
                    networkTable_line.Edge_type = NWedge_type_enum.Dashed_line;
                    networkTable_list.Add(networkTable_line);
                    networkTable_line = new NetworkTable_line_class();
                    networkTable_line.Target = (string)leave_out_line.Left_out_processName.Clone();
                    networkTable_line.Source = (string)leave_out_line.ProcessName.Clone();
                    networkTable_line.Width = leave_out_line.Edge_weight * 2;
                    networkTable_line.Edge_type = NWedge_type_enum.Dashed_line;
                    networkTable_list.Add(networkTable_line);
                }
            }

            List<NetworkTable_line_class> final_networkTable_list = new List<NetworkTable_line_class>();
            #region Remove duplicates
            networkTable_list = networkTable_list.OrderBy(l => l.Source).ThenBy(l => l.Target).ThenByDescending(l => l.Width).ToList();
            int networkTable_count = networkTable_list.Count;
            for (int indexNT = 0; indexNT < networkTable_count; indexNT++)
            {
                networkTable_line = networkTable_list[indexNT];
                if ((indexNT == 0)
                    || (!networkTable_line.Source.Equals(networkTable_list[indexNT - 1].Source))
                    || (!networkTable_line.Target.Equals(networkTable_list[indexNT - 1].Target)))
                {
                    final_networkTable_list.Add(networkTable_line);
                }
            }
            #endregion

            return final_networkTable_list.ToArray();
        }

        private Dictionary<string, bool> Generate_dictionary_with_all_signaling_processes()
        {
            Dictionary<string, bool> signaling_processes_dict = new Dictionary<string, bool>();
            string[] signaling_processes = new string[] { "Cellular communication","Bombesin receptor signaling","Epidermal growth factor family signaling","Gastrointestinal hormone signaling",
            "Interferon signaling","Interleukin receptor signaling","Intracellular common signaling cascades of multiple pathways","Kinin–kallikrein system","Matricellular protein signaling",
            "Neuronal signaling pathways","Pattern recognition signaling","Prostanoid receptor signaling","Purinergic signaling","Signaling by extracellular matrix components","Signaling pathways involved in glucose and lipid homeostasis",
            "Signaling pathways involved in hematopoiesis","Signaling pathways regulating calcium homeostasis","Signaling pathways regulating cardiovascular homeostasis","Signaling pathways regulating steroid and sex hormone synthesis",
            "Signaling pathways regulating water homeostasis","Signaling pathways that control cell proliferation and differentiation","Steroid and sex hormone signaling","TGF-beta superfamily signaling",
            "Thyroid hormone related signaling","TNF superfamily signaling","Gastrin-releasing peptide receptor signaling","Neuromedin B receptor signaling","Epidermal growth factor receptor signaling",
            "Heparin-binding EGF-like growth factor receptor signaling","Macrophage migration inhibitory factor signaling","Neuregulin receptor signaling","Gastric inhibitory polypeptide receptor signaling",
            "Gastrin and cholecystokinin B receptor signaling","Ghrelin receptor signaling","Leptin receptor signaling","Motilin receptor signaling","Vasoactive intestinal peptide receptor signaling",
            "Interferon alpha receptor signaling","Interferon beta receptor signaling","Interferon gamma receptor signaling","Interferon receptor signaling","Interleukin 1 receptor signaling",
            "Interleukin 11 receptor signaling","Interleukin 12 receptor signaling","Interleukin 13 receptor signaling","Interleukin 15 receptor signaling","Interleukin 18 receptor signaling",
            "Interleukin 2 receptor signaling","Interleukin 21 receptor signaling","Interleukin 23 receptor signaling","Interleukin 27 receptor signaling","Interleukin 3 receptor signaling",
            "Interleukin 33 receptor signaling","Interleukin 4 receptor signaling","Interleukin 5 receptor signaling","Interleukin 6 receptor signaling","Interleukin 7 receptor signaling",
            "Interleukin 8 receptor signaling","Interleukin 9 receptor signaling","Adenylyl cyclase signaling pathway","Calcineurin-NFAT signaling pathway","CAM kinase signaling pathway",
            "G-protein coupled receptor signaling pathway","JAK-STAT signaling pathway","Mammalian target of rapamycin signaling pathway","MAPK signaling pathway","NFkB signaling pathway",
            "Nitric oxide signaling pathway","Phospholipase C signaling pathway","PI3 kinase AKT signaling pathway","Protein kinase C signaling pathway","Ras signaling pathway","Bradykinin receptor signaling",
            "Kallidin receptor signaling","CCN intercellular signaling protein family receptor signaling","Fibulin receptor signaling","Osteonectin receptor signaling","Osteopontin receptor signaling",
            "Tenascin receptor signaling","Thrombospondin receptor signaling","Brain derived neurotrophic factor receptor signaling","Cannabinoid receptor signaling","Nerve growth factor receptor signaling",
            "Roundabout signaling","Semaphorin signaling","Somatostatin receptor signaling","C-type lectin receptor signaling","MDA-5 receptor signaling","NOD-like receptor signaling",
            "RIG-I-like receptor signaling","Toll-like receptor signaling","Prostacycline receptor signaling","Prostaglandin D2 receptor signaling","Prostaglandin E2 receptor signaling",
            "Prostaglandin PGF2 alpha receptor signaling","Thromboxane receptor signaling","Purinergic P1 receptor signaling","Purinergic P2X receptor signaling","Purinergic P2Y receptor signaling",
            "Discoidin domain receptor signaling","Hyaluronan receptor CD44 signaling","Hyaluronan-mediated motility receptor signaling","Integrin receptor signaling",
            "Syndecan ectodomain shedding","Syndecan receptor signaling","Versican receptor signaling","Adiponectin receptor signaling","Glucagon receptor signaling",
            "Glucagon-like peptide-1 receptor signaling","Insulin receptor signaling","Peroxisome proliferator-activated receptor alpha signaling","Peroxisome proliferator-activated receptor gamma signaling",
            "Erythropoietin receptor signaling","Granulocyte macrophage colony-stimulating factor receptor signaling","Granulocyte-colony stimulating factor receptor signaling","Leukemia inhibitory factor receptor signaling",
            "Oncostatin-M receptor signaling","Thrombopoietin receptor signaling","Calcitonin receptor signaling","Parathyroid hormone receptor signaling","Vitamin D receptor signaling",
            "Adrenergic receptor signaling","Angiotensin receptor signaling","Endothelin receptor signaling","Muscarinic receptor signaling","Adrenocorticotropic hormone receptor signaling",
            "Corticotropin-releasing hormone receptor signaling","Follicle stimulating hormone receptor signaling","Luteinizing hormone hormone receptor signaling","Antidiuretic hormone receptor signaling",
            "Mineralocorticoid receptor signaling","Natriuretic peptide receptor signaling","Secretin receptor signaling","Ephrin receptor signaling","Fibroblast growth factor receptor signaling","Growth hormone receptor signaling",
            "Hedgehog receptor signaling","Hepatocyte growth factor receptor signaling","HIF-1 receptor signaling pathway","Hippo signaling","Insulin-like growth factor receptor signaling",
            "Notch receptor signaling","Platelet-derived growth factor receptor signaling","Prolactin receptor signaling","Retinoic acid receptor signaling","Vascular endothelial growth factor receptor signaling",
            "WNT-Beta-catenin signaling pathway","Androgen receptor signaling","Estrogen receptor signaling","Glucocorticoid receptor signaling","Progesterone receptor signaling","Activin receptor signaling",
            "Betaglycan signaling","Bone morphogenetic protein receptor signaling","Growth differentiation factor receptor signaling","Inhibin receptor signaling","Nodal growth differentiation factor receptor signaling",
            "Transforming growth factor alpha receptor signaling","Transforming growth factor beta receptor signaling","Thyroid hormone receptor signaling","Thyroid-stimulating hormone receptor signaling","Thyrotropin-releasing hormone receptor signaling",
            "Rank signaling","Tumor necrosis factor alpha receptor signaling","Tumor necrosis factor beta receptor signaling","Erk signaling pathway","SAPK-JNK signaling pathway",
            "Atrial natriuretic peptide receptor signaling","Brain natriuretic peptide receptor signaling","C-type natriuretic peptide receptor signaling" };

            signaling_processes = signaling_processes.OrderBy(l => l).ToArray();
            foreach (string signaling_process in signaling_processes)
            {
                signaling_processes_dict.Add(signaling_process, true);
            }
            return signaling_processes_dict;
        }

        private NetworkTable_line_class[] Generate_networkTable_lines(Leave_out_class leave_out)
        {
            Dictionary<string, bool> scps_that_shall_not_be_connected_with_eachOther_dict = Generate_dictionary_with_all_signaling_processes();

            leave_out.Order_by_processLevel_descending_newSymbolsRankIncreaseDividedByOriginalProcessSize();
            int leave_out_length = leave_out.Leave_out_lines.Length;
            Leave_out_line_class leave_out_line;
            List<Leave_out_line_class> sameLevel_leave_out_list = new List<Leave_out_line_class>();
            NetworkTable_line_class[] new_networkTable_lines;
            List<NetworkTable_line_class> networkTable_list = new List<NetworkTable_line_class>();
            for (int indexL = 0; indexL < leave_out_length; indexL++)
            {
                leave_out_line = leave_out.Leave_out_lines[indexL];
                if ((indexL == 0) || (!leave_out_line.ProcessLevel.Equals(leave_out.Leave_out_lines[indexL - 1].ProcessLevel)))
                {
                    sameLevel_leave_out_list.Clear();
                }
                sameLevel_leave_out_list.Add(leave_out_line);
                if ((indexL == leave_out_length - 1) || (!leave_out_line.ProcessLevel.Equals(leave_out.Leave_out_lines[indexL + 1].ProcessLevel)))
                {
                    new_networkTable_lines = Generate_networkTable_lines_for_one_level(sameLevel_leave_out_list.ToArray(), scps_that_shall_not_be_connected_with_eachOther_dict);
                    networkTable_list.AddRange(new_networkTable_lines);
                }
            }
            return networkTable_list.ToArray();
        }

        public void Generate_scp_scp_network_from_leave_out(Leave_out_class leave_out)
        {
            NetworkTable_line_class[] networkTable_lines = Generate_networkTable_lines(leave_out);
            Scp_nw = new Network_class();
            Scp_nw.Add_from_networkTable_lines(networkTable_lines);
        }
        #endregion

        #region Generate Scp unions for dynamic enrichment analysis
        public string[][] Generate_array_of_scp_unions_between_any_combination_between_two_or_three_neighboring_selected_scps(string[] consideredSCP_names)
        {
            consideredSCP_names = consideredSCP_names.Distinct().OrderBy(l => l).ToArray();
            Network_class scp_nw_considered = this.Scp_nw.Deep_copy();
            scp_nw_considered.Keep_only_input_nodeNames(consideredSCP_names);
            scp_nw_considered.Transform_into_undirected_double_network();
            int nw_length = scp_nw_considered.NW_length;
            Network_line_class nw_line;
            Network_target_line_class target1_line;
            Network_target_line_class target2_line;
            NetworkNode_line_class current_source_node_line;
            NetworkNode_line_class target1_node_line;
            NetworkNode_line_class target2_node_line;
            string current_source_scp;
            string current_target1_scp;
            string current_target2_scp;
            int target_length;
            List<string[]> list_of_scps_of_each_union = new List<string[]>();
            List<string> scps_of_one_union_list = new List<string>();
            for (int indexNW = 0; indexNW < nw_length; indexNW++)
            {
                nw_line = scp_nw_considered.NW[indexNW];
                target_length = nw_line.Targets.Length;
                current_source_node_line = scp_nw_considered.Nodes.Get_indexed_node_line_if_index_is_correct(indexNW);
                current_source_scp = current_source_node_line.Name;
                for (int indexT1 = 0; indexT1 < target_length; indexT1++)
                {
                    target1_line = nw_line.Targets[indexT1];
                    target1_node_line = scp_nw_considered.Nodes.Get_indexed_node_line_if_index_is_correct(target1_line.NW_index);
                    current_target1_scp = target1_node_line.Name;
                    scps_of_one_union_list.Clear();
                    scps_of_one_union_list.Add(current_source_scp);
                    scps_of_one_union_list.Add(current_target1_scp);
                    scps_of_one_union_list.Add("ZZZZZZZZZZZZZZZZZZZZZ-Remove SCP");
                    list_of_scps_of_each_union.Add(scps_of_one_union_list.OrderBy(l => l).ToArray());
                    for (int indexT2 = indexT1 + 1; indexT2 < target_length; indexT2++)
                    {
                        target2_line = nw_line.Targets[indexT2];
                        target2_node_line = scp_nw_considered.Nodes.Get_indexed_node_line_if_index_is_correct(target2_line.NW_index);
                        current_target2_scp = target2_node_line.Name;
                        scps_of_one_union_list.Clear();
                        scps_of_one_union_list.Add(current_source_scp);
                        scps_of_one_union_list.Add(current_target1_scp);
                        scps_of_one_union_list.Add(current_target2_scp);
                        list_of_scps_of_each_union.Add(scps_of_one_union_list.OrderBy(l => l).ToArray());
                    }
                }
            }
            list_of_scps_of_each_union = list_of_scps_of_each_union.OrderBy(l => l.Length).ThenBy(l => l[0]).ThenBy(l => l[1]).ThenBy(l => l[2]).ToList();
            int list_count = list_of_scps_of_each_union.Count;
            List<string[]> unique_list_of_scps_of_each_union = new List<string[]>();
            string[] current_scp_union;
            for (int indexList = 0; indexList < list_count; indexList++)
            {
                current_scp_union = list_of_scps_of_each_union[indexList];
                if ((indexList == 0)
                     || (!current_scp_union[0].Equals(list_of_scps_of_each_union[indexList - 1][0]))
                     || (!current_scp_union[1].Equals(list_of_scps_of_each_union[indexList - 1][1]))
                     || (!current_scp_union[2].Equals(list_of_scps_of_each_union[indexList - 1][2])))
                {
                    if (current_scp_union[2].IndexOf("ZZZZZZZZZZZZZZZZZZZZZ-Remove SCP") == 0)
                    {
                        unique_list_of_scps_of_each_union.Add(new string[] { current_scp_union[0], current_scp_union[1] });
                    }
                    else
                    {
                        unique_list_of_scps_of_each_union.Add(current_scp_union);
                    }
                }
            }
            return unique_list_of_scps_of_each_union.ToArray();
        }
        #endregion

        public Leave_out_scp_scp_network_class Deep_copy_scp_network()
        {
            Leave_out_scp_scp_network_class copy = (Leave_out_scp_scp_network_class)this.MemberwiseClone();
            copy.Options = this.Options.Deep_copy();
            copy.Scp_nw = this.Scp_nw.Deep_copy();
            return copy;
        }
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

}
