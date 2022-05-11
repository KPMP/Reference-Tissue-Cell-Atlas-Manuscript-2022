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
using System.IO;
using Enumerations;
using ReadWrite;
using Statistic;
using Highthroughput_data;

namespace KPMP
{
    enum KPMP_reference_entityClass_enum { E_m_p_t_y, Reference_cluster, Integrated_cluster, Lmd_rnaseq_segment, Lmd_rnaseq_2nd_segment };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class KPMP_singleNucleusCell_randomization_set_line_class
    {
        public int Randomization_no { get; set; }
        public int Considered_patients_count { get; set; }
        public bool Finished { get; set; }
        public bool Generate_figures_for_visualization { get; set; }
        public string[] Considered_patients { get; set; }
        public string[] Considered_centers { get; set; }

        public string ReadWrite_considered_patients
        {
            get { return ReadWriteClass.Get_writeLine_from_array(this.Considered_patients, KPMP_singleNucleusCell_randomization_set_readWriteOptions_class.Array_delimiter); }
            set { this.Considered_patients = ReadWriteClass.Get_array_from_readLine<string>(value, KPMP_singleNucleusCell_randomization_set_readWriteOptions_class.Array_delimiter); }
        }

        public string ReadWrite_considered_centers
        {
            get { return ReadWriteClass.Get_writeLine_from_array(this.Considered_centers, KPMP_singleNucleusCell_randomization_set_readWriteOptions_class.Array_delimiter); }
            set { this.Considered_centers = ReadWriteClass.Get_array_from_readLine<string>(value, KPMP_singleNucleusCell_randomization_set_readWriteOptions_class.Array_delimiter); }
        }

        public KPMP_singleNucleusCell_randomization_set_line_class()
        {
            Randomization_no = -1;
            Considered_patients_count = -1;
            this.Considered_patients = new string[0];
            this.Considered_centers = new string[0];
        }

        public KPMP_singleNucleusCell_randomization_set_line_class Deep_copy()
        {
            KPMP_singleNucleusCell_randomization_set_line_class copy = (KPMP_singleNucleusCell_randomization_set_line_class)this.MemberwiseClone();
            copy.Considered_centers = Array_class.Deep_copy_string_array(this.Considered_centers);
            copy.Considered_patients = Array_class.Deep_copy_string_array(this.Considered_patients);
            return copy;
        }
    }

    class KPMP_singleNucleusCell_randomization_set_readWriteOptions_class : ReadWriteOptions_base
    {
        public static char Array_delimiter { get { return ';'; } }

        public KPMP_singleNucleusCell_randomization_set_readWriteOptions_class(string subdirectory, string fileName)
        {
            File = Global_directory_class.Results_directory + subdirectory + fileName;
            this.Key_propertyNames = new string[] { "Randomization_no", "Considered_patients_count", "ReadWrite_considered_patients", "ReadWrite_considered_centers", "Generate_figures_for_visualization", "Finished" };
            this.Key_columnNames = this.Key_propertyNames;
            this.HeadlineDelimiters = new char[] { Global_class.Tab };
            this.LineDelimiters = new char[] { Global_class.Tab };
            this.Report = ReadWrite_report_enum.Report_main;
            this.File_has_headline = true;
        }
    }

    class KPMP_singleNucleusCell_randomization_set_class
    {
        public KPMP_singleNucleusCell_randomization_set_line_class[] Randomization_sets { get; set; }

        public KPMP_singleNucleusCell_randomization_set_class()
        {
        }

        public void Generate_by_reading_for_singleCellDataset(string fileName)
        {
            Read_randomization_sets_for_singleCellSequencing(fileName);
        }

        private void Read_randomization_sets_for_singleCellSequencing(string fileName)
        {
            KPMP_singleNucleusCell_randomization_set_readWriteOptions_class readWriteOptions = new KPMP_singleNucleusCell_randomization_set_readWriteOptions_class("", fileName);
            readWriteOptions.File = Global_directory_class.Sample_metadata_additional_datasets + fileName;
            this.Randomization_sets = ReadWriteClass.ReadRawData_and_FillArray<KPMP_singleNucleusCell_randomization_set_line_class>(readWriteOptions);
        }

        public KPMP_singleNucleusCell_randomization_set_class Deep_copy()
        {
            KPMP_singleNucleusCell_randomization_set_class copy = new KPMP_singleNucleusCell_randomization_set_class();
            int rando_length = this.Randomization_sets.Length;
            copy.Randomization_sets = new KPMP_singleNucleusCell_randomization_set_line_class[rando_length];
            for (int indexR=0; indexR<rando_length;indexR++)
            {
                copy.Randomization_sets[indexR] = this.Randomization_sets[indexR].Deep_copy();
            }
            return copy;
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class KPMP_leavePatientsOut_cellType_assignment_line_class
    {
        public string Dataset { get; set; }
        public string Cell_type { get; set; }
        public string Cell_type_second_choice { get; set; }
        public string Cluster { get; set; }
        public int RandoNo { get; set; }
        public bool Is_reference_randoNo { get; set; }
        public int Considered_patients_count { get; set; }
        public float Minus_log10_pvalue { get; set; }
        public float Minus_log10_pvalue_of_second_assignement { get; set; }
        public bool Cell_type_matches_rscript_cell_type { get; set; }

        public static KPMP_leavePatientsOut_cellType_assignment_line_class[] Order_by_randoNo(KPMP_leavePatientsOut_cellType_assignment_line_class[] assigned_cellTypes)
        {
            assigned_cellTypes = assigned_cellTypes.OrderBy(l => l.RandoNo).ToArray();
            return assigned_cellTypes;
        }

        public KPMP_leavePatientsOut_cellType_assignment_line_class Deep_copy()
        {
            KPMP_leavePatientsOut_cellType_assignment_line_class copy = (KPMP_leavePatientsOut_cellType_assignment_line_class)this.MemberwiseClone();
            copy.Cell_type = (string)this.Cell_type.Clone();
            copy.Cluster = (string)this.Cluster.Clone();
            copy.Dataset = (string)this.Dataset.Clone();
            return copy;
        }
    }

    class KPMP_leavePatientsOut_cellType_assignment_readWriteOptions_class : ReadWriteOptions_base
    {
        public KPMP_leavePatientsOut_cellType_assignment_readWriteOptions_class(string subdirectory, string fileName)
        {
            string directory = Global_directory_class.Results_directory + subdirectory;
            ReadWriteClass.Create_directory_if_it_does_not_exist(directory);
            this.File = directory + fileName;
            this.Key_propertyNames = new string[] { "Dataset", "Cell_type", "Cell_type_second_choice", "Cluster", "RandoNo", "Is_reference_randoNo", "Considered_patients_count", "Minus_log10_pvalue", "Minus_log10_pvalue_of_second_assignement", "Cell_type_matches_rscript_cell_type" };
            this.Key_columnNames = this.Key_propertyNames;
            this.HeadlineDelimiters = new char[] { Global_class.Tab };
            this.LineDelimiters = new char[] { Global_class.Tab };
            this.File_has_headline = true;
            this.Report = ReadWrite_report_enum.Report_main;
        }
    }

    class KPMP_leavePatientsOut_cellType_assignment_class
    {
        public KPMP_leavePatientsOut_cellType_assignment_line_class[] Assigned_cellTypes { get; set; }
        public KPMP_leavePatientsOut_options_class Options { get; set; }

        public KPMP_leavePatientsOut_cellType_assignment_class()
        {
            this.Assigned_cellTypes = new KPMP_leavePatientsOut_cellType_assignment_line_class[0];
            this.Options = new KPMP_leavePatientsOut_options_class();
        }

        public Dictionary<int, int> Get_consideredPatientsCount_finishedRandomizationsCount_dictionary()
        {
            Dictionary<int, int> consideredPatientsCount_finishedRandoCount_dict = new Dictionary<int, int>();
            this.Assigned_cellTypes = this.Assigned_cellTypes.OrderBy(l => l.Considered_patients_count).ThenBy(l => l.RandoNo).ToArray();
            KPMP_leavePatientsOut_cellType_assignment_line_class assigned_cellType_line;
            int currentConsideredPatientsCount_randoCount = 0;
            int assigned_cellTypes_length = this.Assigned_cellTypes.Length;
            for (int indexAssigned = 0; indexAssigned < assigned_cellTypes_length; indexAssigned++)
            {
                assigned_cellType_line = this.Assigned_cellTypes[indexAssigned];
                if ((indexAssigned == 0)
                    || (!assigned_cellType_line.Considered_patients_count.Equals(this.Assigned_cellTypes[indexAssigned - 1].Considered_patients_count)))
                {
                    currentConsideredPatientsCount_randoCount = 0;
                }
                if ((indexAssigned == 0)
                    || (!assigned_cellType_line.RandoNo.Equals(this.Assigned_cellTypes[indexAssigned - 1].RandoNo))
                    || (!assigned_cellType_line.Considered_patients_count.Equals(this.Assigned_cellTypes[indexAssigned - 1].Considered_patients_count)))
                {
                    currentConsideredPatientsCount_randoCount++;
                }
                if ((indexAssigned == assigned_cellTypes_length - 1)
                    || (!assigned_cellType_line.Considered_patients_count.Equals(this.Assigned_cellTypes[indexAssigned + 1].Considered_patients_count)))
                {
                    consideredPatientsCount_finishedRandoCount_dict.Add(assigned_cellType_line.Considered_patients_count, currentConsideredPatientsCount_randoCount);
                }
            }
            return consideredPatientsCount_finishedRandoCount_dict;
        }

        private void Set_reference_randoNo()
        {
            foreach (KPMP_leavePatientsOut_cellType_assignment_line_class cellType_assignment_line in this.Assigned_cellTypes)
            {
                if (cellType_assignment_line.RandoNo == Options.Reference_randoNO)
                {
                    cellType_assignment_line.Is_reference_randoNo = true;
                }
                else
                {
                    cellType_assignment_line.Is_reference_randoNo = false;
                }
            }
        }

        private void Set_considered_patients_count(KPMP_singleNucleusCell_randomization_set_class randomization_documentation)
        {
            this.Assigned_cellTypes = KPMP_leavePatientsOut_cellType_assignment_line_class.Order_by_randoNo(this.Assigned_cellTypes);
            int assigned_cellTypes_length = this.Assigned_cellTypes.Length;
            KPMP_leavePatientsOut_cellType_assignment_line_class assigned_line;
            int rando_length = randomization_documentation.Randomization_sets.Length;
            KPMP_singleNucleusCell_randomization_set_line_class rando_line = new KPMP_singleNucleusCell_randomization_set_line_class();
            int indexRando = 0;
            int stringCompare = -2;
            for (int indexA = 0; indexA < assigned_cellTypes_length; indexA++)
            {
                assigned_line = this.Assigned_cellTypes[indexA];
                stringCompare = -2;
                while ((indexRando < rando_length) && (stringCompare < 0))
                {
                    rando_line = randomization_documentation.Randomization_sets[indexRando];
                    stringCompare = rando_line.Randomization_no.CompareTo(assigned_line.RandoNo);
                    if (stringCompare < 0)
                    {
                        indexRando++;
                    }
                }
                if (stringCompare != 0) { throw new Exception(); }
                assigned_line.Considered_patients_count = rando_line.Considered_patients_count;
            }
        }

        private void Generate_from_barcode_cellType_instance_private(KPMP_barcode_cellType_class barcode_cellType)
        {
            int barcode_cellTypes_length = barcode_cellType.Barcodes.Length;
            KPMP_barcode_cellType_line_class barcode_celltype_line;
            KPMP_barcode_cellType_line_class last_added_barcode_celltype_line = new KPMP_barcode_cellType_line_class();
            barcode_cellType.Barcodes = KPMP_barcode_cellType_line_class.Order_by_randoNo_cellType_clusterNo(barcode_cellType.Barcodes);
            KPMP_leavePatientsOut_cellType_assignment_line_class cellType_assignement_line;
            List<KPMP_leavePatientsOut_cellType_assignment_line_class> cellType_assignements = new List<KPMP_leavePatientsOut_cellType_assignment_line_class>();
            List<float> minusLogPvalue_first = new List<float>();
            List<float> minusLogPvalue_second = new List<float>();
            for (int indexBarcode = 0; indexBarcode < barcode_cellTypes_length; indexBarcode++)
            {
                barcode_celltype_line = barcode_cellType.Barcodes[indexBarcode];
                if ((indexBarcode == 0)
                    || (!barcode_celltype_line.Rando_no.Equals(barcode_cellType.Barcodes[indexBarcode - 1].Rando_no))
                    || (!barcode_celltype_line.Cell_type.Equals(barcode_cellType.Barcodes[indexBarcode - 1].Cell_type)))
                {
                    minusLogPvalue_first.Clear();
                    minusLogPvalue_second.Clear();
                }
                if ((indexBarcode == 0)
                    || (!barcode_celltype_line.Rando_no.Equals(barcode_cellType.Barcodes[indexBarcode - 1].Rando_no))
                    || (!barcode_celltype_line.Cluster_no.Equals(barcode_cellType.Barcodes[indexBarcode - 1].Cluster_no))
                    || (!barcode_celltype_line.Cell_type.Equals(barcode_cellType.Barcodes[indexBarcode - 1].Cell_type)))
                {
                    last_added_barcode_celltype_line = barcode_celltype_line;
                    minusLogPvalue_first.Add(barcode_celltype_line.Cell_type_1st_minLog10Pvalue);
                    minusLogPvalue_second.Add(barcode_celltype_line.Cell_type_2nd_minLog10Pvalue);
                }
                if (!barcode_celltype_line.Cell_type_2nd_minLog10Pvalue.Equals(last_added_barcode_celltype_line.Cell_type_2nd_minLog10Pvalue)) { throw new Exception(); }
                if (!barcode_celltype_line.Cell_type_1st_minLog10Pvalue.Equals(last_added_barcode_celltype_line.Cell_type_1st_minLog10Pvalue)) { throw new Exception(); }
                if ((indexBarcode == barcode_cellTypes_length - 1)
                    || (!barcode_celltype_line.Rando_no.Equals(barcode_cellType.Barcodes[indexBarcode + 1].Rando_no))
                    || (!barcode_celltype_line.Cell_type.Equals(barcode_cellType.Barcodes[indexBarcode + 1].Cell_type)))
                {
                    cellType_assignement_line = new KPMP_leavePatientsOut_cellType_assignment_line_class();
                    cellType_assignement_line.Cell_type = (string)barcode_celltype_line.Cell_type.Clone();
                    cellType_assignement_line.Cluster = (string)barcode_celltype_line.Cell_type.Clone();
                    cellType_assignement_line.RandoNo = barcode_celltype_line.Rando_no;
                    cellType_assignement_line.Minus_log10_pvalue = Math_class.Get_average(minusLogPvalue_first.ToArray());
                    cellType_assignement_line.Cell_type_second_choice = "One or more";
                    cellType_assignement_line.Minus_log10_pvalue_of_second_assignement = Math_class.Get_average(minusLogPvalue_second.ToArray());
                    cellType_assignements.Add(cellType_assignement_line);
                }
            }
            this.Assigned_cellTypes = cellType_assignements.ToArray();
        }

        private void Remove_not_assigned_cellTypes()
        {
            List<KPMP_leavePatientsOut_cellType_assignment_line_class> keep = new List<KPMP_leavePatientsOut_cellType_assignment_line_class>();
            foreach (KPMP_leavePatientsOut_cellType_assignment_line_class cellType_assignment_line in this.Assigned_cellTypes)
            {
                if (!cellType_assignment_line.Cell_type.Equals(KPMP_data_integration_class.Get_not_assigned_label()))
                {
                    keep.Add(cellType_assignment_line);
                }
            }
            this.Assigned_cellTypes = keep.ToArray();
        }

        public void Generate_from_barcode_cellType_instance(KPMP_barcode_cellType_class barcode_cellType, KPMP_singleNucleusCell_randomization_set_class randomization_documentation)
        {
            Generate_from_barcode_cellType_instance_private(barcode_cellType);
            Set_reference_randoNo();
            Remove_not_assigned_cellTypes();
            Set_considered_patients_count(randomization_documentation);
        }

        public void Write(string subdirectory, string fileName)
        {
            KPMP_leavePatientsOut_cellType_assignment_readWriteOptions_class readWriteOptions = new KPMP_leavePatientsOut_cellType_assignment_readWriteOptions_class(subdirectory, fileName);
            ReadWriteClass.WriteData(this.Assigned_cellTypes, readWriteOptions);
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class KPMP_leavePatientsOut_options_class
    {
        public int Reference_randoNO { get; set; }
        public bool Consider_randoClusters_without_matching_referenceClusters_for_analysis { get; set; }
        public bool Consider_referenceClusters_without_matching_randoClusters_for_analysis { get; set; }
        public int Top_degs_for_cell_type_assignment { get; set; }

        public bool Equals_oter_options(KPMP_leavePatientsOut_options_class other)
        {
            bool equal_options = true;
            if (this.Reference_randoNO!=other.Reference_randoNO) { equal_options = false; }
            if (this.Consider_randoClusters_without_matching_referenceClusters_for_analysis!=other.Consider_randoClusters_without_matching_referenceClusters_for_analysis) { equal_options = false; }
            if (this.Consider_referenceClusters_without_matching_randoClusters_for_analysis!=other.Consider_referenceClusters_without_matching_randoClusters_for_analysis) { equal_options = false; }
            return equal_options;
        }

        public KPMP_leavePatientsOut_options_class()
        {
            this.Reference_randoNO = 0;
            Top_degs_for_cell_type_assignment = 300;
            Consider_randoClusters_without_matching_referenceClusters_for_analysis = false;
            Consider_referenceClusters_without_matching_randoClusters_for_analysis = false;
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class KPMP_randomization_name_class
    {
        const string multiple_assignments_label = "Multiple assignments";
        const string no_assignments_label = "No assignments";

        public static string Multiple_assignments_label { get { return multiple_assignments_label; } }

        public static string No_assignments_label { get { return no_assignments_label; } }

        public static bool Is_unique_assignment(string assignment)
        {
            switch (assignment)
            {
                case multiple_assignments_label:
                case no_assignments_label:
                    return false;
                default:
                    return true;
            }
        }

        public static int[] Get_all_finished_and_unfinished_randomization_nos(string subdirectory, out int[] uncompleted_randoNos)
        {
            string[] mandatory_baseFileNames = new string[] { "Barcode_cellType_associations", "CellCounts_PercentCellType_perPatient", "CellCounts_PercentPatient_perCellType", "CellType_detection_minLog10PDiffs_resolutions_SCT_data", "CellType_detection_resolutions_SCT_data", "CellType_tissueCollection", "CellType_tissueType", "CellTypeDEGsAfterMerging_SCT_data", "ClusterDEGsSelectedRes_SCT_data", "Quality_control" };

            mandatory_baseFileNames = mandatory_baseFileNames.OrderBy(l => l).ToArray();
            string mandatory_baseFileName;
            int mandatrory_baseFileNames_length = mandatory_baseFileNames.Length;
            Dictionary<int, Dictionary<string, List<float>>> randoNo_baseFileName_resolutions_dict = new Dictionary<int, Dictionary<string, List<float>>>();
            Dictionary<string, List<float>> current_baseFileName_resolutions_dict;
            List<float> current_resolutions;
            string completeDirectory = Global_directory_class.Seurat_singleCell_cluster_directory + subdirectory;
            string[] all_completeFileNames = Directory.GetFiles(completeDirectory);
            string completeFileName;
            int completeFileNames_length = all_completeFileNames.Length;
            float resolution =-1;
            int maxPrincipalComponent;
            string baseFileName;
            int randoNo;
            for (int indexC = 0; indexC < completeFileNames_length; indexC++)
            {
                completeFileName = all_completeFileNames[indexC];
                if (completeFileName.IndexOf("randoNo1009")!=-1)
                {
                    string lol = "";
                }
                if (  (completeFileName.IndexOf("AA_Analysis_finished_") != -1)
                    || (completeFileName.IndexOf("AA_Analysis_started_") != -1)
                    || (completeFileName.IndexOf("AA_Analysis_interrupted_") != -1)
                    || (completeFileName.IndexOf("_PCA_genes") != -1))
                { }
                else
                {
                    randoNo = Get_baseFileName_randomization_number_maxPrincipalComponent_and_resolution_from_fileName(completeFileName, out baseFileName, out maxPrincipalComponent);
                    if (!randoNo_baseFileName_resolutions_dict.ContainsKey(randoNo))
                    {
                        randoNo_baseFileName_resolutions_dict.Add(randoNo, new Dictionary<string, List<float>>());
                    }
                    if (!randoNo_baseFileName_resolutions_dict[randoNo].ContainsKey(baseFileName))
                    {
                        randoNo_baseFileName_resolutions_dict[randoNo].Add(baseFileName, new List<float>());
                    }
                    randoNo_baseFileName_resolutions_dict[randoNo][baseFileName].Add(resolution);
                }
            }

            int[] all_randos = randoNo_baseFileName_resolutions_dict.Keys.ToArray();
            int current_rando;
            int all_randos_length = all_randos.Length;
            string[] current_baseFileNames;
            string current_baseFileName = "";
            int current_baseFileNames_length;
            int indexCurrentBaseFileName = 0;
            int baseFileName_compare = -2;

            bool consider_current_rando;

            List<int> completed_randoNos_list = new List<int>();
            List<int> uncompleted_randoNos_list = new List<int>();
            for (int indexRando = 0; indexRando < all_randos_length; indexRando++)
            {
                current_rando = all_randos[indexRando];
                consider_current_rando = true;
                current_baseFileName_resolutions_dict = randoNo_baseFileName_resolutions_dict[current_rando];
                current_baseFileNames = current_baseFileName_resolutions_dict.Keys.OrderBy(l => l).ToArray();
                current_baseFileNames_length = current_baseFileNames.Length;
                for (int indexMandatoryBase = 0; indexMandatoryBase < mandatrory_baseFileNames_length; indexMandatoryBase++)
                {
                    mandatory_baseFileName = mandatory_baseFileNames[indexMandatoryBase];
                    baseFileName_compare = -2;
                    indexCurrentBaseFileName = 0;
                    while ((indexCurrentBaseFileName < current_baseFileNames_length) && (baseFileName_compare < 0))
                    {
                        current_baseFileName = current_baseFileNames[indexCurrentBaseFileName];
                        baseFileName_compare = current_baseFileName.CompareTo(mandatory_baseFileName);
                        if (baseFileName_compare < 0)
                        {
                            indexCurrentBaseFileName++;
                        }
                    }
                    if (baseFileName_compare != 0)
                    {
                        consider_current_rando = false;
                        break;
                    }
                }
                if (consider_current_rando)
                {
                    completed_randoNos_list.Add(current_rando);
                }
                else
                {
                    uncompleted_randoNos_list.Add(current_rando);
                }
            }
            uncompleted_randoNos = uncompleted_randoNos_list.OrderBy(l => l).ToArray();
            return completed_randoNos_list.OrderBy(l => l).ToArray();
        }
 
        public static int Get_baseFileName_randomization_number_maxPrincipalComponent_and_resolution_from_fileName(string fileName, out string baseFileName, out int max_principal_component)
        {
            string fileName_withoutExtension = Path.GetFileNameWithoutExtension(fileName);
            string[] splitStrings = fileName_withoutExtension.Split('_');
            string splitString;
            int splitStrings_length = splitStrings.Length;
            string randoNo_string = "RandoNo";
            string maxPrincipalComponent_string = "maxPC";
            max_principal_component = -1;
            int randoNumber = -1;
            StringBuilder sb = new StringBuilder();
            for (int indexS = 0; indexS < splitStrings_length; indexS++)
            {
                splitString = splitStrings[indexS];
                if (splitString.ToUpper().IndexOf(randoNo_string.ToUpper()) == 0)
                {
                    int indexRandoNo = splitString.ToUpper().IndexOf(randoNo_string.ToUpper()) + randoNo_string.Length;
                    randoNumber = int.Parse(splitString.Substring(indexRandoNo, splitString.Length - indexRandoNo));
                }
                else if (splitString.IndexOf(maxPrincipalComponent_string) == 0)
                {
                    int indexMaxPrincipalComponent_string = splitString.IndexOf(maxPrincipalComponent_string) + maxPrincipalComponent_string.Length;
                    max_principal_component = int.Parse(splitString.Substring(indexMaxPrincipalComponent_string, splitString.Length - indexMaxPrincipalComponent_string));
                }
                else
                {
                    if (sb.Length > 0) { sb.AppendFormat("_"); }
                    sb.AppendFormat("{0}", splitString);
                }
            }
            baseFileName = sb.ToString();
            if (randoNumber == -1) { throw new Exception(); }
            if (max_principal_component == -1) { throw new Exception(); }
            return randoNumber;
        }

        public static string Get_fileName_without_resolution(string fileName)
        {
            string fileName_withoutExtension = Path.GetFileNameWithoutExtension(fileName);
            string[] splitStrings = fileName_withoutExtension.Split('_');
            string splitString;
            int splitStrings_length = splitStrings.Length;
            string resolution_string = "res";
            StringBuilder sb = new StringBuilder();
            for (int indexS = 0; indexS < splitStrings_length; indexS++)
            {
                splitString = splitStrings[indexS];
                if (splitString.IndexOf(resolution_string) == -1)
                {
                    if (sb.Length > 0) { sb.AppendFormat("_"); }
                    sb.AppendFormat(splitString);
                }

            }
            return sb.ToString();
        }

        public static int Get_randomization_number_from_fileName(string fileName)
        {
            string fileName_withoutExtension = Path.GetFileNameWithoutExtension(fileName);
            string pos2;
            int pos3;
            int randoNumber = KPMP_randomization_name_class.Get_baseFileName_randomization_number_maxPrincipalComponent_and_resolution_from_fileName(fileName, out pos2, out pos3);
            return randoNumber;
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class KPMP_barcode_cellType_line_class
    {
        public string Barcode { get; set; }
        public string Cell_type { get; set; }
        public int Cluster_no { get; set; }
        public int Rando_no { get; set; }
        public string Lcm_segment { get; set; }
        public string Integrated_cluster { get; set; }
        public int Considered_patients_count { get; set; }
        public float Resolution { get; set; }
        public float Cell_type_1st_minLog10Pvalue { get; set; }
        public string Cell_type_2nd { get; set; }
        public float Cell_type_2nd_minLog10Pvalue { get; set; }

        public static KPMP_barcode_cellType_line_class[] Order_by_randoNo_cellType_clusterNo(KPMP_barcode_cellType_line_class[] barcode_lines)
        {
            Dictionary<int, Dictionary<string, Dictionary<int, List<KPMP_barcode_cellType_line_class>>>> randoNo_cellType_clusterNo_dict = new Dictionary<int, Dictionary<string, Dictionary<int, List<KPMP_barcode_cellType_line_class>>>>();
            Dictionary<string, Dictionary<int, List<KPMP_barcode_cellType_line_class>>> cellType_clusterNo_dict = new Dictionary<string, Dictionary<int, List<KPMP_barcode_cellType_line_class>>>();
            Dictionary<int, List<KPMP_barcode_cellType_line_class>> clusterNo_dict = new Dictionary<int, List<KPMP_barcode_cellType_line_class>>();
            KPMP_barcode_cellType_line_class barcode_line;
            int barcode_lines_length = barcode_lines.Length;
            for (int indexBarcode = 0; indexBarcode < barcode_lines_length; indexBarcode++)
            {
                barcode_line = barcode_lines[indexBarcode];
                if (!randoNo_cellType_clusterNo_dict.ContainsKey(barcode_line.Rando_no))
                {
                    randoNo_cellType_clusterNo_dict.Add(barcode_line.Rando_no, new Dictionary<string, Dictionary<int, List<KPMP_barcode_cellType_line_class>>>());
                }
                if (!randoNo_cellType_clusterNo_dict[barcode_line.Rando_no].ContainsKey(barcode_line.Cell_type))
                {
                    randoNo_cellType_clusterNo_dict[barcode_line.Rando_no].Add(barcode_line.Cell_type, new Dictionary<int, List<KPMP_barcode_cellType_line_class>>());
                }
                if (!randoNo_cellType_clusterNo_dict[barcode_line.Rando_no][barcode_line.Cell_type].ContainsKey(barcode_line.Cluster_no))
                {
                    randoNo_cellType_clusterNo_dict[barcode_line.Rando_no][barcode_line.Cell_type].Add(barcode_line.Cluster_no, new List<KPMP_barcode_cellType_line_class>());
                }
                randoNo_cellType_clusterNo_dict[barcode_line.Rando_no][barcode_line.Cell_type][barcode_line.Cluster_no].Add(barcode_line);
            }
            int[] randoNos = randoNo_cellType_clusterNo_dict.Keys.ToArray();
            int randoNo;
            int randoNos_length = randoNos.Length;
            string[] cellTypes;
            string cellType;
            int cellTypes_length;
            int[] clusterNos;
            int clusterNo;
            int clusterNos_length;
            randoNos = randoNos.OrderBy(l => l).ToArray();
            List<KPMP_barcode_cellType_line_class> ordered_barcodes = new List<KPMP_barcode_cellType_line_class>();
            for (int indexRandoNos = 0; indexRandoNos < randoNos_length; indexRandoNos++)
            {
                randoNo = randoNos[indexRandoNos];
                cellType_clusterNo_dict = randoNo_cellType_clusterNo_dict[randoNo];
                cellTypes = cellType_clusterNo_dict.Keys.ToArray();
                cellTypes_length = cellTypes.Length;
                cellTypes = cellTypes.OrderBy(l => l).ToArray();
                for (int indexCT = 0; indexCT < cellTypes_length; indexCT++)
                {
                    cellType = cellTypes[indexCT];
                    clusterNo_dict = cellType_clusterNo_dict[cellType];
                    clusterNos = clusterNo_dict.Keys.ToArray();
                    clusterNos = clusterNos.OrderBy(l => l).ToArray();
                    clusterNos_length = clusterNos.Length;
                    for (int indexCNO = 0; indexCNO < clusterNos_length; indexCNO++)
                    {
                        clusterNo = clusterNos[indexCNO];
                        ordered_barcodes.AddRange(clusterNo_dict[clusterNo]);
                    }
                }
            }

            if (Global_class.Check_ordering)
            {
                int ordered_length = ordered_barcodes.Count;
                if (ordered_length != barcode_lines_length) { throw new Exception(); }
                KPMP_barcode_cellType_line_class this_ordered_line;
                KPMP_barcode_cellType_line_class previous_ordered_line;
                for (int indexO = 1; indexO < ordered_length; indexO++)
                {
                    previous_ordered_line = ordered_barcodes[indexO - 1];
                    this_ordered_line = ordered_barcodes[indexO];
                    if ((this_ordered_line.Rando_no.CompareTo(previous_ordered_line.Rando_no)) < 0) { throw new Exception(); }
                    else if ((this_ordered_line.Rando_no.Equals(previous_ordered_line.Rando_no))
                             && (this_ordered_line.Cell_type.CompareTo(previous_ordered_line.Cell_type)) < 0) { throw new Exception(); }
                    else if ((this_ordered_line.Rando_no.Equals(previous_ordered_line.Rando_no))
                             && (this_ordered_line.Cell_type.Equals(previous_ordered_line.Cell_type))
                             && (this_ordered_line.Cluster_no.CompareTo(previous_ordered_line.Cluster_no)) < 0) { throw new Exception(); }
                }
            }
            return ordered_barcodes.ToArray();
        }

        public KPMP_barcode_cellType_line_class Deep_copy()
        {
            KPMP_barcode_cellType_line_class copy = (KPMP_barcode_cellType_line_class)this.MemberwiseClone();
            copy.Barcode = (string)this.Barcode.Clone();
            copy.Cell_type = (string)this.Cell_type.Clone();
            copy.Cell_type_2nd = (string)this.Cell_type_2nd.Clone();
            copy.Lcm_segment = (string)this.Lcm_segment.Clone();
            copy.Integrated_cluster = (string)this.Integrated_cluster.Clone();
            return copy;
        }
    }

    class KPMP_barcode_cellType_inputReadWriteOptions_class : ReadWriteOptions_base
    {
        public KPMP_barcode_cellType_inputReadWriteOptions_class(string completeFileName)
        {
            this.File = completeFileName;
            this.Key_propertyNames = new string[] { "Barcode", "Cluster_no", "Cell_type", "Resolution", "Cell_type_1st_minLog10Pvalue", "Cell_type_2nd", "Cell_type_2nd_minLog10Pvalue" };
            this.Key_columnNames = new string[] { "Barcode", "Cluster", "Cell_type", "Resolution", "Cell_type_1st_minLog10Pvalue", "Cell_type_2nd", "Cell_type_2nd_minLog10Pvalue" };
            this.HeadlineDelimiters = new char[] { Global_class.Tab };
            this.LineDelimiters = new char[] { Global_class.Tab };
            this.File_has_headline = true;
            this.Report = ReadWrite_report_enum.Report_main;
        }
    }

    class KPMP_barcode_cellType_class
    {
        public KPMP_barcode_cellType_line_class[] Barcodes { get; set; }
        public KPMP_leavePatientsOut_options_class Options { get; set; }

        public KPMP_barcode_cellType_class()
        {
            this.Options = new KPMP_leavePatientsOut_options_class();
        }

        private void Rearrange_cellType_order_if_multiple_cellTypes_discovered()
        {
            Dictionary<string, string> oldCellType_newCellType_dict = new Dictionary<string, string>();
            int barcodes_length = this.Barcodes.Length;
            KPMP_barcode_cellType_line_class barcode_line;
            string new_cellType;
            for (int indexB = 0; indexB < barcodes_length; indexB++)
            {
                barcode_line = this.Barcodes[indexB];
                if (!oldCellType_newCellType_dict.ContainsKey(barcode_line.Cell_type))
                {
                    new_cellType = KPMP_data_integration_class.Rearrange_cellType_order_if_multiple_cellTypes_discovered(barcode_line.Cell_type);
                    oldCellType_newCellType_dict.Add((string)barcode_line.Cell_type.Clone(), (string)new_cellType.Clone());
                }
                barcode_line.Cell_type = (string)oldCellType_newCellType_dict[barcode_line.Cell_type].Clone();
            }
        }

        public void Generate(string subdirectory, string add_to_file)
        {
            Read_barcode_clusters(subdirectory, add_to_file);
            Rearrange_cellType_order_if_multiple_cellTypes_discovered();
        }

        private void Read_barcode_clusters(string subdirectory, string add_to_file)
        {
            string complete_directory = Global_directory_class.PostHoc_power_seurat_results_directory + subdirectory;
            string[] completeFileNames = Directory.GetFiles(complete_directory);
            string completeFileName;
            int completeFileNames_length = completeFileNames.Length;
            KPMP_barcode_cellType_line_class[] add_barcodes;
            List<KPMP_barcode_cellType_line_class> barcodes = new List<KPMP_barcode_cellType_line_class>();
            int randomizationNumber;
            int maxPrincipalComponent;
            string baseFileName;
            for (int indexC = 0; indexC < completeFileNames_length; indexC++)
            {
                completeFileName = completeFileNames[indexC];
                if ((completeFileName.IndexOf("Barcode_cellType_associations") != -1)
                    && (completeFileName.IndexOf(add_to_file) != -1))
                {
                    KPMP_barcode_cellType_inputReadWriteOptions_class readWriteOptions = new KPMP_barcode_cellType_inputReadWriteOptions_class(completeFileName);
                    add_barcodes = ReadWriteClass.ReadRawData_and_FillArray<KPMP_barcode_cellType_line_class>(readWriteOptions);
                    randomizationNumber = KPMP_randomization_name_class.Get_baseFileName_randomization_number_maxPrincipalComponent_and_resolution_from_fileName(completeFileName, out baseFileName, out maxPrincipalComponent);
                    foreach (KPMP_barcode_cellType_line_class cluster_line in add_barcodes)
                    {
                        cluster_line.Rando_no = randomizationNumber;
                    }
                    barcodes.AddRange(add_barcodes);
                }
            }
            this.Barcodes = barcodes.ToArray();
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class KPMP_leavePatientsOut_analysis_line_class
    {
        public string EntityClassName { get; set; }
        public int Considered_patients { get; set; }
        public double Normalized_overlap_average { get; set; }
        public double Normalized_overlap_sampleSD { get; set; }
        public double Normalized_notOverlap_average { get; set; }
        public double Normalized_notOverlap_sampleSD { get; set; }
        public double Overlap_average { get; set; }
        public double Overlap_sampleSD { get; set; }
        public double NotOverlap_average { get; set; }
        public double NotOverlap_sampleSD { get; set; }
        public double JaccardIndex_average { get; set; }
        public double JaccardIndex_sampleSD { get; set; }
        public double Correlation_average { get; set; }
        public double Correlation_sampleSD { get; set; }
        public int RandomizationSets_count { get; set; }
        public float Percent_randomizationSets_count { get; set; }
        public int Max_principal_component { get; set; }
        public float Cluster_resolution { get; set; }
        public int Set_no { get; set; }

        public KPMP_leavePatientsOut_analysis_line_class()
        {
            this.Cluster_resolution = 0;
            Considered_patients = 0;
            Normalized_notOverlap_average = 0;
            Normalized_overlap_average = 0;
            Overlap_average = 0;
            NotOverlap_average = 0;
            JaccardIndex_average = 0;
            Correlation_average = 0;
            RandomizationSets_count = 0;
            Percent_randomizationSets_count = 0;
        }

        public KPMP_leavePatientsOut_analysis_line_class Deep_copy()
        {
            KPMP_leavePatientsOut_analysis_line_class copy = (KPMP_leavePatientsOut_analysis_line_class)this.MemberwiseClone();
            copy.EntityClassName = (string)this.EntityClassName.Clone();
            return copy;
        }
    }

    class KPMP_leavePatientsOut_analysis_readWriteOptions : ReadWriteOptions_base
    {
        public KPMP_leavePatientsOut_analysis_readWriteOptions(string subdirectory, string fileName)
        {
            this.File = Global_directory_class.Results_directory + subdirectory + fileName;
            this.Key_propertyNames = new string[] { "EntityClassName","Considered_patients","Normalized_overlap_average","Normalized_overlap_sampleSD","Normalized_notOverlap_average",
                                                    "Normalized_notOverlap_sampleSD", "Overlap_average", "Overlap_sampleSD", "NotOverlap_average", "NotOverlap_sampleSD", "JaccardIndex_average",
                                                    "JaccardIndex_sampleSD", "Correlation_average", "Correlation_sampleSD", "RandomizationSets_count","Percent_randomizationSets_count","Set_no",
                                                    "Max_principal_component","Cluster_resolution" };
            this.Key_columnNames = this.Key_propertyNames;
            this.File_has_headline = true;
            this.HeadlineDelimiters = new char[] { Global_class.Tab };
            this.LineDelimiters = new char[] { Global_class.Tab };
            this.Report = ReadWrite_report_enum.Report_main;
        }
    }

    class KPMP_leavePatientsOut_analysis_class
    {
        public KPMP_leavePatientsOut_analysis_line_class[] Analysis { get; set; }
        public KPMP_leavePatientsOut_options_class Options { get; set; }

        public KPMP_leavePatientsOut_analysis_class()
        {
            this.Analysis = new KPMP_leavePatientsOut_analysis_line_class[0];
            this.Options = new KPMP_leavePatientsOut_options_class();
        }

        private void Add_to_array(KPMP_leavePatientsOut_analysis_line_class[] add_analysis)
        {
            int this_length = this.Analysis.Length;
            int add_length = add_analysis.Length;
            int new_length = this_length + add_length;
            KPMP_leavePatientsOut_analysis_line_class[] new_analysis = new KPMP_leavePatientsOut_analysis_line_class[new_length];
            int indexNew = -1;
            for (int indexThis = 0; indexThis < this_length; indexThis++)
            {
                indexNew++;
                new_analysis[indexNew] = this.Analysis[indexThis];
            }
            for (int indexAdd = 0; indexAdd < add_length; indexAdd++)
            {
                indexNew++;
                new_analysis[indexNew] = add_analysis[indexAdd];
            }
            this.Analysis = new_analysis;
        }

        public void Generate_fisherExact_from_cell_type_assignments(KPMP_leavePatientsOut_cellType_assignment_class cell_type_assignment)
        {
            Dictionary<int, int> consideredPatientsCount_finishedRandoCount = cell_type_assignment.Get_consideredPatientsCount_finishedRandomizationsCount_dictionary();

            KPMP_leavePatientsOut_cellType_assignment_line_class cell_type_assignment_line;
            List<float> currentConsideredPatients_currentCluster_firstMinusLog10pvalues = new List<float>();
            List<float> currentConsideredPatients_currentCluster_secondMinusLog10pvalues = new List<float>();
            int currentConsideredPatients_randoCount = 0;
            KPMP_leavePatientsOut_analysis_line_class analysis_line;
            List<KPMP_leavePatientsOut_analysis_line_class> analysis_list = new List<KPMP_leavePatientsOut_analysis_line_class>();
            int cell_type_assignments_length = cell_type_assignment.Assigned_cellTypes.Length;
            cell_type_assignment.Assigned_cellTypes = cell_type_assignment.Assigned_cellTypes.OrderBy(l => l.Considered_patients_count).ThenBy(l => l.Cell_type).ThenBy(l => l.RandoNo).ToArray();
            float average;
            float sampleSD;
            for (int indexCellTypeAssignment = 0; indexCellTypeAssignment < cell_type_assignments_length; indexCellTypeAssignment++)
            {
                cell_type_assignment_line = cell_type_assignment.Assigned_cellTypes[indexCellTypeAssignment];
                if ((indexCellTypeAssignment == 0)
                    || (!cell_type_assignment_line.Considered_patients_count.Equals(cell_type_assignment.Assigned_cellTypes[indexCellTypeAssignment - 1].Considered_patients_count))
                    || (!cell_type_assignment_line.Cell_type.Equals(cell_type_assignment.Assigned_cellTypes[indexCellTypeAssignment - 1].Cell_type)))
                {
                    currentConsideredPatients_currentCluster_firstMinusLog10pvalues.Clear();
                    currentConsideredPatients_currentCluster_secondMinusLog10pvalues.Clear();
                    currentConsideredPatients_randoCount = 0;
                }
                if ((indexCellTypeAssignment == 0)
                    || (!cell_type_assignment_line.RandoNo.Equals(cell_type_assignment.Assigned_cellTypes[indexCellTypeAssignment - 1].RandoNo))
                    || (!cell_type_assignment_line.Considered_patients_count.Equals(cell_type_assignment.Assigned_cellTypes[indexCellTypeAssignment - 1].Considered_patients_count))
                    || (!cell_type_assignment_line.Cell_type.Equals(cell_type_assignment.Assigned_cellTypes[indexCellTypeAssignment - 1].Cell_type)))
                {
                    currentConsideredPatients_randoCount++;
                }
                currentConsideredPatients_currentCluster_firstMinusLog10pvalues.Add(cell_type_assignment_line.Minus_log10_pvalue);
                currentConsideredPatients_currentCluster_secondMinusLog10pvalues.Add(cell_type_assignment_line.Minus_log10_pvalue_of_second_assignement);
                if ((indexCellTypeAssignment == cell_type_assignments_length - 1)
                    || (!cell_type_assignment_line.Considered_patients_count.Equals(cell_type_assignment.Assigned_cellTypes[indexCellTypeAssignment + 1].Considered_patients_count))
                    || (!cell_type_assignment_line.Cell_type.Equals(cell_type_assignment.Assigned_cellTypes[indexCellTypeAssignment + 1].Cell_type)))
                {
                    analysis_line = new KPMP_leavePatientsOut_analysis_line_class();
                    analysis_line.EntityClassName = "Fisher exact - " + cell_type_assignment_line.Cell_type;
                    analysis_line.Set_no = 0;
                    analysis_line.Considered_patients = cell_type_assignment_line.Considered_patients_count;

                    Math_class.Get_mean_and_sample_sd(currentConsideredPatients_currentCluster_firstMinusLog10pvalues.ToArray(), out average, out sampleSD);
                    if (currentConsideredPatients_currentCluster_firstMinusLog10pvalues.Distinct().ToArray().Length == 1) { sampleSD = 0; }
                    if (Double.IsNaN(sampleSD)) { throw new Exception(); }
                    analysis_line.Overlap_average = average;
                    analysis_line.Overlap_sampleSD = sampleSD;

                    analysis_line.Normalized_overlap_average = 0;
                    analysis_line.Normalized_overlap_sampleSD = 0;

                    analysis_line.NotOverlap_average = 0;
                    analysis_line.NotOverlap_sampleSD = 0;

                    analysis_line.Normalized_notOverlap_average = 0;
                    analysis_line.Normalized_notOverlap_sampleSD = 0;

                    analysis_line.JaccardIndex_average = 0;
                    analysis_line.JaccardIndex_sampleSD = 0;

                    analysis_line.RandomizationSets_count = currentConsideredPatients_currentCluster_firstMinusLog10pvalues.Count;

                    analysis_list.Add(analysis_line);

                    analysis_line = new KPMP_leavePatientsOut_analysis_line_class();
                    analysis_line.EntityClassName = "Fisher exact - " + cell_type_assignment_line.Cell_type;
                    analysis_line.Set_no = 1;
                    analysis_line.Considered_patients = cell_type_assignment_line.Considered_patients_count;

                    Math_class.Get_mean_and_sample_sd(currentConsideredPatients_currentCluster_secondMinusLog10pvalues.ToArray(), out average, out sampleSD);
                    analysis_line.Overlap_average = average;
                    analysis_line.Overlap_sampleSD = sampleSD;

                    analysis_line.Normalized_overlap_average = 0;
                    analysis_line.Normalized_overlap_sampleSD = 0;

                    analysis_line.NotOverlap_average = 0;
                    analysis_line.NotOverlap_sampleSD = 0;

                    analysis_line.Normalized_notOverlap_average = 0;
                    analysis_line.Normalized_notOverlap_sampleSD = 0;

                    analysis_line.JaccardIndex_average = 0;
                    analysis_line.JaccardIndex_sampleSD = 0;

                    analysis_line.RandomizationSets_count = currentConsideredPatients_randoCount;
                    analysis_line.Percent_randomizationSets_count = 100F * (float)analysis_line.RandomizationSets_count / (float)consideredPatientsCount_finishedRandoCount[analysis_line.Considered_patients];
                    if (analysis_line.Percent_randomizationSets_count > 100)
                    {
                        float finishedRando = consideredPatientsCount_finishedRandoCount[analysis_line.Considered_patients];
                        throw new Exception();
                    }

                    analysis_list.Add(analysis_line);
                }
            }
            Add_to_array(analysis_list.ToArray());
        }

        public void Write(string subdirectory, string fileName)
        {
            KPMP_leavePatientsOut_analysis_readWriteOptions readWriteOptions = new KPMP_leavePatientsOut_analysis_readWriteOptions(subdirectory, fileName);
            ReadWriteClass.WriteData(this.Analysis, readWriteOptions);
        }

        public void Read(string subdirectory, string fileName)
        {
            KPMP_leavePatientsOut_analysis_readWriteOptions readWriteOptions = new KPMP_leavePatientsOut_analysis_readWriteOptions(subdirectory, fileName);
            this.Analysis = ReadWriteClass.ReadRawData_and_FillArray<KPMP_leavePatientsOut_analysis_line_class>(readWriteOptions);
        }

    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class KPMP_subsegmentalRNASeq_rawCounts_line_class : IFill_de
    {
        public string Sample_name { get; set; }
        public float Expression_value { get; set; }
        public string Symbol { get; set; }
        public string Patient { get; set; }
        public string Region { get; set; }
        public string Dataset { get; set; }
        public KPMP_value_type_enum Value_type { get; set; }

        public string[] Symbols_for_de { get { return new string[] { Symbol }; } }
        public string[] Names_for_de { get { return new string[] { Region, Patient, Value_type.ToString() }; } }
        public double Value_for_de { get { return Expression_value; } }

        public KPMP_subsegmentalRNASeq_rawCounts_line_class()
        {
            this.Sample_name = "";
            this.Symbol = "";
            this.Patient = "";
            this.Region = "";
            this.Dataset = "";
        }

        public KPMP_subsegmentalRNASeq_rawCounts_line_class Deep_copy()
        {
            KPMP_subsegmentalRNASeq_rawCounts_line_class copy = (KPMP_subsegmentalRNASeq_rawCounts_line_class)this.MemberwiseClone();
            copy.Sample_name = (string)this.Sample_name.Clone();
            copy.Symbol = (string)this.Symbol.Clone();
            copy.Patient = (string)this.Patient.Clone();
            copy.Region = (string)this.Region.Clone();
            copy.Dataset = (string)this.Dataset.Clone();
            return copy;
        }
    }

    class KPMP_subsegmentalRNASeq_rawCounts_class
    {
        public KPMP_subsegmentalRNASeq_rawCounts_line_class[] RawCounts { get; set; }
        public KPMP_integration_paper_metadata_class Dataset_patient { get; set; }
        public string[] Bg_genes { get; set; }
        
        public KPMP_subsegmentalRNASeq_rawCounts_class()
        {
            Dataset_patient = new KPMP_integration_paper_metadata_class();
            RawCounts = new KPMP_subsegmentalRNASeq_rawCounts_line_class[0];
            this.Bg_genes = new string[0];
        }

        private void Add_to_array(KPMP_subsegmentalRNASeq_rawCounts_line_class[] add_rawCounts)
        {
            int this_length = RawCounts.Length;
            int add_length = add_rawCounts.Length;
            int new_length = this_length + add_length;
            KPMP_subsegmentalRNASeq_rawCounts_line_class[] new_rawCounts = new KPMP_subsegmentalRNASeq_rawCounts_line_class[new_length];
            int indexNew = -1;
            for (int indexThis = 0; indexThis < this_length; indexThis++)
            {
                indexNew++;
                new_rawCounts[indexNew] = this.RawCounts[indexThis];
            }
            for (int indexAdd = 0; indexAdd < add_length; indexAdd++)
            {
                indexNew++;
                new_rawCounts[indexNew] = add_rawCounts[indexAdd];
            }
            this.RawCounts = new_rawCounts;
        }

        private void Set_patients_regions()
        {
            int rawCounts_length = this.RawCounts.Length;
            KPMP_subsegmentalRNASeq_rawCounts_line_class rawCounts_line;
            string[] splitStrings;
            for (int indexRC = 0; indexRC < rawCounts_length; indexRC++)
            {
                rawCounts_line = this.RawCounts[indexRC];
                splitStrings = rawCounts_line.Sample_name.Split('_');
                rawCounts_line.Patient = (string)splitStrings[1].Clone();
                rawCounts_line.Region = (string)splitStrings[0].Clone();
            }
        }

        private void Generate_dataset_patient_instance()
        {
            KPMP_integration_paper_metadata_line_class dataset_patient_line;
            List<KPMP_integration_paper_metadata_line_class> dataset_patient_list = new List<KPMP_integration_paper_metadata_line_class>();
            int data_length = this.RawCounts.Length;
            this.RawCounts = this.RawCounts.OrderBy(l => l.Dataset).ThenBy(l => l.Patient).ToArray();
            KPMP_subsegmentalRNASeq_rawCounts_line_class lmd_line;
            for (int indexRawCounts = 0; indexRawCounts < data_length; indexRawCounts++)
            {
                lmd_line = this.RawCounts[indexRawCounts];
                if (  (indexRawCounts == 0)
                    || (!lmd_line.Dataset.Equals(this.RawCounts[indexRawCounts-1].Dataset))
                    || (!lmd_line.Patient.Equals(this.RawCounts[indexRawCounts-1].Patient)))
                {
                    dataset_patient_line = new KPMP_integration_paper_metadata_line_class();
                    dataset_patient_line.Dataset = (string)lmd_line.Dataset.Clone();
                    dataset_patient_line.Libraries = new string[] { (string)lmd_line.Patient.Clone() };
                    dataset_patient_list.Add(dataset_patient_line);
                }
            }
            Dataset_patient.Add_to_array(dataset_patient_list.ToArray());
        }

        private void Generate_bg_genes()
        {
            List<string> bg_genes_list = new List<string>();
            this.RawCounts = this.RawCounts.OrderBy(l => l.Symbol).ToArray();
            int rawCounts_length = this.RawCounts.Length;
            KPMP_subsegmentalRNASeq_rawCounts_line_class rawCounts_line;
            for (int indexRawCounts=0; indexRawCounts<rawCounts_length; indexRawCounts++)
            {
                rawCounts_line = this.RawCounts[indexRawCounts];
                if ((indexRawCounts==0)||(!rawCounts_line.Symbol.Equals(this.RawCounts[indexRawCounts-1].Symbol)))
                {
                    bg_genes_list.Add(rawCounts_line.Symbol);
                }
            }
            this.Bg_genes = bg_genes_list.ToArray();
        }

        public void Generate()
        {
            Read_raw_counts_and_set_as_array();
            Generate_bg_genes();
            Set_patients_regions();
            Generate_dataset_patient_instance();
            //Calculate_subsegmental_averages();
            //Calculate_ratio_and_log2ratio_glomerula_vs_proxTub_and_vice_verse_for_single_values();
            //Keep_only_indicated_subsegments("Glom", "ProxTub", "Ratio Glom vs ProxTub", "Ratio ProxTub vs Glom", "Log2Ratio Glom vs ProxTub", "Log2Ratio ProxTub vs Glom");
            //Add_lines_with_zscore_normalized_expression_of_each_gene();
        }

        public KPMP_integration_paper_metadata_class Get_deep_copy_of_dataset_patient_instance()
        {
            return Dataset_patient.Deep_copy();
        }

        public string[] Get_deep_copy_of_bg_genes()
        {
            int bg_genes_length = this.Bg_genes.Length;
            string[] deep_copy = new string[bg_genes_length];
            for (int indexDC=0; indexDC<bg_genes_length; indexDC++)
            {
                deep_copy[indexDC] = (string)this.Bg_genes[indexDC].Clone();
            }
            return deep_copy;
        }

        public KPMP_standardized_dataset_class Generate_standardized_dataset_instance_with_all_values_filled_in_value_1st()
        {
            KPMP_standardized_dataset_line_class new_standardized_data_line;
            List<KPMP_standardized_dataset_line_class> standardized_data_list = new List<KPMP_standardized_dataset_line_class>();
            int raw_length = this.RawCounts.Length;
            KPMP_subsegmentalRNASeq_rawCounts_line_class subsegmental_rawCounts_line;
            RawCounts = RawCounts.OrderBy(l => l.Symbol).ThenBy(l => l.Region).ThenBy(l => l.Patient).ThenBy(l => l.Value_type).ToArray();
            for (int indexR = 0; indexR < raw_length; indexR++)
            {
                subsegmental_rawCounts_line = this.RawCounts[indexR];
                if ((indexR != 0)
                    && (subsegmental_rawCounts_line.Region.Equals(this.RawCounts[indexR - 1].Region))
                    && (subsegmental_rawCounts_line.Patient.Equals(this.RawCounts[indexR - 1].Patient))
                    && (subsegmental_rawCounts_line.Symbol.Equals(this.RawCounts[indexR - 1].Symbol))
                    && (subsegmental_rawCounts_line.Value_type.Equals(this.RawCounts[indexR - 1].Value_type)))
                {
                    throw new Exception();
                }
                new_standardized_data_line = new KPMP_standardized_dataset_line_class();
                new_standardized_data_line.Cell_segment = (string)subsegmental_rawCounts_line.Region.Clone();
                new_standardized_data_line.Dataset = (string)subsegmental_rawCounts_line.Dataset.Clone();
                new_standardized_data_line.Gene_symbol = (string)subsegmental_rawCounts_line.Symbol.Clone();
                new_standardized_data_line.PatientId = (string)subsegmental_rawCounts_line.Patient.Clone();
                new_standardized_data_line.Value_1st = subsegmental_rawCounts_line.Expression_value;
                new_standardized_data_line.Value_2nd = 0;
                new_standardized_data_line.Value_type_1st = subsegmental_rawCounts_line.Value_type;
                new_standardized_data_line.Value_type_2nd = KPMP_value_type_enum.No_selection;
                standardized_data_list.Add(new_standardized_data_line);
            }
            string dataset = KPMP_dataset_name_class.Lmd_rnaseq_iuosu;
            Dictionary<string, string[]> dataset_bgSymbol_dict = new Dictionary<string, string[]>();
            dataset_bgSymbol_dict.Add(dataset, Get_deep_copy_of_bg_genes());
            KPMP_standardized_dataset_class standard = new KPMP_standardized_dataset_class();
            standard.Add_to_existing_instances(standardized_data_list.ToArray(), dataset_bgSymbol_dict);
            return standard;
        }

        private void Read_raw_counts_and_set_as_array()
        {
            char delimiter = Global_class.Comma;
            string complete_directory = Global_directory_class.Experimental_data_directory + "Integration_paper_SubSegmentalRNASeq_OSUIU_2019-05-21\\Data\\";
            string fileName = "Subsegment_RNAseq_Raw_Counts_QN_add2_5-21-19.csv";
            string complete_fileName = complete_directory + fileName;
            string inputLine;
            StreamReader reader = new StreamReader(complete_fileName);
            string headline = reader.ReadLine();
            string[] columnNames = headline.Split(delimiter);
            string columnName;
            string columnEntry;
            string[] columnEntries;
            int columnEntries_length;
            int columnNames_length = columnNames.Length;
            string current_geneId;
            KPMP_subsegmentalRNASeq_rawCounts_line_class subsegmental_rawCounts_line;
            List<KPMP_subsegmentalRNASeq_rawCounts_line_class> subsegmental_rawCounts_list = new List<KPMP_subsegmentalRNASeq_rawCounts_line_class>();
            while ((inputLine = reader.ReadLine()) != null)
            {
                columnEntries = inputLine.Split(delimiter);
                columnEntries_length = columnEntries.Length;
                current_geneId = columnEntries[0];
                if (columnEntries_length != columnNames_length) { throw new Exception(); }
                for (int indexC = 1; indexC < columnEntries_length; indexC++)
                {
                    columnEntry = columnEntries[indexC];
                    columnName = columnNames[indexC];
                    subsegmental_rawCounts_line = new KPMP_subsegmentalRNASeq_rawCounts_line_class();
                    subsegmental_rawCounts_line.Expression_value = float.Parse(columnEntry);
                    subsegmental_rawCounts_line.Symbol = (string)current_geneId.Clone();
                    subsegmental_rawCounts_line.Sample_name = (string)columnName.Clone();
                    subsegmental_rawCounts_line.Dataset = KPMP_dataset_name_class.Lmd_rnaseq_iuosu;
                    if (columnName.IndexOf("RatioAVG ") == 0)
                    {
                        subsegmental_rawCounts_line.Value_type = KPMP_value_type_enum.Ratioavg;
                    }
                    else if (columnName.IndexOf("P ") == 0)
                    {
                        subsegmental_rawCounts_line.Value_type = KPMP_value_type_enum.Minus_log10_pvalue;
                        subsegmental_rawCounts_line.Expression_value = -(float)Math.Log10(subsegmental_rawCounts_line.Expression_value);
                    }
                    else
                    {
                        subsegmental_rawCounts_line.Value_type = KPMP_value_type_enum.Single_value;
                    }
                    subsegmental_rawCounts_list.Add(subsegmental_rawCounts_line);
                }
            }
            this.RawCounts = subsegmental_rawCounts_list.ToArray();
        }

        public KPMP_subsegmentalRNASeq_rawCounts_class Deep_copy()
        {
            KPMP_subsegmentalRNASeq_rawCounts_class copy = (KPMP_subsegmentalRNASeq_rawCounts_class)this.MemberwiseClone();
            int rawCounts_length = this.RawCounts.Length;
            copy.RawCounts = new KPMP_subsegmentalRNASeq_rawCounts_line_class[rawCounts_length];
            for (int indexRaw = 0; indexRaw < rawCounts_length; indexRaw++)
            {
                copy.RawCounts[indexRaw] = this.RawCounts[indexRaw].Deep_copy();
            }
            return copy;
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class KPMP_singleNucleusCell_averageExpression_line_class
    {
        public string Symbol { get; set; }
        public string Cell_type { get; set; }
        public string Patient { get; set; }
        public KPMP_value_type_enum Value_type { get; set; }
        public float Expression_value { get; set; }
        public string Dataset { get; set; }

        public KPMP_singleNucleusCell_averageExpression_line_class Deep_copy()
        {
            KPMP_singleNucleusCell_averageExpression_line_class copy = (KPMP_singleNucleusCell_averageExpression_line_class)this.MemberwiseClone();
            copy.Symbol = (string)this.Symbol.Clone();
            copy.Cell_type = (string)this.Cell_type.Clone();
            copy.Patient = (string)this.Patient.Clone();
            copy.Dataset = (string)this.Dataset.Clone();
            return copy;
        }
    }

    class KPMP_singleNucleusCell_averageExpression_readWriteOptions_class : ReadWriteOptions_base
    {
        public KPMP_singleNucleusCell_averageExpression_readWriteOptions_class(string directory, string fileName)
        {
            this.File = directory + fileName;
            this.Key_propertyNames = new string[] { "Symbol", "Cell_type", "Patient", "Dataset", "Expression_value" };
            this.Key_columnNames = this.Key_propertyNames;
            this.HeadlineDelimiters = new char[] { Global_class.Tab };
            this.LineDelimiters = new char[] { Global_class.Tab };
            this.File_has_headline = true;
            this.Report = ReadWrite_report_enum.Report_nothing;
        }
    }

    class KPMP_singleNucleusCell_averageExpression_class
    {
        public KPMP_singleNucleusCell_averageExpression_line_class[] Data { get; set; }
        public KPMP_integration_paper_metadata_class Dataset_patient { get; set; }
        public Dictionary<string,string[]> Center_bgGenesInUpperCase_dict { get; set; }

        public KPMP_singleNucleusCell_averageExpression_class()
        {
            this.Data = new KPMP_singleNucleusCell_averageExpression_line_class[0];
            this.Dataset_patient = new KPMP_integration_paper_metadata_class();
            Center_bgGenesInUpperCase_dict = new Dictionary<string, string[]>();
        }

        private void Add_to_array(KPMP_singleNucleusCell_averageExpression_line_class[] add_data)
        {
            int this_length = this.Data.Length;
            int add_length = add_data.Length;
            int new_length = this_length + add_length;
            KPMP_singleNucleusCell_averageExpression_line_class[] new_data = new KPMP_singleNucleusCell_averageExpression_line_class[new_length];
            int indexNew = -1;
            for (int indexThis = 0; indexThis < this_length; indexThis++)
            {
                indexNew++;
                new_data[indexNew] = this.Data[indexThis];
            }
            for (int indexAdd = 0; indexAdd < add_length; indexAdd++)
            {
                indexNew++;
                new_data[indexNew] = add_data[indexAdd];
            }
            this.Data = new_data;
        }

        public void Replace_all_cellType0_by_cellType1(string old_cellType0, string new_cellType1)
        {
            foreach (KPMP_singleNucleusCell_averageExpression_line_class singleNucleusCell_line in this.Data)
            {
                if (singleNucleusCell_line.Cell_type.IndexOf(old_cellType0) == 0)
                {
                    singleNucleusCell_line.Cell_type = singleNucleusCell_line.Cell_type.Replace(old_cellType0, new_cellType1);
                }
            }
        }

        private void Add_center_specific_bg_genesinUpperCase_if_not_already_added(string center)
        {
            if (!Center_bgGenesInUpperCase_dict.ContainsKey(center))
            {
                string bgGenes_completeFileName = Global_kpmp_class.Get_bgGenesInUpperCase_completeFileName("SingleCellNucleus_bgGenes\\", center);
                string[] bgGenes_in_upperCase = ReadWriteClass.Read_string_array(bgGenes_completeFileName);
                Center_bgGenesInUpperCase_dict.Add(center, bgGenes_in_upperCase);
            }
        }

        public void Merge_proxTub_cluster_by_keeping_gene_with_highest_expression()
        {
            int data_length = this.Data.Length;
            KPMP_singleNucleusCell_averageExpression_line_class data_line;
            List<KPMP_singleNucleusCell_averageExpression_line_class> keep_data_list = new List<KPMP_singleNucleusCell_averageExpression_line_class>();
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                data_line = this.Data[indexData];
                if (data_line.Cell_type.IndexOf("PT") == 0)
                {
                    data_line.Cell_type = "PT";
                }
            }

            this.Data = this.Data.OrderBy(l => l.Dataset).ThenBy(l => l.Patient).ThenBy(l => l.Symbol).ThenBy(l => l.Value_type).ThenBy(l => l.Cell_type).ThenByDescending(l => l.Expression_value).ToArray();

            for (int indexData = 0; indexData < data_length; indexData++)
            {
                data_line = this.Data[indexData];
                if (data_line.Cell_type.IndexOf("PT") == 0)
                {
                    if ((indexData == 0)
                        || (!data_line.Dataset.Equals(this.Data[indexData - 1].Dataset))
                        || (!data_line.Patient.Equals(this.Data[indexData - 1].Patient))
                        || (!data_line.Symbol.Equals(this.Data[indexData - 1].Symbol))
                        || (!data_line.Cell_type.Equals(this.Data[indexData - 1].Cell_type))
                        || (!data_line.Value_type.Equals(this.Data[indexData - 1].Value_type)))
                    {
                        keep_data_list.Add(data_line);
                    }
                }
                else
                {
                    keep_data_list.Add(data_line);
                }
            }
            this.Data = keep_data_list.ToArray();
        }

        private void Generate_dataset_patient_instance()
        {
            this.Data = this.Data.OrderBy(l => l.Dataset).ThenBy(l => l.Patient).ToArray();
            int data_length = this.Data.Length;
            KPMP_singleNucleusCell_averageExpression_line_class averageExpression_line;
            KPMP_integration_paper_metadata_line_class dataset_patient_line;
            List<KPMP_integration_paper_metadata_line_class> dataset_patient_list = new List<KPMP_integration_paper_metadata_line_class>();
            for (int indexData=0; indexData<data_length;indexData++)
            {
                averageExpression_line = Data[indexData];
                if (  (indexData==0)
                    || (!averageExpression_line.Dataset.Equals(this.Data[indexData - 1].Dataset))
                    || (!averageExpression_line.Patient.Equals(this.Data[indexData - 1].Patient)))
                {
                    dataset_patient_line = new KPMP_integration_paper_metadata_line_class();
                    dataset_patient_line.Dataset = (string)averageExpression_line.Dataset.Clone();
                    dataset_patient_line.Libraries = new string[] { (string)averageExpression_line.Patient.Clone() };
                    dataset_patient_list.Add(dataset_patient_line);
                }
            }
            this.Dataset_patient.Add_to_array(dataset_patient_list.ToArray());
        }

        private void Set_all_value_types_to_single_value()
        {
            int data_length = this.Data.Length;
            KPMP_singleNucleusCell_averageExpression_line_class averageExpression_line;
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                averageExpression_line = Data[indexData];
                averageExpression_line.Value_type = KPMP_value_type_enum.Single_value;
            }
        }

        public void Generate()
        {
            Read_averageExpression_and_center_specific_bgGenesInUpperCase();
            Set_all_value_types_to_single_value();
            Check_if_datasetNames_match_names_defined_in_KPMP_dataset_name_class();
            Generate_dataset_patient_instance();
        }

        public KPMP_integration_paper_metadata_class Get_deep_copy_of_dataset_patient()
        {
            return Dataset_patient.Deep_copy();
        }

        public KPMP_standardized_dataset_class Generate_standardized_dataset_with_all_valueTypes_filled_in_value_1st()
        {
            int data_length = this.Data.Length;
            KPMP_standardized_dataset_line_class new_standardized_line;
            KPMP_standardized_dataset_line_class new_standardized_line_base = new KPMP_standardized_dataset_line_class();
            List<KPMP_standardized_dataset_line_class> new_standardized_list = new List<KPMP_standardized_dataset_line_class>();

            this.Data = this.Data.OrderBy(l => l.Dataset).ThenBy(l => l.Symbol).ThenBy(l => l.Patient).ThenBy(l => l.Cell_type).ThenBy(l => l.Value_type).ToArray();

            KPMP_singleNucleusCell_averageExpression_line_class singleNucleusCell_line;
            Dictionary<string, string[]> dataset_bgGeneProtein_dict = new Dictionary<string, string[]>();
            List<string> current_dataset_bg_genes = new List<string>();
            for (int indexData = 0; indexData < data_length; indexData++)
            {
                singleNucleusCell_line = this.Data[indexData];
                if ((indexData == 0)
                    || (!singleNucleusCell_line.Dataset.Equals(this.Data[indexData - 1].Dataset)))
                {
                    current_dataset_bg_genes.Clear();
                }
                if ((indexData == 0)
                    || (!singleNucleusCell_line.Dataset.Equals(this.Data[indexData - 1].Dataset))
                    || (!singleNucleusCell_line.Symbol.Equals(this.Data[indexData - 1].Symbol)))
                {
                    current_dataset_bg_genes.Add((string)singleNucleusCell_line.Symbol.Clone());
                }
                if ((indexData == data_length - 1)
                    || (!singleNucleusCell_line.Dataset.Equals(this.Data[indexData + 1].Dataset)))
                {
                    dataset_bgGeneProtein_dict.Add((string)singleNucleusCell_line.Dataset.Clone(), current_dataset_bg_genes.ToArray());
                }
                if ((indexData == 0)
                    || (!singleNucleusCell_line.Dataset.Equals(this.Data[indexData - 1].Dataset))
                    || (!singleNucleusCell_line.Patient.Equals(this.Data[indexData - 1].Patient))
                    || (!singleNucleusCell_line.Cell_type.Equals(this.Data[indexData - 1].Cell_type))
                    || (!singleNucleusCell_line.Symbol.Equals(this.Data[indexData - 1].Symbol)))
                {
                    new_standardized_line_base = new KPMP_standardized_dataset_line_class();
                    new_standardized_line_base.Dataset = (string)singleNucleusCell_line.Dataset.Clone();
                    new_standardized_line_base.PatientId = (string)singleNucleusCell_line.Patient.Clone();
                    new_standardized_line_base.Cell_segment = (string)singleNucleusCell_line.Cell_type.Clone();
                    new_standardized_line_base.Gene_symbol = (string)singleNucleusCell_line.Symbol.Clone();
                }
                if ((indexData != 0)
                    && (singleNucleusCell_line.Dataset.Equals(this.Data[indexData - 1].Dataset))
                    && (singleNucleusCell_line.Patient.Equals(this.Data[indexData - 1].Patient))
                    && (singleNucleusCell_line.Cell_type.Equals(this.Data[indexData - 1].Cell_type))
                    && (singleNucleusCell_line.Symbol.Equals(this.Data[indexData - 1].Symbol))
                    && (singleNucleusCell_line.Value_type.Equals(this.Data[indexData - 1].Value_type)))
                {
                    throw new Exception();
                }
                new_standardized_line = new_standardized_line_base.Deep_copy();
                new_standardized_line.Value_1st = singleNucleusCell_line.Expression_value;
                new_standardized_line.Value_type_1st = singleNucleusCell_line.Value_type;
                new_standardized_line.Value_2nd = 0;
                new_standardized_line.Value_type_2nd = KPMP_value_type_enum.No_selection;
                new_standardized_list.Add(new_standardized_line);
            }
            KPMP_standardized_dataset_class standardized_data = new KPMP_standardized_dataset_class();
            standardized_data.Add_to_existing_instances(new_standardized_list.ToArray(), dataset_bgGeneProtein_dict);
            return standardized_data;
        }

        private void Check_if_datasetNames_match_names_defined_in_KPMP_dataset_name_class()
        {
            int data_length = this.Data.Length;
            KPMP_singleNucleusCell_averageExpression_line_class singleNucleusCell_average_line;
            for (int indexD=0; indexD<data_length;indexD++)
            {
                singleNucleusCell_average_line = this.Data[indexD];
                if (   (!singleNucleusCell_average_line.Dataset.Equals(KPMP_dataset_name_class.SingleCell_premiere))
                    && (!singleNucleusCell_average_line.Dataset.Equals(KPMP_dataset_name_class.SingleCell_ucsf))
                    && (!singleNucleusCell_average_line.Dataset.Equals(KPMP_dataset_name_class.SingleNucleus_ucsd)))
                {
                    throw new Exception();
                }
            }
        }

        private void Read_averageExpression_and_center_specific_bgGenesInUpperCase()
        {
            string directory = Global_directory_class.Seurat_integratedCluster_avgExpression_directory;
            string fileName = "AverageRNAExpression_RNA_counts.txt";
            KPMP_singleNucleusCell_averageExpression_readWriteOptions_class readWriteOptions = new KPMP_singleNucleusCell_averageExpression_readWriteOptions_class(directory,fileName);
            this.Data = ReadWriteClass.ReadRawData_and_FillArray< KPMP_singleNucleusCell_averageExpression_line_class>(readWriteOptions);
            List<KPMP_singleNucleusCell_averageExpression_line_class> keep = new List<KPMP_singleNucleusCell_averageExpression_line_class>();
            foreach (KPMP_singleNucleusCell_averageExpression_line_class data_line in this.Data)
            {
                if (data_line.Expression_value<0) { throw new Exception(); }
                data_line.Value_type = KPMP_value_type_enum.Single_value;
                if (data_line.Dataset.Equals("SC RNASeq I")) { data_line.Dataset = (string)KPMP_dataset_name_class.SingleCell_premiere.Clone(); }
                if (data_line.Expression_value>0)
                {
                    keep.Add(data_line);
                }
            }
            this.Data = keep.ToArray();
            Add_center_specific_bg_genesinUpperCase_if_not_already_added(KPMP_dataset_name_class.SingleCell_premiere);
            Add_center_specific_bg_genesinUpperCase_if_not_already_added(KPMP_dataset_name_class.SingleNucleus_ucsd);
            Add_center_specific_bg_genesinUpperCase_if_not_already_added(KPMP_dataset_name_class.SingleCell_ucsf);
        }

        public KPMP_singleNucleusCell_averageExpression_class Deep_copy()
        {
            KPMP_singleNucleusCell_averageExpression_class copy = (KPMP_singleNucleusCell_averageExpression_class)this.MemberwiseClone();
            int data_length = this.Data.Length;
            copy.Data = new KPMP_singleNucleusCell_averageExpression_line_class[data_length];
            for (int indexD = 0; indexD < data_length; indexD++)
            {
                copy.Data[indexD] = this.Data[indexD].Deep_copy();
            }
            return copy;
        }

    }


}
