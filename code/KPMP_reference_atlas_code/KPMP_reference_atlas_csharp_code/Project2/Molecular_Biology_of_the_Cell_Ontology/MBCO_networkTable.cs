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

namespace MBCO_network_and_nodes
{
    enum NWedge_type_enum { E_m_p_t_y, Arrow, Thick_dotted_line, Dotted_line, Dashed_line, Line };

    class NetworkTable_line_class
    {
        public string Source { get; set; }
        public string Target { get; set; }
        public string Edge_label { get; set; }
        public float Width { get; set; }
        public NWedge_type_enum Edge_type { get; set; }

        public NetworkTable_line_class()
        {
            Width = 1;
            Edge_label = "";
            Edge_type = NWedge_type_enum.Arrow;
        }

        public NetworkTable_line_class Deep_copy()
        {
            NetworkTable_line_class copy = (NetworkTable_line_class)this.MemberwiseClone();
            copy.Source = (string)this.Source.Clone();
            copy.Target = (string)this.Target.Clone();
            copy.Edge_label = (string)this.Edge_label.Clone();
            return copy;
        }
    }
}
