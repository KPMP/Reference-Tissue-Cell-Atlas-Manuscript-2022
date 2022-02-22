#Find minimal marker sets using NS_Forest for sn data

from NSForest_v3 import *
import itertools
import numpy as np
import pandas as pd
import scanpy as sc
import re

mydata=sc.read("GSE121862_UCSD-WU_Single_Nuclei_Cluster_Annotated_Raw_UMI_Matrix.tsv", delimiter='\t',cache=True)
mydata = mydata.transpose()
#mydata.var['mt'] = mydata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
#sc.pp.calculate_qc_metrics(mydata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
sc.pp.normalize_total(mydata, target_sum=1e4)
sc.pp.log1p(mydata)
sc.pp.highly_variable_genes(mydata, min_mean=0.0125, max_mean=3, min_disp=0.5)
mydata.raw = mydata
#sc.pp.regress_out(mydata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(mydata, max_value=10)

clusts = [re.sub("_.*", "", elem) for elem in mydata.obs.index]
clusts = [re.sub("C", "", elem) for elem in clusts]

mapping = pd.read_csv("GSE121862_UCSD-WU_Single_Nuclei_Cluster_Annotations.csv", delimiter=",")

for i in range(0,len(clusts)):
    clusts[i] = mapping['Abbn'][int(clusts[i]) - 1]

clusts = [re.sub("PT-.*", "PT",elem) for elem in clusts]
clusts = [re.sub("ATL-.*", "ATL",elem) for elem in clusts]
clusts = [re.sub("TAL-.*", "TAL",elem) for elem in clusts]
clusts = [re.sub("PC-.*", "PC",elem) for elem in clusts]
clusts = [re.sub("IC-.*", "IC",elem) for elem in clusts]
clusts = [re.sub("EC-.*", "EC",elem) for elem in clusts]


mydata.obs['louvain'] = np.array(clusts)
mydata.obs['louvain'] = pd.Series(mydata.obs['louvain'], dtype="category")

mydata = mydata[~mydata.obs['louvain'].isin(['IMM']),:]

mydata_markers = NS_Forest(mydata) #Runs NS_Forest on scanpy object
