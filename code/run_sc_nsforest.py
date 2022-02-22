#Find minimal marker sets using NS_Forest for sc data

from NSForest_v3 import *
import itertools
import numpy as np
import pandas as pd
import scanpy as sc
import re

mydata=sc.read("Premiere_Raw.data_062220.txt", delimiter='\t',cache=True)
mydata = mydata.transpose()
#mydata.var['mt'] = mydata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
#sc.pp.calculate_qc_metrics(mydata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
sc.pp.normalize_total(mydata, target_sum=1e4)
sc.pp.log1p(mydata)
sc.pp.highly_variable_genes(mydata, min_mean=0.0125, max_mean=3, min_disp=0.5)
mydata.raw = mydata
#sc.pp.regress_out(mydata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(mydata, max_value=10)

barcodes = pd.read_csv("PREMIERE_TIS_JUNE2019_CELLBARCODES_CLUSTER.txt", delimiter='\t', header=None)
tmp = [re.sub("X", "", elem) for elem in mydata.obs.index]
tmp = [re.sub("\\.", "-", elem) for elem in tmp]
tmp = [re.sub("SamplePRE038-1", "SamplePRE038.1", elem) for elem in tmp]
tmp = [re.sub("SamplePRE027-1", "SamplePRE027.1", elem) for elem in tmp]
tmp = pd.DataFrame(tmp)
tmp2 = tmp.merge(barcodes, how="inner", on=0,)

tmp2[2] = [re.sub("PT-.*", "PT",elem) for elem in tmp2[2]]
tmp2[2] = [re.sub("IC-.*", "IC",elem) for elem in tmp2[2]]
tmp2[2] = [re.sub("tPC-IC", "IC",elem) for elem in tmp2[2]]
tmp2[2] = [re.sub("PC-CNT", "PC",elem) for elem in tmp2[2]]
tmp2[2] = [re.sub("EC-.*", "EC",elem) for elem in tmp2[2]]
tmp2[2] = [re.sub("T-.*", "Tcells",elem) for elem in tmp2[2]]


mydata.obs['louvain'] = tmp2[2].to_numpy()
mydata.obs['louvain'] = pd.Series(mydata.obs['louvain'], dtype="category")

mydata_markers = NS_Forest(mydata) #Runs NS_Forest on scanpy object


Markers = list(itertools.chain.from_iterable(mydata_markers['NSForest_Markers'])) #gets list of minimal markers from dataframe for display in scanpy plotting functions
Binary_Markers = list(itertools.chain.from_iterable(mydata_markers['Binary_Genes'])) #gets list of binary markers from dataframe for display in scanpy plotting functions

textfile = open("sc_markers.txt2", "w")
for element in Markers:
    textfile.write(element + "\n")
    
textfile2 = open("sc_binary_markers.txt2", "w")
for element in BinaryMarkers:
    textfile.write(element + "\n")
    
