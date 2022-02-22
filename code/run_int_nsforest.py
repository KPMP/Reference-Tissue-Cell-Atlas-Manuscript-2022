#Find minimal marker sets using NS_Forest for integrated scRNAseq/snRNAseq data

from NSForest_v3 import *
import itertools
import numpy as np
import pandas as pd
import scanpy as sc
import re

adata = sc.read_h5ad("tis.integrated.h5ad")
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata.raw = adata
sc.pp.scale(adata, max_value=10)
adata.obs['louvain'] = adata.obs['celltype'].to_numpy().astype(np.str)
adata.obs['louvain'] = pd.Series(adata.obs['louvain'], dtype="category")
mydata_markers = NS_Forest(adata)

#0: tPC-IC
#1: POD
#2: PT
#3: DTL
#4: PT/PEC
#5: ATL/TAL
#6: TAL
#7: DCT
#8: CNT
#9: PC
#10: IC
#11: EC-AEA
#12: EC-PTC
#13: EC-GC
#14: VSMC/P
#15: MC
#16: FIB
#17: T
#18: NKT/NKC
#19: MYL
#20: MAC
#21: B
