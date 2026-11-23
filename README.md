# ComFlow
ComFlow: an optimal transport framework for deciphering intra-
and intercellular signaling landscapes in 2D and 3D tissue systems

## Installation
Install from PyPI by ```python pip install comflow```or install from source by cd to this directory and pip install .

[Optional] Install tradeSeq in R to analyze the CCC differentially expressed genes.
Currently, tradeSeq version 1.0.1 with R version 3.6.3 has been tested to work.
In order for the R-python interface to work properly, rpy2==3.4.2 and anndata2ri==1.0.6 should be installed. To use this feature, install from PyPI by pip install comflow[tradeSeq] or from source by pip install .[tradeSeq].

## Usage
### Basic usage

Import packages

import comflow as cf
import scanpy as sc
import pandas as pd
import numpy as np
Load a spatial dataset
(e.g., a Visium dataset)

adata = sc.datasets.visium_sge(sample_id='V1_Mouse_Brain_Sagittal_Posterior')
adata.var_names_make_unique()
Basic processing

sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
Specify ligand-receptor pairs

LR=np.array([['Tgfb1', 'Tgfbr1_Tgfbr2', 'Tgfb_pathway'],['Tgfb2', 'Tgfbr1_Tgfbr2', 'Tgfb_pathway'],['Tgfb3', 'Tgfbr1_Tgfbr2', 'Tgfb_pathway']],dtype=str)
df_ligrec = pd.DataFrame(data=LR)
(or use pairs from a ligand-receptor database df_ligrec=ct.pp.ligand_receptor_database(database='CellChat', species='mouse').)

Construct CCC networks
Use collective optimal transport to construct CCC networks for the ligand-receptor pairs with a spatial distance constraint of 200 (coupling between cells with distance greater than 200 is prohibited). For example, the spot-by-spot matrix for the pair Tgfb1 (ligand) and Tgfbr1_Tgfbr2 (receptor)is stored in adata.obsp['comflow-user_database-Tgfb1-Tgfbr1_Tgfbr2']. The total sent or received signal for each pair is stored in adata.obsm['comflow-user_database-sum-sender'] and adata.obsm['comflow-user_database-sum-receiver'].

cf.tl.spatial_communication(adata,
    database_name=database_name, df_ligrec=df_filtered, pathways=pathways, co_exp_ratio=co_exp_ratio, dis_thr=dis_thr, heteromeric=True, pathway_sum=True, copy = True)
## Documentation
See detailed documentation and examples at [https://comflow.readthedocs.io/en/latest/index.html.](https://comflow.readthedocs.io/en/latest/index.html.)

## References
Xie, M. B., Zhu, Y., Hu, Y., Kim, M., Liu, Y. H., Huang, J., Li, B., Rafat, M., Meltzer,S., Hwang, T. H., Cang, Z., Zhao, Y., Zhou, X. M.. ComFlow: an optimal transport framework for deciphering intra- and intercellular signaling landscapes in 2D and 3D tissue systems. 
