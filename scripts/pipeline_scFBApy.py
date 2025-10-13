# -*- coding: utf-8 -*-
"""

@author: bruno.galuzzi
"""

#%% load libraries
import cobra as cb
import numpy as np
import scanpy as sc
from utils_scFBApy import scFBApy,repairNeg,popModel
#%% load cobra model
model = cb.io.read_sbml_model("models/model.xml")
#%% load single cell dataset
adata=sc.read_h5ad("datasets/BC03LNdataset")

#%%
adata,cont=repairNeg(adata,
              "BC03LN_Pooled",
              filter_bulk=True,
              epsilon=1e-4)   
#%% denoising with MAGIC agorithm
#sc.external.pp.magic(adata,knn=3)

#%%Compute optimal fluxes
adata_fluxes_pop=scFBApy(
            model,                                #metabolic model
            adata,                                #scRNA-seq data
            objective="Biomass",                  #objective function to maximize
            cooperation=True,                     #cooperation between cells or not?
            eps=0.001,
            compute_fva=True,                     #compute FVA to integrate RAS 
            npop_fva=5,                        
            type_ras_normalization="max",         #type of normalization (max or sum)   
            and_expression=np.nanmin,             #type of operation for AND expression
            or_expression=np.nansum,              #type of operation for OR expression
            fraction_of_optimum=0,                #fraction of optimum FVA
            processes=1,                          #processes for FVA
            round_c=10                            #digits of approximation
            )
#%%
print("Biomass per cell: ",adata_fluxes_pop.to_df()["Biomass"].mean())
