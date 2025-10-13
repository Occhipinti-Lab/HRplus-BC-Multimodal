# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 11:55:29 2023

@author: bruno.galuzzi
"""
#%% Library
import pandas as pd
import re
from progressbar import ProgressBar,Bar,Percentage
from scanpy import AnnData
from cobra.flux_analysis.variability import find_essential_reactions,find_essential_genes
from cobra import  Model
import cobra as cb
import numpy as np
import copy
from cobra.io import read_sbml_model


#%% function to create a population of models
def popModel(model_orig,
             n_pop,                       #how many models to create
             objective=None,              #objective functions
             c_round=10,
             fraction_of_optimum=0,
             compute_fva=True,
             npop_fva=2,
             processes=1):

    
    model=model_orig.copy()

    #create a model
    
    model_merged=createPopModel(model,n_pop)

    # create objective function
    model_merged=newObjfun(model_merged,objective,n_pop)
    

    if compute_fva:
        
        #create a model

        model_merged2=createPopModel(model_orig.copy(),npop_fva)
        
        # create objective function
        model_merged2=newObjfun(model_merged2,objective,npop_fva)
        

        
        reactions=[el.id for el in model.reactions]
        reactions_ex=[el.id for el in model.exchanges]
        reaction_noex=[el for el in reactions if el not in reactions_ex]
        reaction_list=[el+"_cell0" for el in reaction_noex]
        reaction_list.extend([el+"_#" for el in reactions_ex])
        
        dfFVA=cb.flux_analysis.flux_variability_analysis(model_merged2,
                                                         fraction_of_optimum=fraction_of_optimum,
                                                         processes=processes,
                                                         reaction_list=reaction_list                           
                                                                                    ).round(c_round)
    
        for reaction in reaction_noex:

            for i in range(n_pop):
                  model_merged.reactions.get_by_id(reaction+"_cell"+str(i)).bounds=(dfFVA.loc[reaction+"_cell0","minimum"],
                                                                                   dfFVA.loc[reaction+"_cell0","maximum"])
        for reaction in reactions_ex:
                  model_merged.reactions.get_by_id(reaction+"_#").bounds=(dfFVA.loc[reaction+"_#","minimum"],
                                                                          dfFVA.loc[reaction+"_#","maximum"])                                                                 
    return model_merged


def createPopModel(model,n_pop):
    #
    # Create a population of networks
    #
    #
    exchanges=[el.id for el in model.exchanges]
    reactions_noexchange=[el for el in model.reactions if el.id not in exchanges]
    reactions_exchange=[el for el in model.reactions if el.id  in exchanges]
    
    #considero i metaboliti delle exchanges
    metabolite_exchange=[]
    for exchange in reactions_exchange:
        for metabolite in exchange.metabolites:
            metabolite_exchange.append(metabolite.id)
    metabolite_exchange=list(set(metabolite_exchange))
    
    for reaction in reactions_exchange:
        reaction.id=reaction.id+"_#"              #gli do un nome a parte
    
    #remove exchange reactions from the model
    model.remove_reactions(reactions_exchange)  
    

    model_merged=Model("model_merged")
    reactions_to_add=[]

    for i in range(n_pop):
        reactions2=copy.deepcopy(reactions_noexchange)
        for reaction in reactions2:
            reaction.id=reaction.id+"_cell"+str(i)
            for el in reaction.metabolites: #change also the metabolite names
                if "_cell" not in el.id and el.id not in metabolite_exchange:
                    el.id=el.id+"_cell"+str(i)
                
        reactions_to_add.extend(reactions2)
    
    #add reactions of the super-network
    model_merged.add_reactions(reactions_to_add)
    
    #add the exchange
    model_merged.add_reactions(reactions_exchange)

    return model_merged

def newObjfun(model,objective,n_pop):
    #come somma di tutte
    reactions_ids=[el.id for el in model.reactions]
    coefficients = dict()
    for reaction in model.reactions: 
        coefficients[model.reactions.get_by_id(reaction.id)] = 0
    if objective+"_#" in reactions_ids:
            coefficients[model.reactions.get_by_id(objective+"_#")] = 1
    else:
        for i in range(n_pop):
            coefficients[model.reactions.get_by_id(objective+"_cell"+str(i))] = 1                    
        
    model.objective=coefficients

    return model

def scPopFBA(model_orig,
          npop,
          ras_adata,
          compute_fva=True,
          type_ras_normalization="sum",
          eps=0,
          return_adata=True
          ):
    
    rasMatrix=ras_adata.to_df()   
    #single FBA con modello unico
    model=model_orig.copy()


    if type_ras_normalization=="sum":
        rasMatrix=rasMatrix.div(rasMatrix.sum())
    else:
        rasMatrix=rasMatrix.div(rasMatrix.max())
    rasMatrix=rasMatrix.fillna(0)

    for i in range(npop):
        for reaction in rasMatrix.columns:
            bounds_original=model.reactions.get_by_id(reaction+"_cell"+str(i)).bounds
            valRas=rasMatrix.iloc[i].loc[reaction]

            #rimappo i ras
            bounds=(valRas*bounds_original[0]-eps,
                    valRas*bounds_original[1]+eps)
            
            model.reactions.get_by_id(reaction+"_cell"+str(i)).bounds=bounds
            
        
    if return_adata:
        obj_res=model.optimize()
        dfTot=df2matrix(obj_res.fluxes.round(10),npop)
        adata=AnnData(dfTot.T)
        return adata,model
    else:
        return dfTot,model


"""
Function to compute scFBA
"""
def scFBA(model_orig,
          objective,
          npop,
          ras_adata,
          type_ras_normalization="sum",
          compute_fva=False,
          eps=0,
          round_c=10,
          processes=1,
          fraction_of_optimum=0,
          return_adata=True
          ):
    
    #%% normalize ras values
    rasMatrix=ras_adata.to_df()   
    if type_ras_normalization=="sum":
        rasMatrix=rasMatrix.div(rasMatrix.sum())
    else:
        rasMatrix=rasMatrix.div(rasMatrix.max())
    rasMatrix=rasMatrix.fillna(0)    
    #%% Set up the bounds using FVA
    model_orig2=model_orig.copy()
    reactions=[el.id for el in model_orig2.reactions]
    
    if compute_fva:
        dfFVA=cb.flux_analysis.flux_variability_analysis(model_orig2,fraction_of_optimum=fraction_of_optimum,processes=processes).round(round_c)
        for reaction in reactions:
            model_orig2.reactions.get_by_id(reaction).bounds=(dfFVA.loc[reaction,"minimum"],dfFVA.loc[reaction,"maximum"])

    

    dfTot=pd.DataFrame(index=[reaction.id for reaction in model_orig.reactions],columns=["cell"+str(i) for i in range(npop)])


    for i in range(npop):

        model=model_orig2.copy()
        for reaction in rasMatrix.columns:
            bounds_original=model.reactions.get_by_id(reaction).bounds
            valRas=rasMatrix.iloc[i].loc[reaction]

            bounds=(valRas*bounds_original[0]-eps,
                    valRas*bounds_original[1]+eps)
            
            model.reactions.get_by_id(reaction).bounds=bounds
            
        dfTot.loc[:,"cell"+str(i)]=model.optimize().fluxes.round(round_c)
            
    
    if return_adata:
        adata=AnnData(dfTot.T)
        return adata
    else:
        return dfTot


def scFBApy(model_orig,
            adata,
            objective,
            cooperation=True,
            compute_fva=True,
            eps=0,
            npop_fva=2,
            type_ras_normalization="max",
            and_expression=np.nanmin,
            or_expression=np.nansum,
            fraction_of_optimum=0,
            processes=1,
            round_c=10
            ):
    
    
    npop=adata.shape[0]
    
    model = read_sbml_model("models/model.xml")   
    exchanges=[el.id for el in model.exchanges]

    #%% Compute RAS

    ras_object=RAS_computation(adata,model)
    ras_adata=ras_object.compute(drop_duplicates=False,
                                  drop_na_rows=True,
                                  and_expression=and_expression,
                                  or_expression=or_expression,
                                  add_essential_reactions=False,
                                  add_essential_genes=False
                                 )

    #%% Compute fluxes 
    if cooperation:
        #cooperation between cells
        model_x_pop=model.copy()
    
        #modify exchanges
        for reaction in exchanges: 
           bounds=model_x_pop.reactions.get_by_id(reaction).bounds
           bounds=(bounds[0]*npop,bounds[1]*npop)
    
           model_x_pop.reactions.get_by_id(reaction).bounds=bounds
        
        #create super-network
        model_merged=popModel(model_x_pop,
                              npop,
                              objective,
                              compute_fva=compute_fva,
                              npop_fva=npop_fva
                              )

        #compute fluxes
        adata_fluxes_pop,model_integrated=scPopFBA(model_merged,
                            npop,
                            ras_adata,
                            eps=eps,
                            type_ras_normalization=type_ras_normalization,
                            return_adata=True)
    else:
        #no cooperation between cells
        adata_fluxes_pop=scFBA(model,
                          objective,
                          npop,
                          ras_adata,
                          eps=eps,
                          compute_fva=compute_fva,
                          type_ras_normalization=type_ras_normalization,
                          return_adata=True,
                          fraction_of_optimum=fraction_of_optimum,
                          processes=processes,
                          round_c=round_c
                          )    
    
    return adata_fluxes_pop
    
   
"""
Class to compute the RAS values

"""

class RAS_computation:

    def __init__(self,adata,model):
                                                       
        self._logic_operators = ['and', 'or', '(', ')']
        self.val_nan = np.nan

        # Build the dictionary for the GPRs
        df_reactions = pd.DataFrame(index=[reaction.id for reaction in model.reactions])
        gene_rules=[reaction.gene_reaction_rule for reaction in model.reactions]   
        
        gene_rules=[el.replace("OR","or").replace("AND","and").replace("(","( ").replace(")"," )") for el in gene_rules]        
        df_reactions['rule'] = gene_rules
        df_reactions = df_reactions.reset_index()
        df_reactions = df_reactions.groupby('rule').agg(lambda x: sorted(list(x)))
        
        self.dict_rule_reactions = df_reactions.to_dict()['index']

        # build useful structures for RAS computation
        self.model = model
        self.count_adata = adata.copy()
        self.genes = self.count_adata.var.index.intersection([gene.id for gene in model.genes])
        
        #check if there is one gene at least 
        if len(self.genes)==0:
            if len(self.genes) == 0:
                raise ValueError(
                    "❌ ERROR: No gene of the count matrix is in the metabolic model!\n"
                    "🔍 Please check if the gene IDs in adata.var_names match the metabolic model (likely Ensembl IDs)."
                )
        
        self.cell_ids = list(self.count_adata.obs.index.values)
        self.count_df_filtered = self.count_adata.to_df().T.loc[self.genes]
 
    def compute(self,
                or_expression=np.nansum,    # type of operation to do in case of an or expression (max, sum, mean)
                and_expression=np.nanmin,   # type of operation to do in case of an and expression(min, sum)
                drop_na_rows=True,          # if True remove the nan rows of the ras  matrix
                drop_duplicates=False,      # if true, remove duplicates rows
                regexp=re.compile(r"\([a-zA-Z0-9-.:\s]+\)"),  # regular expression inside a parenthesis
                print_progressbar=True,     # if True, print the progress bar
                add_count_metadata=True,    # if True add metadata of cells in the ras adata
                add_met_metadata=True,      # if True add metadata from the metabolic model (gpr and compartments of reactions)
                add_essential_reactions=False,
                add_essential_genes=False
                ):

        self.or_function = or_expression
        self.and_function = and_expression
        
        ras_df = pd.DataFrame(index=range(len(self.dict_rule_reactions)), columns=self.cell_ids)
        ras_df[:][:] = self.val_nan
        
        if print_progressbar:
            pbar = ProgressBar(widgets=[Percentage(), Bar()], maxval=len(self.dict_rule_reactions)).start()
            i = 0
        
        # for loop on reactions
        ind = 0       
        for rule, reaction_ids in self.dict_rule_reactions.items():
            if len(rule) != 0:
                # there is one gene at least in the formula
                rule_split = rule.split()
                rule_split_elements = list(filter(lambda x: x not in self._logic_operators, rule_split))  # remove of all logical operators
                rule_split_elements = list(np.unique(rule_split_elements))                                # genes in formula
                
                # which genes are in the count matrix?                
                genes_in_count_matrix = list(set([el for el in rule_split_elements if el in self.genes]))
                genes_notin_count_matrix = list(set([el for el in rule_split_elements if el not in self.genes]))


                if len(genes_in_count_matrix) > 0: #there is at least one gene in the count matrix
                     if len(rule_split) == 1:
                         #one gene --> one reaction
                         ras_df.iloc[ind] = self.count_df_filtered.loc[genes_in_count_matrix]
                     else:                        
                        # more genes in the formula
                        lista = re.findall(regexp, rule)
                        if len(lista) == 0:
                             #or/and sequence
                             matrix = self.count_df_filtered.loc[genes_in_count_matrix].values
                             if len(genes_notin_count_matrix) > 0:
                                matrix = np.vstack([matrix, [self.val_nan for el in self.cell_ids]])

                             if 'or' in rule_split: 
                                ras_df.iloc[ind] = self.or_function(matrix, axis=0)
                             else:
                                ras_df.iloc[ind] = self.and_function(matrix, axis=0)
                        else:
                            # ho almeno una tonda
                            data = self.count_df_filtered.loc[genes_in_count_matrix]  # dataframe of genes in the GPRs
                            genes = data.index
                            j = 0
                             
                            for cellid in self.cell_ids:    #for loop on the cells
                                lista_cell = lista.copy()
                                rule_cell = rule
                                 
                                while len(lista_cell) > 0:
                                    #
                                    for el in lista_cell:
                                        #print(el[1:-1])
                                        value = self._evaluate_expression(el[1:-1].split(), data[cellid], genes)
                                        rule_cell = rule_cell.replace(el, str(value))   
                                    lista_cell = re.findall(regexp, rule_cell)      
         
                                ras_df.iloc[ind, j] = self._evaluate_expression(rule_cell.split(), data[cellid], genes)
                                j=j+1
      
            ind = ind+1
            #update percentage
            if print_progressbar:
                pbar.update(i+1)
                i = i+1
        
        if print_progressbar:
            pbar.finish()
        
        ras_df=ras_df.astype("float")    
        ras_df['REACTIONS'] = [reaction_ids for rule,reaction_ids in self.dict_rule_reactions.items()]
        
        reactions_common = pd.DataFrame()
        reactions_common["REACTIONS"] = ras_df['REACTIONS']
        reactions_common["proof2"] = ras_df['REACTIONS']
        reactions_common = reactions_common.explode('REACTIONS')
        reactions_common = reactions_common.set_index("REACTIONS")

        ras_df = ras_df.explode("REACTIONS")
        ras_df = ras_df.set_index("REACTIONS")

        if drop_na_rows:
            ras_df = ras_df.dropna(how="all")
            
        if drop_duplicates:
            ras_df = ras_df.drop_duplicates()
        
        #create AnnData structure for RAS
        ras_adata = AnnData(ras_df.T)

        #add metadata
        if add_count_metadata:
            ras_adata.var["common_gprs"] = reactions_common.loc[ras_df.index]
            ras_adata.var["common_gprs"] = ras_adata.var["common_gprs"].apply(lambda x: ",".join(x))
            for el in self.count_adata.obs.columns:
                ras_adata.obs["countmatrix_"+el]=self.count_adata.obs[el]

        if add_met_metadata:
            if len(self.model.compartments)>0:
                  ras_adata.var['compartments']=[list(self.model.reactions.get_by_id(reaction).compartments) for reaction in ras_adata.var.index]  
                  ras_adata.var['compartments']=ras_adata.var["compartments"].apply(lambda x: ",".join(x))
            
            ras_adata.var['GPR rule'] = [self.model.reactions.get_by_id(reaction).gene_reaction_rule for reaction in ras_adata.var.index]

        if add_essential_reactions:            
            essential_reactions=find_essential_reactions(self.model)
            essential_reactions=[el.id for el in essential_reactions]            
            ras_adata.var['essential reactions']=["yes" if el in essential_reactions else "no" for el in ras_adata.var.index]
        
        if add_essential_genes:
            essential_genes=find_essential_genes(self.model)
            essential_genes=[el.id for el in essential_genes]
            ras_adata.var['essential genes']=[" ".join([gene for gene in genes.split()  if gene in essential_genes]) for genes in ras_adata.var["GPR rule"]]
        
        return ras_adata


    def _check_number(self,value):
      try:
        float(value)
        return True
      except ValueError:
        return False

    def _evaluate_expression(self, rule_split, values_cellid, genes):
        
        #ci sono per forza solo or
        rule_split2 = list(filter(lambda x: x != "or" and x!="and", rule_split))   

        values = list()
        i=0
        for el in rule_split2:
             if self._check_number(el):
                 values.append(float(el))
             elif el in genes:
                 values.append(values_cellid[el])
             else:
                 values.append(self.val_nan)
                 i=i+1
                 
        if i==len(rule_split2):
            return self.val_nan
        if "or" in rule_split:
            #or sequence
            return self.or_function(values)
        else:
            #and sequence
            return self.and_function(values)
        
        
        
def df2matrix(df,npop):
    
    index=df.index
    index=[el.split("_cell")[0] for el in index]
    reactions=list(set(index))

    dfCells=pd.DataFrame(index=reactions,columns=["cell"+str(i) for i in range(npop)])


    for reaction in reactions:
        if "_#" in reaction:
            dfCells.loc[reaction,:]=df.loc[reaction]
        else:
            valori=[df.loc[reaction+"_cell"+str(i)] for i in range(npop)]
            dfCells.loc[reaction,:]=valori

    return dfCells


def repairNeg(adata,
              bulk_el,
              filter_bulk=True,
              epsilon=1e-4):
    
    adata2=adata.copy()
    
    
    genes=adata2.var.index
    cells=[cell for cell in adata2.obs.index if cell!=bulk_el]
    df=adata2.to_df().loc[cells,:]
    dfBulk=adata2.to_df().loc[bulk_el,:]
    adata2=adata2[cells,:]

    cont=0
    
    valori=df.loc[cells,:].mean()<2*epsilon
    
    
    for gene in genes:
        if valori.loc[gene] and dfBulk.loc[gene]>epsilon:
            cont=cont+1
            df.loc[:,gene]=dfBulk.loc[gene]
            
    print("Total genes changes:",str(cont))
    adata2.X=df.values
       
    adata2=adata2[adata2.obs.index!=bulk_el,:]
        
    return adata2,cont