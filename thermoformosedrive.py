#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np


# In[ ]:


import pandas as pd


# In[ ]:


from equilibrator_api import ComponentContribution, Q_
CC = ComponentContribution()


# In[ ]:


from equilibrator_assets.generate_compound import create_compound, get_or_create_compound


# In[ ]:


from ast import literal_eval


# In[ ]:


def CacheGen(filepath):
    df = pd.read_csv(filepath, sep='\t') 
    
    compounds = []
    Reagents = []
    Products = []
    
    for i in range(len(df['Index'])):
        reagents = literal_eval(df['Reagents'][i])
        Reagents.append(reagents)
        products = literal_eval(df['Products'][i])
        Products.append(products)

        for j in range(len(reagents)):
            if reagents[j] not in compounds:
                compounds.append(reagents[j])
        for k in range(len(products)):
            if products[k] not in compounds:
                compounds.append(products[k])
        print(i)
                    
    compound_cache = get_or_create_compound(CC.ccache, compounds, mol_format="smiles", error_log='./ErrorLog.tsv')
    
    return(compound_cache, compounds)


# In[ ]:


def ThermoGen(filepath, compound_cache, compounds, name):
    rels = pd.read_csv(filepath, sep='\t')
    indexes = []
    rules = []
    Reagents = []
    Products = []
    for i in range(len(rels['Index'])):
        indexes.append(rels['Index'][i])
        rules.append(rels['Rule'][i])
        reagents = literal_eval(rels['Reagents'][i])
        Reagents.append(reagents)
        products = literal_eval(rels['Products'][i])
        Products.append(products)
        
    mus = []
    sigma_vecs = []
    for c in compound_cache:
        mu = (CC.predictor.preprocess.get_compound_prediction(c))[0]
        sigma_vec = (CC.predictor.preprocess.get_compound_prediction(c))[1]
        mus.append(mu)
        sigma_vecs.append(sigma_vec)
    
    error_log = pd.read_csv('./ErrorLog.tsv', sep='\t')
    final_compounds = []
    for i in range(len(compounds)):
        if error_log['status'][i] == 'valid':
            final_compounds.append(compounds[i])

    print(len(mus))
    print(len(final_compounds))
    
    EnergyChanges = []
    for i in range(len(rels['Index'])):
        dummy_mus = []
        dummy_sigma_vecs = []
        dummy_compounds = []
        dummy_coefficients = []
        reagents = literal_eval(rels['Reagents'][i])
        products = literal_eval(rels['Products'][i])
        for j in range(len(reagents)):
            dummy_compounds.append(reagents[j])
            dummy_coefficients.append(-1)
        for k in range(len(products)):
            dummy_compounds.append(products[k])
            dummy_coefficients.append(1)
        valid_reaction = True
        for m in range(len(dummy_compounds)):
            if dummy_compounds[m] not in final_compounds:
                valid_reaction = False
                break
            else: 
                dummy_mus.append(mus[final_compounds.index(dummy_compounds[m])])
                dummy_sigma_vecs.append(sigma_vecs[final_compounds.index(dummy_compounds[m])])
        if valid_reaction == True:
            S = np.zeros(len(dummy_compounds))
            for n in range(len(dummy_coefficients)):
                S[n] = dummy_coefficients[n]
            dummy_mus = Q_(dummy_mus, "kJ/mol")
            dummy_sigma_vecs = Q_(dummy_sigma_vecs, "kJ/mol")
            standard_dgs = S.T @ dummy_mus
            U = S.T @ dummy_sigma_vecs
            EnergyChanges.append(standard_dgs._magnitude.round(2))
        else:
            EnergyChanges.append('NaN')
    
    outputdata = {'Index':indexes, 'Reagents':Reagents, 'Products':Products, 'Rule':rules, 'Energy Change':EnergyChanges}
    outputdf = pd.DataFrame(outputdata)
    outputdf.to_csv(f'{name}RelsThermo.tsv', header=None, index=None, sep='\t', mode='a')
    return(outputdf)


# In[ ]:


get_ipython().run_cell_magic('time', '', "a, b = CacheGen('formosedriveProcessedRels.tsv')\n")


# In[ ]:


get_ipython().run_cell_magic('time', '', "a = ThermoGen('formosedriveProcessedRels.tsv', a, b, 'Formose_drive')\n")


# In[ ]:





# In[ ]:





# In[ ]:




