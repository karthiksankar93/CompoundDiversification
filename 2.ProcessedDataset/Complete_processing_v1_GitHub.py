# -*- coding: utf-8 -*-
"""
Created on Fri Oct 15 12:54:53 2021

@author: karthiksankar2
"""

import pandas as pd
from joblib import Parallel, delayed
from helper_reaction_processing_v2_GitHub import rxn_process
from helper_template_extraction_v1_GitHub import get_templates
from helper_template_application_check_v2_GitHub import template_application
import os
from pathlib import Path

# read the original PANDAS dataframe from Jin et al.
def get_data_df (fpath = 'Complete_dataset.csv'):
    return pd.read_csv(fpath, index_col = 0)

#get the root path!
root_path = Path(os.path.dirname(os.path.abspath(__file__))).parent.absolute()

#read in the dataset
datasub = get_data_df (os.path.join(root_path, '1.RawDataset', 'Complete_dataset.csv'))


#parallelized reaction processing: pick a principal reactant based on SCScore, then process such that you get  to
#principal reactant >> products
processed_output = Parallel (n_jobs=20, verbose=1)(delayed(rxn_process)(rxn_smi) for rxn_smi in datasub['Reaction SMILES'])

react_smiles = []
rxn_smiles = []

count = 0

for item in processed_output:
    try:    
        react_smiles.append (item[0])
        rxn_smiles.append (item[1])
    except:
        react_smiles.append (None)
        rxn_smiles.append (None)
 
datasub ['react_smiles'] = react_smiles
datasub ['rxn_smiles'] = rxn_smiles

#get rid of null react_smiles (and rxn_smiles) values
bool_series = pd.notnull (datasub['react_smiles'])
datasub = datasub[bool_series]

#Extract all templates:

#a print statement to get started
print("Extracting all templates")

#parallelized template extraction
templates = Parallel (n_jobs=20, verbose=1)(delayed(get_templates)(rxn_smi) for rxn_smi in datasub['rxn_smiles'])

#store the template
datasub["template"]=templates

inputs = []

for ix in list(datasub.index):
    inputs.append ((datasub['react_smiles'][ix], datasub['template'][ix], datasub['rxn_smiles'][ix]))
    
print ('Parallel analysis')
verify = Parallel (n_jobs=20, verbose=5)(delayed(template_application)(i) for i in inputs)

#store the results
print ('Store the results')
datasub['Template Verification (True/False)'] = verify

#save the file
datasub.to_csv ('Complete_dataset_results.csv')
