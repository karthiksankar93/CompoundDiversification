# -*- coding: utf-8 -*-
"""
Created on Sat Sep  4 08:27:43 2021

@author: karthiksankar2
"""

import rdkit.Chem as Chem
from rdchiral.initialization import rdchiralReaction, rdchiralReactants
from rdchiral.main import rdchiralRun

def template_application (inputs):
    
    react_smiles = inputs[0]
    template = inputs[1]
    rxn_smiles = inputs[2]
    
    try:
        #load the reactant
        rct = rdchiralReactants (react_smiles)
        
        #load the reaction
        rxn = rdchiralReaction (template)
        
        #run the reaction with the template
        outcomes = rdchiralRun (rxn ,rct , combine_enantiomers=False)
        
        #get product goal
        prod_goal = Chem.MolFromSmiles (rxn_smiles.split('>')[2])
        [a.ClearProp('molAtomMapNumber') for a in prod_goal.GetAtoms()]
        prod_goal = Chem.MolToSmiles(prod_goal, True)
        #sometimes stereochem takes another canonicalization
        prod_goal = Chem.MolToSmiles(Chem.MolFromSmiles(prod_goal), True)
        
        if prod_goal in outcomes:
            result = True    
        else:
            result = False
            
    except:
        result = False
    
    return result