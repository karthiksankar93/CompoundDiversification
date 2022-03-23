# -*- coding: utf-8 -*-
"""
Created on Fri Sep  3 15:22:13 2021

@author: karthiksankar2

from helper_template_extraction_v1.py

created by me on Tue Jan 26 15:09:40 2021
"""

# Might be good to provide this file manually!
from template_extractor_newSettings import extract_from_reaction

def get_templates(rxn_smi):
    
    # extracts the template    
    try:
        #convert reaction into a dictionary, ensure it is setup to obtain forward templates!
        reaction = {}
        rct_0, rea_0, prd_0 = rxn_smi.split('>')
        reaction['reactants'] = prd_0
        reaction['products'] = rct_0
        reaction['_id'] = 0
                
        #extract the template
        template = extract_from_reaction(reaction)['reaction_smarts']
    
    # fails to extract template
    except:
        template=None
        
    return template
