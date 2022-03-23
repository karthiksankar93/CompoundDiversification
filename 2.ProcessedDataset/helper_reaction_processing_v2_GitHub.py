# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 15:19:02 2021

@author: karthiksankar2
"""

import rdkit.Chem as Chem
import re
import operator


'''
This is a standalone, importable SCScorer model. It does not have tensorflow as a
dependency and is a more attractive option for deployment. The calculations are
fast enough that there is no real reason to use GPUs (via tf) instead of CPUs (via np)
'''

import math, sys, random, os
import numpy as np
import time
import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
import json
import gzip
import six

import os

score_scale = 5.0
min_separation = 0.25

FP_len = 1024
FP_rad = 2

def sigmoid(x):
  return 1 / (1 + math.exp(-x))

class SCScorer():
    def __init__(self, score_scale=score_scale):
        self.vars = []
        self.score_scale = score_scale
        self._restored = False

    def restore(self, weight_path='model.ckpt-10654.as_numpy.json.gz', FP_rad=FP_rad, FP_len=FP_len):
        self.FP_len = FP_len; self.FP_rad = FP_rad
        self._load_vars(weight_path)
        #print('Restored variables from {}'.format(weight_path))

        if 'uint8' in weight_path or 'counts' in weight_path:
            def mol_to_fp(self, mol):
                if mol is None:
                    return np.array((self.FP_len,), dtype=np.uint8)
                fp = AllChem.GetMorganFingerprint(mol, self.FP_rad, useChirality=True) # uitnsparsevect
                fp_folded = np.zeros((self.FP_len,), dtype=np.uint8)
                for k, v in six.iteritems(fp.GetNonzeroElements()):
                    fp_folded[k % self.FP_len] += v
                return np.array(fp_folded)
        else:
            def mol_to_fp(self, mol):
                if mol is None:
                    return np.zeros((self.FP_len,), dtype=np.float32)
                return np.array(AllChem.GetMorganFingerprintAsBitVect(mol, self.FP_rad, nBits=self.FP_len,
                    useChirality=True), dtype=np.bool)
        self.mol_to_fp = mol_to_fp

        self._restored = True
        return self

    def smi_to_fp(self, smi):
        if not smi:
            return np.zeros((self.FP_len,), dtype=np.float32)
        return self.mol_to_fp(self, Chem.MolFromSmiles(smi))

    def apply(self, x):
        if not self._restored:
            raise ValueError('Must restore model weights!')
        # Each pair of vars is a weight and bias term
        for i in range(0, len(self.vars), 2):
            last_layer = (i == len(self.vars)-2)
            W = self.vars[i]
            b = self.vars[i+1]
            x = np.matmul(x, W) + b
            if not last_layer:
                x = x * (x > 0) # ReLU
        x = 1 + (score_scale - 1) * sigmoid(x)
        return x

    def get_score_from_smi(self, smi='', v=False):
        if not smi:
            return ('', 0.)
        fp = np.array((self.smi_to_fp(smi)), dtype=np.float32)
        if sum(fp) == 0:
            if v: print('Could not get fingerprint?')
            cur_score = 0.
        else:
            # Run
            cur_score = self.apply(fp)
            if v: print('Score: {}'.format(cur_score))
        mol = Chem.MolFromSmiles(smi)
        if mol:
            smi = Chem.MolToSmiles(mol, isomericSmiles=True, kekuleSmiles=True)
        else:
            smi = ''
        return (smi, cur_score)

    def _load_vars(self, weight_path):
        if weight_path.endswith('pickle'):
            import cPickle as pickle
            with open(weight_path, 'rb') as fid:
                self.vars = pickle.load(fid)
                self.vars = [x.tolist() for x in self.vars]
        elif weight_path.endswith('json.gz'):
            with gzip.GzipFile(weight_path, 'r') as fin:    # 4. gzip
                json_bytes = fin.read()                      # 3. bytes (i.e. UTF-8)
                json_str = json_bytes.decode('utf-8')            # 2. string (i.e. JSON)
                self.vars = json.loads(json_str)
                self.vars = [np.array(x) for x in self.vars]

#load the model
model = SCScorer()
model.restore('model.ckpt-10654.as_numpy.json.gz')

def SCScore_Calculator (smiles_input):
    '''This module takes SMILES string as input, and provides as output the SCScore. The original algorithm
    used the SCScorer in ASKCOS https://askcos.mit.edu/. Here, I use a standalone version of SCScorer; expect
    minor differences in results'''
    
    #get SCScore from model
    (smi, sco) = model.get_score_from_smi(smiles_input)
    
    #return statement with results
    return (sco)
    
def rxn_process (smiles):
    
    # splits the reactants and products
    all_reactants, all_products = smiles.split('>>')
    
    reactants = [Chem.MolFromSmiles(smi) for smi in all_reactants.split('.')]
    
    #reactSCScore dictionary initialization
    reactSCScore = {}
    
    #get information about the atom map numbers in the product
    prod_maps_list_0 = []
    for prod in all_products.split('.'):
        prod_maps_list_0.extend (re.findall('\:([[0-9]+)\]', prod))
    prod_maps_0 = set (prod_maps_list_0)
    
    for react in reactants:

        # Make sure the react molecule contributes to the product
        react_maps_0 = set (re.findall('\:([[0-9]+)\]', Chem.MolToSmiles(react)))
        
        if len (react_maps_0.intersection (prod_maps_0))==0:
            continue
    
        react_SCScore_calc_temp = SCScore_Calculator (Chem.MolToSmiles (react))
        reactSCScore [Chem.MolToSmiles (react)] = react_SCScore_calc_temp
        
    #print reactant with highest SCScore
    principal_reactant = max (reactSCScore.items(), key = operator.itemgetter(1))[0]
    
    # Multiple products = enumerate
    for react in reactants:
    
        # Make sure all have atom mapping
        if Chem.MolToSmiles (react) != principal_reactant:
            continue
        
        #get reactant SMILES
        react_smi = Chem.MolToSmiles (react, True)
    
        #Re-parse products
        products = [Chem.MolFromSmiles(smi) for smi in all_products.split('.')]
    
        # Get rid of products when they don'e contribute to this reactant
        react_maps = set (re.findall('\:([[0-9]+)\]', react_smi))
        product_smi_list = []
        prod_list = []
        for mol in products:
            used = False
            for a in mol.GetAtoms():
                if a.HasProp('molAtomMapNumber'):
                    if a.GetProp('molAtomMapNumber') in react_maps:
                        used = True
                    else:
                        a.ClearProp('molAtomMapNumber')
            if used:
                product_smi_list.append (Chem.MolToSmiles(mol,True))
                prod_list.append (mol)
        
        # Overall goal: Remove atom map numbers present in the reactants but not in products (step 0)
        
        # Step 1: What are the atom map numbers in products
        prod_maps_list = []
        for prod in prod_list:
            prod_maps_list.extend (re.findall('\:([[0-9]+)\]', Chem.MolToSmiles (prod, True)))
        
        prod_maps = set (prod_maps_list)
        
        # Step 2: Remove atom map numbers in the reactant that don't contribute to product atoms
        for a in react.GetAtoms():
            if a.HasProp('molAtomMapNumber'):
                if a.GetProp('molAtomMapNumber') not in prod_maps:
                    a.ClearProp('molAtomMapNumber')
            
        products_smi = '.'.join(product_smi_list)
        
        react_smi_2 = Chem.MolToSmiles (react, True)
        
        if products_smi == react_smi:
            continue
            
        #print unmapped reactants too
        [a.ClearProp('molAtomMapNumber') for a in react.GetAtoms()]
        react_smiles = Chem.MolToSmiles(react, True)
        rxn_smiles = '{}>>{}'.format (react_smi_2, products_smi)
        return (react_smiles, rxn_smiles)

if __name__ == "__main__":
    output = rxn_process ('[C:1]([CH3:2])([CH3:3])=[N:4][NH:5][C:6]([c:7]1[c:8]([CH3:15])[c:9]([O:13][CH3:14])[cH:10][cH:11][cH:12]1)=[O:16].[Cl:17][O:18][N:19]=[CH:20][c:21]1[cH:22][c:23]([CH3:28])[cH:24][c:25]([CH3:27])[cH:26]1.[Cl:29][CH:30]([Cl:31])[Cl:32].[K+:33].[K+:34].[O-:35][C:36]([O-:37])=[O:38].[OH2:39]>>[C:1]1([CH3:2])([CH3:3])[N:4]([NH:5][C:6]([c:7]2[c:8]([CH3:15])[c:9]([O:13][CH3:14])[cH:10][cH:11][cH:12]2)=[O:16])[C:20]([c:21]2[cH:22][c:23]([CH3:28])[cH:24][c:25]([CH3:27])[cH:26]2)=[N:19][O:18]1')
    print ('The react_smiles is: {}'.format(output [0]))
    print ('The rxn_smiles is: {}'.format (output [1]))