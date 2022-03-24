# -*- coding: utf-8 -*-
"""
Created on Wed Oct 20 19:03:15 2021

@author: karthiksankar2
"""
# all import statements
import pandas as pd 
import os
import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
from rdchiral.initialization import rdchiralReaction, rdchiralReactants
from rdchiral.main import rdchiralRun
from rdkit import DataStructs
import numpy as np
from joblib import Parallel, delayed

def get_fp_parallel (smi, getfp_label):
    
    #set fingerprint type    
    if getfp_label == 'Morgan2noFeat':
        getfp = lambda smi: AllChem.GetMorganFingerprint(Chem.MolFromSmiles(smi), 2, useFeatures=False)
    elif getfp_label == 'Morgan3noFeat':
        getfp = lambda smi: AllChem.GetMorganFingerprint(Chem.MolFromSmiles(smi), 3, useFeatures=False)
    elif getfp_label == 'Morgan2Feat':
        getfp = lambda smi: AllChem.GetMorganFingerprint(Chem.MolFromSmiles(smi), 2, useFeatures=True)
    elif getfp_label == 'Morgan3Feat':
        getfp = lambda smi: AllChem.GetMorganFingerprint(Chem.MolFromSmiles(smi), 3, useFeatures=True)
    else:
        raise ValueError('Unknown getfp label')
                
    return getfp (smi)

def do_one_bulk (similarity_label, getfp_label, val_index_list):
    
    '''This function will take as input a list of validation indeces. Then, it will run the search algorithm (in this case do_one_v4) one at a time for all
    validation indeces. The function output is also a list of all do_one_v4 outputs corresponding to the validation indeces of the input list'''
    
    output_list = []
    
    for index in val_index_list:
        output_list.append (do_one_v4(similarity_label, getfp_label, index))
    
    return output_list

def do_one_v4 (similarity_label, getfp_label, ix,  max_prec=100, testlimit = 100):
    ''' This takes as input the similarity label, fingerprint label, and the test/validation set index.
    It runs similarity-based diversification for the test/validation compound of the given index.
    It returns (found_rank, return_rxn, sim_to_idealProdPrevious, sim_reactant_to_idealProd). 'found_rank' is the rank at which the correct solution was found.
    'return_rxn' is the reaction SMILES associated with the correct reaction. 'sim_to_idealProdPrevious' is the similarity between the proposed product and 
    the recorded product. 'sim_reactant_to_idealProd' is the similarity between the recorded reactant and the recorded product.
    If the algorithm runs fine mechanically, but is simply unable to find a good solution --> it will return (9999,9999,sim_to_idealProdPrevious, sim_reactant_to_idealProd).
    If the algorithm is unable to even run, it will return (9998,9998,9998,9998)'''
    
    
    #setup a similarity metric
    if similarity_label == 'Tanimoto':
        similarity_metric = DataStructs.BulkTanimotoSimilarity
    elif similarity_label == 'Dice':
        similarity_metric = DataStructs.BulkDiceSimilarity
    elif similarity_label == 'TverskyA': # weighted towards punishing onlyA
        def similarity_metric(x, y):
            return DataStructs.BulkTverskySimilarity(x, y, 1.5, 1.0)
    elif similarity_label == 'TverskyB': # weighted towards punishing onlyB
        def similarity_metric(x, y):
            return DataStructs.BulkTverskySimilarity(x, y, 1.0, 1.5)
    else:
        raise ValueError('Unknown similarity label')

    #set fingerprint type    
    if getfp_label == 'Morgan2noFeat':
        getfp = lambda smi: AllChem.GetMorganFingerprint(Chem.MolFromSmiles(smi), 2, useFeatures=False)
    elif getfp_label == 'Morgan3noFeat':
        getfp = lambda smi: AllChem.GetMorganFingerprint(Chem.MolFromSmiles(smi), 3, useFeatures=False)
    elif getfp_label == 'Morgan2Feat':
        getfp = lambda smi: AllChem.GetMorganFingerprint(Chem.MolFromSmiles(smi), 2, useFeatures=True)
    elif getfp_label == 'Morgan3Feat':
        getfp = lambda smi: AllChem.GetMorganFingerprint(Chem.MolFromSmiles(smi), 3, useFeatures=True)
    else:
        raise ValueError('Unknown getfp label')
    
    # overall try/ except capture for reactions that fail in the algorithm
    try:
    
        #get reactant smiles
        react_smiles = []
        react_smiles.append (datasub_val ['react_smiles'][ix])

        #loads reactant SMILES into RDChiral object
        rct = rdchiralReactants(datasub_val ['react_smiles'][ix])

        #get the fingerprint of the product
        fp = datasub_val['react_fp'][ix]

        #calculates similarity metric between fingerprint 
        # and all fingerprints in the database
        sims = similarity_metric (fp, [fp_ for fp_ in datasub['react_fp']])

        #sort the similarity metric in reverse order
        js = np.argsort(sims) [::-1]

        #get the product goal from the reaction smiles string
        prod_goal = Chem.MolFromSmiles (datasub_val['rxn_smiles'][ix].split('>')[2])
        [a.ClearProp('molAtomMapNumber') for a in prod_goal.GetAtoms()]
        prod_goal = Chem.MolToSmiles(prod_goal, True)

        #sometimes stereochem takes another canonicalization
        prod_goal = Chem.MolToSmiles(Chem.MolFromSmiles(prod_goal), True)

        # Get probability of precursors
        probs = {}
        coreact = {}

        for ji,j in enumerate (js[:max_prec]):

            #get the index in the similarity search
            jx = datasub.index[j]

            try:  
                #extract the template
                template = datasub['template'][jx]

                #get product reference fingerprint
                prods_ref_fp = getfp(datasub['rxn_smiles'][jx].split('>')[2])

                #load the reaction into an rdchiral object
                rxn = rdchiralReaction(template)


            except:

                #if it fails, move to the next reaction
                continue 

            try:
                # try running the rdchiral reaction
                outcomes = rdchiralRun (rxn, rct, combine_enantiomers=False)

            #if running the rdchiral reaction is not possible go into this exception routine
            except Exception as e:
                print(e)
                outcomes = []
                continue

            for products in outcomes:

                #proposed product fingerprint
                products_fp = getfp(products)

                #compare proposed and precedent product fp
                products_sim = similarity_metric(products_fp, [prods_ref_fp])[0]
                
                #compute overall similarity
                if products in probs:
                    if (products_sim*sims[j]) > probs[products]:
                        probs[products] = products_sim*sims[j]
                        coreact[products] = datasub ['Reaction Partner (s)'][jx]
                else:
                    probs[products] = products_sim*sims[j]
                    coreact[products] = datasub ['Reaction Partner (s)'][jx]
        
        #initialize found_rank            
        found_rank = 9999
        
        #initialize variable
        sim_to_idealProdPrevious = -1

        #convert prod_goal to its fingerprint
        prod_goal_fp = getfp(prod_goal)


        for r, (prod, prob) in enumerate(sorted(probs.items(), key=lambda x:x[1], reverse=True)[:testlimit]):         

            #convert prod to its fingeprint
            prod_fp = getfp (prod)

            #compute similarity to ideal product
            sim_to_idealProd = similarity_metric (prod_fp, [prod_goal_fp])[0]
            
            #if it is better than the previous identified solution, save!
            if sim_to_idealProd > sim_to_idealProdPrevious:
                sim_to_idealProdPrevious = sim_to_idealProd
                found_rank = r + 1
                recoveredProduct = prod
                coreactant = coreact[prod]

        #get reactant fingerprint
        reactant_fp = getfp(datasub_val ['react_smiles'][ix])
        sim_reactant_to_idealProd = similarity_metric (reactant_fp, [prod_goal_fp])[0]

        #assemble the reaction
        return_rxn = ''
        
        #add the reactants to return_rxn
        return_rxn += datasub_val ['react_smiles'][ix]

        #add coreactant information to return_rxn
        for reactant in coreactant:
            return_rxn += '.' + reactant

        #add product information to return_rxn
        return_rxn += '>>' + recoveredProduct

        #if it is not a great suggestion, return 9999
        if sim_to_idealProdPrevious <= sim_reactant_to_idealProd:
            found_rank = 9999
            return_rxn = 9999

        #return results
        return (found_rank, return_rxn, sim_to_idealProdPrevious, sim_reactant_to_idealProd)
    
    #unexpected error during the search- return 9998
    except:
        return (9998, 9998, 9998, 9998)

    
# Evaluate the results!!
def ranks_to_acc(found_at_rank, fid=None):
    
    tot = float(len(found_at_rank))
    accs = []
    for n in [1, 3, 5, 10, 20, 50,9997, 9998, 9999]:
        accs.append(sum([r <= n for r in found_at_rank]) / tot)
        
    return accs

# read the  dataset
def get_data_df (fpath = 'dataset_input.csv'):
    return pd.read_csv(fpath, index_col = 0)
       
if __name__ == '__main__':
    
    import math
    from ast import literal_eval
    from pathlib import Path
    
    #validation parameters to test
    all_getfp_labels = ['Morgan2Feat']
    all_similarity_labels = ['Tanimoto']
    
    #enter the number of desired cores
    n_cores_desired = 18    
    
    #get the root path!
    root_path = Path(os.path.dirname(os.path.abspath(__file__))).parent.absolute()
    
    #read in the dataset
    data = get_data_df (os.path.join(root_path, '3.DataSplit', 'dataset_input.csv'))
    
    #convert the string in the dataframe into list (result of storing the dataframe as csv!)
    data['Reaction Partner (s)'] =  data.str_Reaction_Partner.apply(lambda x: literal_eval(str(x)))
    
    #drop the string column str_Reaction_Partner
    data.drop(columns = ['str_Reaction_Partner'], inplace = True)
    
    #split the dataset
    datasub = data[data['dataset']== 'train']
    
    #loop through all fingerprint labels
    for fp in all_getfp_labels:
        
        # loop through similarity labels
        for sim in all_similarity_labels:
            
            #set up fingerprint label
            getfp_label = fp
            
            #set up similarity label
            similarity_label = sim
            
            #set the test dataframe. It is just called datasub_val for convenience- the code is from validation runs!
            datasub_val = data[data['dataset']== 'test']
            
            #compute the fingerprints in a parallelized fashion
            print ('Getting fingerprint information for the training set: {}'.format (getfp_label))
            train_fp = Parallel (n_jobs=20, verbose=3)(delayed(get_fp_parallel)(smi, getfp_label) for smi in datasub['react_smiles'])
            datasub['react_fp'] = train_fp
    
            print ('Getting fingerprint information for the test set: {}'.format (getfp_label))
            val_fp = Parallel (n_jobs=20, verbose=3)(delayed(get_fp_parallel)(smi, getfp_label) for smi in datasub_val ['react_smiles'])
            datasub_val['react_fp'] = val_fp   
            
            #get a list of validation indeces
            val_list = list(datasub_val.index)
        
            #split validation list into n_cores_desired separate lists, batches the large job into smaller parts!
            start = 0
            end = len (val_list)
            n_equal_split = math.ceil (len (val_list)/n_cores_desired)
            
            #it is best if each core had atleast 1000 examples
            n = max (1000, n_equal_split)
            input_list = []
            for i in range (start,end,n):
                x = i
                input_list.append(val_list[x:x+n])
            
            #number of cores must equal the length of the input_list
            n_cores = len (input_list)
            
            #parallelize
            output = Parallel(n_jobs=n_cores, verbose=10)(delayed(do_one_bulk)(similarity_label, getfp_label, i) for i in input_list)
            
            #combine the output list
            output_cum = []
            for item in output:
                output_cum.extend(item)
            
            #parse the output
            rank = []
            rxn = []
            sim_prod_prod = []
            sim_reac_prod = []
        
            for item in output_cum:
                rank.append (item[0])
                rxn.append (item[1])
                sim_prod_prod.append (item[2])
                sim_reac_prod.append (item[3])
            
            #add the list to the pandas dataframe
            datasub_val ['Rank'] = rank
            datasub_val ['Reaction'] = rxn
            datasub_val ['Similarity (proposed product vs recorded product)'] = sim_prod_prod
            datasub_val ['Similarity (recorded reactant vs recorded product)'] = sim_reac_prod            
            
            #drop the react_fp
            datasub_val.drop (columns = 'react_fp', inplace = True)
            
            #change to out folder
            os.chdir (os.path.join(root_path, 'TopK_Accuracy_Evaluation', 'out_test'))
            
            #save the results
            datasub_val.to_pickle ('datasub_test_result_{}_{}.pkl'.format(getfp_label,similarity_label))
            datasub_val.to_csv ('datasub_test_result_{}_{}.csv'.format(getfp_label,similarity_label))
