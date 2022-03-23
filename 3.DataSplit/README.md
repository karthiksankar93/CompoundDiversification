# Data Splitting Strategy
The dataset was split randomly into training, validation, and testing splits. This code takes as input the processed dataset. The output from this code is available in this folder as dataset_input.csv

# dataset_input.csv details
* Data type: PANDAS dataframe 
* Column: Reaction SMILES - This is the original reaction smiles associated with the reaction 
* Column: reactionID - This is the reaction ID associated with every reaction 
* Column: react_smiles - After reaction ennumeration, this is the principal reactant under consideration 
* Column: rxn_smiles - This is the reaction SMILES that results from reaction processing/ ennumeration 
* Column: template - This template was extracted from the reaction 
* Column: Template Verification (True/False) - True means that the extracted template, when applied to the reactant, generates the product of interest!
* dataset - This is the dataset split into training, validation, and test splits
* Reaction Partner (s) - This is a list containing potential reaction partners for a given reaction. This information is implicitly available in the template. This is originally a list that was converted to string while saving to csv. Look into 'from ast import literal_eval' to reconvert the string to list.
