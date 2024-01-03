# Final output of processing the dataset

The final output from all the processing in this folder is: Complete_dataset_results.csv. Its description is as follows: \
Complete_dataset_results.csv description: 
* Data type: PANDAS dataframe 
* Column: Reaction SMILES - This is the original reaction smiles associated with the reaction 
* Column: reactionID - This is the reaction ID associated with every reaction 
* Column: react_smiles - After reaction ennumeration, this is the principal reactant under consideration 
* Column: rxn_smiles - This is the reaction SMILES that results from reaction processing/ ennumeration 
* Column: template - This template was extracted from the reaction 
* Column: Template Verification (True/False) - True means that the extracted template, when applied to the reactant, generates the product of interest! 


# Other file details

model.ckpt-10654.as_numpy.json.gz - This is the SCScore model developed by Coley et al. It is just here for convenience. \
Complete_processing_v1_GitHub.py - This takes the complete_dataset.csv as input. Then it processes the reactions, extracts the templates, and checks that the extracted templates are appropriate. \
helper_reaction_processing_v2_GitHub.py - This processes the reaction SMILES so that there is a single principal reactant. \
helper_template_extraction_v1_GitHub.py - This extracts the template for every reaction. \
helper_template_application_check_v2_GitHub.py - It applies every template to make sure it is okay. \
Complete_dataset_results.csv - This is the final dataset after processing efforts. \
template_extractor_newSettings.py - This file contains the change in settings to RDChiral. Important to use this for template extraction for maximum accuracy!.
