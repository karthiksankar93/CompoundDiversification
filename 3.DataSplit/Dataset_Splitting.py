#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  2 11:07:16 2022

@author: karthiksankar2
"""
#import pandas
import pandas as pd
import os
from pathlib import Path

# read the processed dataset
def get_data_df (fpath = 'Complete_dataset_results.csv'):
    return pd.read_csv(fpath, index_col = 0)

#get the root path!
root_path = Path(os.path.dirname(os.path.abspath(__file__))).parent.absolute()

#read in the dataset
dataset_results = get_data_df (os.path.join(root_path, 'ProcessedDataset', 'Complete_dataset_results.csv'))

# keep entries where you are able to extract and apply templates
dataset_results_v2 = dataset_results [dataset_results['Template Verification (True/False)'] == True]

# import numpy
import numpy as np
import time

# split the datset
def split_data_df(data, val_frac=0.1, test_frac=0.1, shuffle=True, seed=123):
    # Define shuffling
    if shuffle:
        if seed is None:
            np.random.seed(int(time.time()))
        else:
            np.random.seed(seed)
        def shuffle_func(x):
            np.random.shuffle(x)
    else:
        def shuffle_func(x):
            pass 
    
    indeces = data.index.tolist()
    N = len(indeces)
    print('{} rows total'.format(N))

    shuffle_func(indeces)
    train_end = int((1.0 - val_frac - test_frac) * N)
    val_end = int((1.0 - test_frac) * N)

    for i in indeces[:train_end]:
        data.set_value(i, 'dataset', 'train')
    for i in indeces[train_end:val_end]:
        data.set_value(i, 'dataset', 'val')
    for i in indeces[val_end:]:
        data.set_value(i, 'dataset', 'test')
    print(data['dataset'].value_counts())
    
#split the dataset
split_data_df (dataset_results_v2)

#save the dataset
dataset_results_v2.to_pickle ('dataset_input_1000.pkl')