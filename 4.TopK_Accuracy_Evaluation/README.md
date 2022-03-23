# TopK_Accuracy_Evaluation

* Top_k_accuracy_v7_parallel_val_GitHub.py - This is the code I used for validation/ learning the appropriate fingerprint + similarity metric to use. It uses the validation dataset + training dataset. It tests a range of fingerprints (different radii, feature settings, etc.) and similarity metrics (Tanimoto, Dice, etc.).
* Top_k_accuracy_v7_parallel_test_GitHub.py - This is the code I used for testing. It uses the test dataset + training dataset. It uses Morgan Fingerprint of radius 2, with features. It uses tanimoto similarity metric.
* Folder: out_test - This folder contains the results of the test set in both .csv and .pkl format. The top-k accuracy analysis, and associated plots (fast filter plot, similarity of proposed product vs. recorded product plot) are available in Test_Set_Analysis.ipynb. api_client.py and FastFilter.ipynb both are used to compute fast filter scores for the proposed reactions.
* Folder: out_val - This folder contains the results of the validation set in both .csv and .pkl format. 

