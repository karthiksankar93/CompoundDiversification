# RawDataset

Original dataset comes from 'Predicting Organic Reaction Outcomes with Weisfeiler-Lehman Network' by Wengong Jin, Connor W. Coley, Regina Barzilay, Tommi Jaakkola.  The dataset can be accessed through https://github.com/connorcoley/rexgen_direct/tree/master/rexgen_direct/data.

The 'train.txt.proc', 'valid.txt.proc', and 'test.txt.proc' were all combined together in 'Complete_dataset.csv'. Further, every reaction has a 'reactionID' assigned to it. This keeps track of the reaction IDs; just in case one needs to come back to the reactions later on!

One might be able to open the files using pandas, and the command pandas.read_csv('Complete_dataset.csv', index_col = 0). In case you are having challenges, the dataset is also available at https://www.dropbox.com/sh/ggsgy4u0pt3xgy2/AAAphnqY-XhnHBDC8S03UX7Aa?dl=0. The file name in the dropbox folder is 'Complete_dataset.csv'.
