#!/usr/bin/env python
# coding: utf-8

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns 
from sklearn.model_selection import train_test_split

def main():
    data_dir="/Users/laurasansc/github/algorithms_bioinformatics/data/"

    raw_data_file = data_dir + "raw/predictions_smm.txt"

    raw_data = pd.read_csv(raw_data_file, sep="\t")

    # Group by human
    human_data = raw_data[raw_data['species'].isin(['human'])]

    alleles = list(human_data['allele'].unique())
    #len(alleles)

    allele_to_take = 19
    allele_no = 0
    print("Selecting the alleles ... ")
    for allele in alleles:
        if allele_no <= 19:
            print(allele_no, allele)
            allele_data = human_data.loc[(human_data['allele'] == allele) &
                                         (human_data['ic50'] < 500) &
                                         (human_data['length'] == 9) ].sort_values(by='ic50', ascending=True)
            train, test = train_test_split(allele_data, test_size = 0.2)
            train = pd.DataFrame(train)
            test = pd.DataFrame(test)
            train.to_csv(data_dir + '/train/' + allele.replace('*','_').replace(' ','_')+'_train.csv', index=False)
            test.to_csv(data_dir + '/test/' + allele.replace('*','_').replace(' ','_')+'_test.csv', index=False)
        allele_no = allele_no +1


if __name__ == "__main__":
    main()