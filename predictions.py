#!/usr/bin/env python
# coding: utf-8

# In[64]:


import numpy as np
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import pandas as pd
import os, os.path
from collections import defaultdict


def initialize_matrix(peptide_len, alphabet):
    init_matrix = [0]*peptide_len

    for i in range(0, peptide_len):
        row = {}
        for letter in alphabet: 
            row[letter] = 0.0
        init_matrix[i] = row
        
    return init_matrix


def from_psi_blast(file_name, peptide_len):

    f = open(file_name, "r")
    nline = 0
    for line in f:
        sline = str.split( line )
        if nline == 0:
        # recover alphabet
            alphabet = [str]*len(sline)
            for i in range(0, len(sline)):
                alphabet[i] = sline[i] 
            matrix = initialize_matrix(peptide_len, alphabet)
        else:
            i = int(sline[0])
            for j in range(2,len(sline)):
                matrix[i-1][alphabet[j-2]] = float(sline[j])       
        nline+= 1
            
    return matrix


def score_peptide(peptide, matrix):
    acum = 0
    for i in range(0, len(peptide)):
        acum += matrix[i][peptide[i]]
    return acum


def main():
    
    print("Calculating the predictive values...")
    
    data_dir = "/Users/laurasansc/github/algorithms_bioinformatics/data/"
    
    # Groups of data, splitted before
    groups = ["train", "test", "train_cluster_centroids"]
    methods = ["hobohm","heuristics"]
    original_df = pd.read_csv(data_dir+"filtered_dataset.csv")
    pcc_df = original_df["allele"].unique()

    # Read evaluation files  (test)
    pcc_dict = defaultdict(list)
    for group in groups:
        for method in methods:
            if group == "train_cluster_centroids" and method == "heuristics":
                break
            for evaluation_file in os.listdir(data_dir+group+'/'):
                if evaluation_file.endswith("_"+group+".csv"): 
                    allele = evaluation_file.replace("_"+group,"")
                    evaluation_df = pd.read_csv(data_dir+group+'/'+evaluation_file)
                    # Read sequence and transform ic50
                    raw_evaluation_peptides = evaluation_df['sequence'].to_list()
                    #  x = 1 - log(IC50)/log(50000)
                    raw_evaluation_targets = (1 - (np.log(evaluation_df.ic50)/np.log(50000))).to_list()

                    #evaluation_peptides, evaluation_targets

                    peptide_len = 9
                    evaluation_peptides = []
                    evaluation_targets = []
                    for i in range(0, len(raw_evaluation_peptides)):
                        if len(raw_evaluation_peptides[i]) == peptide_len:
                            evaluation_peptides.append(raw_evaluation_peptides[i])
                            evaluation_targets.append(raw_evaluation_targets[i])

                    # Define which PSSM file to use (file save from pep2mat)

                    pssm_file = data_dir + method + "_PSSM/w_matrix_" + allele 

                    w_matrix = from_psi_blast(pssm_file, peptide_len)

                    evaluation_predictions = []
                    for i in range(len(evaluation_peptides)):
                        score = score_peptide(evaluation_peptides[i], w_matrix)
                        evaluation_predictions.append(score)
                        #print(evaluation_peptides[i], score, evaluation_targets[i])

                    pcc = pearsonr(evaluation_targets, evaluation_predictions)
                    #print(allele + "\t")
                    #print("PCC: ", pcc[0])


                    pcc_dict[allele[:-4]].append(round(pcc[0], 4))

                    #print(pcc_dict)

                    s = "PCC: "+ str(round(pcc[0], 4))

                    plt.clf()
                    plt.close('all')
                    plt.figure(figsize=(7,5))
                    plt.scatter(evaluation_targets, evaluation_predictions, color="gray")
                    plt.title('Pearson Correlation ' + allele[:-4])
                    plt.ylim(0)
                    plt.xlabel('Evaluation targets transformed IC50')
                    plt.ylabel('Evaluation predictions')
                    plt.text(max(evaluation_targets)/(max(evaluation_targets)+(max(evaluation_targets)/2)), max(evaluation_predictions)/(max(evaluation_predictions)*0.40), s)
                    fname = data_dir+ "pcc_plots/"+ allele[:-4] + "_"+ method + "_"+ group + '.png'  
                    plt.savefig(fname, format='png', dpi=100)

    pcc_df = pd.DataFrame.from_dict(pcc_dict, orient='index', columns=['Hobohm_train', 'Heuristics_train', 'Hobohom_test', 'Heuristics_test', 'Hobohm_train_cluster_centroids'])
    pcc_df['allele'] = pcc_df.index
    pcc_df= pcc_df[['allele','Hobohm_train', 'Heuristics_train', 'Hobohom_test', 'Heuristics_test', 'Hobohm_train_cluster_centroids']]


    pcc_df.to_csv(data_dir +"pcc_results.csv", index=False)


if __name__ == "__main__":
    main()