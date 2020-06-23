#!/usr/bin/env python
# coding: utf-8

# # Clustering sequence weighting

import numpy as np
import pandas as pd
from time import time
import math
import os


# ### Functions

def load_data_utils():
    # Load the alphabet file
    alphabet_file = data_dir + "Matrices/alphabet"

    alphabet = np.loadtxt(alphabet_file, dtype=str)

    # Load the matrices file 
    blosum62_file = data_dir + "Matrices/blosum62.freq_rownorm"
    _blosum62 = np.loadtxt(blosum62_file, dtype=float).T

    blosum62 = {}

    for i, letter_1 in enumerate(alphabet):

        blosum62[letter_1] = {}

        for j, letter_2 in enumerate(alphabet):

            blosum62[letter_1][letter_2] = _blosum62[i, j]
            
    # Load the bg information
    bg_file = data_dir + "Matrices/bg.freq.fmt"
    _bg = np.loadtxt(bg_file, dtype=float)

    bg = {}
    for i in range(0, len(alphabet)):
        bg[alphabet[i]] = _bg[i]

    return alphabet, blosum62, bg
    
def initialize_matrix(peptide_len, alphabet):
    init_matrix = [0]*peptide_len

    for i in range(0, peptide_len):
        row = {}
        for letter in alphabet: 
            row[letter] = 0.0
        init_matrix[i] = row
        
    return init_matrix


def c_matrix(peptides, peptide_len, alphabet):
    
    cmatrix = initialize_matrix(peptide_len, alphabet)

    for position in range(0, peptide_len):
        for peptide in peptides:
            cmatrix[position][peptide[position]] += 1
    return cmatrix


def clustering_weighting(peptides, cmatrix):

    weights = {}
    for peptide in peptides:
        # apply sequence weighting
        weights[peptide] = 1
        
    return weights


def f_matrix(peptides, peptide_len, alphabet):
    
    fmatrix = initialize_matrix(peptide_len, alphabet)

    for position in range(0, peptide_len):
        n = 0;

        for peptide in peptides:
            fmatrix[position][peptide[position]] += weights[peptide]
            n += weights[peptide]

        for letter in alphabet: 
            fmatrix[position][letter] = fmatrix[position][letter]/n

    return fmatrix


def g_matrix(f_matrix, peptide_len, alphabet, blosum62):
    
    gmatrix = initialize_matrix(peptide_len, alphabet)

    for position in range(0, peptide_len):
        for letter_1 in alphabet:
            for letter_2 in alphabet:
              gmatrix[position][letter_1] += fmatrix[position][letter_2] * blosum62[letter_1][letter_2]

    return gmatrix

def p_matrix(f_matrix, g_matrix, peptide_len, alphabet, beta, neff):
    
    pmatrix = initialize_matrix(peptide_len, alphabet)
    alpha = neff - 1

    for position in range(0, peptide_len):
        for a in alphabet:
            pmatrix[position][a] = (alpha*fmatrix[position][a] + beta*gmatrix[position][a]) / (alpha + beta)

    return pmatrix


def w_matrix(p_matrix, bg, peptide_len, alphabet):
    wmatrix = initialize_matrix(peptide_len, alphabet)

    for position in range(0, peptide_len):

        for letter in alphabet:
            if pmatrix[position][letter] > 0:
                wmatrix[position][letter] = 2 * math.log(pmatrix[position][letter]/bg[letter])/math.log(2)
            else:
                wmatrix[position][letter] = -999.9

    return wmatrix


def to_psi_blast(matrix):

    header = ["", "A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
    print ('{:>4} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8}'.format(*header)) 
    letter_order = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

    for i, row in enumerate(matrix):
        scores = []
        scores.append(str(i+1) + " A")

        for letter in letter_order:
            score = row[letter]
            scores.append(round(score, 4))

        print('{:>4} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8}'.format(*scores)) 




def to_psi_blast_file(matrix, file_name):
    
    with open(file_name, 'w') as file:
        header = ["", "A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
        file.write ('{:>4} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8}\n'.format(*header)) 
        letter_order = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

        for i, row in enumerate(matrix):
            scores = []
            scores.append(str(i+1) + " A")

            for letter in letter_order:
                score = row[letter]
                scores.append(round(score, 4))

            file.write('{:>4} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8}\n'.format(*scores)) 


# ### Data directory


data_dir = "/Users/laurasansc/github/algorithms_bioinformatics/data/"


# ### Main

def main():
    print("Getting PSSM from Hobohm clustering algorithm...")
    
    alphabet, blosum62, bg = load_data_utils()

    #sequence_weighting = False
    beta = 50
    peptide_len = 9
    train_dir = data_dir + 'train/'

    # Iterate through data directory to find train files 
    for peptides_file in os.listdir(train_dir):
        if peptides_file.endswith("_train_cluster_centroids.csv"): 
            allele = peptides_file.replace("_train_cluster_centroids.csv","")
            # Read file 
            peptides_df = pd.read_csv(train_dir+peptides_file)
            raw_peptides = peptides_df["sequence"].to_list()

            peptides = []
            for i in range(0, len(raw_peptides)):
                if len(raw_peptides[i]) == peptide_len:          
                    peptides.append(raw_peptides[i])
                #else:
                 #   print("Peptide length too short discard", raw_peptides[i])
            neff = len(peptides)
            ### get count matrices
            cmatrix = c_matrix(peptides, peptide_len, alphabet)
            ### weights 
            weights = clustering_weighting(peptides, cmatrix)
            ### f, g, p and w matrices
            fmatrix = f_matrix(peptides, peptide_len, alphabet)
            gmatrix = g_matrix(fmatrix, peptide_len, alphabet, blosum62)
            pmatrix = p_matrix(fmatrix, g_matrix, peptide_len, alphabet, beta, neff)
            wmatrix = w_matrix(pmatrix, bg, peptide_len, alphabet)

            # Write out PSSM in Psi-Blast format to file
            file_name = data_dir + "HobohmI_PSSM/w_matrix_" + allele + '.csv'
            to_psi_blast_file(wmatrix, file_name)


if __name__ == "__main__":
    main()



