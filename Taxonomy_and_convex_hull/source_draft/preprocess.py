from __future__ import division 
import sys

import numpy as np
import pickle
import time
import random
import os
from tqdm import tqdm 
import pandas as pd

from scipy.spatial import distance
from scipy.stats import scoreatpercentile 

from mpl_toolkits.axes_grid1 import make_axes_locatable

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.gridspec as gridspec

from Bio import SeqIO

from sklearn.model_selection import train_test_split
from sklearn.decomposition import IncrementalPCA






# Read the fasta file
def get_seq(filename, get_header = False):
    assert filename.endswith('.fasta'), 'Not a fasta file.'
    
    records = list(SeqIO.parse(filename, "fasta"))
    records_seq = [i.seq for i in records]
    headers = [i.description for i in records]
    if get_header == True:
        return records_seq, headers
    else:
        return records_seq
    
# Output indices of AAs at each position
def potts_index(sequence): # Output all possible AAs on all positions
    for i in range(len(sequence)):
        assert len(sequence[i]) == len(sequence[0])
        
        aa_at_pos = []
        for i in range(len(sequence[0])):
            tmp_list = []
            for j in range(len(sequence)):
                if (sequence[j][i] in tmp_list) == False:
                    tmp_list.append(sequence[j][i])
            aa_at_pos.append(tmp_list)
        return [''.join(i) for i in aa_at_pos]    


# Convert alphabet sequences to one-hot Potts representation.
def convert_potts(sequence, aa_at_pos):
    N = len(sequence) # samples
    n = len(aa_at_pos) # positions
    q_n = np.array([len(i) for i in aa_at_pos]).astype(int)
    length = sum(q_n)

    v_traj_onehot = np.zeros((N,length)).astype(int)
    for j in range(n):
        sumq = sum(q_n[:j])
        for i in range(N):
            pos_to_1 = sumq + aa_at_pos[j].find(sequence[i][j])
            v_traj_onehot[i,pos_to_1] = 1 # set a position as 1
    return v_traj_onehot, q_n


# append sequences into a list container 
def create_MSA(Seq_samples, samples, positions):
    single_seq, seq_list = [], []
    for jj in range(samples):
        single_seq = "".join((str(ii) for ii in Seq_samples[jj]))
        seq_list.append(single_seq)
    return seq_list


def pandas_list_to_array(df):
    return np.transpose(np.array(df.values.tolist()), (0, 2, 1))

def preprocess_inputs(df, token2int, cols=['Sequence']):
    return pandas_list_to_array(df[cols].applymap(lambda seq: [token2int[x] for x in seq]))

def OneHot_encode(data,sample):
    b = np.eye(21)[data[sample]].transpose(0,2,1)[:,:,0]
    return b















