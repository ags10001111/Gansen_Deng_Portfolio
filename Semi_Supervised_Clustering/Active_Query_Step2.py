"""
@author: Gansen Deng
"""
import numpy as np
from numpy import dot, arccos, sin
from numpy.linalg import norm
import pandas as pd
import scipy
from FCM import cmeans
from sklearn.metrics import adjusted_rand_score
from MPCKmeans_SR_notabs import *
from MPCKmeansmf_SR_notabs import *
from copy import deepcopy
from concurrent.futures import ThreadPoolExecutor, as_completed
from random import sample
import itertools

import warnings
warnings.filterwarnings('error')


def prob_similar(i, j, u):
    """
    The function for calculating the probability of two instances will be similar after query
    i: The index of first instance
    j: The index of second instance
    u: The clustering membership vector obtained from Fuzzy c-means
    """
    return u[:,i].dot(u[:,j])

def Entropy(i, j, u):
    """
    The function for obtaining the entropy of a candidate instance pair
    i: The index of first instance
    j: The index of second instance
    u: The clustering membership vector obtained from Fuzzy c-means
    """

    prob = prob_similar(i,j,u)
    if prob == 0:
        prob = prob + 1/u.shape[1]
    if prob == 1:
        prob = prob - 1/u.shape[1]
    emc = - prob * np.log2(prob) - (1-prob) * np.log2(1-prob)

    return emc

def dist_x2y(x, y, num_col, cat_col, sr_col, cat_dist_mat):
    """
    The function for calculating the distances between two data points
    """
    p = x.shape[0]
    sr_col_a = [val for sr_col_sub in sr_col for val in sr_col_sub]

    if len(sr_col_a)>0:
        x2y = np.zeros(p + len(sr_col))
        x2y[sr_col_a] = x[sr_col_a] - y[sr_col_a]
        x2y[p:(p + len(sr_col))] = [sin(arccos(np.nanmin([pearsonr(x[sr_col_sub],  y[sr_col_sub])[0],1]))/2) for sr_col_sub in sr_col]
    else:
        x2y = np.zeros(p)
   
    x2y[num_col] = x[num_col] - y[num_col]
    if cat_dist_mat == None: # If cat_dist_mat is None, then use the Hamming distance
        x2y[cat_col] = [int(int(x[cat_col[i]]) != int(y[cat_col[i]])) for i in range(len(cat_col))]
    else:   
        x2y[cat_col] = [np.array(cat_dist_mat[i])[int(x[cat_col[i]]), int(y[cat_col[i]])] for i in range(len(cat_col))]
    return x2y

def DIV0(pair1, pair2, data, A, num_col, cat_col, sr_col, cat_dist_mat):
    """
    The function for calculating the dissimiarity between two query pairs
    """

    i, j = pair1
    k, l = pair2
    G = np.linalg.cholesky(A).T
    dist_i2j = dist_x2y(data[i], data[j], num_col, cat_col, sr_col, cat_dist_mat)
    dist_k2l = dist_x2y(data[k], data[l], num_col, cat_col, sr_col, cat_dist_mat)
    v1 = dot(G, dist_i2j)
    v2 = dot(G, dist_k2l)
    cs = np.clip(dot(v1, v2)/np.maximum((norm(v1) * norm(v2)), 1e-6), -1 + 1e-10, 1 - 1e-10)
    return sin(arccos(cs))   
    

def DIV(pair, B, data, A, num_col, cat_col, sr_col, cat_dist_mat):
    """
    The function for calculating the diversity between the candidate pair and the pairs that are already in the batch
    pair: The indeces of the candidate pair
    B: The collection of the pairs already in the query batch
    data: The numpy data array
    A: The distance matrix
    num_col: The list containing the indeces of numerical columns
    cat_col: The list containing the indeces of categorical columns
    sr_col: The list containing groups of self-reported columns (Each group is a sub-list)
    cat_dist_mat: The list of matries defining the distance between categorical values (One element for one variable). If None, the Hamming distance will be used
    """
    
    # The minimum dissimilarity
    div_min = 2

    for pair_B in B:
        div = DIV0(pair, pair_B, data, A, num_col, cat_col, sr_col, cat_dist_mat)
        if div < div_min:
            div_min = div
    return div_min

def find_query_batch(AP, size, u, data, A, num_col, cat_col, sr_col, cat_dist_mat):
    """
    AP (array): The index set of the anchor points obtained from Step I
    size (int): The size of the query batch 
    u: The clustering membership vector obtained from Fuzzy c-means
    data: The numpy data array
    A: The learned distance matrix
    num_col: The list containing the indeces of numerical columns
    cat_col: The list containing the indeces of categorical columns
    sr_col: The list containing groups of self-reported columns (Each group is a sub-list)
    cat_dist_mat: The list of matries defining the distance between categorical values (One element for one variable). If None, the Hamming distance will be used
    """

    if size < 1:
        raise ValueError('The budget set must be at least 1')

    B = []
    n = data.shape[0]
    NAP = np.setdiff1d(np.arange(n), AP)  # The index set of the points that are not anchor points

    # Precompute EMC matrix for NAP vs. AP in parallel
    with ThreadPoolExecutor() as executor:
        futures = {(nap, ap): executor.submit(Entropy, nap, ap, u)
                   for nap in NAP for ap in AP}
        EMC_mat = np.zeros((len(NAP), len(AP)))
        for (nap, ap), future in futures.items():
            i, j = NAP.tolist().index(nap), AP.tolist().index(ap)
            EMC_mat[i, j] = future.result()

    # Normalize the EMC matrix
    EMC_mat_s = (EMC_mat - EMC_mat.min()) / EMC_mat.max()

    # Find the first query pair with the largest EMC
    fp_index = np.unravel_index(np.argmax(EMC_mat_s), EMC_mat_s.shape)
    B.append((NAP[fp_index[0]], AP[fp_index[1]]))

    mask = np.zeros_like(EMC_mat, dtype=bool)
    mask[fp_index[0], fp_index[1]] = True

    while len(B) < size:
        # Precompute DIV matrix for all pairs (NAP, AP) not in B in parallel
        with ThreadPoolExecutor() as executor:
            futures = {(nap, ap): executor.submit(DIV, (nap, ap), B, data, A, num_col, cat_col, sr_col, cat_dist_mat)
                       for nap in NAP for ap in AP}
            div_mat = np.zeros((len(NAP), len(AP)))
            for (nap, ap), future in futures.items():
                i, j = NAP.tolist().index(nap), AP.tolist().index(ap)
                div_mat[i, j] = future.result()
        
        # Normalize the DIV matrix
        div_mat_s = (div_mat - div_mat.min()) / div_mat.max()

        # Aggregate the normalized EMC and DIV matrices
        agg_mat_s = EMC_mat_s + div_mat_s
        agg_mat_sm = np.ma.array(agg_mat_s, mask=mask)
        
        # Find the next query pair
        qp_index = np.unravel_index(np.argmax(agg_mat_sm), agg_mat_sm.shape)
        B.append((NAP[qp_index[0]], AP[qp_index[1]]))
        mask[qp_index[0], qp_index[1]] = True

    return B

def query2cons2(pairs, labels):
    ml = []  # Set of Must Link Constraints
    cl = []  # Set of Cannot Link Constraints
    for pair in pairs:
            if labels[pair[0]] == labels[pair[1]]:
                ml.append(pair)
            else:
                cl.append(pair)
    return ml, cl    

def find_random_batch(AP, size, n):
    """
    AP (array): The index set of the anchor points obtained from Step I
    size (int): The size of the query batch 
    u: The clustering membership vector obtained from Fuzzy c-means
    n: The size of the whole dataset
    """

    if size < 1:
        raise ValueError('The budget set must be at least 1')

    B = []
    NAP = np.setdiff1d(np.arange(n), AP)  # The index set of the points that are not anchor points

    all_pairs = list(itertools.product(AP, NAP))
    B = sample(all_pairs, size)

    return B
    




