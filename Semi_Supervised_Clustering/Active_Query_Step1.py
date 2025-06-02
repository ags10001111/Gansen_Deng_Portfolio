"""
@author: Gansen Deng
"""
import numpy as np
import pandas as pd
from sklearn.metrics import pairwise_distances
import scipy
from scipy.stats import pearsonr
from numpy import dot, arccos, sin

def find_anchor_points(dist_mat, labels):
    ''' 
    Find the set of anchor points based on the pairwise disance matric of each cluster
    dist_mat: The list of pairwise disance matrices of each cluster
    labels: The label vector of each point
    '''
    
    K = len(dist_mat)
    Anc = [-1] * K
    Anc_den = [-1] * K

    for i in range(len(dist_mat)):
        ni = np.sum(np.array(labels) == i)
        den = dist_mat[i].sum(axis = 1)/ni
        Anc[i] = np.where(np.array(labels) == i)[0][np.argmin(den)]
        Anc_den[i] = np.min(den)

    return {'Index': Anc, 'Density':Anc_den}

def find_anchor_center(Anc, X, A, num_col, cat_col, sr_col, cat_dist_mat):
    ''' 
    Find the center of the anchor point set
    Anc: The indices of anchor points to query
    '''
    pdm = pairwise_distances(X.iloc[Anc], metric= SR_dist, A = A, num_col = num_col, cat_col = cat_col, 
                                 sr_col = sr_col, cat_dist_mat = cat_dist_mat)
    den = pdm.sum(axis = 1)
    ctr_idx = np.argmin(den)

    return ctr_idx

def SR_dist(x, y, A, num_col, cat_col, sr_col, cat_dist_mat):
    """
    The function for calculting the distance between two data points
    x: The first data point array
    y: The second data point array
    num_col: The list containing the indeces of numerical columns
    cat_col: The list containing the indeces of categorical columns
    sr_col: The list containing groups of self-reported columns (Each group is a sub-list)
    cat_dist_mat: The list of matries defining the distance between categorical values (One element for one variable). If None, the Hamming distance will be used
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
    dist_xy = np.dot(np.dot(x2y, A), x2y)
    return dist_xy

def dist_xu(x, u, A, num_col, cat_col, sr_col, cat_dist_mat):
        """
        The function for calculating the distance between a data point and a cluster center, i.e., (x - mu)^T A (x - mu)
        """

        p = x.shape[0]
        
        # Extract numerical and specific range columns for u
        u_num = np.array([u[i] for i in num_col])
        u_cat = [u[i] for i in cat_col]
        sr_col_a = [val for sr_col_sub in sr_col for val in sr_col_sub]

        if len(sr_col_a)>0:
            # Initialize x2u array
            x2u = np.zeros(p + len(sr_col))
            u_sr = np.array([u[i] for i in sr_col_a])
            x2u[sr_col_a] = abs(x[sr_col_a] - u_sr)
            
            # Calculate Pearson correlation and arccos
            x2u[p:(p + len(sr_col))] = [sin(arccos(np.nanmin([pearsonr(x[sr_col_sub],  np.array([u[i] for i in sr_col_sub]))[0],0]))/2) for sr_col_sub in sr_col]
        else:
            x2u = np.zeros(p)
        
        # Calculate differences for numerical and specific range columns
        x2u[num_col] = abs(x[num_col] - u_num)
        
        # Calculate differences for categorical columns
        if cat_dist_mat is None:  # Use Hamming distance if cat_dist_mat is None
            x2u[cat_col] = np.array([sum(int(x[cat_col[i]]) != int(k) * v for k, v in u_cat[i].items()) for i in range(len(cat_col))])
        else:
            x2u[cat_col] = np.array([sum(np.array(cat_dist_mat[i])[int(x[cat_col[i]]), int(k)] * v for k, v in u_cat[i].items()) for i in range(len(cat_col))])
        
        # Compute the final distance
        dist_xu = np.dot(np.dot(x2u, A), x2u)
        return np.sqrt(dist_xu)

def get_dist_mat(X, labels, A, num_col, cat_col, sr_col, cat_dist_mat):
    """
    The function for getting the pairwise distance matrix
    """
    K = len(np.unique(labels))
    dist_mat_list = []
    for i in range(K):
        pdm = pairwise_distances(X.loc[np.array(labels) == i], metric= SR_dist, A = A, num_col = num_col, cat_col = cat_col, 
                                 sr_col = sr_col, cat_dist_mat = cat_dist_mat)
        dist_mat_list.append(pdm)
    return dist_mat_list

def get_dist_mat1(X, labels, clu_centers, A, num_col, cat_col, sr_col, cat_dist_mat):
    """
    The function for getting the pairwise distance matrix
    """
    K = len(np.unique(labels))
    dist_mat_list = []
    for i in range(K):
        ni = np.sum(np.array(labels) == i)
        pdm = np.array([
                dist_xu(X.loc[np.array(labels) == i].iloc[j], clu_centers[i], A, num_col, cat_col, sr_col, cat_dist_mat)
                for j in range(ni)
            ])
        dist_mat_list.append(pdm)
    return dist_mat_list

def query2cons1(anc, labels):
    ml = []  # Set of Must Link Constraints
    cl = []  # Set of Cannot Link Constraints
    for i, a in enumerate(anc):
        for j in range(i+1, len(anc)):
            if labels[a] == labels[anc[j]]:
                ml.append((a, anc[j]))
            else:
                cl.append((a, anc[j]))
    return ml, cl

def query2cons1p(anc, labels, X, A, num_col, cat_col, sr_col, cat_dist_mat):
    ml = []  # Set of Must Link Constraints
    cl = []  # Set of Cannot Link Constraints
    Q = 0 # Record the number of queries
    while len(anc) > 1:
        ctr_idx = find_anchor_center(anc, X, A, num_col, cat_col, sr_col, cat_dist_mat)
        ctr_point = anc[ctr_idx]  # The anchor point at the center
        remaining_anchors = []  # List to hold remaining anchors after processing

        for j in range(len(anc)):
            if j != ctr_idx:
                if labels[ctr_point] == labels[anc[j]]:
                    ml.append((ctr_point, anc[j]))
                else:
                    cl.append((ctr_point, anc[j]))
                    remaining_anchors.append(anc[j])
                Q += 1

        anc = remaining_anchors
    return ml, cl, Q

def query2cons1p_c(anc, labels, X, A, num_col, cat_col, sr_col, cat_dist_mat):
    """
    Optimized function to generate Must-Link (ML) and Cannot-Link (CL) constraints
    with reduced comparisons using a center-based approach.
    
    Parameters:
    - anc: List of anchor indices.
    - labels: List or array of labels corresponding to each index.
    - X: Data matrix.
    - A: Metric matrix.
    - num_col, cat_col, sr_col: Columns for numeric, categorical, and specific SR variables.
    - cat_dist_mat: Precomputed categorical distance matrix.
    
    Returns:
    - ml: List of Must-Link constraints as pairs of indices.
    - cl: List of Cannot-Link constraints as pairs of indices.
    - Q: Total number of queries performed.
    """
    Q = 0    # Record the number of queries

    while len(anc) > 1:
        # Find the center anchor point
        ctr_idx = find_anchor_center(anc, X, A, num_col, cat_col, sr_col, cat_dist_mat)
        ctr_point = anc.pop(ctr_idx)  # Remove and retrieve the center anchor point

        # Partition anchors based on whether they Must-Link or Cannot-Link
        same_label = [x for x in anc if labels[ctr_point] == labels[x]]
        diff_label = [x for x in anc if labels[ctr_point] != labels[x]]

        # Add constraints and update queries
        Q += len(anc)

        # Update remaining anchors for the next iteration
        anc = diff_label

    return Q
