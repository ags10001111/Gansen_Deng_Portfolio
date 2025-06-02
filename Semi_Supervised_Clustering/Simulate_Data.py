import numpy as np
from math import *
from scipy.stats import random_correlation
from sklearn.preprocessing import StandardScaler
from itertools import combinations
from random import sample
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn.cluster import KMeans
from kmodes.kprototypes import KPrototypes
from sklearn.metrics import adjusted_rand_score

from active_semi_clustering.semi_supervised.pairwise_constraints import PCKMeans, MPCKMeansMF, COPKMeans
from cobras_ts.cobras_kmeans import COBRAS_kmeans
from cobras_ts.querier.labelquerier import LabelQuerier

import metric_learn
import os
from sklearn.datasets import make_circles, make_moons, load_wine, load_digits, load_iris, load_breast_cancer, fetch_olivetti_faces, fetch_covtype
from active_semi_clustering.exceptions import EmptyClustersException
from active_semi_clustering.active.pairwise_constraints import ExampleOracle, MinMax, NPU, ExploreConsolidate

from Active_Query_Step1 import *
from MPCKmeans_SR_notabs import *
from MPCKmeans1 import *
from MPCKmeansmf_SR_notabs import *
import FCM
from Active_Query_Step2 import *
from Categorical_Distance import *
from NPU_SR import *


def generate_points_on_sphere(r, K, p):
    """
    Generate K points uniformly on a p-dimensional sphere with radius r.
    p >= 2
    return: list of list
    """
    def angle_to_coord(angles):
        x = []
        cum = 1 # record cumulitive product of sines
        for p in range(len(angles)):
            x += [r * cum * cos(angles[p])]
            cum  *= sin(angles[p])
        x.append( r * cum)
        return x
    
    return [angle_to_coord([i/K*pi] * (p-2) + [i/K*2*pi]) for i in range(K)]

def load_sphere(P1, L_P2, P3, P4, K, n_per_cluster, r = None, b=0.5, nCatLevels=3, sigma=1, seed=1, random_scale=False):
    """
    Generate mixed data with K clusters and n_per_cluster number of data points per cluster. The data is from Gaussian mixture model with K mixtures and the centers are uniformly located on a (P1+P2+P3)-dimensional sphere with radius r.
    
    P1: The number of continuous variables
    L_P2: The list consisting of the number of self-reported variables in each group
    P3: The number of categorical variables
    P4: The number of irrelevant variables
    K: The number of clusters
    r: The radius of the (P1+P2+P3)-dimensional sphere
    b: The standard deviation of the subjective bias
    n_per_cluster: The number of data points in each cluster
    nCatLevels: The number of levels for each categorical variable
    sigma: The standard deviation of each variable
    random_scale: A boolean variable indicating if the scale of each variable is random or not
    correlation: The correlation value between some variables within blocks of the covariance matrix
    """
    rng = np.random.RandomState(seed)
    P2 = sum(L_P2)
    P = P1 + P2 + P3
    X = np.zeros((1, P))
    y = np.zeros((1,))

    if r is None:
        r = 2 * sigma * sqrt(P/K)
        #r = 5
    
    # Generate K cluster centers on a sphere
    centers = generate_points_on_sphere(r, K, P)
    
    # Create a base covariance matrix with a block structure
    Sigma_base0 = sigma ** 2 * np.eye(P)  # Start with an identity matrix scaled by sigma^2
    Sigma_base = Sigma_base0.copy()
    eigenvalues = np.random.dirichlet(np.ones(P)) * P
    Sigma_base = random_correlation.rvs(eigenvalues)

    # Generate data for each cluster
    for k in range(K):
        mu0 = centers[k]
        #Sigma_base[P1:(P1+P2), P1:(P1+P2)] = random_correlation.rvs(eigenvalues)
        X0 = rng.multivariate_normal(mu0, Sigma_base, n_per_cluster)
        y0 = np.repeat(k, n_per_cluster)
        X = np.row_stack((X, X0))
        y = np.concatenate((y, y0))
    
    X = X[1:, :]  # Remove the initial row of zeros
    y = y[1:]
    
    # Apply scaling if random_scale is True
    if random_scale:
        scale = np.square(rng.normal(1, 1, P))
    else:
        scale = np.ones(P)
    X = transform(X, np.diag(scale))

    # Add subjective bias to self-reported (P2) variables
    if P2 > 0:
        k = P1
        i = 1
        for P2_0 in L_P2:
            #e = rng.normal((-1)**(i+1), b, n_per_cluster * K)
            #e = np.random.weibull(3, n_per_cluster * K) * 4
            e = np.random.gamma(2, 2, n_per_cluster * K) * ((-1)**(i+1))
            e_arr = np.tile(e, (P2_0, 1))
            noise = np.random.normal(0, 0.2, e_arr.shape)
            X[:, k:(k + P2_0)] = X[:, k:(k + P2_0)] + (e_arr + noise).T
            k = k + P2_0
            i = i + 1

    # Generate categorical variables from continuous data
    for p in range(P1 + P2, P1 + P2 + P3):
        X[:, p] = pd.cut(X[:, p], bins=nCatLevels, labels=range(nCatLevels))

    # Add irrelevant variables if P4 > 0
    if P4 > 0:
        X1, _, _, _ = load_opposite(P4, n_per_cluster * K, 5, seed+1, random_scale)
        X1 = rng.permutation(X1)
        X = np.column_stack((X, X1[:X.shape[0],:]))
    
    return X, y, scale

def load_opposite(P, N, mu=10, seed=1, random_scale=False):
    rng = np.random.RandomState(seed)
    X = np.zeros((1,P))
    y = np.zeros((1,))
    num_class = 0
    for p in range(P):
        sigma = 1
        mu0 = np.zeros((P,))
        mu0[p] = -mu
        mu1 = np.zeros((P,))
        mu1[p] = mu
        Sigma = sigma**2*np.diag(np.repeat(1, P))
        if N%2 == 0:
            X0 = rng.multivariate_normal(mu0, Sigma, int(N/2))
            X1 = rng.multivariate_normal(mu1, Sigma, int(N/2))
        else:
            X0 = rng.multivariate_normal(mu0, Sigma, int(N/2))
            X1 = rng.multivariate_normal(mu1, Sigma, int(N/2) + 1)
        y0 = np.repeat(p, N)
        num_class += 1
        X = np.row_stack((X,np.row_stack((X1, X0))))
        y = np.concatenate((y, y0)) 
    X = X[1:,:]
    y = y[1:]
    if random_scale:
        scale = np.square(rng.normal(1,1,P))
    else:
        scale = np.ones(P)
    X = transform(X, np.diag(scale))
    return X, y, num_class, scale

def ARI_clustering(X, y, C, method='kmeans', S=None, D=None):
    if method.lower() == 'kmeans':
        model = KMeans(n_clusters=C)
        model.fit(X)
    elif method.lower() == 'pckmeans':
        model = PCKMeans(n_clusters=C)
        model.fit(X, ml=list(S), cl=list(D))
    elif method.lower() == 'mpckmeansmf':
        # multiple full rank metric
        model = MPCKMeansMF(n_clusters=C)
        model.fit(X, ml=list(S), cl=list(D))    
    elif method.lower() == 'mpckmeans':
        model = MPCKMeans(n_clusters=C)
        model.fit(X, ml=list(S), cl=list(D))      
    return adjusted_rand_score(y, model.labels_)

def sim30(P1, L_P2, P3, P4, K, n_per_cluster = 20, b=0.5, nCatLevels=4, sigma=1, seed0=1, random_scale=False, H = False, B = 10, B1 = 6):
    r = None
    P2 = sum(L_P2)
    sim_df = pd.DataFrame()
    Prop_df = pd.DataFrame()
    Prop1_df = pd.DataFrame()
    MF_df = pd.DataFrame()
    MF1_df = pd.DataFrame()
    #kp_df = pd.DataFrame()
    for seed in range(seed0, seed0 + 30):
        sim_result = pd.DataFrame([sim0(P1, L_P2, P3, P4, K, n_per_cluster, r, b, nCatLevels, sigma, seed, random_scale, H, B, B1)])
        sim_df = pd.concat([sim_df, sim_result], ignore_index = True)
    
    Prop_df = sim_df.iloc[:, 0].reset_index(drop=True).to_frame(name="ARI")
    Prop_df["SR_Prop"] = P2 / (P1 + P2)
    Prop_df["P4"] = P4
    Prop_df["Method"] = "Diagonal"

    Prop1_df = sim_df.iloc[:, 1].reset_index(drop=True).to_frame(name="ARI")
    Prop1_df["SR_Prop"] = P2 / (P1 + P2)
    Prop1_df["P4"] = P4
    Prop1_df["Method"] = "Diagonal (NO COR)"

    MF_df = sim_df.iloc[:, 2].reset_index(drop=True).to_frame(name="ARI")
    MF_df["SR_Prop"] = P2 / (P1 + P2)
    MF_df["P4"] = P4
    MF_df["Method"] = "Dense"

    MF1_df = sim_df.iloc[:, 3].reset_index(drop=True).to_frame(name="ARI")
    MF1_df["SR_Prop"] = P2 / (P1 + P2)
    MF1_df["P4"] = P4
    MF1_df["Method"] = "Dense (NO COR)"

    #kp_df = sim_df.iloc[:, 4].reset_index(drop=True).to_frame(name="ARI")
    #kp_df["SR_Prop"] = P2 / (P1 + P2)
    #kp_df["P4"] = P4
    #kp_df["Method"] = "K-Prototype"

    result_df = pd.concat([Prop_df, Prop1_df, MF_df, MF1_df])
    return result_df

def sim_cat30(P1, L_P2, P3, P4, K, cat_dist=True, n_per_cluster=20, b=0.5, nCatLevels=4, sigma=1, seed0=1, random_scale=False, H=False, B=10, B1=8):
    sim_df = pd.DataFrame()
    Prop_df = pd.DataFrame()
    MF_df = pd.DataFrame()

    # Loop through 30 seeds
    for seed in range(seed0, seed0 + 30):
        sim_result = pd.DataFrame([sim_cat0(P1, L_P2, P3, P4, K, cat_dist, n_per_cluster, None, b, nCatLevels, sigma, seed, random_scale, H, B, B1)])
        sim_df = pd.concat([sim_df, sim_result], ignore_index=True)

    # Create proportion and method-specific dataframes
    Prop_df = sim_df.iloc[:, 0].reset_index(drop=True).to_frame(name="ARI")
    Prop_df["CAT_Prop"] = P3 / (P1 + P3)
    Prop_df["K"] = K
    Prop_df["CAT"] = cat_dist
    Prop_df["Method"] = "Diagonal"

    MF_df = sim_df.iloc[:, 1].reset_index(drop=True).to_frame(name="ARI")
    MF_df["CAT_Prop"] = P3 / (P1 + P3)
    MF_df["K"] = K
    MF_df["CAT"] = cat_dist
    MF_df["Method"] = "Dense"

    # Combine Prop_df and MF_df
    result_df = pd.concat([Prop_df, MF_df], ignore_index=True)

    return result_df

def sim30_AQ1_test(P1, L_P2, P3, P4, K, n_per_cluster = 20, b=0.5, nCatLevels=4, sigma=1, seed0=1, random_scale = False, B1 = 8):
    P2 = sum(L_P2)
    result_diag = pd.DataFrame()
    result_mf = pd.DataFrame()
    result_km = pd.DataFrame()
    
    if B1 >= K:
        r = K
        sim_df = pd.DataFrame()
        for seed in range(seed0, seed0 + 30):
            X = pd.DataFrame(load_sphere(P1, L_P2, P3, P4, K, n_per_cluster, r, b, nCatLevels, sigma, seed, random_scale)[0])
            y = load_sphere(P1, L_P2, P3, P4, K, n_per_cluster, r, b, nCatLevels, sigma, seed, random_scale)[1]
            P = P1 + P2 + P3 + P4
            num_col = list(range(0, P1)) + list(range(P1 + P2 + P3, P))
            k = P1
            sr_col = []
            for P2_0 in L_P2:
                sr_col0 = list(range(k, k + P2_0))
                k = k + P2_0
                sr_col.append(sr_col0)
            cat_col = list(range(P1+P2, P1 + P2 + P3))
            cat_dist_mat = cooccur2(X, cat_col)
            sim_result = pd.DataFrame([sim_AQ1_test(X, y, num_col, sr_col, cat_col, cat_dist_mat, B1)])
            sim_df = pd.concat([sim_df, sim_result], ignore_index = True)
        result_diag = pd.DataFrame([{'SR_Prop': P2/(P1+P2), 'K': K, 'b': b, 'B': B1, 'Prop': np.mean(sim_df.iloc[:,0] == K), 'Mean': np.mean(sim_df.iloc[:,0])}])
        result_mf = pd.DataFrame([{'SR_Prop': P2/(P1+P2), 'K': K, 'b': b, 'B': B1, 'Prop': np.mean(sim_df.iloc[:,1] == K), 'Mean': np.mean(sim_df.iloc[:,1])}])
        result_km = pd.DataFrame([{'SR_Prop': P2/(P1+P2), 'K': K, 'b': b, 'B': B1, 'Prop': np.mean(sim_df.iloc[:,2] == K), 'Mean': np.mean(sim_df.iloc[:,2])}])
    return result_diag, result_mf, result_km

def sim_con30(P1, L_P2, P4, K, n_per_cluster = 20, b=0.5, nCatLevels=4, sigma=1, seed0=1, random_scale=False, H = False, B = 10, B1 = 6):
    r = None
    P2 = sum(L_P2)
    sim_df = pd.DataFrame()
    for seed in range(seed0, seed0 + 30):
        sim_result = pd.DataFrame([sim_con0(P1, L_P2, P4, K, n_per_cluster, r, b, nCatLevels, sigma, seed, random_scale, H, B, B1)])
        sim_df = pd.concat([sim_df, sim_result], ignore_index = True)
    
    # Define method mappings
    method_names = [
        "Diagonal", "Dense", "PCKM", "MPCKM",
        "KM", "COBRAS"
    ]

    # Extract results for each method and add meta information
    result_dfs = []
    for idx, method in enumerate(method_names):
        method_df = sim_df.iloc[:, idx].reset_index(drop=True).to_frame(name="ARI")
        method_df["SR_Prop"] = P2 / (P1 + P2)  # SR proportion
        method_df["P4"] = P4  # Number of irrelevant features
        method_df["B"] = B  # Budget
        method_df["K"] = K  # Number of clusters
        method_df["Method"] = method  # Method name
        result_dfs.append(method_df)

    # Combine all method-specific DataFrames into a single DataFrame
    result_df = pd.concat(result_dfs, ignore_index=True)
    return result_df

def sim0(P1, L_P2, P3, P4, K, n_per_cluster = 20, r = None, b=0.5, nCatLevels=4, sigma=1, seed=1, random_scale=False, H = False, B = 10, B1 = 6):
    X = pd.DataFrame(load_sphere(P1, L_P2, P3, P4, K, n_per_cluster, r, b, nCatLevels, sigma, seed, random_scale)[0])
    y = load_sphere(P1, L_P2, P3, P4, K, n_per_cluster, r, b, nCatLevels, sigma, seed, random_scale)[1]
    P2 = sum(L_P2)
    P = P1 + P2 + P3 + P4
    num_col = list(range(0, P1)) + list(range(P1 + P2 + P3, P))
    k = P1
    sr_col = []
    for P2_0 in L_P2:
        sr_col0 = list(range(k, k + P2_0))
        k = k + P2_0
        sr_col.append(sr_col0)
    cat_col = list(range(P1+P2, P1 + P2 + P3))

    num_col1 = list(range(0, P1+P2)) + list(range(P1 + P2 + P3, P))

    #cat_dist_mat = cooccur2(X, cat_col)
    cat_dist_mat = None

    ARI_prop = sim_sr(X, y, num_col, sr_col, cat_col, cat_dist_mat, H, B, B1)

    ARI_prop1 = sim_sr(X, y, num_col1, [[]], cat_col, cat_dist_mat, H, B, B1)

    ARI_prop_MF = sim_mfsr(X, y, num_col, sr_col, cat_col, cat_dist_mat, H, B, B1)

    ARI_prop_MF1 = sim_mfsr(X, y, num_col, [[]], cat_col, cat_dist_mat, H, B, B1)

    #ARI_km = sim_kp(X, y, num_col, cat_col)

    ARI_dict = {'Proposed': ARI_prop, 'Proposed_NOCOR': ARI_prop1, 'Prop_MF': ARI_prop_MF, 'Prop_MF_NOCOR': ARI_prop_MF1}

    return ARI_dict

def sim_cat0(P1, L_P2, P3, P4, K, cat_dist, n_per_cluster = 20, r = None, b=0.5, nCatLevels=4, sigma=1, seed=1, random_scale=False, H = False, B = 10, B1 = 8):
    X = pd.DataFrame(load_sphere(P1, L_P2, P3, P4, K, n_per_cluster, r, b, nCatLevels, sigma, seed, random_scale)[0])
    y = load_sphere(P1, L_P2, P3, P4, K, n_per_cluster, r, b, nCatLevels, sigma, seed, random_scale)[1]
    P2 = sum(L_P2)
    P = P1 + P2 + P3 + P4
    num_col = list(range(0, P1)) + list(range(P1 + P2 + P3, P))
    k = P1
    sr_col = []
    for P2_0 in L_P2:
        sr_col0 = list(range(k, k + P2_0))
        k = k + P2_0
        sr_col.append(sr_col0)
    cat_col = list(range(P1+P2, P1 + P2 + P3))

    if cat_dist:
        cat_dist_mat = cooccur2(X, cat_col)

    else:
        cat_dist_mat = None

    ARI_prop = sim_sr(X, y, num_col, sr_col, cat_col, cat_dist_mat, H, B, B1)

    ARI_prop_MF = sim_mfsr(X, y, num_col, sr_col, cat_col, cat_dist_mat, H, B, B1)

    ARI_dict = {'Diagonal': ARI_prop, 'Dense': ARI_prop_MF}

    return ARI_dict

def sim_con0(P1, L_P2, P4, K, n_per_cluster = 20, r = None, b=0.5, nCatLevels=4, sigma=1, seed=1, random_scale=False, H = False, B = 10, B1 = 6):
    X = pd.DataFrame(load_sphere(P1, L_P2, 0, P4, K, n_per_cluster, r, b, nCatLevels, sigma, seed, random_scale)[0])
    y = load_sphere(P1, L_P2, 0, P4, K, n_per_cluster, r, b, nCatLevels, sigma, seed, random_scale)[1]
    P2 = sum(L_P2)
    P = P1 + P2 + P4
    num_col = list(range(0, P1)) + list(range(P1 + P2, P))
    k = P1
    sr_col = []
    for P2_0 in L_P2:
        sr_col0 = list(range(k, k + P2_0))
        k = k + P2_0
        sr_col.append(sr_col0)

    num_col1 = list(range(0, P1+P2)) + list(range(P1 + P2, P))

    ARI_prop = sim_sr(X, y, num_col, sr_col, [], None, H, B, B1)

    ARI_prop_MF = sim_mfsr(X, y, num_col, sr_col, [], None, H, B, B1)

    ARI_pckm = sim_pckm(X, y, list(range(0, P)))

    ARI_mpckm = sim_mpckm(X, y, list(range(0, P)))

    ARI_km = sim_km(X, y)

    ARI_COBRAS = sim_COBRAS(np.array(X), y, B)

    ARI_dict = {'Proposed': ARI_prop, 'Prop_MF': ARI_prop_MF, 'PCKM': ARI_pckm, 'MPCKM': ARI_mpckm, 'K-Means': ARI_km, 'COBRAS':ARI_COBRAS}

    return ARI_dict

def sim_lam0(P1, L_P2, P3, P4, K, n_per_cluster = 20, r = None, b=0.5, nCatLevels=4, sigma=1, seed=1, random_scale=False, H = True, B = 10, B1 = 8, lam = 1, lam_nd = 0):
    X = pd.DataFrame(load_sphere(P1, L_P2, P3, P4, K, n_per_cluster, r, b, nCatLevels, sigma, seed, random_scale)[0])
    y = load_sphere(P1, L_P2, P3, P4, K, n_per_cluster, r, b, nCatLevels, sigma, seed, random_scale)[1]
    P2 = sum(L_P2)
    P = P1 + P2 + P3 + P4
    num_col = list(range(0, P1)) + list(range(P1 + P2 + P3, P))
    k = P1
    sr_col = []
    for P2_0 in L_P2:
        sr_col0 = list(range(k, k + P2_0))
        k = k + P2_0
        sr_col.append(sr_col0)
    cat_col = list(range(P1+P2, P1 + P2 + P3))

    cat_dist_mat = cooccur2(X, cat_col)
    #cat_dist_mat = None

    ARI_prop = sim_sr(X, y, num_col, sr_col, cat_col, cat_dist_mat, H, B, B1, lam)

    ARI_prop_MF = sim_mfsr(X, y, num_col, sr_col, cat_col, cat_dist_mat, H, B, B1, lam, lam_nd)

    ARI_dict = {'Proposed': ARI_prop, 'Prop_MF': ARI_prop_MF}

    return ARI_dict

def sim_lam30(P1, L_P2, P3, P4, K, n_per_cluster = 20, b=0.5, nCatLevels=4, sigma=1, seed0=1, random_scale=False, H = False, B = 10, B1 = 8, lam = 1, lam_nd = 0):
    r = None
    P2 = sum(L_P2)
    sim_df = pd.DataFrame()
    Prop_df = pd.DataFrame()
    MF_df = pd.DataFrame()
    for seed in range(seed0, seed0 + 30):
        sim_result = pd.DataFrame([sim_lam0(P1, L_P2, P3, P4, K, n_per_cluster, r, b, nCatLevels, sigma, seed, random_scale, H, B, B1, lam, lam_nd)])
        sim_df = pd.concat([sim_df, sim_result], ignore_index = True)
    sim_df = pd.concat([sim_df, sim_df.describe().loc[["mean", "std"]]])
    sim_df = sim_df.round(3)
    Prop_df = pd.DataFrame([{'SR_Prop': P2/(P1+P2), 'P3' : P3, 'P4': P4, 'K': K, 'b': b, 'lam': lam, 'Method': 'Proposed', 'Mean': sim_df.iloc[-2,0], 'SD': sim_df.iloc[-1,0]}])
    MF_df = pd.DataFrame([{'SR_Prop': P2/(P1+P2), 'P3' : P3, 'P4': P4, 'K': K, 'b': b, 'lam': lam, 'Method': 'MF', 'Mean': sim_df.iloc[-2,1], 'SD': sim_df.iloc[-1,1]}])
    result_df = pd.concat([Prop_df, MF_df])
    return result_df

def sim_lam_nd0(P1, L_P2, P3, P4, K, n_per_cluster = 20, r = 3, b=0.5, nCatLevels=4, sigma=1, seed=1, random_scale=False, H = True, B = 10, B1 = 8, lam = 8.7, lam_nd = 0):
    X = pd.DataFrame(load_sphere(P1, L_P2, P3, P4, K, n_per_cluster, r, b, nCatLevels, sigma, seed, random_scale)[0])
    y = load_sphere(P1, L_P2, P3, P4, K, n_per_cluster, r, b, nCatLevels, sigma, seed, random_scale)[1]
    P2 = sum(L_P2)
    P = P1 + P2 + P3 + P4
    num_col = list(range(0, P1)) + list(range(P1 + P2 + P3, P))
    k = P1
    sr_col = []
    for P2_0 in L_P2:
        sr_col0 = list(range(k, k + P2_0))
        k = k + P2_0
        sr_col.append(sr_col0)
    cat_col = list(range(P1+P2, P1 + P2 + P3))

    cat_dist_mat = cooccur2(X, cat_col)
    #cat_dist_mat = None

    ARI_prop_MF = sim_mfsr(X, y, num_col, sr_col, cat_col, cat_dist_mat, H, B, B1, lam, lam_nd)

    ARI_dict = {'Prop_MF': ARI_prop_MF}

    return ARI_dict

def sim_lam_nd30(P1, L_P2, P3, P4, K, n_per_cluster = 20, b=0.5, nCatLevels=4, sigma=1, seed0=1, random_scale=False, H = False, B = 10, B1 = 8, lam = 8.7, lam_nd = 0):
    r = K
    P2 = sum(L_P2)
    sim_df = pd.DataFrame()
    Prop_df = pd.DataFrame()
    MF_df = pd.DataFrame()
    for seed in range(seed0, seed0 + 30):
        sim_result = pd.DataFrame([sim_lam_nd0(P1, L_P2, P3, P4, K, n_per_cluster, r, b, nCatLevels, sigma, seed, random_scale, H, B, B1, lam, lam_nd)])
        sim_df = pd.concat([sim_df, sim_result], ignore_index = True)
    sim_df = pd.concat([sim_df, sim_df.describe().loc[["mean", "std"]]])
    sim_df = sim_df.round(3)
    MF_df = pd.DataFrame([{'SR_Prop': P2/(P1+P2), 'P3' : P3, 'P4': P4, 'K': K, 'b': b, 'lam_nd': lam_nd, 'Method': 'MF', 'Mean': sim_df.iloc[-2,0], 'SD': sim_df.iloc[-1,0]}])
    return MF_df

def sim_sr(X, y, num_col, sr_col, cat_col, cat_dist_mat, H = False, B = 10, B1 = 8, lam = 4.1):
    sr_col_a = [val for sr_col_sub in sr_col for val in sr_col_sub]
    p2 = len(sr_col) if len(sr_col_a) > 0 else 0

    if len(cat_col)>0:
        model1 = KPrototypes(n_clusters=B1, n_init=200)
        model1.fit(np.array(X), categorical = cat_col)
        labels1 = model1.labels_

    else:
        model1 = KMeans(n_clusters=B1)
        model1.fit(np.array(X))
        labels1 = model1.labels_

    dist_mat = get_dist_mat(X, labels1, A = np.identity(len(X.columns)), num_col = num_col + sr_col_a, cat_col = cat_col, sr_col = [[]], cat_dist_mat = None)
    Anc_Ind = find_anchor_points(dist_mat, labels1)['Index']
    
    K1 = max(len(np.unique(y[Anc_Ind])),2)
    model2 = MPCKMeans_sr(n_clusters=K1, lam = lam)
    model2.fit(np.array(X), ml=query2cons1(Anc_Ind, y)[0], cl=query2cons1(Anc_Ind, y)[1], num_col = num_col, cat_col = cat_col, sr_col = sr_col, cat_dist_mat = cat_dist_mat, H = False)
    labels2 = model2.labels_
    dist_mat = get_dist_mat(X, labels2, A = model2.A_, num_col = num_col, cat_col = cat_col, sr_col = sr_col, cat_dist_mat = cat_dist_mat)
    Anc_Ind = find_anchor_points(dist_mat, labels2)['Index']

    _, u, _, _, _, _, _ = FCM.cmeans(
	np.array(X).T, K1, 2, A=model2.A_, num_col = num_col, cat_col = cat_col, sr_col = sr_col, cat_dist_mat = cat_dist_mat, error=0.005, maxiter=1000, init=None
)
    
    query_batch = find_query_batch(np.array(Anc_Ind), B, u, np.array(X), model2.A_, num_col, cat_col, sr_col, cat_dist_mat)

    ML = query2cons1(Anc_Ind, y)[0] + query2cons2(query_batch, y)[0]
    CL = query2cons1(Anc_Ind, y)[1] + query2cons2(query_batch, y)[1]

    model3 = MPCKMeans_sr(n_clusters=K1, lam = lam)
    model3.fit(np.array(X), ml=ML, cl=CL, num_col = num_col, cat_col = cat_col, sr_col = sr_col, cat_dist_mat = cat_dist_mat, H = H)
    labels3 = model3.labels_

    ARI = adjusted_rand_score(y, labels3)

    return ARI

def sim_mfsr(X, y, num_col, sr_col, cat_col, cat_dist_mat, H = False, B = 10, B1 = 8, lam = 11.1, lam_nd = 0):
    sr_col_a = [val for sr_col_sub in sr_col for val in sr_col_sub]
    p2 = len(sr_col) if len(sr_col_a) > 0 else 0
    if len(cat_col)>0:
        model1 = KPrototypes(n_clusters=B1, n_init=50)
        model1.fit(np.array(X), categorical = cat_col)
        labels1 = model1.labels_

    else:
        model1 = KMeans(n_clusters=B1)
        model1.fit(np.array(X))
        labels1 = model1.labels_

    dist_mat = get_dist_mat(X, labels1, A = np.identity(len(X.columns)), num_col = num_col + sr_col_a, cat_col = cat_col, sr_col = [[]], cat_dist_mat = None)
    Anc_Ind = find_anchor_points(dist_mat, labels1)['Index']
    
    K1 = max(len(np.unique(y[Anc_Ind])),2)
    model2 = MPCKMeansMF_sr(n_clusters=K1, lam = lam, lam_nd = lam_nd)
    model2.fit(np.array(X), ml=query2cons1(Anc_Ind, y)[0], cl=query2cons1(Anc_Ind, y)[1], num_col = num_col, cat_col = cat_col, sr_col = sr_col, cat_dist_mat = cat_dist_mat, H = False)
    labels2 = model2.labels_
    dist_mat = get_dist_mat(X, labels2, A = model2.A_, num_col = num_col, cat_col = cat_col, sr_col = sr_col, cat_dist_mat = cat_dist_mat)
    Anc_Ind = find_anchor_points(dist_mat, labels2)['Index']

    _, u, _, _, _, _, _ = FCM.cmeans(
	np.array(X).T, K1, 2, A=model2.A_, num_col = num_col, cat_col = cat_col, sr_col = sr_col, cat_dist_mat = cat_dist_mat, error=0.005, maxiter=1000, init=None
)
    
    query_batch = find_query_batch(np.array(Anc_Ind), B, u, np.array(X), model2.A_, num_col, cat_col, sr_col, cat_dist_mat)

    ML = query2cons1(Anc_Ind, y)[0] + query2cons2(query_batch, y)[0]
    CL = query2cons1(Anc_Ind, y)[1] + query2cons2(query_batch, y)[1]

    model3 = MPCKMeansMF_sr(n_clusters=K1, lam = lam, lam_nd = lam_nd)
    model3.fit(np.array(X), ml=ML, cl=CL, num_col = num_col, cat_col = cat_col, sr_col = sr_col, cat_dist_mat = cat_dist_mat, H = H)
    labels3 = model3.labels_

    ARI = adjusted_rand_score(y, labels3)

    return ARI

def sim_pckm(X, y, num_col, B = 10, B1 = 8):
    model1 = PCKMeans(n_clusters=B1)
    model1.fit(np.array(X), ml=[], cl=[])
    labels1 = model1.labels_

    dist_mat = get_dist_mat(X, labels1, A = np.identity(len(X.columns)), num_col = num_col, cat_col = [], sr_col = [[]], cat_dist_mat = None)
    Anc_Ind = find_anchor_points(dist_mat, labels1)['Index']
    
    K1 = max(len(np.unique(y[Anc_Ind])),2)
    model2 = PCKMeans(n_clusters=K1)
    model2.fit(np.array(X), ml=query2cons1(Anc_Ind, y)[0], cl=query2cons1(Anc_Ind, y)[1])
    labels2 = model2.labels_
    dist_mat = get_dist_mat(X, labels2, A = np.identity(len(X.columns)), num_col = num_col, cat_col = [], sr_col = [[]], cat_dist_mat = None)
    Anc_Ind = find_anchor_points(dist_mat, labels2)['Index']

    _, u, _, _, _, _, _ = FCM.cmeans(
	np.array(X).T, K1, 2, A=np.identity(len(X.columns)), num_col = num_col, cat_col = [], sr_col = [[]], cat_dist_mat = None, error=0.005, maxiter=1000, init=None
)
    
    query_batch = find_query_batch(np.array(Anc_Ind), B, u, np.array(X), np.identity(len(X.columns)), num_col, [], [[]], None)

    ML = query2cons1(Anc_Ind, y)[0] + query2cons2(query_batch, y)[0]
    CL = query2cons1(Anc_Ind, y)[1] + query2cons2(query_batch, y)[1]

    model3 = PCKMeans(n_clusters=K1)
    model3.fit(np.array(X), ml=ML, cl=CL)
    labels3 = model3.labels_

    ARI = adjusted_rand_score(y, labels3)

    return ARI

def sim_mpckm(X, y, num_col, B = 10, B1 = 8):
    model1 = MPCKMeans(n_clusters=B1)
    model1.fit(np.array(X), ml=[], cl=[])
    labels1 = model1.labels_

    dist_mat = get_dist_mat(X, labels1, A = np.identity(len(X.columns)), num_col = num_col, cat_col = [], sr_col = [[]], cat_dist_mat = None)
    Anc_Ind = find_anchor_points(dist_mat, labels1)['Index']
    
    K1 = max(len(np.unique(y[Anc_Ind])),2)
    model2 = MPCKMeans(n_clusters=K1)
    model2.fit(np.array(X), ml=query2cons1(Anc_Ind, y)[0], cl=query2cons1(Anc_Ind, y)[1])
    labels2 = model2.labels_
    dist_mat = get_dist_mat(X, labels2, A = model2.A_, num_col = num_col, cat_col = [], sr_col = [[]], cat_dist_mat = None)
    Anc_Ind = find_anchor_points(dist_mat, labels2)['Index']

    _, u, _, _, _, _, _ = FCM.cmeans(
	np.array(X).T, K1, 2, A=model2.A_, num_col = num_col, cat_col = [], sr_col = [[]], cat_dist_mat = None, error=0.005, maxiter=1000, init=None
)
    
    query_batch = find_query_batch(np.array(Anc_Ind), B, u, np.array(X), model2.A_, num_col, [], [[]], None)

    ML = query2cons1(Anc_Ind, y)[0] + query2cons2(query_batch, y)[0]
    CL = query2cons1(Anc_Ind, y)[1] + query2cons2(query_batch, y)[1]

    model3 = MPCKMeans(n_clusters=K1)
    model3.fit(np.array(X), ml=ML, cl=CL)
    labels3 = model3.labels_

    ARI = adjusted_rand_score(y, labels3)

    return ARI

def sim_kp(X, y, num_col, cat_col, B1 = 8):
    model1 = KPrototypes(n_clusters=B1)
    model1.fit(np.array(X), categorical = cat_col)
    labels1 = model1.labels_

    dist_mat = get_dist_mat(X, labels1, A = np.identity(len(X.columns)), num_col = num_col, cat_col = cat_col, sr_col = [[]], cat_dist_mat = None)
    Anc_Ind = find_anchor_points(dist_mat, labels1)['Index']
    
    K1 = max(len(np.unique(y[Anc_Ind])),2)
    model2 = KPrototypes(n_clusters=K1)
    model2.fit(np.array(X), categorical = cat_col)
    labels2 = model2.labels_

    ARI = adjusted_rand_score(y, labels2)

    return ARI

def sim_km(X, y, B1 = 8):
    model1 = KMeans(n_clusters=B1)
    model1.fit(np.array(X))
    labels1 = model1.labels_

    dist_mat = get_dist_mat(X, labels1, A = np.identity(len(X.columns)), num_col = list(range(len(X.columns))), cat_col = [], sr_col = [[]], cat_dist_mat = None)
    Anc_Ind = find_anchor_points(dist_mat, labels1)['Index']
    
    K1 = max(len(np.unique(y[Anc_Ind])),2)
    model2 = KMeans(n_clusters=K1)
    model2.fit(np.array(X))
    labels2 = model2.labels_

    ARI = adjusted_rand_score(y, labels2)

    return ARI

def sim_AQ1_test(X, y, num_col, sr_col, cat_col, cat_dist_mat, B1 = 8):
    sr_col_a = [val for sr_col_sub in sr_col for val in sr_col_sub]
    if False:
        model1 = MPCKMeans_sr(n_clusters=B1)
        model1.fit(np.array(X), ml=[], cl=[], num_col = num_col, cat_col = cat_col, sr_col = sr_col, cat_dist_mat = cat_dist_mat, H = False)
        labels1 = model1.labels_
        #centers1 = model1.cluster_centers_

        
        p2 = len(sr_col) if len(sr_col_a) > 0 else 0
        dist_mat = get_dist_mat(X, labels1, A= np.identity(len(X.columns) + p2), num_col = num_col, cat_col = cat_col, sr_col = sr_col, cat_dist_mat = cat_dist_mat)
        Anc_Ind = find_anchor_points(dist_mat, labels1)['Index']
        
        K1 = len(np.unique(y[Anc_Ind]))

        model2 = MPCKMeansMF_sr(n_clusters=B1)
        model2.fit(np.array(X), ml=[], cl=[], num_col = num_col, cat_col = cat_col, sr_col = sr_col, cat_dist_mat = cat_dist_mat, H = False)
        labels2 = model2.labels_
        #centers2 = model2.cluster_centers_

        dist_mat = get_dist_mat(X, labels2, A = np.identity(len(X.columns) + p2), num_col = num_col, cat_col = cat_col, sr_col = sr_col, cat_dist_mat = cat_dist_mat)
        Anc_Ind = find_anchor_points(dist_mat, labels2)['Index']
        
        K2 = len(np.unique(y[Anc_Ind]))

    K1 = 2
    K2 = 2
    model3 = KPrototypes(n_clusters=B1, n_init=200)
    model3.fit(np.array(X), categorical = cat_col)
    labels3 = model3.labels_
    #centers3 = model3.cluster_centers_

    dist_mat = get_dist_mat(X, labels3, A = np.identity(len(X.columns)), num_col = num_col + sr_col_a, cat_col = cat_col, sr_col = [[]], cat_dist_mat = None)
    Anc_Ind = find_anchor_points(dist_mat, labels3)['Index']

    K3 = len(np.unique(y[Anc_Ind]))

    return {'Diag': K1, 'MF': K2, 'KM': K3}

def sim_query30(P1, L_P2, P4, K, n_per_cluster=20, b=0.5, nCatLevels=4, sigma=1, seed0=1, random_scale=False, H=False, B=10):
    """
    Simulates clustering performance over 30 runs and aggregates results.

    Parameters:
        P1 (int): Number of numerical features.
        L_P2 (list): List of sizes of SR feature groups.
        P4 (int): Number of irrelevant features.
        K (int): Number of clusters.
        n_per_cluster (int): Number of points per cluster.
        nCatLevels (int): Number of levels in categorical variables.
        sigma (float): Noise level for data generation.
        seed0 (int): Starting seed for random number generator.
        random_scale (bool): Whether to randomize scaling.
        H (bool): Additional option for control.
        B (int): Total number of budgeted queries.
        B1 (int): Number of budgeted queries for must-link constraints.

    Returns:
        pd.DataFrame: Aggregated clustering results across 30 runs.
    """
    P2 = sum(L_P2)  # Total SR feature count
    B1 = B + 3
    sim_df = pd.DataFrame()  # DataFrame to store results across seeds

    # Run simulations for 30 seeds
    for seed in range(seed0, seed0 + 30):
        # Collect results from sim_query0
        sim_result = pd.DataFrame([sim_query0(P1, L_P2, P4, K, n_per_cluster, None, b, nCatLevels, sigma, seed, random_scale, H, B, B1)])
        sim_df = pd.concat([sim_df, sim_result], ignore_index=True)

    # Define method mappings
    method_names = [
        "Diagonal", "Diagonal (1R)", "Diagonal (2R)", "Diagonal (12R)",
        "Dense", "Dense (1R)", "Dense (2R)", "Dense (12R)",
        "MM_DIAG", "MM_DENSE", "NPU_DIAG", "NPU_DENSE"
    ]

    # Extract results for each method and add meta information
    result_dfs = []
    for idx, method in enumerate(method_names):
        method_df = sim_df.iloc[:, idx].reset_index(drop=True).to_frame(name="ARI")
        method_df["SR_Prop"] = P2 / (P1 + P2)  # SR proportion
        method_df["P4"] = P4  # Number of irrelevant features
        method_df["B"] = B  # Budget
        method_df["K"] = K  # Number of clusters
        method_df["Method"] = method  # Method name
        result_dfs.append(method_df)

    # Combine all method-specific DataFrames into a single DataFrame
    result_df = pd.concat(result_dfs, ignore_index=True)
    return result_df


def sim_query0(P1, L_P2, P4, K, n_per_cluster = 20, r = None, b=0.5, nCatLevels=4, sigma=1, seed=1, random_scale=False, H = False, B = 10, B1 = 6):
    X = load_sphere(P1, L_P2, 0, P4, K, n_per_cluster, r, b, nCatLevels, sigma, seed, random_scale)[0]
    y = load_sphere(P1, L_P2, 0, P4, K, n_per_cluster, r, b, nCatLevels, sigma, seed, random_scale)[1]
    P2 = sum(L_P2)
    P = P1 + P2 + P4
    num_col = list(range(0, P1)) + list(range(P1 + P2, P))
    k = P1
    sr_col = []
    for P2_0 in L_P2:
        sr_col0 = list(range(k, k + P2_0))
        k = k + P2_0
        sr_col.append(sr_col0)

    num_col1 = list(range(0, P1+P2)) + list(range(P1 + P2, P))

    NPU_SR = ARI_semi_active(X, y, K, B, 'sr', 'npu', num_col, [], sr_col)

    NPU_MFSR = ARI_semi_active(X, y, K, B, 'mfsr', 'npu', num_col, [], sr_col)

    ARI_prop = sim_sr_query(pd.DataFrame(X), y, num_col, sr_col, B, B1, False, False)

    ARI_prop_1r = sim_sr_query(pd.DataFrame(X), y, num_col, sr_col, B, B1, True, False)

    ARI_prop_2r = sim_sr_query(pd.DataFrame(X), y, num_col, sr_col, B, B1, False, True)

    ARI_prop_12r = sim_sr_query(pd.DataFrame(X), y, num_col, sr_col, B, B1, True, True)

    ARI_prop_MF = sim_mfsr_query(pd.DataFrame(X), y, num_col, sr_col, B, B1, False, False)

    ARI_prop_MF_1r = sim_mfsr_query(pd.DataFrame(X), y, num_col, sr_col, B, B1, True, False)

    ARI_prop_MF_2r = sim_mfsr_query(pd.DataFrame(X), y, num_col, sr_col, B, B1, False, True)

    ARI_prop_MF_12r = sim_mfsr_query(pd.DataFrame(X), y, num_col, sr_col, B, B1, True, True)

    MM_SR = ARI_semi_active(X, y, K, B, 'sr', 'minmax', num_col, [], sr_col)

    MM_MFSR = ARI_semi_active(X, y, K, B, 'mfsr', 'minmax', num_col, [], sr_col)

    #ARI_COBRAS = sim_COBRAS(X, y, B)

    ARI_dict = {'Diagonal': ARI_prop, 'Diagonal_1r': ARI_prop_1r, 'Diagonal_2r': ARI_prop_2r, 'Diagonal_12r': ARI_prop_12r,
    'Dense': ARI_prop_MF, 'Dense_1r': ARI_prop_MF_1r, 'Dense_2r': ARI_prop_MF_2r, 'Dense_12r': ARI_prop_MF_12r,
    'MM_DIAG': MM_SR, 'MM_DENSE': MM_MFSR, 'NPU_DIAG': NPU_SR, 'NPU_DENSE': NPU_MFSR
    }

    return ARI_dict

def sim_sr_query(X, y, num_col, sr_col, B = 10, B1 = 6, S1_Random = False, S2_Random = False):
    sr_col_a = [val for sr_col_sub in sr_col for val in sr_col_sub]
    p2 = len(sr_col) if len(sr_col_a) > 0 else 0

    if S1_Random:
        Anc_Ind = sample(range(X.shape[0]), B1)

    else:
        model1 = KMeans(n_clusters=B1)
        model1.fit(np.array(X))
        labels1 = model1.labels_

        dist_mat = get_dist_mat(X, labels1, A = np.identity(len(X.columns)), num_col = num_col + sr_col_a, cat_col = [], sr_col = [[]], cat_dist_mat = None)
        Anc_Ind = find_anchor_points(dist_mat, labels1)['Index']
    
    K1 = max(len(np.unique(y[Anc_Ind])),2)
    B1_used = query2cons1p_c(Anc_Ind, y, X, np.identity(len(X.columns)+p2), num_col, [], sr_col, None)
    ML = query2cons1(Anc_Ind, y)[0]
    CL = query2cons1(Anc_Ind, y)[1]

    if B - B1_used > 1:
        if S2_Random:
            query_batch = find_random_batch(Anc_Ind, B - B1_used, X.shape[0])
        else:
            model2 = MPCKMeans_sr(n_clusters=K1)
            model2.fit(np.array(X), ml=ML, cl=CL, num_col = num_col, cat_col = [], sr_col = sr_col, cat_dist_mat = None, H = False)
            labels2 = model2.labels_
            dist_mat = get_dist_mat(X, labels2, A = model2.A_, num_col = num_col, cat_col = [], sr_col = sr_col, cat_dist_mat = None)
            Anc_Ind = find_anchor_points(dist_mat, labels2)['Index']

            _, u, _, _, _, _, _ = FCM.cmeans(
            np.array(X).T, K1, 2, A=model2.A_, num_col = num_col, cat_col = [], sr_col = sr_col, cat_dist_mat = None, error=0.005, maxiter=1000, init=None
        )
            
            query_batch = find_query_batch(np.array(Anc_Ind), B - B1_used, u, np.array(X), model2.A_, num_col, [], sr_col, None)

        ML = query2cons1(Anc_Ind, y)[0] + query2cons2(query_batch, y)[0]
        CL = query2cons1(Anc_Ind, y)[1] + query2cons2(query_batch, y)[1]

    model3 = MPCKMeans_sr(n_clusters=K1)
    model3.fit(np.array(X), ml=ML, cl=CL, num_col = num_col, cat_col = [], sr_col = sr_col, cat_dist_mat = None, H = False)
    labels3 = model3.labels_

    ARI = adjusted_rand_score(y, labels3)

    return ARI

def sim_mfsr_query(X, y, num_col, sr_col, B = 10, B1 = 6, S1_Random = False, S2_Random = False):
    sr_col_a = [val for sr_col_sub in sr_col for val in sr_col_sub]
    p2 = len(sr_col) if len(sr_col_a) > 0 else 0

    if S1_Random:
        Anc_Ind = sample(range(X.shape[0]), B1)

    else:
        model1 = KMeans(n_clusters=B1)
        model1.fit(np.array(X))
        labels1 = model1.labels_

        dist_mat = get_dist_mat(X, labels1, A = np.identity(len(X.columns)), num_col = num_col + sr_col_a, cat_col = [], sr_col = [[]], cat_dist_mat = None)
        Anc_Ind = find_anchor_points(dist_mat, labels1)['Index']
    
    K1 = max(len(np.unique(y[Anc_Ind])),2)
    B1_used = query2cons1p_c(Anc_Ind, y, X, np.identity(len(X.columns)+p2), num_col, [], sr_col, None)

    ML = query2cons1(Anc_Ind, y)[0]
    CL = query2cons1(Anc_Ind, y)[1]

    if B - B1_used > 1:
        if S2_Random:
            query_batch = find_random_batch(Anc_Ind, B - B1_used, X.shape[0])
        else:
            model2 = MPCKMeansMF_sr(n_clusters=K1)
            model2.fit(np.array(X), ml=ML, cl=CL, num_col = num_col, cat_col = [], sr_col = sr_col, cat_dist_mat = None, H = False)
            labels2 = model2.labels_
            dist_mat = get_dist_mat(X, labels2, A = model2.A_, num_col = num_col, cat_col = [], sr_col = sr_col, cat_dist_mat = None)
            Anc_Ind = find_anchor_points(dist_mat, labels2)['Index']

            _, u, _, _, _, _, _ = FCM.cmeans(
            np.array(X).T, K1, 2, A=model2.A_, num_col = num_col, cat_col = [], sr_col = sr_col, cat_dist_mat = None, error=0.005, maxiter=1000, init=None
        )
            
            query_batch = find_query_batch(np.array(Anc_Ind), B - B1_used, u, np.array(X), model2.A_, num_col, [], sr_col, None)

        ML = query2cons1(Anc_Ind, y)[0] + query2cons2(query_batch, y)[0]
        CL = query2cons1(Anc_Ind, y)[1] + query2cons2(query_batch, y)[1]

    model3 = MPCKMeansMF_sr(n_clusters=K1)
    model3.fit(np.array(X), ml=ML, cl=CL, num_col = num_col, cat_col = [], sr_col = sr_col, cat_dist_mat = None, H = False)
    labels3 = model3.labels_

    ARI = adjusted_rand_score(y, labels3)

    return ARI

def sim_COBRAS(X, y, budget):
    clusterer = COBRAS_kmeans(X, LabelQuerier(y), budget)
    clustering, intermediate_clusterings, runtimes, ml, cl = clusterer.cluster()
    return adjusted_rand_score(y, clustering.construct_cluster_labeling())

def ARI_semi_active(X, y, K, nc, semi, active, num_col, cat_col = [], sr_col = []):
    oracle = ExampleOracle(y, max_queries_cnt=nc)
    
    if semi.lower() == 'pckmeans':
        clusterer = PCKMeans(n_clusters = K)
    elif semi.lower() == 'mpckmeans':
        clusterer = MPCKMeans(n_clusters = K)
    elif semi.lower() == 'copkmeans':
        clusterer = COPKMeans(n_clusters = K)
    elif semi.lower() == 'sr':
        clusterer = MPCKMeans_sr(n_clusters = K)
    elif semi.lower() == 'mfsr':
        clusterer = MPCKMeansMF_sr(n_clusters = K)
        
        
    if active.lower() == 'minmax':
        active_learner = MinMax(n_clusters=K)
        if (semi.lower() == 'sr') or (semi.lower() == 'mfsr'):
            active_learner.fit(X, oracle)
            pairwise_constraints = active_learner.pairwise_constraints_
            ML = [tuple(sublist) for sublist in pairwise_constraints[0]]
            CL = [tuple(sublist) for sublist in pairwise_constraints[1]]
            clusterer.fit(X, ml=ML, cl=CL, num_col = num_col, cat_col = cat_col, sr_col = sr_col)
        else:
            active_learner.fit(X, oracle)
            pairwise_constraints = active_learner.pairwise_constraints_
            ML = pairwise_constraints[0]
            CL = pairwise_constraints[1]
            clusterer.fit(X, ml=ML, cl=CL)
    
    elif active.lower() == 'npu':
        if (semi.lower() == 'sr') or (semi.lower() == 'mfsr'):
            active_learner = NPU_SR(clusterer=clusterer)
            active_learner.fit(X, oracle, num_col = num_col, cat_col = cat_col, sr_col = sr_col)
            pairwise_constraints = active_learner.pairwise_constraints_
            ML = [tuple(sublist) for sublist in pairwise_constraints[0]]
            CL = [tuple(sublist) for sublist in pairwise_constraints[1]]
            clusterer.fit(X, ml=ML, cl=CL, num_col = num_col, cat_col = cat_col, sr_col = sr_col)
        else:
            active_learner = NPU(clusterer=clusterer)
            active_learner.fit(X, oracle)
            pairwise_constraints = active_learner.pairwise_constraints_
            ML = pairwise_constraints[0]
            CL = pairwise_constraints[1]
            clusterer.fit(X, ml=ML, cl=CL)
    return adjusted_rand_score(y, clusterer.labels_)    

def cp_app(X, y, y_all, labels0, num_col, sr_col, cat_col, cat_dist_mat, B = 10, B1 = 6):
    sr_col_a = [val for sr_col_sub in sr_col for val in sr_col_sub]
    p2 = len(sr_col) if len(sr_col_a) > 0 else 0

    model1 = KMeans(n_clusters=B1)
    model1.fit(np.array(X))
    labels1 = model1.labels_

    dist_mat = get_dist_mat(X, labels1, A = np.identity(len(X.columns)), num_col = num_col + sr_col_a, 
                            cat_col = cat_col, sr_col = [[]], cat_dist_mat = cat_dist_mat)
    Anc_Ind = find_anchor_points(dist_mat, labels1)['Index']
    
    K1 = max(len(np.unique(y[Anc_Ind])),2)
    B1_used = query2cons1p_c(Anc_Ind, y, X, np.identity(len(X.columns)+p2), num_col, cat_col, sr_col, cat_dist_mat)

    ML = query2cons1(Anc_Ind, y)[0]
    CL = query2cons1(Anc_Ind, y)[1]

    if B - B1_used > 1:
        model2 = MPCKMeans_sr(n_clusters=K1)
        model2.fit(np.array(X), ml=ML, cl=CL, num_col = num_col, cat_col = cat_col, sr_col = sr_col, cat_dist_mat = cat_dist_mat, H = False)
        labels2 = model2.labels_
        dist_mat = get_dist_mat(X, labels2, A = model2.A_, num_col = num_col, cat_col = [], sr_col = sr_col, cat_dist_mat = None)
        Anc_Ind = find_anchor_points(dist_mat, labels2)['Index']

        _, u, _, _, _, _, _ = FCM.cmeans(
        np.array(X).T, K1, 2, A=model2.A_, num_col = num_col, cat_col = [], sr_col = sr_col, cat_dist_mat = None, error=0.005, maxiter=1000, init=None
    )
            
        query_batch_diag = find_query_batch(np.array(Anc_Ind), B - B1_used, u, np.array(X), model2.A_, num_col, [], sr_col, None)


        model2 = MPCKMeansMF_sr(n_clusters=K1)
        model2.fit(np.array(X), ml=ML, cl=CL, num_col = num_col, cat_col = cat_col, sr_col = sr_col, cat_dist_mat = cat_dist_mat, H = False)
        labels2 = model2.labels_
        dist_mat = get_dist_mat(X, labels2, A = model2.A_, num_col = num_col, cat_col = [], sr_col = sr_col, cat_dist_mat = None)
        Anc_Ind = find_anchor_points(dist_mat, labels2)['Index']

        _, u, _, _, _, _, _ = FCM.cmeans(
        np.array(X).T, K1, 2, A=model2.A_, num_col = num_col, cat_col = [], sr_col = sr_col, cat_dist_mat = None, error=0.005, maxiter=1000, init=None
    )
            
        query_batch_dense = find_query_batch(np.array(Anc_Ind), B - B1_used, u, np.array(X), model2.A_, num_col, [], sr_col, None)

        model2 = MPCKMeans(n_clusters=K1)
        model2.fit(np.array(X), ml=ML, cl=CL)
        labels2 = model2.labels_
        dist_mat = get_dist_mat(X, labels2, A = model2.A_, num_col = num_col + sr_col_a, cat_col = [], sr_col = [[]], cat_dist_mat = None)
        Anc_Ind = find_anchor_points(dist_mat, labels2)['Index']

        _, u, _, _, _, _, _ = FCM.cmeans(
        np.array(X).T, K1, 2, A=model2.A_, num_col = num_col + sr_col_a, cat_col = [], sr_col = [[]], cat_dist_mat = None, error=0.005, maxiter=1000, init=None
    )
            
        query_batch_mpckm = find_query_batch(np.array(Anc_Ind), B - B1_used, u, np.array(X), model2.A_, num_col + sr_col_a, [], [[]], None)

        model2 = PCKMeans(n_clusters=K1)
        model2.fit(np.array(X), ml=ML, cl=CL)
        labels2 = model2.labels_
        dist_mat = get_dist_mat(X, labels2, A = np.identity(len(X.columns)), num_col = num_col + sr_col_a, cat_col = [], sr_col = [[]], cat_dist_mat = None)
        Anc_Ind = find_anchor_points(dist_mat, labels2)['Index']

        _, u, _, _, _, _, _ = FCM.cmeans(
        np.array(X).T, K1, 2, A=np.identity(len(X.columns)), num_col = num_col + sr_col_a, cat_col = [], sr_col = [[]], cat_dist_mat = None, error=0.005, maxiter=1000, init=None
    )
            
        query_batch_pckm = find_query_batch(np.array(Anc_Ind), B - B1_used, u, np.array(X), np.identity(len(X.columns)), num_col + sr_col_a, [], [[]], None)

    ARI_df = pd.DataFrame()

    for b in range(1, B - B1_used):
        model_diag = MPCKMeans_sr(n_clusters=K1)
        model_dense = MPCKMeansMF_sr(n_clusters=K1)
        model_mpckm = MPCKMeans(n_clusters=K1)
        model_pckm = PCKMeans(n_clusters=K1)
        
        ML = query2cons1(Anc_Ind, y)[0] + query2cons2(query_batch_diag[:b], y)[0]
        CL = query2cons1(Anc_Ind, y)[1] + query2cons2(query_batch_diag[:b], y)[1]
        model_diag.fit(np.array(X), ml=ML, cl=CL, num_col = num_col, cat_col = [], sr_col = sr_col, cat_dist_mat = None, H = False)
        labels3 = model_diag.labels_ # Cluster labels of superinstaces
        labels = labels0.copy()
        for i in range(len(labels3)):
            labels[labels0 == i] =  labels3[i]

        ARI = adjusted_rand_score(y_all, labels)
        new_row = {'ARI': ARI, 'B': B1_used + b, 'Method': 'Diagonal'}
        ARI_df = pd.concat([ARI_df, pd.DataFrame([new_row])], ignore_index=True)

        ML = query2cons1(Anc_Ind, y)[0] + query2cons2(query_batch_dense[:b], y)[0]
        CL = query2cons1(Anc_Ind, y)[1] + query2cons2(query_batch_dense[:b], y)[1]
        model_dense.fit(np.array(X), ml=ML, cl=CL, num_col = num_col, cat_col = [], sr_col = sr_col, cat_dist_mat = None, H = False)
        labels3 = model_dense.labels_
        labels = labels0.copy()
        for i in range(len(labels3)):
            labels[labels0 == i] =  labels3[i]
        ARI = adjusted_rand_score(y_all, labels)
        new_row = {'ARI': ARI, 'B': B1_used + b, 'Method': 'Dense'}
        ARI_df = pd.concat([ARI_df, pd.DataFrame([new_row])], ignore_index=True)

        ML = query2cons1(Anc_Ind, y)[0] + query2cons2(query_batch_mpckm[:b], y)[0]
        CL = query2cons1(Anc_Ind, y)[1] + query2cons2(query_batch_mpckm[:b], y)[1]
        model_mpckm.fit(np.array(X), ml=ML, cl=CL)
        labels3 = model_mpckm.labels_
        labels = labels0.copy()
        for i in range(len(labels3)):
            labels[labels0 == i] =  labels3[i]
        ARI = adjusted_rand_score(y_all, labels)
        new_row = {'ARI': ARI, 'B': B1_used + b, 'Method': 'MPCKM'}
        ARI_df = pd.concat([ARI_df, pd.DataFrame([new_row])], ignore_index=True)

        ML = query2cons1(Anc_Ind, y)[0] + query2cons2(query_batch_pckm[:b], y)[0]
        CL = query2cons1(Anc_Ind, y)[1] + query2cons2(query_batch_pckm[:b], y)[1]
        model_pckm.fit(np.array(X), ml=ML, cl=CL)
        labels3 = model_pckm.labels_
        labels = labels0.copy()
        for i in range(len(labels3)):
            labels[labels0 == i] =  labels3[i]
        ARI = adjusted_rand_score(y_all, labels)
        new_row = {'ARI': ARI, 'B': B1_used + b, 'Method': 'PCKM'}
        ARI_df = pd.concat([ARI_df, pd.DataFrame([new_row])], ignore_index=True)

        clusterer = COBRAS_kmeans(np.array(X), LabelQuerier(y), b + B1_used)
        clustering, intermediate_clusterings, runtimes, ml, cl = clusterer.cluster()
        labels = labels0.copy()
        for i in range(len(labels3)):
            labels[labels0 == i] =  clustering.construct_cluster_labeling()[i]
        ARI = adjusted_rand_score(y_all, labels)

        new_row = {'ARI': ARI, 'B': B1_used + b, 'Method': 'COBRAS'}
        ARI_df = pd.concat([ARI_df, pd.DataFrame([new_row])], ignore_index=True)

        ARI_df.to_csv('CP_APP.csv')

    return ARI_df

def ga_app(X, y, y_all, labels0, num_col, sr_col, cat_col, cat_dist_mat, B = 10, B1 = 6):
    sr_col_a = [val for sr_col_sub in sr_col for val in sr_col_sub]
    p2 = len(sr_col) if len(sr_col_a) > 0 else 0

    model1 = KMeans(n_clusters=B1)
    model1.fit(np.array(X))
    labels1 = model1.labels_

    dist_mat = get_dist_mat(X, labels1, A = np.identity(len(X.columns)), num_col = num_col + sr_col_a, 
                            cat_col = cat_col, sr_col = [[]], cat_dist_mat = cat_dist_mat)
    Anc_Ind = find_anchor_points(dist_mat, labels1)['Index']
    
    K1 = max(len(np.unique(y[Anc_Ind])),2)
    B1_used = query2cons1p_c(Anc_Ind, y, X, np.identity(len(X.columns)+p2), num_col, cat_col, sr_col, cat_dist_mat)

    ML = query2cons1(Anc_Ind, y)[0]
    CL = query2cons1(Anc_Ind, y)[1]

    if B - B1_used > 1:
        model2 = MPCKMeans_sr(n_clusters=K1)
        model2.fit(np.array(X), ml=ML, cl=CL, num_col = num_col, cat_col = cat_col, sr_col = sr_col, cat_dist_mat = cat_dist_mat, H = False)
        labels2 = model2.labels_
        dist_mat = get_dist_mat(X, labels2, A = model2.A_, num_col = num_col, cat_col = [], sr_col = sr_col, cat_dist_mat = None)
        Anc_Ind = find_anchor_points(dist_mat, labels2)['Index']

        _, u, _, _, _, _, _ = FCM.cmeans(
        np.array(X).T, K1, 2, A=model2.A_, num_col = num_col, cat_col = [], sr_col = sr_col, cat_dist_mat = None, error=0.005, maxiter=1000, init=None
    )
            
        query_batch_diag = find_query_batch(np.array(Anc_Ind), B - B1_used, u, np.array(X), model2.A_, num_col, [], sr_col, None)


        model2 = MPCKMeansMF_sr(n_clusters=K1)
        model2.fit(np.array(X), ml=ML, cl=CL, num_col = num_col, cat_col = cat_col, sr_col = sr_col, cat_dist_mat = cat_dist_mat, H = False)
        labels2 = model2.labels_
        dist_mat = get_dist_mat(X, labels2, A = model2.A_, num_col = num_col, cat_col = [], sr_col = sr_col, cat_dist_mat = None)
        Anc_Ind = find_anchor_points(dist_mat, labels2)['Index']

        _, u, _, _, _, _, _ = FCM.cmeans(
        np.array(X).T, K1, 2, A=model2.A_, num_col = num_col, cat_col = [], sr_col = sr_col, cat_dist_mat = None, error=0.005, maxiter=1000, init=None
    )
            
        query_batch_dense = find_query_batch(np.array(Anc_Ind), B - B1_used, u, np.array(X), model2.A_, num_col, [], sr_col, None)

        model2 = MPCKMeans(n_clusters=K1)
        model2.fit(np.array(X), ml=ML, cl=CL)
        labels2 = model2.labels_
        dist_mat = get_dist_mat(X, labels2, A = model2.A_, num_col = num_col + sr_col_a, cat_col = [], sr_col = [[]], cat_dist_mat = None)
        Anc_Ind = find_anchor_points(dist_mat, labels2)['Index']

        _, u, _, _, _, _, _ = FCM.cmeans(
        np.array(X).T, K1, 2, A=model2.A_, num_col = num_col + sr_col_a, cat_col = [], sr_col = [[]], cat_dist_mat = None, error=0.005, maxiter=1000, init=None
    )
            
        query_batch_mpckm = find_query_batch(np.array(Anc_Ind), B - B1_used, u, np.array(X), model2.A_, num_col + sr_col_a, [], [[]], None)

        model2 = PCKMeans(n_clusters=K1)
        model2.fit(np.array(X), ml=ML, cl=CL)
        labels2 = model2.labels_
        dist_mat = get_dist_mat(X, labels2, A = np.identity(len(X.columns)), num_col = num_col + sr_col_a, cat_col = [], sr_col = [[]], cat_dist_mat = None)
        Anc_Ind = find_anchor_points(dist_mat, labels2)['Index']

        _, u, _, _, _, _, _ = FCM.cmeans(
        np.array(X).T, K1, 2, A=np.identity(len(X.columns)), num_col = num_col + sr_col_a, cat_col = [], sr_col = [[]], cat_dist_mat = None, error=0.005, maxiter=1000, init=None
    )
            
        query_batch_pckm = find_query_batch(np.array(Anc_Ind), B - B1_used, u, np.array(X), np.identity(len(X.columns)), num_col + sr_col_a, [], [[]], None)

    ARI_df = pd.DataFrame()

    for b in range(1, B - B1_used):
        model_diag = MPCKMeans_sr(n_clusters=K1)
        model_dense = MPCKMeansMF_sr(n_clusters=K1)
        model_mpckm = MPCKMeans(n_clusters=K1)
        model_pckm = PCKMeans(n_clusters=K1)
        
        ML = query2cons1(Anc_Ind, y)[0] + query2cons2(query_batch_diag[:b], y)[0]
        CL = query2cons1(Anc_Ind, y)[1] + query2cons2(query_batch_diag[:b], y)[1]
        model_diag.fit(np.array(X), ml=ML, cl=CL, num_col = num_col, cat_col = [], sr_col = sr_col, cat_dist_mat = None, H = False)
        labels3 = model_diag.labels_ # Cluster labels of superinstaces
        labels = labels0.copy()
        for i in range(len(labels3)):
            labels[labels0 == i] =  labels3[i]

        ARI = adjusted_rand_score(y_all, labels)
        new_row = {'ARI': ARI, 'B': B1_used + b, 'Method': 'Diagonal'}
        ARI_df = pd.concat([ARI_df, pd.DataFrame([new_row])], ignore_index=True)

        ML = query2cons1(Anc_Ind, y)[0] + query2cons2(query_batch_dense[:b], y)[0]
        CL = query2cons1(Anc_Ind, y)[1] + query2cons2(query_batch_dense[:b], y)[1]
        model_dense.fit(np.array(X), ml=ML, cl=CL, num_col = num_col, cat_col = [], sr_col = sr_col, cat_dist_mat = None, H = False)
        labels3 = model_dense.labels_
        labels = labels0.copy()
        for i in range(len(labels3)):
            labels[labels0 == i] =  labels3[i]
        ARI = adjusted_rand_score(y_all, labels)
        new_row = {'ARI': ARI, 'B': B1_used + b, 'Method': 'Dense'}
        ARI_df = pd.concat([ARI_df, pd.DataFrame([new_row])], ignore_index=True)

        ML = query2cons1(Anc_Ind, y)[0] + query2cons2(query_batch_mpckm[:b], y)[0]
        CL = query2cons1(Anc_Ind, y)[1] + query2cons2(query_batch_mpckm[:b], y)[1]
        model_mpckm.fit(np.array(X), ml=ML, cl=CL)
        labels3 = model_mpckm.labels_
        labels = labels0.copy()
        for i in range(len(labels3)):
            labels[labels0 == i] =  labels3[i]
        ARI = adjusted_rand_score(y_all, labels)
        new_row = {'ARI': ARI, 'B': B1_used + b, 'Method': 'MPCKM'}
        ARI_df = pd.concat([ARI_df, pd.DataFrame([new_row])], ignore_index=True)

        ML = query2cons1(Anc_Ind, y)[0] + query2cons2(query_batch_pckm[:b], y)[0]
        CL = query2cons1(Anc_Ind, y)[1] + query2cons2(query_batch_pckm[:b], y)[1]
        model_pckm.fit(np.array(X), ml=ML, cl=CL)
        labels3 = model_pckm.labels_
        labels = labels0.copy()
        for i in range(len(labels3)):
            labels[labels0 == i] =  labels3[i]
        ARI = adjusted_rand_score(y_all, labels)
        new_row = {'ARI': ARI, 'B': B1_used + b, 'Method': 'PCKM'}
        ARI_df = pd.concat([ARI_df, pd.DataFrame([new_row])], ignore_index=True)

        clusterer = COBRAS_kmeans(np.array(X), LabelQuerier(y), b + B1_used)
        clustering, intermediate_clusterings, runtimes, ml, cl = clusterer.cluster()
        labels = labels0.copy()
        for i in range(len(labels3)):
            labels[labels0 == i] =  clustering.construct_cluster_labeling()[i]
        ARI = adjusted_rand_score(y_all, labels)

        new_row = {'ARI': ARI, 'B': B1_used + b, 'Method': 'COBRAS'}
        ARI_df = pd.concat([ARI_df, pd.DataFrame([new_row])], ignore_index=True)

        ARI_df.to_csv('GA_APP.csv')

    return ARI_df


def transform(X, A):
    return np.dot(X, mat_sqrt(A))

def mat_sqrt(A):
    """
    return V st A=V*V'
    """
    w, v = np.linalg.eigh(A)
    eps = 1e-8
    w = np.array([x if x>eps else 0 for x in w])
    assert np.all(w>=0), "A is not PSD"
    return np.dot(v, np.sqrt(np.diag(np.abs(w)))).dot(v.T)