import numpy as np
import gurobipy as gp
from gurobipy import GRB
import matplotlib.ticker as ticker
import plotly.express as px
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import pandas as pd
pd.set_option('display.float_format', lambda x: '%.3f' % x)


class InterpretableClusteringOptimalTree:
    
    def __init__(self, max_depth=2, min_leaf_samples=1, min_cluster_number=2, sil_lb=-1, sil_ub=1):
        self.max_depth = max_depth
        self.min_leaf_samples = min_leaf_samples
        self.min_cluster_number = min_cluster_number
        self.sil_lb = sil_lb
        self.sil_ub = sil_ub

        # index for all nodes in a tree
        self.n_ind = [i + 1 for i in range(2 ** (self.max_depth + 1) - 1)]
        # index for branch nodes
        self.b_ind = self.n_ind[:-2**self.max_depth]
        # index for leaf nodes
        self.l_ind = self.n_ind[-2**self.max_depth:]
        # index left nodes
        self.left_nodes = self.n_ind[1::2]
        
    def fit(self, x, dist):
        self.n = x.shape[0]
        self.p = x.shape[1]
        self.dist = dist
        # scale data
        if np.min(x) >= 0:
            if np.max(x) <= 1:
                x = x
               
            else:
                self.scales = np.max(x, axis=0)
                self.scales[self.scales == 0] = 1
                x = x / self.scales
        else:
            x = (x - np.min(x)) / (np.max(x) - np.min(x))
            
        # solve MIP
        m, a, b, c, d, l, z, M, q, r, s, K, w, eps = self._constructMIO(x)
        m.setParam('NonConvex', 2)
        m.setParam('Timelimit', 2000)

        if self.max_depth >= 3:
            m.setParam('MIPFocus', 3)
            m.setParam('NoRelHeurTime', 200)
            m.setParam('MIPGap', 0.05)
            m.setParam('Cuts', 0)
            m.setParam('GomoryPasses', 1)
            
        m.optimize()
        m.write('model.lp')

        # get parameters
        self.a_ = {ind:a[ind].x for ind in a}
        self.b_ = {ind:b[ind].x for ind in b}
        self.c_ = {ind:c[ind].x for ind in c}
        self.d_ = {ind:d[ind].x for ind in d}
        self.l_ = {ind:l[ind].x for ind in l}
        self.z_ = {ind:z[ind].x for ind in z}
        self.M_ = {ind:M[ind].x for ind in M}
        self.q_ = {ind:q[ind].x for ind in q}
        self.r_ = {ind:r[ind].x for ind in r}
        self.s_ = {ind:s[ind].x for ind in s}
        self.K_ = {ind:K[ind].x for ind in K}
        self.w_ = {ind:w[ind].x for ind in w}
        self.eps_ = eps
         
    def _constructMIO(self, x):
        """MIO formulation for ICOT"""
    
        # compute all pairwise distances
        dist = self.dist

        # create a model
        m = gp.Model('m')

        # create variables
        a = m.addVars(self.p, self.b_ind, vtype=GRB.BINARY, name='a') # splitting feature
        b = m.addVars(self.b_ind, lb=0, ub=1, vtype=GRB.CONTINUOUS, name='b') # splitting threshold
        c = m.addVars(self.n, self.l_ind, lb=0, ub=np.sqrt(self.p), vtype=GRB.CONTINUOUS, name='c') # average distance of i from cluster t
        d = m.addVars(self.b_ind, vtype=GRB.BINARY, name='d') # splitting option
        z = m.addVars(self.n, self.n_ind, vtype=GRB.BINARY, name='z') # leaf node assignment
        l = m.addVars(self.l_ind, vtype=GRB.BINARY, name='l') # leaf node activation
        r = m.addVars(self.n, lb=0, ub=np.max(dist), vtype=GRB.CONTINUOUS, name='r') # average distance to all points in own cluster
        q = m.addVars(self.n, lb=0, ub=np.max(dist), vtype=GRB.CONTINUOUS, name='q') # average distance to all points in nearest cluster
        M = m.addVars(self.n, lb=0, ub=np.max(dist), vtype=GRB.CONTINUOUS, name='M') # max(r,q)
        s = m.addVars(self.n, lb=self.sil_lb, ub=self.sil_ub, vtype=GRB.CONTINUOUS, name='s') # silhoutte score
        K = m.addVars(self.l_ind, lb=0, ub=self.n-1, vtype=GRB.INTEGER, name='K') # number of points in cluster
        w = m.addVars(self.n, self.l_ind, vtype=GRB.BINARY, name='w') # auxillary variable for the second minimum

        # objective function
        m.setObjective(s.sum() / self.n, GRB.MAXIMIZE)
        
        # constraints
        m.addConstrs(s[i] * M[i] == (q[i] - r[i])  for i in range(self.n)) # 1
        m.addConstrs(q[i] <= M[i] for i in range(self.n)) # 2
        m.addConstrs(r[i] <= M[i] for i in range(self.n)) # 3
        m.addConstrs(gp.quicksum(c[i,t] * z[i,t] for t in self.l_ind) == r[i] for i in range(self.n)) # 7
        m.addConstrs(K[t] == gp.quicksum(z[i,t] for i in range(self.n)) for t in self.l_ind) # 9
        m.addConstrs(gp.quicksum(dist[i,j] * z[j,t] for j in range(self.n)) == c[i,t] * (K[t]) for t in self.l_ind for i in range(self.n)) # 8
        m.addConstrs(d[t] == gp.quicksum(a[j,t] for j in range(self.p)) for t in self.b_ind) # 10
        m.addConstrs(b[t] <= d[t] for t in self.b_ind) # 11
        m.addConstrs(d[t] <= d[t//2] for t in self.b_ind if t != 1) # 12
        m.addConstrs(z[i,t] <= l[t] for t in self.l_ind for i in range(self.n)) # 13
        m.addConstrs(gp.quicksum(z[i,t] for i in range(self.n)) >= self.min_leaf_samples * l[t] for t in self.l_ind) # 14
        m.addConstrs(gp.quicksum(z[i,t] for t in self.l_ind) == 1  for i in range(self.n)) # 15
        m.addConstrs(gp.quicksum(z[i,t] for i in range(self.n)) >= d[t] for t in self.b_ind) 
        m.addConstrs(gp.quicksum(z[i,2*t] for i in range(self.n)) >= d[t] for t in self.b_ind)
        m.addConstrs(gp.quicksum(z[i,2*t+1] for i in range(self.n)) >= d[t] for t in self.b_ind)
        for i in range(self.n):
            for t in  self.l_ind:
                m.addGenConstrIndicator(z[i,t], False, q[i] <= c[i,t])
                m.addGenConstrIndicator(z[i,t], False, q[i] >= c[i,t] - np.sqrt(self.p)*(1 - w[i,t]))
        m.addConstrs(gp.quicksum(w[i,t] * (1 - z[i,t]) for t in self.l_ind) == 1 for i in range(self.n))
        m.addConstrs(w[i,t] * z[i,t] == 0 for t in self.l_ind for i in range(self.n))
        m.addConstr(gp.quicksum(l[t] for t in self.l_ind) >= self.min_cluster_number)  
        
        # compute epsilon
        eps = self._comp_eps(x)
        
        # find ancestors and set constraints for hierarchy 
        for t in self.n_ind:
            ancestors = []
            anc = t // 2
            if t > 1:
                ancestors.append(t)
                while anc != 0:
                    ancestors.append(anc)
                    anc = anc // 2
                for k in range(len(ancestors) - 1):
                    if ancestors[k] in self.left_nodes:
                        m.addConstrs(gp.quicksum(a[j,ancestors[k+1]] * (x[i,j] + eps[j]) for j in range(self.p))
                                     +
                                     (1 + np.max(eps)) * (1 - d[ancestors[k+1]])
                                     <=
                                     b[ancestors[k+1]] + (1 + np.max(eps)) * (1 - z[i,t]) for i in range(self.n))
                    else:
                        m.addConstrs(gp.quicksum(a[j,ancestors[k+1]] * x[i,j] for j in range(self.p))
                        >=
                        b[ancestors[k+1]] - (1 - z[i,t]) for i in range(self.n))
        return m, a, b, c, d, l, z, M, q, r, s, K, w, eps

    @staticmethod
    def _comp_eps(x):
        """compute the minimum distance among all observations within each feature"""
        eps = []
        for j in range(x.shape[1]):
            xj = x[:,j]
            # drop duplicates
            xj = np.unique(xj)
            # sort
            xj = np.sort(xj)[::-1]
            # distance
            e = [1]
            for i in range(len(xj)-1):
                e.append(xj[i] - xj[i+1])
            # min distance
            eps.append(np.min(e) if np.min(e) else 1)
        return eps


    
def tree_picture(model_fit, depth, x_train):
    node_per_level = [2**i for i in range(depth+1)]
    nodes_ind = [i + 1 for i in range(2 ** (depth + 1) - 1)]
    leaf_ind = nodes_ind[-2**depth:]
    x_rectangles = list(np.arange(1,sum(node_per_level)+1,1))
    x_curr = [1+3*i for i in range(node_per_level[-1])]
    y_rectangles = []
    x_start_edges = []
    for i in range(len(node_per_level)-1,-1,-1):
        node_indx = x_rectangles[node_per_level[i]-1:node_per_level[i]*2 - 1]
        for j in range(len(node_indx)):
            x_rectangles[node_indx[j]-1] = x_curr[j]
        x_curr = [(x_curr[::2][k] + x_curr[1::2][k] + 2)/2 - 1 for k in range(int(len(x_curr)/2))]  
    for i in range(len(node_per_level)):
        y_curr = [1+2*i] * node_per_level[::-1][i] 
        y_rectangles += y_curr
    y_rectangles = y_rectangles[::-1]
    
    # draw arrows
    x_left_start = x_rectangles[:node_per_level[-1]-1]
    x_right_start = [x + 2 for x in x_left_start]
    y_left_start = y_rectangles[:node_per_level[-1]-1]
    y_right_start = y_rectangles[:node_per_level[-1]-1]
    x_start = sum([[x_left_start[i], x_right_start[i]] for i in range(len(x_left_start))], [])
    y_start = sum([[y_left_start[i], y_right_start[i]] for i in range(len(y_left_start))], [])
    x_end = [x + 1 for x in x_rectangles[1:]]
    y_end = sum([[2+2*i] * node_per_level[::-1][i] for i in range(len(node_per_level)-1)],[])[::-1]
    
    # draw recrangles
    plt.figure(figsize=(20, 10))
    ax = plt.gca()
    ax.set_ylim(0, 2*(depth+1) + 1)
    ax.set_xlim(0, 3*(2**depth) + 1)
    for i in range(len(x_rectangles)):
        ax.add_patch(Rectangle((x_rectangles[i], y_rectangles[i]), 2, 1, fill=None, alpha=1, lw=2))
        
    # create lists of a, b and d
    # d
    d = [v for k,v in model_fit.d_.items()]
    # b
    b = [round(v,5) for k,v in model_fit.b_.items()]

    # a
    a = [[0] * x_train.shape[1] for _ in range(sum(node_per_level) - node_per_level[-1])]
    for k, v in model_fit.a_.items():
        if int(v) == 1:
            a[k[1]-1][k[0]] = 1
            
    # c 
    c = [-1] * node_per_level[-1]
    for k, v in model_fit.c_.items():
        if v == 1:
            c[k[1] - node_per_level[-1]] = k[0]
            
    # K
    K = [np.round(v) for k,v in model_fit.K_.items()]
    
    # epsilon
    epsilon = [np.round(e,5) for e in model_fit.eps_]

    for i in range(len(x_rectangles[:node_per_level[-1]-1])):
        #a_ind = [list(np.where(a[i])[0] + 1)[0] if int(np.max(a[i])) > 0 else None]
        a_ind = [int(list(np.where(a[i])[0] + 1)[0]) if int(np.max(a[i])) > 0 else -1]
        if a_ind[0] > 0:
            eps_j = epsilon[a_ind[0]-1]
            b[i] = b[i] - eps_j
        else:
            a_ind[0] = None
        ax.text([x + 0.1 for x in x_rectangles[:node_per_level[-1]-1]][i], [y + 0.6 for y in y_rectangles[:node_per_level[-1]-1]][i], f'split = {bool(int(abs(d[i])))}', fontsize=36/depth)
        ax.text([x + 0.1 for x in x_rectangles[:node_per_level[-1]-1]][i], [y + 0.2 for y in y_rectangles[:node_per_level[-1]-1]][i], f'x_{a_ind[0]} <= {np.round(b[i],3)}', fontsize=36/depth)
        
    # add info for leaf nodes
    for i in range(len(x_rectangles[node_per_level[-1]-1:])):
        ax.text([x + 0.2 for x in x_rectangles[node_per_level[-1]-1:]][i], [y + 0.4 for y in y_rectangles[node_per_level[-1]-1:]][i], f'samples={int(K[i])}', fontsize=36/depth)
    # plot edges
    for i in range(len(x_start)):
        plt.plot([x_end[i],x_start[i]], [y_end[i],y_start[i]], '-', color='black')

    ax.xaxis.set_major_locator(ticker.NullLocator())
    ax.yaxis.set_major_locator(ticker.NullLocator())
    plt.show()