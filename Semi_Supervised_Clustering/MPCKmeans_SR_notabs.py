import numpy as np
import scipy

from active_semi_clustering.exceptions import EmptyClustersException
from active_semi_clustering.semi_supervised.pairwise_constraints.constraints import preprocess_constraints
from Step1_Impute import infer_membership_from_label
from itertools import combinations
from constrain_A import preprocess_constraints1
from math import *
from scipy.stats import pearsonr
from numpy import dot, arccos, sin

np.seterr('raise')

class MPCKMeans_sr:
    "MPCK-Means-S-D that learns only a single (S) diagonal (D) matrix"
    def __init__(self, n_clusters=3, max_iter=50, w=0.5, lam=4.1):
        """
        n_clusters: The number of clusters
        max_iter: The maximum number of iteration
        w: The weight of Constraint Violation Penalty
        lam: The weight of the regularizer
        """

        self.n_clusters = n_clusters
        self.max_iter = max_iter
        self.w = w
        self.lam = lam #4.1

    def fit(self, X, ml=[], cl=[], num_col = [], cat_col = [], sr_col = [], cat_dist_mat = None, H = False):
        """
        X: The 2d numpy data array to cluster
        ml: Must Link Constraints (A list of tuples)
        cl: Cannot Link Constraints (A list of tuples)
        num_col: The list containing the indeces of numerical columns
        cat_col: The list containing the indeces of categorical columns
        sr_col: The list containing groups of self-reported columns (Each group is a sub-list)
        cat_dist_mat: The list of matries defining the distance between categorical values (One element for one variable). If None, the Hamming distance will be used
        H: A boolen value indicating whether query argumenting should be used
        """

        # The function for creating transitive closure for the ML set
        def ML_transitive_closure(a):
            """
            a: The list of tuples to create the transitive closure
            """

            closure = set(a)
            while True:
                new_relations = set((x,w) for x,y in closure for q,w in closure if q == y)

                closure_until_now = closure | new_relations

                if closure_until_now == closure:
                    break

                closure = closure_until_now

            return closure
        
        # The function for creating transitive closure for the CL set
        def CL_transitive_closure(ml, cl):
            """
            ml: The list of must-link tuples
            cl: The list of cannot-link tuples
            """

            closure = set(cl)
            while True:
                new_relations = set((x,w) for x,y in closure for q,w in ml if y == q)
                new_relations1 = set((w,y) for x,y in closure for q,w in ml if x == q)

                closure_until_now = closure | new_relations | new_relations1

                if closure_until_now == closure:
                    break

                closure = closure_until_now

            return closure
                
        def abs_diff(d1, d2):
            # d2 is the previous center
            return sum([abs(d1[key] - d2[key]) for key in d2.keys()])
        
        # Preprocess constraints
        ml_graph, cl_graph, neighborhoods = preprocess_constraints(ml, cl, X.shape[0])

        all_constraints = ml + cl
        #list1, list2 = zip(*all_constraints)
        query_instances = None

        # Checking whehter to do query argumenting
        if H:
            # Augmenting constraints
            H = infer_membership_from_label(ml, cl, X.shape[0], self.n_clusters) # Infer fuzzy membership matrix H
            ml_a = set() # Augmented Similar Set
            cl_a = set() # Augmented Dissimilar Set
            
            # The unlabeled set
            U = set(combinations(range(X.shape[0]),2)) - ML_transitive_closure(ml) - CL_transitive_closure(ML_transitive_closure(ml), cl)
            # Add pairs to the augmented set
            for i, j in U:
                inner_product = H[i].dot(H[j])
                if inner_product > 1/self.n_clusters:
                    ml_a.add((i,j))
                if inner_product < 1/self.n_clusters:
                    cl_a.add((i,j))
            
            mla_graph, cla_graph = preprocess_constraints1(list(ml_a), list(cl_a), X.shape[0])

        else:
            H , mla_graph, cla_graph = None, None, None
            
        # Initialize metrics
        sr_col_a = [val for sr_col_sub in sr_col for val in sr_col_sub]
        if len(sr_col_a)>0:
            A = np.identity(X.shape[1] + len(sr_col))
        else:
            A = np.identity(X.shape[1])

        # Initialize cluster centers
        cluster_centers = self._initialize_cluster_centers(X, neighborhoods, A, num_col, cat_col, sr_col, cat_dist_mat)

        # Repeat until convergence
        for iteration in range(self.max_iter):
            prev_cluster_centers = cluster_centers.copy()
            prev_cluster_centers_nc = [prev_cluster_center[i] for i in num_col + sr_col_a for prev_cluster_center in prev_cluster_centers] # Non-categorical Cluster Center
            prev_cluster_centers_c = [prev_cluster_center[i] for i in cat_col for prev_cluster_center in prev_cluster_centers] # Categorical Cluster Center

            # Find farthest pair of points according to each metric
            farthest = self._find_farthest_pairs_of_points(X, A, num_col, cat_col, sr_col, cat_dist_mat)

            # Assign clusters
            labels = self._assign_clusters(X, cluster_centers, A, farthest, ml_graph, cl_graph, query_instances, H, self.w, self.lam, num_col, cat_col, sr_col, cat_dist_mat)
   
            # Estimate means
            cluster_centers = self._get_cluster_centers(X, labels, cat_col)
            cluster_centers_nc = [cluster_center[i] for i in num_col + sr_col_a for cluster_center in cluster_centers] # Non-categorical Cluster Center
            cluster_centers_c = [cluster_center[i] for i in cat_col for cluster_center in cluster_centers] # Categorical Cluster Center

            # Update metrics
            A = self._update_metrics(X, labels, cluster_centers, farthest, ml_graph, cl_graph, query_instances, H, self.w, self.lam, num_col, cat_col, sr_col, cat_dist_mat)

            # Check for convergence
            cluster_centers_shift_nc = [a - b for a, b in zip(prev_cluster_centers_nc, cluster_centers_nc)]
            converged_nc = np.allclose(np.array(cluster_centers_shift_nc), np.zeros(np.array(cluster_centers_shift_nc).shape), atol=1e-6, rtol=0)

            cluster_centers_shift_c = [abs_diff(cluster_centers_c[i], prev_cluster_centers_c[i]) for i in range(len(cat_col))]
            converged_c = np.allclose(np.array(cluster_centers_shift_c), np.zeros(np.array(cluster_centers_shift_c).shape), atol=1e-6, rtol=0)

            if converged_nc and converged_c: # If both categorical centers & non-categorical centers are unchanged
                print('\t', iteration, 'CONVERGED!')
                break

            elif iteration == self.max_iter - 1:
                print('\t', 'NO CONVERGED!')

        self.cluster_centers_, self.labels_, self.A_, self.ml_, self.cl_ = cluster_centers, labels, A, ml, cl

        return self

    def partial_fit(self, X, new_ml=[], new_cl=[], num_col = [], cat_col = [], sr_col = [], cat_dist_mat = None):
        """
        The function for adding new constraints to the current model
        """
        def abs_diff(d1, d2):
            # d2 is the previous center
            return sum([abs(d1[key] - d2[key]) for key in d2.keys()])
        
        ml = self.ml_ + new_ml
        cl = self.cl_ + new_cl

        # Preprocess constraints
        ml_graph, cl_graph, neighborhoods = preprocess_constraints(ml, cl, X.shape[0])

        # Initialize metrics
        A = self.A_

        # Initialize cluster centers
        cluster_centers = self.cluster_centers_
        sr_col_a = [val for sr_col_sub in sr_col for val in sr_col_sub]

        # Repeat until convergence
        for iteration in range(self.max_iter):
            prev_cluster_centers = cluster_centers.copy()
            prev_cluster_centers_nc = [prev_cluster_center[i] for i in num_col + sr_col_a for prev_cluster_center in prev_cluster_centers] # Non-categorical Cluster Center
            prev_cluster_centers_c = [prev_cluster_center[i] for i in cat_col for prev_cluster_center in prev_cluster_centers] # Categorical Cluster Center

            # Find farthest pair of points according to each metric
            farthest = self._find_farthest_pairs_of_points(X, A, num_col, cat_col, sr_col, cat_dist_mat)

            # Assign clusters
            labels = self._assign_clusters(X, cluster_centers, A, farthest, ml_graph, cl_graph, None, None, None, self.w, self.lam, num_col, cat_col, sr_col, cat_dist_mat)

            # Estimate means
            cluster_centers = self._get_cluster_centers(X, labels, cat_col)
            cluster_centers_nc = [cluster_center[i] for i in num_col + sr_col_a for cluster_center in cluster_centers] # Non-categorical Cluster Center
            cluster_centers_c = [cluster_center[i] for i in cat_col for cluster_center in cluster_centers] # Categorical Cluster Center

            # Update metrics
            A = self._update_metrics(X, labels, cluster_centers, farthest, ml_graph, cl_graph, None, None, None, self.w, self.lam, num_col, cat_col, sr_col, cat_dist_mat)

            # Check for convergence
            cluster_centers_shift_nc = [a - b for a, b in zip(prev_cluster_centers_nc, cluster_centers_nc)]
            converged_nc = np.allclose(np.array(cluster_centers_shift_nc), np.zeros(np.array(cluster_centers_shift_nc).shape), atol=1e-6, rtol=0)

            cluster_centers_shift_c = [abs_diff(cluster_centers_c[i], prev_cluster_centers_c[i]) for i in range(len(cat_col))]
            converged_c = np.allclose(np.array(cluster_centers_shift_c), np.zeros(np.array(cluster_centers_shift_c).shape), atol=1e-6, rtol=0)

            if converged_nc and converged_c: # If both categorical centers & non-categorical centers are unchanged
                print('\t', iteration, 'CONVERGED!')
                break

        self.cluster_centers_, self.labels_, self.A_, self.ml_, self.cl_ = cluster_centers, labels, A, ml, cl

        return self

    def _find_farthest_pairs_of_points(self, X, A, num_col, cat_col, sr_col, cat_dist_mat):
        """
        The function for finidng the farthest pairs of points based on current distance metric
        """

        farthest = None
        n = X.shape[0]
        max_distance = 0

        for i in range(n):
            for j in range(i + 1, n):
                distance = self._dist(X[i], X[j], A, num_col, cat_col, sr_col, cat_dist_mat)
                if distance > max_distance:
                    max_distance = distance
                    farthest = (i, j, distance)

        assert farthest is not None

        return farthest

    def _initialize_cluster_centers(self, X, neighborhoods, A, num_col, cat_col, sr_col, cat_dist_mat):
        """
        The function for initializing cluster centers based on the given constraints
        """
        neighborhood_centers = [X[neighborhood].mean(axis=0).tolist() for neighborhood in neighborhoods]
        for i in range(len(neighborhoods)):
            ni = len(neighborhoods[i])
            for j in cat_col:
                unique, counts = np.unique(X[np.ix_([neighborhoods[i]][0], [j])], return_counts=True)
                results = dict(zip(unique, counts/ni))
                all_keys = np.unique(X[:,j]).astype(int)
                for key in all_keys:
                    if key not in results.keys():
                            results[key] = 0
                neighborhood_centers[i][j] = results
        neighborhood_sizes = np.array([len(neighborhood) for neighborhood in neighborhoods])
        neighborhood_weights = neighborhood_sizes / neighborhood_sizes.sum()

        if len(neighborhoods) > self.n_clusters:
            cluster_centers = neighborhood_centers[self._weighted_farthest_first_traversal(np.array(neighborhood_centers), neighborhood_weights, self.n_clusters, 
                                                                                           A, num_col, cat_col, sr_col)]
        else:
            if len(neighborhoods) > 0:
                cluster_centers = neighborhood_centers
            else:
                cluster_centers = np.empty((0, X.shape[1]))

            if len(neighborhoods) < self.n_clusters:
                remaining_cluster_centers = X[np.random.choice(X.shape[0], self.n_clusters - len(neighborhoods), replace=False), :].tolist()
                for i in range(len(remaining_cluster_centers)):
                    for j in cat_col:
                            results = {int(remaining_cluster_centers[i][j]): 1}
                            all_keys = np.unique(X[:,j]).astype(int)
                            for key in all_keys:
                                    if key not in results.keys():
                                            results[key] = 0
                            remaining_cluster_centers[i][j] = results
                cluster_centers = np.concatenate([cluster_centers, remaining_cluster_centers])

        return cluster_centers
    
    def _dist(self, x, y, A, num_col, cat_col, sr_col, cat_dist_mat):
        """
        The function for calculating the distance between two data points based on the current distance matrix, i.e., (x - y)^T A (x - y)
        """

        p = x.shape[0]
        sr_col_a = [val for sr_col_sub in sr_col for val in sr_col_sub]

        if len(sr_col_a)>0:
            x2y = np.zeros(p + len(sr_col))
            x2y[sr_col_a] = x[sr_col_a] - y[sr_col_a]
            x2y[p:(p + len(sr_col))] = [sin(arccos(np.nanmin([pearsonr(x[sr_col_sub],  y[sr_col_sub])[0],0]))/2) for sr_col_sub in sr_col]
        else:
            x2y = np.zeros(p)
        
        x2y[num_col] = x[num_col] - y[num_col]
        
        if cat_dist_mat is None: # If cat_dist_mat is None, then use the Hamming distance
            x2y[cat_col] = (x[cat_col] != y[cat_col]).astype(int)
        else:   
            x2y[cat_col] = [np.array(cat_dist_mat[i])[int(x[cat_col[i]]), int(y[cat_col[i]])] for i in range(len(cat_col))]

        dist_xy = np.dot(np.dot(x2y, A), x2y)
        return sqrt(dist_xy)
    
    def _dist_x2y(self, x, y, num_col, cat_col, sr_col, cat_dist_mat):
        """
        The function for calculating the difference vector between two data points, i.e., x - y
        """
        p = x.shape[0]

        sr_col_a = [val for sr_col_sub in sr_col for val in sr_col_sub]

        if len(sr_col_a)>0:
            x2y = np.zeros(p + len(sr_col))
            x2y[sr_col_a] = x[sr_col_a] - y[sr_col_a]
            x2y[p:(p + len(sr_col))] = [sin(arccos(np.nanmin([pearsonr(x[sr_col_sub],  y[sr_col_sub])[0],0]))/2) for sr_col_sub in sr_col]
        else:
            x2y = np.zeros(p)

        x2y[num_col] = x[num_col] - y[num_col]
        
        if cat_dist_mat is None: # If cat_dist_mat is None, then use the Hamming distance
            x2y[cat_col] = (x[cat_col] != y[cat_col]).astype(int)
        else:   
            x2y[cat_col] = [np.array(cat_dist_mat[i])[int(x[cat_col[i]]), int(y[cat_col[i]])] for i in range(len(cat_col))]
        
        return x2y
    
    def _dist_xu(self, x, u, A, num_col, cat_col, sr_col, cat_dist_mat):
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
            x2u[sr_col_a] = x[sr_col_a] - u_sr
            
            # Calculate Pearson correlation and arccos
            x2u[p:(p + len(sr_col))] = [sin(arccos(np.nanmin([pearsonr(x[sr_col_sub],  np.array([u[i] for i in sr_col_sub]))[0],0]))/2) for sr_col_sub in sr_col]
        else:
            x2u = np.zeros(p)
        
        # Calculate differences for numerical and specific range columns
        x2u[num_col] = x[num_col] - u_num
        
        # Calculate differences for categorical columns
        if cat_dist_mat is None:  # Use Hamming distance if cat_dist_mat is None
            try:
                x2u[cat_col] = np.array([sum(int(x[cat_col[i]]) != int(k) * v for k, v in u_cat[i].items()) for i in range(len(cat_col))])
            except Exception as e:
                x2u[cat_col] = np.array([(len(cat_col) - 1)/len(cat_col) for _ in range(len(cat_col))])
        else:
            x2u[cat_col] = np.array([sum(np.array(cat_dist_mat[i])[int(x[cat_col[i]]), int(k)] * v for k, v in u_cat[i].items()) for i in range(len(cat_col))])
        
        # Compute the final distance
        dist_xu = np.dot(np.dot(x2u, A), x2u)
        return sqrt(dist_xu)
    
    def _dist_x2u(self, x, u, num_col, cat_col, sr_col, cat_dist_mat):
        """
        The function for calculating the difference vector between a data point and a cluster center, i.e., x - mu
        """

        p = x.shape[0]
        
        # Extract numerical, specific range, and categorical columns for u
        u_num = np.array([u[i] for i in num_col])
        u_cat = [u[i] for i in cat_col]
        sr_col_a = [val for sr_col_sub in sr_col for val in sr_col_sub]

        if len(sr_col_a) > 0:
            # Initialize x2u array
            x2u = np.zeros(p + len(sr_col))
            u_sr = np.array([u[i] for i in sr_col_a])
            x2u[sr_col_a] = x[sr_col_a] - u_sr
            # Calculate Pearson correlation and arccos
            x2u[p:(p + len(sr_col))] = [sin(arccos(np.nanmin([pearsonr(x[sr_col_sub],  np.array([u[i] for i in sr_col_sub]))[0],0]))/2) for sr_col_sub in sr_col]
        else:
            # Initialize x2u array
            x2u = np.zeros(p)
            
        # Calculate differences for numerical and specific range columns
        x2u[num_col] = x[num_col] - u_num
        # Calculate differences for categorical columns
        if cat_dist_mat is None:  # Use Hamming distance if cat_dist_mat is None
            x2u[cat_col] = np.array([
                sum(int(x[cat_col[i]]) != int(k) * v for k, v in u_cat[i].items())
                for i in range(len(cat_col))
            ])
        else:
            x2u[cat_col] = np.array([
                sum(np.array(cat_dist_mat[i])[int(x[cat_col[i]]), int(k)] * v for k, v in u_cat[i].items())
                for i in range(len(cat_col))
            ])
        
        return x2u

    def _dist_uu(self, u1, u2, A, num_col, cat_col, sr_col):
        """
        The function for calculating the distance between a data point and a cluster center, i.e., (x - mu)^T A (x - mu)
        """
        p = len(u1)

        # Extract numerical and specific range columns for u
        u1_num = np.array([u1[i] for i in num_col])
        u2_num = np.array([u2[i] for i in num_col])
        u1_cat = [u1[i] for i in cat_col]
        u2_cat = [u2[i] for i in cat_col]
        sr_col_a = [val for sr_col_sub in sr_col for val in sr_col_sub]

        if len(sr_col_a)>0:
            # Initialize x2u array
            u2u = np.zeros(p + len(sr_col))
            u1_sr = np.array([u1[i] for i in sr_col_a])
            u2_sr = np.array([u2[i] for i in sr_col_a])
            u2u[sr_col_a] = u1_sr - u2_sr
            
            # Calculate Pearson correlation and arccos
            u2u[p:(p + len(sr_col))] = [sin(arccos(np.nanmin([pearsonr(np.array([u1[i] for i in sr_col_sub]),  np.array([u2[i] for i in sr_col_sub]))[0],0]))/2) for sr_col_sub in sr_col]
        else:
            u2u = np.zeros(p)
        
        # Calculate differences for numerical and specific range columns
        u2u[num_col] = u1_num - u2_num
        
        # Calculate differences for categorical columns
        u2u[cat_col] = np.array([sum([u1_cat[i][k] - u2_cat[i][k] for k in u1_cat[i]]) for i in range(len(cat_col))])

        # Compute the final distance
        dist_uu = np.dot(np.dot(u2u, A), u2u)
        return sqrt(dist_uu)

    def _objective_fn(self, X, i, labels, cluster_centers, cluster_id, A, farthest, ml_graph, cl_graph, query_instances, H, w, lam, num_col, cat_col, sr_col, cat_dist_mat):
        """
        The function for calculating the objective function to minimize
        """

        def f_m(x_i, x_j):
            return self._dist(x_i, x_j, A, num_col, cat_col, sr_col, cat_dist_mat)

        def f_c(x_i, x_j):
            return farthest[2] - self._dist(x_i, x_j, A, num_col, cat_col, sr_col, cat_dist_mat)

        x_i = X[i]
        cluster_center = cluster_centers[cluster_id]
        
        # Calculate term_d
        term_d = self._dist_xu(x_i, cluster_center, A, num_col, cat_col, sr_col, cat_dist_mat) + lam * (np.trace(A) - np.log(np.linalg.det(A) + 1e-9))
        
        # Calculate term_m
        term_m = 0
        for j in ml_graph[i]:
            if labels[j] >= 0 and labels[j] != labels[i]:
                term_m += 2 * w * f_m(x_i, X[j]) if i > j else 2 * w * f_m(X[j], x_i)

        if H is not None:
            H_i = H[i]
            for j in range(i+1, X.shape[0]):
                if i in query_instances and j in query_instances:
                    if j not in ml_graph[i] and j not in cl_graph[i] and labels[j] != labels[i]:
                        coef = self.n_clusters / (self.n_clusters - 1) * (H_i.dot(H[j]) - 1 / self.n_clusters)
                        coef_S = max(coef, 0)
                        term_m += 1 * w * coef_S * f_m(X[j], x_i)

        # Calculate term_c
        term_c = 0
        for j in cl_graph[i]:
            if labels[j] == labels[i]:
                term_c += 2 * w * f_c(x_i, X[j]) if i > j else 2 * w * f_c(X[j], x_i)

        if H is not None:
            for j in range(i+1, X.shape[0]):
                if i in query_instances and j in query_instances:
                    if j not in ml_graph[i] and j not in cl_graph[i] and labels[j] == labels[i]:
                        coef = self.n_clusters / (self.n_clusters - 1) * (H_i.dot(H[j]) - 1 / self.n_clusters)
                        coef_D = max(-coef, 0)
                        term_c += 1 * w * coef_D * f_c(x_i, X[j]) if i > j else 1 * w * coef_D * f_c(X[j], x_i)

        return term_d + term_m + term_c


    def _assign_clusters(self, X, cluster_centers, A, farthest, ml_graph, cl_graph, query_instances, H, w, lam, num_col, cat_col, sr_col, cat_dist_mat):
        """
        The function for assigning a data point to a cluster
        """

        n_samples, n_clusters = X.shape[0], len(cluster_centers)
        labels = np.full(n_samples, fill_value=-1)

        best_labels = None
        best_objective_value = float('inf')

        for repeat in range(5):  # Repeat the assignment process five times
            current_labels = np.full(n_samples, fill_value=-1)  # Temporary labels for this iteration
            index = np.random.permutation(n_samples)  # Randomize order deterministically for this iteration

            total_objective = 0  # Track the total objective for this run
            
            for i in index:             
                # Compute objectives for each cluster
                objectives = [
                    self._objective_fn(X, i, current_labels, cluster_centers, cluster_id, A, farthest, ml_graph, cl_graph, query_instances, H, w, lam, num_col, cat_col, sr_col, cat_dist_mat)
                    for cluster_id in range(n_clusters)
                ]
                
                # Assign to the cluster with the lowest objective
                current_labels[i] = np.argmin(objectives + np.random.uniform(0, 1e-6, size=len(objectives)))  # Tie-breaking
                total_objective += min(objectives)  # Add the minimum objective to the total
            
            # Update the best assignment if this run has a lower total objective
            if total_objective < best_objective_value:
                best_labels = current_labels.copy()
                best_objective_value = total_objective

            # After all iterations, use the best assignment
            labels = best_labels

        n_samples_in_cluster = np.bincount(labels, minlength=n_clusters)
        empty_clusters = np.where(n_samples_in_cluster == 0)[0]

        while len(empty_clusters) > 0:
            # Calculate distances of all points to their assigned cluster centers
            distances = np.array([
                self._dist_xu(X[i], cluster_centers[labels[i]], A, num_col, cat_col, sr_col, cat_dist_mat)
                for i in range(n_samples)
            ])
            far_from_centers = distances.argsort()[::-1]

            for i, cluster_id in enumerate(empty_clusters):
                # Reassign the farthest points to the empty clusters
                far_point_idx = far_from_centers[i]
                new_center = X[far_point_idx]
                cluster_centers[cluster_id] = new_center
                labels[far_point_idx] = cluster_id
                n_samples_in_cluster[cluster_id] = 1

            n_samples_in_cluster = np.bincount(labels, minlength=n_clusters)
            empty_clusters = np.where(n_samples_in_cluster == 0)[0]

        return labels


    def _update_metrics(self, X, labels, cluster_centers, farthest, ml_graph, cl_graph, query_instances, H, w, lam, num_col, cat_col, sr_col, cat_dist_mat):
        # This function assumes that A is a diagonal matrix (More efficient, especially for high-dimensional data)
        """
        The function for calculating the updated metric
        """

        N, D = X.shape
        sr_col_a = [val for sr_col_sub in sr_col for val in sr_col_sub]
        if len(sr_col_a) > 0:
            p2 = len(sr_col)
        else:
            p2 = 0

        A = np.zeros((D + p2, D + p2))
        

        for d in range(D + p2):
            if d in num_col + sr_col_a:
                term_x = np.sum([(x[d] - cluster_centers[labels[i]][d]) ** 2 for i, x in enumerate(X)]) # The numeric & SR element 
            elif d >= D:
                term_x = np.sum([sin(arccos(np.nanmin([pearsonr(x[sr_col[d-D]],  [cluster_centers[labels[i]][k] for k in sr_col[d-D]])[0],0]))/2) for i, x in enumerate(X)])# The correlation element
            else:
                if cat_dist_mat is None:
                    term_x = np.sum([(x[cat_col[cat_col.index(d)]] != k).astype(int) * v for i, x in enumerate(X) for k,v in cluster_centers[labels[i]][d].items()])
                else:
                    term_x = np.sum([np.array(cat_dist_mat[cat_col.index(d)])[int(x[cat_col[cat_col.index(d)]]), int(k)] * v for i, x in enumerate(X) for k,v in cluster_centers[labels[i]][d].items()]) # The categorical element
            
            term_m = 0
            for i in range(N):
                for j in ml_graph[i]:
                    if labels[i] != labels[j]:
                        if d in num_col + sr_col_a:
                            term_m += w * (X[i,d] - X[j,d]) ** 2
                        elif d >= D:
                            term_m += w * sin(arccos(np.nanmin([pearsonr(X[i,sr_col[d-D]],  X[j,sr_col[d-D]])[0],0]))/2)
                        else:
                            if cat_dist_mat is None:
                                term_m += w * (X[i,d] != X[j,d]).astype(int)
                            else:
                                term_m += w * np.array(cat_dist_mat[cat_col.index(d)])[int(X[i,d]), int(X[j,d])]
                
                if (H is not None):
                    for j in range(i+1, X.shape[0]):
                        if i in query_instances and j in query_instances:
                            if j not in ml_graph[i] and j not in cl_graph[i] and labels[j] != labels[i]:
                                coef = self.n_clusters/(self.n_clusters - 1) * ((H[i].dot(H[j])) - 1/self.n_clusters)
                                coef_S = max(coef, 0)
                                if d in num_col + sr_col_a:
                                    term_m += w/2 * coef_S * (X[i,d] - X[j,d]) ** 2
                                elif d >= D:
                                    term_m += w/2 * coef_S * sin(arccos(np.nanmin([pearsonr(X[i,sr_col[d-D]],  X[j,sr_col[d-D]])[0], 0]))/2)
                                else:
                                    if cat_dist_mat is None:
                                        term_m += w/2 * coef_S * (X[i,d] != X[j,d]).astype(int)
                                    else:
                                        term_m += w/2 * coef_S * np.array(cat_dist_mat[cat_col.index(d)])[int(X[i,d]), int(X[j,d])]

            term_c = 0
            for i in range(N):
                for j in cl_graph[i]:
                    if labels[i] == labels[j]:
                        if d in num_col + sr_col_a:
                            tmp = ((X[farthest[0], d] - X[farthest[1], d]) ** 2 - (X[i, d] - X[j, d]) ** 2)
                        elif d >= D:
                            tmp =  sin(arccos(np.nanmin([pearsonr(X[farthest[0],sr_col[d-D]],  X[farthest[1],sr_col[d-D]])[0],0]))/2) - sin(arccos(np.nanmin([pearsonr(X[i,sr_col[d-D]],  X[j,sr_col[d-D]])[0],0]))/2)
                        else:
                            if cat_dist_mat is None:
                                tmp = (X[farthest[0],d] != X[farthest[1],d]).astype(int) -(X[i,d] != X[j,d]).astype(int)
                            else:
                                tmp = np.array(cat_dist_mat[cat_col.index(d)])[int(X[farthest[0],d]), int(X[farthest[1],d])] - np.array(cat_dist_mat[cat_col.index(d)])[int(X[i,d]), int(X[j,d])]
                        
                        term_c += w * max(tmp, 0)

                if (H is not None):
                    for j in range(i+1, X.shape[0]):
                        if i in query_instances and j in query_instances:
                            if j not in ml_graph[i] and j not in cl_graph[i] and labels[j] == labels[i]:
                                coef = self.n_clusters * ((H[i].dot(H[j])) - 1/self.n_clusters)
                                coef_D = max(-coef, 0)
                                if d in num_col + sr_col_a:
                                    tmp = ((X[farthest[0], d] - X[farthest[1], d]) ** 2 - (X[i, d] - X[j, d]) ** 2)
                                elif d >= D:
                                    tmp =  sin(arccos(np.nanmin([pearsonr(X[farthest[0],sr_col[d-D]],  X[farthest[1],sr_col[d-D]])[0],0]))/2) - sin(arccos(np.nanmin([pearsonr(X[i,sr_col[d-D]],  X[j,sr_col[d-D]])[0],0]))/2)
                                else:
                                    if cat_dist_mat is None:
                                        tmp = (X[farthest[0],d] != X[farthest[1],d]).astype(int) - (X[i,d] != X[j,d]).astype(int)
                                    else:
                                        tmp = np.array(cat_dist_mat[cat_col.index(d)])[int(X[farthest[0],d]), int(X[farthest[1],d])] - np.array(cat_dist_mat[cat_col.index(d)])[int(X[i,d]), int(X[j,d])]
                                
                                term_c += w/2 * coef_D * max(tmp, 0)

            A[d, d] = N * 1 / max(term_x/lam + term_m/lam + term_c/lam + 1, 1e-9)
        return A
    
    def _get_cluster_centers(self, X, labels, cat_col):
        """
        The function for obtaining the cluster centers based on the label assignment
        """

        clu_center = [X[labels == i].mean(axis=0).tolist() for i in range(self.n_clusters)]
        for i in range(self.n_clusters):
            ni = sum(labels == i)
            for j in cat_col:
                unique, counts = np.unique(X[np.ix_([labels == i][0], [j])], return_counts=True)
                results = dict(zip(unique, counts/ni))
                all_keys = np.unique(X[:,j]).astype(int)
                for key in all_keys:
                    if key not in results.keys():
                        results[key] = 0
                clu_center[i][j] = results
        return clu_center
    
    def _weighted_farthest_first_traversal(self, points, weights, k, A, num_col, cat_col, sr_col):
        """
        The function for finding the fartherest points to initialize cluster centers
        """
        traversed = []

        # Choose the first point randomly (weighted)
        i = np.random.choice(len(points), size=1, p=weights)[0]
        traversed.append(i)

        # Find remaining n - 1 maximally separated points
        for _ in range(k - 1):
            max_dst, max_dst_index = 0, None

            for i in range(len(points)):
                if i not in traversed:
                    dst = np.array([(self._dist_uu(points[i], points[j], A, num_col, cat_col, sr_col)) for j in traversed]).min()
                    weighted_dst = weights[i] * dst

                    if weighted_dst > max_dst:
                        max_dst = weighted_dst
                        max_dst_index = i

            traversed.append(max_dst_index)

        return traversed

    def _weighted_farthest_first_traversal1(self, X, centers, weights, A, num_col, cat_col, sr_col, cat_dist_mat):
        """
        Selects the next cluster center by finding the farthest weighted point.

        Parameters:
            X: np.ndarray
                The dataset (n_samples x n_features).
            centers: list
                List of current cluster centers.
            weights: np.ndarray
                Weights for each data point in X.
            A: np.ndarray
                Metric matrix.
            num_col: list
                Indices of numerical columns in X.
            cat_col: list
                Indices of categorical columns in X.
            sr_col: list
                Indices of special constraint features in X.
            cat_dist_mat: dict
                Precomputed distance matrix for categorical variables.

        Returns:
            int:
                Index of the farthest weighted point.
        """
        max_dst, max_dst_index = 0, None

        # Iterate through all data points
        for i in range(X.shape[0]):
            # Compute the minimum distance from the current point to any center
            weighted_dst = np.array([weights[j] * self._dist_xu(X[i,:], centers[j], A, num_col, cat_col, sr_col, cat_dist_mat)
            for j in range(len(centers))]).min()

            # Update the farthest point if necessary
            if weighted_dst > max_dst:
                max_dst = weighted_dst
                max_dst_index = i

        return max_dst_index
    
