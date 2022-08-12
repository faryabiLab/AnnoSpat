#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
from sklearn.base import BaseEstimator
from sklearn.metrics.pairwise import euclidean_distances
class ssKMeans(BaseEstimator):

    def __init__(self, k, max_iter=10000, random_state=0, tol=1e-10):
        self.k = k
        self.max_iter = max_iter
        self.random_state = random_state
        self.tol = tol

    def _e_step(self, X):
        
        temp = euclidean_distances(X, self.cluster_centers_,
                                 squared=True).argmin(axis=1)
        self.labels_ = np.copy(self.labels_initial_)
        self.labels_[self.mask] = temp[self.mask]

               

    def _average(self, X):
        return X.mean(axis=0)

    def _m_step(self, X):
        X_center = None
        for center_id in range(self.k):
            center_mask = self.labels_ == center_id
            if not np.any(center_mask):
                # The centroid of empty clusters is set to the center of
                # everything
                if X_center is None:
                    X_center = self._average(X)
                self.cluster_centers_[center_id] = X_center
            else:
                self.cluster_centers_[center_id] =                     self._average(X[center_mask])

    def fit(self, X, centroids, assigned_centroids):
        n_samples = X.shape[0]
        vdata = np.mean(np.var(X, 0))

        #random_state = check_random_state(self.random_state)
        #self.labels_ = random_state.permutation(n_samples)[:self.k] #not actualy label but a temp variable used to initialize centers
        #self.cluster_centers_ = X[self.labels_]
        
        self.labels_initial_=assigned_centroids
        self.mask = (assigned_centroids==-1)
        self.labels_=[]
        #print(self.labels_initial_, type(self.labels_initial_))
        #print(self.mask, type(self.mask))
        
        
        self.cluster_centers_ = centroids
        

        for i in range(self.max_iter):
            #print(i)
            centers_old = self.cluster_centers_.copy()

            self._e_step(X)
            self._m_step(X)

            if np.sum((centers_old - self.cluster_centers_) ** 2) < self.tol * vdata:
                break

        return self

