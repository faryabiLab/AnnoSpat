#!/usr/bin/env python
# coding: utf-8

# In[21]:


import warnings
import pandas as pd
import numpy as np
import os
import time
##import matplotlib.pyplot as plt
from sklearn import preprocessing
from sklearn.model_selection import train_test_split
import pickle
from sklearn.metrics.pairwise import euclidean_distances
from AnnoSpat.dependencies.cluster_cells import ssKMeans
from AnnoSpat.dependencies.helper_functions import *
import subprocess

warnings.filterwarnings("ignore", category=FutureWarning)


# In[26]:


def anno(data_path, signature_path, op_dir, p1, p2,toDrop, roi_name_col,disease_status_col,classifier,th,adaptive_th_vec,sampling_ratio,suffix,fileseparator):
    print('\nStarting cell type annotation.....\n')

    #toDrop=''
    
    raw_data_file = pd.read_csv(data_path, index_col=0, sep=fileseparator) # contains "cells" on rows and "proteins", "ROI_name","disease status" (optional), cell coordinates (optional) on columns
    signatures=pd.read_csv(signature_path , index_col=0) 

    # # IMC-T1D
    logtr=1 
    un=1 

    raw_data = raw_data_file.loc[:, p1:p2]# data_imct1d.transpose()
    if (len(toDrop) !=0):
        raw_data.drop(toDrop, axis =1, inplace =True)

    # get ROI name for each cell
    roi_names_allcells = raw_data_file.loc[:, roi_name_col]
    if( len(disease_status_col) == 0):
        roi_names_status_allcells = raw_data_file.loc[:, [roi_name_col]]
    else:
        roi_names_status_allcells = raw_data_file.loc[:, [roi_name_col,disease_status_col]]


    unique_rois = roi_names_status_allcells.drop_duplicates()

    unique_rois = roi_names_status_allcells.drop_duplicates()

    if(len(disease_status_col) != 0):
        noDiseaseConditions= len(unique_rois[disease_status_col ].unique() )
    else:
        noDiseaseConditions=1

    if( (len(disease_status_col) == 0) | (noDiseaseConditions ==1) ):
        tr_rois, te_rois = train_test_split(unique_rois[roi_name_col].values,test_size = 1-sampling_ratio, random_state=0)
    else:
        tr_rois, te_rois, tr_status, te_status = train_test_split(unique_rois[roi_name_col].values,unique_rois[disease_status_col ].values,test_size = 1-sampling_ratio,  stratify = unique_rois[disease_status_col ].values, random_state=0)

    raw_data_tr = raw_data.loc[ roi_names_allcells.isin(tr_rois) ]
    raw_data_te = raw_data.loc[ roi_names_allcells.isin(te_rois) ]

    dict_celltypes= { ii : signatures.columns[ii] for ii in range(len(signatures.columns)) }
    dict_celltypes.update({ len(signatures.columns): 'Unknown'})

    #np.save('../Python_objects/CTmap'+suffix+'.npy', dict_celltypes) 
    #dict_celltypes = np.load('../Python_objects/CTmap'+suffix+'.npy',allow_pickle='TRUE').item()

    imc_forSSC_raw=raw_data_tr.copy()
    imc_forSSC_raw.index=np.arange(len(imc_forSSC_raw) ) #to avoid descripency in initial labelling

    #### pre-processing data

    # #Scale input vectors individually to unit norm (vector length) so that euclidean dist in high dim space correlates with cosine distance (cells/samples lie on a unit sphere in the feature/proteomic space; direction of cells in proteomic space is more important than the magnitutde; take and example of 4 cells in 2D proteomic space)
    imc_forSSC_log =  imc_forSSC_raw
    if (logtr==1):
        imc_forSSC_log = pd.DataFrame( np.log(1+imc_forSSC_raw), index=imc_forSSC_raw.index,columns= imc_forSSC_raw.columns )


    imc_forSSC=pd.DataFrame(imc_forSSC_log, index=imc_forSSC_raw.index,columns= imc_forSSC_raw.columns )
    if(un==1):
        #log tr data (reduce the effect of outliers)
        #imc_forSSC=pd.DataFrame(  imc_forSSC, index=imc_forSSC.index,columns= imc_forSSC.columns ) #less accurate than log tr in terms of visulaization
        imc_forSSC=pd.DataFrame( preprocessing.normalize(imc_forSSC_log, norm='l2'), index=imc_forSSC_raw.index,columns= imc_forSSC_raw.columns )


    #### initialize

    k=signatures.shape[1] +1 #(1 extra cluster for unknown cell type)

    # Initialise centroids  
    centroids = np.zeros((k , imc_forSSC.shape[1]) , dtype=object)

    # Create a list to store which class/cell-type-centroid is assigned to each sample/label vector
    assigned_centroids = np.zeros(len(imc_forSSC), dtype = np.int32)-1

    #### get order of cell type abundance

    dec_abudance_ind=get_cell_type_abundance(imc_forSSC_raw , signatures)

    #### Initial labelling to get the supervised samples

    (assigned_centroids,centroids ) = initial_labelling(assigned_centroids, centroids, imc_forSSC, signatures, dec_abudance_ind, th,adaptive_th_vec)

    write_initial_labels(assigned_centroids , raw_data_tr.index,raw_data_te.index,op_dir, suffix)

    tic = time.time()
    
    print("-----Estimating cell types using Semi supervised clustering-----\n")
    
    kmeans = ssKMeans(k=k )

    kmeans.fit( imc_forSSC.values , centroids, assigned_centroids)

#     with open('models/model_kmeans.pkl', 'wb') as dictionary_file:
#         pickle.dump(kmeans, dictionary_file)
#     with open('AnnoSpat/models/model_kmeans.pkl', 'rb') as dictionary_file:
#         kmeans = pickle.load(dictionary_file)

    print("-------------------------Done-------------------------\n")
    toc=time.time()
    print("Time taken: "+str((toc-tic)/60) +"seconds\n")

    suffix

    pd.DataFrame( kmeans.cluster_centers_ , columns = imc_forSSC.columns).to_csv(op_dir+'/centroids_'+suffix+'.csv', index=False)

    assigned_labels=pd.DataFrame(kmeans.labels_)#.to_csv()
    assigned_labels.index=raw_data_tr.index

    (tr_labels,tr_labels_numericLabels) = process_cell_labels(assigned_labels,dict_celltypes,raw_data_tr.index)

    ### Train classifer

    use_trained_model =0 # 0 if predicting first time, 1 if not

    data_tr = raw_data_tr.copy() #required for training classifiers
    data_te =raw_data_te.copy()

    # classify test cells (may take about a minute)
    te_labels = classify_test_cells(classifier,data_tr, tr_labels_numericLabels['label'], data_te )

    predictNsave(classifier,suffix,op_dir, te_labels, raw_data_te.index ,tr_labels, tr_labels_numericLabels, raw_data_tr.index, dict_celltypes )

    # saved labels/cell type annotations in 
    print("Saving Annotations in: "+op_dir+'/trte_labels_numericLabels_'+classifier+suffix+'.csv')

    print("-------------------------Done-------------------------\n")
# In[27]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[19]:


### spatial analysis

import pathlib
def spatial(marks, labels, out, ip_roi):
    
    #Rscript pp.R --marks "Alpha_Cells,B_Cells" --labels /mnt/data1/gw/research/sc_integration/labels/labels_aanchal_imc.csv --out out_dir/ --input /mnt/data1/gw/research/sc_integration/data/imc/split_wide_data/20181130_HPAP003_1_C_ROI4.tif_mat.csv
    
    
    #print('pwd: ', subprocess.check_output('pwd') )
    command = 'sudo Rscript'
    SCRIPT_PATH = pathlib.Path(__file__).parent
    path2script = str(SCRIPT_PATH)+'/pp.R' #'/sample_Rscript.R'
    #print("R script at: ",path2script)
    args = [marks, labels, out, ip_roi]
    cmd = [path2script] + args#[command, path2script] + args
    cmd
    print("printing the sytem command:",cmd)
    os.system("chmod +x "+path2script)
    os.system("sudo "+path2script+" --marks "+marks+" --labels "+labels+" --out "+ out+" --input "+ip_roi)
    #subprocess.check_output(cmd, universal_newlines=True)
    #x = subprocess.check_output(cmd, universal_newlines=True)
    #print('The #rows and cols is (from R script):', x)


# In[ ]:




