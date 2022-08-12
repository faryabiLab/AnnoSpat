#!/usr/bin/env python
# coding: utf-8

# In[1]:

def get_data(dataname, suffix='NA'):

    gt=[]
    data_dir = '/mnt/data2/aanchal/data/'+dataname+'/'
    if (dataname=='IMC_T1D'):
        subfolder_name ='raw_data' # sub-folder containing raw data

        raw_data_file = pd.read_csv(data_dir+'/'+subfolder_name+'/mgDF.csv', index_col=0) # contains "cells" on rows and "proteins", "ROI_name" and "disease status" on columns (can also contain spatial info)

        if(suffix =='NA'):
            signatures=[]
        else:
            #signatures=pd.read_csv(data_dir+"signatures_T1D_withImmuneCelltypes_withNegMarkers"+suffix+".csv" , index_col=0)
            signatures=pd.read_csv(data_dir+"signatures_T1D"+suffix+".csv" , index_col=0) 
            

        roi_name_col = 'TIFFfilename'
        disease_status_col = 'Status'
        p1 = 'HLA.ABC' # first protein/column in tehe protein expression file
        p2 = 'Ghrelin' # last protein/column in tehe protein expression file
        toDrop=''

    return raw_data_file, signatures, roi_name_col,disease_status_col,p1,p2,toDrop, suffix, data_dir,subfolder_name,gt



#Convert ROI names to standard ROI names which can be easily programtically manipulated
def map_roi_names(old_names):
    df=old_names #raw_data_file[['cell_id','TIFFfilename']]
    df['Donor']=old_names['TIFFfilename'].str.replace(r'-', '')         .str.replace('HPAP38', 'HPAP038')         .str.replace('ICRH91', 'ICRH091')         .str.replace('ICRH99', 'ICRH099')         .str.extract(pat='([H][P][A][P]...|[I][C][R][H]...)')#.to_csv('aa.csv') # doesnt work with ROI names having '-' in donor name, but we only care abt the following 6 donors 
    #df = df[ (df['Donor']== 'HPAP002')  |  (df['Donor']== 'HPAP015' ) | (df['Donor']== 'HPAP020' )| (df['Donor']== 'HPAP005')| (df['Donor']== 'HPAP006')| (df['Donor']== 'HPAP014' )  ]
    df['Part']=old_names['TIFFfilename'].str.replace('Head', 'head')         .str.replace('HEAD', 'head')         .str.replace('HPAP003_1_', 'HPAP003_head_')         .str.replace('HPAP003_2_', 'HPAP003_body_')         .str.replace('HPAP003_3_', 'HPAP003_tail_')         .str.replace('HPAP008_1_', 'HPAP008_head_')         .str.replace('HPAP008_2_', 'HPAP008_body_')         .str.replace('HPAP008_3_', 'HPAP008_tail_')         .str.replace('HPAP019_1_', 'HPAP019_head_')         .str.replace('HPAP019_2_', 'HPAP019_body_')         .str.replace('HPAP019_3_', 'HPAP019_tail_')         .str.replace('HPAP035_1_', 'HPAP035_head_')         .str.replace('HPAP035_2_', 'HPAP035_body_')         .str.replace('HPAP035_3_', 'HPAP035_tail_')         .str.replace('HPAP036_1_', 'HPAP036_head_')         .str.replace('HPAP036_2_', 'HPAP036_body_')         .str.replace('HPAP036_3_', 'HPAP036_tail_')         .str.replace('HPAP014_1_', 'HPAP014_head_')         .str.replace('HPAP014_2_', 'HPAP014_body_')         .str.replace('HPAP014_3_', 'HPAP014_tail_')         .str.replace('HPAP006_1_', 'HPAP006_head_')         .str.replace('HPAP006_2_', 'HPAP006_body_')         .str.replace('HPAP006_3_', 'HPAP006_tail_')         .str.replace('ICRH100_4_', 'ICRH100_None_')         .str.replace('ICRH91_4_', 'ICRH91_None_')         .str.replace('ICRH99_2_', 'ICRH99_None_')         .str.extract(pat='([hbtBTN]...)')
    df['Part']=df['Part'].str.capitalize() 
    df['ROI']=old_names['TIFFfilename'].str.replace(r'-', '')         .str.replace(r' 0', '')         .str.replace(r' 1', '1')         .str.replace(r' 2', '2')         .str.replace(r' 3', '3')         .str.extract(pat='([R][O][I].)')#.to_csv('aa.csv')
        #.unique()#.str.replace(r'_', '').unique()
    df['Donor_Part_ROI']= df['Donor']+'_'+df['Part']+'_'+df['ROI']
    return df


import numpy as np
def get_cell_type_abundance(imc_forSSC_raw, signatures):
    i=0
    cta=np.zeros(len(signatures.columns))
    #cell types which are scattered all around (high expression in majority of cellls, like immune cells) might overwrite the labels of sparsely occuring cell types (like beta)
    # so traverse cell types in the decreasing order of abundance (sum of each protein over all cells)
    for celltype in (signatures.columns):

        markers=signatures[celltype].drop_duplicates().dropna().values.tolist()

        protein_exp=imc_forSSC_raw[markers[0]]#only 1 marker used for ordering abundance of cell types 
        #protein_exp=data_slice.iloc[:,0]
        celltype_specific_marker_expression=( protein_exp-np.percentile(protein_exp,0.01) )/np.percentile(protein_exp,99.99)#normalize the expression scale of marker protein 
        celltype_specific_marker_expression[celltype_specific_marker_expression>1]=1

        cell_type_abundance=  (celltype_specific_marker_expression.sum(axis=0)) #.mean(axis=1) #firs column of data_slice has cell type specific marker #proxy to #NoOfCellsInWhichProteinIsExpressed

        cta[i]=cell_type_abundance
        i=i+1
        #print(cell_type_abundance)
        
    return np.argsort(-cta)
    
import matplotlib.pyplot as plt
def initial_labelling(assigned_centroids, centroids, imc_forSSC, signatures, dec_abudance_ind, th,adaptive_th_vec):
    
    #print("-------------->",type(th),type(float(th) ))
    # imc_forSSC = scale_normalize(imc_forSSC)
    
    k = centroids.shape[0]
    # assigning unknown cells and their centroid
    inter=np.where ( imc_forSSC.index != -100000 )[0] #since index has 0, use a value which is surely not there in the index
    for i in range( len(signatures.columns) ):
        ##print (i)
        markers=signatures[ signatures.columns[i] ].drop_duplicates().dropna().values.tolist() 
        markers_neg_= [char for char in markers if  char.endswith('-')] 
        markers_neg = [x[:-1] for x in markers_neg_] #remove negative sign '-' from the name of negative markers
        markers_pos = list(set(markers) - set(markers_neg_ ))
        #print(markers_pos, markers_neg, markers_neg_)
        
        data_slice = pd.concat([imc_forSSC[markers_pos],   imc_forSSC[markers_neg].max()-imc_forSSC[markers_neg] ], axis=1)     
        
        #data_slice = scale_normalize(data_slice)###################### scale normalization ########################
        c = data_slice.prod(axis=1) #data_slice.sum(axis=1)
        
        #c = imc_forSSC[  signatures[ signatures.columns[i] ].drop_duplicates().dropna().values.tolist()  ].sum(axis=1) 
        th2 = np.percentile(c ,   float(th[1]) )#np.mean(c)- (th[1]*np.std(c)) #
        protein_negative_cells_ind =  np.where( c<= th2) [0]
        inter= set(inter).intersection ( set(protein_negative_cells_ind) )
    assigned_centroids[ list( inter ) ] = k-1
    centroids[ k -1]= imc_forSSC.loc[inter] .mean(axis=0)
    
    good_cells_ind_allcts_list=[]   
    for i in dec_abudance_ind: 

        celltype=signatures.columns[i]
        markers=signatures[celltype].drop_duplicates().dropna().values.tolist()
        #print(celltype,': ',markers)

        
        markers_neg_= [char for char in markers if char.endswith('-')] #notice the _ (this variable stores the marker names with - sign at end as given in marker file)
        markers_neg= [x[:-1] for x in markers_neg_] #marker name without - sign , to fetch proteins from imc_forSSC
        markers_pos  = list(set(markers) - set(markers_neg_ ))
        data_slice= pd.concat([imc_forSSC[markers_pos],   imc_forSSC[markers_neg].max()-imc_forSSC[markers_neg] ], axis=1) 
        #data_slice=imc_forSSC[markers]
        
        #data_slice = scale_normalize(data_slice)###################### scale normalization ########################
        celltype_scores=data_slice.prod(axis=1) #labelling with all markers

        adaptive_th=adaptive_th_vec[i]
        th0=np.percentile(celltype_scores, adaptive_th)
        #th0=np.percentile(celltype_scores, th[0])#np.mean(celltype_scores)+  (th[0]*np.std(celltype_scores) )#
        th1=np.percentile(celltype_scores, float(th[0]) )#np.mean(celltype_scores)+  (th[0]*np.std(celltype_scores) )#
        good_cells_ind= ( ( celltype_scores > th0 ) & (celltype_scores < th1 ) )
                
#         plt.figure()
#         (celltype_scores).plot.kde()
#         plt.axvline(x=  th0 )
#         plt.axvline(x=  th1 )
#         plt.title(str( signatures.columns[i]) + ', th=' + str(round(th0,2))+ ',' + str(round(th1,2)) )
#         plt.savefig('figures/cts_distribution/'+str( signatures.columns[i] )+'.png')
        
        #print( len(np.where(good_cells_ind==1)[0]) )
        indices=good_cells_ind[good_cells_ind].index

        assigned_centroids[ good_cells_ind ]=i

        good_cells=imc_forSSC.loc[  indices ]
        centroids[i] = good_cells.mean(axis=0)

        good_cells_ind_allcts_list=good_cells_ind_allcts_list+ list(  indices  )

    return (assigned_centroids,centroids)


import pandas as pd
def process_cell_labels(assigned_labels,dict_celltypes,original_index):
    labels_ssc= pd.DataFrame(assigned_labels).replace( dict_celltypes ).assign(item=original_index)
    labels_ssc.rename( columns={0:'label'}, inplace=True)
    labels_ssc = labels_ssc[['item','label']]

    labels_numericLabels= pd.DataFrame(assigned_labels)
    labels_numericLabels=labels_numericLabels.assign(item=original_index)
    labels_numericLabels.rename( columns={0:'label'}, inplace=True)
    labels_numericLabels = labels_numericLabels[['item','label']]

    return (labels_ssc,labels_numericLabels)


def predictNsave(classifier,suffix, op_dir, te_labels, te_original_index ,tr_labels_ssc,tr_labels_numericLabels, tr_original_index, dict_celltypes):
    (te_labels_ssc,te_labels_numericLabels) = process_cell_labels(te_labels,dict_celltypes,te_original_index)
    #(tr_labels_ssc,tr_labels_numericLabels) = process_cell_labels(tr_labels,dict_celltypes,tr_original_index)
    
    trte_labels = pd.concat([tr_labels_ssc, te_labels_ssc], axis= 0)
    trte_labels.sort_values(by=['item']).to_csv(op_dir+'/trte_labels_'+classifier+suffix+'.csv', index=False)
    trte_labels_numericLabels = pd.concat([tr_labels_numericLabels, te_labels_numericLabels], axis= 0)
    trte_labels_numericLabels.sort_values(by=['item']).to_csv(op_dir+'/trte_labels_numericLabels_'+classifier+suffix+'.csv', index=False)
# In[ ]:


import time
import os
import pickle
from sklearn import preprocessing
from sklearn.metrics.pairwise import euclidean_distances
def classify_test_cells(classifier,data_tr, labels_tr, data_te,suffix=[],data_dir=[] ): #, use_trained_model 
    print("Estimating cell types...")
#     if not os.path.exists('models'):
#         os.makedirs('models')
    
    if(classifier =='ELM'):
        from scipy import linalg
        from sklearn.preprocessing import OneHotEncoder
        X_train = data_tr.copy()
        y_train = labels_tr.copy()
        n_classes = len(set( y_train ))
        n_samples, n_features = X_train.shape 

        X_test = data_te.copy()

        input_size = X_train.shape[1]
        hidden_size = 1000

        np.random.seed(0)
        input_weights = np.random.normal(size=[input_size,hidden_size])
        biases = np.random.normal(size=[hidden_size])

        def relu(x):
            return np.maximum(x, 0, x)
        def hidden_nodes(X):
            G = np.dot(X, input_weights)
            G = G + biases
            H = relu(G)
            return H
        def predict(X):
            out = hidden_nodes(X)
            out = np.dot(out, output_weights)
            return out

        onehotencoder = OneHotEncoder(categories='auto')
        #scaler = MinMaxScaler()

        #X_train = scaler.fit_transform(train.values[:,1:])
        y_train = onehotencoder.fit_transform(y_train.values.reshape(-1, 1)).toarray()
        y_train
 

        st=time.time()
        print("-----Estimating cell types using ELM-----\n")
        output_weights = np.dot(linalg.pinv2(hidden_nodes(X_train)), y_train)
        end= time.time()
        print("-------------------------Done-------------------------\n")
        print("Time taken: "+str((end-st)/60) +"seconds\n")

#         with open('models/model_'+classifier+'.pkl', 'wb') as dictionary_file:
#             pickle.dump(output_weights, dictionary_file)
#         with open('/mnt/data0/apps/anaconda/anaconda2_2019_10/envs/aanchalPy3/lib/python3.7/site-packages/AnnoSpat/models/model_'+classifier+'.pkl', 'rb') as dictionary_file:
#             output_weights = pickle.load(dictionary_file)
            
        
        #PREDICT
        prediction = predict(X_test)
        #correct = 0
        te_labels = prediction.argmax(1)
        

    elif(classifier=='ssc_usingSavedCentroids'):
        X_test = data_te.copy()
        X_test_processed = preprocessing.normalize(np.log(1+X_test), norm='l2')
        saved_centroids = pd.read_csv(str(data_dir)+'labels_ssc/centroids_ssc'+str(suffix)+'.csv')

        te_labels = euclidean_distances(X_test , saved_centroids, squared=True).argmin(axis=1)


    elif(classifier=='SVM_SGD_elasticPenalty'):
        from sklearn.linear_model import SGDClassifier
        st=time.time()
        clf_svmsgd = SGDClassifier(penalty='elasticnet')
        clf_svmsgd.fit(data_tr, labels_tr )
        end=time.time()
        end-st

        # with open('models/model_'+classifier, 'wb') as dictionary_file:
        #     pickle.dump(clf_svmsgd, dictionary_file)

        # with open('models/model_'+classifier, 'rb') as dictionary_file:
        #     clf_svmsgd = pickle.load(dictionary_file)

        te_labels = clf_svmsgd.predict(data_te)
    elif(classifier == 'MLP'):
        from sklearn.neural_network import MLPClassifier
        st=time.time()
        clf_mlp = MLPClassifier(random_state=1, max_iter=300).fit(data_tr, labels_tr )
        #clf_mlp.predict_proba(X_test[:1])
        end=time.time()
        end-st

        # with open('models/model_'+classifier, 'wb') as dictionary_file:
        #     pickle.dump(clf_mlp, dictionary_file)

        # with open('models/model_'+classifier, 'rb') as dictionary_file:
        #     clf_mlp = pickle.load(dictionary_file)

        te_labels = clf_mlp.predict(data_te)
    elif (classifier =='RF'):

        from sklearn.ensemble import RandomForestClassifier
        st=time.time()
        clf_rf = RandomForestClassifier(max_depth=2, random_state=0)
        clf_rf.fit(data_tr, labels_tr )
        end=time.time()
        end-st

        classifier='RF'
        # with open('models/model_'+classifier, 'wb') as dictionary_file:
        #     pickle.dump(clf_rf, dictionary_file)
        #     #RF

        # with open('models/model_'+classifier, 'rb') as dictionary_file:
        #     clf_rf = pickle.load(dictionary_file)

        te_labels = clf_rf.predict(data_te)
        
    return te_labels






def write_initial_labels(assigned_centroids, raw_data_tr_index,raw_data_te_index,op_dir, suffix):
    initial_labels=pd.DataFrame(assigned_centroids, index=raw_data_tr_index)
    initial_labels=initial_labels.assign(item=raw_data_tr_index)
    initial_labels.rename( columns={0:'label'}, inplace=True)
    initial_labels = initial_labels[['item','label']]
    
    df_te = pd.DataFrame( np.zeros((len(raw_data_te_index)))-1   ).astype(int)
    df_te.rename( columns={0:'label'}, inplace=True)
    df_te['item'] = raw_data_te_index
    initial_labels_trte = pd.concat([initial_labels ,df_te])
    
    
    initial_labels_trte.to_csv(op_dir+'/initiallabels_numericLabels'+suffix+'.csv', index=False)
    
    
    
def scale_normalize(df):
    result = df.copy()
    for feature_name in df.columns:
        max_value = np.percentile(df[feature_name],99.99)
        min_value = np.percentile(df[feature_name],0.01)
        result[feature_name] = (df[feature_name] - min_value) / (max_value - min_value)
        result[result > 1]=1
    return result


        
        
        




