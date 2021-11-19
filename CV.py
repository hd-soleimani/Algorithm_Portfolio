#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 10 12:26:36 2021

@author: Hamed Soleimani
The University of Melbourne

"""



import numpy as np
import pandas as pd
import os
import glob
import time
start = time.time()

fold =10    #number of folds


metadata_csv=glob.glob(os.path.join('*.csv'))   #Reading metadata csv file in the folder

for f in metadata_csv:
    

    
    metadata = pd.read_csv(f)


    n_instance=len(metadata)  #number of the instances
    
    algorithms_df=pd.DataFrame(metadata.loc[:,['algo' in i for i in metadata.columns]])  #the columns contain "algo" inside
    
    algorithms_arr=algorithms_df.columns
    n_algorithms=len(algorithms_arr)  #number of the algorithms
    

        
    p=algorithms_df.copy()
    
    p = p.dropna(axis = 0, how = 'all')     #if all of the value of the algo are NaN, that instance has to be removed

    n_CV=len(p)      #number of Cross Validation we need for the metadata
    
    

    n_elements=int(np.round(n_CV/fold))  #the number of elements in each fold
    
    CV_final=pd.DataFrame()  #an empty dataframe

    CV=np.zeros(n_CV)    
    CV=pd.DataFrame(CV)
    CV=CV.rename(columns={0: 'CV'})  #Crearing dataframe "CV" with the lenghts of "len(p)" and column "CV"
    
    remain_CV=CV.copy()  #Having a copy of dataframe "CV"
    
    for i in range(1,fold):
        fold_tmp=remain_CV.sample(n=n_elements,random_state=200)   #random dample to the number=n_elements based on the remaining indecies
        fold_tmp['CV']=i   #assign fold number
        CV_final=pd.concat([CV_final, fold_tmp], axis=0)  #mergging
        remain_CV=remain_CV.drop(fold_tmp.index) #drop the selected ones from the original dataframe

    
    #last fold will be made based on the remainings:
    fold_tmp=remain_CV
    fold_tmp['CV']=fold
    CV_final=pd.concat([CV_final, fold_tmp], axis=0).sort_index()
            

    CV_final.to_csv('CV_metadata_'+str(f))  #save csv file

    
    
process_time=time.time()-start
print("Time in min=",process_time/3660)
