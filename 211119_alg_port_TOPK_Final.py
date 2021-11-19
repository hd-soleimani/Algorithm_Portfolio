#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 10:55:36 2021

@author: Hamed Soleimani
The University of Melbourne

"""



import numpy as np
import pandas as pd
import xlwt
import os
import glob
import time
from collections import Counter
start = time.time()


style0 = xlwt.easyxf('font: name Times New Roman, color-index red, bold on',
num_format_str='#,##0.00')
style1 = xlwt.easyxf(num_format_str='D-MMM-YY')
wb = xlwt.Workbook()


epsilon=0 #becasue of Shanon calculation


def shanon_ideal(data,epsilon,n_algorithms):

    performance=np.zeros((len(data),n_algorithms))   #we want to have performance of all of the algorithms
    performance_algs=np.zeros(n_algorithms)
    
    for i in range(len(data)):
        tmp_min=np.min(data.iloc[i])

        for j in range(n_algorithms):
            if data.iloc[i,j]==tmp_min:
                performance[i,j]=1
            elif tmp_min!=0 and ((data.iloc[i,j]-tmp_min)/tmp_min)<=epsilon:
                performance[i,j]=1


    performance_algs=np.sum(performance,axis=0)  #the good performance of each algorithm

    shanon_ideal=0
    for k in range(len(performance_algs)):
        if performance_algs[k]!=0:
            shanon_ideal=shanon_ideal-(performance_algs[k]/len(data)*np.log(performance_algs[k]/len(data)))

    return performance_algs,shanon_ideal


def shanon (index,perf,data,shanon_ideal):
    shanon_port=0
    for alg in index:

        if perf[alg]!=0:
            shanon_port=shanon_port-(perf[alg]/len(data)*np.log(perf[alg]/len(data)))


    shanon_port=shanon_port/shanon_ideal
        
    return shanon_port

                      

#TopN function
def TopN(Data,n_Alg,n_Instance):
    Topi=np.zeros(n_Instance)
    for i in range(0,n_Instance):
        Topi[i]=np.argmin(Data.iloc[i,:-1])
        
    Topi_best=Counter(Topi).most_common()
    Topi_N=[t[0] for t in Topi_best]
    return Topi_N[0:n_Alg]

#importing metadata to a DataFrame:
#metadata=pd.read_csv("metadata_ASP-POTASSCO.csv")


metadata_csv=glob.glob(os.path.join('*.csv'))   #Reading metadata csv file in the folder (root folder)

for f in metadata_csv:
    
    ws = wb.add_sheet(str(f))
    ws.write(0, 0, 'K', style0)
    ws.write(0, 1, 'CV_score_mean', style0)
    ws.write(0, 2, 'CV_score_std', style0)
    ws.write(0, 3, 'Shanon_eps_'+str(epsilon), style0)




    
    
    metadata = pd.read_csv(f)


    n_instance=len(metadata)  #number of the instances
    
    algorithms_df=pd.DataFrame(metadata.loc[:,['algo' in i for i in metadata.columns]])  #the columns contain "algo" inside
    
    algorithms_arr=algorithms_df.columns
    n_algorithms=len(algorithms_arr)  #number of the algorithms
    
    
    l=4
    for i in algorithms_arr:
        ws.write(0, l, i)
        l+=1
        
    p=algorithms_df.copy()
    
    p = p.dropna(axis = 0, how = 'all')     #if all of the value of the algo are NaN, that instance has to be removed
    #p = p.dropna(axis = 0, how = 'all') #if all of the value of the algo are NaN, that instance has to be removed
    p=p.reset_index(drop=True)
    
    
    #Since there is an issue with NaN: Model select NaN as decision Variables
    #We have to replace NaNs with BigM:
    bigm=1000000000
    p = p.fillna(bigm)
    
    
    
    
    #cross validation loop:
        
    df_cv=pd.read_csv('/Users/hsoleimani/Desktop/Andres/OPT/CV/CV_metadata_'+str(f))
    metadata_cv=p.copy()
    metadata_cv['CV']=df_cv['CV']
    
    fold=10
    small_number=0.0001    #in some scenarios, there were erros due to np.min(iloc[i,:-1])=0, so in those cases, I use small_number instead

    

    #calculating Shanon_ideal:
    
    performances,shanon_best=shanon_ideal(p,epsilon,n_algorithms)
    

    for iter in range (1,n_algorithms+1):  #with iter, we can consider all cardinality constraints (K) from 1 algorithm
    #in portfolio to the max (n_algorithms)
    #for iter in range (7,8):
        K=iter
        print("K",K)
        cv_prediction=np.zeros(fold)
        CV_score=np.zeros(fold)
        obj_true=np.zeros(fold)
        shanon_TOPK=np.zeros(fold)

        
        
        for c in range (1,fold+1): #in cross validation, we need to run model for each fold once
        #for c in range (1,2):
        
            #Creating Test and Train Sets:
            test=metadata_cv.loc[(metadata_cv['CV'] == c)]
            train=metadata_cv.loc[(metadata_cv['CV'] != c)]
            
            n_train=len(train)
            
            n_test=len(test)
            
            train_no_CV=train.drop('CV',1)
            test_no_CV=test.drop('CV',1)

            
            
            selected_algorithm_train=TopN(train_no_CV,K,n_train)
            
            if len(selected_algorithm_train)<K:#since TOPK does not cover all algs:
                for i in range(0,n_algorithms): 
                    if len(selected_algorithm_train)!=K and i not in selected_algorithm_train:
                        selected_algorithm_train.append(i)
                    

            #selected_algorithm_train=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,16,17,18,19,20,21,22,23,24,25,26,27,28,29]

            #Prediction: The results of the training model are now considered for predicting the objective function of test set:
            prediction_objective=0
            
            #columns_s = [int(s) for s in selected_algorithm_train]
            
            #Regret:
            for i in range(0,n_test):
                tmp_prediction_objective=0
                tmp_prediction_objective=np.min(test.iloc[i,:-1][(int(s) for s in selected_algorithm_train)])
                #print("here:",i,"i",np.min(test_no_CV.iloc[i,:-1][(int(s) for s in selected_algorithm_train)]),(int(s) for s in selected_algorithm_train))
                #tmp_prediction_objective=np.min(test_no_CV.iloc[i , columns_s])
                
                if np.min(test.iloc[i,:-1])==0:
                    tmp_prediction_objective=tmp_prediction_objective/small_number
                    #print("zero",i,"value",np.min(test.iloc[i,:-1][(int(s) for s in selected_algorithm_train)]))
                else:
                    #print("non zero",i,"value",np.min(test.iloc[i,:-1][(int(s) for s in selected_algorithm_train)]))

                    tmp_prediction_objective=tmp_prediction_objective/np.min(test.iloc[i,:-1]) #OR test_no_CV
                #print(c,prediction_objective,tmp_prediction_objective,"and value",np.min(test.iloc[i,:-1]))
                prediction_objective=max(prediction_objective,tmp_prediction_objective)
                


            tmp_index_TOPK=[]
            for i in selected_algorithm_train:
                tmp_index_TOPK.append(int(i))

            shanon_TOPK[c-1]=shanon (tmp_index_TOPK,performances,p,shanon_best)
            cv_prediction[c-1]=prediction_objective
            
        mean_cv_prediction=np.mean(cv_prediction)   
        STD_cv_prediction=np.std(cv_prediction) 
        shanon_TOPK_mean=np.mean(shanon_TOPK)


            
        
        j=0
        for i in range(4,n_algorithms+4):
            if j in selected_algorithm_train:
                idx =[index for index in range(len(selected_algorithm_train)) if selected_algorithm_train[index] == j]
                idx=idx[0]+1  #index +1: not to start from zero
                ws.write(iter, i, idx)
            else:
                ws.write(iter, i, 0)
            j+=1
    
        
        ws.write(iter, 0, K, style0)
        ws.write(iter, 1, mean_cv_prediction, style0)
        ws.write(iter, 2, STD_cv_prediction, style0)
        ws.write(iter, 3, shanon_TOPK_mean, style0)




    
wb.save('Results_TopK.xls')

    
process_time=time.time()-start
print("Time in min=",process_time/3660)
