#!/usr/bin/env python2
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
start = time.time()


style0 = xlwt.easyxf('font: name Times New Roman, color-index red, bold on',
num_format_str='#,##0.00')
style1 = xlwt.easyxf(num_format_str='D-MMM-YY')
wb = xlwt.Workbook()


#importing metadata to a DataFrame:
#metadata_csv=pd.read_csv("ASP-POTASSCO.csv")


metadata_csv=glob.glob(os.path.join('*.csv'))   #Reading metadata csv file in the folder
                      
for f in metadata_csv:
    
    ws = wb.add_sheet(str(f))
    ws.write(0, 0, 'K', style0)
    ws.write(0, 1, "Number of Var", style0)
    ws.write(0, 2, "Regret_train", style0)
    ws.write(0, 3, 'MIP Gap', style0)
    ws.write(0, 4, 'Time (s)', style0)
    ws.write(0, 5, 'CV_prediction_mean', style0)
    ws.write(0, 6, 'CV_prediction_STD', style0)


    
    
    metadata = pd.read_csv(f)


    n_instance=len(metadata)  #number of the instances
    
    algorithms_df=pd.DataFrame(metadata.loc[:,['algo' in i for i in metadata.columns]])  #the columns contain "algo" inside
    
    algorithms_arr=algorithms_df.columns
    n_algorithms=len(algorithms_arr)  #number of the algorithms
    
    
    l=7
    for i in algorithms_arr:
        ws.write(0, l, i)
        l+=1
        
    p=algorithms_df.copy()
    
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
    CV_all_mean=np.zeros(n_algorithms)

    
    from gurobi import Model, GRB, quicksum

    for iter in range (1,n_algorithms+1):  #with iter, we can consider all cardinality constraints (K) from 1 algorithm
    #in portfolio to the max (n_algorithms)
    #for iter in range (7,8):
        K=iter
        print("K",K)
        cv_prediction=np.zeros(fold)
        CV_score=np.zeros(fold)
        #obj_true=np.zeros(fold)
        objective=np.zeros(fold)

        
        
        # #finding the optimum values for that "K" considering all metadata to find metalagorithms: these are true labels in all cases of CV (fixed)
        # #and calculating regret portfolio         
        # I_true=[i for i in range (0,n_instance)]
        # A_true=[a for a in range (0,n_algorithms)]
        # N_true=[(i,a) for a in A_true for i in I_true]
        # mdl_true=Model('alg_sel_port_true')
        # x_true=mdl_true.addVars(A_true, vtype=GRB.BINARY)
        # u_true=mdl_true.addVars(N_true,vtype=GRB.BINARY)
        # G_true= mdl_true.addVar()
        # mdl_true.modelSense=GRB.MINIMIZE
        # mdl_true.setObjective(G_true)
        # mdl_true.addConstrs((p.iloc[i,a]/max(small_number,np.min(p.iloc[i,:-1])))*u_true[i,a]<=G_true for i in I_true for a in A_true);
        # mdl_true.addConstrs((quicksum(u_true[i,a] for a in A_true) >=1 ) for i in I_true );
        # mdl_true.addConstrs( x_true[a]>=u_true[i,a]  for i in I_true for a in A_true);
        # mdl_true.addConstr(quicksum(x_true[a] for a in A_true) <= K);
        # mdl_true.Params.TimeLimit=10000
        # mdl_true.optimize()
        # selected_algorithms_true=[a for a in A_true if x_true[a].x>0.99]
        # objective_opt_true=mdl_true.objVal
        # print("OBJECTIVE_true",objective_opt_true)
        # objective=mdl_true.objVal
        # gap=mdl_true.MIPGap
        # variables=mdl_true.NumVars
        # run_time=mdl_true.Runtime
        
        
        for c in range (1,fold+1): #in cross validation, we need to run model for each fold once
        #for c in range (1,3):
        
            #Creating Test and Train Sets:
            test=metadata_cv.loc[(metadata_cv['CV'] == c)]
            train=metadata_cv.loc[(metadata_cv['CV'] != c)]
            
            n_train=len(train)
            n_test=len(test)
            
            I=[i for i in range (0,n_train)]
            A=[a for a in range (0,n_algorithms)]
            N=[(i,a) for a in A for i in I]
        
            mdl=Model('alg_sel_port')
            
            x=mdl.addVars(A, vtype=GRB.BINARY)
            u=mdl.addVars(N,vtype=GRB.BINARY)
            G= mdl.addVar()
            
            mdl.modelSense=GRB.MINIMIZE
            mdl.setObjective(G)
        
            #since we do not need to have column "CV" in our portfolio calculations, we have to ignore it: iloc[i,:-1]
            mdl.addConstrs((train.iloc[i,a]/max(small_number,np.min(train.iloc[i,:-1])))*u[i,a]<=G for i in I for a in A);
            #mdl.addConstrs((train.iloc[i,a]/np.min(train.iloc[i,:-1]))*u[i,a]<=G for i in I for a in A);

            
            mdl.addConstrs((quicksum(u[i,a] for a in A) >=1 ) for i in I );
            mdl.addConstrs( x[a]>=u[i,a]  for i in I for a in A);
            mdl.addConstr(quicksum(x[a] for a in A) <= K);
            
            #termination rules:
                #https://www.gurobi.com/documentation/9.0/refman/mip_models.html     
                #the number of discovered feasible integer solutions exceeds the specified value
                #mdl.Params.SolutionLimit=1
                
                #mdl.Params.MIPGap=0.2
            mdl.Params.TimeLimit=10000
            mdl.optimize()
            
            selected_algorithms=[a for a in A if x[a].x>0.99]
            
                #Attributes:
                #https://www.gurobi.com/documentation/9.0/refman/attributes.html
            
            
            objective[c-1]=mdl.objVal
            print("OBJECTIVE",objective[c-1])

            
        
        
            #Prediction: The results of the training model are now considered for predicting the objective function of test set:
            prediction_objective=0
            for i in range(0,n_test):
                tmp_prediction_objective=np.min(test.iloc[i][(s for s in selected_algorithms)])
                if np.min(test.iloc[i,:-1])==0:
                    tmp_prediction_objective=tmp_prediction_objective/small_number
                else:
                    tmp_prediction_objective=tmp_prediction_objective/np.min(test.iloc[i,:-1])

                prediction_objective=max(prediction_objective,tmp_prediction_objective)


            #objective=mdl.objVal
            gap=mdl.MIPGap
            variables=mdl.NumVars
            run_time=mdl.Runtime
            cv_prediction[c-1]=prediction_objective
            
            
            #Label: The results of the true model are now considered for finding the objective function of the test set:
            # label_objective=0
            # for i in range(0,n_test):
            #     tmp_label_objective=np.min(test.iloc[i][(s for s in selected_algorithms_true)])
            #     if np.min(test.iloc[i,:-1])==0:
            #         tmp_label_objective=tmp_label_objective/small_number
            #     else:
            #         tmp_label_objective=tmp_label_objective/np.min(test.iloc[i,:-1])
            #     label_objective=max(label_objective,tmp_label_objective)   
            # obj_true[c-1]=label_objective
            # print("OBJECTIVE_true",label_objective)
            
        
            # if prediction_objective ==0:
            #     CV_score[c-1]=1
            # else:
            #     CV_score[c-1]=label_objective/prediction_objective   #the ratio:  objective_of_optimise_model with metaalgorithms/objective_of_prediction_test_set
            # CV_score_mean=np.mean(CV_score)
            
        
        
        mean_cv_prediction=np.mean(cv_prediction)
        STD_cv_prediction=np.std(cv_prediction)

        CV_all_mean[K-1]=mean_cv_prediction    #record all of the mean_cv_prediction in one array to have a chart

        j=0
        for i in range(7,n_algorithms+7):
            ws.write(iter, i, x[j].x)
            j+=1
    
        
        ws.write(iter, 0, K, style0)
        ws.write(iter, 1, variables, style0)
        ws.write(iter, 2, objective[c-1], style0)
        ws.write(iter, 3, gap, style0)
        ws.write(iter, 4, run_time, style0)
        ws.write(iter, 5, mean_cv_prediction, style0)
        ws.write(iter, 6, STD_cv_prediction, style0)


    
wb.save('Results_OPT.xls')

    
process_time=time.time()-start
print("Time in min=",process_time/3660)
