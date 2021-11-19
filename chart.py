#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  3 15:52:02 2021

@author: hsoleimani
The University of Melbourne
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import xlrd

from collections import OrderedDict
linestyles_dict = OrderedDict(
    [('solid',               (0, ())),
     ('loosely dotted',      (0, (1, 10))),
     ('dotted',              (0, (1, 5))),
     ('densely dotted',      (0, (1, 1))),

     ('loosely dashed',      (0, (5, 10))),
     ('dashed',              (0, (5, 5))),
     ('densely dashed',      (0, (5, 1))),

     ('loosely dashdotted',  (0, (3, 10, 1, 10))),
     ('dashdotted',          (0, (3, 5, 1, 5))),
     ('densely dashdotted',  (0, (3, 1, 1, 1))),

     ('loosely dashdotdotted', (0, (3, 10, 1, 10, 1, 10))),
     ('dashdotdotted',         (0, (3, 5, 1, 5, 1, 5))),
     ('densely dashdotdotted', (0, (3, 1, 1, 1, 1, 1)))])


#input:epsilon: 0 or 0.05 or 0.1 or 0.2 or 0.5 or 1 or 2 
#epsilon = float(input('Enter epsilon value out of 0 or 0.05 or 0.1 or 0.2 or 0.5 or 1 or 2: '))

epsilon=[0,0.05,0.1,0.2,0.5,1,2]
#epsilon=[0]

for e in epsilon:
    
    #Reading the name of all sheets and then use that in a loop:
    main_excel = xlrd.open_workbook('Results_OPT_211118_Final.xls')
    sheet_array = main_excel.sheet_names()


    for i in sheet_array:

        #name of the xls files
        result_opt='Results_OPT_211118_Final.xls'
        result_topk='Results_TopK_211119_Final_eps0.xls'
        result_psfs='Results_PSFS_211111_dropallNaN_newCV.xls'
        result_rnd='Results_RND_211111_dropallNaN_newCV.xls'
        result_icarus='Results_ICARUS_211111_dropallNaN_newCV.xls'





        #read dataframe from excel files
        df_result_opt=pd.read_excel(result_opt,sheet_name=i)
        df_result_topk=pd.read_excel(result_topk,sheet_name=i)

        i=i.replace('.csv','')

        df_result_psfs=pd.read_excel(result_psfs,sheet_name=i)
        df_result_rnd=pd.read_excel(result_rnd,sheet_name=i)
        df_result_icarus=pd.read_excel(result_icarus,sheet_name=i)


        #the columns contain "algo" inside
        algorithms_df=pd.DataFrame(df_result_opt.loc[:,['algo' in i for i in df_result_opt.columns]])  

        algorithms_arr=algorithms_df.columns
        n_algorithms=len(algorithms_arr)  #number of the algorithms


        #cardinality in ICARUS:
        cardinality_e=df_result_icarus.loc[df_result_icarus['Epsilon']==e,"Cardinality"]
        cardinality_h=df_result_icarus.loc[df_result_icarus['Epsilon']==0.5,"Cardinality"]
        cardinality_k=df_result_icarus.loc[df_result_icarus['Epsilon']==2,"Cardinality"]



        #make those columns we want in arrays
        CV_mean_opt=df_result_opt['CV_prediction_mean']
        CV_mean_opt=np.log10(CV_mean_opt)

        CV_mean_topk=df_result_topk['CV_score_mean']
        CV_mean_topk=np.log10(CV_mean_topk)

        CV_mean_psfs=df_result_psfs.loc[df_result_psfs['Epsilon']==e,"CV_prediction_mean"]
        CV_mean_psfs=np.log10(CV_mean_psfs)

        CV_mean_rnd=df_result_rnd.loc[df_result_rnd['Epsilon']==e,"CV_prediction_mean"]
        CV_mean_rnd=np.log10(CV_mean_rnd)

        CV_mean_icarus=df_result_icarus.loc[df_result_icarus['Epsilon']==e,"CV_prediction_mean"]
        CV_mean_icarus=np.log10(CV_mean_icarus)

        CV_mean_psfs2=df_result_psfs.loc[df_result_psfs['Epsilon']==0.5,"CV_prediction_mean"]
        CV_mean_psfs2=np.log10(CV_mean_psfs2)

        CV_mean_rnd2=df_result_rnd.loc[df_result_rnd['Epsilon']==0.5,"CV_prediction_mean"]
        CV_mean_rnd2=np.log10(CV_mean_rnd2)

        CV_mean_icarus2=df_result_icarus.loc[df_result_icarus['Epsilon']==0.5,"CV_prediction_mean"]
        CV_mean_icarus2=np.log10(CV_mean_icarus2)

        CV_mean_psfs3=df_result_psfs.loc[df_result_psfs['Epsilon']==2,"CV_prediction_mean"]
        CV_mean_psfs3=np.log10(CV_mean_psfs3)

        CV_mean_rnd3=df_result_rnd.loc[df_result_rnd['Epsilon']==2,"CV_prediction_mean"]
        CV_mean_rnd3=np.log10(CV_mean_rnd3)

        CV_mean_icarus3=df_result_icarus.loc[df_result_icarus['Epsilon']==2,"CV_prediction_mean"]
        CV_mean_icarus3=np.log10(CV_mean_icarus3)


        X = np.linspace(1, n_algorithms, n_algorithms)
        Ya = CV_mean_opt
        Yb = CV_mean_topk

        Yc=CV_mean_psfs
        Yd=CV_mean_rnd
        Ye=CV_mean_icarus

        Yf=CV_mean_psfs2
        Yg=CV_mean_rnd2
        Yh=CV_mean_icarus2

        Yi=CV_mean_psfs3
        Yj=CV_mean_rnd3
        Yk=CV_mean_icarus3


        plt.plot(X, Ya,linestyle='-',label="MIP") # solid
        plt.plot(X, Yb,linestyle='--',label="TOPK")# dashed
        plt.plot(X, Yc,linestyle='-.',label=("PSFS, e="+str(e)))# dashdot
        plt.plot(X, Yd,linestyle=':',label=("RND, e="+str(e)))# dotted
        #plt.plot(X, Ye,'g^',label=("ICARUS, e="+str(epsilon))) # OR you can use: 'r--', 'bs' , 'g^' , 'bo'
        #plt.scatter(cardinality,Ye,'g^',label=("ICARUS, e="+str(epsilon))) # OR you can use: 'r--', 'bs' , 'g^' , 'bo'
        plt.scatter(cardinality_e,Ye,color='brown',label=("ICARUS, e="+str(e)))


        #plt.plot(X, Yf,linestyle=linestyles_dict['loosely dashdotdotted'],label="PSFS, e=0.5") 
        #plt.plot(X, Yg,linestyle=linestyles_dict['densely dotted'],label="RND, e=0.5") 
        #plt.scatter(cardinality_h,Yh,color='red',label="ICARUS, e=0.5")


        #plt.plot(X, Yi,linestyle=linestyles_dict['loosely dotted'],label="PSFS, e=2")
        #plt.plot(X, Yj,linestyle=linestyles_dict['densely dashed'],label="RND, e=2") 
        #plt.scatter(cardinality_k,Yk,color='purple',label="ICARUS, e=2")





        plt.xlabel("number of algorithms in portfolio")
        plt.ylabel("mean Regret of CV repetitions")
        #plt.legend()
        lgd=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5)) #legend outside

        plt.grid(True)


        plt.title(i)
        #plt.savefig(i+'.eps')
        plt.savefig(i+'_e'+str(e)+'.eps', bbox_extra_artists=(lgd,), bbox_inches='tight') #legend show in save fig
        #plt.savefig(i+'_eselected'+'.pdf', bbox_extra_artists=(lgd,), bbox_inches='tight') #legend show in save fig



        plt.show()
