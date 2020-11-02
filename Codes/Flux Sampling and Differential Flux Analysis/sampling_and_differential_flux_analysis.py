# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 20:47:09 2020

@author: Piyush
"""

from cobra.sampling import ACHRSampler
from cobra.sampling import OptGPSampler
import cobra
import matplotlib.pyplot as plt
import scipy as sp
import numpy as np
import pandas as pd
from sklearn import decomposition
from sklearn import datasets
from sklearn.preprocessing import scale
from sklearn.linear_model import LinearRegression
import statsmodels
from scipy.stats import sem, t
from scipy import mean
from statsmodels.sandbox.stats.multicomp import multipletests
import seaborn as sns
from scipy.stats import hypergeom

modelIN=cobra.io.load_matlab_model('HumanGEMNHBESARS_cons.mat')
modelUIN=cobra.io.load_matlab_model('HumanGEMNHBEMock_cons.mat')
modelIB=cobra.io.load_matlab_model('HumanGEMBiopsySARS.mat')
modelUIB=cobra.io.load_matlab_model('HumanGEMBiopsyMock.mat')

achr_IN = ACHRSampler(modelIN, thinning=100)
achr_UIN= ACHRSampler(modelUIN,thinning=100)
achr_IB = ACHRSampler(modelIB, thinning=100)
achr_UIB= ACHRSampler(modelUIB,thinning=100)

samples_IB=achr_IB.sample(10000)
samples_UIB=achr_UIB.sample(10000)
samples_IN=achr_IN.sample(10000)
samples_UIN=achr_UIN.sample(10000)

#Mann Whitney U Test

def kstest(samplesI,samplesUI,file_name):

    rxns1=set(samples_IB.columns)
    rxns2=set(samples_UIB.columns)

    rxn_c=rxns1.intersection(rxns2)

    pvals=[]
    rxnid=[]
    fc=[]

    for rxn in rxn_c:
        data1=samples_IB[rxn].round(decimals=4)
        data2=samples_UIB[rxn].round(decimals=4)
        data1=data1.sample(n=1000)
        data2=data2.sample(n=1000)
        if((data1.std()!=0 and data1.mean()!=0) or (data2.std()!=0 and data2.mean()!=0)):
            kstat,pval=sp.stats.ks_2samp(data1,data2)
            foldc=(data1.mean()-data2.mean())/(data1.mean()+data2.mean())
            pvals.append(pval)
            rxnid.append(rxn)
            fc.append(foldc)
    data_mwu=pd.DataFrame({'Reaction':rxnid,'Pvalue':pvals})
    data_mwu=data_mwu.set_index('Reaction')
#    plt.hist(data_mwu['Pvalue'],100)
#    plt.xlabel('P-value')
#    plt.ylabel('Frequency')
#    plt.title('P-value distribution')
    reject,padj,_,_=statsmodels.stats.multitest.multipletests(data_mwu['Pvalue'], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
    data_mwu['Padj']=padj
    data_mwu['Reject']=reject
    data_mwu['FC']=fc
    data_sigFC=data_mwu.loc[(abs(data_mwu['logFC'])>0.82) & (data_mwu['Padj']<0.05),:]
    file=file_name+"_Enriched.csv"
    data_sigFC.to_csv(file)

kstest(samples_IB,samples_UIB,'Biopsy')
kstest(samples_IN,samples_UIN,'NHBE')