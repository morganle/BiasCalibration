# -*- coding: utf-8 -*-
"""
Created on Wed Oct  6 10:22:59 2021

@author: rhodesle
Aim is to run the simulation under the truth, mles and bias-corrected parameters
to look at the distributions of these
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import os
os.environ["KMP_DUPLICATE_LIB_OK"]="TRUE"

#from facloc import facloc
from experiment_functions import experiment

from tqdm import tqdm
from numpy.random import SeedSequence, default_rng

seed = 12345
seed = SeedSequence(seed)

m = 100
n = 100

# reps - number of replications
reps = 5000

# true parameters
true_params = np.array([1/3, 1.26, 0.83, 1.96, 0.83])

settings = [[1,3],[1,2,3],[1]]

n_trucks = 5

#%%
"""Read in the data"""
results = pd.read_csv('facloc_experiment_truth.csv')
results = results.loc[results['n0'] == n]
results = results.loc[results['Observations']==m]

mle_cols = ['mle_'+str(i) for i in range(5)]
theta_cols = ['theta_'+str(i) for i in range(5)]

for setting in settings:

    # warehouse locations
    warehouses = pd.read_csv('locations.csv')
    warehouses = warehouses.loc[warehouses['setting'].isin(setting)]
    warehouse_coords = []
    for i in warehouses.index:
        warehouse_coords.extend(warehouses.loc[i,['x','y']].tolist())
    warehouse_coords = np.array(warehouse_coords)
        
    num_warehouses = warehouses.shape[0]
    
    # truck allocation
    num_trucks = [n_trucks for _ in range(warehouses.shape[0])]

    
    all_sim_obs = pd.DataFrame()
    """run at the true parameters"""
    total_kpi, rep_kpi, trace = facloc(true_params, warehouse_coords, num_trucks, reps, seed=seed)

    all_sim_obs['true'] = rep_kpi


    """Perform replications"""
    for i in tqdm(results.index):
        # check the calibration was applied
        if results.loc[i,'iterations']>0:
            # run at mle
            mles = results.loc[i,mle_cols]
            total_kpi, rep_kpi, trace = facloc(mles, warehouse_coords, num_trucks, reps, seed=seed)
    
            all_sim_obs['mle'+str(i)] = rep_kpi
            
            # run at bias corrected
            corrected = results.loc[i,theta_cols]
            total_kpi, rep_kpi, trace = facloc(corrected, warehouse_coords, num_trucks, reps, seed=seed)
    
            all_sim_obs['corrected'+str(i)] = rep_kpi
        
            all_sim_obs[['true','mle'+str(i),'corrected'+str(i)]].hist()
            plt.show()


    all_sim_obs.to_csv('facloc_hists_'+str(num_warehouses)+'_'+str(n)+'_'+str(m)+'.csv',index=False)
    
    
#%% Now look at the resulting bias

for setting in settings:
    # warehouse locations
    warehouses = pd.read_csv('locations.csv')
    warehouses = warehouses.loc[warehouses['setting'].isin(setting)]
    num_warehouses = warehouses.shape[0]
    print("Setting with "+str(num_warehouses)+" warehouses.")
    
    sim_obs = pd.read_csv('facloc_hists_'+str(num_warehouses)+'_'+str(n)+'_'+str(m)+'.csv')
    
    mle_cols = [c for c in sim_obs.columns if 'mle' in c]
    corrected_cols = [c for c in sim_obs.columns if 'corrected' in c]
    
    means = sim_obs.mean()
    
    bias = means - means['true']
    
    mle_bias = bias[mle_cols]
    corrected_bias = bias[corrected_cols]
    
    step = (round(max([mle_bias.max(),corrected_bias.max()]),2) - round(min([mle_bias.min(),corrected_bias.min()]),2))/20
                     
    bins = np.arange(round(min([mle_bias.min(),corrected_bias.min()])-step,2),
                     round(max([mle_bias.max(),corrected_bias.max()])+step,2),
                     step = step)
    
    print("Average bias at MLE:\t%.4f" % mle_bias.mean())
    #mle_bias.hist()
    ax1 = mle_bias.plot.hist(bins=bins,alpha=0.8,title=str(num_warehouses)+' warehouses: bias',label='MLE')
    
    print("Average bias at calibrated:\t%.4f" % corrected_bias.mean())
    #corrected_bias.hist(alpha=0.5)
    ax2 = corrected_bias.plot.hist(bins=bins,alpha=0.8,label='Recalibrated')
    
    # make the plots more attractive
    ax1.set_xlabel('Bias')
    ax1.legend([ax1,ax2],['1','2'])
    
    plt.show()
    
    
    
