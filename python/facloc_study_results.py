# -*- coding: utf-8 -*-
"""
Created on Mon Oct  4 10:00:22 2021

@author: rhodesle

Aim is to run the simulation under the truth, mles and bias-corrected parameters
to look at the distributions of these
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import os
os.environ["KMP_DUPLICATE_LIB_OK"]="TRUE"

from facloc import facloc

from tqdm import tqdm
from numpy.random import SeedSequence

seed = 12345
seed = SeedSequence(seed)

m = 100
n = 100

# true parameters
true_params = np.array([1/3, 1.26, 0.83, 1.96, 0.83])

# warehouse locations
setting = [1,3] # initially go for the 5 warehouses
n_trucks = 5

warehouses = pd.read_csv('locations.csv')
warehouses = warehouses.loc[warehouses['setting'].isin(setting)]
warehouse_coords = []
for i in warehouses.index:
    warehouse_coords.extend(warehouses.loc[i,['x','y']].tolist())
warehouse_coords = np.array(warehouse_coords)
num_trucks = [n_trucks for _ in range(warehouses.shape[0])]

"""Simulation function and properties."""
# reps - number of replications
reps = 10000


all_sim_obs = pd.DataFrame()
"""run at the true parameters"""
total_kpi, rep_kpi, trace = facloc(true_params, warehouse_coords, num_trucks, reps, seed=seed)

all_sim_obs['true'] = rep_kpi



"""Read in the data"""
results = pd.read_csv('facloc_experiment.csv')

results = results.loc[results['n0'] == n]
results = results.loc[results['Observations']==m]

mle_cols = ['mle_'+str(i) for i in range(5)]
theta_cols = ['theta_'+str(i) for i in range(5)]

for i in tqdm(results.index):
    # check the calibration was applied
    if results.loc[i,'iterations']==0:
        continue
    
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


all_sim_obs.to_csv('facloc_hists_'+str(n)+'_'+str(m)+'.csv',index=False)
