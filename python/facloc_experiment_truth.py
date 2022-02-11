# -*- coding: utf-8 -*-
"""
Experiment based on the facloc simulation model in which we use an approximation 
of the true simulation output as the target for the bias correction.
"""
import numpy as np
import pandas as pd
import os
os.environ["KMP_DUPLICATE_LIB_OK"]="TRUE"

from facloc import facloc_nominal, facloc_mles
from experiment_functions_true import experiment_truth

from numpy.random import SeedSequence, default_rng

# experiment settings
G = 100
dim = 5
# number of bootstraps
B = 400

# ratio of reduction in predictive variance
ratio = 100

seed = 12345
seeds = SeedSequence(seed).spawn(G)

ms = [100]
ns = [100]


# true parameters
true_params = np.array([1/3, 1.26, 0.83, 1.96, 0.83])
# bound on parameters
lower_bound = np.array([0.01, -10000, 0.01, -10000, 0.01])

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

kpi_limits = [0,1]

"""Simulation function and properties."""


def facloc_data(params, number, seed = None):
    """
    A function to generate the data
    """
    num_dists = 3

    rng = default_rng(seed)
    if isinstance(number, int):
        number = [number for _ in range(num_dists)]

    observations = {'Interarrivals' : rng.exponential(scale = 1/params[0], size=number[0]),
                    'PickupTimes' : rng.lognormal(mean=params[1], sigma = params[2], size=number[1]),
                    'DeliveryTimes' : rng.lognormal(mean=params[3], sigma = params[4], size=number[2])
                    }
    return(observations)


def sim_function(params, n, est_n = False, seed=None):
    """
    Wrapper function to pass to the simulation. By default, do not produce the
    extra information needed to increase the size of n.
    """
    return(facloc_nominal(params, warehouse_coords, num_trucks, n, est_n=est_n, seed=seed))


""" The experiment. """
total_output = []

for m in ms:
    for n in ns:
        print(m,n,ratio)
        
        output, hit_boundary = experiment_truth(sim_function, facloc_data, 
                                                facloc_mles, true_params, 
                                                lower_bound, n, m, G, B, seeds, 
                                                ratio=ratio, 
                                                performance_limits=kpi_limits)
    
        output['percentage_reduction'] = (abs(output['bootstrap_bias_est']) - abs(output['final_bias']))/abs(output['bootstrap_bias_est'])
        # Overall mean
        print(output.mean())
            
        output['Observations'] = m
        output['n0'] = n
        output['ratio'] = ratio

        total_output.append(output)
            
        used = output.loc[output['iterations']>0]
        #used[['bootstrap_bias_est','final_bias']].hist()
        # mean when procedure used
        print(used.mean())


# compile the tables together
total_output = pd.concat(total_output)

file_name = 'facloc_experiment_truth.csv'
exists = os.path.exists(file_name)
total_output.to_csv(file_name,header=(not exists),index=False,mode='a')
