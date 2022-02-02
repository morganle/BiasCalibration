# -*- coding: utf-8 -*-
"""
Functions required to run the Facility Location simulation code

@author: rhodesle
"""

from __future__ import print_function
import facloc_sim
import matlab

import pandas as pd
import numpy as np

from numpy.random import SeedSequence
from statsmodels.api import OLS

def facloc(theta, warehouse_coords, num_trucks, reps, seed):
    """
    Function to call the simulation model for the facility location problem. 
    It links with the simulation model that is built in MatLab using the facloc_sim
    package.

    Parameters
    ----------
    theta : array-like
        Input parameters for the simulation.
    warehouse_coords : array-like
        The coordinates of the warehouses.
    num_trucks : array-like
        The number of trucks at each warehouse.
    reps : int
        Number of replications of the simulation to perform.
    seed : int or SeedSequence
        The starting seed for the RNG in the simulation.

    Returns
    -------
    total_output : float
        Overall KPI
    rep_output : np.array
        KPI on each day of the simulation
    trace : np.array
        Details of the random variables generated within the simulation model
    """
    # open the MATLAB connection
    my_facloc_sim = facloc_sim.initialize()

    # convert the inputs to matlab variables
    if isinstance(theta, list):
        thetaIn = matlab.double(theta, size=(1, len(theta)))
    else:
        thetaIn = matlab.double(theta.tolist(), size=(1, 5))
    if isinstance(warehouse_coords, list):
        xIn = matlab.double(warehouse_coords, size=(1, len(warehouse_coords)))
    else:
        xIn = matlab.double(warehouse_coords.tolist(), size=(1, len(warehouse_coords)))
    if isinstance(num_trucks, list):
        numTrucksIn = matlab.double(num_trucks, size=(1, len(warehouse_coords)/2))
    else:
        numTrucksIn = matlab.double(num_trucks.tolist(), size=(1, len(warehouse_coords)/2))
    
    runlengthIn = matlab.double([reps], size=(1, 1))

    # convert the seed (by generating a state)
    if isinstance(seed,SeedSequence):
        seedIn = matlab.double([seed.generate_state(1)[0]], size=(1, 1))
    elif isinstance(seed, int):
        seedIn = matlab.double([seed], size=(1, 1))
    else:
        # terminate the MATLAB connection
        my_facloc_sim.terminate()
        raise Exception('Seed must be an integer or numpy.random.SeedSequence object.')

    # run the simulation
    total_output, rep_output, trace = my_facloc_sim.FACLOC(thetaIn, xIn, numTrucksIn, 
                                                           runlengthIn, seedIn, nargout=3)

    # terminate the MATLAB connection
    my_facloc_sim.terminate()

    # convert the output data into numpy arrays
    rep_output = np.array(rep_output._data)
    trace = np.array(trace._data).reshape(trace.size[::-1]).T

    return(total_output, rep_output, trace)

def facloc_nominal(theta, warehouse_coords, num_trucks, reps, est_n=False, seed=None):
    """
    The nominal experiment.
    """
    if not isinstance(seed, SeedSequence):
        seed = SeedSequence(seed)

    """Run the simulation."""
    kpi, rep_kpi, rvs = facloc(theta, warehouse_coords, num_trucks, reps, seed)
    
    """post simulation"""
    rvs = pd.DataFrame(rvs, columns = ['rep', 'Interarrivals', 'PickupTimes', 'DeliveryTimes'])

    data = rvs.groupby('rep').apply(facloc_mles)
    # attach the performance
    data['Y'] = rep_kpi[np.logical_not(np.isnan(rep_kpi))]

    #print(data)

    # calculate the gradient using the W&S gradient estimator
    lin_model = OLS(data['Y'], data.drop('Y', axis=1)).fit()
    grad = lin_model.params.to_numpy()

    extra = {}
    if est_n:
        # if we are wanting to estimate n, need to calculate a few more things
        # variance of simulation output
        extra['var'] = data['Y'].var()
        # variance of the residuals
        extra['s2_epsilon'] = (lin_model.resid**2).sum()/(data.shape[0] - len(theta) - 1)
        # variance covariance of internal mles
        extra['cov'] = data.drop('Y', axis=1).cov()
        extra['batch_size'] = 1

    output = {'KPI' : data['Y'].mean(), 'gradient' : grad, 'extra' : extra, 'data':data}
    return(output)


def facloc_mles(observations):
    """
    A function to calculate the MLEs of a set of observations.

    observations is a dictionary
    """
    lambda_hat = interarrival_mles(observations['Interarrivals'])
    p_mean_log_hat, p_std_log_hat = lognormal_mles(observations['PickupTimes'])
    d_mean_log_hat, d_std_log_hat = lognormal_mles(observations['DeliveryTimes'])

    return(pd.Series([lambda_hat, p_mean_log_hat, p_std_log_hat, d_mean_log_hat, d_std_log_hat]))


def lognormal_mles(data):
    log_data = np.log(data)
    mu_hat = np.mean(log_data)
    sigma_hat = np.sqrt(np.var(log_data,ddof=0))
    return(mu_hat,sigma_hat)

def interarrival_mles(data):
    return(1/np.mean(data))
