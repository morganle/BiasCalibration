# -*- coding: utf-8 -*-
"""
A modified framework for the experiments that uses an estimate of the true 
performance measure as the bias corrected output to remove the dependence upon 
the bias estimator used.
"""
from experiment_functions import nominal_experiment, select_reps, bias_correction, collect_data

import pandas as pd
import numpy as np
import time

from scipy.stats import t
from numpy.random import SeedSequence

n_true = 10000

def experiment_truth(sim, data_generation, max_likelihood_estimation, 
                     true_params, bound, reps, m, G, B, seeds, ratio, 
                     performance_limits = [-np.inf,np.inf]):
    """
    A function that will perform macro-replications of the recalibration 
    procedure to correct for bias due to input modelling. It uses the Bootstrap
    bias estimator. Any simulation can be used, as long as the inputs and 
    outputs match those used in the code.

    Parameters
    ----------
    sim : function
        A function to run the simulation.
    data_generation : function
        A function to generate real world observations.
    max_likelihood_estimation : function
        A function to calculate the MLE of the observed data.
    true_params : array-like/list
        The true parameters for the system.
    bound : : array-like
        The lower bounds on input parameters.v
    reps : int
        Number of replications of the simulation.
    m : int or list
        The number of observations to generate from each distribution.
    G : int
        Number of macro-replications to perform.
    B : int
        Number of bootstrap samples to perform.
    seeds : list
        List of seeds to use for each of the macro-replications.
    ratio : int/float
        The desired comparative reduction in the variance. (In the paper, use 100.)
    performance_limits : list, optional
        A list of the upper and lower bounds for the perfomrance measure. The 
        default is [-np.inf,np.inf], i.e. no bounds.

    Returns
    -------
    output_data : pandas.DataFrame
        Contains most of the summary statistics for each macro-replication.
    boundary_hits : dictionary
        Each element is the number of times the boundary was hit in each 
        iteration of the macro-replication.
    """

    """Set up the output data structure."""
    output_data = {'macro_rep' : [i for i in range(G)],
                   'sim_reps' : [0 for i in range(G)],
                   'iterations' : [0 for i in range(G)],
                   'Bootstrap' : [0 for i in range(G)],
                   'h_width' : [0 for i in range(G)],
                   'bootstrap_bias_est' : [0 for i in range(G)],
                   'final_bias' : [0 for i in range(G)],
                   'Y_mle' : [0 for i in range(G)],
                   'Y_bar_bar' : [0 for i in range(G)],
                   'Bias_adjusted_Y' : [0 for i in range(G)],
                   'final_KPI' : [0 for i in range(G)],
                   'hit_boundary' : [0 for i in range(G)],
                   'corr_time' : [0 for i in range(G)],
                   'total_time' : [0 for i in range(G)]
                  }

    # add keys for the parameter values in each macro-replication
    for i in range(len(true_params)):
        output_data['mle_'+str(i)] = [0 for i in range(G)]
        output_data['theta_'+str(i)] = [0 for i in range(G)]

    # boundary hits dictionary
    boundary_hits = {}
    
    """Simulate at the truth first to provide a target."""
    true_results = nominal_experiment(true_params, n_true, sim, True, seeds[0])
    true_var = true_results['extra']['var']

    """Iterate over the macro-replications."""
    for g in range(G):
        begin = time.time()
        print(g)

        # Sort out how we handle the seeds
        seed = seeds[g] if isinstance(seeds[g], SeedSequence) else SeedSequence(seeds[g])

        # Generate the observed data
        observations = data_generation(true_params, m, seed.spawn(1)[0])

        # Caluclate the MLEs for the observed data
        mle = max_likelihood_estimation(observations)

        # estimate at the MLE
        mle_results = nominal_experiment(mle, reps, sim, True, seed.spawn(1)[0])
        Y_mle = mle_results['KPI']

        # originial bias
        bootstrap_bias_est = Y_mle - true_results['KPI']

        # bias adjusted value
        bias_adjusted_Y = true_results['KPI']
        print('Target:\t',bias_adjusted_Y)
        
        h_width = t.ppf(0.975,df = reps-1)*np.sqrt(mle_results['extra']['var']/reps + true_var/n_true)

        B_iter = None
        Y_bar_bar = None

        # calculate an appropriate n
        n = select_reps(mle,mle_results,reps,bias_adjusted_Y,bound,ratio=ratio)
        print('new n:\t',n)
        # if n is greater than reps, redo the nominal experiment
        if n>reps:
            mle_results = nominal_experiment(mle, n, sim, True, seed.spawn(1)[0])
	        # update the values of the mle, bias and target
            Y_mle = mle_results['KPI']
            bootstrap_bias_est = Y_mle - true_results['KPI']
            bias_adjusted_Y = Y_mle - bootstrap_bias_est
            print("New target:\t",bias_adjusted_Y)

        # only continue if target is positive
        if bias_adjusted_Y > performance_limits[0] and bias_adjusted_Y < performance_limits[1]:
            # perform the correction
            begin_correction = time.time()
            results = bias_correction(mle,bias_adjusted_Y,mle_results,bound,h_width,n,sim,seed=seed.spawn(1)[0])
            end_correction = time.time()

            print('Final output:\t',results['final_KPI'])

            # update the output data dictionary
            collect_data(output_data, results, mle, B_iter, h_width,
                         bootstrap_bias_est, Y_mle, Y_bar_bar, bias_adjusted_Y,
                         g, n)
            output_data['corr_time'][g] = end_correction-begin_correction

            # record the sequence of hitting the boundary
            boundary_hits[g] = results['pg_iterations'] if results['boundary'] else 0

        else:
            print('Bias adjusted KPI is not physically attainable')
            # update the output data dictionary
            results = {'theta_new' : [None for _ in mle]}
            collect_data(output_data, results, mle, B_iter, h_width,
                         bootstrap_bias_est, Y_mle, Y_bar_bar, bias_adjusted_Y,
                         g, n)
            boundary_hits[g] = 0

        output_data['total_time'][g] = time.time() - begin
        
        # write results to a temporary file
        pd.DataFrame(output_data).to_csv('temp_true_results.csv',index=False,header=False)
    # convert output_data into a pandas DataFrame
    output_data = pd.DataFrame(output_data)

    return(output_data, boundary_hits)