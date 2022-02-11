"""
Framework for the experiments.
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time

from numpy import dot
from numpy.random import SeedSequence, default_rng
from scipy.stats import t
from tqdm import trange
from copy import deepcopy

VB = True
ITER_MAX = 100

def experiment(sim, data_generation, max_likelihood_estimation, 
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

        # estimate the true performance
        Y_bar_bar, h_width = bootstrap_bias(observations, B, reps, sim, seed.spawn(1)[0],
                                            max_likelihood_estimation)

        # originial bias
        bootstrap_bias_est = Y_bar_bar-Y_mle

        # bias adjusted value
        bias_adjusted_Y = 2*Y_mle - Y_bar_bar
        print('Target:\t',bias_adjusted_Y)

        B_iter = B
        # ensure that the target is not negative
        while bias_adjusted_Y <= 0 and B_iter < 600:
            B_iter <- B_iter+100
            # estimate the true performance
            Y_bar_bar, h_width = bootstrap_bias(observations, B_iter, reps, sim,
                                                seed.spawn(1)[0],
                                                max_likelihood_estimation)

            # originial bias
            bootstrap_bias_est = Y_bar_bar-Y_mle

            # bias adjusted value
            bias_adjusted_Y = 2*Y_mle - Y_bar_bar
            print('Target:\t',bias_adjusted_Y)

        # calculate an appropriate n
        n = select_reps(mle,mle_results,reps,bias_adjusted_Y,bound,ratio=ratio)
        print('new n:\t',n)
        # if n is greater than reps, redo the nominal experiment
        if n>reps:
            mle_results = nominal_experiment(mle, n, sim, True, seed.spawn(1)[0])
	        # update the values of the mle, bias and target
            Y_mle = mle_results['KPI']
            bootstrap_bias_est = Y_bar_bar-Y_mle
            bias_adjusted_Y = 2*Y_mle - Y_bar_bar
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
        pd.DataFrame(output_data).to_csv('temp_results.csv',index=False,header=False)
    # convert output_data into a pandas DataFrame
    output_data = pd.DataFrame(output_data)

    return(output_data, boundary_hits)


"""
Update the output data.
"""
def collect_data(data, corr_results, mle, B, b_h_width,
                 bias_est, Y_mle, Y_bar_bar, Y_target, g, reps):
    """
    A function to update the summary output data.
    """
    data['sim_reps'][g] = reps
    data['Bootstrap'][g] = B
    data['h_width'][g] = b_h_width
    data['bootstrap_bias_est'][g] = bias_est
    data['Y_mle'][g] = Y_mle
    data['Y_bar_bar'][g] = Y_bar_bar
    data['Bias_adjusted_Y'][g] = Y_target

    try:
        data['iterations'][g] = corr_results['iterations']
        data['final_bias'][g] = corr_results['final_KPI'] - Y_target
        data['final_KPI'][g] = corr_results['final_KPI']
        data['hit_boundary'][g] = corr_results['boundary']
    except:
        # if these are not available, just add a None
        data['iterations'][g] = None
        data['final_bias'][g] = None
        data['final_KPI'][g] = None
        data['hit_boundary'][g] = None

    for i in range(len(mle)):
        data['mle_'+str(i)][g] = mle[i]
        data['theta_'+str(i)][g] = corr_results['theta_new'][i]


"""
Bootstrap function
"""
def bootstrap_bias(data, B, reps, sim, seed, mle_function):
    """
    A function to perform a boostrap bias estimation, first by resampling from
    the data, then calculating MLEs and running the simulation for these
    bottstrap samples.

    Parameters
    ----------
    data : dictionary
        Each element is a list of observations from the 'real- world'.
    B : int
        Number of bootstrap samples to perform.
    reps : int
        Number of replications of the simulation.
    sim : function
        A function to run the simulation.
    seed : numpy.random.SeedSequence
        The starting seed for the procedure.
    mle_function : function
        A function to calculate the MLE of the observed data.

    Returns
    -------
    Y_bar_bar : float
        The mean of the bootstrap sampled simulation outputs.
    h_width : float
        Half width of a 95% CI for Y_bar_bar
    """
    # RNG
    rng = default_rng(seed)
    Y_bar = []

    for b in trange(B):
        params_star = []
        
        bootstrap_sample = {}
        for label, obs in data.items():
            m = len(obs)
            bootstrap_sample[label] = bootstrapping(m, obs, rng)
        params_star = mle_function(bootstrap_sample)
            
        # run the simulation
        Y_bar.append(sim(params_star, reps, seed=seed.spawn(1)[0])['KPI'])
        
    # calculate the bootstrap mean
    Y_bar_bar = np.mean(Y_bar)

    ## what would a 95% CI for mean(Y_bar) be?
    h_width = t.ppf(0.975,df = B-1)*np.sqrt(np.var(Y_bar,ddof=1)/B)

    return(Y_bar_bar, h_width)  # Y_bar_bar being returned


def bootstrapping(m, data, rng):
    return(rng.choice(data,m,replace=True))


"""
Data functions
"""
def data_generation(params, number, seed = None):
    """
    A funtion to generate the data - this based on exponential only
    """
    num_dists = len(params)

    rng = default_rng(seed)
    if isinstance(number, int):
        number = [number for _ in range(num_dists)]

    observations = {i : rng.exponential(scale = params[i], size=number[i])
                    for i in range(num_dists)}
    return(observations)


"""
Nominal experiment
"""
def nominal_experiment(params, n, sim, est_n, seed):
    output = sim(params, n, est_n=est_n, seed=seed)
    return(output)

"""
Selecting the number of observations to use, n
"""
def select_reps(mle,mle_results,n0,target,limits,ratio):
    """
    A method to choose the number of replications to use within the algorithm.
    It does this by finding n such that predictive variance of the linear 
    approximation to the response surface at the proposed step is less than the 
    simulation variance divided by ratio (see equation (18) of the paper).

    Parameters
    ----------
    mle : numpy.array
        The Maximum Likelihood Estimate of the parameters (could use other estimators).
    mle_results : dictionary
        Contains the results of the nominal experiment.
    n0 : int
        Number of replications used in the nominal experiment.
    target : float
        bias corrected output.
    limits : array-like
        The lower bounds on input parameters.
    ratio : int/float
        The desired comparative reduction in the variance. (In the paper, use 100.)

    Returns
    -------
    reps : int
        Suggested number of replications to use.
    """
    # first sort results from the nominal experiment
    est = mle_results['KPI']
    grad = mle_results['gradient']
    sigma2 = mle_results['extra']['var']
    s2_epsilon = mle_results['extra']['s2_epsilon']
    var_cov = mle_results['extra']['cov']
    # invert the variance covariance matrix
    inv_var_cov = np.linalg.inv(var_cov)
    batch_size = mle_results['extra']['batch_size']

    # what would be the step in this circumstance
    take_step = produce_step(mle,est,grad,mle,target,limits)
    step = take_step['parameters'] - mle

    reps = n0
    while approx_pred_var(step,sigma2,s2_epsilon,batch_size,inv_var_cov,reps)> sigma2/ratio:
        reps += batch_size
    
    if VB:
        print('variance ratio')
        print([sigma2/reps, 
               s2_epsilon/(reps/batch_size - len(step) - 2) * step@inv_var_cov@step,
               sigma2/ratio])
    return(reps)


def approx_pred_var(step,sim_var,res_var,batch_size,inv_mle_cov,n):
    """
    Calculate the approximate predictive variance at the step.

    Parameters
    ----------
    step : numpy.array
        The step taken by the SQP algorithm.
    sim_var : float
        Estimated variance of the simulation output.
    res_var : float
        Estimated variance of the residuals of the linear model to estimate the 
        gradient.
    batch_size : int
        The size of batches if replications are batched to estimate the gradient.
    inv_mle_cov : np.array/matrix
        The inverse of the variance-covariance matrix of the MLEs.
    n : int
        Number of replications.

    Returns
    -------
    pred_var : float
        Approximate predictive variance.
    """
    
    var_r = sim_var/n
    var_grad = res_var/(n/batch_size - len(step) - 2) * (step@inv_mle_cov@step)
    pred_var = var_r + var_grad + np.sqrt(var_r * var_grad)
    return(pred_var)

"""
Bias Correction functions
"""
def bias_correction(start,target,initial_nominal,limits,epsilon,reps,sim,seed=None):
    """
    The actual function for bias correction and recalibration.

    Parameters
    ----------
    start : array-like
        Initial parameters (normally MLE).
    target : float
        bias corrected output.
    initial_nominal : dictionary
        Results from the nominal experiment.
    limits : array-like
        The bounds on variables.
    epsilon : float
        Tolerance on bias.
    reps : int
        The number of replications to be run.
    sim : function
        The simulation function.
    seed : int or SeedSequence, optional
        Seed for the random number generator. The default is None.

    Returns
    -------
    results : dictionary
        A dictionary with keys: 
            'final_kpi' - final estimated KPI
            'theta_new' - recalibrated parameters (in numpy.array)
            'iterations' - number of iterations of the algorithm
            'boundary' - was the boundary hit at any point
            'pg_iterations' - list of number of boundary hits in each iteration
    """
    print('Target:\t',target)
    print('Tolerance:\t',epsilon)
    
    seed = seed if isinstance(seed, SeedSequence) else SeedSequence(seed)
    
    mle = start
    params = start
  
    # initial nominal experiment
    est = [initial_nominal['KPI']]
    grad = initial_nominal['gradient']
    constraint = False
    boundaryHits = []
    
    # vector of bias
    bias = [est[0] - target]
    
    # counter
    iter = 0
    
    # loop whilst the bias is too big
    while abs(est[iter] - target) > epsilon and iter<ITER_MAX:   #keep going whilst the estimate is outsideof the CI for Y_bar_bar
        ## step to the new point
        take_step = produce_step(params,est[iter],grad,mle,target,limits)
        params = take_step['parameters']
        constraint = take_step['hit_boundary']
        boundaryHits.append(take_step['boundary_hits'])
        
        iter += 1
        
        # nominal experiment at current position (do not estimate n)
        output = nominal_experiment(params, reps, sim, False, seed.spawn(1)[0])
        est.append(output['KPI'])
        #if VB:
        #    print(est)
        #    print("Target:\t%.4f" %target)
        grad = output['gradient']
        # update bias
        bias.append(est[iter] - target)
        
    plot_bias_trace(bias,epsilon)
    results = {'final_KPI' : est[iter], 'theta_new' : params, 'iterations' : iter,
               'boundary' : constraint, 'pg_iterations' : boundaryHits}
    return(results)

def produce_step(params,sim_est,grad,mle,target,limits):
    """
    A method that calculates the optimal solution to the SQP sub-problem.

    Parameters
    ----------
    params : numpy.array
        The incumbent set of input parameters in the algorithm.
    sim_est : float
        The estimated KPI at this set of parameters.
    grad : numpy.array
        The estimated gradient of the response surface at this set of parameters.
    mle : numpy.array
        The Maximum Likelihood Estimate of the parameters (could use other estimators).
    target : float
        The target value for the simulation output.
    limits : list/array-like
        The lower bounds on the parameters.

    Returns
    -------
    output : dictionary
        'parameters' - numpy.array of the updated set of parameters (after step)
        'hit_boundary' - boolean indicator of whether the projected gradient was used
        'boundary_hits' - how many times the boundary was hit
    """
    # step to the new point
    step = ((target-sim_est) + dot(grad,params-mle))/dot(grad,grad) * grad - \
            (params-mle)
    
    # update parameters
    new_params = params + step
    print(new_params) if VB else None
  
    #check for feasibility
    boundaryHits = 0
    constraint = False
    proj_grad = deepcopy(grad)
    while sum(new_params<limits)>0:
        constraint = True
        boundaryHits += 1
        #if infeasible, project to feasible and set gradient component to 0
        correction = np.zeros(shape=len(params))
        for i in range(len(params)):
            if new_params[i] <= limits[i]:
                correction[i] = limits[i] - new_params[i]
                new_params[i] = limits[i]
                proj_grad[i] = 0
                
        print('Correction',correction) if VB else None
        print('proj grad', proj_grad) if VB else None

        # approximate value at new position
        proj_est = target + dot(grad,correction)
    
        # now add the projected gradient step:
        if np.linalg.norm(proj_grad) > 0:
            step = ((target-proj_est)/dot(proj_grad,proj_grad))*proj_grad
        else:
            step = np.zeros(shape=len(params))
    
        # update parameters
        new_params = new_params + step
        print("PG step:\t",new_params) if VB else None
  
    return({'parameters' : new_params, 'hit_boundary' : constraint, 
            'boundary_hits' : boundaryHits})

def plot_bias_trace(bias, hwidth):
    """
    A function to plot the trace of how the estimated bias evolves through the 
    algorithm.

    Parameters
    ----------
    bias : array-like
        An array or list of the estimated bias in each iteration of the 
        algorithm.
    hwidth : float
        The tolerance of the algorithm.
    """
    bias_series = pd.Series(bias)
    # plot
    ax = bias_series.plot(ylabel='Bias',xlabel='Iteration')
    ax.axhline(hwidth, linestyle='--',color='black')
    ax.axhline(-hwidth, linestyle='--',color='black')
    
    plt.show()
