# BiasCalibration
Repository to hold code from publication "Reducing and Calibrating for Input Model Bias in Computer Simulation", INFORMS Journal on Computing.

The experiments include a Stochastic Activity Network and M/M/1/k queue (written in R) and an order delivery simulation (written in Python). These are in the two separate folders.

The generic algorithm exists in both languages:

* in R, this is contained within `experiment_functions_projected_gradient.R`

* in Python, this is contained within `experiment_functions.py`

Within the paper, the algorithm is written for both upper and lower bounds on parameter values. In this code, only lower bounds were considered.

# Example 1: M/M/1/k queue
The simulation model and associated data generation and bootstrapping functions are in `GG1ksimEdits.R`. This was not developed the authors, but by Dr Wei Xie. The script `MM1k_experiment.R` runs the experiment, currently set-up for the first experiment.

For evaluating the change in the total Mean Squared Error, use the `MSE_MM1k.R` script.

# Example 2: Stochastic Activity Network (SAN)
This example can be found in the supplementary material of the publication named above. For the Performance measures, we use both the expected project duration and the probability of the duration exceeding some threshold.

The SAN simulation model, data generation and bootstrapping functions are contained within the `SAN.R` script. The `SAN_experiments.R` and `ResultsAnalysis_SAN.R` scripts perform the experiments and analysis.

# Example 3: Order Delivery Simulation
The simulation model (`FACLOC.m`) is written in MATLAB and adapted from an facility location problem example within the SimOpt Library:

[SimOpt Library](http://www.simopt.org) (2021) Simulation optimization (SimOpt) library. ([Github repository](http://github.com/simopt-admin/simopt/wiki))

The `facloc.py` file contains details of calling the simulation model from Python. This requires installation of MATLAB Runtime. Instructions for this are in:
`BiasCalibration/python/facloc_sim/for_redistribution_files_only/GettingStarted.html`

* The necessary input file is `locations.csv`, which contains the locations of the warehouses.

* The experiment in the main paper can be run using the `facloc_experiment.py` script, with results analysis performed by `facloc_results.R`.

* The experiments in the supplementary material can be found using `facloc_experiment_truth.py`, which runs the algorithm but using an estimate of the true performance measure at the true input parameters as the target for the calibration. 

* Further simulation replications and analysis of the resulting input parameter settings is performed by `facloc_study_results.py` and `facloc_experiment_true_results.py`. This examines whether less bias is passed through the simulation when we change the system settings (in this case the number of warehouses). 

* The `facloc_ks_stat.R` also considers whether the distribution of the outputs are closer to the truth following the calibration of the input parameters.
