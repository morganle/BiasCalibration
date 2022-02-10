# source the generic functions used to run the experiment
source("experiment_functions_projected_gradient.R")
# source the simulation model, the nominal experiment function,
# the data.generation function and the bootstrapping function
source("GG1ksimEdits.R")

values_of_m <- c(100,500,5000)
values_of_lam <- c(0.6,0.9)

# number of bootstraps
B <- 400
# nT - number of replications
nT<-500
# runlength - number of customers simulated
runlength<-1000

set.seed(54321)
# number of macro-replications
G <- 100
seeds <- floor(runif(G, 1,10000))

batch.G <- G/10  

for(g in 1:batch.G) {
  batch.seeds <- seeds[1:10+(g-1)*10]
  
  # true parameters
  true.params <- c(values_of_lam[1],1)
  # number of data points
  m <- values_of_m[2]

  start_time <- Sys.time()
  results <- experiment(true.params, m, G, B, seeds, 0.001, sim=MM1K_nominal, cap=100, nT=nT , runlength=runlength )
  end_time <- Sys.time()

  print(end_time - start_time)

  write.csv(results,paste0("Results/Experiment1_",true.params[1],"_",m,"_",g,".csv"),
            row.names = F)
}
