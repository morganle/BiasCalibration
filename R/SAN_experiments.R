# source the generic functions used to run the experiment
source("experiment_functions_projected_gradient.R")
# source the simulation model, the nominal experiment function,
# the data.generation function and the bootstrapping function
source("SAN.R")

set.seed(54321)
# number of macro-replications
G <- 1000
seeds <- floor(runif(G, 1,1000000))

# true parameters
true.params <- c(1,2,3,4,5,6,7,8,9,10,11,12,13)
# bound on parameters
lower_bound <- rep(0.01,length(true.params))
# number of bootstraps
B <- 400
# nT - number of replications
nT <- 50
nT <- max(nT, length(true.params)+21)
# late threshold
late <- 60
times<-c()

for(ratio in c(4,10,25,50,100)){
  for(m in c(5000,500,50)) {
  print(paste('Ratio:',ratio,' m',m))
  start_time <- Sys.time()
  # we want the duration
  output<-experiment(true.params, lower_bound, nT, m, G, B, seeds, ratio=ratio, 0.1, sim=SAN_nominal, late_thres=late , duration = T , est_n=T)
  end_time <- Sys.time()
  print(end_time - start_time)
  
  hist(output[[1]]$iterations, main=paste(c("Iterations",m)))
  hist(output[[1]]$sim_reps, main=paste(c("Number of Simulation replications",m)))
  
  write.table(output[[1]],file=paste0("SAN_",m,"_",nT,"_",ratio,".csv"),append = FALSE, sep=",",row.names = FALSE,col.names = TRUE)
  
  # look at the number of iterations of the projected gradient approach
  boundaryHits <- output[[2]]
  iterationsNeeded <- c()
  for (i in 1:length(boundaryHits)) {
    # only bother if it was used
    iterationsNeeded <- c(iterationsNeeded,boundaryHits[[i]])
  }
  O <- c(m,ratio,"D",hist(iterationsNeeded,seq(-1,length(true.params),by=1))$counts)
  write.table(t(O),"hit_boundary.csv",sep=",",row.names = FALSE, col.names = FALSE, append = TRUE,quote = F)
  }
}


######try for the probability of being "late".#########
set.seed(54321)
# number of macro-replications
G <- 1000
seeds <- floor(runif(G, 1,1000000))

# true parameters
true.params <- c(1,2,3,4,5,6,7,8,9,10,11,12,13)
# bound on parameters
lower_bound <- rep(0.01,length(true.params))
# number of bootstraps
B <- 400
# nT - number of replications
nT <- 50

for(m in c(5000,500,50)) {
  # late threshold
  for(late in c(70,80)){
    start_time <- Sys.time()
    output<-experiment(true.params, lower_bound, nT, m, G, B, seeds, 0.1, sim=SAN_nominal, late_thres=late, duration = F, est_n=T )
    end_time <- Sys.time()
    
    print(end_time - start_time)
    hist(output[[1]]$iterations)
    
    write.table(output[[1]],file=paste0("SAN_late_",late,"_",m,"_",nT,".csv"),append = FALSE, sep=",",row.names = FALSE,col.names = TRUE) 

    # look at the number of iterations of the projected gradient approach
    boundaryHits <- output[[2]]
    iterationsNeeded <- c()
    for (i in 1:length(boundaryHits)) {
      # only bother if it was used
      iterationsNeeded <- c(iterationsNeeded,boundaryHits[[i]])
    }
    O <- c(m,late,hist(iterationsNeeded,seq(-1,length(true.params),by=1))$counts)
    write.table(t(O),"hit_boundary.csv",sep=",",row.names = FALSE, col.names = FALSE, append = TRUE,quote = F)
  }
}
