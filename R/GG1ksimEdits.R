# estimate the steady-state number of customers in G/G/1/cap queue
# Date : 11/10/2012            Author : Wei Xie

# Input:
#    arrival time ~ Gamma(parA),  parA matrix dim(1,2) 
#    service time ~ Gamma(parS),
#    cap - system capacity
#    nT - number of replications at design points
#    runlength - number of customers finished per run

# Ouput:
#    statN - average number customer in system
#    param_bars - internal estimates of the arrival and service time distribution parameters

GG1ksim = function(parA,parS,cap,nT,runlength)
{
  k = nrow(parA);             # number of design points
  n0 = max(nT);
  
  param_bars<-matrix(nrow=n0,ncol=2)

  statN = matrix(0,k,n0);   # average number of customers in system
  warmup = 300;
  runlength = runlength+warmup;
  
  
  for (m in 1:k){               # loop over design points
    for (j in 1:nT[m]){         # loop over replications
      temp = 0;                 # average number of customers in system
      arrN = 1;                 # record arrivals
      depN = 0;                 # record departures
      stat = 1;                 # number of customers in current system
      arrivals<-c()
      sers<-c()
      
      timer = rgamma(1,shape=parA[m,1],scale=parA[m,2]);
      arrivalT = timer + rgamma(1,shape=parA[m,1],scale=parA[m,2]);
      serviceT = timer + rgamma(1,shape=parS[m,1],scale=parS[m,2]);
      sers<-c(sers,serviceT-timer)
      arrivals<-c(arrivals,arrivalT-timer)
      
      flag = 0;
      while (depN <= runlength){
        if ((depN == warmup) && (flag == 0)){
          startT = timer;
          temp = 0;
          flag = 1;
        }
        
        if (arrivalT < serviceT){
          temp = temp + stat*(arrivalT-timer);
          timer = arrivalT;
          arrN = arrN+1;
          if (stat == 0){
            serviceT = timer+rgamma(1,shape=parS[m,1],scale=parS[m,2]);
            sers<-c(sers,serviceT-timer)
          }
          if(stat < cap){
            stat = stat+1;
          }
          arrivalT = timer+rgamma(1,shape=parA[m,1],scale=parA[m,2]);
          arrivals<-c(arrivals,arrivalT-timer)
        }else{
          temp = temp + stat*(serviceT-timer);
          stat = stat-1;
          timer = serviceT;
          depN = depN+1;
          if(stat>0){
            serviceT = timer+rgamma(1,shape=parS[m,1],scale=parS[m,2]);
            sers<-c(sers,serviceT-timer)
          }else{
            serviceT = Inf;
          }
        }
      }
      statN[m,j] = temp/(timer-startT);
      param_bars[j,1:2]<-c(1/mean(arrivals),1/mean(sers))
    }                             # end of loop over replications
  }                               # end of loop over design points
  
  OO<-list(statN,param_bars)
  return(OO);
}



MM1K_nominal <- function(theta,cap,nT,runlength){
  # function to run the M/M/1/K queue and output mean and gradient
  # theta - input parameters
  # cap - capacity of the queue
  # nT - number of replications
  # runlength - number of customers simulated
  sim_output <- GG1ksim(matrix(c(1,1/theta[1]),nrow=1),matrix(c(1,1/theta[2]),nrow=1),cap,c(nT),runlength)
  data <- data.frame(lambda=sim_output[[2]][,1],mu=sim_output[[2]][,2],Y=t(sim_output[[1]]))
  # calculate the gradient
  grad <- as.numeric(lm(Y ~ lambda + mu,data)$coefficients[2:3])
  #N3<-sim_output[[3]]
  output <- list(KPI=mean(data$Y), gradient = grad)#, N3=N3)
  return(output)
}

######### Data
data.generation <- function(params,number) {
  # params is vector of parameters for the exponential distribution
  # number is the number of observations wanted
  if (length(number)==1) {
    number <- rep(number,length(params))
  }
  observations <- list()
  for(i in 1:length(params)) {
    observations[[i]] <- rexp(number[i],rate = params[i])
  }
  return(observations)
}

######### Bootstrapping
bootstrap.bias <- function(data, B, reps, sim,...) {
  # function to estimate the bias
  # data - a list of observation sets
  # B - number of bootstraps
  simu<-match.fun(sim)
  Y_bar<-c()
  
  for(b in 1:B){
    params.star <- rep(0,length(data))
    for(i in 1:length(data)){
      m <- length(data[[i]])
      ## fit the MLE, the params are rates
      params.star[i] <-1/ mean(bootstrapping(m,data[[i]])) 
    }
    Y_bar[b] <- simu(theta=params.star,n=reps,...)$KPI
  }
  ## what would a 95% CI for mean(Y_bar) be?
  h_width<-qt(0.975,B-1)*sqrt(var(Y_bar))/sqrt(B)
  out<-list( Y_bar_bar = mean(Y_bar), h_width = h_width )
  return(out)  # Y_bar_bar being returned
}

# lam arrival rate
# mu arrival rate
## example use GG1ksim(matrix(c(1,1/lam),nrow=1),matrix(c(1,1/mu),nrow=1),1000,c(100),500)