### SAN simulation ###

# think of a network of 13 arcs with 9 nodes (a to i). All activities must be completed
# during the project. The longest path of the network is thus the project duration.
# each arc an be thought of as an activity with a stochastic duration

# all activities follow a exponential distribution with some mean value for example

# mean activity times
#theta<-c(1,2,3,4,5,6,7,8,9,10,11,12,13)

#path_means<-c(theta[2]+theta[6]+theta[11], theta[1]+theta[3]+theta[6]+theta[11],
#                theta[1]+theta[5]+theta[11], theta[1]+theta[4]+theta[7]+theta[9]+theta[11],
#                  theta[1]+theta[4]+theta[7]+theta[10]+theta[13],
#                      theta[1]+theta[4]+theta[8]+theta[12]+theta[13])

#E_p_dur<-max(path_means) 

## lets say over 40 weeks project duration would be late and we're interested in the
## probability of the project finishing late, expected finish time and probability that 
## arc e-h is on the longest path

# n - number replications
# theta - mean activity durations
# late_thres - what time unit does the project become late?

### both functions need some editing if we want to change the performance measures they output
SAN<-function(theta,n,late_thres,duration = T){
  
  param_bars<-matrix(nrow=n,ncol=length(theta))
  stat_N<-matrix(nrow=n,ncol=1)
  
  Y<-numeric(length(n))  ## project durations
  W<-rep(0,n)  ## indicator of late or not
  Z<-rep(0,n)  ## indicator if e-h on the longest path
  for (i in 1:n){
    X<-rexp(13, rate=1/theta)
    param_bars[i,]<-X
    paths<-c(X[2]+X[6]+X[11], 
             X[1]+X[3]+X[6]+X[11],
             X[1]+X[5]+X[11], 
             X[1]+X[4]+X[7]+X[9]+X[11],
             X[1]+X[4]+X[7]+X[10]+X[13],
             X[1]+X[4]+X[8]+X[12]+X[13])
    Y[i]<- max(paths)
    W[i] <- (max(paths) > late_thres)
    if (duration) {
      stat_N[i,1]<-Y[i]
    }
    else {
      stat_N[i,1]<-W[i]
    }
  }
  OO<-list(stat_N,param_bars)
  return(OO);
}
#TEST<-SAN(theta,100,50,T)

SAN_nominal <- function(theta,n,late_thres,duration,est_n = F){
  # function to run the M/M/1/K queue and output mean and gradient
  # n - number replications
  # theta - mean activity durations
  # late_thres - what time unit does the project become late?
  sim_output <- SAN(theta,n,late_thres,duration)
  data <- data.frame( t1 = sim_output[[2]][,1], t2 = sim_output[[2]][,2], t3 = sim_output[[2]][,3],
                      t4 = sim_output[[2]][,4], t5 = sim_output[[2]][,5], t6 = sim_output[[2]][,6],
                      t7 = sim_output[[2]][,7], t8 = sim_output[[2]][,8], t9 = sim_output[[2]][,9],
                      t10 = sim_output[[2]][,10], t11 = sim_output[[2]][,11], t12 = sim_output[[2]][,12],
                      t13 = sim_output[[2]][,13], Y=sim_output[[1]])
  # first batch it
  no.batches <- max(length(theta) + 21,n/10)
  batch.size <- floor(n/no.batches)
  temp_data <- data[1:no.batches,]
  for(i in 1:no.batches) {
    batch <- 1:batch.size + (i-1)*batch.size
    temp_data[i,] <- sapply(data[batch,], mean)
  }
  
  # calculate the gradient
  lin_model = lm(Y ~ t1 + t2 + t3 + t4 + t5 + t6 + t7 + t8 + t9 + t10 + t11 + t12 + t13,temp_data)
  grad <- as.numeric(lin_model$coefficients[2:14])
  
  extra = list()
  if(est_n){
    # if we are wanting to estimate n, need to calculate a few more things
    # variance of simulation output
    extra$var <- var(data$Y)
    # variance of the residuals
    extra$s2_epsilon <- sum(lin_model$residuals^2)/(no.batches - length(theta) - 1)
    # variance covariance of internal mles
    extra$cov <- cov(temp_data[,1:length(theta)])
    extra$batch.size <- batch.size
  }
  
  output <- list(KPI=mean(data$Y), gradient = grad, extra = extra)
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
    observations[[i]] <- rexp(number[i],rate = 1/params[i])  ## edited here for the SAN example as the params are means not rates
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
      ## fit the MLE, the params are means not rates
      params.star[i] <- mean(bootstrapping(m,data[[i]])) 
    }
    Y_bar[b] <- simu(theta=params.star,n=reps,...)$KPI
  }
  ## what would a 95% CI for mean(Y_bar) be?
  h_width<-qt(0.975,B-1)*sqrt(var(Y_bar))/sqrt(B)
  out<-list( Y_bar_bar = mean(Y_bar), h_width = h_width )
  return(out)  # Y_bar_bar being returned
}

#TEST<-SAN_nominal(theta,100,50,T,T)
#TEST$extra
