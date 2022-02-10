# Framework for the experiments.

experiment <- function(true.params, bound, reps, m, G, B, seeds, ratio, sim, ...) {
  # true.params - the true parameters for the input distributions
  # bound       - the lower bounds on the parameter values
  # reps        - the number of replications of the simulation to run in the nominal 
  #               experiment
  # m           - the number of data points for input distributions (may be a vector)
  # G           - number of macro-replications
  # B           - number of bootstraps
  # seeds       - list of starting seeds for the experiment
  # ratio       - the ratio used for choosing the number of replications to use
  # sim         - the simulation function
  # ...         - additional argument for the simulation model

  Output.data <- data.frame(macro.rep=1:G,sim_reps=rep(reps,G),iterations=rep(0,G), Bootstrap=rep(0,G), h_width=rep(0,G),
                            bootstrap_bias_est=rep(0,G), final_bias=rep(0,G) ,
                            Y_mle = rep(0,G), Y_bar_bar=rep(0,G),
                            Bias_adjusted_Y=rep(0,G), final.KPI=rep(0,G), hit_boundary=rep(F,G) )
  # add columns for the parameters - mle and end

  # rename the data.frame
  new.names <- c(names(Output.data),rep(0,2*length(true.params)))
  for (i in 1:length(true.params)) {
    new.names[dim(Output.data)[2]+i] <- paste0("mle.",i)
    new.names[dim(Output.data)[2]+length(true.params)+i] <- paste0("theta.",i)
  }
  Output.data <- cbind(Output.data,matrix(0,nrow = G, ncol=2*length(true.params)))
  names(Output.data) <- new.names

  # boundary hits list
  boundaryHits <- list()

  # iterate over the macro-replications
  corr_times<-c()
  times<-c()
  for (g in 1:G) {
    time_in<-as.numeric(proc.time()[3])
    print(g)
    set.seed(seeds[g])

    # generate the observed data
    observations <- data.generation(true.params,m)

    # calculate the MLEs
    mle <- rep(0,length(true.params))
    for (i in 1:length(observations)) {
      mle[i] <- mean(observations[[i]]) ## edited here for the SAN example as the params are means not rates
    }
    #print(mle)

    # estimate at the MLE
    mle_results <- Nominal.experiment(mle,reps,sim,...)
    Y_mle<-mle_results[[1]]

    # estimate the true performance
    BBias <-bootstrap.bias(observations,B,reps,sim,...)
    Y_bar_bar <- BBias$Y_bar_bar
    h_width <- BBias$h_width

    # originial bias
    bootstrap_bias_est<-Y_bar_bar-Y_mle

    # bias adjusted value
    bias_adjusted_Y <- 2*Y_mle - Y_bar_bar
    print(bias_adjusted_Y)

    B_iter <- B
    # ensure that the target is not negative
    while (bias_adjusted_Y <= 0 & B_iter < 600) {
      B_iter <- B_iter+100
      # estimate the true performance
      BBias <-bootstrap.bias(observations,B_iter,sim,...)
      Y_bar_bar <- BBias$Y_bar_bar
      h_width <- BBias$h_width

      # originial bias
      bootstrap_bias_est<-Y_bar_bar-Y_mle

      # bias adjusted value
      bias_adjusted_Y <- 2*Y_mle - Y_bar_bar
      print(c(B_iter,bias_adjusted_Y))
    }

    # calculate an appropriate n
    n <- select.reps(mle,mle_results,reps,bias_adjusted_Y,bound,ratio=ratio)
    print(c('new n:',n))
    if(n>reps){
	    mle_results <- Nominal.experiment(mle,reps,sim,...)
	    # update the values of the mle, bias and target
	    Y_mle<-mle_results[[1]]
	    bootstrap_bias_est <- Y_bar_bar-Y_mle
	    bias_adjusted_Y <- 2*Y_mle - Y_bar_bar
	    print("New target")
    }

    if (bias_adjusted_Y > 0) {
      # perform the correction
      timer_in<-as.numeric(proc.time()[3])
      results <- bias.correction(start = mle, target = bias_adjusted_Y, initial_nominal = mle_results, limits = bound, epsilon = h_width, reps=n, sim, ...)
      timer<- as.numeric(proc.time()[3]) - timer_in
      print(c('Final output:',results$final.KPI))
      corr_times<-c(corr_times,timer)
      print(c("correction time",corr_times[g]))

      theta.new <- results$theta.new

      # report the results in the data.frame
      Output.data[g,] <- c(g, n, results$iterations, B_iter, h_width, bootstrap_bias_est, bias_adjusted_Y-results$final.KPI ,
                           Y_mle, Y_bar_bar, bias_adjusted_Y, results$final.KPI, results$boundary,
                           mle,theta.new)


      #record the sequence of hitting the boundary
      if (results$boundary) {
        boundaryHits[[g]] <- results$pg.iterations
      } else {
        boundaryHits[[g]] <- 0
      }
    } else {
      theta.new <- rep(NA,length(true.params))

      # report the results in the data.frame
      Output.data[g,] <- c(g, n,NA, B_iter, h_width, bootstrap_bias_est, NA,
                           Y_mle, Y_bar_bar, bias_adjusted_Y, NA, NA,
                           mle,theta.new)

      #record the sequence of hitting the boundary
      boundaryHits[[g]] <- 0
    }
    #write.table(Output.data[g,],outputFile,sep=",",row.names = FALSE, col.names = FALSE, append = TRUE)
    time_out<-as.numeric(proc.time()[3])
    times<-c(times,time_out-time_in)
  }
  print(mean(times))
  print(mean(corr_times))

  return(list(Output.data,boundaryHits,mean(times),mean(corr_times)))
}

######### Bias functions
# the function for the bootstrap bias estimator is simulation specific
# as it requires fitting the Maximum Likelihood Estimators which depend on the 
# distribution
bootstrapping<-function(m,data){
  return(sample(data,m,replace=TRUE))
}

######### Data
# the function for the data generation is simulation specific
# as it requires sampling from the true input distributions

######### Nominal experiment
Nominal.experiment <- function(params,n,sim,...) {
  # run n replications of the simulation with input parameters params
  # output has a list, gradient vector and list of performance measures
  simu<-match.fun(sim)
  output <- simu(theta=params,n=n,...)
  return(output)
}

######### method to select n
select.reps <- function(mle,mle_results,n0,target,limits,ratio){
  # a function to choose the number of replications to reduce the variance, by 
  # considering the ratio between the simulation variance and the predictive 
  # variance of the linear model at the proposed step
  #
  # mle         - the MLE for the input parameters
  # mle_results - the results from the nominal experiment
  # n0          - the initial number of replications used
  # target      - the desired performance measure from the simulation
  # limits      - lower bounds for the input parameter values
  # ratio       - the acceptable ratio between the simulation variance and the predictive variance
  
  # first sort results from the nominal experiment
  est <- mle_results$KPI
  grad <- mle_results$gradient
  sigma2 <- mle_results[[3]]$var
  s2_epsilon <- mle_results[[3]]$s2_epsilon
  var_cov <- mle_results[[3]]$cov
  inv.var.cov <- solve(var_cov)
  batch.size <- mle_results[[3]]$batch.size

  # what would be the step in this circumstance
  take_step <- produce_step(mle,est,grad,mle,target,limits)
  step <- take_step$parameters - mle

  # increase the reps until condition on predictive variance is met
  reps <- n0
  while(approx.pred.var(step,sigma2,s2_epsilon,batch.size,inv.var.cov,reps)>sigma2/ratio){
    reps <- reps+batch.size
  }
  print('variance ratio')
  print(c(sigma2/reps, s2_epsilon/(reps/batch.size - length(step) - 2) * as.numeric(t(step)%*%inv.var.cov%*%step),sigma2/ratio))
  return(reps)
}
approx.pred.var <- function(step,sim.var,res.var,batch.size,inv.mle.cov,n){
  # calculate the predictive variance of the linear model at the step
  #
  # step        - the proposed step in the algorithm
  # sim.var     - the variance of the simulation output
  # res.var     - the variance of the residuals from the linear model
  # batch.size  - the number of replications in a batch
  # inv.mle.cov - the inverse of the internal MLE covariance matrix from the 
  #               Weiland and Schmeiser gradient estimator
  # n           - the number of replications performed
  
  var.r <- sim.var/n
  var.grad <- res.var/(n/batch.size - length(step) - 2) * as.numeric(t(step)%*%inv.mle.cov%*%step)
  pred.var <- var.r + var.grad + sqrt(var.r * var.grad)
  return(pred.var)
}

######### Correction
bias.correction <- function(start,target,initial_nominal,limits,epsilon,reps,sim,...) {
  # function to perform the bias correction algorithm
  #
  # start           - the current values of the parameters
  # target          - the desired performance measure from the simulation
  # initial_nominal - results from the nominal experiment at start
  # limits          - lower bounds for the input parameter values
  # epsilon         - the tolerated distance from target used for the stopping rule
  # reps            - the number of replications of the simulation to perform
  # sim             - the simulation function
  # ...             - additional arguments for the simulation model

  # vector of parameters
  mle <- start
  params <- start

  # initial nominal experiment
  est <- c(initial_nominal$KPI)
  grad <- initial_nominal$gradient
  constraint <- F
  boundaryHits <- c()

  # vector of bias
  bias <- c(target - est[1])

  # counter
  iter <- 1

  # loop whilst the bias is too big
  while (est[iter] > (target+epsilon) | est[iter] < (target-epsilon) ) {   #keep going whilst the estimate is outside of the CI for Y_bar_bar
    # step to the new point
    take_step <- produce_step(params,est[iter],grad,mle,target,limits)
    params <- take_step$parameters
    constraint <- take_step$hit_boundary
    boundaryHits[iter] <- take_step$boundaryHits

    iter <- iter+1

    # nominal experiment at current position
    output <- Nominal.experiment(params,reps,sim,...)
    est[iter] <- output$KPI
    
    grad <- output$gradient
    
    # update bias
    bias[iter] <- target - est[iter]
  }
  #plot.bias.trace(bias,epsilon)
  results <- list(final.KPI = est[iter], theta.new = params, 
                  iterations = iter-1, boundary = constraint, 
                  pg.iterations = boundaryHits)
  return(results)
}

produce_step <- function(params,sim_est,grad,mle,target,limits){
  # Function to calculate the next step of the iterative recalibration process
  #
  # params  - current set of input parameters
  # sim_est - estimate of the current performance measure at params
  # grad    - estimated gradient of the response
  # mle     - the MLE for the input parameters
  # target  - the estimate of the unbiased performance measure
  # limits  - lower bounds on the input parameters
  
  # step to the new point
  step<-((target-sim_est)/as.numeric(t(grad)%*%grad))*grad
  + (as.numeric(t(grad)%*%(params-mle))/as.numeric(t(grad)%*%grad))*grad
  - (params-mle)

  # update parameters
  params <- params + step

  #check for feasibility
  boundaryHits <- 0
  constraint <- F
  proj_grad <- grad
  while (sum(params<limits)>0) {
    constraint <- T
    boundaryHits <- boundaryHits + 1
    #if infeasible, project to feasible and set gradient component to 0:
    correction <- rep(0,length(params))
    for (i in 1:length(params)) {
      if (params[i] <= limits[i]) {
        correction[i] <- limits[i] - params[i]
        params[i] <- max(params[i],limits[i])
        proj_grad[i] <- 0
      }
    }
    # approximate value at new position
    proj_est <- target + as.numeric(t(grad)%*%correction)

    # now add the projected gradient step:
    step<-((target-proj_est)/as.numeric(t(proj_grad)%*%proj_grad))*proj_grad

    # update parameters
    params <- params + step
    #print(c("PG step",params))
  }
  return(list(parameters=params,hit_boundary=constraint,boundaryHits = boundaryHits))
}

plot.bias.trace <- function(bias,h.width) {
  # plot of the estimated bias over the duration of the algorithm
  # include an indication of the tolerated region
  plot(1:length(bias)-1,bias,type = 'l',xlab="Iteration",ylab="Bias")
  points(c(0,length(bias)),c(0,0),type='l',col=2)
  points(c(0,length(bias)),rep(h.width,2),type='l',lty=2,col=2)
  points(c(0,length(bias)),rep(-h.width,2),type='l',lty=2,col=2)
}
