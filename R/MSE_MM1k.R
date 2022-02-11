
source("GG1ksimEdits.R")
set.seed(12345)

lam = 0.8
mu = 1
rho = lam/mu
B = 400
m = 500
G = 1000  ## macro reps

# nT - number of replications
nT<-100
# runlength - number of customers simulated
runlength<-1000


MM1k<-function(rho){
  return( rho/(1-rho) - (101*(rho)^101)/(1-rho^101))
}

eta_gradient<-function(lam,mu){
  k<-100
  grad1<- mu/((mu-lam)^2)   - ( (k+1)^2*mu^(k+1)*lam^k ) / (( mu^(k+1) - lam^(k+1) )^2)  # diff by lam
  grad2<- -lam/((mu-lam)^2) + ( (k+1)^2*mu^(k)*lam^(k+1))/ (( mu^(k+1) - lam^(k+1) )^2)  # diff by mu
  return(c(grad1,grad2))
}
IU_calc<-function(arrs,sers){
  Y<-numeric(B)
  for(i in 1:B){
    lam_star<-1/mean(sample(arrs,m,replace = TRUE))
    mu_star<-1/mean(sample(sers,m,replace = TRUE))
    #if(lam_star/mu_star<1){
    Y[i]<-MM1k(lam_star/mu_star)  # no simulation noise
    #}
    # else{
    #   Y[i] = MM1K_nominal(c(lam_star,mu_star), cap=100, nT=nT , runlength=runlength)$KPI 
    # }
  }
  return( (1/(B-1))*sum( (mean(Y) - Y[i])^2 ) )
}


## when sim response known
bias.correction <- function(start,target,limits,epsilon) {
  #print(target)
  #print(epsilon)
  
  # function that actually performs the correction
  # vector of parameters
  mle <- start
  params <- start
  
  # initial nominal experiment
  est <- MM1k(params[1]/params[2]) #c(initial_nominal$KPI)
  grad <- eta_gradient(params[1],params[2]) #initial_nominal$gradient
  constraint <- F
  boundaryHits <- c()
  
  # vector of bias
  bias <- c(target - est[1])
  
  # counter
  iter <- 1
  
  # loop whilst the bias is too big
  while (est[iter] > (target+epsilon) | est[iter] < (target-epsilon) ) {   #keep going whilst the estimate is outsideof the CI for Y_bar_bar
    # step to the new point
    step<-((target-est[iter])/as.numeric(t(grad)%*%grad))*grad+ (as.numeric(t(grad)%*%(params-mle))/as.numeric(t(grad)%*%grad))*grad- (params-mle)
    
    # update parameters
    params <- params + step
    #print(params)
    
    #check for feasibility
    boundaryHits[iter] <- 0
    proj_grad <- grad
    while (sum(params<limits)>0) {
      constraint <- T
      boundaryHits[iter] <- boundaryHits[iter] + 1
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
    
    iter <- iter+1
    
    # nominal experiment at current position
    #output <- MM1k(params[1]/params[2])#Nominal.experiment(params,sim,...)

    est[iter] <- MM1k(params[1]/params[2]) #output$KPI
    grad <- eta_gradient(params[1],params[2]) #output$gradient

   # update bias
    bias[iter] <- target - est[iter]
    if(iter > 30){
      iter<-10000000
      break
    }
  }
  #plot.bias.trace(bias,epsilon)
  results <- list(final.KPI = est[iter], theta.new = params, iterations = iter-1, boundary = constraint, pg.iterations = boundaryHits)
  return(results)
}

## parametric bootstrap
out = numeric(G)
rho_hat = numeric(G)
lam_hat = numeric(G)
mu_hat = numeric(G)
b_est = numeric(G)
IU_est = numeric(G)
MSE_est = numeric(G)


eta_c = MM1k(lam/mu)
for( i in 1:G){
  arrs<-rexp(m,lam)
  sers<-rexp(m,mu)
  lam_hat[i]<-1/mean(arrs)
  mu_hat[i]<-1/mean(sers)
  rho_hat[i] = lam_hat[i]/mu_hat[i]
  out[i]<-MM1k(rho_hat[i])
  b_est[i]<- abs(out[i] - eta_c)
  IU_est[i]<-IU_calc(arrs,sers)
  MSE_est[i]<- IU_est[i] + b_est[i]^2
}
h_width = 0.025*MM1k(rho)


## bias correction gets stuck when the mle is really far from the truth.
## This is because when the truth is known we don't have a sensible CI 
## from which to choose h_width.


#inds<-which(rho_hat > 0.82 & rho_hat  < 1.05) # & (rho_h < 0.97 | rho_h > 1.03))

#lam_hat<-lam_hat[inds]
#mu_hat<-mu_hat[inds]
#rho_hat<-rho_hat[inds]#
#b_est<-b_est[inds]
#IU_est<-IU_est[inds]#
#MSE_est<-MSE_est[inds]
#G<-length(rho_hat)

# before bias correction

plot(MSE_est)
points(IU_est,col="blue")
points(b_est,col="green")


#pos_bias_inds = which(rho_hat > 0.8)
#neg_bias_inds = which(rho_hat < 0.8)
#b_pos<-mean(abs( out[pos_bias_inds] - rep(eta_c,length(pos_bias_inds))))
#IU_pos<-var(out[pos_bias_inds])
#b_neg<-mean(abs( out[neg_bias_inds] - rep(eta_c,length(neg_bias_inds))))
#IU_neg<-var(out[neg_bias_inds])

## now with the bias correction

corrected_params<-matrix(rep(0,2*G),ncol=2)
iterations<-numeric(G)
remaining_bias<-numeric(G)
for(i in 1:G){
  print(i)
bias_adjusted_Y <- MM1k(rho)
results<-bias.correction(start = c(lam_hat[i],mu_hat[i]), target = bias_adjusted_Y, limits = c(0.01,0.01), epsilon = h_width)
iterations[i]<-results$iterations
corrected_params[i,]<- results$theta.new
remaining_bias[i]<- abs(results$final.KPI - MM1k(rho))
}

G<-500
remaining_bias<-remaining_bias[1:G]
iterations<-iterations[1:G]

mean(b_est)
length(which(iterations>0))/G
mean(iterations)
max(iterations)
mean((b_est[1:G] - remaining_bias)/ b_est[1:G])


lam_corrected<-corrected_params[(1:G),1]
mu_corrected<-corrected_params[(1:G),2]

mean(abs(lam_corrected - lam_hat[1:G]))
mean(abs(mu_corrected - mu_hat[1:G]))
h_width

rho_corrected<-lam_corrected/mu_corrected

plot(rho_hat)
abline(h=lam,col="green")
points(rho_corrected, col="blue")

out2<-numeric(G)
b_corr_est = numeric(G)
IU_corr_est = numeric(G)
MSE_corr_est = numeric(G)
for( i in 1:G){
  out2[i]<-MM1k(rho_corrected[i])
  b_corr_est[i]<- abs(out2[i] - eta_c)
  arrs<-rexp(m,lam_corrected[i])
  sers<-rexp(m,mu_corrected[i])
  IU_corr_est[i]<-IU_calc(arrs,sers)
  MSE_corr_est[i]<- IU_corr_est[i] + b_corr_est[i]^2
}

par(mfrow=c(1,2))

plot(b_est,col="green",ylab="MSE, IU and Bias",main="Before bias-correction",ylim=c(0,max(MSE_est)))
points(IU_est,pch=18,col="cyan4")
points(b_est,pch=18,col="darkolivegreen2")
points(MSE_est,pch=20,col="purple")

plot(b_corr_est,col="green",ylab="MSE, IU and Bias",main="After bias-correction",ylim=c(0,max(MSE_corr_est)))

points(IU_corr_est,pch=18,col="cyan4")#"darkolivegreen2")
points(b_corr_est,pch=18,col="darkolivegreen2")#="cyan4")
points(MSE_corr_est,pch=20,col="purple")

par(mfrow=c(1,3))

s_MSE_est<-sort(MSE_est)
s_MSE_corr_est<-sort(MSE_corr_est)

cap<-floor(G*0.9)
inds1<-which(MSE_est<s_MSE_est[cap])
inds2<-which(MSE_corr_est<s_MSE_corr_est[cap])

boxplot(MSE_est[inds1],MSE_corr_est,main="MSE",names=c("Before", "After"),col=c("purple","purple"))   

boxplot(b_est[inds1],b_corr_est,main="Bias",names=c("Before", "After"),col=c("green","green"))

boxplot(IU_est[inds1],IU_corr_est,main="IU",names=c("Before", "After"),col=c("blue","blue"))

length(which( b_est - b_corr_est >0))/G
length(which( IU_est - IU_corr_est >0))/G
length(which( MSE_est - MSE_corr_est >0))/G

max(MSE_est)

#which(rho_hat>lam) == which( IU_est - IU_corr_est >0 )

#s_MSE_est<-sort(MSE_est)[1:950]
#s_MSE_corr_est<-sort(MSE_corr_est)[1:950]

#r_inds<-which(MSE_est>5)
#boxplot(MSE_est[-r_inds],IU_est[-r_inds],b_est[-r_inds],col=c("black","blue","green"))
#boxplot(MSE_corr_est[-r_inds],IU_corr_est[-r_inds],b_corr_est[-r_inds],col=c("black","blue","green"))

#b_corr_pos<-mean(abs( out2[pos_bias_inds] - rep(eta_c,length(pos_bias_inds))))
#IU_corr_pos<-var(out2[pos_bias_inds])
#b_corr_neg<-mean(abs( out2[neg_bias_inds] - rep(eta_c,length(neg_bias_inds))))
#IU_corr_neg<-var(out2[neg_bias_inds])


