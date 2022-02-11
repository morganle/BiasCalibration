## data - file name read SAN, number of data points, number of replications

## columns of data
# g                           - the macro replication number (range from 1-10 as I had to run in parallel)
# results iterations          - how many iterations did the SQP do before stopping
# h_width                     - the half width of the 95% confidence interval around Y_bar_bar (this informs the stopping rule)
# bootstrap_bias_est          - the estimate of the bias caused by input modelling 
# Y_bar_bar.results.final.KPI - how far we end up from the target after SQP
# Y_mle                       - the output of a run of the nominal experiment at the MLE
# Y_bar_bar                   - the bootstrap estimate
# Bias_adjusted_Y             - the target we're trying to calibrate to
# results.final.KPI           - the output of a run of the nominal experiment at the calibrated input parameters
# MLE i                       - MLE of input i
# theta.new.1                 - bias calibrated input i

mean_norm_distance <- function(data) {
  distances <- c()
  col.start <- which(names(data)=="mle_1")
  col.end <- which(names(data)=="theta_1") - 1
  dimensions <- col.end + 1 - col.start
  for (i in 1:dim(data)[1]) {
    squared_dist <- 0
    for (j in col.start:col.end) {
      squared_dist <- squared_dist + (data[i,j] - data[i,j+dimensions])^2
    }
    distances[i] <- sqrt(squared_dist)
  }
  return(mean(distances))
}

total_sim_reps <- function(data_row,n0){
  N <- n0
  # number of replications for the bootstrapping
  N <- N + data_row[,4]*n0
  # count the additional reps from rerunning the nominal
  if(data_row[,2]>n0){
    N <- N + data_row[,2]
  }
  # reps from the iteration
  N <- N + data_row[,2]*data_row[,3]
  return(N)
}

################################ FACLOC Experiment ##############


ms <- c(100)
variance_ratios <- c(100)
ns <- c(100,500)
par(mfrow=c(2,1))
for (n0 in ns) {
  filename <- paste0("facloc_experiment_",n0)
  for (m in ms) {
    for (r in variance_ratios){
      data<-read.csv(paste0(filename,'.csv'),header=TRUE)
      G<-dim(data)[1]
      
      n_bar <- mean(data$sim_reps)
      n_bar_h_width <- qnorm(0.975)*sqrt(var(data$sim_reps)/length(data$sim_reps))
      
      b_hat<-mean(data$bootstrap_bias_est)
      b_hat_h_width <- qnorm(0.975)*sqrt(var(data$bootstrap_bias_est)/length(data$bootstrap_bias_est))
      SQPs<-which(data$iterations>0)
      delta<-length(which(data$iterations>0))/G
      K_bar<-mean(data$iterations)
      m_K<-as.numeric(max(data$iterations))
      psi <- sum(data$hit_boundary=="True")/G * 100
      
      gamma<-mean((abs(data$bootstrap_bias_est[SQPs])-abs(data$final_bias[SQPs]))/abs(data$bootstrap_bias_est[SQPs]))*100
      gamma_h_width<- qnorm(0.975)*sqrt(var((abs(data$bootstrap_bias_est[SQPs])-abs(data$final_bias[SQPs]))/abs(data$bootstrap_bias_est[SQPs]))/(delta*G))*100
      tau<-mean_norm_distance(data)
      h_wid<-mean(data$h_width)
      
      hist(data$iterations,main=paste('iterations',m,r))
      hist(data$sim_reps,main=paste('new n:',m,r))
      
      # count total number of simulation reps
      total_reps <- rep(0,G)
      for(i in 1:G){
        total_reps[i] <- total_sim_reps(data[i,],n0)
      }
      hist(total_reps[data$iterations>0 & data$iterations<100])
      print(c(m,r,mean(total_reps),mean(total_reps[data$iterations>0]),qnorm(0.975)*sd(total_reps[data$iterations>0])/sqrt(delta*G)))
      
      OO<-c(m,n0,r,n_bar,n_bar_h_width,b_hat,b_hat_h_width,delta,K_bar,m_K,psi,gamma,gamma_h_width,tau,h_wid)
      #print(OO)
      write.table(t(OO),paste0(filename,"_KPI.csv"),sep=",",row.names = FALSE, col.names = FALSE, append = TRUE)
    }
  }
}
par(mfrow=c(1,1))
plot(data$bootstrap_bias_est, data$sim_reps)



############# Changing the system settings. #######
warehouses <- c(3,5,7)
m <- 100
n <- 100

for(w in warehouses){
  sim_obs <- read.csv(paste0('facloc_hists_',w,'_',n,'_',m,'.csv'))
  
  G <- (dim(sim_obs)[2]-1)/2
  
  mle_cols <- (1:79)*2
  corrected_cols <- mle_cols+1
  
  means <- apply(sim_obs,2,mean)
  
  bias <- means - means[1]
  
  mle_bias <- bias[mle_cols]
  corrected_bias <- bias[corrected_cols]
  
  step <- (round(max(mle_bias,corrected_bias),2) - round(min(mle_bias,corrected_bias),2))/10
  bins <- seq(floor(min(mle_bias,corrected_bias)/step)*step,ceiling(max(mle_bias,corrected_bias)/step)*step,step)
  
  print(paste(c("Average bias at MLE:",mean(mle_bias))))
  
  # produce a plot
  plot.name <- paste0("plots/Bias_",w,".jpeg")
  jpeg(plot.name, width = 4.8, height = 4.5, units = 'in', res = 600)
  
  col1 <- "red"
  col2 <- "blue"
  y.max <- max(hist(mle_bias,plot=F)$counts,hist(corrected_bias,plot=F)$counts)
  hist(mle_bias,bins, ylim=c(0,y.max),col=col1, main=paste0("Bias: ",w,' Warehouses'), xlab='Bias estimate')
  hist(corrected_bias, col=col2, add=T)
  hist(mle_bias,bins, col=col1, add=T)
  hist(corrected_bias, col=rgb(0,0,1,0.5), add=T)
  if(w==3){
    legend('topright',c("MLE","Recalibrated"),col=c(col1,col2),pch=15)
  }
  else{
    legend('topleft',c("MLE","Recalibrated"),col=c(col1,col2),pch=15)
  }
  
  dev.off()
}


##################
t_reps <- matrix(nrow=100,ncol=2)
for (n0 in ns) {
  j <- which(ns==n0)
  filename <- paste0("facloc_experiment_",n0)
  for (m in ms) {
    for (r in variance_ratios){
      data<-read.csv(paste0(filename,'.csv'),header=TRUE)
      G<-dim(data)[1]
      
      delta<-length(which(data$iterations>0))/G
      
      # count total number of simulation reps
      total_reps <- rep(0,G)
      for(i in 1:G){
        t_reps[i,j] <- total_sim_reps(data[i,],n0)
      }
      
    }
  }
}
plot(t_reps,xlim=c(0,5000),ylim=c(0,5000))
hist(t_reps[,1]-t_reps[,2])
