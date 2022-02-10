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
  col.start <- which(names(data)=="mle.1")
  col.end <- which(names(data)=="theta.1") - 1
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
  if(data_row[,2]>n0){
    N <- N + data_row[,2]
  }
  N <- N + data_row[,2]*data_row[,3]
  return(N)
}


ms <- c(50,500,5000)
variance_ratios <- c(100)
ns <- c(1000,100)
par(mfrow=c(2,1))
for (n0 in ns) {
  for (m in ms) {
    for (r in variance_ratios){
      data<-read.csv(paste0("SAN_",m,"_",n0,"_",r,".csv"),header=TRUE)
      G<-dim(data)[1]
      
      n_bar <- mean(data$sim_reps)
      n_bar_h_width <- 1.96*sqrt(var(data$sim_reps)/length(data$sim_reps))
    
      b_hat<-mean(data$bootstrap_bias_est)
      b_hat_h_width <- 1.96*sqrt(var(data$bootstrap_bias_est)/length(data$bootstrap_bias_est))
      SQPs<-which(data$iterations>0)
      delta<-length(which(data$iterations>0))/G
      K_bar<-mean(data$iterations)
      m_K<-as.numeric(max(data$iterations))
      psi <- sum(data$hit_boundary)/G * 100
      
      gamma<-mean((abs(data$bootstrap_bias_est[SQPs])-abs(data$final_bias[SQPs]))/abs(data$bootstrap_bias_est[SQPs]))*100
      gamma_h_width<- 1.95*sqrt(var((abs(data$bootstrap_bias_est[SQPs])-abs(data$final_bias[SQPs]))/abs(data$bootstrap_bias_est[SQPs]))/G)*100
      tau<-mean_norm_distance(data)
      h_wid<-mean(data$h_width)
      
      hist(data$iterations,main=paste('iterations',m,r))
      hist(data$sim_reps,main=paste('new n:',m,r))
      
      # count total number of simulation reps
      total_reps <- rep(0,G)
      for(i in 1:G){
        total_reps[i] <- total_sim_reps(data[i,],n0)
      }
      print(c(m,r,sum(total_reps)))
      
      OO<-c(m,n0,r,n_bar,n_bar_h_width,b_hat,b_hat_h_width,delta,K_bar,m_K,psi,gamma,gamma_h_width,tau,h_wid)
      write.table(t(OO),"SAN_KPIs_updated.csv",sep=",",row.names = FALSE, col.names = FALSE, append = TRUE)
    }
  }
}

for (m in ms) {
  for (r in variance_ratios){
    data<-read.csv(paste0("SAN_",m,"_50_",r,".csv"),header=TRUE)
    G<-dim(data)[1]
    
    b_hat<-mean(data$bootstrap_bias_est)
    b_hat_h_width <- 1.96*sqrt(var(data$bootstrap_bias_est)/length(data$bootstrap_bias_est))
    SQPs<-which(data$iterations>0)
    delta<-length(which(data$iterations>0))/G
    K_bar<-mean(data$iterations)
    m_K<-as.numeric(max(data$iterations))
    psi <- sum(data$hit_boundary)/G * 100
    
    gamma<-mean((abs(data$bootstrap_bias_est[SQPs])-abs(data$final_bias[SQPs]))/abs(data$bootstrap_bias_est[SQPs]))*100
    gamma_h_width<- 1.95*sqrt(var((abs(data$bootstrap_bias_est[SQPs])-abs(data$final_bias[SQPs]))/abs(data$bootstrap_bias_est[SQPs]))/G)*100
    tau<-mean_norm_distance(data)
    h_wid<-mean(data$h_width)
    
    hist(data$sim_reps,main=paste(m,r))
    
    OO<-c(m,r,b_hat,b_hat_h_width,delta,K_bar,m_K,psi,gamma,gamma_h_width,tau,h_wid)
    write.table(t(OO),"SAN_KPIs_updated.csv",sep=",",row.names = FALSE, col.names = FALSE, append = TRUE)
  }
}



# late
for (n0 in ns) {
  for (m in ms) {
    for (late in c(80,70)) {
      data<-read.csv(paste0("SAN_late_",late,"_",m,"_1000.csv"),header=TRUE)
      G<-dim(data)[1]
      
      n_bar <- mean(data$sim_reps)
      n_bar_h_width <- 1.96*sqrt(var(data$sim_reps)/length(data$sim_reps))
      
      b_hat<-mean(data$bootstrap_bias_est)
      b_hat_h_width <- 1.96*sqrt(var(data$bootstrap_bias_est)/length(data$bootstrap_bias_est))
      SQPs<-which(data$iterations>0)
      delta<-length(which(data$iterations>0))/G
      K_bar<-mean(data$iterations)
      m_K<-as.numeric(max(data$iterations))
      psi <- sum(data$hit_boundary)/G * 100
      
      gamma<-mean((abs(data$bootstrap_bias_est[SQPs])-abs(data$final_bias[SQPs]))/abs(data$bootstrap_bias_est[SQPs]))*100
      gamma_h_width<- 1.96*sqrt(var((abs(data$bootstrap_bias_est[SQPs])-abs(data$final_bias[SQPs]))/abs(data$bootstrap_bias_est[SQPs]))/G)*100
      tau<-mean_norm_distance(data)
      h_wid<-mean(data$h_width)
      
      OO<-c(m,late,n0,n_bar,n_bar_h_width,b_hat,b_hat_h_width,delta,K_bar,m_K,psi,gamma,gamma_h_width,tau,h_wid)
      
      write.table(t(OO),"SAN_late_KPIs_updated.csv",sep=",",row.names = FALSE, col.names = FALSE, append = TRUE)
    }
  }
}

