library(kSamples)
settings <- c(5,3,7)

AD <- T

data <- read.csv(paste0('facloc_hists_',settings[1],'_100_100.csv'))
# macro-reps
M <- (dim(data)[2] - 1)/2

mle_ks <- matrix(nrow=M,ncol=length(settings))
mle_ad <- matrix(nrow=M,ncol=length(settings))
cal_ks <- matrix(nrow=M,ncol=length(settings))
cal_ad <- matrix(nrow=M,ncol=length(settings))
corr_mle_ks <- matrix(nrow=M,ncol=length(settings))
corr_mle_ad <- matrix(nrow=M,ncol=length(settings))

mle_bias <- matrix(nrow=M,ncol=length(settings))
cal_bias <- matrix(nrow=M,ncol=length(settings))
corr_mle_bias <- matrix(nrow=M,ncol=length(settings))

for( s in settings ){
  i <- which(settings==s)
  data <- read.csv(paste0('facloc_hists_',s,'_100_100.csv'))
  
  true_reps <- data[,1]
  
  for (m in 1:M){
    # perform test on mle
    mle_reps <- data[,2*m]
    test <- ks.test(mle_reps,true_reps)
    mle_ks[m,i] <- test$statistic
    
    mle_bias[m,i] <- mean(mle_reps)-mean(true_reps)
    
    # perform test on debiased output
    corr_mle_reps <- mle_reps - mle_bias[m,i]
    test <- ks.test(corr_mle_reps,true_reps)
    corr_mle_ks[m,i] <- test$statistic
    
    # perform test on recalibrated inputs
    cal_reps <- data[,2*m+1]
    test <- ks.test(cal_reps,true_reps)
    cal_ks[m,i] <- test$statistic
    
    if(AD){
      # now do AD test
      test <- ad.test(mle_reps,true_reps)
      mle_ad[m,i] <- test$ad[1,1]
      # now do AD test
      test <- ad.test(corr_mle_reps,true_reps)
      corr_mle_ad[m,i] <- test$ad[1,1]
      # perform test on recalibrated inputs
      test <- ad.test(cal_reps,true_reps)
      cal_ad[m,i] <- test$ad[1,1]
    }
    
    # calculate the biases
    mle_bias[m,i] <- mean(mle_reps)-mean(true_reps)
    corr_mle_bias[m,i] <- mean(corr_mle_reps)-mean(true_reps)
    cal_bias[m,i] <- mean(cal_reps)-mean(true_reps)
  }
}

# now worry about plots
step <- 0.05
col1 <- rgb(255, 0, 0, max = 255, alpha = 125)
col2 <- rgb(0, 0, 255, max = 255, alpha = 125)
col3 <- rgb(0, 255,255, max = 255, alpha = 125)
for( i in 1:length(settings)){
  s <- settings[i]
  # KS plots
  breaks <- seq(0,ceiling(max(mle_ks[,i],cal_ks[,i])/step)*step,step)
  # look at max frequency
  y_max <- max(hist(mle_ks[,i],breaks=breaks,plot=F)$counts,hist(cal_ks[,i],breaks=breaks,plot=F)$counts,hist(corr_mle_ks[,i],breaks=breaks,plot=F)$counts)
  
  plot.name <- paste0("KS_3hist_",s,".jpeg")
  jpeg(plot.name, width = 6, height = 5, units = 'in', res = 600)
  hist(mle_ks[,i], col=col1, main=paste(s,'Warehouses: K-S Statistic'),
           breaks=breaks,ylim=c(0,y_max),
           xlab = "Kolmogorov-Smirnov Statistic")
  
  hist(cal_ks[,i],add=T,col = col2, breaks = breaks)
  hist(corr_mle_ks[,i],add=T,col = col3, breaks = breaks)
  legend('topright',c('MLE','Bias corrected','Recalibrated'),col=c(col1,col3,col2),pch=15)
  
  dev.off()
  
  print(mean(mle_ks))
  print(mean(cal_ks))
  
  hist(mle_ks[,i]-cal_ks[,i])
  
  # AD plots
  if(AD){
    h = hist(mle_ad[,i],plot=F)
    x_max <- max(c(mle_ad[,i],cal_ad[,i]))
    
    if(max(mle_ad[,i])<2*max(cal_ad[,i]) && max(mle_ad[,i])>0.5*max(cal_ad[,i])){
      ad_step <- h$breaks[2] - h$breaks[1]
      bins <- seq(0,ceiling(x_max/ad_step)*ad_step,ad_step)
      
      y_max <- max(h$counts,hist(cal_ad[,i],plot=F,breaks=bins)$counts)
      
      hist(mle_ad[,i], col=col1, main=paste(s,'warehouses: AD stat'),
           xlab = "Anderson-Darling Statistic",ylim=c(0,y_max))
      hist(cal_ad[,i],add=T,col = col2, breaks = bins)
    }
    else{
      y_max <- max(h$counts,hist(cal_ad[,i],plot=F)$counts)
      hist(mle_ad[,i], col=col1, main=paste(s,'warehouses: AD stat'),
           xlab = "Anderson-Darling Statistic",ylim=c(0,y_max))
      hist(cal_ad[,i],add=T,col = col2)
      hist(corr_mle_ad[,i],add=T,col = col3)
    }
    
    legend('topright',c('MLE','Bias corrected','Recalibrated'),col=c(col1,col3,col2),pch=15)
    
    hist(mle_ad[,i]-cal_ad[,i])
  }
  
}

# look at changes in GoF stats as a function of change in bias
for(i in 1:length(settings)){
  
  plot(mle_bias[,i]-cal_bias[,i], mle_ks[,i]-cal_ks[,i])
  plot(abs(mle_bias[,i])-abs(cal_bias[,i]), mle_ks[,i]-cal_ks[,i])
  
  plot(mle_bias[,i]-cal_bias[,i], mle_ad[,i]-cal_ad[,i])
  plot(abs(mle_bias[,i])-abs(cal_bias[,i]), mle_ad[,i]-cal_ad[,i])
  
}