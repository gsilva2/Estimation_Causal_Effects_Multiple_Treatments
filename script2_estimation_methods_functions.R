############################################ Functions: Estimation of Causal Effects ############################################ 

#################################################################################################################################

#Loading required packages
library(varhandle)
library(Matching)
library(data.table)
library(rlist)
library(nnet)
library(gtools)
library(survey)
library(mgcv)
library(MASS)
library(BART)

#################################################################################################################################


#################### Section 2: Functions to Estimate Causal Effects for Different Estimation Methods ####################

# All of these functions take as input a dataset generated using the functions in Script 1
# We assume Z=3 for these functions

# The functions here estimate the ATTs (ATT_2_1, ATT_3_1) among those who received treatment 1


########## Multiple Imputation Functions (Based on GAMs) ##########

#spline="tp" is MI-TP
#spline="cubic" is MI-Cub
#M is number of imputed datasets
mi_calculator<-function(dataset, spline="tp", M=25, discard=TRUE){
  
  if(discard==TRUE){
    dataset<-dataset[complete.cases(dataset), ]
    dataset<-as.data.frame(dataset)
    y<-which(names(dataset)=="gps_1")
    z<-which(names(dataset)=="gps_3")
    dataset<-dataset[ ,-(y:z)]
    
    w<-which(names(dataset)=="gps_disc_gps_1")
    names(dataset)[w:(w+2)]<-c("gps_1", "gps_2", "gps_3")
  }
  
  
  if(discard==FALSE){
    dataset<-dataset
    dataset<-as.data.frame(dataset)
  }  
  
  
  ntreat<-length(table(dataset$t_z))  
  
  ########## Model Fitting Procedure ########## 
  #So now we want to use the data matrix and the observed outcome, Y, to generate "ntreat" models
  #Predictions are saved in an array with M rectangular matrices, each with # of rows equal to # of patients and # of columns equal to # of treatments 
  
  
  predictions<-array(NA, dim=c(sum(dataset$t_z==1), ntreat, M))
  for(j in 2:ntreat){
    data_fit<-dataset[dataset$t_z==j, ]
    #Number of covariates
    a<-which(names(data_fit)=="t_z")-1
    
    
    #Fit model using observed data
    #This model includes all covariate entries-2 and ntreat-1 generalized propensity scores
    # We consider different types of splines here
    if(spline=="tp2"){
      if(a==5){
        
        model<-mgcv::gam(Y_obs~(V1+V2+V3+ts(gps_2, gps_3)), data=data_fit)
      }
      if(a==10){
        
        model<-mgcv::gam(Y_obs~(V1+V2+V3+V4+V5+V6+V7+V8+ts(gps_2, gps_3)), data=data_fit)
      }
    }
    
    if(spline=="cubic"){
      if(a==5){
        
        model<-mgcv::gam(Y_obs~(V1+V2+V3+s(gps_2, bs="cr", k=8)+s(gps_3, bs="cr", k=8)), data=data_fit)
      }
      if(a==10){
        
        model<-mgcv::gam(Y_obs~(V1+V2+V3+V4+V5+V6+V7+V8+s(gps_2, bs="cr", k=8)+s(gps_3, bs="cr", k=8)), data=data_fit)
      }
    }
    
    if(spline=="tp"){
      if(a==5){
        
        model<-mgcv::gam(Y_obs~(V1+V2+V3+s(gps_2, bs="tp")+s(gps_3, bs="tp")), data=data_fit)
      }
      if(a==10){
        
       model<-mgcv::gam(Y_obs~(V1+V2+V3+V4+V5+V6+V7+V8+s(gps_2, bs="tp")+s(gps_3, bs="tp")), data=data_fit)
        
      }
    }
    
    #Multiply Impute missing potential outcomes 
    data_pred<-dataset[dataset$t_z==1, ]
    new<-predict(model, newdata=data_pred, type="lpmatrix")
    betas<-mvrnorm(n=M, coef(model), model$Vp)
    
    y_new<-c()
    for(m in 1:M){
      y_new<-cbind(y_new, new%*%betas[m, ])
    }
    
    predictions[ ,j , ]<-y_new
    
  }
  
  #Potential outcomes for treatment =1 have been observed for the group receiving the base treatment, treatment 1
  predictions[ , 1, ]<-dataset$Y_obs[dataset$t_z==1]
  
  #Calculate differences in means and SEs
  res<-matrix(NA, nrow=M, ncol=4)
  for(m in 1:M){
    est<-predictions[, , m]
    ######NOTE: Change this if Z !=3
  
    diff_1<-est[,2]-est[,1]
    res_1_est<-mean(diff_1)
    res_1_se<-sd(diff_1)/sqrt(length(est[,1]))
    
    diff_2<-est[,3]-est[,1]
    res_2_est<-mean(diff_2)
    res_2_se<-sd(diff_2)/sqrt(length(est[,1]))
    
    res[m, ]<-c(res_1_est, res_2_est, res_1_se, res_2_se)
  }
  
  #Output of res is SE - square third and fourth columns to get within-imputation variance
  res[,3:4]<-res[,3:4]^2
  
  #Combining M estimates and within-imputation variances to get imputation estimate of ATT + SEs
  att_est<-colMeans(res[,1:2])
  att_se<-sqrt(colMeans(res[,3:4])+(1+(1/M))*c(var(res[,1]), var(res[,2])))
  #Obtained from Yuan (2016) but also available in Rubin (1987)
  df<-(M-1)*(1+(colMeans(res[,3:4])/((1+(1/M))*c(var(res[,1]), var(res[,2])))))^2
  
  
  
  output<-matrix(NA, nrow=2, ncol=4)
  output[,1]<-att_est
  output[,2]<-att_se
  output[1, 3:4]<-c(att_est[1]-qt(0.975, df[1])*att_se[1], att_est[1]+qt(0.975, df[1])*att_se[1])
  output[2, 3:4]<-c(att_est[2]-qt(0.975, df[2])*att_se[2], att_est[2]+qt(0.975, df[2])*att_se[2])
  
  output<-as.data.frame(output)
  
  colnames(output)<-c("Estimate", "SE", "Low_95", "High_95")
  rownames(output)<-c("ATT_2_1", "ATT_3_1")
  
  
  return(output)
  
}


########## Multiple Imputation Functions (Based on BARTs) ##########

#covs=FALSE is MI-BART
#covs=TRUE is MI-BARTC
#M is number of imputed datasets
bart_calculator<-function(dataset, covs=TRUE, M=25, discard=TRUE){
  
  if(discard==TRUE){
    dataset<-dataset[complete.cases(dataset), ]
    dataset<-as.data.frame(dataset)
    y<-which(names(dataset)=="gps_1")
    z<-which(names(dataset)=="gps_3")
    dataset<-dataset[ ,-(y:z)]
    
    w<-which(names(dataset)=="gps_disc_gps_1")
    names(dataset)[w:(w+2)]<-c("gps_1", "gps_2", "gps_3")
  }
  
  
  if(discard==FALSE){
    dataset<-dataset
    dataset<-as.data.frame(dataset)
  }  
  
  
  ntreat<-length(table(dataset$t_z))  
  
  ########## Model Fitting Procedure ########## 
  #So now we want to use the data matrix and the observed outcome, Y, to generate "ntreat" models
  #Predictions are saved in an array with M rectangular matrices, each with # of rows equal to # of patients and # of columns equal to # of treatments 
  
  predictions<-array(NA, dim=c(sum(dataset$t_z==1), ntreat, M))
  for(j in 2:ntreat){
    data_fit<-dataset[dataset$t_z==j, ]
    #Number of covariates
    a<-which(names(data_fit)=="t_z")-1
    
    
    #Fit model using observed data
    #This model includes all covariate entries-2 and ntreat-1 generalized propensity scores
    
    if(covs==TRUE){
      x<-as.matrix(data_fit[ ,1:(a-2)])
      x<-cbind(x, data_fit$gps_2, data_fit$gps_3)
    }
    else{x<-as.matrix(cbind(data_fit$gps_2, data_fit$gps_3))}
    y<-data_fit$Y_obs
    
    hi<-capture.output({bar<-wbart(x, y, ndpost=M)})
    
    
    #Multiply Impute missing potential outcomes 
    data_pred<-dataset[dataset$t_z==1, ]
    
    if(covs==TRUE){
      x_pred<-as.matrix(data_pred[ ,1:(a-2)])
      x_pred<-cbind(x_pred, data_pred$gps_2, data_pred$gps_3)
    }
    else{x_pred<-as.matrix(cbind(data_pred$gps_2, data_pred$gps_3))}
    
    hi2<-capture.output({pred<-predict(bar, x_pred)})
    
    predictions[ ,j , ]<-t(pred)
    
  }
  
  #Potential outcomes for treatment =1 have been observed for the group receiving the base treatment, treatment 1
  predictions[ , 1, ]<-dataset$Y_obs[dataset$t_z==1]
  
  #Calculate differences in means and SEs
  res<-matrix(NA, nrow=M, ncol=4)
  for(m in 1:M){
    est<-predictions[, , m]
    ######NOTE: Change this if Z !=3
    
    diff_1<-est[,2]-est[,1]
    res_1_est<-mean(diff_1)
    res_1_se<-sd(diff_1)/sqrt(length(est[,1]))
    
    diff_2<-est[,3]-est[,1]
    res_2_est<-mean(diff_2)
    res_2_se<-sd(diff_2)/sqrt(length(est[,1]))
    
    res[m, ]<-c(res_1_est, res_2_est, res_1_se, res_2_se)
  }
  
  
  #Output of res is SE - square third and fourth columns to get within-imputation variance
  res[,3:4]<-res[,3:4]^2
  
  #Combining M estimates and within-imputation variances to get imputation estimate of ATT + SEs
  att_est<-colMeans(res[,1:2])
  att_se<-sqrt(colMeans(res[,3:4])+(1+(1/M))*c(var(res[,1]), var(res[,2])))
  #obtained from Yuan (2016) but also available in Rubin (1987)
  df<-(M-1)*(1+(colMeans(res[,3:4])/((1+(1/M))*c(var(res[,1]), var(res[,2])))))^2
  
  
  
  output<-matrix(NA, nrow=2, ncol=4)
  output[,1]<-att_est
  output[,2]<-att_se
  output[1, 3:4]<-c(att_est[1]-qt(0.975, df[1])*att_se[1], att_est[1]+qt(0.975, df[1])*att_se[1])
  output[2, 3:4]<-c(att_est[2]-qt(0.975, df[2])*att_se[2], att_est[2]+qt(0.975, df[2])*att_se[2])
  
  output<-as.data.frame(output)
  
  colnames(output)<-c("Estimate", "SE", "Low_95", "High_95")
  rownames(output)<-c("ATT_2_1", "ATT_3_1")
  
  return(output)
  
}

