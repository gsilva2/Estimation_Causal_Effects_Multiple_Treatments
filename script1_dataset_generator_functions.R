############################################ Functions: Estimation of Causal Effects ############################################ 

#################################################################################################################################

#Loading required packages
library(sn)
library(nnet)

#################################################################################################################################


#################### Section 1: Functions to Generate Dataset, Propensity Score Matrix, Potential Outcomes, Observed Outomes ####################

#These functions use as inputs the same simulation parameters as in the manuscript. Please refer to the manuscript for more details


############################### X_data simulator function ###############################

#This function will simulate the covariate data under various simulation configurations
#Default arguments are given, these can easily be changed
#Function returns simulated X data and treatment indicator, t_z

#Z = number of treatments -> input can either be 3, 5, 10
#n1 = number of individuals in treatment 1
#gamma = controls number of people in each treatment group
#When Z=3, n2=gamma*n1,  n3=gamma^2*n1
#When Z=5, n2=n4=gamma*n1,  n3=n5=gamma^2*n1
#When Z=10, n6=n1, n2=n4=n7=n9=gamma*n1,  n3=n5= n8=n10=gamma^2*n1

#Each subject has an L-dimensional vector X_i of covariates
#X_i varies depending on treatment received
#X_i|Z_i=z is multivariate skewed t distribution with degrees of freedom df, location=mu_t, Omega=Sigma_t, skewness=eta

X_data_simulator<-function(Z=3, n1=300, gamma=1, L=5, d=0, df=Inf, eta=0, lambda=0, sigmasq2=1, sigmasq3=1){
  
  #Let's simulate the vector of counts for each treatment group
  n<-rep(NA, times=Z)
  n[1]<-n1
  for(j in 2:Z){
    if(j%in%c(2, 4, 7, 9)){n[j]=gamma*n1}
    if(j%in%c(3, 5, 8, 10)){n[j]=gamma^2*n1}
    if(j%in%c(6)){n[j]=n1}
  }
  
  #Treatment indicator vector, t_z
  if(Z==3){t_z<-c(rep(1, n[1]), rep(2, n[2]), rep(3, n[3]))}
  if(Z==5){t_z<-c(rep(1, n[1]), rep(2, n[2]), rep(3, n[3]), rep(4, n[4]), rep(5, n[5]))}
  if(Z==10){t_z<-c(rep(1, n[1]), rep(2, n[2]), rep(3, n[3]), rep(4, n[4]), rep(5, n[5]), rep(6, n[6]), rep(7, n[7]), rep(8, n[8]), rep(9, n[9]), rep(10, n[10]))}
  
  
  #Generating mu_z vectors
  oneL<-as.matrix(rep(1, L))
  #d_mat contains d_t vectors in its rows
  d_mat<-diag(Z)
  
  if(Z==3){mu<-d*sqrt((1+sigmasq2+sigmasq3)/3)}
  if(Z==5){mu<-d*sqrt((1+2*sigmasq2+2*sigmasq3)/5)}
  if(Z==10){mu<-d*sqrt((1+9)/10)}
  
  d_mat<-mu*d_mat
  #Each row of mu_mat represents a mu_t
  mu_mat<-matrix(NA, nrow=Z, ncol=L)
  for(j in 1:Z){
    mu_mat[j, ]<-kronecker(oneL, d_mat[, j])[1:L]
  }
  
  #LXL Sigma matrix - off diagonal entries are lambda for all 
  
  Sigma_I<-diag(L)
  Sigma_A<-matrix(lambda, nrow=L, ncol=L)
  diag(Sigma_A)<-1
  Sigma_B<-matrix(lambda, nrow=L, ncol=L)
  diag(Sigma_B)<-sigmasq2
  Sigma_C<-matrix(lambda, nrow=L, ncol=L)
  diag(Sigma_C)<-sigmasq3
  
  
  Sigma_array<-array(data=NA, dim=c(L, L, Z))
  for(j in 1:Z){
    if(Z%in%c(3, 5)){
      if(j==1){Sigma_array[ , , j]<-Sigma_A}
      if(j%in%c(2, 4)){Sigma_array[ , , j]<-Sigma_B}
      if(j%in%c(3, 5)){Sigma_array[ , , j]<-Sigma_C}
    }
    else{Sigma_array[ , , j]<-1*Sigma_I}
  }
  
  #L covariates and N=sum(n) individuals in the dataset
  #So X has N rows and L columns
  #Now we generate X_i's for everyone in our dataset
  X<-matrix(NA, nrow=length(t_z), ncol=L)
  
  for(j in 1:Z){
    if(j==1){X[1:n[j], ]<-rmst(n=n[j], xi=mu_mat[j, ], Omega=Sigma_array[ , , j], alpha=rep(eta, L), nu=df)}
    if(j!=1){
      prev<-sum(n[1:(j-1)])
      X[(prev+1):(prev+n[j]), ]<-rmst(n=n[j], xi=mu_mat[j, ], Omega=Sigma_array[ , , j], alpha=rep(eta, L), nu=df)}
  }
  
  X<-as.data.frame(cbind(X, t_z))
  X$t_z<-as.factor(X$t_z)
  return(X)
  
}


############################### gps_matrix function ###############################

#This function uses the dataset generated by X_data_simulator and uses it to estimate the vector of GPS for each unit
#This functions then discards units using the discarding rule in Gutman & Lopez and detailed in the manuscript
#Then GPS is reestimated for all retained units
#This function returns a list of two sets of GPS matrices (one with discarding, the other without)
#The first element of the list [the GPS matrix before discarding] has an indicator function (ind) denoting whether the observation was discarded or kept


gps_matrix<-function(data_X){
  ntreat<-length(table(data_X$t_z))
  
  
  #### Original GPS estimation -> before discarding
  #Generate propensity score using multinomial logistic regression model with al second-level interactions
  ps_mod<-nnet::multinom(t_z~.*., data=data_X, trace=FALSE)
  #Now we have a data frame of propensity scores - each row represents a person
  prop_scores<-fitted(ps_mod)
  
  ind_rule<-c()
  
  #### Discarding step - discard units with PS outside of range suggested in Gutman & Lopez paper
  disc_rule<-matrix(NA, nrow=ntreat, ncol=2)
  for(i in 1:ntreat){
    sub<-prop_scores[,i]
    disc_rule[i, ]<-c(max(tapply(sub, data_X$t_z, min)), min(tapply(sub, data_X$t_z, max)))
  }
  
  for(i in 1:dim(data_X)[1]){
    vec<-rep(NA, ntreat)
    for(j in 1:ntreat){
      ind<-prop_scores[i,j]>disc_rule[j, 1]&prop_scores[i,j]<disc_rule[j, 2]
      vec[j]<-ind
    }
    if(sum(vec)==ntreat){ind_rule[i]<-1}else{ind_rule[i]<-0}
  }
  
  data_X_keep<-data_X[ind_rule==1, ]
  
  #### Fit GPS model again -> after discarding
  ps_mod<-nnet::multinom(t_z~.*., data=data_X_keep, trace=FALSE)
  #Now we have an updated data frame of propensity scores - each row represents a person
  prop_scores_keep<-fitted(ps_mod)
  prop_scores_2<-matrix(NA, nrow=dim(data_X)[1], ncol=ntreat)
  for(i in 1:dim(data_X)[1]){
    row_name<-rownames(data_X)[i]
    if(row_name%in%rownames(prop_scores_keep)){
      j<-which(rownames(prop_scores_keep)==row_name)
      prop_scores_2[i, ]<-prop_scores_keep[j, ]}
  }
  
  prop_scores<-as.data.frame(prop_scores)
  prop_scores_2<-as.data.frame(prop_scores_2)
  
  
  colnames(prop_scores)<-paste("gps_", colnames(prop_scores), sep="")
  colnames(prop_scores_2)<-paste("gps_disc_", colnames(prop_scores), sep="")
  
  prop_scores$ind<-ind_rule
  return(list(prop_scores,prop_scores_2))
}


############################### Y_data simulator function ###############################

#This function will take the dataset generated by X_data_simulator and uses it to simulate the full set of potential outcomes
#Function returns a matrix with Z columns and the number of rows corresponds to the number of individuals in the dataset

#The simulation factors unknown to the investigator that can be modified while simulating the data: theta and the form m(x) for the relationship between X and Y
#For this demonstration, we assume m=0
#Y_i(t)=(X_i)beta_t+e_{t, i}
#e_{t, i} are N(0, 1)
#Beta is simulated as Uniform(-theta, theta)


Y_data_simulator<-function(X_simulator, m=0, theta=1){
  
  #Generating a matrix of betas for the different treatment groups
  #each row represents betas for a different treatment group
  #each column represents betas for a different covariate (total # of covariates: L)
  n<-which(names(X_simulator)=="t_z")
  X_data<-as.matrix(X_simulator[ ,-n])
  treatment_ind<-X_simulator[ ,n]
  
  
  
  Z<-length(table(treatment_ind))
  
  beta<-runif(Z*dim(X_data)[2], min=-theta, max=theta)
  beta<-matrix(beta, nrow=Z, ncol=dim(X_data)[2])
  
  
  #Matrix of all potential outcomes - we know all the entries in this
  Y_mat<-matrix(NA, nrow=length(treatment_ind), ncol=Z)
  for(i in 1:length(treatment_ind)){
    for(j in 1:Z){
      prod<-X_data[i, ]%*%beta[j, ]
      
      if(m==0){Y_mat[i, j]<-prod+rnorm(1)}
      
    }
  }
  
  Y_mat<-as.data.frame(Y_mat)
  colnames(Y_mat)<-paste("Y_", colnames(Y_mat), sep="")
  return(Y_mat)
}



############################### dataset_generator function ###############################

#This function inputs X_data_simulator, Y_data_simulator and results from gps_matrix to return full data matrix
#Full data matrix has all simulated and estimated components, as well as a Y_obs column

dataset_generator<-function(X_mat, Y_mat, gps){
  dataset<-as.data.frame(cbind(X_mat, Y_mat, gps[[1]], gps[[2]]))
  
  n<-which(names(X_mat)=="t_z")
  treatment_ind<-X_mat[ ,n]
  
  Y_obs<-c()
  
  for(i in 1:dim(X_mat)[1]){
    Y_obs[i]<-Y_mat[i, treatment_ind[i]]
  }
  
  
  dataset$Y_obs<-Y_obs
  return(dataset)
}

#################################################################################################################################
