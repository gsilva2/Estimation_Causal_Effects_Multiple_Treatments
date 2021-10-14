############################################ Functions: Estimation of Causal Effects ############################################ 

#################################################################################################################################

#### Read in scripts with functions ####

#Set working directory to match folder where script functions (script1_dataset_generator_function, script2_estimation_methods_functions) are saved 

source("script1_dataset_generator_functions.R")
source("script2_estimation_methods_functions.R")

#################### Section 3: Simulating a Dataset, Using this Dataset to Estimate Causal Effects with Different Estimation Methods ####################

############################### Simulate a Dataset ###############################

#Simulate covariates + treatment indicator
#Keep most parameters to default
X_dat<-X_data_simulator(Z=3, n1=600, gamma=2, L=5, d=0.5)

#Estimate GPS matrix for X_dat (with discarding, without discarding)
gps_mat<-gps_matrix(X_dat)

#Generate full set of potential outcomes - this full set will not be used to estimate causal effects since we don't get to observe this full set
#In analysis, only Y_obs will be used, we are simply generating the full set of potential outcomes to show how it was done for the simulation
potential_out_mat<-Y_data_simulator(X_dat, m=0, theta=1)

#Compile all of the generated data to obtain a full dataset and a column with observed outcome
#Observed outcome is a function of potential outcome matrix, treatment indicator
data<-dataset_generator(X_dat, potential_out_mat, gps_mat)

#Drop columns with potential outcomes since we don't get to observe this full set
#Only keep Y_obs
#Note: names vector will differ depending on value of Z
names<-c("Y_V1", "Y_V2", "Y_V3")
col_num<-which(colnames(data)==names)
data<-data[,-col_num]

############################### Estimate Causal Effects ###############################

# The functions here estimate the ATTs (ATT_2_1, ATT_3_1) among those who received treatment 1

#Here we will estimate causal effects with discarding
#To estimate causal effects without discarding, set discard=FALSE

########## MI-Cub ##########

mi_cub_res<-mi_calculator(data, spline="cubic", M=25, discard=TRUE)
mi_cub_res

########## MI-TP ##########

mi_tp_res<-mi_calculator(data, spline="tp", M=25, discard=TRUE)
mi_tp_res

########## MI-BART ##########

mi_bart_res<-bart_calculator(data, covs="FALSE", M=25, discard=TRUE)
mi_bart_res

########## MI-BARTC ##########

mi_bartc_res<-bart_calculator(data, covs="TRUE", M=25, discard=TRUE)
mi_bartc_res
