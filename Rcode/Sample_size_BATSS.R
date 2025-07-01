########################################################################################################
########################################################################################################
# 8888888 .d8888b. 88888888888 888     888                                     
#   888  d88P  Y88b    888     888     888                                     
#   888  888    888    888     888     888                                     
#   888  888           888     888     888                                     
#   888  888           888     888     888                                     
#   888  888    888    888     888     888                                     
#   888  Y88b  d88P    888     Y88b. .d88P                                     
# 8888888 "Y8888P"     888      "Y88888P"                                      
# 
# 8888888b.     d8888 888b    888 88888888888 888    888 8888888888 8888888b.  
# 888   Y88b   d88888 8888b   888     888     888    888 888        888   Y88b 
# 888    888  d88P888 88888b  888     888     888    888 888        888    888 
# 888   d88P d88P 888 888Y88b 888     888     8888888888 8888888    888   d88P 
# 8888888P" d88P  888 888 Y88b888     888     888    888 888        8888888P"  
# 888      d88P   888 888  Y88888     888     888    888 888        888 T88b   
# 888     d8888888888 888   Y8888     888     888    888 888        888  T88b  
# 888    d88P     888 888    Y888     888     888    888 8888888888 888   T88b 
########################################################################################################
########################################################################################################
########################################################################################################


# File title: PANTHER_in_BATSS.R 
# 
# Input: N/A
# 
# Output: summary information of PANTHER study - tbc if any files saved
# 
# Purpose: Test a trial of PANTHER style to see if can be studied the sample size calculation of the trial.


# ------------------ Author ---------------- Date ---------------- Update:
#            | Manel Benlloch Guajardo | 01MAY2025 |   Initial version of the file.

########################################################################################################
########################################################################################################
########################################################################################################

library(BATSS)
library(tidyr)
library(magrittr)
library(dplyr)
library(INLA)
#  Location of package : "C:/Users/nb221/AppData/Local/R/win-library/4.4"

#######################################################################################################
################### Study - preloading
################### Set the Probability file with all info         ####################################
################### Probability distribution of incidence          ####################################

setwd(paste0(rstudioapi::getSourceEditorContext()$path,"/.."))


#  Distributions in columns, each row corresponding to an outcome value.
primOutDist_panth<- read.csv(paste0(getwd(),"/../","excel_distributions/VFDdistributions_logodds.csv"),header=TRUE) 

#######################################################################################################

# END of setup - Start of the coding of the trial:

# Simulation
#BATSS using INLA's default normal priors, N(0,1000)

colnames(primOutDist_panth) #select the most appropiate ones

multinom_rand_dset1<- data.frame(primOutDist_panth[,6]) # hyper Control
multinom_rand_dset2<- data.frame(primOutDist_panth[,15]) # hyper Treatment
Trials<-25

source(paste0(getwd(),"/batss_glm_breakdown.R"))

scenario2 = batss.glm.pom(   
    model           = y ~ treatment,
    var             = list(y = multinomial_random,
                           treatment = treatalloc.fun),
    var.control     = list(y = list(size= 1)),
    family          = "pom",
    link            = 'identity',
    beta            = list(multinom_rand_dset1,# this is the control
                           multinom_rand_dset2 # this is the treatment 
                           ), 
    which           = c(2),   # Select which groups are treatments
    R               = Trials,
    control.fixed = list(mean = list( treatment = 0), prec = 0.1),
    control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, config = TRUE),
    alternative     = c("greater")
    , # One or two sided hypothesis
    map_probabilities = TRUE, # This variable will apply new maps to the days of
                              # day free support- estimation.
                             
        # RAR:
        # RAR option not used as we have uniform sampling 
        #   - continuously sample to Max sample size/ Futitility/ Efficacy
    #RAR             = prob.trippa,
    #RAR.control     = list("gamma"=3, "eta"=1.4,"nu"=0.1),
    #delta.RAR       = 0,
    
    prob0           = c("UC"=1,"Simvastatin"=1),#,"Baricitinib"=1),
    N               = 529*2, # Assume the maximum cap of hypoinflammatory is reached
    interim         = list(recruited=list(m0 = 89*2  #89*3 # Trigger interim at 89 patients per arm
                                         ,m  = 49*2  # As per the recruitment expected Do interims at 49/ arm
                                         )),
    eff.arm         = efficacy.arm.fun, # Efficiency function of posteriors
    delta.eff       = log(1.1), # Select which interims select efficiency beta P(beta > delta.fut)
    eff.arm.control = list(b.eff = 0.84), # select the probability of the posterior > beta  
    fut.arm         = futility.arm.fun,
    delta.fut       = log(1.075), # select the analysed efficiency beta P(beta > delta.fut)
    fut.arm.control = list(b.fut = 1-0.78), # select the probability of the posterior > beta  
    delta.RAR       = 0,
    computation     = "sequential",
    #mc.cores        = parallel::detectCores()-1,
    H0              = FALSE,
    eff.trial=efficacy.arm.fun,
    fut.trial=futility.arm.fun,
    RAR = NULL,
    extended = 2)


saveRDS(scenario2,paste0(getwd(),"/../Results/Results_out.rds") )

scenario2
summary(scenario2)
scenario2$H1

for (i in 1:40){
    print(scenario2$H1$trial[[i]]) }

matrix(((scenario2$H1$estimate)), 
       ncol = 3, 
       byrow = TRUE)

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

scenario1<-scenario2

print(scenario1)
summary(scenario1)
plot(scenario1)

scenario1 |> summary()   #pe - maginal using the log efficacy OR // #pf - marginal using the log futility OR
scenario1$H1$trial[[1]]

scenario1$H1$estimate
matrix(((scenario1$H1$estimate)), 
       ncol = 3, 
       byrow = TRUE) # THIRD COLUMN IS THE OR

means_dist<-exp(mean(matrix(unlist(unlist(scenario1$H1$estimate)), ncol = 3, byrow = TRUE)[,3]))
hist(exp(matrix(unlist(unlist(scenario1$H1$estimate)), ncol = 3, byrow = TRUE))[,3],breaks=15)

sd(exp(matrix(unlist(unlist(scenario1$H1$estimate)), ncol = 3, byrow = TRUE)[,3]))


INLA::inla(OSFD ~  treatment , 
     family='pom',
     data = data, 
     control.fixed = list(mean = list( treat = 0), prec = 0.1),
     control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, config = TRUE))



