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

Trials <- 10   #number of trials 

#  Distributions in columns, each row corresponding to an outcome value.
primOutDist_panth<- read.csv(paste0(getwd(),"/../","excel_distributions/VFDdistributions_logodds.csv"),header=TRUE) |> 
    #make distribution of active the same as the passive
    dplyr::select(p_hypo_c, p_hypo_t, p_hypo_pooled, p_hyper_c, p_hyper_t )

# samples_per_trial
# Generate multinomial outcomes - of the 5 distributions
sapp_prob<-sapply(primOutDist_panth, function(x)data.frame(rmultinom(Trials*3
                                     ,size=1
                                     ,prob=x)))

sapp_prob[,'p_hypo_c'] # cAN GET OUT A SPECIFIC DISTRIBUTION
sapp_prob[1,] # CAN GET OUT A SPECIFIC number WITHIN THE SAMPLING (all dists)


#######################################################################################################
################### Study simulation - set parameters                 ####################################
################### CONSIDER ONLY FOR HYPOIMMFLAMATORY:                 ####################################

 
treatment <- factor(rep(c("UC","Simvastatin","Baricitinib"),Trials*3)) # IN our model we still have 3 treatments
levels(treatment)

# Data from generated into rows - need to turn to columns
df_generated <-do.call(rbind.data.frame, sapp_prob[,1])
y<-apply(df_generated[,1:30], 1,function(x)which( x==1 ))-2

tmp <- file.path(tempdir(), "test.log")


# Multinomial function - bespoke to get the value on 30 days
# Information i need - 
# samp_n - is set by the m (ie how many individuals)
# prob_dist - is set by the beta - so need to fit a list
# as no prob_dist - also not num_dist

multinomial_generation <- 
    function(n = nrow(prob_dist), # default to be nrows of prob dist- if want diferent, specify
             prob_dist,
             size){
        
        sapp_prob<-apply(prob_dist  # mu
                         ,1
                         ,function(x)
                             t(rmultinom(n=n ,   # m information is in X
                                         size=size,
                                         prob=x)) 
        )|> t()
        
        return(c(matrix(apply(sapp_prob, 1,function(x)which( x==1)))))
    }

# different version - defined in battss_glm_breakdown - not using n (n is number of wors in prob_dist)
multinomial_generation(
                    prob_dist   = t(primOutDist_panth[,1]) #primOutDist_panth[1:4,]
                   ,size        = 1)

apply(t(primOutDist_panth[,1])  # mu
     ,1
     ,function(x)
         t(rmultinom(n=1 ,   # m information is in X
                     size=1,
                     prob=x)) 
        )|> t()

a <- multinomial_random(prob_dist=t(primOutDist_panth[,1]),
                        #num_dist =1,
                        size=1) # test how it would work

a


# function for the allocation (this is the randomisation which should be minimisation but fairly equal to this)
treatalloc.fun  = function(m,prob){
    prob = abs(prob)/sum(abs(prob)) 
    m0.g = floor(prob*m)
    m0   = sum(m0.g)
    factor(rep(names(prob),m0.g+rmultinom(1,m-m0,prob)),
           levels=names(prob))
}

# test on m = 60 patients and equal allocation per group
table(treatalloc.fun(m=60,prob=c(UC=1,Simvastatin=1,Baricitinib=1)))
table(treatalloc.fun(m=61,prob=c(UC=1,Simvastatin=1,Baricitinib=1))) # test on 61, where last patient at random

# function For the efficacy check 
efficacy.arm.fun = function(posterior,b.eff){
    posterior > b.eff 
}

# function For the futility check 
futility.arm.fun = function(posterior,b.fut){
    posterior < b.fut
}
# test Efficacy and Test futility
efficacy.arm.fun(0.6, b.eff = 0.84) # set > 0.84
futility.arm.fun(0.9, b.fut=0.22)  # set < 0.22

# function
futility.trial.fun = function(fut.target){
    all(fut.target)
}

# test 
futility.trial.fun(c(B=TRUE,C=TRUE,D=TRUE,E=TRUE,F=TRUE)) # Is any arm futilitydeclared?
futility.trial.fun(c(B=TRUE,C=TRUE,D=TRUE,E=TRUE,F=FALSE))

#######################################################################################################

# END of setup - Start of the coding of the trial:

# simulation
#BATSS using INLA's default normal priors, N(0,1000)
Trials<-25




## Now try the code from inla above but with the data to be input in batss.glm

Pheno<-1
multinom_rand_dset1<-data.frame(control_dist  = c(primOutDist_panth[1:3,Pheno]+0.03,
                                       1-sum(primOutDist_panth[1:3,Pheno]+0.03))) # to be deleted - add random measure so it's not 0 everywhere 
Pheno<-2
multinom_rand_dset2<-data.frame(treatment_dist= c(primOutDist_panth[1:3,Pheno]+0.03,
                                       1-sum(primOutDist_panth[1:3,Pheno]+0.03))) # to be deleted - add random measure so it's not 0 everywhere 


multinom_rand_dset1<- data.frame(primOutDist_panth[,1])
multinom_rand_dset2<- data.frame(primOutDist_panth[,2])

source(paste0(getwd(),"/batss_glm_breakdown.R"))
scenario1 = batss.glm.pom(   
    model           = y ~ treatment,
    var             = list(y = multinomial_random,
                           treatment = treatalloc.fun),
    var.control     = list(y = list(size= 1)),
    family          = "pom",
    link            = 'identity',
    beta            = list(multinom_rand_dset2,
                           multinom_rand_dset1  # this is the control
                           ), # this is the treatment
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
    N               = 504*3, # Assume the maximum cap of hypoinflammatory is reached
    interim         = list(recruited=list(m0=80*3 #89*3 # Trigger interim at 89 patients per arm
                                         ,m = 49*3  # As per the recruitment expected Do interims at 49/ arm
                                         )),

    eff.arm         = efficacy.arm.fun, # Efficiency function of posteriors
    delta.eff       = 1.1, # Select which interims select efficiency beta P(beta > delta.fut)
    eff.arm.control = list(b.eff = 0.78), # select the probability of the posterior > beta  
    fut.arm         = futility.arm.fun,
    delta.fut       = 1.075, # select the analysed efficiency beta P(beta > delta.fut)
    fut.arm.control = list(b.fut = 0.22), # select the probability of the posterior > beta  
    delta.RAR=0,
    computation     = "sequential",
    #mc.cores        = parallel::detectCores()-1,
    H0              = FALSE,
    eff.trial=NULL,
    fut.trial=NULL,
    RAR = NULL,
    extended = 1)

print(scenario1)
summary(scenario1)
plot(scenario1)

scenario1$H0




INLA::inla(OSFD ~  treatment , 
     family='pom',
     data = data, 
     control.fixed = list(mean = list( treat = 0), prec = 0.1),
     control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, config = TRUE))



