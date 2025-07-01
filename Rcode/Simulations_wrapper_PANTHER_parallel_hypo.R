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


# File title: Simulations_wrapper_PANTHER
# 
# Input: N/A
# 
# Output: Output a single table with informaition on eff/fut/non-conv/neither results in many sims
# 
# Purpose: Model PANTHER using BATTS and a lot of simulations


# ------------------ Author ---------------- Date ---------------- Update:
#            | Manel Benlloch Guajardo | 01MAY2025 |   Initial version of the file.

########################################################################################################
########################################################################################################
########################################################################################################
rm(list=ls())

library(BATSS)
library(tidyr)
library(magrittr)
library(dplyr)
library(INLA)

# Load the distribution of probabilities - and choose hypo/hyper samples - 
# CHANGE PATHS TO APPROPIATE ONE IF getwd() not selecting the right location
primOutDist_panth<- read.csv(paste0(getwd(),"../excel_distributions/VFDdistributions_logodds.csv"),header=TRUE) 
source(paste0(getwd(),"../Rcode/batss_glm_breakdown.R"))


# beta_0_select<-primOutDist_panth %>% dplyr::select(p_hyper_c,
#                                                    p_hyper_minus30, p_hyper_minus20, p_hyper_minus10,
#                                                    p_hyper_05, p_hyper_10, p_hyper_20, p_hyper_30,
#                                                    p_hyper_40, p_hyper_50,p_hyper_60) %>%
#     mutate(p_hyper_0=p_hyper_c)

# below the hypo distributions
beta_0_select<-
    primOutDist_panth %>% dplyr::select(p_hypo_pooled,
                                        p_hypo_minus10, p_hypo_05, p_hypo_10,
                                        p_hypo_20, p_hypo_30, p_hypo_40, p_hypo_50) %>%
    mutate(p_hypo_0=p_hypo_pooled)

number_betas<-ncol(beta_0_select)-1

if (length(commandArgs(trailingOnly=TRUE))!=0){
    args=commandArgs(trailingOnly=TRUE)
    print(args) #print information of underlying run
    Trials = as.numeric(args[1])
    pres_wd = as.character(args[2])} else {
    setwd(paste0(rstudioapi::getSourceEditorContext()$path,"/../.."))
    pres_wd = getwd()}
if(length(commandArgs(trailingOnly=TRUE)) <= 2) {
    number_cores   = 1 # setting the m variable 
    unique_core_id = 1
    core_id        = 1
    Trials_alloc   = Trials
} else {
    
    unique_core_id=as.numeric(args[3])
    number_cores = as.numeric(args[4])
    
    # Transformations following
    if(as.numeric(args[3])%%number_betas==0){core_id=number_betas}else{core_id=as.numeric(args[3])%%number_betas}
    Trials_alloc = ceiling(Trials/(number_cores/number_betas))
    message(glue::glue('This is the number of HYPER trials calculated in this core: {Trials_alloc} and these are the number of cores {number_cores} 
                     - but the id of this one is {unique_core_id}'))
} 


mc_cores=18
print(glue::glue('trying with {mc_cores}'))


Wrapper<- function(
        number_node =1,
        beta_list,
        model,var,var.control=NULL,family="gaussian",link="identity",
        which,alternative = "greater",R=1e+4,N,interim,prob0,
        delta.eff=0,delta.fut=delta.eff, delta.RAR=0,
        eff.arm,eff.arm.control=NULL,
        eff.trial=NULL,eff.trial.control=NULL,
        fut.arm,fut.arm.control=NULL,
        fut.trial=NULL,fut.trial.control=NULL,
        RAR=NULL,RAR.control=NULL,
        H0=TRUE,
        seed_set=1,
        computation="parallel",
        mc.cores=getOption("mc.cores", 3L),
        map_probabilities =FALSE,
        #linux.os = NA,
        extended=0, ...){
    
    data <- list()
    summary<-as.data.frame(matrix(NA,length(beta_list)-1,4,
                                  dimnames=list(paste0(colnames(beta_list[,-1])),
                                                c('Efficacy','Futility','Unresolved','Non-convergence'))))
    
    
    beta_baseline<-beta_list[1]
    beta<-list(beta_baseline,# this is the control - always
               beta_list[number_node+1])  # this is the treatment - changes by iteration
    
    
    glm_pom_run <- batss.glm.pom(seed_set= number_node,
                                 model=model,var=var,var.control=var.control,family=family,link=link,
                                 beta=beta,which,alternative ,R=R,N=N,
                                 interim=interim,
                                 prob0=prob0,
                                 delta.eff=delta.eff,delta.fut=delta.fut, delta.RAR=delta.RAR,
                                 eff.arm=eff.arm,eff.arm.control=eff.arm.control,
                                 eff.trial=eff.trial,eff.trial.control=eff.trial.control,
                                 fut.arm=fut.arm,fut.arm.control=fut.arm.control,
                                 fut.trial=fut.trial,fut.trial.control=fut.trial.control,
                                 RAR=RAR,RAR.control=RAR.control,
                                 H0=H0,
                                 computation='parallel',
                                 mc.cores=mc.cores,
                                 map_probabilities=map_probabilities,
                                 extended=extended,...)
    
    
    data[[number_node]]<-matrix(glm_pom_run$H1$estimate, ncol = 3, byrow = TRUE,
                                  dimnames=list(paste0('Trial-',colnames(beta_list[number_node]),'-',1:R),
                                                c('Looks','Result','posterior Mode'))) 
    
    summary[number_node,1]<-glm_pom_run$H1$target$global$efficacy[2] # efficacy
    summary[number_node,2]<-glm_pom_run$H1$target$global$futility[2] # futility
    summary[number_node,3]<-1-(glm_pom_run$H1$target$global$futility[2]+glm_pom_run$H1$target$global$efficacy[2]) # unresolved
    summary[number_node,4]<-glm_pom_run$H1$target$global$nonconverg[2] # issues
    
    out <-list(data = data, beta_list=beta_list, summary=summary,all_details=glm_pom_run,
               data_inputs=list(seed_set= number_node,
                              model=model,var=var,var.control=var.control,family=family,link=link,
                              beta=beta,which,alternative ,R=R,N=N,
                              interim=interim,
                              prob0=prob0,
                              delta.eff=delta.eff,delta.fut=delta.fut, delta.RAR=delta.RAR,
                              eff.arm=eff.arm,eff.arm.control=eff.arm.control,
                              eff.trial=eff.trial,eff.trial.control=eff.trial.control,
                              fut.arm=fut.arm,fut.arm.control=fut.arm.control,
                              fut.trial=fut.trial,fut.trial.control=fut.trial.control,
                              RAR=RAR,RAR.control=RAR.control,
                              H0=H0,
                              computation=computation,
                              mc.cores=mc.cores,
                              map_probabilities=map_probabilities,
                              extended=extended))
}


t2=Sys.time()

Trials<-8

results_wrap<-Wrapper(   
    number_node =m,
    beta_list =beta_0_select,
    model           = y ~ treatment,
    var             = list(y = multinomial_random,
                           treatment = treatalloc.fun),
    var.control     = list(y = list(size= 1)),
    family          = "pom",
    link            = 'identity',
    #beta            = list(multinom_rand_dset1,# this is the control
    #                       multinom_rand_dset2 # this is the treatment 
    #), 
    which           = c(2),   # Select which groups are treatments
    R               = Trials,
    control.fixed = list(mean = list( treatment = 0), prec = 0.1),
    control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, config = TRUE),
    alternative     = c("greater"), # One or two sided hypothesis
    map_probabilities = TRUE, # This variable will apply new maps to the days of
    
    prob0           = c("UC"=1,"Simvastatin"=1),#,"Baricitinib"=1),
    N               = 504*2, # Assume the maximum cap of hypoinflammatory is reached
    interim         = list(recruited=list(m0 = 89*2    #89*3 # Trigger interim at 89 patients per arm
                                          ,m  = 49*2)), # As per the recruitment expected Do interims at 49/ arm
    eff.arm         = efficacy.arm.fun, # Efficiency function of posteriors
    delta.eff       = log(1.1), # Select which interims select efficiency beta P(beta > delta.fut)
    eff.arm.control = list(b.eff = 0.84), # select the probability of the posterior > beta  
    fut.arm         = futility.arm.fun,
    delta.fut       = log(1.075), # select the analysed efficiency beta P(beta > delta.fut)
    fut.arm.control = list(b.fut = 1-0.78), # select the probability of the posterior > beta  
    delta.RAR       = 0,
    computation     = "parallel",
    mc.cores        = 4,
    H0              = FALSE,
    eff.trial=efficacy.arm.fun,
    fut.trial=futility.arm.fun,
    RAR = NULL,
    extended = 2)


t3<-Sys.time()
print(t3-t2)

results_wrap$data_inputs

saveRDS(results_wrap,
        paste0(pres_wd,'/Results/simulation',m,'.rds'))

# View this in the local R data session
results_wrap_100

