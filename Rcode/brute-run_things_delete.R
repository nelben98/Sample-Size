# generate data

model           = y ~ treatment
var             = list(y = multinomial_random,
                       treatment = treatalloc.fun)
var.control     = list(y = list(size= 1))
family          = "pom"
link            = 'identity'
beta            = list(multinom_rand_dset1,# this is the treatment
                       multinom_rand_dset2  # this is the control
                       )
which           = c(2)
   # Select which groups are treatments
R               = Trials
control.fixed = list(mean = list( treatment = 0), prec = 0.1)
control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, config = TRUE)
alternative     = c("greater")
 # One or two sided hypothesis
map_probabilities = TRUE
# This variable will apply new maps to the days of
# day free support- estimation.

# RAR:
# RAR option not used as we have uniform sampling 
#   - continuously sample to Max sample size/ Futitility/ Efficacy
#RAR             = prob.trippa,
#RAR.control     = list("gamma"=3, "eta"=1.4,"nu"=0.1),
#delta.RAR       = 0,
prob0           = c("UC"=1,"Simvastatin"=1) #,"Baricitinib"=1),
N               = 529*2 # Assume the maximum cap of hypoinflammatory is reached
interim         = list(recruited=list(m0=89*2 #89*3 # Trigger interim at 89 patients per arm
                                     ,m = 49*2))  # As per the recruitment expected Do interims at 49/ arm
eff.arm         = efficacy.arm.fun # Efficiency function of posteriors
delta.eff       = log(1.1)  # Select which interims select efficiency beta P(beta > delta.fut)
eff.arm.control = list(b.eff = 0.84)  # select the probability of the posterior > beta  
fut.arm         = futility.arm.fun
delta.fut       = log(1.075)  # select the analysed efficiency beta P(beta > delta.fut)
fut.arm.control = list(b.fut = 1-0.78)  # select the probability of the posterior > beta  
delta.RAR=0
computation     = "sequential"
#mc.cores        = parallel::detectCores()-1,
H0              = FALSE
eff.trial=efficacy.arm.fun
fut.trial=futility.arm.fun
RAR = NULL
extended = 1
delta.RAR=0
N = 529*2


R=3
    
#
#
seed_setting=10
lapply(R,
       batss.trial.pom,
       data=data,
       seed_setting=seed_setting,
       model=model,
       link=link,
       family=family,
       beta=beta,
       RAR=RAR,RAR.control=RAR.control,twodelta=twodelta,
       delta.eff=delta.eff,delta.fut=delta.fut,
       eff.arm=eff.arm,eff.trial=eff.trial,delta.RAR=delta.RAR,
       eff.arm.control=eff.arm.control,eff.trial.control=eff.trial.control,
       fut.arm=fut.arm,fut.trial=fut.trial,
       fut.arm.control=fut.arm.control,fut.trial.control=fut.trial.control,
       id.target=id.target,n.target=n.target,
       id.look=id.look,n.look=n.look,prob0=prob0,
       id.group=id.group,n.group=n.group,groupvar=groupvar,
       var=var,var.control=var.control,id.var=id.var,n.var=n.var,
       #linux.os=linux.os,
       map_probabilities =map_probabilities,
       extended = extended)
#  Now try a dataset until it breaks
#
#
pb <- txtProgressBar(min = 0, max = 1000, style = 3) 

#data_batts<-data
for (i in 1:1000){
    setTxtProgressBar(pb, i-1)
    POMPOMPOM<-batss.trial.pom(int = i,
                    data = data_batts,
                    model = model,
                    link = 'identity',
                    family,
                    beta,
                    prob0,
             RAR=NULL,RAR.control=NULL,
             eff.arm,eff.trial,
             eff.arm.control,
             eff.trial.control = NULL,
             fut.arm,
             fut.trial,
             fut.arm.control,
             fut.trial.control = NULL,
             id.target,
             n.target,
             id.look,
             n.look,
             id.group,
             n.group,
             groupvar,
             twodelta,delta.eff,
             delta.fut,delta.RAR,
             var,
             var.control,
             id.var,
             n.var,
             map_probabilities,
             extended)
    
    estimate = batss.res.e(list(POMPOMPOM),id.target)
    if(estimate[,"type",,drop=FALSE]==0 & is.na(estimate[,"mid",,drop=FALSE])){
        break
        close(pb)}
}


tar.p    = batss.res.tp(estimate,id.target)

tar.g    = batss.res.tg(estimate,id.target)
eff.p    = batss.res.ep(estimate,id.target,n.look)
eff.g    = batss.res.eg(estimate,id.target,n.look)
fut.p    = batss.res.fp(estimate,id.target,n.look)
fut.g    = batss.res.fg(estimate,id.target,n.look)
#
##
###
####
#####
######
#######
#######
data<-data_batts

set.seed((n.look+1)*int)  

# generate data for initial panel
n = m = N = prob = ref = target = ref = active = mu = posterior = NULL
env = new.env()
assign("m",id.look[1,"m"], envir = env)
assign("n",id.look[1,"n"], envir = env)
assign("prob",prob0 , envir = env)
assign("var",var , envir = env)
assign("var.control",var.control , envir = env)
assign("N", id.look$n[n.look], envir = env)  
assign("ref",id.group$ref, envir = env)

#call functions from 'var'
pos.col <- 2                                                                             # initialize the column indicator (starting at two, column one is the response)
for(vw in 2:length(var)) {                                                               # cycle through all variables in 'var' - starting from 2 because the first variable here is the response
    tmp_nam <- names(var)[vw]                                                              # store current variable name
    args_ <- plyr::.(n=m,m=m,prob=prob)                                                    # set function arguments, these are preset in the trial function loop
    if (tmp_nam %in% names(var.control)) args_ <- c(args_, var.control[[tmp_nam]])         # add additional arguments if specified in 'var.control'
    tmp_var <- R.utils::doCall(var[[tmp_nam]], envir = env, args = args_)                  # call variable generating function
    if (!is.matrix(tmp_var)) {                                                             # check if the generated data is NOT a matrix 
        data[, pos.col] <- tmp_var                                                           # fill column in 'data' according to position indication
        pos.col <- pos.col+1                                                                 # increase position indicator
    } else {
        colnames(tmp_var) <- paste0(tmp_nam,1:dim(tmp_var)[2])                               # name columns of matrix 'name'1,'name'2, etc
        for (jj in 1:dim(tmp_var)[2]) {
            data[, pos.col] <- tmp_var[,jj]                                                    # cycle through columns and fill 'data' accordingly
            pos.col <- pos.col+1                                                               # keep track of position in 'data'
        }
    } 
}
if (family !='pom'){
    #X = model.matrix(as.formula(paste0("~",strsplit(model,"~")[[1]][2])),data=data)
    X <- model.matrix(model[-2], data = data)                                            
    #---
    XB = X%*%beta
    assign("mu",switch(link,
                       "identity" = XB,
                       "log"      = exp(XB),
                       "logit"    = INLA::inla.link.logit(XB, inverse=TRUE),
                       "probit"   = INLA::inla.link.probit(XB, inverse=TRUE),
                       "robit"    = INLA::inla.link.robit(XB, inverse=TRUE),
                       "cauchit"  = INLA::inla.link.cauchit(XB, inverse=TRUE),
                       "loglog"   = INLA::inla.link.loglog(XB, inverse=TRUE),
                       "cloglog"  = INLA::inla.link.cloglog(XB, inverse=TRUE)),envir=env)
    
    tmp_nam <- names(var)[1] 
    args_ <- plyr::.(n=m,mu=mu)                                                             # create a quoted(!) list of available 'ingredients' 
    if (!(identical(var[[1]],rbinom) || identical(var[[1]],rnbinom))) {
        names(args_)[1:2] <- formalArgs(var[[1]])[1:2]  
    } else {
        if (identical(var[[1]],rbinom)) {
            names(args_)[1:2] <- c("n","prob")
        }
    }                                                                                       # rename the list objects to the names required by specified formula (NOTE: the order of the items is set, if a function requires a different order this will not work, clever rearranging may be needed)
    if (tmp_nam %in% names(var.control)) args_ <- c(args_, var.control[[tmp_nam]])          # add extra arguments if provided
    data[, id.var[1]] = R.utils::doCall(var[[1]], args = args_, envir = env)                # execute in 'env' environment with unused arguments allowed
    
} else if (family =='pom'){
    # n doesn't make sense to be used - as rmultinom is high dimensional so will not work as rbinom
    # and therefore feeding the number of rows substitutes the n=m
    X <-  model.matrix(model[-2], data = data) #data.matrix(qdapTools::mtabulate(as.data.frame(t(data))))
    unlisted_beta <- matrix(unlist(beta), ncol = 2, byrow = TRUE)
    XB = X%*%t(unlisted_beta)
    assign("mu",switch(link,
                       "identity" = XB),envir=env) # only one type of link function
    
    tmp_nam <- names(var)[1] 
    args_ <- plyr::.(#n=m, Not select n - as it will generate a SINGLE RMULTINOM per subject - so n=1 (numb subject implicit in number of rows)
        mu=mu)                                                             # create a quoted(!) list of available 'ingredients' 
    if (!(identical(var[[1]],rbinom) || identical(var[[1]],rnbinom))) {
        names(args_)[1] <- formalArgs(var[[1]])[1]  
    } else {
        if (identical(var[[1]],rbinom)) {
            names(args_)[1] <- c("prob")
        }
    }                                                                                       # rename the list objects to the names required by specified formula (NOTE: the order of the items is set, if a function requires a different order this will not work, clever rearranging may be needed)
    if (tmp_nam %in% names(var.control)) args_ <- c(args_, var.control[[tmp_nam]])          # add extra arguments if provided
    
    data[, id.var[1]] = R.utils::doCall(var[[1]], args = args_, envir = env) -2 # minus two to make the mapping {-1,28}
    
    if(map_probabilities){
        data[, id.var[1]]<- 
            dplyr::case_when(data[, id.var[1]] == -1  ~ 1,
                             data[, id.var[1]] >=  0  & data[, id.var[1]] <= 2  ~ 2,
                             data[, id.var[1]] >=  3  & data[, id.var[1]] <= 9 ~ 3,
                             data[, id.var[1]] >=  10 & data[, id.var[1]] <= 13 ~ 4,
                             data[, id.var[1]] >=  14 & data[, id.var[1]] <= 17 ~ 5,
                             data[, id.var[1]] >=  18 & data[, id.var[1]] <= 19 ~ 6,
                             data[, id.var[1]] >=  20 & data[, id.var[1]] <= 21 ~ 7,
                             data[, id.var[1]] >=  22 & data[, id.var[1]] <= 23 ~ 8,
                             data[, id.var[1]] >=  24 & data[, id.var[1]] <= 26 ~ 9,
                             data[, id.var[1]] >=  27 ~ 10)
        # if data is in factor - need to remove empty levels
        if(class(data[, id.var[1]]) =='factor'){
            data[, id.var[1]]<-as.numeric(droplevels(as.factor(data[, id.var[1]])))
        }
    }
}


inla(formula=model, 
     family='pom',
     data = data, 
     control.fixed = list(mean = list( treat = 0), prec = 0.1),
     verbose=TRUE,
     control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, config = TRUE))

#
##
###
####
#####
######
#######
#######

#### Example with code from 'Ed_Code_Bayesian_seqdes.R'
# master[[g]][[p]][[i]]
# Phenotype 1 
# profile 2 (0.1 OR)
# i (interim) -1

head(master[[1]][[2]][[1]])

inla(OSFD2 ~  treat , 
     family='pom',
     data = master[[1]][[2]][[1]], 
     control.fixed = list(mean = list( treat = 0), prec = 0.1),
     control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, config = TRUE))


model = y ~ treatment

X <- data.matrix(qdapTools::mtabulate(as.data.frame(t(data))))
model.matrix(model[-2], data = (qdapTools::mtabulate(as.data.frame(t(data)))))      


######################################    INVESTIGATE WHICH ONE IS BETTER WITH CODE
#
#
#


fit_hypo <-fit # fit_hypo statistic as per Ed's code

fit_hypo$marginals.fixed$treat	
summary(fit_hypo)$fixed[2,1] # mean
summary(fit_hypo)$fixed[2,2] # Sd
inla.pmarginal(log(efficacyOR),fit_hypo$marginals.fixed$treat)>=efficacyPThresh
inla.pmarginal(log(futilityOR),fit_hypo$marginals.fixed$treat)>=futilityPThresh


plot(fit_hypo$marginals.fixed$treat[,1]  ,fit_hypo$marginals.fixed$treat[,2])+ abline(v=log(1.1),col='green') + abline(v=log(1.075),col='red')

plot(fit_hypo$marginals.fixed$`(Intercept`[,1]  ,fit_hypo$marginals.fixed$`(Intercept`[,2])




setwd(paste0(rstudioapi::getSourceEditorContext()$path,"/.."))
#  Distributions in columns, each row corresponding to an outcome value.
primOutDist_panth<- read.csv(paste0(getwd(),"/../","excel_distributions/VFDdistributions_logodds.csv"),header=TRUE) 


ggplot(data=primOutDist_panth, aes(x=days)) +
    geom_line(aes(y=p_hypo_c,col='red'))+
    geom_line(aes(y=p_hypo_80,col='blue'))+
    #geom_line(aes(y=p_hypo_05))+
    theme_minimal() +
    theme(axis.text.x=element_text(size=rel(1.8)),
          axis.text.y=element_text(size=rel(1.8)),
          plot.title = element_text(size = rel(1.5)))+
    scale_x_continuous(name="Days", limits=c(-1, 28), n.breaks = 29) +
    scale_y_continuous(name="Probability")+
    ggtitle("Probability distribution organ support")






scenario2 = batss.glm.pom(   
    model           = y ~ treatment,  # POM regression specification
    var             = list(y = multinomial_random,
                           treatment = treatalloc.fun),
    var.control     = list(y = list(size= 1)), # Characteristics of Multinomial dist
    family          = "pom",
    link            = 'identity',
    beta            = list(multinom_rand_dset1,  # this is the dist of control
                           multinom_rand_dset2 ),# this is the dist of treatment 
    which           = c(2),   # Select which groups are treatments
    R               = Trials, # Number of trials to Simulate

    alternative     = c("greater"), # One or two sided hypothesis
    prob0           = c("UC"=1,"Simvastatin"=1),  ## Keep at 2 treatments for now
    
    N               = 504*2, # Assume the maximum cap of hypoinflammatory is reached
    interim         = list(recruited=list(m0 = 89*2   # Trigger interim at 89 patients per arm
                                          ,m  = 49*2  # As per the recruitment expected Do interims at 49/ arm
    )),
    
    eff.arm         = efficacy.arm.fun, # Efficiency function of posteriors
    delta.eff       = log(1.1), # Select which interims select efficiency beta P(beta > delta.fut)
    eff.arm.control = list(b.eff = 0.84), # select the probability of the posterior > beta  
    
    fut.arm         = futility.arm.fun,
    delta.fut       = log(1.075), # select the analysed efficiency beta P(beta > delta.fut)
    fut.arm.control = list(b.fut = 1-0.78), # select the probability of the posterior > beta  
    
    delta.RAR       = 0,
    eff.trial=efficacy.arm.fun,
    fut.trial=futility.arm.fun,
    
    # Information on Prior
    control.fixed = list(mean = list( treatment = 0), prec = 0.1),
    control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, config = TRUE), 
    
    # Extra specifications:
    computation     = "parallel",
    mc.cores        = parallel::detectCores()-1,
    H0              = FALSE,
    RAR = NULL,
    map_probabilities = TRUE, # This variable will apply new maps to the days of
                              #        day free support- estimation.
    extended = 2)

    