# generate data

model           = y ~ treatment
var             = list(y = multinomial_random,
                       treatment = treatalloc.fun)
var.control     = list(y = list(size= 1))
family          = "pom"
link            = 'identity'
beta            = list(multinom_rand_dset2,# this is the treatment
                       multinom_rand_dset1  # this is the control
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
N               = 504*2 # Assume the maximum cap of hypoinflammatory is reached
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
N = 504*2
int=1
    
if(sum(!is.na(match(c("m0","m"),names(interim.recruited))))==2){
    size_look = seq(interim.recruited$m0,N,interim.recruited$m)  
    size_look[length(size_look)] = N
}else{
    interim.recruited <- sort(interim.recruited)        #sort interim recruited to make sure the values are ordered
    size_look = c(interim.recruited,N)  
}
n.look  = length(size_look)   

id.look = data.frame(pos = 1:n.look,
                     id  = paste0("n=",size_look),
                     n   = size_look,
                     m   = c(size_look[1],size_look[-1]-size_look[-n.look]))
# generate predictors
env0 = new.env()
assign("m", id.look$m[1], envir = env0)
assign("n", id.look$n[1], envir = env0)
assign("prob",prob0,envir = env0)

assign("var.control", var.control, envir = env0)  
covar <- vector("list",length(var[-1]))
for (ii in 1:length(var[-1])) {
    tmp_nam <- names(var)[ii+1]
    args_ <- plyr::.(n=n,m=m,prob=prob)
    if (tmp_nam %in% names(var.control)) args_ <- c(args_, var.control[[tmp_nam]])
    covar[[ii]] <- R.utils::doCall(var[[ii+1]], envir = env0, args = args_)   #call functions directly
    if (is.matrix(covar[[ii]])) colnames(covar[[ii]]) <- paste0(tmp_nam,1:dim(covar[[ii]])[2])
}

n.var <- length(all.vars(model))
id.var <- names(var)
if (any(sapply(covar,is.matrix))) {
    where.mat <- which(sapply(covar,is.matrix))
    tmp_var <- covar[where.mat]
    id.var <- names(var)[1]
    for (ii in 1:length(covar)) {
        if (ii %in% (where.mat)) {
            id.var <- c(id.var,colnames(covar[[ii]]))
        } else {
            id.var <- c(id.var,names(var[ii+1]))
        }
    }
}

m0     = length(covar[[1]])
if(length(m0)>1|m0[1]!=id.look$m[1]){stop("different predictor length")}

data <- as.data.frame(matrix(NA,m0,n.var,
                             dimnames=list(paste0("1-",1:m0),id.var)))
pos.col <- 2
for (var.count in 1:(length(var)-1)){
    if (!is.matrix(covar[[var.count]])) {
        data[,pos.col] = covar[[var.count]]
        pos.col <- pos.col+1
    } else {
        for (jj in 1:dim(covar[[var.count]])[2]) {
            data[,pos.col] = covar[[var.count]][,jj]
            pos.col <- pos.col+1
        }
    }
}
#
#
#
#


#  Now try a dataset until it breaks
#
#

    int<-2
    # cat(paste0("\t start:",int,"\n"))
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
                           "log" = exp(XB),
                           "logit" = INLA::inla.link.logit(XB, inverse=TRUE),
                           "probit" = INLA::inla.link.probit(XB, inverse=TRUE),
                           "robit" = INLA::inla.link.robit(XB, inverse=TRUE),
                           "cauchit" = INLA::inla.link.cauchit(XB, inverse=TRUE),
                           "loglog" = INLA::inla.link.loglog(XB, inverse=TRUE),
                           "cloglog" = INLA::inla.link.cloglog(XB, inverse=TRUE)),envir=env)
        
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
                           "identity" = XB),envir=env)
        
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
                                 data[, id.var[1]] ==  0  ~ 2,
                                 data[, id.var[1]] >=  1  & data[, id.var[1]] <= 9 ~ 3,
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
    
    # open the dots! - INLA function without ... reading
    # 
    
    
    fit =NA
    fit =try( inla(formula=model, data=data, family=family,
                   control.family=list(control.link=list(model=link)),
                   control.fixed = list(mean = list( treat = 0), prec = 0.1),
                   control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, config = TRUE),
                   verbose=FALSE)
              ,silent = TRUE)
    if(!grepl('Error',fit[1])){print(i)} # if find error in the computation of INLA



fit<-NULL
3 - 5 - 11 # these are the ones which count





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


model           = y ~ treatment

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
