

# Load the following functions to input the program:
# 


#######################################################
#######################################################
#######################################################
#######################################################

# Multinomial function - bespoke to get the value on 30 days
# prob_dist - is set by the beta - so need to fit a list
# Size - Number of subjects to be randomize (since it's one observation per subject, and n-subjects is captured in the size of prob_dist) set to =1
    # Information needed for obtaining random samples- 
    # samp_n - is set by the m (ie how many individuals) - [obtained by rows in prob_dist]

multinomial_random <- 
    function(prob_dist,
             size){
        
        
        sapp_prob<-apply(prob_dist  # mu
                         ,1
                         ,function(x)
                             t(rmultinom(n=1 ,   # m information is in X
                                         size=size,
                                         prob=x)) 
        )|> t()
        
        return(c(matrix(apply(sapp_prob, 1,function(x)which( x==1)))))
    }

# function for the allocation (this is the randomisation which should be minimisation but fairly equal to this)
treatalloc.fun  = function(m,prob){
    prob = abs(prob)/sum(abs(prob)) 
    m0.g = floor(prob*m)
    m0   = sum(m0.g)
    factor(rep(names(prob),m0.g+rmultinom(1,m-m0,prob)),
           levels=names(prob))
}

        # test on m = 60 patients and equal allocation per group
#table(treatalloc.fun(m=60,prob=c(UC=1,Simvastatin=1,Baricitinib=1)))
#table(treatalloc.fun(m=61,prob=c(UC=1,Simvastatin=1,Baricitinib=1))) # test on 61, where last patient at random


# function For the efficacy check 
efficacy.arm.fun = function(posterior,b.eff){
    posterior > b.eff 
}

# function For the futility check 
futility.arm.fun = function(posterior,b.fut){
    posterior < b.fut
}

# test Efficacy and Test futility
# efficacy.arm.fun(0.6, b.eff = 0.84) # set > 0.84
# futility.arm.fun(0.9, b.fut=0.22)  # set < 0.22

# function
futility.trial.fun = function(fut.target){
    all(fut.target)
}

# test 
# futility.trial.fun(c(B=TRUE,C=TRUE,D=TRUE,E=TRUE,F=TRUE)) # Is any arm futilitydeclared?
# futility.trial.fun(c(B=TRUE,C=TRUE,D=TRUE,E=TRUE,F=FALSE))

#######################################################
#######################################################
#######################################################
#######################################################

# BATTS GLM FUNCTION - Modified for the POM analysis:
# Input - as the usual batts.glm, with the addition of map_probabilities 
#   map_probabilities - true/false: defines if the probabilities are mapped to less categories (algo converges)
# 
# Output - the same as original.

batss.glm.pom = function(
        model,var,var.control=NULL,family="gaussian",link="identity",
        beta,which,alternative = "greater",R=1e+4,N,interim,prob0,
        delta.eff=0,delta.fut=delta.eff, delta.RAR=0,
        eff.arm,eff.arm.control=NULL,
        eff.trial=NULL,eff.trial.control=NULL,
        fut.arm,fut.arm.control=NULL,
        fut.trial=NULL,fut.trial.control=NULL,
        RAR=NULL,RAR.control=NULL,
        H0=TRUE,computation="parallel",
        mc.cores=getOption("mc.cores", 3L),
        map_probabilities =FALSE,
        #linux.os = NA,
        extended=0, ...){
    
    #---    
    call <- match.call()                       # save call
    model <- as.formula(model)                 # allow for string and formula input
    #--- 
    
    ##
    ## dataset structure and useful definitions 
    ##
    
    message("    Initialisation")    
    
    #some checks
    n = m = prob = NULL
    #error messages
    if (!is.character(model) && ! plyr::is.formula(model)) 
        stop("invalid 'model' argument")
    if (!is.list(var) || !all(sapply(var,is.function)))
        stop("'var' must be a list of functions")
    if (is.null(intersect(names(var),intersect(setdiff(unlist(strsplit(all.vars(model),"1")),all.vars(model)),
                                               setdiff(unlist(strsplit(all.vars(model),"2")),all.vars(model))))) && 
        !setequal(all.vars(model),names(var)))
        stop("all variables in the model formula must have a generating function in the 'var' list")
    if (!(family %in% names(INLA::inla.models()$likelihood))){
        stop("invalid 'family' argument, see help files and inla documentation for available families")
    } else {
        if (!(family %in% c("gaussian","binomial","nbinomial","poisson")))
            warning("functionality only tested for gaussian, binomial, negative binomial and poisson distributions")
    }
    if (!(link %in% INLA::inla.models()$likelihood[[family]]$link) || !(link %in% c("identity","log","logit","probit","robit","cauchit","loglog","cloglog")))
        stop("'link' not supported, see help files and inla documentation for available link functions")
    if (!is.null(interim)){
        if(!inherits(interim,"list")){stop("'interim' should be a list")}
        interim.recruited = interim$recruited
    }else{
        stop("'interim' should be provided")
    } 
    if (!is.null(interim.recruited) && !is.numeric(unlist(interim.recruited))) 
        stop("'interim.recruited' must be a (list of) numeric vector(s)")
    if (!is.null(interim.recruited) && any(interim.recruited < 0)) 
        stop("negative interim recruitment numbers not allowed")
    if ((N < 0) || length(N)>1) {
        stop("total sample size 'N' must be a positive scalar")
        N <-  floor(N)
    }
    if (length(which)>length(beta) )
        stop("number of targets greater than number of parameters")
    if ((!is.null(RAR) && !(is.function(RAR))) || 
        (!is.null(eff.arm) && !(is.function(eff.arm))) || (!is.null(eff.trial) && !(is.function(eff.trial))) ||
        (!is.null(fut.arm) && !(is.function(fut.arm))) || (!is.null(fut.trial) && !(is.function(fut.trial))))
        stop("'RAR', 'eff.arm', 'eff.trial', 'fut.arm', 'fut.trial' must be functions or NULL")
    
    #warnings
    if (!is.null(interim.recruited) && any(interim.recruited > N)) {
        warning("some interim analyses are outside the maximum sample size and will be ignored")
        interim.recruited <- interim.recruited[interim.recruited<N]
    }
    if (all(prob0<0) && sum(prob0)!=1)
        warning("sum of 'prob0' not equal to 1")
    
    
    
    # size per look
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
    # group
    groupvar <- names(var)[2]
    n.group  = nlevels(data[,groupvar])
    id.group = data.frame(pos = 1:n.group,
                          id = levels(data[,groupvar]),
                          reference = levels(data[,groupvar])==levels(data[,groupvar])[1],
                          active = TRUE,
                          row.names = levels(data[,groupvar]))
    # define covariate name corresponding to group
    tmp = labels(terms(model))
    whichw = rep(FALSE,length(tmp))
    for(pw in 1:length(tmp)){
        if(is.factor(data[,tmp[pw]])){
            whichw[pw] = all(!is.na(match(id.group$id,levels(data[,tmp[pw]]))))&
                all(!is.na(match(levels(data[,tmp[pw]]),id.group$id)))
        }
    }
    if(sum(whichw)==1){
        groupvar = tmp[whichw]
    }else{
        if(sum(whichw)==0){
            stop("the variable corresponding to the treatment isn't the expected factor")
        }else{
            stop("2 factors share the same levels")
        }
    }
    # look
    id.look = cbind(id.look, 
                    matrix(NA,n.look,n.group,
                           dimnames=list(id.look$id,id.group$id)))
    # generate X matrix for names
    X = model.matrix(model[-2], data = data)                               #create model matrix straight from formula
    
    if(ncol(X)!=length(beta)){stop("length of 'beta' not compatible with X matrix")}
    names(beta) = colnames(X)
    # targets (only defined once)
    n.target  = length(which)
    if(length(alternative)==1){alternative=rep(alternative,n.target)
    }else{if(length(alternative)!=n.target){stop("length(alternative)!=n.target")}}
    id.target = data.frame(pos = NA,
                           id  = colnames(X)[which],
                           alternative = alternative,
                           group = NA, active = TRUE, 
                           look = NA, efficacy = NA, futility = NA, nonconver=NA, # added
                           low = NA, mid = NA, high = NA,
                           row.names =  colnames(X)[which])
    id.target$group = sapply(id.target$id,function(x){
        levels(data[,groupvar])[which(sapply(split(X[,x]!=0,data[,groupvar]),any))]
    })
    id.target = id.target[order(id.target$group),]
    id.target$pos = 1:n.target
    
    # delta vector(s)
    # mw = c(match("delta",names(eff.arm.control)),
    #        match("delta",names(fut.arm.control)))
    
    # none
    if(identical(delta.fut,delta.eff)){
        twodelta = FALSE 
        if (length(delta.eff)==1) delta.eff = delta.fut = rep(delta.eff,n.look)
        if (length(delta.eff)!=n.look) stop("length of delta not equal to number of looks")
    }else{
        # both
        if(!(is.null(eff.arm) || is.null(fut.arm))){  
            twodelta = TRUE        
            # efficacy
            if(length(delta.eff)==1){
                delta.eff = rep(delta.eff,n.look)
            }else{if(length(delta.eff)!=n.look){
                stop("length of delta not equal to number of looks")
            }}
            # futility
            if(length(delta.fut)==1){
                delta.fut = rep(delta.fut,n.look)
            }else{if(length(delta.fut)!=n.look){
                stop("length of delta not equal to number of looks")
            }}
            # one
        }else{
            twodelta = FALSE
            # unique
            #tmp = par[[mw[!is.na(mw)]]]
            tmp <- if (is.null(fut.arm)) delta.eff else delta.fut
            if(length(tmp)==1){
                tmp = rep(tmp,n.look)
            }else{if(length(tmp)!=n.look){
                stop("length of delta not equal to number of looks")
            }}
            # assign
            delta.eff = delta.fut = tmp
        }}  
    
    if (length(delta.RAR==1)) delta.RAR = rep(delta.RAR,n.look)
    if (length(delta.RAR)!=n.look) stop("length of delta.RAR not equal to number of looks")
    
    # trial stopping rules
    if(is.null(eff.trial) && !is.null(eff.arm)){
        eff.trial = function(eff.target){all(eff.target)}             #use function directly
        #---
    }
    if(is.null(fut.trial) && !is.null(fut.arm)){
        fut.trial = function(fut.target){all(fut.target)}             #use function directly
        #---
    }
    
    # seeds
    if(length(R)==1){id.seed=1:R}else{id.seed=R}    
    
    #############################
    # H1
    #############################
    
    if(if(is.list(beta)) {!all(beta[which] |> unlist() ==0)}
       else {!all(beta[which]==0)}){
        message("Evaluation of H1")           
        H1 = TRUE
        # sequential
        # 
        # 
        if(computation!="parallel"){
            trial_r = lapply(id.seed,
                             batss.trial.pom,
                             data=data,
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
                             extended = extended, ...)
            
            
            # parallel
        }else{if(computation=="parallel"){           
            # unix via forking
            if(Sys.info()[[1]]!="Windows"){
                trial_r = parallel::mclapply(id.seed,batss.trial.pom,
                                             data=data,model=model,link=link,family=family,beta=beta,
                                             RAR=RAR,RAR.control=RAR.control,twodelta=twodelta,delta.eff=delta.eff,delta.fut=delta.fut,
                                             eff.arm=eff.arm,eff.trial=eff.trial,delta.RAR=delta.RAR,
                                             eff.arm.control=eff.arm.control,eff.trial.control=eff.trial.control,
                                             fut.arm=fut.arm,fut.trial=fut.trial,
                                             fut.arm.control=fut.arm.control,fut.trial.control=fut.trial.control,
                                             id.target=id.target,n.target=n.target,
                                             id.look=id.look,n.look=n.look,prob0=prob0,
                                             id.group=id.group,n.group=n.group,groupvar=groupvar,
                                             var=var,var.control=var.control,id.var=id.var,n.var=n.var,
                                             # linux.os=linux.os,
                                             extended = extended, mc.cores=mc.cores,
                                             map_probabilities =map_probabilities,mc.set.seed = FALSE,...)
                # windows without forking
            }else{
                cl = parallel::makeCluster(getOption("cl.cores", mc.cores))
                parallel::clusterEvalQ(cl, c(library(INLA)))
                #parallel::clusterExport(cl, transfer, envir = .GlobalEnv)           
                #        parallel::clusterExport(cl, c(".expit"), envir = environment())           
                trial_r = parallel::parLapply(cl=cl,id.seed,batss.trial.pom,
                                              data=data,model=model,link=link,family=family,beta=beta,
                                              RAR=RAR,RAR.control=RAR.control,twodelta=twodelta,delta.eff=delta.eff,delta.fut=delta.fut,
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
                                              extended = extended,...)
                parallel::stopCluster(cl)    
            }
        }}     
        
        
        ##
        ## results
        ##
        estimate = batss.res.e(trial_r,id.target)
        tar.p    = batss.res.tp(estimate,id.target)
        tar.g    = batss.res.tg(estimate,id.target)
        eff.p    = batss.res.ep(estimate,id.target,n.look)
        eff.g    = batss.res.eg(estimate,id.target,n.look)
        fut.p    = batss.res.fp(estimate,id.target,n.look)
        fut.g    = batss.res.fg(estimate,id.target,n.look)
        sample   = batss.res.s1(trial_r,group=id.group$id,
                                type=c(apply(estimate[,"type",,drop=FALSE],2:3,paste0,collapse="")),
                                early=c(apply(estimate[,"look",,drop=FALSE]<n.look,2:3,all)))
        scenario = batss.res.s2(sample,target=id.target$id)
        res_H1   = list(estimate = estimate,
                        target   = list(par=tar.p,global=tar.g),
                        efficacy = list(par=eff.p,global=eff.g),
                        futility = list(par=fut.p,global=fut.g),
                        sample=sample,scenario=scenario)
        trial_H1    = trial_r
    }else{
        H1 = FALSE
    }
    #############################
    # H0
    #############################
    
    if( if(is.list(beta)) {all(beta[which] |> unlist() ==0)} else{all(beta[which]==0)} | H0==TRUE ){
        
        message("    Evaluation of H0")
        H0 = TRUE
        beta0 = beta
        beta0[which] = 0
        # sequential
        # 
        # 
        if(computation!="parallel"){
            trial_r = lapply(id.seed,batss.trial.pom,
                             data=data,model=model,link=link,family=family,beta=beta0,
                             RAR=RAR,RAR.control=RAR.control,twodelta=twodelta,delta.eff=delta.eff,delta.fut=delta.fut,
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
                             extended = extended, ...)
            
            # parallel    
        }else{if(computation=="parallel"){           
            # unix via forking
            if(Sys.info()[[1]]!="Windows"){
                trial_r = parallel::mclapply(id.seed,batss.trial.pom,
                                             data=data,model=model,link=link,family=family,beta=beta0,
                                             RAR=RAR,RAR.control=RAR.control,twodelta=twodelta,delta.eff=delta.eff,delta.fut=delta.fut,
                                             eff.arm=eff.arm,eff.trial=eff.trial,delta.RAR=delta.RAR,
                                             eff.arm.control=eff.arm.control,eff.trial.control=eff.trial.control,
                                             fut.arm=fut.arm,fut.trial=fut.trial,
                                             fut.arm.control=fut.arm.control,fut.trial.control=fut.trial.control,
                                             id.target=id.target,n.target=n.target,
                                             id.look=id.look,n.look=n.look,prob0=prob0,
                                             id.group=id.group,n.group=n.group,groupvar=groupvar,
                                             var=var,var.control=var.control,id.var=id.var,n.var=n.var,
                                             #linux.os=linux.os,
                                             extended = extended, 
                                             map_probabilities =map_probabilities,
                                             mc.cores=mc.cores,mc.set.seed = FALSE,...)
                # windows without forking
            }else{
                cl = parallel::makeCluster(getOption("cl.cores", mc.cores))
                parallel::clusterEvalQ(cl, c(library(INLA)))
                #parallel::clusterExport(cl, transfer, envir = .GlobalEnv)           
                #        parallel::clusterExport(cl, c(".expit"), envir = environment())           
                trial_r = parallel::parLapply(cl=cl,id.seed,batss.trial.pom,
                                              data=data,model=model,link=link,family=family,beta=beta0,
                                              RAR=RAR,RAR.control=RAR.control,twodelta=twodelta,delta.eff=delta.eff,delta.fut=delta.fut,
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
                                              extended=extended, ...)
                parallel::stopCluster(cl)    
            }
        }}     
        
        ##
        ## 
        ## 
        ##
        #browser()
        estimate = batss.res.e(trial_r,id.target)
        tar.p    = batss.res.tp(estimate,id.target)
        tar.g    = batss.res.tg(estimate,id.target)
        eff.p    = batss.res.ep(estimate,id.target,n.look)
        eff.g    = batss.res.eg(estimate,id.target,n.look)
        fut.p    = batss.res.fp(estimate,id.target,n.look)
        fut.g    = batss.res.fg(estimate,id.target,n.look)
        sample   = batss.res.s1(trial_r,group=id.group$id,
                                type=c(apply(estimate[,"type",,drop=FALSE],2:3,paste0,collapse="")),
                                early=c(apply(estimate[,"look",,drop=FALSE]<n.look,2:3,all)))
        scenario = batss.res.s2(sample,target=id.target$id)
        res_H0   = list(estimate = estimate,
                        target   = list(par=tar.p,global=tar.g),
                        efficacy = list(par=eff.p,global=eff.g),
                        futility = list(par=fut.p,global=fut.g),
                        sample=sample,scenario=scenario)
        trial_H0    = trial_r 
        
    }else{
        H0 = FALSE
    }
    
    ##
    ## output
    ##
    
    message("    Results")    
    look        = id.look[,c("pos","id","n","m")]
    FE  = data.frame(pos=1:ncol(X),id=colnames(X),target=FALSE,
                     row.names = colnames(X))
    FE[id.target$id,"target"] = TRUE
    if(H0  && family != 'pom'){FE[,'Beta (H0)'] = beta0}
    else if(H0 && family == 'pom') {FE[,'Beta (H0)'] = sapply(beta ,paste0)}
    if(H1 && family != 'pom'){FE[,'Beta (H1)'] = beta}
    else if(H1 && family == 'pom') {FE[,'Beta (H1)'] = sapply(beta ,paste0)}
    # 
    par = list(RAR=RAR, group=id.group[,c("pos","id","reference")],
               seed=id.seed, H0=H0, H1=H1, version=utils::packageVersion("BATSS"))
    out = list(beta = FE, look = look, par=par)    
    if(H0){
        out$H0 = res_H0
        if(extended>0){out$H0$trial = trial_H0}
    }
    if(H1){
        out$H1 = res_H1
        if(extended>0){out$H1$trial = trial_H1}
    }
    #---
    out$call <- call
    out$type <- "glm"
    #---
    class(out) = "batss"
    out
}






#######################################################
#######################################################
#######################################################
#######################################################

# BATSS trial -POM
# Same as the original batts.trial - with a few changes on:
#       1- How the random sample for subjects is generated if fam = 'POM'
#       2- adding map probabilities 
#       3- Adding some extra tests to smooth out some issues in the INLA convergence - and keep iterating


batss.trial.pom = function(int,data,model,link,family,beta,prob0,
                       RAR,RAR.control,
                       eff.arm,eff.trial,
                       eff.arm.control,eff.trial.control,
                       fut.arm,fut.trial,
                       fut.arm.control,fut.trial.control,
                       id.target,n.target,
                       id.look,n.look,
                       id.group,n.group,groupvar,
                       twodelta,delta.eff,delta.fut,delta.RAR,
                       var,var.control,id.var,n.var,
                       map_probabilities,
                       #linux.os=linux.os,
                       extended,...){
    # int=2
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
    #---
    #cat("B")
    # prepare
    
    posterior.fun = function(inf,fit,delta){
        prob = INLA::inla.pmarginal(delta, fit$marginals.fixed[[unlist(inf[1])]])
        ifelse(inf[2]=="greater",1-prob,prob)
    }    
    mx.posterior_eff.lt = mx.posterior_fut.lt = matrix(NA,nrow=n.look,ncol=n.target,dimnames=list(id.look$id,id.target$id))
    if (!is.null(RAR))  mx.posterior_RAR.lt = mx.posterior_eff.lt
    mx.futility.lt = mx.efficacy.lt = matrix(FALSE,nrow=n.look,ncol=n.target,
                                             dimnames=list(id.look$id,id.target$id))
    mx.rprob.lt = matrix(NA,nrow=n.look,ncol=n.group,
                         dimnames=list(id.look$id,id.group$id))   
    #cat("C")    
    
    dots <- rlang::dots_list(...,.named=TRUE) # dots=NULL
    
    # loop
    #if(INLA::inla.os.type()=="linux"&!is.na(linux.os)){
    #    INLA::inla.binary.install(os=linux.os,verbose=TRUE,md5.check=FALSE)       
    #    }
    
    for(lw in 1:n.look){# lw=0; lw=lw+1
        # size
        #cat(.p("look:",lw,"\n"))
        INLA_fail = FALSE # set a new variable to explore whether the INLA connection failed or not
        fit =NA
        
        temp = table(data[,groupvar])
        id.look[lw,names(temp)] = temp 
        assign("n",temp, envir = env)  
        assign("ref",id.group$ref, envir = env) 
        #cat("D")
        # fit 

        
        if ("control.family" %in% names(dots)) {
            control.link <- list(control.link=list(model=link))
            dots$control.family <- c(dots$control.family,control.link)
            
            fit = do.call(INLA::inla,c(list(formula=model, data=data, family=family,
                                            verbose=FALSE),dots))
        } else {
            # Added the try function - will attempt to resolve -if not just save as failed function
            fit = try(
                    do.call(INLA::inla,
                            c(list(formula=model, data=data, family=family,
                                   control.family=list(control.link=list(model=link)),
                                   verbose=FALSE)
                            ,dots)
                            )
                    ,silent = TRUE)
            if (grepl('Error',fit[1])) INLA_fail=TRUE
        }
        #cat("E")           
        # posteriors, efficacy and futility
        aw = id.target$active
        if (!is.null(eff.arm) & INLA_fail==FALSE) {
            mx.posterior_eff.lt[lw,aw] = apply(id.target[aw,c("id","alternative"),drop=FALSE],1,
                                               posterior.fun,fit=fit,delta=delta.eff[lw]) 
        } else if (!is.null(eff.arm) & INLA_fail==TRUE) {
            mx.posterior_eff.lt[lw,aw] = NA
        }  else if (is.null(eff.arm)) {
            mx.posterior_eff.lt[lw,aw] = NA
        }
        
        if ((twodelta || (is.null(eff.arm) && !is.null(fut.arm))) && INLA_fail==FALSE){ # NOT ADD ISSE - BUT DO IF ITS THE CASE
            mx.posterior_fut.lt[lw,aw] = apply(id.target[aw,c("id","alternative"),drop=FALSE],1,
                                               posterior.fun,fit=fit,delta=delta.fut[lw])               
        }else{
            if (!is.null(fut.arm) & INLA_fail==FALSE) {
                mx.posterior_fut.lt[lw,aw] = mx.posterior_eff.lt[lw,aw]   
            } else  if (!is.null(fut.arm) & INLA_fail==TRUE) {
                mx.posterior_fut.lt[lw,aw] = NA
            } else  if (is.null(fut.arm) ){
                mx.posterior_fut.lt[lw,aw] = NA
            }
        }
        if (!is.null(RAR)& INLA_fail==FALSE) {
            mx.posterior_RAR.lt[lw,aw] = apply(id.target[aw,c("id","alternative"),drop=FALSE],1,
                                               posterior.fun,fit=fit,delta=delta.RAR[lw])
        } else if (!is.null(RAR)& INLA_fail==TRUE) {
            mx.posterior_RAR.lt[lw,aw] = NA
        }
        
        #cat("F")           
        # update mx.futility.lt and mx.efficacy.lt
        for(tw in 1:n.target){
            if(aw[tw]){
                
                # efficacy
                assign("posterior",mx.posterior_eff.lt[lw,tw], envir = env)
                # assign("target",names(id.look[lw,names(temp)])==id.target[tw,"group"],envir = env)
                if (is.null(eff.arm) || is.na(delta.eff[lw])) {
                    mx.efficacy.lt[lw,tw] = FALSE
                } else {
                    mx.efficacy.lt[lw, tw] = R.utils::doCall(eff.arm, args = c(plyr::.(posterior=posterior,n=n,N=N,target=target,ref=ref),eff.arm.control), envir = env)        #call function instead of parsing and evaluating string
                }
                
                # futility
                if(twodelta || (is.null(eff.arm) && !is.null(fut.arm))){
                    assign("posterior",mx.posterior_fut.lt[lw,tw], envir = env)
                }
                if (is.null(fut.arm) || is.na(delta.fut[lw])) {
                    mx.futility.lt[lw,tw] = FALSE
                } else {
                    mx.futility.lt[lw, tw] = R.utils::doCall(fut.arm, args = c(plyr::.(posterior=posterior,n=n,N=N,target=target,ref=ref),fut.arm.control), envir = env)        #call function instead of parsing and evaluating string
                }
                #---
            }else{
                mx.efficacy.lt[lw,tw] = FALSE
                mx.futility.lt[lw,tw] = FALSE
            }
        }
        #cat("G")           
        eff.target = apply(mx.efficacy.lt[1:lw,,drop=FALSE],2,any)
        fut.target = apply(mx.futility.lt[1:lw,,drop=FALSE],2,any)
        
        #Futility and Efficacy - Exist / INLA has not failed to converge --> save as true/false to stop study
        if (!is.null(eff.arm) & INLA_fail==FALSE) {
            if(!is.na(eff.target)) eff.stop = eff.target else eff.stop =FALSE
        } else if (!is.null(eff.arm) & INLA_fail==TRUE) {eff.stop = FALSE 
        } else if (is.null(eff.arm)) {eff.stop = FALSE }
        
        if (INLA_fail==FALSE & !is.null(fut.arm)){ 
            if(!is.na(fut.target)) fut.stop =fut.target else fut.stop =FALSE
        } else if ( INLA_fail==TRUE & !is.null(fut.arm)){ fut.stop = FALSE
        } else if (is.null(fut.arm)) {fut.stop = FALSE}
        
        #---
        # efficacy   - Save results in vector to display    
        if(INLA_fail==FALSE & any(mx.efficacy.lt[lw,aw])){
            # identify arms
            ew = which(mx.efficacy.lt[lw,]&aw)
            # inactive arms according to eff.trial
            id.target$active[ew] = FALSE
            id.group[id.target$group[ew],"active"] = FALSE
            # save estimate and adapt list of target 
            id.target$look[ew]     = lw
            id.target$efficacy[ew] = TRUE                
            id.target[ew,c("low","mid","high")] = fit$summary.fixed[id.target$id[ew],
                                                                    c("0.025quant","mean","0.975quant")]
        }
        # futility   - Save results in vector to display
        if(INLA_fail==FALSE & any(mx.futility.lt[lw,aw])){
            # identify arms
            fw = which(mx.futility.lt[lw,]&aw)
            # inactive arms according to fut.trial
            id.target$active[fw] = FALSE
            id.group[id.target$group[fw],"active"] = FALSE
            # save estimate and adapt list of target 
            id.target$look[fw]     = lw
            id.target$futility[fw] = TRUE                
            id.target[fw,c("low","mid","high")] = fit$summary.fixed[id.target$id[fw],
                                                                    c("0.025quant","mean","0.975quant")]
        }        
        # stop trial due to no active parameters or last look
        all.stop = (eff.stop|fut.stop)| all(!id.target$active)| lw==n.look        

        if(all.stop){
            if(any(id.target$active)){
                aw = which(id.target$active)
                id.target$look[aw]     = lw
                if (INLA_fail==TRUE){id.target[aw,c("low","mid","high")] = c('NA','NA','NA')
                                    id.target$nonconver = TRUE # added
                                    message('Error in the whole running')} else{
                id.target[aw,c("low","mid","high")] = fit$summary.fixed[id.target$id[aw],
                                                                        c("0.025quant","mean","0.975quant")]}
            }
            break
            # continue
        }else{
            
            # prob per group
            if(!is.null(RAR)){
                # prob per group
                assign("posterior",mx.posterior_RAR.lt[lw,id.target$active], envir = env)
                assign("active",id.group$active, envir = env)
                #prob = .eval(RAR,envir=env) 
                assign("n",unlist(id.look[lw, id.group$id]),envir = env)  
                #assign ingredients to environment 'env' 
                assign("ref",id.group$ref,envir = env)
                assign("N",id.look$n[n.look],envir = env)
                assign("RAR.control", RAR.control, envir = env)
                prob = R.utils::doCall(RAR, args = c(plyr::.(posterior=posterior,n=n,N=N,ref=ref,active=active), RAR.control) ,envir = env)      #call function RAR in environment 'env'
                #---
            }else{
                prob = prob0[id.group$active]
            }
            names(prob) = id.group$id[id.group$active]
            mx.rprob.lt[lw,names(prob)] = prob/sum(prob)
            
            # predictors
            assign("n",id.look[lw+1,"n"],envir=env)
            assign("m",id.look[lw+1,"m"],envir=env)
            assign("prob",prob,envir=env)
            set.seed(lw+(n.look+1)*int)
            
            assign("var.control",var.control,envir=env)
            covar <- vector("list",length(var[-1]))
            for (ii in 1:length(var[-1])) {
                tmp_nam <- names(var)[ii+1]
                args_ <- plyr::.(n=m,m=m,prob=prob)
                if (tmp_nam %in% names(var.control)) args_ <- c(args_, var.control[[tmp_nam]])
                covar[[ii]] <- R.utils::doCall(var[[ii+1]], envir = env, args = args_)    
                if (is.matrix(covar[[ii]])) colnames(covar[[ii]]) <- paste0(tmp_nam,1:dim(covar[[ii]])[2])
            }
            
            #--- Now iterate on the next look - save the extra batch of patients as 'NEW' and add to the 
            #       previous batch of patients at the interim.
            new = as.data.frame(matrix(NA,id.look[lw+1,"m"],n.var,
                                       dimnames=list(paste0(lw+1,"-",1:id.look[lw+1,"m"]),id.var)))
            
            pos.col <- 2
            for (var.count in 1:(length(var)-1)){
                if (!is.matrix(covar[[var.count]])) {
                    new[,pos.col] = covar[[var.count]]
                    pos.col <- pos.col+1
                } else {
                    for (jj in 1:dim(covar[[var.count]])[2]) {
                        new[,pos.col] = covar[[var.count]][,jj]
                        pos.col <- pos.col+1
                    }
                }
            }
            
            if (family !='pom'){
                
                X <- model.matrix(model[-2], data = new)                                            
                XB = X%*%beta[colnames(X)]
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
                new[, id.var[1]] = R.utils::doCall(var[[1]], args = args_, envir = env)                # execute in 'env' environment with unused arguments allowed
                
                # data appending
                data = rbind(data,new)
            } else if (family =='pom'){
                # n doesn't make sense to be used - as rmultinom is high dimensional so will not work as rbinom
                # and therefore feeding the number of rows substitutes the n=m
                
                X <- data.matrix(qdapTools::mtabulate(as.data.frame(t(new)))) # if fails try using this  X<-model.matrix(model[-2], data = new)
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
                
                new[, id.var[1]] = R.utils::doCall(var[[1]], args = args_, envir = env) -2 # minus two to make the mapping {-1,28}
            
                if(map_probabilities){
                    new[, id.var[1]]<- 
                        dplyr::case_when(new[, id.var[1]] == -1  ~ 1,
                                         new[, id.var[1]] ==  0  ~ 2,
                                         new[, id.var[1]] >=  1  & new[, id.var[1]] <= 9 ~ 3,
                                         new[, id.var[1]] >=  10 & new[, id.var[1]] <= 13 ~ 4,
                                         new[, id.var[1]] >=  14 & new[, id.var[1]] <= 17 ~ 5,
                                         new[, id.var[1]] >=  18 & new[, id.var[1]] <= 19 ~ 6,
                                         new[, id.var[1]] >=  20 & new[, id.var[1]] <= 21 ~ 7,
                                         new[, id.var[1]] >=  22 & new[, id.var[1]] <= 23 ~ 8,
                                         new[, id.var[1]] >=  24 & new[, id.var[1]] <= 26 ~ 9,
                                         new[, id.var[1]] >=  27 ~ 10)
                }
                # data appending
                data = rbind(data,new)
            }
        }# end continue
        #cat(".")
    }# end loop
    
    # output
    # cat(paste0("\t end:",int,"\n"))
    
    colnames(mx.rprob.lt)      = paste0("r(",colnames(mx.rprob.lt),")")
    colnames(id.look)[-c(1:4)] = paste0("n(",colnames(id.look)[-c(1:4)],")")
    colnames(mx.posterior_eff.lt)  = paste0("pe(",colnames(mx.posterior_eff.lt),")")
    colnames(mx.posterior_fut.lt)  = paste0("pf(",colnames(mx.posterior_fut.lt),")")
    list(target = id.target, look = cbind(id.look,mx.posterior_eff.lt,mx.posterior_fut.lt,mx.rprob.lt),
         data   = if(extended==2){data}else{NULL})
}




############                                            ############
############                                            ############
############   -    OTHER NON EXPORTED FUNCTIONS        ############
############        FROM BATSS  - NO CHANGES ON THESE   ############
############                                            ############


# estimate matrix
batss.res.e = function(trial_r,id.target){       
    out = array(unlist(lapply(trial_r,function(x){
        type = rep(NA,nrow(x$target))
        for(i in 1:length(type)){
            eff = x$target$efficacy[i]&!is.na(x$target$efficacy[i])
            fut = x$target$futility[i]&!is.na(x$target$futility[i])
            type[i] = ifelse(eff&fut,3,ifelse(eff,1,ifelse(fut,2,0)))
        }
        c(x$target$look,type,x$target$mid)
    })),
    dim=c(nrow(id.target),3,length(trial_r)),
    dimnames=list(id.target$id,c("look","type","mid")))
    out[,"type",][is.na(out[,"type",])] = 0
    out
}
# target per parameter
batss.res.tp = function(estimate,id.target){
    out = id.target[,c("pos","id","alternative","group")] 
    out$efficacy   = apply(estimate[,"type",,drop=FALSE]==1,1,mean)
    out$futility   = apply(estimate[,"type",,drop=FALSE]==2,1,mean)
    out$both       = apply(estimate[,"type",,drop=FALSE]==3,1,mean) 
    out$nonconverg = apply(estimate[,"type",,drop=FALSE]==0,1,mean) # added
    colnames(out)[1] = ""    
    out
}

# target global   
batss.res.tg = function(estimate,id.target){
    out = data.frame(pos=1:2,id=c("At least one","All"),
                     alternative="",group="",efficacy=NA,futility=NA,both=NA,nonconverg=NA )# added
    out$nonconverg[1]= mean(apply(estimate[,"type",,drop=FALSE]==0,3,sum)>0)# added
    out$efficacy[1]  = mean(apply(estimate[,"type",,drop=FALSE]==1,3,sum)>0)
    out$futility[1]  = mean(apply(estimate[,"type",,drop=FALSE]==2,3,sum)>0)
    out$both[1]      = mean(apply(estimate[,"type",,drop=FALSE]==3,3,sum)>0)
    out$nonconverg[2]= mean(apply(estimate[,"type",,drop=FALSE]==0,3,all)>0)# added
    out$efficacy[2]  = mean(apply(estimate[,"type",,drop=FALSE]==1,3,all)>0)
    out$futility[2]  = mean(apply(estimate[,"type",,drop=FALSE]==2,3,all)>0)
    out$both[2]      = mean(apply(estimate[,"type",,drop=FALSE]==3,3,all)>0)
    colnames(out)[1] = ""    
    out
}  
# efficacy per target parameter
batss.res.ep = function(estimate,id.target,n.look){
    out = id.target[,c("pos","id","alternative","group")]
    out$early   = apply((estimate[,"type",,drop=FALSE]==1)*(estimate[,"look",,drop=FALSE]<n.look),1,mean)
    out$last    = apply((estimate[,"type",,drop=FALSE]==1)*(estimate[,"look",,drop=FALSE]==n.look),1,mean)
    out$overall = apply((estimate[,"type",,drop=FALSE]==1),1,mean)
    colnames(out)[1] = ""    
    out
}
# efficacy global
batss.res.eg = function(estimate,id.target,n.look){
    out = data.frame(pos=1:2,id=c("At least one","All"),
                     alternative="",group="",early=NA,last=NA,overall=NA)
    # at least one
    out$early[1] = mean(apply(estimate[,"type",,drop=FALSE]==1&(estimate[,"look",,drop=FALSE]<n.look),3,sum)>0)
    out$last[1]  = mean(apply(estimate[,"type",,drop=FALSE]==1&(estimate[,"look",,drop=FALSE]==n.look),3,sum)>0)
    out$overall[1]  = mean(apply(estimate[,"type",,drop=FALSE]==1,3,sum)>0)
    # all
    out$early[2] = mean(apply(estimate[,"type",,drop=FALSE]==1&(estimate[,"look",,drop=FALSE]<n.look),3,all)>0)
    out$last[2]  = mean(apply(estimate[,"type",,drop=FALSE]==1&(estimate[,"look",,drop=FALSE]==n.look),3,all)>0)
    out$overall[2]  = mean(apply(estimate[,"type",,drop=FALSE]==1,3,all)>0)
    # out
    colnames(out)[1] = ""    
    out
}  
# futility per target parameter
batss.res.fp = function(estimate,id.target,n.look){
    out = id.target[,c("pos","id","alternative","group")]
    out$early   = apply((estimate[,"type",,drop=FALSE]==2)*(estimate[,"look",,drop=FALSE]<n.look),1,mean)
    out$last    = apply((estimate[,"type",,drop=FALSE]==2)*(estimate[,"look",,drop=FALSE]==n.look),1,mean)
    out$overall = apply((estimate[,"type",,drop=FALSE]==2),1,mean)
    colnames(out)[1] = ""    
    out
}
# efficacy global
batss.res.fg = function(estimate,id.target,n.look){
    out = data.frame(pos=1:2,id=c("At least one","All"),
                     alternative="",group="",early=NA,last=NA,overall=NA)
    # at least one
    out$early[1] = mean(apply(estimate[,"type",,drop=FALSE]==2&(estimate[,"look",,drop=FALSE]<n.look),3,sum)>0)
    out$last[1]  = mean(apply(estimate[,"type",,drop=FALSE]==2&(estimate[,"look",,drop=FALSE]==n.look),3,sum)>0)
    out$overall[1]  = mean(apply(estimate[,"type",,drop=FALSE]==2,3,sum)>0)
    # all
    out$early[2] = mean(apply(estimate[,"type",,drop=FALSE]==2&(estimate[,"look",,drop=FALSE]<n.look),3,all)>0)
    out$last[2]  = mean(apply(estimate[,"type",,drop=FALSE]==2&(estimate[,"look",,drop=FALSE]==n.look),3,all)>0)
    out$overall[2]  = mean(apply(estimate[,"type",,drop=FALSE]==2,3,all)>0)
    # out
    colnames(out)[1] = ""    
    out
}  
batss.res.s1 = function(trial_r,group,type,early){
    size = as.data.frame(matrix(unlist(lapply(trial_r,function(x,group){
        x$look[max(x$target$look,na.rm=TRUE),paste0("n(",group,")")]
    },group=group)),byrow=TRUE,ncol=length(group),
    dimnames=list(names(trial_r),group)))
    cbind(size,type,early)    
}
batss.res.s2 = function(sample,target){
    tablew = table(sample$type)
    tablew = tablew[order(tablew,decreasing=TRUE)]
    out    = data.frame(pos=1:length(tablew),id=names(tablew),
                        overall=c(tablew)/sum(tablew),early=NA)
    early  = round(tapply(sample$early,sample$type,mean),2)
    out[names(early),"early"] = early
    for(i in 1:length(target)){
        out = cbind(out,as.numeric(substr(out$id,i,i)))
        colnames(out)[ncol(out)] = target[i]
    }
    out = cbind(out[,-c(3:4)],out[,3:4])    
    colnames(out)[1] = ""
    out
}



