### Bayesian group sequential design simulation code
### Ordinal primary outcome, proportional odds model
### Author: Ed Waddingham, Imperial Clinical Trials Unit
### February 2025


setwd(paste0(rstudioapi::getSourceEditorContext()$path,"/.."))


## Load libraries
library(remotes)
library(MASS)
library(extraDistr)
library(dplyr)
library(brms)
library(INLA)
library(brinla)
library(ggplot2)
library(gridExtra)
library(stringr)

## Edit this line to read in the distrbutions for the primary outcome.
#  Distributions in columns, each row corresponding to an outcome value.
#  Code assumes names will be, eg, "pcontrol_subgrp1" for control in subphenotype 1
#  and  "pactive_OR1_1_subgrp2" for active treatment with OR=1.1 in subphenotype 2
primOutDist<-read.csv(paste0(getwd(),"/../","excel_distributions/prob_prim_out_VANISH.csv"),header=TRUE) %>%
        #make distribution of active the same as the passivel
        mutate(pactive_OR1_2_subgrp1 = pcontrol_subgrp1, 
               pactive_OR1_2_subgrp2 = pcontrol_subgrp2,
               pactive_OR1_4_subgrp1 = pcontrol_subgrp1, 
               pactive_OR1_4_subgrp2 = pcontrol_subgrp2)

#### KEY INPUT PARAMETERS TO SET START HERE

# number of simulations
nSims<-10

# set random number seed
set.seed(1)


interimN<-list()
## Edit the lines below to read in the sample size at each interim.
## The input files should have column 1 for control and column 2 for active treatment.
## Row 1 gives the sample size at interim analysis 1, etc.
## The final row should contain the final max sample size.

# subphenotype 1 - Hypo (70%)
interimN[[1]]<-read.csv(paste0(getwd(),"/../","excel_distributions/interimN1.csv"),header=TRUE)

# Subpheno 2 - Hyper (30%)
interimN[[2]]<-read.csv(paste0(getwd(),"/../","excel_distributions/interimN2.csv"),header=TRUE)

# total SS (both arms) in each subphenotype and overall
phenoN<-list()
phenoN[[1]]<-interimN[[1]][,1]+interimN[[1]][,2]
phenoN[[2]]<-interimN[[2]][,1]+interimN[[2]][,2]
totalN<-phenoN[[1]]+phenoN[[2]]

# number of interims
nInterims<-nrow(interimN[[1]])


# indicators showing which phenotypes are included in each interim
phenoInt<-list()
# (this is based on the PANTHER design where the smaller phenotype is not always included)

# phenotype 1
phenoInt[[1]]<-rep(1,nInterims)
# phenotype 2
phenoInt[[2]]<-c(0,rep(0:1, times=nInterims/2))

# number and names of phenotypes
nPhenotypes<-2
phenotypes<-c("1","2")


# active arm sample size cap in each phenotype (again, part of PANTHER design)
# Recruitment to an arm is stopped if the cap (based on equivalent frequentist SS) is reached
sampleCap<-c(504,529)  # (phenotype 1, phenotype 2)


## Stopping triggers
# odds ratio to be exceeded (and minimum posterior probability thereof) in order to demonstrate efficacy
efficacyOR<-1.1
efficacyPThresh<-0.84
# odds ratio to be undershot (and minimum posterior probability thereof) in order to demonstrate futility
futilityOR<-1.075
futilityPThresh<-0.78

## active treatment arm profiles to simulate for each treatment
## operating characteristics are calculated for each of these profiles
## and these are put together to construct power curves.
## Profile names should correspond to the names of distributions in primOutDist but without the "pactive_OR"nprefix and "_subgrpx" suffix
profiles<-c("1_2","1_4")
nProfiles<-length(profiles)

#### INPUT PARAMETERS TO SET END HERE


# initialise dataframes for simulated data
    # and results dataframes/vectors


simData<-list()
postCoefmean<-list()
postCoefse<-list()

stopEfficacy<-list()
stopFutility<-list()
trialResult<-list()
stoppingPoint<-list()
postOR<-list()
postOREff<-list()
postORFut<-list()
percentSuccess<-list()
percentFutile<-list()
percentExceededCap<-list()
perDrugPhenoN<-list()
power<-list()
alpha<-list()
percentDecided<-list()
percentUndecided<-list()
correctFut<-list()
incorrectFut<-list()
pheno_N<-list()
actualN<-list()
expectedN<-list()
q80N<-vector()
controlStoppingPoint<-list()
overallStoppingPoint<-rep(0,nSims)
nStoppedTreatments<-list()

# Setup of the dataframes/lists that will hold all the information during the simulations:
for (g in 1:nPhenotypes){
    simData[[g]]<-list()
    stopEfficacy[[g]]<-list()
    stopFutility[[g]]<-list()
    postCoefmean[[g]]<-list()
    postCoefse[[g]]<-list()
    trialResult[[g]]<-list()
    controlStoppingPoint[[g]]<-rep(0,nSims)
    stoppingPoint[[g]]<-list()
    pheno_N[[g]]<-list()
    postOR[[g]]<-list()
    postOREff[[g]]<-list()
    postORFut[[g]]<-list()
    percentSuccess[[g]]<-list()
    percentFutile[[g]]<-list()
    percentExceededCap[[g]]<-list()
    perDrugPhenoN[[g]]<-list()
    power[[g]]<-list()
    alpha[[g]]<-list()
    percentDecided[[g]]<-list()
    percentUndecided[[g]]<-list()
    correctFut[[g]]<-list()
    incorrectFut[[g]]<-list()
    expectedN[[g]]<-list()
    
    for (p in 1:nProfiles) {
        
        simData[[g]][[p]]<-data.frame(matrix(nrow=interimN[[g]][nInterims,1]+interimN[[g]][nInterims,2],ncol=2))
        names(simData[[g]][[p]])<-c("treat","OSFD")
        simData[[g]][[p]]$treat<-c(rep(0,interimN[[g]][nInterims,1]),rep(1,interimN[[g]][nInterims,2]))
        simData[[g]][[p]]$pheno<-phenotypes[g]
        postCoefmean[[g]][[p]]<-data.frame(matrix(nrow=nSims,ncol=nInterims))
        postCoefse[[g]][[p]]<-data.frame(matrix(nrow=nSims,ncol=nInterims))
        stopEfficacy[[g]][[p]]<-data.frame(matrix(nrow=nSims,ncol=nInterims))
        stopFutility[[g]][[p]]<-data.frame(matrix(nrow=nSims,ncol=nInterims))
        trialResult[[g]][[p]]<-data.frame(matrix(nrow=nSims,ncol=nInterims))
        perDrugPhenoN[[g]][[p]]<-list()
        stoppingPoint[[g]][[p]]<-rep(nInterims,nSims)
        trialResult[[g]][[p]][,]<-"Inconclusive"
        postOR[[g]][[p]]<-vector()
        postOREff[[g]][[p]]<-list()
        postORFut[[g]][[p]]<-list()
        expectedN[[g]][[p]]<-vector()
        percentSuccess[[g]][[p]]<-vector()
        percentFutile[[g]][[p]]<-vector()
        percentExceededCap[[g]][[p]]<-vector()
        percentDecided[[g]][[p]]<-vector()
        percentUndecided[[g]][[p]]<-vector()
        for (i in 1:nInterims){
            postOREff[[g]][[p]][[i]]<-vector()
            postORFut[[g]][[p]][[i]]<-vector()
        }
    }
    power[[g]]<-vector()
    alpha[[g]]<-vector()
    correctFut[[g]]<-vector()
    incorrectFut[[g]]<-vector()
}


# Simulate outcome data:
# Outcome is simulated for all participants up to max sample size
# The data is later subsetted to give a dataset for each INTERIM analysis

for (m in 1:nSims) { #Repeat the study for each simulation
    rawSimData_c<-list()
    rawSimData_t<-list()
    simDataAll<-list()
    simDataRand<-list()
    
    for (g in 1:nPhenotypes){
        
        # control arm: simulate draws from multinomial distribution - for each of their phenotypes
        rawSimData_c[[g]]<-
            t(rmultinom(interimN[[g]][nInterims,1]
                       ,size=1
                       ,prob=primOutDist[,names(primOutDist)==paste0("pcontrol_subgrp",phenotypes[g])]
                       )
              )
        
        # treatment arm: simulate draws from multinomial distribution for each treatment profile (two treatments)
        rawSimData_t[[g]]<-list()
        
        for (p in 1:nProfiles){ # FOR EACH OF TWO TREATMENTS
            rawSimData_t[[g]][[p]]<-
                t(rmultinom(interimN[[g]][nInterims,2]
                           ,size=1
                           ,prob=primOutDist[,names(primOutDist)==paste0("pactive_OR",profiles[p],"_subgrp",phenotypes[g])]
                           )
                  )
            
            # + The above yields data with a column for each possible outcome value and 1s indicating the outcome for each observation
            #    next few lines convert this to a single column giving the outcome value
            # + '2' is subtracted because the lowest category for our organ support outcome is not 1 but -1 (representing death) 
            #    here the control and active data are also stacked into a single dataset (one for each profile/phenotype)
            for (i in 1:nrow(simData[[g]][[p]])) {
                if(simData[[g]][[p]]$treat[i]==0){
                    simData[[g]][[p]]$OSFD[i]<-which(rawSimData_c[[g]][i,]==1)-2 #shift to account for death ('-1')
                }			
                if(simData[[g]][[p]]$treat[i]==1){
                    simData[[g]][[p]]$OSFD[i]<-which(rawSimData_t[[g]][[p]][i-interimN[[g]][nInterims,1],]==1)-2		
                }
            }
        }
        
    }# Now we have simulated the subject's survival days '-1' --> '28' days
    
    
    
    interimDataAll<-list()
    simDataAll<-list()
    simDataRand<-list()
    
    ## + Append the simulated data in both phenotypes together to create a 
    ##   single dataset simDataAll[[p]] for each profile p
    ## + Subset this to create interimDataAll[[p]][[i]] for each interim i
    
    for (p in 1:nProfiles){
        interimDataAll[[p]]<-list()
        
        simDataAll[[p]]<-data.frame(matrix(nrow=0,ncol=2)) # create empty dframe to hold both cols
        names(simDataAll[[p]])<-c("treat","OSFD")
        
        for (g in 1:nPhenotypes){
            simDataAll[[p]]<-rbind(simDataAll[[p]],simData[[g]][[p]])
        }
        # randomise order of simulated data
        simDataRand[[p]]<-simDataAll[[p]][order(sample(1:nrow(simDataAll[[p]]))),]
        rownames(simDataRand[[p]])<-NULL
        
        # subset first N_int patients for each interim
        for (i in 1:(nInterims)){
            
            interimDataAll[[p]][[i]]<-simDataRand[[p]][1:totalN[i],]
        }	
    }
    
    ### Set up analysis datasets 
    ###     -INLA matrix analysis
    ###     -Posterior dataset
    ###     -Master dataset by phenotype
    
    modelINLA<-list()
    master<-list()
    post<-list()
    
    for (g in 1:nPhenotypes){
        modelINLA[[g]]<-list()
        master[[g]]<-list()
        post[[g]]<-list()
        
        
        for (p in 1:nProfiles){
            modelINLA[[g]][[p]]<-list()
            post[[g]][[p]]<-list()
            master[[g]][[p]]<-list()
        }
        
    }
    
    
    for (p in 1:nProfiles){
        for (i in 1:(nInterims)){
            
            
            # Recoding the organ support free days into 10 categories 
            #  (Since INLA can only fit proportional odds model for up to 
            #  10 categories)
            interimDataAll[[p]][[i]]$OSFD2 <- case_when(interimDataAll[[p]][[i]]$OSFD == -1  ~ 1,
                                                        interimDataAll[[p]][[i]]$OSFD ==  0  ~ 2,
                                                        interimDataAll[[p]][[i]]$OSFD >=  1  & interimDataAll[[p]][[i]]$OSFD <= 9 ~ 3,
                                                        interimDataAll[[p]][[i]]$OSFD >=  10 & interimDataAll[[p]][[i]]$OSFD <= 13 ~ 4,
                                                        interimDataAll[[p]][[i]]$OSFD >=  14 & interimDataAll[[p]][[i]]$OSFD <= 17 ~ 5,
                                                        interimDataAll[[p]][[i]]$OSFD >=  18 & interimDataAll[[p]][[i]]$OSFD <= 19 ~ 6,
                                                        interimDataAll[[p]][[i]]$OSFD >=  20 & interimDataAll[[p]][[i]]$OSFD <= 21 ~ 7,
                                                        interimDataAll[[p]][[i]]$OSFD >=  22 & interimDataAll[[p]][[i]]$OSFD <= 23 ~ 8,
                                                        interimDataAll[[p]][[i]]$OSFD >=  24 & interimDataAll[[p]][[i]]$OSFD <= 26 ~ 9,
                                                        interimDataAll[[p]][[i]]$OSFD >=  27 ~ 10)
            
            # split the data by phenotype again
            for (g in 1:nPhenotypes){
                master[[g]][[p]][[i]]<-interimDataAll[[p]][[i]][interimDataAll[[p]][[i]]$pheno==phenotypes[g],]
            }
        }
    }
    
    
    
    ### Carry out analyses - INLA and Markov chains
    
    for (g in 1:nPhenotypes){
        
        for (p in 1:nProfiles){
            for (i in 1:(nInterims)){ if (phenoInt[[g]][i]==1) {
                
                # only carry out analysis if treatment efficacy undecided in this phenotype
                if (trialResult[[g]][[p]][m,i]=="Inconclusive"){
                    
                    # remove empty levels to avoid crashing the model fit
                    master[[g]][[p]][[i]]$OSFD2<-as.numeric(droplevels(as.factor(master[[g]][[p]][[i]]$OSFD2)))
                    ## INLA analysis of ordinal outcome
                    modelINLA[[g]][[p]][[i]] <- inla(OSFD2 ~  treat , family='pom',
                                                     data = master[[g]][[p]][[i]], 
                                                     control.fixed = list(mean = list( treat = 0), prec = 0.1),
                                                     control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, config = TRUE))
                    summary(modelINLA[[g]][[p]][[i]])
                    
                    
                    
                    
                    
                    # Store OR estimates and evaluate stopping rules
                    
                    post[[g]][[p]][[i]]<-modelINLA[[g]][[p]][[i]]$marginals.fixed$treat	
                    postCoefmean[[g]][[p]][m,i]<-summary(modelINLA[[g]][[p]][[i]])$fixed[2,1]
                    postCoefse[[g]][[p]][m,i]<-summary(modelINLA[[g]][[p]][[i]])$fixed[2,2]
                    stopEfficacy[[g]][[p]][m,i]<-1-inla.pmarginal(log(efficacyOR),post[[g]][[p]][[i]])>=efficacyPThresh
                    stopFutility[[g]][[p]][m,i]<-inla.pmarginal(log(futilityOR),post[[g]][[p]][[i]])>=futilityPThresh
                    
                    
                    # update trial result for this simulation 
                    
                    if(stopEfficacy[[g]][[p]][m,i]){
                        trialResult[[g]][[p]][m,i:(nInterims)]<-"Effective"
                        stoppingPoint[[g]][[p]][m]<-i
                        postOREff[[g]][[p]][[i]]<-c(postOREff[[g]][[p]][[i]],postCoefmean[[g]][[p]][m,i])
                    }
                    if(stopFutility[[g]][[p]][m,i]){
                        trialResult[[g]][[p]][m,i:(nInterims)]<-"Futile"
                        stoppingPoint[[g]][[p]][m]<-i
                        postORFut[[g]][[p]][[i]]<-c(postORFut[[g]][[p]][[i]],postCoefmean[[g]][[p]][m,i])
                    }
                    if(trialResult[[g]][[p]][m,i]=="Inconclusive" & sum(master[[g]][[p]][[i]]$treat==1)>=sampleCap[g]){
                        trialResult[[g]][[p]][m,i:(nInterims)]<-"Exceeded max"
                        stoppingPoint[[g]][[p]][m]<-i
                        postOREff[[g]][[p]][[i]]<-c(postOREff[[g]][[p]][[i]],"")
                        postORFut[[g]][[p]][[i]]<-c(postORFut[[g]][[p]][[i]],"")
                    }
                    if(i==stoppingPoint[[g]][[p]][m]){postOR[[g]][[p]]<-c(postOR[[g]][[p]],postCoefmean[[g]][[p]][m,i])}
                    
                }
                
                
            }}
            
        }
        
    }
}


#head(modelINLA)


# percentage of simulations that are successful / futile by ith interim
for (i in 1:(nInterims)){
    for (g in 1:nPhenotypes){
        for (p in 1:nProfiles){
            percentSuccess[[g]][[p]][i]<-sum(trialResult[[g]][[p]][,i]=="Effective")*100/nSims
            percentFutile[[g]][[p]][i]<-sum(trialResult[[g]][[p]][,i]=="Futile")*100/nSims
            percentExceededCap[[g]][[p]][i]<-sum(trialResult[[g]][[p]][,i]=="Exceeded max")*100/nSims
            percentDecided[[g]][[p]][i]<-percentSuccess[[g]][[p]][i] + percentFutile[[g]][[p]][i] + percentExceededCap[[g]][[p]][i] 
            percentUndecided[[g]][[p]][[i]]<-100-percentDecided[[g]][[p]][i]
            
            perDrugPhenoN[[g]][[p]][[i]]<-interimN[[g]][,2][pmin(i,stoppingPoint[[g]][[p]])]
            expectedN[[g]][[p]][i]<-mean(perDrugPhenoN[[g]][[p]][[i]])
        }
    }
}





# reshaping results for table/graph
# powerStats.csv is a table showing on a cumulative by-interim basis:
# probabilities of success/futility/exceeding cap, posterior OR values, and sample size summary stats
# for each phenotype/odds ratio.

# powerCurves.png shows graphically the probabilities of graduation or rejection (the latter includes exceeding the cap)
# for each phenotype/odds ratio.

# powerBars.png is similar to powerCurves but in a bar chart format and with probabilities of graduation / rejection / exceeding cap / remaining undecided.

ORLongVector<-list()
preORLongVector<-list()
ORShortVector<-list()
preORLongVectorAll<-list()
ORLongVectorAll<-list()

# profileVector<-list()
ORMatrix<-list()
powerMatrix<-list()
futilityMatrix<-list()
phenoVector<-list()
powerVector<-list()
futilityVector<-list()
capVector<-list()
undecidedVector<-list()


powerVectorAll<-vector()
futilityVectorAll<-vector()
capVectorAll<-vector()
undecidedVectorAll<-vector()
phenoVectorAll<-vector()

for (g in 1:nPhenotypes){
    ORShortVector[[g]]<-vector()
    ORLongVector[[g]]<-vector()
    preORLongVector[[g]]<-vector()
    ORLongVectorAll[[g]]<-vector()
    ORMatrix[[g]]<-matrix(nrow=0,ncol=nProfiles)
    powerVector[[g]]<-vector()
    powerMatrix[[g]]<-vector()
    futilityVector[[g]]<-vector()
    futilityMatrix[[g]]<-vector()
    capVector[[g]]<-vector()
    undecidedVector[[g]]<-vector()
    phenoVector[[g]]<-rep(phenotypes[g],nProfiles*nInterims)
    
    for (p in 1:nProfiles){
        ORShortVector[[g]]<-c(ORShortVector[[g]],exp(mean(postOR[[g]][[p]])))
        ORLongVector[[g]]<-c(ORLongVector[[g]],rep(exp(mean(postOR[[g]][[p]])),nInterims))
        preORLongVector[[g]]<-c(preORLongVector[[g]],rep(profileVector[p],nInterims))
        
        powerMatrix[[g]]<-cbind(powerMatrix[[g]],percentSuccess[[g]][[p]])
        powerVector[[g]]<-c(powerVector[[g]],percentSuccess[[g]][[p]])
        futilityMatrix[[g]]<-cbind(futilityMatrix[[g]],percentFutile[[g]][[p]]+percentExceededCap[[g]][[p]])
        
        futilityVector[[g]]<-c(futilityVector[[g]],percentFutile[[g]][[p]])
        capVector[[g]]<-c(capVector[[g]],percentExceededCap[[g]][[p]])
        undecidedVector[[g]]<-c(undecidedVector[[g]],percentUndecided[[g]][[p]])
    }
    for (i in 1:nInterims){
        ORMatrix[[g]]<-rbind(ORMatrix[[g]],t(ORShortVector[[g]]))
    }
    ORLongVectorAll<-c(ORLongVectorAll,ORLongVector[[g]])
    preORLongVectorAll<-c(preORLongVectorAll,preORLongVector[[g]])
    powerVectorAll<-c(powerVectorAll,powerVector[[g]])
    futilityVectorAll<-c(futilityVectorAll,futilityVector[[g]])	
    capVectorAll<-c(capVectorAll,capVector[[g]])
    undecidedVectorAll<-c(undecidedVectorAll,undecidedVector[[g]])
    phenoVectorAll<-c(phenoVectorAll,phenoVector[[g]])
}

ORLongVectorAll<-c(ORLongVectorAll,ORLongVectorAll)
intNumberVectorAll<-rep(c(1:nInterims),nPhenotypes*nProfiles)
phenoIntNumberVectorAll<-rep(paste0(phenoVectorAll,intNumberVectorAll),4)

resultVectorAll<-c(rep("graduated",length(powerVectorAll)),rep("rejected due to futility",length(powerVectorAll)),rep("exceeded cap",length(powerVectorAll)),rep("undecided",length(powerVectorAll)))
powerVectorAll<-c(powerVectorAll,futilityVectorAll,capVectorAll,undecidedVectorAll)
phenoVectorAll<-rep(phenoVectorAll,4)

powerCurves<-data.frame(cbind(ORLongVectorAll,powerVectorAll,resultVectorAll,phenoVectorAll,intNumberVectorAll,phenoIntNumberVectorAll))
names(powerCurves)<-c("OR","Probability","Result","Subphenotype","InterimNumber","phenoInterimNumber")
powerCurves$OR<-as.numeric(str_replace(powerCurves$OR,"_","."))
powerCurves$Probability<-as.numeric(powerCurves$Probability)

powerStats<-data.frame(cbind(rep(nSims,nInterims),rep(efficacyOR,nInterims),rep(efficacyPThresh,nInterims),rep(futilityOR,nInterims),rep(futilityPThresh,nInterims),phenoInt[[1]]*phenoN[[1]],phenoInt[[2]]*phenoN[[2]],ORMatrix[[1]],powerMatrix[[1]],futilityMatrix[[1]], ORMatrix[[2]],powerMatrix[[2]],futilityMatrix[[2]]))
names(powerStats)<-c("nSims","efficacyOR","efficacyPThresh","futilityOR","futilityPThresh",paste0("maxN_",phenotypes[1]),paste0("maxN_",phenotypes[2]),rep(paste0("OR_",phenotypes[1]),nProfiles),rep(paste0("power_",phenotypes[1]),nProfiles),rep(paste0("futility_",phenotypes[1]),nProfiles),rep(paste0("OR_",phenotypes[2]),nProfiles),rep(paste0("power_",phenotypes[2]),nProfiles),rep(paste0("futility_",phenotypes[2]),nProfiles))

write.csv(powerStats,"powerStats.csv")

ggplot(powerCurves[powerCurves$InterimNumber==nInterims,])+
    geom_line(aes(x = OR, y = Probability, color = Subphenotype,linetype = Result) )+
    geom_point(aes(x = OR, y = Probability, color = Subphenotype,shape = Result))+
    labs(x="Proportional odds ratio",y="Probability of result (percent)")

ggsave("powerCurves.png")

get_legend<-function(myggplot){
    tmp <- ggplot_gtable(ggplot_build(myggplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
}

powerCurves$Result<-factor(powerCurves$Result,levels=c("graduated","rejected due to futility","exceeded cap","undecided"))

palette<-c("green","red","orange","cyan")

plot1<-ggplot(powerCurves[powerCurves$InterimNumber==nInterims&powerCurves$Subphenotype=="1",], aes(fill=Result, y=Probability, x=OR)) + 
    geom_bar(position="fill", stat="identity")+theme(legend.position="none")+labs(x="OR in subgroup 1")+scale_fill_manual(values=palette)+ geom_text(position=position_fill(vjust=0.5),aes(y=Probability,x=OR,label = ifelse(Probability==0,"",round(Probability/100,2))), size=3)

plot2<-ggplot(powerCurves[powerCurves$InterimNumber==nInterims&powerCurves$Subphenotype=="2",], aes(fill=Result, y=Probability, x=OR)) + 
    geom_bar(position="fill", stat="identity")+labs(x="OR in subgroup 2")+scale_fill_manual(values=palette)+ geom_text(position=position_fill(vjust=0.5),aes(y=Probability,x=OR,label = ifelse(Probability==0,"",round(Probability/100,2))), size=3)



legend<-get_legend(plot2)

plot2<-plot2 +theme(legend.position="none")



gg <- arrangeGrob(plot1,plot2, legend, ncol=3, widths=c(1.5,1.5,1)) #generates graph
ggsave("powerBars.png",gg)
