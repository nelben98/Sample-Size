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
seed_set =1

#
#
## fit the regression:
## 
## 

fit = try(INLA::inla(formula=model, data=data, family=family, control.compute=list(openmp.strategy="small"),num.threads=1,
                   control.family=list(control.link=list(model=link)),
                   verbose=FALSE)
    )

#  Now try a dataset until it breaks
#
#
















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



#####
######
######
###### Sequence to analyse the treatment and OSFD2 distributions
######
######




bind_rows(scenario2$H1$trial, .id = "trial")


info_inter<-matrix(NA, nrow = 962*10,ncol=4)
colnames(info_inter) <-c('treat','OSFD','pheno','OSFD2','trial')

info_inter_BATSS<-data.frame(matrix(NA, nrow = 50*10,ncol=4))
colnames(info_inter_BATSS) <-c("UC", "Simvastatin",'trial','cat')

for (i in 1:50){
    if (i == 1) {
        info_inter_BATSS[1:10,c("UC")] <- data.frame(scenario2$H1$trial[[1]]$data) |> filter(treatment=='UC') |> select(Freq)
        info_inter_BATSS[1:10,c( "Simvastatin")] <- data.frame(scenario2$H1$trial[[1]]$data) |> filter(treatment=='Simvastatin') |> select(Freq)
        info_inter_BATSS[1:10,'trial']<-i
        info_inter_BATSS[1:10,'cat']<-data.frame(scenario2$H1$trial[[1]]$data) |> filter(treatment=='Simvastatin') |> select(y)
        
        }
    else if (i > 1) {
        info_inter_BATSS[c((i-1)*10+1):c((i*10)),c("UC")] <-           data.frame(scenario2$H1$trial[[i]]$data) |> filter(treatment=='UC') |> select(Freq)
        info_inter_BATSS[c((i-1)*10+1):c((i*10)),c( "Simvastatin")] <- data.frame(scenario2$H1$trial[[i]]$data) |> filter(treatment=='Simvastatin') |> select(Freq)
        info_inter_BATSS[c((i-1)*10+1):c((i*10)),c( "cat")] <- data.frame(scenario2$H1$trial[[i]]$data) |> filter(treatment=='Simvastatin') |> select(y)
        info_inter_BATSS[c((i-1)*10+1):c((i*10)),'trial']<-i
        
        }
}

#if 0 then control if 1 then treatment
dist_issues<-
        cbind(
        info_inter_BATSS |> 
            mutate(nrow = row_number(),
                   Simva= as.factor(Simvastatin)) |>
            select(cat,Simvastatin)|>
           group_by(cat) %>% 
           summarise(trt_batss = sum(Simvastatin)),
        bind_rows(interimDataAll[[1]], .id = "trial") |>
            filter(treat==1) |>
            mutate(nrow = row_number()) |> 
            select(trial,treat,OSFD2,nrow) |> 
            pivot_wider(names_from = treat
                       ,values_from = OSFD2)|> 
            mutate(trt = as.factor(`1`)) |>
            count(trt) |> 
            rename(trt_algo=n) |> select(!trt)
        )|>
    left_join(
        cbind(
        info_inter_BATSS |> 
            mutate(nrow = row_number()) |>
            select(cat,UC)|>
            group_by(cat) %>% 
            summarise(UC_batss = sum(UC)),
        bind_rows(interimDataAll[[1]], .id = "trial") |>
            filter(treat==0) |>
            mutate(nrow = row_number()) |> 
            select(trial,treat,OSFD2,nrow) |> 
            pivot_wider(names_from = treat
                       ,values_from = OSFD2)|> 
            mutate(UC = as.factor(`0`)) |>
            count(UC) |> 
            rename(UC_algo=n) |>select(!UC)
        ), by = join_by(cat)
    ) 



ggplot(data=dist_issues, aes(x=cat)) +
    geom_line(aes(y=trt_batss ,col='Batss'))+
    geom_line(aes(y=trt_algo  ,col='Algo'))+
    theme_minimal() +
    scale_x_continuous(name="Days", n.breaks = 10)


ggplot(data=dist_issues, aes(x=cat)) +
    geom_line(aes(y=UC_batss  ,col='Batss'))+
    geom_line(aes(y=UC_algo  ,col='Algo'))+
    theme_minimal() +
    scale_x_continuous(name="Days", n.breaks = 10)
