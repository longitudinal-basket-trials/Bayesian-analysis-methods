library(R2jags)

######################
### Implementation ###
######################

######################
### Load Data ########
######################

# "Example_Simulated_Data.Rdata" = example of a (single) simulated data set with treatment effect parameters: 
#  true_theta = c(0.59, 1.17, 1.02, 0.95, 0.13, 0.75)
#  Note: dataset simulated according to the procedure specified in 'Bayesian borrowing for basket trials with longitudinal outcomes'

load(file="Example_Simulated_Data.Rdata") 
head(df)
##################
## Stratified-L ##
##################

# clinically relevant threshold for treatment effect
cutoff.theta = c(0.25, 0.30)

# L2 mu prior hyperparameters
Prior.theta = 0
Prior.sd.theta = 1

Prior.g00 = c(20,3)
Prior.g01 = c(0, 1) 
Prior.g10 = c(0, 1) 

# L2 SD prior hyperparameters

zero=0
zero[1] <- 0
zero[2] <- 0
R=matrix(nrow=2,ncol=2) # N.B. if using Wishart prior for covariance, R is hyperparameters for inverse covariance matrix (see model file)
R[1,1] <- 1
R[1,2] <- 0 # prior hyperparameters for (inverse) covariance between intercept and slope (if using Wishart prior)
R[2,1] <- 0 # prior hyperparameters for (inverse) covariance between intercept and slope (if using Wishart prior)
R[2,2] <- 1

# L1 SD prior hyperparameters
PriorL1sigHN=1

# set initial values for MCMC
inits <- function(){
  list(
    sig1 = 1,
    Tau00 = 1, Tau10 = 1
  )
}

# parameters to return
parameters <- c("theta", "pCat")

pCat = numeric(6)
nthetaM = list()

set.seed(123)

kMod=length(unique(df$Module.long))
Module.long=df$Module.long

# run Stratified-L model
for(k in 1:kMod){
  
  y = df$y.true[Module.long == k]
  NObs = length(y)
  totaln = NObs/4
  patientID = c(1:totaln)
  patientID = (rep(patientID,each=4))
  t1 = df$timet[Module.long == k]
  Trt1 = df$true.trt.long[Module.long == k]
  Trt = Trt1[seq(1, length(Trt1), 4)]
  
  data <- list("NObs", "y", "t1", "Trt", 
               "totaln", "patientID",
               "Prior.theta", "Prior.sd.theta", 
               "Prior.g00", "Prior.g01", "Prior.g10",
               "zero", "R",
               "PriorL1sigHN", "cutoff.theta")
  
  Strat_L <- jags(data=data, inits, parameters.to.save=parameters, 
                  model.file= "Stratified-L.txt", 
                  n.chains = 4,
                  n.burnin = 13000, 
                  n.iter = 30000)
  
  nthetaM[[k]] = Strat_L$BUGSoutput$summary
  pCat[k] = Strat_L$BUGSoutput$summary[2,1] # p(theta>0.25)
  Strat_L_list<-list("summary"=nthetaM, "pCat"=pCat)
  
}

##############
### BHM-L ####
##############
y = df$y.true

kMod=length(unique(df$Module.long))
totaln=length(unique(df$patientID))
NObs=length(y)
Module.long=df$Module.long
Module=Module.long[seq(1, length(Module.long), 4)]
patientID=c(rep((1:totaln), each=4))

t1=df$timet
Trt.long= df$true.trt.long
Trt = Trt.long[seq(1, length(Trt.long), 4)]

cutoff.theta = c(0.25, 0.30)

# prior hyperparameters for L3 mus 
Prior.g00 = c(20,3)
Prior.g01 = c(0, 1) 
Prior.g10 = c(0, 1) 
Prior.g11 = c(0, 1) 

# prior hyperparameters for L3 SDs
PriorL3sigHN = c(1, 1, 1, 0.125)

# L2 prior SD hyperparameters
zero=0
zero[1] <- 0
zero[2] <- 0
R=matrix(nrow=2,ncol=2) 
R[1,1] <- 1
R[1,2] <- 0
R[2,1] <- 0 
R[2,2] <- 1

# L1 prior SD hyperparameters
PriorL1sigHN=1

# data list
data <- list("y", 
             "Module", 
             "t1", 
             "Trt", 
             "patientID",
             "kMod",
             "totaln",
             "NObs",
             "Prior.g00", "Prior.g01", "Prior.g10", "Prior.g11", "PriorL3sigHN",
             "zero", "R",
             "PriorL1sigHN",
             "cutoff.theta")

# set initial values for MCMC
inits <- function(){
  list(
    sigma00 = 1, sigma01 = 1, sigma10 = 1, sigma11 = 1,
    Tau00=1, Tau10=1,
    sig1=1
  )
}

# parameters to return
parameters <- c("theta", "pCat")

# run BHM-L model
set.seed(123)
BHM_L <- jags(data=data, 
                 inits, 
                 parameters.to.save=parameters, 
                 model.file="BHM-L.txt",
                 n.chains = 4,
                 n.burnin = 13000, 
                 n.iter = 30000) 

##############
### EXNEX-L ##
##############
y = df$y.true

kMod=length(unique(df$Module.long))
totaln=length(unique(df$patientID))
NObs=length(y)
Module.long=df$Module.long
Module=Module.long[seq(1, length(Module.long), 4)]
patientID=c(rep((1:totaln), each=4))

t1=df$timet
Trt.long= df$true.trt.long
Trt = Trt.long[seq(1, length(Trt.long), 4)]

nex.theta = 0
nex.sig = 1

pMix = c(0.5, 0.5)

Prior.mu.theta = c(0, 1)
Prior.tau.HN = 0.125

cutoff.theta = c(0.25, 0.30)

# prior hyperparameters for L3 mus 
Prior.g00 = c(20,3)
Prior.g01 = c(0, 1) 
Prior.g10 = c(0, 1)  

# prior hyperparameters for L3 SD
PriorL3sigHN = c(1, 1, 1)  

# L2 prior SD hyperparameters
zero=0
zero[1] <- 0
zero[2] <- 0
R=matrix(nrow=2,ncol=2) 
R[1,1] <- 1
R[1,2] <- 0
R[2,1] <- 0 
R[2,2] <- 1

# L1 prior SD hyperparameters
PriorL1sigHN=1

# data list
data <- list("y", 
             "Module", 
             "t1", 
             "Trt", 
             "patientID",
             "kMod",
             "totaln",
             "NObs",
             "Prior.mu.theta", "Prior.tau.HN",
             "Prior.g00", "Prior.g01", "Prior.g10", "PriorL3sigHN",
             "zero", "R",
             "PriorL1sigHN",
             "cutoff.theta",
             "nex.theta", "nex.sig", "pMix")

# set initial values for MCMC
inits <- function(){
  list(
    sigma00 = 1, sigma01 = 1, sigma10 = 1, 
    Tau00=1, Tau10=1,
    sig1=1
  )
}

# parameters to return
parameters <- c("theta", "pCat")

# run EXNEX-L model
set.seed(123)
EXNEX_L <- jags(data=data, inits, parameters.to.save=parameters, 
                model.file="EXNEX-L.txt", 
                n.chains = 4,
                n.burnin = 13000, 
                n.iter = 30000)


#############
### HD-L ####
#############
NObs=length(df$y.true)
y1=df$y.true
y = array(0, dim = c(2, NObs))
y[1,] = y1
y[2,] = y1

kMod=length(unique(df$Module.long))

Module.long=df$Module.long
Module=Module.long[seq(1, length(Module.long), 4)]

totaln=length(unique(df$patientID))
patientID=c(rep((1:totaln), each=4))

t1 = array(0, dim = c(2, NObs))
t1[1, ] = df$timet
t1[2, ] = df$timet

Trt.long= df$true.trt.long
Trt = array(dim = c(2, totaln))
Trt[1, ] = Trt.long[seq(1, length(Trt.long), 4)]
Trt[2, ] = Trt.long[seq(1, length(Trt.long), 4)]

s0 = 0.15

lSlab = 0.01
uSlab = 1

spike = 100

Prior.theta = rep(0, kMod)
Prior.sd.theta = rep(1, kMod)

cutoff.theta = c(0.25, 0.30)

# prior hyperparameters for L3 mus

Prior.g00 = rep(20, kMod)
Prior.g01 = rep(0, kMod)
Prior.g10 = rep(0, kMod)

prior.sd.g00=rep(3, kMod)
prior.sd.g01=rep(1, kMod)
prior.sd.g10=rep(1, kMod)

# prior hyperparameters for L2 SD

zero=0
zero[1] <- 0
zero[2] <- 0
R=matrix(nrow=2,ncol=2) 
R[1,1] <- 1
R[1,2] <- 0
R[2,1] <- 0 
R[2,2] <- 1

# prior hyperparameters for L1 SD

PriorL1sigHN=rep(1, kMod)

#data list
data <- list("Module.long",
              "prior.sd.g00","prior.sd.g01","prior.sd.g10",
              "kMod", "NObs", "totaln", "y", "patientID", "Module", "t1", "Trt",
              "Prior.theta", "Prior.sd.theta",
              "s0", "lSlab", "uSlab", "spike",
              "Prior.g00", "Prior.g01", "Prior.g10",
              "zero", "R",
              "PriorL1sigHN",
              "cutoff.theta")

# set initial values for MCMC
inits <- function(){
  list(
    Tau00=1, Tau10=1,
    sig1=rep(1, kMod)
  )}

# parameters to return
parameters <- c("theta", "pCat")

# run HD-L model
set.seed(123)
HD_L <- jags(data=data,
                inits,
                parameters.to.save=parameters,
                model.file="HD-L.txt",
                n.chains = 4,
                n.burnin = 13000,
                n.iter = 30000)
