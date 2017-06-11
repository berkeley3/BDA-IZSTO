#  model  1 (fixed effects+no baseline)----------------
crossover <- read.csv('crossover.csv')

## import library----------
library(BRugs)

# write model----------
cat("model{ 
    for (i in 1:N){
    # --------------------
    # model for oed
    # --------------------
    oed[i]~dnorm( mu.oed[i], tau.oed)
    mu.oed[i] <- b.oed[1] + b.oed[2] *(period[i]-1)+ b.oed[3] *(treatment[i]-1) + a.oed.fixed[ patient[i] ] 
    # --------------------
    # model for dbp
    # --------------------
    dbp[i]~dnorm( mu.dbp[i], tau.dbp)
    mu.dbp[i] <- b.dbp[1] + b.dbp[2] *(period[i]-1)+ b.dbp[3] *(treatment[i]-1) + a.dbp.fixed[ patient[i] ] 
    }
    for (i in 1:n){ 
    # 
    # Non-informative priors for individual/patients fixed effects		
    a.oed.fixed[i]~dnorm( 0.0, 0.001)
    a.dbp.fixed[i]~dnorm( 0.0, 0.001)
    }
    for (i in 1:p){
    b.oed[i]~dnorm( 0.0, 0.001)
    b.dbp[i]~dnorm( 0.0, 0.001)
    }
    tau.oed~dgamma( 0.001,0.001)
    tau.dbp~dgamma( 0.001,0.001)
    #
    #
    s2[1]<-1/tau.oed
    s2[2]<-1/tau.dbp
    
    for( i in 1:N ){
    res1[i] <- oed[i] - mu.oed[i]
    res2[i] <- dbp[i] - mu.dbp[i]
    }
    R[1] <- 1 - pow( sd(res1[1:N])/sd(oed[1:N]), 2)
    R[2] <- 1 - pow( sd(res2[1:N])/sd(dbp[1:N]), 2)
    }", file ='model1.txt')

# model Check------------

modelCheck('model1.txt')

# load data------

datalist <- list(N=nrow(crossover),
                 n = unique(length(crossover$patient)),
                 p = 3,
                 oed = crossover$oed,
                 dbp = crossover$dbp,
                 period = crossover$period,
                 treatment = crossover$treatment,
                 patient = crossover$patient)

modelData(bugsData( datalist))

# compile model----
modelCompile()

#inits-------------
inits1 = list( b.oed=c(0,0,0), b.dbp=c(0,0,0), tau.oed=1, tau.dbp=1)
BRugs::bugsInits(list(inits1), numChains=1, fileName='inits1.txt')
modelInits('inits1.txt')
modelGenInits()

# bur-in------
modelUpdate(20000)     
# set parameters----
parameters <- c('tau.oed', 'tau.dbp', 'b.oed','b.dbp')
samplesSet(parameters)     
dicSet()
modelUpdate(50000)       

dicStats()
samplesStats('*')

### JAGS
library(rjags)
# write model----------
cat("model{ 
    for (i in 1:N){
    # --------------------
    # model for oed
    # --------------------
    oed[i]~dnorm( mu.oed[i], tau.oed)
    mu.oed[i] <- b.oed[1] + b.oed[2] *(period[i]-1)+ b.oed[3] *(treatment[i]-1) + a.oed.fixed[ patient[i] ] 
    # --------------------
    # model for dbp
    # --------------------
    dbp[i]~dnorm( mu.dbp[i], tau.dbp)
    mu.dbp[i] <- b.dbp[1] + b.dbp[2] *(period[i]-1)+ b.dbp[3] *(treatment[i]-1) + a.dbp.fixed[ patient[i] ] 
    }
    for (i in 1:n){ 
    # 
    # Non-informative priors for individual/patients fixed effects		
    a.oed.fixed[i]~dnorm( 0.0, 0.001)
    a.dbp.fixed[i]~dnorm( 0.0, 0.001)
    }
    for (i in 1:p){
    b.oed[i]~dnorm( 0.0, 0.001)
    b.dbp[i]~dnorm( 0.0, 0.001)
    }
    tau.oed~dgamma( 0.001,0.001)
    tau.dbp~dgamma( 0.001,0.001)
    #
    #
    s2[1]<-1/tau.oed
    s2[2]<-1/tau.dbp
    
    for( i in 1:N ){
    res1[i] <- oed[i] - mu.oed[i]
    res2[i] <- dbp[i] - mu.dbp[i]
    }
    R[1] <- 1 - pow( sd(res1[1:N])/sd(oed[1:N]), 2)
    R[2] <- 1 - pow( sd(res2[1:N])/sd(dbp[1:N]), 2)
    }", file ='model1jags.txt')

## write data as list------
datalist <- list(N=nrow(crossover),
                 n = unique(length(crossover$patient)),
                 p = 3,
                 oed = crossover$oed,
                 dbp = crossover$dbp,
                 period = crossover$period,
                 treatment = crossover$treatment,
                 patient = crossover$patient)

## initialize values-----

inits1 = list( b.oed=c(0,0,0), b.dbp=c(0,0,0), tau.oed=1, tau.dbp=1)
inits2 = list( b.oed=c(1,1,1), b.dbp=c(1,1,1), tau.oed=2, tau.dbp=2)

# set parameter(s) to be monitored-----
parameters <- c('tau.oed', 'tau.dbp', 'b.oed','b.dbp')

adaptSteps = 500              # Number of steps to "tune" the samplers.
burnInSteps = 5000            # Number of steps to "burn-in" the samplers.
nChains = 2                  # Number of chains to run.
numSavedSteps=50000           # Total number of steps in chains to save.
thinSteps=1                   # Number of steps to "thin" (1=keep every step).
nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.
# Create, initialize, and adapt the model:
jagsModel = jags.model( "model1jags.txt" , data=datalist ,  inits=list(inits1,inits2),
                        n.chains=nChains , n.adapt=adaptSteps )


cat( "Burning in the MCMC chain...\n" )
update( jagsModel , n.iter=burnInSteps )
# The saved MCMC chain:
cat( "Sampling final MCMC chain...\n" )
codaSamples = coda.samples( jagsModel , variable.names=parameters ,
                            n.iter=nIter , thin=thinSteps )

dic.samples(jagsModel , variable.names=parameters ,
            n.iter=nIter , thin=thinSteps )