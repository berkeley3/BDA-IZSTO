#  model  2 (random effects+no baseline)----------------
crossover <- read.csv('crossover.csv')

## import library----------
library(BRugs)

# write model----------
cat("model{ 
    for (i in 1:N){
    # --------------------
    # model for oed
    # --------------------
    oed[i]~dnorm( mu.oed[i], tau.oed[1])
    mu.oed[i] <- b.oed[1] + b.oed[2] *(period[i]-1)+ b.oed[3] *(treatment[i]-1) + a.oed[ patient[i] ] 
    # --------------------
    # model for dbp
    # --------------------
    dbp[i]~dnorm( mu.dbp[i], tau.dbp[1])
    mu.dbp[i] <- b.dbp[1] + b.dbp[2] *(period[i]-1)+ b.dbp[3] *(treatment[i]-1) + a.dbp[ patient[i] ] 
    }
    for (i in 1:n){ 
    # 
    # Non-informative priors for individual/patients fixed effects		
    a.oed[i]~dnorm( 0.0, tau.oed[2])
    a.dbp[i]~dnorm( 0.0, tau.dbp[2])
    }
    for (i in 1:p){
    b.oed[i]~dnorm( 0.0, 0.001)
    b.dbp[i]~dnorm( 0.0, 0.001)
    }
    tau.oed[1]~dgamma( 0.001,0.001)
    tau.oed[2]~dgamma( 0.001,0.001)
    tau.dbp[1]~dgamma( 0.001,0.001)
    tau.dbp[2]~dgamma( 0.001,0.001)
    #
    #
    s2[1]<-1/tau.oed[1]
    s2[2]<-1/tau.oed[2]
    s2[3]<-1/tau.dbp[1]
    s2[4]<-1/tau.dbp[2]
    r[1] <- s2[2]/(s2[2]+s2[1])
    r[2] <- s2[4]/(s2[4]+s2[3])
    
    for( i in 1:N ){
    res1[i] <- oed[i] - mu.oed[i]
    res2[i] <- dbp[i] - mu.dbp[i]
    }
    R[1] <- 1 - pow( sd(res1[1:N])/sd(oed[1:N]), 2)
    R[2] <- 1 - pow( sd(res2[1:N])/sd(dbp[1:N]), 2)
    }", file ='model2.txt')

# model Check------------

modelCheck('model2.txt')

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
inits1 = list( b.oed=c(0,0,0), b.dbp=c(0,0,0), tau.oed=c(1,1), tau.dbp=c(1,1))
BRugs::bugsInits(list(inits1), numChains=1, fileName='inits2.txt')
modelInits('inits2.txt')
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

