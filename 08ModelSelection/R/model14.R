#  model  2 (random effects+ baseline)----------------
crossover.base <- read.csv('R\\crossover.base.csv')

## import library----------
library(BRugs)

# write model----------
cat("model{ 
    for (i in 1:N){
    # --------------------
    # model for oed
    # --------------------
    oed[i]~dnorm( mu.oed[i], tau.oed[1])
    mu.oed[i] <- b.oed[1] + b.oed[2] *(period[i]-1)+ b.oed[3] *(treatment[i]-1) + g[1]*a.oed[ patient[i] ] + (1-g[1])*a.oed.fixed[ patient[i] ] + g[2]*b.oed[4]*oedbase[i]
    # --------------------
    # model for dbp
    # --------------------
    dbp[i]~dnorm( mu.dbp[i], tau.dbp[1])
    mu.dbp[i] <- b.dbp[1] + b.dbp[2] *(period[i]-1)+ b.dbp[3] *(treatment[i]-1) + g[1]*a.dbp[ patient[i] ]+(1-g[1])*a.dbp.fixed[ patient[i] ] + g[2]*b.dbp[4]*dbpbase[i]
    }
    for (i in 1:n){ 
    # 
    # Hyper priors for individual/patients random effects
    a.oed[i]~dnorm( 0.0, tau.oed[2])
    a.dbp[i]~dnorm( 0.0, tau.dbp[2])
    # 
    # Non-informative priors for individual/patients fixed effects		
    a.oed.fixed[i]~dnorm( 0.0, 0.001)
    a.dbp.fixed[i]~dnorm( 0.0, 0.001)
    }
    for (i in 1:p){
    b.oed[i]~dnorm( 0.0, 0.001)
    b.dbp[i]~dnorm( 0.0, 0.001)
    }
    tau.oed[1]~dgamma( 0.001,0.001)
    tau.oed[2]~dgamma( a, a)
    tau.dbp[1]~dgamma( 0.001,0.001)
    tau.dbp[2]~dgamma( a, a)
    a <- g[1]*0.001 + (1-g[1])*1
    #
    #
    s2[1]<-1/tau.oed[1]
    s2[2]<-1/tau.oed[2]
    s2[3]<-1/tau.dbp[1]
    s2[4]<-1/tau.dbp[2]
    
    for( i in 1:N ){
    res1[i] <- oed[i] - mu.oed[i]
    res2[i] <- dbp[i] - mu.dbp[i]
    }
    R[1] <- 1 - pow( sd(res1[1:N])/sd(oed[1:N]), 2)
    R[2] <- 1 - pow( sd(res2[1:N])/sd(dbp[1:N]), 2)
    }", file ='model14.txt')

# model Check------------

modelCheck('model14.txt')

# load data------

datalist <- list(N=nrow(crossover.base),
                 n = unique(length(crossover.base$patient)),
                 p = 4, 
                 g=c(1,1),
                 oed = crossover.base$oed,
                 dbp = crossover.base$dbp,
                 period = crossover.base$period,
                 treatment = crossover.base$treatment,
                 patient = crossover.base$patient,
                 oedbase = crossover.base$oedbase,
                 dbpbase = crossover.base$dbpbase)


modelData(bugsData( datalist))

# compile model----
modelCompile()

#inits-------------
inits1 = list( b.oed=c(0,0,0,0), b.dbp=c(0,0,0,0), tau.oed=c(1,1), tau.dbp=c(1,1))
BRugs::bugsInits(list(inits1), numChains=1, fileName='inits14.txt')
modelInits('inits14.txt')
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

