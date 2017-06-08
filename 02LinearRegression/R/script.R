
#-- Hierarchical Modeling: example of meta-analysis: BRugs------------------

### specify model
modelstring="
model
{
  for (k in 1:K1){ 
         logor[k] <- log(or[k])
         selogor[k] <- log(U[k]/L[k])/(2*1.96)
		 precision.logor[k] <- 1/pow(  selogor[k], 2) 
         logor[k] ~ dnorm( theta[k],  precision.logor[k] ) 
         theta[k]~dnorm( mu.theta, tau.theta )
         OR[k] <- exp(theta[k]) 
    }

    for (k in 1:K2){ 
        for (i in 1:2){ 
        	N[i,k] <- Y[i,1,k]+Y[i,2,k]
	        Y[i,1,k] ~ dbin( p[i,k], N[i,k] )
    	    logit(p[i,k]) <- a[k] + theta[K1+k] * equals(i,1)
		}
        theta[K1+k] ~ dnorm( mu.theta, tau.theta )
        OR[K1+k] <- exp(theta[K1+k]) 
    }
    # priors 
    for( k in 1:3) { a[k] ~ dnorm( 0.0, 0.001) } 
    mu.theta  ~ dnorm( 0.0, 0.001)
    tau.theta ~ dgamma( 0.001, 0.001)

    s2.theta <- 1/tau.theta

    s.theta <- sqrt( s2.theta )
    OR.total <- exp(mu.theta)
  
}
"
.temp = file("myModel.txt","w")  
writeLines(modelstring,con=.temp) 
close(.temp)


###  Specify the data as a list
library(BRugs)
datalist <- list( K1=7, K2=3, 
              or=c(3.89, 3.97, 3.88, 17.47, 5.35, 9.1, 3.41), 
              L=c(0.92, 2.2, 2.47, 14.24, 2.44, 5.57, 2.94),
              U=c(16.3, 7.16, 6.08, 21.43, 11.74, 14.86, 3.96), 
              Y=structure(.Data=c(49, 33, 29958, 70186, 12, 89, 0, 118, 29, 171,  4, 81), .Dim = c(2, 2, 3))
)

modelCheck('myModel.txt')
modelData(bugsData( datalist))

modelCompile()
inits1 <- list(mu.theta=3.5, tau.theta=1, a=c(0,0,0), theta=rep(1,10))
#BRugs::bugsInits(list(inits1), numChains=1, fileName='inits1.txt')
#modelInits('inits1.txt')
modelInits(bugsInits(list(inits1), numChains=1))
# eventually for non-parents parameters: modelGenInits()
#modelGenInits()
# burn-in
modelUpdate(20000)     
# set parameters and update model
parameters = c( "mu.theta" , "OR.total", "tau.theta", "OR" )     # The parameter(s) to be monitored.

samplesSet(parameters)     
dicSet()
modelUpdate(50000)       

#sample stats
samplesStats("*") 
samplesStats("OR.total", thin=10) 

samplesDensity("*")              # plot the densities,
samplesHistory("*", mfrow = c(2, 2)) # plot the chain,

samplesDensity("tau")              # plot the densities,
samplesDensity("k", thin=50)        plot the densities,

samplesBgr("*")             # plot the bgr statistics
samplesAutoC("*", 1)     # plot autocorrelation

dicStats()


#-- Hierarchical Modeling: example of meta-analysis: JAGS------------------

### specify model
modelstring="
model
{
  for (k in 1:K1){ 
  selogor[k] <- log(U[k]/L[k])/(2*1.96)
  precision.logor[k] <- 1/pow(  selogor[k], 2) 
  logor[k] ~ dnorm( theta[k],  precision.logor[k] ) 
  theta[k]~dnorm( mu.theta, tau.theta )
  OR[k] <- exp(theta[k]) 
  }
  
  for (k in 1:K2){ 
  for (i in 1:2){ 
  Y[i,1,k] ~ dbin( p[i,k], N[i,k] )
  logit(p[i,k]) <- a[k] + theta[K1+k] * equals(i,1)
  }
  theta[K1+k] ~ dnorm( mu.theta, tau.theta )
  OR[K1+k] <- exp(theta[K1+k]) 
  }
  # priors 
  for( k in 1:3) { a[k] ~ dnorm( 0.0, 0.001) } 
  mu.theta  ~ dnorm( 0.0, 0.001)
  tau.theta ~ dgamma( 0.001, 0.001)
  
  s2.theta <- 1/tau.theta
  
  s.theta <- sqrt( s2.theta )
  OR.total <- exp(mu.theta)
  
}
"
.temp = file("myModel.txt","w")  
writeLines(modelstring,con=.temp) 
close(.temp)


###  Specify the data as a list
library(rjags)
datalist <- list( K1=7, K2=3, 
                  logor=log(c(3.89, 3.97, 3.88, 17.47, 5.35, 9.1, 3.41)), 
                  L=c(0.92, 2.2, 2.47, 14.24, 2.44, 5.57, 2.94),
                  U=c(16.3, 7.16, 6.08, 21.43, 11.74, 14.86, 3.96), 
                  Y=structure(.Data=c(49, 33, 29958, 70186, 12, 89, 0, 118, 29, 171,  4, 81), .Dim = c(2, 2, 3)),
                   N=structure(.Data=c(49+29958,33+70186,12,89+118,29+4,171+81), .Dim=c(2,3)))            


inits1 <- list(mu.theta=3.5, tau.theta=1, a=c(0,0,0), theta=rep(1,10))
parameters <- c( "mu.theta" , "OR.total", "tau.theta", "OR" )     # The parameter(s) to be monitored.


adaptSteps <- 500              # Number of steps to "tune" the samplers.
burnInSteps <- 5000            # Number of steps to "burn-in" the samplers.
nChains <- 1                  # Number of chains to run.
numSavedSteps <- 50000           # Total number of steps in chains to save.
thinSteps <- 1                   # Number of steps to "thin" (1=keep every step).
nIter <- ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.
# Create, initialize, and adapt the model:
jagsModel <- jags.model( "myModel.txt" , data=datalist ,  inits=inits1 ,
                        n.chains=nChains , n.adapt=adaptSteps )

cat( "Burning in the MCMC chain...\n" )
update( jagsModel , n.iter=burnInSteps )
# The saved MCMC chain:
cat( "Sampling final MCMC chain...\n" )
codaSamples = coda.samples( jagsModel , variable.names=parameters ,
                            n.iter=nIter , thin=thinSteps )

summary(codaSamples)
plot(codaSamples)
mcmcChain <- as.matrix(codaSamples)
OR.total <- mcmcChain[,'OR.total']
plot(OR.total)

