model_string='
data{
int N1;
int y1[N1, 4]; 
}
parameters{
real<lower=0,upper=1> SeELISA;
real<lower=0,upper=1> SpELISA;
real<lower=0,upper=1> Seifc;
real<lower=0,upper=1> Spifc;
real<lower=0,upper=1> phi1;
real<lower=0> sigma;
real alpha;
real U1[N1];
real<lower=0,upper=1> z1[N1];
}
transformed parameters{
real pi1[N1];
simplex[4] p1[N1];

for (i in 1:N1){
pi1[i] = inv_logit((alpha + U1[i]))*z1[i];
p1[i,1] = pi1[i]*SeELISA*Seifc + (1-pi1[i])*(1-SpELISA)*(1-Spifc);
p1[i,2] = pi1[i]*SeELISA*(1-Seifc) + (1-pi1[i])*(1-SpELISA)*Spifc;
p1[i,3] = pi1[i]*(1-SeELISA)*Seifc + (1-pi1[i])*SpELISA*(1-Spifc);
p1[i,4] = pi1[i]*(1-SeELISA)*(1-Seifc) + (1-pi1[i])*SpELISA*Spifc;
}
}
model{
//
SeELISA ~ beta(15, 2);
Seifc ~ beta(15, 2);
SpELISA ~ beta(15, 2);
Spifc ~ beta(15, 2);
//
phi1 ~ beta(30,100);
alpha ~ normal(0,5);
sigma ~ cauchy(0, 5);
//
for (j in 1:N1){
U1[j] ~ normal( 0 , sigma);
z1[j] ~ beta(phi1, 1-phi1);
}

for(i in 1:N1){
y1[i] ~ multinomial(p1[i]);
}
}
'

library(rstan)


sta_dso2 <- stan_model(model_code = model_string)
load('output\\datNZ.RData')
dataList <- list(y1 = as.matrix(dat), N1 = 18)


stanFit <- sampling(object = sta_dso2, 
                    data = dataList, 
                    chains = 1, iter = 9000, 
                    warmup = 1000, thin = 1)
