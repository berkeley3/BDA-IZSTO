## ----stan2, cache=T------------------------------------------
library(rstan)
modelString <- " 
data{
// First declare all variables in the data block
int<lower=0> N1;
int<lower=0> y1;
int<lower=0> N2;
int<lower=0> y2;
}
parameters{
real<lower=0, upper=1> theta1;
real<lower=0, upper=1> theta2;
}
model{
theta1 ~ beta(1,1);
y1 ~ binomial(N1, theta1);
theta2 ~ beta(1,1);
y2 ~ binomial(N2, theta2);
}
generated quantities{
real theta;
theta = theta2 - theta1;
}
" 


stanDso <- stan_model(model_code = modelString_2)


dataList = list(N1 =52, N2 =48, y1=43, y2 = 44)

stanFit <- sampling(object = stanDso, 
                    data = dataList, 
                    chains = 1, iter = 9000, warmup = 1000, thin = 1)

samples <- extract(stanFit)
thetas_post <- samples[['theta']]


## ----dependson="stan2", fig.height=3, fig.width=3-------------------------

plot(density(thetas_post), main= "posterior probability", 
     xlab=expression(theta[2]-theta[1]) )


## ----dependson="stan2", fig.height=3, fig.width=3-------------------------

quantile(thetas_post, probs=c(.025, .5, .975))

