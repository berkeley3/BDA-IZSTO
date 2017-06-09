load('R\\cortbowl.RData')
cortbowl$totCortlog <- log(cortbowl$totCort)
cortbowl$Ring <- as.integer(as.factor(cortbowl$Ring)) 
cortbowl$Age_z <- scale(cortbowl$Age)

mod_mat <- as.data.frame(model.matrix(~ Implant + days + Implant:days, data=cortbowl))
data_list <- list(N = nrow(cortbowl),
                  totCortlog = cortbowl$totCortlog,
                  ImplantP = mod_mat$ImplantP,
                  days20 = mod_mat$days20,
                  daysbefore = mod_mat$daysbefore,
                  Age_z = as.numeric(cortbowl$Age_z),
                  Ring = cortbowl$Ring,
                  ImplantP_X_days20 = mod_mat$`ImplantP:days20`,
                  ImplantP_X_daysbefore = mod_mat$`ImplantP:daysbefore`,
                  N_Ring = 151,
                  N_Brood = 54,
                  Brood = as.numeric((cortbowl$Brood))
)


model_string <- '
data{
int N;
real totCortlog[N];
real ImplantP[N];
real days20[N];
real daysbefore[N];
real Age_z[N];
int Brood[N];	
int Ring[N];
real ImplantP_X_days20[N];
real ImplantP_X_daysbefore[N];
int N_Brood;
int N_Ring;
}
parameters{
real Intercept;
real beta_ImplantP;
real beta_days20;
real beta_daysbefore;
real beta_Age_z;
real beta_ImplantP_X_days20;
real beta_ImplantP_X_daysbefore;
real<lower=0> sigma;
real vary_Ring[N_Ring];
real vary_Brood[N_Brood];
real<lower=0> sigma_Ring;
real<lower=0> sigma_Brood;
}
model{
real vary[N];
real glm[N];
// Priors
Intercept ~ normal( 0 , 100 );
beta_ImplantP ~ normal( 0 , 100 );
beta_days20 ~ normal( 0 , 100 );
beta_daysbefore ~ normal( 0 , 100 );
beta_Age_z ~ normal( 0 , 100 );
beta_ImplantP_X_days20 ~ normal( 0 , 100 );
beta_ImplantP_X_daysbefore ~ normal( 0 , 100 );
sigma_Ring ~ uniform( 0 , 100 );
sigma_Brood ~ uniform( 0 , 100 );    
sigma ~ uniform( 0 , 100 );
// Varying effects
for ( j in 1:N_Ring ) vary_Ring[j] ~ normal( 0 , sigma_Ring );
for ( j in 1:N_Brood ) vary_Brood[j] ~ normal( 0 , sigma_Brood );
// Fixed effects
for ( i in 1:N ) {
vary[i] <- vary_Ring[Ring[i]] + vary_Brood[Brood[i]];
glm[i] <- vary[i] + Intercept
+ beta_ImplantP * ImplantP[i]
+ beta_days20 * days20[i]
+ beta_daysbefore * daysbefore[i]
+ beta_Age_z * Age_z[i]
+ beta_ImplantP_X_days20 * ImplantP_X_days20[i]
+ beta_ImplantP_X_daysbefore * ImplantP_X_daysbefore[i];
}
totCortlog ~ normal( glm , sigma );
}
generated quantities{
real dev;
real vary[N];
real glm[N];
dev <- 0;
for ( i in 1:N ) {
vary[i] <- vary_Ring[Ring[i]];
glm[i] <- vary[i] + Intercept
+ beta_ImplantP * ImplantP[i]
+ beta_days20 * days20[i]
+ beta_daysbefore * daysbefore[i]
+ beta_Age_z * Age_z[i]
+ beta_ImplantP_X_days20 * ImplantP_X_days20[i]
+ beta_ImplantP_X_daysbefore * ImplantP_X_daysbefore[i];
dev <- dev + (-2) * normal_log( totCortlog[i] , glm[i] , sigma );
}
}
'

library(rstan)
stan_dso <- stan_model(model_code=model_string)

stanFit <- sampling(object = stan_dso, 
                    data = data_list, 
                    chains = 1, iter = 9000, 
                    warmup = 1000, thin = 1)

source('R\\stanmer2.R')
stanmer2(stanFit)

library(shinystan)
launch_shinystan(stanFit)
