## ---- echo=FALSE, results='asis' or jump to load(workspace) ----------------------------------
load('R\\datac.RData')
library(knitr)
my_model <- '
data{
    int N;
    int outcome[N];
    real treatment[N];
    real condition[N];
}

parameters{
    real Intercept;
    real beta_treatment;
    real beta_condition;
}

model{
    real glm[N];
    // Priors
    Intercept ~ normal( 0 , 100 );
    beta_treatment ~ normal( 0 , 100 );
    beta_condition ~ normal( 0 , 100 );
    // Fixed effects
    for ( i in 1:N ) {
        glm[i] <- Intercept
                + beta_treatment * treatment[i]
                + beta_condition * condition[i];
        glm[i] <- inv_logit( glm[i] );
    }
    outcome ~ binomial(1 , glm );
}

generated quantities{
    real dev;
    real glm_rep[N];
    int y_rep[N]; // A new vector!
    dev <- 0;
    for ( i in 1:N ) {
        glm_rep[i] <- Intercept
                + beta_treatment * treatment[i]
                + beta_condition * condition[i];
        y_rep[i] <- bernoulli_rng(inv_logit(Intercept
                + beta_treatment * treatment[i]
                + beta_condition * condition[i])); // replication!
        dev <- dev + (-2) * binomial_log( outcome[i] , 1 , inv_logit(glm_rep[i]) );
    }
}
' 

datalist <- list(outcome = datac$outcome,
              treatment = datac$treatment,
              condition = datac$condition,
              N = nrow(datac))
mymod<- stan( model_code = my_model,
               data = datalist)

## ---- echo=FALSE---------------------------------------------------------
load('output\\workspace.RData')

samples <- extract(mymod, pars='y_rep')


library(shinystan)
launch_shinystan(mymod)


