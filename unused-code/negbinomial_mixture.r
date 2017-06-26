full_ngb_mixed_model = """
data {
  int<lower=0> N; #samples
  real<lower=0> mean1; # mean counts inferred
  real<lower=0> tau1;
  int<lower=0> y[N];
}

parameters{
  real<lower=0> mean2;
  real<lower=0>  tau2;
  real<lower=0> alpha;
  real<lower=0> beta;
  real<lower=0,upper=1>  pi;
}


model {
  pi ~ Beta(alpha, beta);
  for( i in 1:N){
    target += log_mix(pi,neg_binomial_2_lpmf(y[i] | mean1,tau1),neg_binomial_2_lpmf(y[i] | mean2,tau2));
  }
}
"""
library("rstan") # observe startup messages
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


data = list(N=N,mean1=mean,tau1=dispersion,y=mixed)

#compiles the model
stan_model <- stan_model(model_code =full_ngb_mixed_model)

#fit via HMC
fit <- sampling(stan_model, data=data, par=c("pi","mean2","tau2"), iter = 1000, chains = 3, thin=1)

summary(fit)
traceplot(fit)


library("bayesplot")
library("ggplot2")

posterior <-as.matrix(extract(fit)$pi)
colnames(posterior)<-paste("pi",1:N)
plot_title <- ggtitle("Posterior distributions with medians and 80% intervals")
mcmc_areas(posterior, point_est="mean",prob = 0.8) + plot_title
