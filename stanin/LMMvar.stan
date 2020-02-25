/*
  LMM with individuals means for each component
  as well as random intercepts for subject and items
  Smaller variance component for Glob attachment
*/ 

data {
	int<lower=1> N;                    // Number of observations
	real y[N];  		            //outcome
	int<lower=1> K; // N of conditions
	int<lower=1, upper=K> attachment[N];  //predictor

	int<lower=1> S;                  //number of subjects
	int<lower=1, upper=S> subj[N];   //subject id
	int<lower=1> I;                  //number of items
	int<lower=1, upper=I> item[N];   //item id
}

parameters {
 // For random effects
	vector[S] u; //subject intercepts
	real<lower=0> sigma_u;//subj sd

	vector[I] w; //items intercepts
	real<lower=0> sigma_w;//items sd

  
  // Non-centre parameters
  real mu_beta;
  real<lower = 0> sigma_beta;
  vector[K] beta_raw;
  
  real mu_sigma;
  real<lower=0> sigma_sigma;
  ordered[2] sigma_raw;
}

transformed parameters{
  vector[K] beta;
  vector[N] mu;
	positive_ordered[2] sigma;

  for(n in 1:N){
    beta[attachment[n]] = mu_beta + sigma_beta * beta_raw[attachment[n]];
    mu[n] = beta[attachment[n]] + u[subj[n]] + w[item[n]];
    if(attachment[n] == 1){
      sigma[1] =  mu_sigma + sigma_sigma * sigma_raw[1];
    }
    else{
      sigma[2] =  mu_sigma + sigma_sigma * sigma_raw[2];
    }
  }
}

model {

  // Priors
  mu_beta ~ cauchy(5, 2.5);
  sigma_beta ~ normal(0, 2.5);
  beta_raw ~ normal(0, 1);
  
  mu_sigma ~ cauchy(0, 2.5);
  sigma_sigma ~ cauchy(0, 2.5);
	sigma_raw ~ normal(0, 1);

	// REs priors
	sigma_u ~ normal(0,2.5);
	u ~ normal(0, sigma_u); //subj random effects

	sigma_w ~ normal(0,2.5);
	w ~ normal(0, sigma_w); //items random effects

  // likelihood
  for(n in 1:N){
    if(attachment[n] == 1){
      y[n] ~ lognormal(mu[n], sigma[1]);
    }
    else{
      y[n] ~ lognormal(mu[n], sigma[2]);
    }
  }
}

generated quantities{
	real log_lik[N];
	vector[N] y_tilde;
	for (n in 1:N){
    if(attachment[n] == 1){
  		log_lik[n] = lognormal_lpdf(y[n] | mu[n], sigma[1]);
    	y_tilde[n] = lognormal_rng(mu[n], sigma[1]);
    }
    else{
  		log_lik[n] = lognormal_lpdf(y[n] | mu[n], sigma[2]);
    	y_tilde[n] = lognormal_rng(mu[n], sigma[2]);
    }
	}
}
