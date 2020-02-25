/* 
 LMM with individuals means for each attachment type
 as well as random intercepts for subject and items
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
	real<lower=0> sigma;		// residual sd
	
    // For random effects
	vector[S] u; //subject intercepts
	real<lower=0> sigma_u;//subj sd

	vector[I] w; //items intercepts
	real<lower=0> sigma_w;//items sd

	
  // Non-centre parameters
  real mu_beta;
  real<lower=0> sigma_beta;
  vector[K] beta_raw;

}

transformed parameters{
  vector[K] beta;
  vector[N] mu;

  for(n in 1:N){
    beta[attachment[n]] = mu_beta + sigma_beta * beta_raw[attachment[n]];
    mu[n] = beta[attachment[n]] + u[subj[n]] + w[item[n]];
  }
}

model {
  // Priors
  mu_beta ~ normal(6, 4);
  sigma_beta ~ normal(0, 1);
  beta_raw ~ normal(0, 1);
	sigma ~ cauchy(0, 2.5);
	
	// REs priors
	sigma_u ~ normal(0,2.5);
	u ~ normal(0, sigma_u); //subj random effects

	sigma_w ~ normal(0,2.5);
	w ~ normal(0, sigma_w); //items random effects

  // likelihood
  y ~ lognormal(mu, sigma);
}

generated quantities{
	real log_lik[N];
	vector[N] y_tilde;
	for (n in 1:N){
		log_lik[n] = lognormal_lpdf(y[n] | mu[n], sigma);
		y_tilde[n] = lognormal_rng(mu[n], sigma);
	}
}
