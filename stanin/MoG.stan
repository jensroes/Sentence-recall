/*
 This model is adapted from Vasishth et al. 2007 
 Means and variance same across condition but the mixing proportion is allowed to vary
 Variance for longer latencies is larger
*/


data {
	int<lower=1> N;                    // Number of observations
	real y[N];  		            //outcome
	int<lower=1> K;             // N of conditions 
	int<lower=1, upper=K> attachment[N];  //predictor
	
	int<lower=1> S;                  //number of subjects
	int<lower=1, upper=S> subj[N];   //subject id

	int<lower=1> I;                  //number of items
	int<lower=1, upper=I> item[N];   //items id
	
}


parameters {
	real<lower=0> delta;			// distribution + extra component
	real<lower = 0.001, upper=.999> theta[K]; //probability of extreme values
	
  // For random effects
  vector[S] u; //subject intercepts
	real<lower=0> sigma_u;//subj sd

	vector[I] w; //items intercepts
	real<lower=0> sigma_w;//items sd

  // Non-centre parameters
  real mu_beta;
  real<lower = 0> sigma_beta;
  real beta_raw;

  real mu_sigma;
  real<lower=0> sigma_sigma;
  ordered[2] sigma_raw;

}

transformed parameters{
  matrix[N,2] mu;
	matrix[2,3] log_theta;
	real beta;
	positive_ordered[2] sigma;
  
  sigma = mu_sigma + sigma_sigma * sigma_raw;

  for(i in 1:K){
    log_theta[1,i] = log(theta[i]);
    log_theta[2,i] = log1m(theta[i]);
  }
  
  for(n in 1:N){
      beta = mu_beta + sigma_beta * beta_raw;
      mu[n,1] = beta + u[subj[n]] + w[item[n]];
      mu[n,2] = beta + delta + u[subj[n]] + w[item[n]];
	}
}

model {
  // Priors
  theta ~ beta(2, 2);

  mu_beta ~ normal(6, 4);
  sigma_beta ~ normal(0, 1);
  beta_raw ~ normal(0, 1);

	delta ~ normal(0, 1);

  mu_sigma ~ normal(0, 2.5);
  sigma_sigma ~ normal(0, 1);
  sigma_raw ~ normal(0, 1);
    
	// REs priors
	sigma_u ~ normal(0,2.5);
	u ~ normal(0, sigma_u); //subj random effects

	sigma_w ~ normal(0,2.5);
	w ~ normal(0, sigma_w); //items random effects
	
  
  // Likelihood	
  for (n in 1:N){
    target += log_sum_exp(
 		     log_theta[1,attachment[n]] + lognormal_lpdf(y[n] | mu[n,1], sigma[1]), 
 		     log_theta[2,attachment[n]] + lognormal_lpdf(y[n] | mu[n,2], sigma[2]) 
    ); 
  }
}

generated quantities{
	real log_lik[N];
	vector[N] y_tilde;
  real<lower=0,upper=1> theta_tilde[K]; 

  // likelihood: 
  for(n in 1:N){ 
	    log_lik[n] = log_sum_exp(
 		      log_theta[1,attachment[n]] + lognormal_lpdf(y[n] | mu[n,1], sigma[1]), 
 		      log_theta[2,attachment[n]] + lognormal_lpdf(y[n] | mu[n,2], sigma[2])
	    );
       theta_tilde[attachment[n]] = bernoulli_rng(theta[attachment[n]]); 
          if(theta_tilde[attachment[n]]) { 
              y_tilde[n] = lognormal_rng(mu[n,1], sigma[1]);
          }
          else{
              y_tilde[n] = lognormal_rng(mu[n,2], sigma[2]);
          }
    }
}



