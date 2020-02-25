/*
 This model is adapted from Vasishth et al. 2007 
 This model is testing the prediction that high attachment increases the onset latency.
 For global ambiguities this is capture by the mixing proportion theta that requires longer latencies.
*/


data {
	int<lower=1> N;                    // Number of observations
  	real y[N];  		            //outcome
	int<lower=1> K;             // N of copy task components 
	int<lower=1, upper=K> attachment[N];  //predictor
	
	int<lower=1> S;                  //number of subjects
	int<lower=1, upper=S> subj[N];   //subject id

	int<lower=1> I;                  //number of items
	int<lower=1, upper=I> item[N];   //items id
	
}


parameters {
//  real beta;
	real<lower=0> delta;			// distribution + extra component
	real<lower = 0.001, upper=.999> theta; //probability of extreme values
	
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
	vector[2] log_theta;
  vector[N] RE;
	real beta = mu_beta + sigma_beta * beta_raw;
	positive_ordered[2] sigma;

  sigma = mu_sigma + sigma_sigma * sigma_raw;

  log_theta[1] =  log(theta);
  log_theta[2] =  log1m(theta);

  for(n in 1:N){
      RE[n] = u[subj[n]] + w[item[n]];
	}
}

model {
  // Priors
  mu_beta ~ normal(6, 4);
  sigma_beta ~ normal(0, 1);
  beta_raw ~ normal(0, 1);
  
  //beta ~ normal(6, 1);
	delta ~ normal(0, 1);
  theta ~ beta(1, 1);

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
    if(attachment[n] == 1){ // Global ambiguous
      target += log_sum_exp(
   		     log_theta[2] + lognormal_lpdf(y[n] | beta + delta + RE[n], sigma[2]), 
   		     log_theta[1] + lognormal_lpdf(y[n] | beta + RE[n], sigma[1]) 
      );
    }
    if(attachment[n] == 2){
      target += lognormal_lpdf(y[n] | beta + RE[n], sigma[1]);  
    }
    if(attachment[n] == 3){
      target += lognormal_lpdf(y[n] | beta + delta + RE[n], sigma[2]);  
    }
  }
}

generated quantities{
	real log_lik[N];
	vector[N] y_tilde;
  real<lower=0,upper=1> theta_attachment; 
  
  // likelihood: 
  for(n in 1:N){ 
    if(attachment[n] == 1){
	    log_lik[n] = log_sum_exp(
 		      log_theta[2] + lognormal_lpdf(y[n] | beta + delta + RE[n], sigma[2]), 
 		      log_theta[1] + lognormal_lpdf(y[n] | beta + RE[n], sigma[1]));
      theta_attachment = bernoulli_rng(theta); 
      if(theta_attachment) { 
          y_tilde[n] = lognormal_rng(beta + RE[n], sigma[1]);
      }
      else{
        y_tilde[n] = lognormal_rng(beta + delta + RE[n], sigma[2]);
      }
    }
    if(attachment[n] == 2){
      log_lik[n] = lognormal_lpdf(y[n] | beta + RE[n], sigma[1]);
      y_tilde[n] = lognormal_rng(beta + RE[n], sigma[1]);
    }
    if(attachment[n] == 3){
      log_lik[n] = lognormal_lpdf(y[n] | beta + delta + RE[n], sigma[2]);
      y_tilde[n] = lognormal_rng(beta + delta + RE[n], sigma[2]);
    }
  }
}
