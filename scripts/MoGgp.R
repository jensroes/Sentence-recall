# Load packages
#source("functions/functions.R")
library(tidyverse)
library(rstan)
library(rethinking)
library(magrittr)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Sampling parameters
set.seed(125)
n_chain = n_core = 3 # number of cores/chains
iterations = 20000

d <- read_csv("dfs/attachment.csv") %>%
  filter(
    onset > 200,
    ld == 0,
    n_edits < 11) %>%
  mutate(subj = as.integer(factor(subj)))

d$attachment <- d %$% factor(np, levels = c("Glob",  "NP", "VP"), ordered = TRUE)

# Data as list for stan input
attachment <- d %$% as.numeric(plyr::mapvalues(attachment, from = levels(attachment), to= 1:length(unique(attachment))))
dat <-  within( list(), {
                  N <- nrow(d)
                  y <- d$onset # DV
                  attachment <- attachment
                  K <- max(attachment)
                  subj <- d$subj
                  S <- max(d$subj)
                  item <- as.integer(factor(d$item))
                  I <- length(unique(d$item))
                } );str(dat)


# Initialise start values
start <- function(chain_id = 1){
    list(beta_raw = 0, 
         mu_beta = 6,
         sigma_beta = .01,
         sigma_raw = c(.1,.2), 
         mu_sigma = 0,
         sigma_sigma = .1,
         delta = 0.1,
         theta = c(.5,.5),
         u = rep(0, max(d$subj)),
         w = rep(0, length(unique(d$item))),
         sigma_u = .1,
         sigma_w = .1
    )
  }

start_ll <- lapply(1:n_chain, function(id) start(chain_id = id) )

# --------------
# Stan models ##
# --------------
#---- 
# Load modelp
mog <- stan_model(file = "stanin/MoGgp.stan")

# Check model
m <- sampling(mog, chain = 1, iter = 1, data = dat) 

# Fit model
m <- sampling(mog, 
              data = dat,
              init = start_ll,
              iter = iterations,
              warmup = iterations/2,
              chains = n_chain, 
              cores = n_core, 
              include = FALSE,
              pars = c("RE", "log_theta"), #"beta_raw", "mu_beta", "sigma_beta" , "z_w", "z_u", 
          #    "L_u", "L_w"
              save_warmup = FALSE, # Don't save the warmup
              refresh = 250,
              thin = 1,
              seed = 81,
              control = list(max_treedepth = 16,
                             adapt_delta = 0.99,
                             stepsize = 2)
)

# Save posterior
saveRDS(m, 
        file = "stanout/MoGgp.rda",
        compress = "xz")

# Traceplots
param <- c("beta",  "delta", "theta", "sigma") 
traceplot(m, param, inc_warmup = FALSE )
#include = TRUE, unconstrain = FALSE, inc_warmup = FALSE, window = NULL
summary(print(m, pars = param, probs = c(.025,.975)))

