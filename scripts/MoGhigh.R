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
iterations = 8000

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
         theta = c(.5),
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
mog <- stan_model(file = "stanin/MoGhigh.stan")

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
              pars = c("RE", "log_theta"), 
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
        file = "stanout/MoGhigh.rda",
        compress = "xz")

# Traceplots
param <- c("beta",  "delta", "theta", "sigma") 
traceplot(m, param, inc_warmup = FALSE )
#include = TRUE, unconstrain = FALSE, inc_warmup = FALSE, window = NULL
summary(print(m, pars = param, probs = c(.025,.975)))

as.data.frame(m, pars = param) %>%
  gather(Parameter, value) %>%
  #  mutate(value = exp(value)) %>%
  dplyr::group_by(Parameter) %>%
  dplyr::summarise_all(funs(
    M = mean(.),
    Lower = HPDI(., prob = .95)[1],
    Upper = HPDI(., prob = .95)[2])) %>% 
  print(n=25)


as.data.frame(m, pars = c("beta", "beta2")) %>%
  gather(Parameter, value) %>%
  mutate(value = exp(value),
         Attachment = ifelse(Parameter == "beta", "correct", 
                             ifelse(Parameter == "beta2", "incorrect", NA))) %>%
  filter(value < 1500) %>%
  ggplot(aes(x=value, color = Attachment, fill=Attachment)) +
  geom_histogram(aes(y=..density..),  position = "identity", alpha = .1, bins = 60) +
  geom_density(aes(y=..density..), fill = "transparent") +
  ggtitle("Posterior distributions: MoG") +
  theme_minimal() +
  scale_color_hue("Attachment: ", c = 80, l=20) +
  scale_fill_hue("Attachment: ", c = 80, l=20) +
#  scale_x_continuous(breaks = seq(200,500,50)) +
  xlab("Posterior onset latency (in ms)") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "top",
        legend.justification = "left",
        legend.title = element_text(face = "bold"),
        legend.key.width = unit(.75, "cm"),
        plot.title = element_text(size = 14, face = "bold"))

as.data.frame(m, pars = c("theta")) %>%
  as.tibble() %>%
  gather(Parameter, value) %>%
  separate(Parameter, into =c("drop", "attachment"),  extra = "merge") %>%
  mutate(attachment = ifelse(attachment == "1]","correct","incorrect" ) ) %>%
  select(-drop) %>%
  dplyr::rename(theta = value) %>%
  ggplot(aes(x=theta, color = attachment, fill = attachment)) +
  geom_histogram(aes(y=..density..), position = "identity", 
                 alpha = .4, bins = 40) +
  geom_density(aes(y=..density..), fill = "transparent") +
  ggtitle("Posterior distributions: MoG") +
  theme_minimal() +
  scale_color_hue("Attachment: ", c = 80, l=20) +
  scale_fill_hue("Attachment: ", c = 80, l=20) +
  xlab(bquote(theta~"(probability of incorrect attachment)")) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "top",
        legend.justification = "left",
        legend.title = element_text(face = "bold"),
        legend.key.width = unit(.75, "cm"),
        plot.title = element_text(size = 14, face = "bold"))

#param <- c("beta", "theta", "sigma", "sigmap_e", "sigma_e") 


  
