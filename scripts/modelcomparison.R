rm(list =ls())
library(loo)
library(tidyverse)
library(rstan)

path <- "stanout/"
(files <- dir(path, pattern = ".rda"))
ms <- gsub(files, pattern = ".rda", replacement = "")

for(i in 1:length(files)){
  m <- readRDS(paste0(path, files[i] ))
  log_lik <- extract_log_lik(m, merge_chains = F) 
  r_eff <- relative_eff(exp(log_lik)) 
  assign(paste0("loo_",ms[i]), loo(log_lik, r_eff = r_eff, cores = 2))
  print(ms[i]); rm(list = "m")
}

(loos <- ls(pattern = "loo_"))
mc <- do.call(loo::loo_compare, lapply(loos, as.name))

tibble(loos) %>%
  mutate(Model = paste0("model", 1:n())) -> modelnames

mc %>% as.data.frame() %>%
  mutate(Model=row.names(.)) %>%
  select(Model, elpd_diff:se_elpd_loo) %>%
  remove_rownames() %>%
  left_join(modelnames) %>%
  mutate(Model = gsub("loo_", "", loos)) %>%
  select(-loos) -> mc;mc

file_out <- paste0("loos/modelcomparisons.csv")
write_csv(mc, file_out)

# Model weights
(loos_list <- do.call(list, lapply(loos, as.name)))
loo_model_weights(loos_list, method = "pseudobma")
loo_model_weights(loos_list, method = "stacking")
