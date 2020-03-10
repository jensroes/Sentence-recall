library(xtable)
library(tidyverse)

d <- read_csv("loos/modelcomparisons.csv") %>%
  mutate(model = c("M6", "M7", "M5", "M4", "M3", "M2", "M1")) %>%
  select(model, elpd_diff, se_diff) 

print(xtable(d, digits = 1), include.rownames = F)
