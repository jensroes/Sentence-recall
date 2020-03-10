library(tidyverse)
library(ggthemes)
library(ggeffects)
library(tidybayes)

m <- readRDS("stanout/LMM.rda")

as.data.frame(m, pars = c("beta") ) %>%
  as_tibble() %>%
  rename(glob = `beta[1]`,
         low = `beta[2]`,
         high = `beta[3]`) %>%
  gather(Parameter, value) %>%
  mutate(value = exp(value)) %>%
  mutate(Parameter = recode_factor(Parameter, 
                                   glob = "Global ambiguity",
                                   low = "Low attachment",
                                   high = "High attachment",
                                   .ordered = T)) -> ps

ps %>%
  pivot_wider(names_from = Parameter, 
              values_from = value, 
              values_fn = list(value = list)) %>%
  unnest() %>%
  mutate(`Temporary ambiguity` = (`Low attachment` + `High attachment`)/2) %>%
  transmute(diff = `Temporary ambiguity`-`Global ambiguity`) %>%
  summarise(M = mean(diff), 
            lo = quantile(diff, probs = .025),
            up = quantile(diff, probs = .975)) %>%
  pivot_longer(cols = M:up, names_to = "Parameter", values_to = "value") %>%
  pull(value) %>% round(0) -> diff; diff


ps %>%
  ggplot(aes(y = .25, x = value, fill = Parameter)) +
  stat_halfeyeh(alpha = .5) +
  scale_fill_ptol("") +
  scale_x_continuous(breaks = seq(100, 1100, 50), limits = c(390, 570)) +
  labs(x = "Onset latency [in msecs]",
       y = "",
       title = "Experiment 2: M2", #subtitle = bquote("Posterior onset latency"),
       subtitle = bquote("Difference:"~Delta*hat(beta)==.(diff[1])*" msecs; 95% PI["*.(diff[2])*","~.(diff[3])*"]")) +
  theme_bw(base_size = 12) +
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), 
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(.1,0.2, 0, -.5), "cm"),
        axis.text.y = element_blank(),
        legend.position = "bottom",
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.justification = c(1,0),
        legend.margin = margin(unit(c(-5,0,0,0), "cm")),
        axis.ticks.y = element_blank()) 

ggsave(filename = "/home/jens/Documents/Latexsource/posters/cuny2020/gfx/Exp2LMMv2.pdf",  width = 5.75, height = 3.75,  device = cairo_pdf, bg = "transparent")
  


