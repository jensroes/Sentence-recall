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
                                   .ordered = T)) %>%
  ggplot(aes(y = .25, x = value, fill = Parameter)) +
  stat_halfeyeh(alpha = .5) +
  scale_fill_ptol("") +
  labs(x = "Onset latency [in msecs]",
       y = "") +
  theme_bw(base_size = 8) +
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        legend.position = "bottom",
        axis.ticks.y = element_blank()) -> beta_plot;beta_plot

ggsave(filename = "gfx/Exp3LMM.png", plot = beta_plot, 
       device = png(),
       dpi = 320,
       width = 6, height = 2.5, bg = "transparent")



