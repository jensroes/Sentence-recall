library(tidyverse)
library(ggthemes)
library(ggeffects)
library(gridExtra)
library(tidybayes)

m <- readRDS("stanout/MoG.rda")

ps_theta <- as.data.frame(m, pars = c("theta") ) %>%
  gather(Parameter, value) %>%
  mutate(value = 1-value) %>%
  mutate(Parameter = recode_factor(Parameter, `theta[1]` =  "glob",
                                   `theta[2]` =  "low",
                                   `theta[3]` =  "high")) %>%
  mutate(Parameter = factor(Parameter, levels = c("low",  "glob",  "high"), ordered = TRUE)) %>%
  as_tibble();ps_theta

ps_theta %>%
  mutate(Parameter = recode_factor(Parameter, low = "Low\nattachment",
                                   glob = "Global\nambiguity",
                                              high = "High\nattachment")) %>% 
  ggplot(aes(y = .01, x = value,fill = Parameter, alpha = stat(abs(x) > .5))) +
  geom_vline(xintercept = c(.5), linetype = "dashed", alpha = .3) +
  stat_halfeyeh() +
  scale_fill_wsj(name ="") +
  scale_alpha_discrete(range=c(.25,.65) ) +
  scale_x_continuous(limits = c(0,1)) +
  guides(alpha = F) +
  labs(y = "Posterior probability", x = bquote(hat(theta))) +
  theme_bw(base_size = 12) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "top",
        legend.justification = "right",
        legend.key.width = unit(1, "cm"),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), 
        panel.grid.minor = element_blank()) 

filename <- "~/Documents/Latexsource/slides/slides-socsci-NTU-2020/gfx/recall_mog.pdf"
ggsave(filename, width = 6.5, height = 5, device = cairo_pdf)



