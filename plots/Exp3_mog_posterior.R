library(tidyverse)
library(ggthemes)
library(ggeffects)
library(gridExtra)
library(tidybayes)

m <- readRDS("stanout/MoG.rda")

as.data.frame(m, pars = c("beta", "delta") ) %>%
  transmute(delta1 = exp(beta + delta),
#            delta2 = exp(beta) * exp(delta),
            beta = exp(beta)) %>%
  gather(Parameter, value) %>%
  #  mutate(Parameter = factor(Parameter, labels = c("Global ambiguous",  "Low attachment", "High attachment"))) %>%
  as_tibble() -> ps_beta;ps_beta

ps_theta <- as.data.frame(m, pars = c("theta") ) %>%
  gather(Parameter, value) %>%
  mutate(value = 1-value) %>%
  mutate(Parameter = recode_factor(Parameter, `theta[1]` =  "glob",
                                   `theta[2]` =  "low",
                                   `theta[3]` =  "high")) %>%
  mutate(Parameter = factor(Parameter, levels = c("low",  "glob",  "high"), ordered = TRUE)) %>%
  as_tibble();ps_theta

ps_beta %>%
  mutate(Parameter = factor(Parameter, levels = c("beta", "delta"), ordered = T)) %>%
  ggplot(aes(y = .0, x = value, fill = Parameter)) +
  stat_halfeyeh(alpha = .45) +
  scale_fill_manual("Mixture components", values = c("violetred4", "olivedrab4")) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = bquote(atop("Onset latency [in msecs]")),
       y = "",
       title = "a. Mixture components") +
  theme_bw(base_size = 8) +
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), 
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(vjust = -.25),
        axis.text.y = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") -> beta_plot;beta_plot

ps_theta %>%
  ggplot(aes(y = Parameter, x = value, fill = stat(abs(x) < .5))) +
  geom_vline(xintercept = c(.5), linetype = "dashed", alpha = .3) +
  stat_halfeyeh(alpha = .45) + 
  scale_fill_manual(values = c("olivedrab4", "olivedrab3")) +
  scale_x_continuous(limits = c(0,1)) +
  scale_y_discrete(labels = c('glob'= "Global\nambiguity",
                              'low' = "Low\nattachment",
                              'high' = "High\nattachment"), position = "right") +
  labs(y = "", x = bquote("Probability of long latencies"~hat(theta)),title = "b. Mixing proportion", caption = "") +
  theme_bw(base_size = 8) +
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), 
        panel.grid.minor = element_blank(),
        legend.position = "none") -> theta_plot; theta_plot

p <- grid.arrange(beta_plot, theta_plot, ncol = 2, widths =c(.35,.65))

ggsave(filename = "../manuscript/gfx/Exp3MoG.png", plot = p, 
       width = 6, height = 3,
       device = png(),
       dpi = 320,
       bg = "transparent")



