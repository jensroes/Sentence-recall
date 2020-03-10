library(tidyverse)
library(ggthemes)
library(ggeffects)
library(rethinking)
library(gridExtra)
library(tidybayes)
library(magrittr)

m <- readRDS("stanout/MoG.rda")

as.data.frame(m, pars = c("beta", "delta") ) %>%
  transmute(delta = exp(beta + delta),
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
       title = "Experiment 2: M6", subtitle = bquote("a. Mixture components"), caption = "") +
  theme_bw(base_size = 12) +
  #  annotate("text", y = c(.15,.15), x = c(450, 500), label = c("beta", "beta + delta"), parse = T, size = 3) +
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), 
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        plot.caption = element_text(vjust = 6),
        plot.margin = unit(c(.25,0, -1, -.5), "cm"),
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
  labs(y = "", x = bquote("Probability of long latencies"),title = "", subtitle = bquote("b. Mixing proportion"~hat(theta)), caption = "") +
  theme_bw(base_size = 12) +
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), 
        plot.subtitle = element_text(vjust = -.1),
#        plot.caption = element_text(vjust = 2),
        plot.margin = unit(c(-.10,-.1, -.27, .5), "cm"), #   trbl
        panel.grid.minor = element_blank(),
        legend.position = "none") -> theta_plot; theta_plot

p <- grid.arrange(beta_plot, theta_plot, ncol = 2, widths =c(.4, .6))

ggsave(filename = "/home/jens/Documents/Latexsource/posters/cuny2020/gfx/Exp2v2.pdf", plot = p, width = 5.75, height = 3.95,  device = cairo_pdf, bg = "transparent")



  