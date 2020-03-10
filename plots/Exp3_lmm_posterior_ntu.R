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
                                   glob = "Global\nambiguity",
                                   low = "Low\nattachment",
                                   high = "High\nattachment",
                                   .ordered = T)) -> ps

ggplot(ps, aes(y = .1, x = value, fill = Parameter)) +
  stat_halfeyeh(alpha = .45) + 
  scale_fill_wsj(name ="") +
  labs(y = "Posterior probability", x = bquote("Onset latency [in msecs]"),
       title = "") +
  theme_bw(base_size = 12) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "top",
        legend.justification = "right",
        legend.key.width = unit(1, "cm"),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), 
        panel.grid.minor = element_blank()) 

filename <- "~/Documents/Latexsource/slides/slides-socsci-NTU-2020/gfx/recall_lmm.pdf"
ggsave(filename, width = 6.5, height = 5, device = cairo_pdf)



