
# load/source stuff ----
library(tidyverse)
source(here::here('R', 'functions.R'))


# set up experiment parameters ----

## simulation/bootstrap parameters ----

# number of simulation/bootstrap iterations
iters <- 1000

## sampling parameters ----

# number of sampling units (e.g., hauls)
su_num <- 500

# number of samples per sampling unit (e.g., ages/lengths)
# vector to test differences in number of samples across sampling unit
su_samp <- c(20, 10)

# vector pf probabilities to select options of number of samples in sampling unit 
p_su_samp <- c(0.9, 0.1)

## population unit parameters ----

# number of population units sampled (e.g., number of schools/year classes)
pu <- 5

# number of population categories (e.g., ages or length bins)
pc <- 15

# CV around mean within pop'n unit (e.g., CV in mean age of school)
pu_cv <- 0.2

# pop'n exponential decay (e.g., lnM, tied to inverse of number of pop'n categories so that pop'n = 0.01 at largest category)
d <- log(0.01) / (1 - pc)


# experiment 1: ----
# do we get unbiased estimates of pop'n composition sampling different schools?
# what is effect of expansion complexity?

## get simulated pop'n ----
sim_popn <- get_popn(d, pu, pc, pu_cv)

# need to run to get example desired for sim_popn
sim_popn_multi <- sim_popn
sim_popn_uni <- sim_popn
sim_popn_rec <- sim_popn

# plot example
sim_popn_multi$p_true %>% 
  tidytable::mutate(popn_strctr = sim_popn_multi$popn_strctr) %>% 
  tidytable::bind_rows(sim_popn_uni$p_true %>% 
                         tidytable::mutate(popn_strctr = sim_popn_uni$popn_strctr)) %>% 
  tidytable::bind_rows(sim_popn_rec$p_true %>% 
                         tidytable::mutate(popn_strctr = sim_popn_rec$popn_strctr)) %>% 
  tidytable::mutate(popn_strctr = factor(popn_strctr, levels = c('recruitment pulse', 'multimodal', 'unimodal'))) -> plot_dat

ex_plot <- ggplot(data = plot_dat, aes(x = as.factor(cat), y = p_true, fill = popn_strctr)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~popn_strctr, ncol = 1, scales = 'free_y') +
  scico::scale_fill_scico_d(palette = 'roma')  +
  theme_bw() +
  xlab('category') +
  ylab('comp')

ggsave(filename = "ex_popn_strctr.png",
       plot = ex_plot,
       path = here::here("figs"),
       width = 6.5,
       height = 5,
       units = "in")


