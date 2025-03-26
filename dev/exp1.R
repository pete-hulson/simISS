
# load/source stuff ----
library(tidyverse)
source(here::here('R', 'functions.R'))


# set up experiment parameters ----

## simulation/bootstrap parameters ----

# number of simulation/bootstrap iterations
iters <- 1000

## sampling parameters ----

# number of sampling units (e.g., hauls)
su_num <- 250

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
pu_cv <- 0.25

# pop'n exponential decay (e.g., lnM, tied to inverse of number of pop'n categories so that pop'n = 0.01 at largest category)
d <- log(0.01) / (1 - pc)


# experiment 1: ----
# do we get unbiased estimates of pop'n composition sampling different schools?
# what is effect of expansion complexity?

## get simulated pop'n ----
sim_popn <- get_popn(d, pu, pc, pu_cv)

## run sim loop ----
tictoc::tic()
rr_sim <- purrr::map(1:iters, ~sim_comp(su_num, sim_popn, su_samp, p_su_samp))
runtime_sim <- tictoc::toc()

 # unlist results
do.call(mapply, c(list, rr_sim, SIMPLIFY = FALSE))$comp %>% 
  tidytable::map_df(., ~as.data.frame(.x), .id = "sim") -> res_sim

examp_sim <- list(sim_popn = sim_popn, res_sim = res_sim)

# save results
saveRDS(examp_sim,
        file = here::here('output', 'examp_sim.rds'))

## plot results ----

### simulation results ----
res_sim %>% 
  tidytable::pivot_longer(cols = c(samp_p_wtd, samp_p_unwtd)) %>% 
  tidytable::select(cat, comp_type = name, comp = value, selex_type) %>% 
  tidytable::bind_rows(sim_popn$p_true %>% 
                         tidytable::mutate(comp_type = 'true') %>% 
                         tidytable::rename(comp = p_true)) -> plot_dat
# add a random comp
res_sim %>% 
  tidytable::filter(sim == sample(1:iters, 1)) %>% 
  tidytable::select(-sim) %>% 
  tidytable::rename(rand_wtd = samp_p_wtd,
                    rand_unwtd = samp_p_unwtd) %>% 
  tidytable::pivot_longer(cols = c(rand_wtd, rand_unwtd)) %>% 
  tidytable::rename(comp_type = name, comp = value) -> rand_sim

sim_plot <- ggplot(data = plot_dat, aes(x = as.factor(cat), y = comp)) +
  geom_bar(data = plot_dat %>% tidytable::filter(comp_type %in% c('true')), stat = 'identity', 
           aes(fill = comp_type), alpha = 0.5) +
  geom_boxplot(data = plot_dat %>% tidytable::filter(comp_type %in% c('samp_p_wtd', 'samp_p_unwtd')), 
               aes(fill = comp_type)) +
  geom_line(data = rand_sim, aes(x = cat, y = comp, col = comp_type), linetype = 'dashed') +
  geom_point(data = rand_sim, size = 1.5, col = 'black', aes(x = cat, y = comp, shape = comp_type)) +
  facet_wrap(~selex_type, ncol = 1) +
  scico::scale_fill_scico_d(palette = 'roma')  +
  scale_color_manual(values = c(scico::scico(3, palette = 'roma')[1], scico::scico(3, palette = 'roma')[2])) +
  theme_bw() +
  xlab('category')

ggsave(filename = "ex_sim_popn.png",
       plot = sim_plot,
       path = here::here("figs"),
       width = 6.5,
       height = 5,
       units = "in")

## notes: ----
# so long as the schools are represented across the haul samples in proportion to their relative abundance, 
# the mean across iterations in unbiased, whether weighted or not, so,
# you don't need to weight by haul's catch even if the number of samples within a haul aren't proportional to catch
# however, for any given iteration, the composition resulting from weighted or unweighted comps can be different.
# the uncertainty around weighted comps is larger than unweighted comps
