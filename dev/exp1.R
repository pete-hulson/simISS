
# load/source stuff ----
library(tidyverse)
source(here::here('R', 'base_functions.R'))


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
# note: sample pop'n category sampling without replacement for example plots
sim_popn <- get_popn(d, pu, pc, pu_cv, replace_pc = FALSE, plot_name = 'exp1_popn')

## run sim loop ----
tictoc::tic()
rr_exp1 <- purrr::map(1:iters, ~sim_comp(su_num, sim_popn, su_samp, p_su_samp))
runtime_exp1 <- tictoc::toc()

 # unlist results
do.call(mapply, c(list, rr_exp1, SIMPLIFY = FALSE))$comp %>% 
  tidytable::map_df(., ~as.data.frame(.x), .id = "sim") -> res_sim

exp1_res <- list(sim_popn = sim_popn, res_sim = res_sim)

# save results
saveRDS(exp1_res,
        file = here::here('output', 'exp1_res.rds'))

## plot results ----

### simulation results ----
res_sim %>% 
  tidytable::rename(Wtd = samp_p_wtd, Unwtd = samp_p_unwtd) %>% 
  tidytable::pivot_longer(cols = c(Wtd, Unwtd)) %>% 
  tidytable::select(cat, 'Expansion type' = name, comp = value, selex_type) %>% 
  tidytable::bind_rows(sim_popn$p_true %>% 
                         tidytable::mutate(comp_type = 'True') %>% 
                         tidytable::rename(comp = p_true)) -> plot_dat
# add a random comp
res_sim %>% 
  tidytable::filter(sim == sample(1:iters, 1)) %>% 
  tidytable::select(-sim) %>% 
  tidytable::pivot_longer(cols = c(samp_p_wtd, samp_p_unwtd)) %>% 
  tidytable::mutate(shape = case_when(name == 'samp_p_wtd' ~ 24,
                                      name == 'samp_p_unwtd' ~ 25),
                    name = case_when(name == 'samp_p_wtd' ~ 'Random Wtd',
                                     name == 'samp_p_unwtd' ~ 'Random Unwtd')) %>% 
  tidytable::select(cat, selex_type, 'Expansion type' = name, comp = value, shape) -> rand_sim

plot <- ggplot(data = plot_dat, aes(x = as.factor(cat), y = comp)) +
  geom_bar(data = plot_dat %>% tidytable::filter(comp_type %in% c('True')), stat = 'identity', 
           aes(col = comp_type), alpha = 0) +
  geom_boxplot(data = plot_dat %>% tidytable::filter(`Expansion type` %in% c('Wtd', 'Unwtd')), 
               aes(fill = `Expansion type`), position = position_dodge(0.8), width = 0.5, alpha = 0.7) +
  geom_line(data = rand_sim, linewidth = 0.25, aes(x = cat, y = comp, col = `Expansion type`)) +
  # geom_point(data = rand_sim, fill = 'white', shape = rand_sim$shape, aes(x = cat, y = comp)) +
  facet_wrap(~selex_type, ncol = 1, scales = 'free_y') +
  scale_color_manual(values = c(scico::scico(3, palette = 'roma')[1], scico::scico(3, palette = 'roma')[2], scico::scico(3, palette = 'roma')[3]))  +
  scale_fill_manual(values = c(scico::scico(3, palette = 'roma')[1], scico::scico(3, palette = 'roma')[2])) +
  theme_bw() +
  xlab('Category') +
  ylab('Composition') +
  guides(col = 'none')

ggsave(filename = "exp1_sim.png",
       plot = plot,
       path = here::here("figs"),
       width = 6.5,
       height = 5,
       units = "in")

