
# load/source stuff ----
library(tidyverse)
library(future)
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
pu <- 25

# number of population categories (e.g., ages or length bins)
pc <- 15

# CV around mean within pop'n unit (e.g., CV in mean age of school)
pu_cv <- 0.25

# pop'n exponential decay (e.g., lnM, tied to inverse of number of pop'n categories so that pop'n = 0.01 at largest category)
d <- log(0.01) / (1 - pc)


# experiment 3: can a given realization of a sample give the 'true' iss back? ----

# number of bootstrap replicates
bs_iters <- 10


#start timer
tictoc::tic()
## run simulation ----
rr_bs <- purrr::map(1:bs_iters, ~bs_sim(d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters))
# end timer
runtime <- tictoc::toc()

(runtime$toc - runtime$tic) / (60 * bs_iters) * 1000 / 60

## unlist & save results ----
do.call(mapply, c(list, rr_bs, SIMPLIFY = FALSE))$rss_se %>% 
  tidytable::map_df(., ~as.data.frame(.x), .id = "sim") -> rss_se
do.call(mapply, c(list, rr_bs, SIMPLIFY = FALSE))$iss_bs %>% 
  tidytable::map_df(., ~as.data.frame(.x), .id = "sim") -> iss_bs
do.call(mapply, c(list, rr_bs, SIMPLIFY = FALSE))$popn_strctr %>% 
  tidytable::map_df(., ~as.data.frame(.x), .id = "sim") %>% 
  tidytable::rename(popn_strctr = .x) -> popn_strctr

res_bs <- list(rss_se = rss_se, iss_bs = iss_bs, popn_strctr = popn_strctr)

# save results
saveRDS(res_bs,
        file = here::here('output', paste0('res_bs.rds')))

## plot results ----

rss_se %>%  
  tidytable::left_join(popn_strctr) %>% 
  tidytable::summarise(iss_wtd = psych::harmonic.mean(rss_wtd, zero = FALSE),
                       iss_unwtd = psych::harmonic.mean(rss_unwtd, zero = FALSE),
                       .by = popn_strctr) %>% 
  tidytable::bind_rows(rss_se %>%  
                         tidytable::left_join(popn_strctr) %>% 
                         tidytable::mutate(popn_strctr = 'combined') %>% 
                         tidytable::summarise(iss_wtd = psych::harmonic.mean(rss_wtd, zero = FALSE),
                                              iss_unwtd = psych::harmonic.mean(rss_unwtd, zero = FALSE),
                                              .by = popn_strctr)) %>% 
  tidytable::pivot_longer(cols = c(iss_wtd, iss_unwtd)) %>% 
  tidytable::rename(type = name, iss = value) %>% 
  tidytable::mutate(popn_strctr = factor(popn_strctr, levels = c('recruitment pulse', 'multimodal', 'unimodal', 'combined'))) -> true_iss
  
iss_bs %>%  
  tidytable::left_join(popn_strctr) %>%
  tidytable::rename(iss_wtd = bs_iss_wtd, iss_unwtd = bs_iss_unwtd) %>% 
  tidytable::bind_rows(iss_bs %>%  
                         tidytable::left_join(popn_strctr) %>%
                         tidytable::mutate(popn_strctr = 'combined') %>% 
                         tidytable::rename(iss_wtd = bs_iss_wtd, iss_unwtd = bs_iss_unwtd)) %>% 
  tidytable::pivot_longer(cols = c(iss_wtd, iss_unwtd)) %>% 
  tidytable::rename(type = name, iss = value) %>% 
  tidytable::mutate(popn_strctr = factor(popn_strctr, levels = c('recruitment pulse', 'multimodal', 'unimodal', 'combined'))) -> bs_iss

# plot for all pop'n structures
bs_plot <- ggplot(data = bs_iss, aes(x = type, y = iss, fill = type)) + 
  geom_boxplot() +
  geom_point(data = true_iss, shape = 24, size = 3.777, fill = 'white', aes(x = type, y = iss)) +
  facet_wrap(~popn_strctr, nrow = 1) +
  scico::scale_color_scico_d(palette = 'roma') +
  scico::scale_fill_scico_d(palette = 'roma') +
  theme_bw()

ggsave(filename = "bs_test_popn.png",
       plot = bs_plot_popn,
       path = here::here("figs"),
       width = 6.5,
       height = 5,
       units = "in")

# plot for combined pop'n structures
bs_plot_comb <- ggplot(data = bs_iss %>% tidytable::filter(popn_strctr == 'combined'), aes(x = type, y = iss, fill = type)) + 
  geom_boxplot() +
  geom_point(data = true_iss %>% tidytable::filter(popn_strctr == 'combined'), shape = 24, size = 3.777, fill = 'white', aes(x = type, y = iss)) +
  scico::scale_color_scico_d(palette = 'roma') +
  scico::scale_fill_scico_d(palette = 'roma') +
  theme_bw()

ggsave(filename = "bs_test_comb.png",
       plot = bs_plot_comb,
       path = here::here("figs"),
       width = 6.5,
       height = 5,
       units = "in")
