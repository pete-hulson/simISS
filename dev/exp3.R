
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
# run simulation
rr_bs <- purrr::map(1:bs_iters, ~bs_sim(d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters))
# end timer
runtime <- tictoc::toc()

(runtime$toc - runtime$tic) / (60 * bs_iters) * 1000 / 60




# unlist results
do.call(mapply, c(list, rr_bs, SIMPLIFY = FALSE))$rss_se %>% 
  tidytable::map_df(., ~as.data.frame(.x), .id = "sim") -> rss_se
do.call(mapply, c(list, rr_bs, SIMPLIFY = FALSE))$iss_bs %>% 
  tidytable::map_df(., ~as.data.frame(.x), .id = "sim") -> iss_bs

res_bs <- list(rss_se = rss_se, iss_bs = iss_bs)

# save results
saveRDS(res_bs,
        file = here::here('output', paste0('res_bs.rds')))

# plot results
rss_se %>% 
  tidytable::pivot_longer(cols = c(rss_wtd, rss_unwtd)) %>% 
  tidytable::rename(sim_rss = value) %>% 
  tidytable::mutate(weighting = case_when(name == 'rss_wtd' ~ TRUE, .default = FALSE)) %>% 
  tidytable::select(-name) %>% 
  tidytable::left_join(iss_bs %>% 
                         tidytable::pivot_longer(cols = c(bs_iss_wtd, bs_iss_unwtd)) %>% 
                         tidytable::rename(bs_iss = value) %>% 
                         tidytable::mutate(weighting = case_when(name == 'bs_iss_wtd' ~ TRUE, .default = FALSE)) %>% 
                         tidytable::select(-name)) -> plot_data


bs_plot <- ggplot(data = plot_data, aes(x = sim_rss, y = bs_iss, col = weighting)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  scico::scale_color_scico_d(palette = 'roma') +
  xlab("RSS for sampling event of generated population") +
  ylab("ISS for bootstrap of sampling event") +
  theme_bw()

ggsave(filename = "bs_test.png",
       plot = bs_plot,
       path = here::here("figs"),
       width = 6.5,
       height = 5,
       units = "in")








