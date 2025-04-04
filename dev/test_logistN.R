
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

# get simulated pop'n ----
# note: sample pop'n category sampling without replacement for example plots
sim_popn <- get_popn(d, pu, pc, pu_cv, plot = FALSE)

# run sim loop ----
rr_test <- purrr::map(1:iters, ~sim_comp(su_num, sim_popn, su_samp, p_su_samp))

# unlist results
res_sim <- do.call(mapply, c(list, rr_test, SIMPLIFY = FALSE))$comp %>% 
  tidytable::map_df(., ~as.data.frame(.x), .id = "sim") %>% 
  tidytable::pivot_longer(cols = c('samp_p_wtd', 'samp_p_unwtd')) %>% 
  tidytable::rename(comp_type = name, p_obs = value) %>% 
  tidytable::mutate(comp_type = case_when(comp_type == 'samp_p_wtd' ~ 'wtd',
                                          .default = 'unwtd'))
# compute input sample size
iss_sim <- res_sim %>% 
  tidytable::left_join(sim_popn$p_true) %>% 
  tidytable::summarise(rss = sum(p_true * (1- p_true)) / sum((p_obs - p_true) ^ 2),
                       .by = c(sim, selex_type, comp_type)) %>% 
  tidytable::summarise(iss = psych::harmonic.mean(rss, zero = FALSE),
                       .by = c(selex_type, comp_type))


# test logistic-normal estimation ----

# set up data for logistic-normal
data <- list(exp = sim_popn$p_true %>% 
               tidytable::select(-N_c), 
             obs = res_sim,
             iss = iss_sim)

# remove 0's
if(any(data$obs == 0) || any(data$exp == 0)) {
  # small constant
  eps <- 1e-6
  
  # expected
  data$exp <- data$exp %>% 
    # add constant in
    tidytable::mutate(p_true = case_when(p_true == 0 ~ eps,
                                         .default = p_true)) %>% 
    # renormalize
    tidytable::mutate(p_true = p_true / sum(p_true), 
                      .by = c(selex_type))
  
  # observed
  data$obs <- data$obs %>% 
    # add constant in
    tidytable::mutate(p_obs = case_when(p_obs == 0 ~ eps,
                                        .default = p_obs)) %>% 
    # renormalize
    tidytable::mutate(p_obs = p_obs / sum(p_obs), 
                      .by = c(sim, selex_type, comp_type))
}



# vector of selectivity types tested
selex_t <- unique(data$iss$selex_type)

# run for wtd expansion, iid
rr_wtd_iid <- purrr::map(1:length(selex_t),
                         ~est_logistic_normal(start_sigma = log(10),
                                              cov_strc = 'iid',
                                              data = data, 
                                              selex_t = selex_t[.],
                                              comp_t = 'wtd'))

# run for unwtd expansion, iid
rr_unwtd_iid <- purrr::map(1:length(selex_t),
                           ~est_logistic_normal(start_sigma = log(10),
                                                cov_strc = 'iid',
                                                data = data, 
                                                selex_t = selex_t[.],
                                                comp_t = 'unwtd'))

# run for wtd expansion, 1DAR1
rr_wtd_1DAR1 <- purrr::map(1:length(selex_t),
                           ~est_logistic_normal(start_sigma = log(10),
                                                start_rho = 0.1,
                                                cov_strc = '1DAR1',
                                                data = data, 
                                                selex_t = selex_t[.],
                                                comp_t = 'wtd'))

# run for unwtd expansion, 1DAR1
rr_unwtd_1DAR1 <- purrr::map(1:length(selex_t),
                             ~est_logistic_normal(start_sigma = log(10),
                                                  start_rho = 0.1,
                                                  cov_strc = '1DAR1',
                                                  data = data, 
                                                  selex_t = selex_t[.],
                                                  comp_t = 'unwtd'))







# restructure results

iss_sim %>% 
  tidytable::left_join(
    # wtd iid
    do.call(mapply, c(list, rr_wtd_iid, SIMPLIFY = FALSE))$res %>% 
      tidytable::map_df(., ~as.data.frame(.x), .id = "sel_t") %>% 
      tidytable::select(selex_type, comp_type, sigma_iid = sigma) %>% 
      tidytable::bind_rows(
        # unwtd iid
        do.call(mapply, c(list, rr_unwtd_iid, SIMPLIFY = FALSE))$res %>% 
          tidytable::map_df(., ~as.data.frame(.x), .id = "sel_t") %>% 
          tidytable::select(selex_type, comp_type, sigma_iid = sigma)) %>% 
      tidytable::left_join(
        # wtd 1DAR1
        do.call(mapply, c(list, rr_wtd_1DAR1, SIMPLIFY = FALSE))$res %>% 
          tidytable::map_df(., ~as.data.frame(.x), .id = "sel_t") %>% 
          tidytable::select(selex_type, comp_type, sigma_1DAR1 = sigma, rho_1DAR1 = rho) %>% 
          tidytable::bind_rows(
            # unwtd 1DAR1
            do.call(mapply, c(list, rr_unwtd_1DAR1, SIMPLIFY = FALSE))$res %>% 
              tidytable::map_df(., ~as.data.frame(.x), .id = "sel_t") %>% 
              tidytable::select(selex_type, comp_type, sigma_1DAR1 = sigma, rho_1DAR1 = rho))))





