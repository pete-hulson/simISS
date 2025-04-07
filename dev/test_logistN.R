
# load/source stuff ----
library(tidyverse)
library(future)
source(here::here('R', 'base_functions.R'))
source(here::here('R', 'stats_functions.R'))


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




# notes ----
# want to test whether taking mean of sigmas estimated within an iteration gives similar
# value to estimate of sigma across combined iterations
#
# mean across iterations tends to be smaller than the sigma using all iterations
# that is, uncertainty is larger estimating sigma across combined iteration dataframe

# estimate logsitN params for each iter ----
combs_it <- tidytable::expand_grid(iter = 1:iters,
                                   selex = unique(data$iss$selex_type), 
                                   comp = unique(data$iss$comp_type))

# iid
tictoc::tic()
future::plan(multisession, workers = 6)
res_iis_iter <- par_run_iid(data, combs_it)
future::plan(sequential)
runtime_test_iid_iter <- tictoc::toc()

# 1DAR1
tictoc::tic()
future::plan(multisession, workers = 6)
res_1DAR1_iter <- par_run_1DAR1(data, combs_it)
future::plan(sequential)
runtime_test_iid_iter <- tictoc::toc()


res_iis_iter %>% 
  tidytable::summarise(mu_sigma = median(sigma),
                       sd_sigma = sd(sigma), .by = c(selex_type, comp_type))


# estimate logistN params overall ----
# dataframe of selectivity/composition expansion types tested

combs <- tidytable::expand_grid(selex = unique(data$iss$selex_type), 
                                comp = unique(data$iss$comp_type))

# run for iid
rr_iid <- purrr::map(1:dim(combs)[1],
                     ~est_logistic_normal(start_sigma = log(10),
                                          cov_strc = 'iid',
                                          data = data, 
                                          selex_t = combs$selex[.],
                                          comp_t = combs$comp[.]))

# run for 1DAR1
rr_1DAR1 <- purrr::map(1:dim(combs)[1],
                       ~est_logistic_normal(start_sigma = log(10),
                                            start_rho = 0.1,
                                            cov_strc = '1DAR1',
                                            data = data, 
                                            selex_t = combs$selex[.],
                                            comp_t = combs$comp[.]))

# unlist results
iss_sim %>% 
  tidytable::left_join(
    # iid
    do.call(mapply, c(list, rr_iid, SIMPLIFY = FALSE))$res %>% 
      tidytable::map_df(., ~as.data.frame(.x), .id = "comb") %>% 
      tidytable::select(selex_type, comp_type, sigma_iid = sigma) %>% 
      tidytable::left_join(
        # 1DAR1
        do.call(mapply, c(list, rr_1DAR1, SIMPLIFY = FALSE))$res %>% 
          tidytable::map_df(., ~as.data.frame(.x), .id = "comb") %>% 
          tidytable::select(selex_type, comp_type, sigma_1DAR1 = sigma, rho_1DAR1 = rho)))

res_test %>% 
  tidytable::summarise(sigma = mean(sigma),
                       sd_sigma = sd(sigma), .by = c(selex_type, comp_type))




#' function to test iid for each iter in parallel (6 cores)
#' 
par_run_iid <- function(data, combs_it){
  run1 %<-% purrr::map(1:1000, 
                       ~est_logistic_normal(start_sigma = log(10),
                                            cov_strc = 'iid',
                                            data = list(exp = data$exp,
                                                        obs = data$obs[sim == combs_it$iter[.]],
                                                        iss = data$iss), 
                                            selex_t = combs_it$selex[.],
                                            comp_t = combs_it$comp[.])) %seed% TRUE
  run2 %<-% purrr::map(1001:2000, 
                       ~est_logistic_normal(start_sigma = log(10),
                                            cov_strc = 'iid',
                                            data = list(exp = data$exp,
                                                        obs = data$obs[sim == combs_it$iter[.]],
                                                        iss = data$iss), 
                                            selex_t = combs_it$selex[.],
                                            comp_t = combs_it$comp[.])) %seed% TRUE
  run3 %<-% purrr::map(2001:3000, 
                       ~est_logistic_normal(start_sigma = log(10),
                                            cov_strc = 'iid',
                                            data = list(exp = data$exp,
                                                        obs = data$obs[sim == combs_it$iter[.]],
                                                        iss = data$iss), 
                                            selex_t = combs_it$selex[.],
                                            comp_t = combs_it$comp[.])) %seed% TRUE
  run4 %<-% purrr::map(3001:4000, 
                       ~est_logistic_normal(start_sigma = log(10),
                                            cov_strc = 'iid',
                                            data = list(exp = data$exp,
                                                        obs = data$obs[sim == combs_it$iter[.]],
                                                        iss = data$iss), 
                                            selex_t = combs_it$selex[.],
                                            comp_t = combs_it$comp[.])) %seed% TRUE
  run5 %<-% purrr::map(4001:5000, 
                       ~est_logistic_normal(start_sigma = log(10),
                                            cov_strc = 'iid',
                                            data = list(exp = data$exp,
                                                        obs = data$obs[sim == combs_it$iter[.]],
                                                        iss = data$iss), 
                                            selex_t = combs_it$selex[.],
                                            comp_t = combs_it$comp[.])) %seed% TRUE
  run6 %<-% purrr::map(5001:6000, 
                       ~est_logistic_normal(start_sigma = log(10),
                                            cov_strc = 'iid',
                                            data = list(exp = data$exp,
                                                        obs = data$obs[sim == combs_it$iter[.]],
                                                        iss = data$iss), 
                                            selex_t = combs_it$selex[.],
                                            comp_t = combs_it$comp[.])) %seed% TRUE
  
  
  res <- do.call(mapply, c(list, run1, SIMPLIFY = FALSE))$res %>% 
    tidytable::map_df(., ~as.data.frame(.x), .id = "sim") %>% 
    tidytable::bind_rows(do.call(mapply, c(list, run2, SIMPLIFY = FALSE))$res %>% 
                           tidytable::map_df(., ~as.data.frame(.x), .id = "sim")) %>% 
    tidytable::bind_rows(do.call(mapply, c(list, run3, SIMPLIFY = FALSE))$res %>% 
                           tidytable::map_df(., ~as.data.frame(.x), .id = "sim")) %>% 
    tidytable::bind_rows(do.call(mapply, c(list, run4, SIMPLIFY = FALSE))$res %>% 
                           tidytable::map_df(., ~as.data.frame(.x), .id = "sim")) %>% 
    tidytable::bind_rows(do.call(mapply, c(list, run5, SIMPLIFY = FALSE))$res %>% 
                           tidytable::map_df(., ~as.data.frame(.x), .id = "sim")) %>% 
    tidytable::bind_rows(do.call(mapply, c(list, run6, SIMPLIFY = FALSE))$res %>% 
                           tidytable::map_df(., ~as.data.frame(.x), .id = "sim")) %>% 
    tidytable::select(-sim)
  
  res
}

#' function to test 1DAR1 for each iter in parallel (6 cores)
#' 
par_run_1DAR1 <- function(data, combs_it){
  run1 %<-% purrr::map(1:1000, 
                       ~est_logistic_normal(start_sigma = log(10),
                                            cov_strc = '1DAR1',
                                            data = list(exp = data$exp,
                                                        obs = data$obs[sim == combs_it$iter[.]],
                                                        iss = data$iss), 
                                            selex_t = combs_it$selex[.],
                                            comp_t = combs_it$comp[.])) %seed% TRUE
  run2 %<-% purrr::map(1001:2000, 
                       ~est_logistic_normal(start_sigma = log(10),
                                            cov_strc = '1DAR1',
                                            data = list(exp = data$exp,
                                                        obs = data$obs[sim == combs_it$iter[.]],
                                                        iss = data$iss), 
                                            selex_t = combs_it$selex[.],
                                            comp_t = combs_it$comp[.])) %seed% TRUE
  run3 %<-% purrr::map(2001:3000, 
                       ~est_logistic_normal(start_sigma = log(10),
                                            cov_strc = '1DAR1',
                                            data = list(exp = data$exp,
                                                        obs = data$obs[sim == combs_it$iter[.]],
                                                        iss = data$iss), 
                                            selex_t = combs_it$selex[.],
                                            comp_t = combs_it$comp[.])) %seed% TRUE
  run4 %<-% purrr::map(3001:4000, 
                       ~est_logistic_normal(start_sigma = log(10),
                                            cov_strc = '1DAR1',
                                            data = list(exp = data$exp,
                                                        obs = data$obs[sim == combs_it$iter[.]],
                                                        iss = data$iss), 
                                            selex_t = combs_it$selex[.],
                                            comp_t = combs_it$comp[.])) %seed% TRUE
  run5 %<-% purrr::map(4001:5000, 
                       ~est_logistic_normal(start_sigma = log(10),
                                            cov_strc = '1DAR1',
                                            data = list(exp = data$exp,
                                                        obs = data$obs[sim == combs_it$iter[.]],
                                                        iss = data$iss), 
                                            selex_t = combs_it$selex[.],
                                            comp_t = combs_it$comp[.])) %seed% TRUE
  run6 %<-% purrr::map(5001:6000, 
                       ~est_logistic_normal(start_sigma = log(10),
                                            cov_strc = '1DAR1',
                                            data = list(exp = data$exp,
                                                        obs = data$obs[sim == combs_it$iter[.]],
                                                        iss = data$iss), 
                                            selex_t = combs_it$selex[.],
                                            comp_t = combs_it$comp[.])) %seed% TRUE
  
  
  res <- do.call(mapply, c(list, run1, SIMPLIFY = FALSE))$res %>% 
    tidytable::map_df(., ~as.data.frame(.x), .id = "sim") %>% 
    tidytable::bind_rows(do.call(mapply, c(list, run2, SIMPLIFY = FALSE))$res %>% 
                           tidytable::map_df(., ~as.data.frame(.x), .id = "sim")) %>% 
    tidytable::bind_rows(do.call(mapply, c(list, run3, SIMPLIFY = FALSE))$res %>% 
                           tidytable::map_df(., ~as.data.frame(.x), .id = "sim")) %>% 
    tidytable::bind_rows(do.call(mapply, c(list, run4, SIMPLIFY = FALSE))$res %>% 
                           tidytable::map_df(., ~as.data.frame(.x), .id = "sim")) %>% 
    tidytable::bind_rows(do.call(mapply, c(list, run5, SIMPLIFY = FALSE))$res %>% 
                           tidytable::map_df(., ~as.data.frame(.x), .id = "sim")) %>% 
    tidytable::bind_rows(do.call(mapply, c(list, run6, SIMPLIFY = FALSE))$res %>% 
                           tidytable::map_df(., ~as.data.frame(.x), .id = "sim")) %>% 
    tidytable::select(-sim)
  
  res
}
