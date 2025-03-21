
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


# experiment 3: can a given realization of a sample give the 'true' iss back? ----

# number of bootstrap replicates
bs_iters <- 10

## get a sample realization ----

# generate log-normal catch
c_realize <- data.frame(haul = 1:haul_num, c_h = exp(stats::rnorm(haul_num, 0, 1))) 

# determine which school is being sampled (based on relative abundance)
sc_smp <- colSums(rmultinom(haul_num, 1, p_sc) * 1:3)

# generate a realization of samples across hauls
rr <- purrr::map(1:haul_num, ~gen_haul_samp(sc_smp[.]))

do.call(mapply, c(list, rr, SIMPLIFY = FALSE))$samp_h %>% 
  tidytable::map_df(., ~as.data.frame(.x), .id = "haul") %>% 
  tidytable::rename(cat = V1,
                    samp_h = V2) -> haul_realize

# compute the realization composition
haul_realize %>% 
  tidytable::left_join(haul_realize %>% 
                         tidytable::summarise(haul_samp = sum(samp_h), .by = haul) %>% 
                         tidytable::mutate(prop_s = haul_samp / sum(haul_samp))) %>% 
  tidytable::left_join(c_realize %>% 
                         tidytable::mutate(prop_c = c_h / sum(c_h))) %>% 
  tidytable::summarise(samp_wtd = sum(samp_h * prop_s * prop_c), # weight samples by sample and catch proportion
                       samp_unwtd = sum(samp_h), # do not weight samples
                       .by = cat) %>% 
  tidytable::mutate(real_p_wtd = samp_wtd / sum(samp_wtd),
                    real_p_unwtd = samp_unwtd / sum(samp_unwtd)) %>% 
  tidytable::select(cat, real_p_wtd, real_p_unwtd) -> p_realize


## define bootstrap functions ----

bs_haul_samp <- function(haul_num,
                         haul_realize,
                         c_realize){
  
  # bootstrap hauls
  bs_h <- data.frame(haul = sample(1:haul_num, haul_num, replace = TRUE)) %>% 
    tidytable::mutate(id = .I)
  
  # get bootstrap haul sample size results
  bs_h %>% 
    tidytable::left_join(haul_realize) %>% 
    tidytable::summarise(haul_samp = sum(samp_h), .by = id) %>% 
    tidytable::mutate(prop_s = haul_samp / sum(haul_samp)) %>% 
    tidytable::rename(haul = id) -> bs_prop_s
  
  bs_prop_s %>% 
    tidytable::summarise(nss = sum(haul_samp)) -> bs_nss
  
  # get catch proportions
  bs_h %>% 
    tidytable::left_join(c_realize) %>% 
    tidytable::mutate(prop_c = c_h / sum(c_h)) %>% 
    tidytable::select(haul = id, c_h, prop_c) -> bs_prop_c
  
  # get bootstrap comp results
  bs_h %>% 
    tidytable::left_join(haul_realize) %>% 
    tidytable::select(haul = id, cat, samp_h) %>% 
    tidytable::left_join(bs_prop_s) %>% 
    tidytable::left_join(bs_prop_c) %>% 
    tidytable::summarise(samp_wtd = sum(samp_h * prop_s * prop_c), # weight samples by sample and catch proportion
                         samp_unwtd = sum(samp_h), # do not weight samples
                         .by = cat) %>% 
    tidytable::mutate(samp_p_wtd = samp_wtd / sum(samp_wtd),
                      samp_p_unwtd = samp_unwtd / sum(samp_unwtd)) %>% 
    tidytable::select(cat, samp_p_wtd, samp_p_unwtd) -> bs_comp_it
  
  list(bs_comp_it = bs_comp_it, bs_nss = bs_nss)
}

sim_bs <- function(haul_num, iters){
  
  # get a sample realization
  
  # generate log-normal catch
  c_realize <- data.frame(haul = 1:haul_num, c_h = exp(stats::rnorm(haul_num, 0, 1))) 
  
  # determine which school is being sampled (based on relative abundance)
  sc_smp <- colSums(rmultinom(haul_num, 1, p_sc) * 1:3)
  
  # generate a realization of samples across hauls
  rr <- purrr::map(1:haul_num, ~gen_haul_samp(sc_smp[.]))
  
  do.call(mapply, c(list, rr, SIMPLIFY = FALSE))$samp_h %>% 
    tidytable::map_df(., ~as.data.frame(.x), .id = "haul") %>% 
    tidytable::rename(cat = V1,
                      samp_h = V2) -> haul_realize
  
  # compute the realization composition
  haul_realize %>% 
    tidytable::left_join(haul_realize %>% 
                           tidytable::summarise(haul_samp = sum(samp_h), .by = haul) %>% 
                           tidytable::mutate(prop_s = haul_samp / sum(haul_samp))) %>% 
    tidytable::left_join(c_realize %>% 
                           tidytable::mutate(prop_c = c_h / sum(c_h))) %>% 
    tidytable::summarise(samp_wtd = sum(samp_h * prop_s * prop_c), # weight samples by sample and catch proportion
                         samp_unwtd = sum(samp_h), # do not weight samples
                         .by = cat) %>% 
    tidytable::mutate(real_p_wtd = samp_wtd / sum(samp_wtd),
                      real_p_unwtd = samp_unwtd / sum(samp_unwtd)) %>% 
    tidytable::select(cat, real_p_wtd, real_p_unwtd) -> p_realize
  
  # bootstrap the realization
  rr <- purrr::map(1:iters, ~bs_haul_samp(haul_num, haul_realize, c_realize))
  
  
  do.call(mapply, c(list, rr, SIMPLIFY = FALSE))$bs_comp_it %>% 
    tidytable::map_df(., ~as.data.frame(.x), .id = "sim") -> res_bs
  
  res_bs %>% 
    tidytable::left_join(p_realize) %>% 
    tidytable::summarise(rss_wtd = sum(real_p_wtd * (1- real_p_wtd)) / sum((samp_p_wtd - real_p_wtd) ^ 2),
                         rss_unwtd = sum(real_p_unwtd * (1- real_p_unwtd)) / sum((samp_p_unwtd - real_p_unwtd) ^ 2),
                         .by = sim) -> rss_res_bs
  
  rss_res_bs %>% 
    tidytable::left_join(do.call(mapply, c(list, rr, SIMPLIFY = FALSE))$bs_nss %>% 
                           tidytable::map_df(., ~as.data.frame(.x), .id = "sim")) %>% 
    tidytable::summarise(iss_wtd = psych::harmonic.mean(rss_wtd, zero = FALSE),
                         iss_unwtd = psych::harmonic.mean(rss_unwtd, zero = FALSE),
                         mean_nss = mean(nss)) -> iss_res_bs
  
  list(iss_res_bs = iss_res_bs)
}

## run bootstrap ----

purrr::map(1:bs_iters, ~sim_bs(haul_num, iters))
















