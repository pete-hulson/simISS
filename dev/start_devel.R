
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

## run sim loop ----
tictoc::tic()
rr_sim <- purrr::map(1:iters, ~sim_comp(su_num, sim_popn, su_samp, p_su_samp))
runtime_sim <- tictoc::toc()

 # unlist results
do.call(mapply, c(list, rr_sim, SIMPLIFY = FALSE))$comp %>% 
  tidytable::map_df(., ~as.data.frame(.x), .id = "sim") -> res_sim

# save results
saveRDS(res_sim,
        file = here::here('output', 'res_sim.rds'))

## plot results ----
res_sim %>% 
  tidytable::pivot_longer(cols = c(samp_p_wtd, samp_p_unwtd)) %>% 
  tidytable::select(cat, comp_type = name, comp = value) %>% 
  tidytable::bind_rows(res_sim %>% 
                         tidytable::filter(sim == sample(1:iters, 1)) %>% 
                         tidytable::select(-sim) %>% 
                         tidytable::rename(rand_wtd = samp_p_wtd,
                                           rand_unwtd = samp_p_unwtd) %>% 
                         tidytable::pivot_longer(cols = c(rand_wtd, rand_unwtd)) %>% 
                         tidytable::rename(comp_type = name, comp = value)) %>% 
  tidytable::bind_rows(sim_popn$p_true %>% 
                         tidytable::mutate(comp_type = 'true') %>% 
                         tidytable::rename(comp = p_true)) -> plot_dat

sim_plot <- ggplot(data = plot_dat, aes(x = as.factor(cat), y = comp)) +
  geom_bar(data = plot_dat %>% tidytable::filter(comp_type %in% c('true')), stat = 'identity', 
           aes(fill = comp_type), alpha = 0.5) +
  geom_boxplot(data = plot_dat %>% tidytable::filter(comp_type %in% c('samp_p_wtd', 'samp_p_unwtd')), 
               aes(fill = comp_type)) +
  geom_line(data = plot_dat %>% tidytable::filter(comp_type %in% c('rand_wtd', 'rand_unwtd')) %>% tidytable::left_join(data.frame(comp_type = c('rand_wtd', 'rand_unwtd'), col = c(scico::scico(3, palette = 'roma')[1], scico::scico(3, palette = 'roma')[2]))), 
            aes(x = cat, col = comp_type), linetype = 'dashed') +
  scico::scale_fill_scico_d(palette = 'roma')  +
  scale_color_manual(values = c(scico::scico(3, palette = 'roma')[1], scico::scico(3, palette = 'roma')[2])) +
  theme_bw() +
  xlab('category')

ggsave(filename = "sim_popn.png",
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

# experiment 2: ----
# what are the factors that influence input sample size?

# number of simulation replicates for testing iss axes of influence
sim_reps <- 2


## test expansion weighting & pop'n structure ----

# run simulation
tictoc::tic()
rr_exp <- purrr::map(1:sim_reps, ~rep_sim(d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters))
runtime_exp <- tictoc::toc()

# plot results
plot_exp(rr_exp)


## test pop'n unit structure (spread around mean category) ----

# set levels of cv
pu_cv_test <- c(0.05, 0.1, 0.25, 1)

# run simulation
tictoc::tic()
rr_cv <- purrr::map(1:sim_reps, ~purrr::map(1:length(pu_cv_test), ~rep_sim(d, pu, pc, pu_cv = pu_cv_test[.], su_num, su_samp, p_su_samp, iters)))
runtime_cv <- tictoc::toc()

# save & plot results
plot_sim(rr_cv, 'cv', pu_cv_test, "CV", fact_perc = TRUE)


## test number of pop'n units ----

# set numbers of pop'n units
npu_test <- c(5, 10, 25, 100)

# run simulation
tictoc::tic()
rr_npu <- purrr::map(1:sim_reps, ~purrr::map(1:length(npu_test), ~rep_sim(d, pu = npu_test[.], pc, pu_cv, su_num, su_samp, p_su_samp, iters)))
runtime_npu <- tictoc::toc()

# save & plot results
plot_sim(rr_npu, 'npu', npu_test, "N_PU")


## test number of categories (i.e., longevity, growth) ----

# set numbers of categories
cat_test <- c(10, 15, 25, 50)

# run simulation
tictoc::tic()
rr_cat <- purrr::map(1:sim_reps, ~purrr::map(1:length(cat_test), ~rep_sim(d, pu, pc = cat_test[.], pu_cv, su_num, su_samp, p_su_samp, iters)))
runtime_cat <- tictoc::toc()

# save & plot results
plot_sim(rr_cat, 'cat', cat_test, "N_Cat")


## test number of sampling units (i.e., number of hauls) ----

# set numbers of categories
su_test <- c(50, 100, 250, 500)

# run simulation
tictoc::tic()
rr_su <- purrr::map(1:sim_reps, ~purrr::map(1:length(su_test), ~rep_sim(d, pu, pc, pu_cv, su_num = su_test[.], su_samp, p_su_samp, iters)))
runtime_su <- tictoc::toc()

# save & plot results
plot_sim(rr_su, 'su', su_test, "N_SU")


## test sample size within sampling units ----

# set numbers of categories
samp_test <- c(20, 50, 100, 250)

# run simulation
tictoc::tic()
rr_samp <- purrr::map(1:sim_reps, ~purrr::map(1:length(samp_test), ~rep_sim(d, pu, pc, pu_cv, su_num, su_samp = c(samp_test[.], 10), p_su_samp = c(1, 0), iters)))
runtime_samp <- tictoc::toc()

# save & plot results
plot_sim(rr_samp, 'samp', samp_test, "N_Samp")

# calc runtime for 500 simulations
((runtime_exp$toc - runtime_exp$tic) + 
    (runtime_cv$toc - runtime_cv$tic) + 
    (runtime_npu$toc - runtime_npu$tic) + 
    (runtime_cat$toc - runtime_cat$tic) + 
    (runtime_su$toc - runtime_su$tic) + 
    (runtime_samp$toc - runtime_samp$tic)) / (60 * sim_reps) * 250 / 60










# weighting reduces the iss, both result in a sim that results in reduction from nss









# (runtime$toc - runtime$tic) / (60 * sim_reps) * 500 / 60 # time for 500 reps in hours













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
















