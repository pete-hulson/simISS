
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


# experiment 2: ----
# what are the factors that influence input sample size?

# number of simulation replicates for testing iss axes of influence
sim_reps <- 10

## test expansion weighting & pop'n structure ----

# run simulation
tictoc::tic()
rr_exp <- purrr::map(1:sim_reps, ~rep_sim(d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters))
runtime_exp <- tictoc::toc()

# calc runtime test
(runtime_exp$toc - runtime_exp$tic) / (60 * sim_reps) * 1000 / 60

# save & plot results
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











