# experiment 2: ----
# what are the factors that influence input sample size?


# load libraries/source fcns ----
library(tidyverse)
library(future)
source(here::here('R', 'base_functions.R'))
source(here::here('R', 'stats_functions.R'))
source(here::here('R', 'exp2_functions.R'))


# set up experiment parameters ----

## simulation/bootstrap parameters ----

# full run?
full_run = FALSE

# number of simulation replicates for testing iss axes of influence
sim_reps <- 5

# number of desired simulation replicates
X <- 500

# number of bootstrap iterations
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


# run exp2 tests in parallel (on 7 cores) ----
if(isTRUE(full_run)){
  tictoc::tic()
  future::plan(multisession, workers = 7)
  run_exp2_tests(X, d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters, cov_strc = 'iid')
  future::plan(sequential)
  runtime_test_exp2 <- tictoc::toc()
}else{
  tictoc::tic()
  future::plan(multisession, workers = 7)
  run_exp2_tests(sim_reps, d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters, cov_strc = 'iid')
  future::plan(sequential)
  runtime_test_exp2 <- tictoc::toc()
}


# calc runtime for X simulations
if(isTRUE(full_run)){
  cat("Run took", round((runtime_test_exp2$toc - runtime_test_exp2$tic) / 60 / 60, digits = 1), "hrs")
} else{
  cat("Run is estimated to take", round((runtime_test_exp2$toc - runtime_test_exp2$tic) / (60 * sim_reps) * X / 60, digits = 1), "hrs")
}
