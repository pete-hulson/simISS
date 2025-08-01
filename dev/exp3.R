# experiment 3: ----
# can bootstrapping a given realization of a sample give the 'true' iss back?


# load libraries/source fcns ----
library(tidyverse)
library(future)
library(forcats)
source(here::here('R', 'base_functions.R'))
source(here::here('R', 'stats_functions.R'))
source(here::here('R', 'exp3_functions.R'))

# set up experiment parameters ----

## simulation/bootstrap parameters ----

# full run?
full_run = FALSE

# number of bootstrap replicates
bs_iters <- 5

# number of desired bootstrap replicates
X <- 500

# number of bootstrap iterations
iters <- 500

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


# run exp3 in parallel ----

# get number of cores
numCore <- parallel::detectCores()

# run on 10 cores
if(numCore > 10){
  if(isTRUE(full_run)){
    tictoc::tic()
    future::plan(multisession, workers = 10)
    run_bs_test(X, d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters, boot_2nd = FALSE, numCore)
    future::plan(sequential)
    runtime_test_exp3 <- tictoc::toc()
  } else{
    tictoc::tic()
    future::plan(multisession, workers = 10)
    run_bs_test(bs_iters, d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters, boot_2nd = FALSE, numCore)
    future::plan(sequential)
    runtime_test_exp3 <- tictoc::toc()
  }
}

# run on 7 cores
if(numCore < 10){
  if(isTRUE(full_run)){
    tictoc::tic()
    future::plan(multisession, workers = 7)
    run_bs_test(bs_iters = round(10 * X / 7, digits = 0), d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters, boot_2nd = FALSE, numCore)
    future::plan(sequential)
    runtime_test_exp3 <- tictoc::toc()
  } else{
    tictoc::tic()
    future::plan(multisession, workers = 7)
    run_bs_test(bs_iters, d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters, boot_2nd = FALSE, numCore)
    future::plan(sequential)
    runtime_test_exp3 <- tictoc::toc()
  }
}

# calc runtime for X replicates
if(isTRUE(full_run)){
  cat("Run took", round((runtime_test_exp3$toc - runtime_test_exp3$tic) / 60 / 60, digits = 1), "hrs")
} else{
  if(numCore > 10){
    cat("Run is estimated to take", round((runtime_test_exp3$toc - runtime_test_exp3$tic) / (60 * bs_iters) * X / 60, digits = 1), "hrs")
  } else{
    cat("Run is estimated to take", round((runtime_test_exp3$toc - runtime_test_exp3$tic) / (60 * bs_iters) * round(10 * X / 7, digits = 0) / 60, digits = 1), "hrs")
  }
}
