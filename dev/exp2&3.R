# run both experiment 2 & 3 ----


# load libraries/source functions ----
library(tidyverse)
library(future)
source(here::here('R', 'base_functions.R'))
source(here::here('R', 'exp2_functions.R'))
source(here::here('R', 'exp3_functions.R'))


# set up experiment parameters ----

## simulation/bootstrap parameters ----

# full run?
full_run = FALSE

# number of simulation replicates for testing iss axes of influence (exp2)
sim_reps <- 5

# number of bootstrap replicates (exp3)
bs_iters <- 5

# number of desired bootstrap replicates
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


# experiment 3: ----

# get number of cores
numCore <- parallel::detectCores()

# run on 10 cores
if(numCore > 10){
  if(isTRUE(full_run)){
    tictoc::tic()
    future::plan(multisession, workers = 10)
    run_bs_test(X, d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters, numCore)
    future::plan(sequential)
    runtime_test_exp3 <- tictoc::toc()
  } else{
    tictoc::tic()
    future::plan(multisession, workers = 10)
    run_bs_test(bs_iters, d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters, numCore)
    future::plan(sequential)
    runtime_test_exp3 <- tictoc::toc()
  }
}

# run on 7 cores
if(numCore < 10){
  if(isTRUE(full_run)){
    tictoc::tic()
    future::plan(multisession, workers = 7)
    run_bs_test(bs_iters = round(10 * X / 7, digits = 0), d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters, numCore)
    future::plan(sequential)
    runtime_test_exp3 <- tictoc::toc()
  } else{
    tictoc::tic()
    future::plan(multisession, workers = 7)
    run_bs_test(bs_iters, d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters, numCore)
    future::plan(sequential)
    runtime_test_exp3 <- tictoc::toc()
  }
}


# experiment 2: ----

if(isTRUE(full_run)){
  tictoc::tic()
  future::plan(multisession, workers = 7)
  run_exp2_tests(X, d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters)
  future::plan(sequential)
  runtime_test_exp2 <- tictoc::toc()
}else{
  tictoc::tic()
  future::plan(multisession, workers = 7)
  run_exp2_tests(sim_reps, d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters)
  future::plan(sequential)
  runtime_test_exp2 <- tictoc::toc()
}


# calc runtime test for X runs ----
if(isTRUE(full_run)){
  cat("Exp2 run took", (runtime_test_exp2$toc - runtime_test_exp2$tic) / 60 / 60, "hrs\n")
  cat("Exp3 run took", (runtime_test_exp3$toc - runtime_test_exp3$tic) / 60 / 60, "hrs")
} else{
  runtime_exp2 <- round((runtime_test_exp2$toc - runtime_test_exp2$tic) / (60 * sim_reps) * X / 60, digits = 1)
  cat("Exp2 run is estimated to take", runtime_exp2, "hrs\n")
  if(numCore > 10){
    runtime_exp3 <- round((runtime_test_exp3$toc - runtime_test_exp3$tic) / (60 * bs_iters) * X / 60, digits = 1)
    cat("Exp3 run is estimated to take", runtime_exp3, "hrs\n")
  } else{
    runtime_exp3 <- round((runtime_test_exp3$toc - runtime_test_exp3$tic) / (60 * bs_iters) * round(10 * X / 7, digits = 0) / 60, digits = 1)
    cat("Exp3 run is estimated to take", runtime_exp3, "hrs\n")
  }
  cat("Exp 2&3 total run time estimated to take", round((runtime_exp2 + runtime_exp2) / 24, digits = 1), "days")
}
