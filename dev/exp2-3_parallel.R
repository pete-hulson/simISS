
# load/source stuff ----
library(tidyverse)
library(future)
source(here::here('R', 'functions.R'))


# set up experiment parameters ----

# number of simulation replicates for testing iss axes of influence
sim_reps <- 1000

# number of bootstrap replicates
bs_iters <- 1000

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


# experiment 2: ----
# what are the factors that influence input sample size?

# run exp2 tests in parallel
tictoc::tic()
future::plan(multisession, workers = 6)
run_exp2_tests(sim_reps, d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters)
future::plan(sequential)
runtime_test_exp2 <- tictoc::toc()


# experiment 3: ----
# can a given realization of a sample give the 'true' iss back?

# run exp3 in parallel
tictoc::tic()
future::plan(multisession, workers = 10)
run_bs_test(bs_iters, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters)
future::plan(sequential)
runtime_test_exp3 <- tictoc::toc()


# calc runtime for 500 runs ----
# experiment 2
(runtime_test_exp2$toc - runtime_test_exp2$tic) / (60 * sim_reps) * 1000 / 60
# experiment 3
(runtime_test_exp3$toc - runtime_test_exp3$tic) / (60 * bs_iters) * 1000 / 60


