
# load/source stuff ----
library(tidyverse)
source(here::here('R', 'base_functions.R'))
source(here::here('R', 'stats_functions.R'))
source(here::here('R', 'exp2_functions.R'))



# set up experiment parameters ----

## simulation/bootstrap parameters ----

# full run?
full_run = FALSE

# number of simulation replicates for testing iss axes of influence
sim_reps <- 25

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
pu <- 5

# number of population categories (e.g., ages or length bins)
pc <- 15

# CV around mean within pop'n unit (e.g., CV in mean age of school)
pu_cv <- 0.2

# pop'n exponential decay (e.g., lnM, tied to inverse of number of pop'n categories so that pop'n = 0.01 at largest category)
d <- log(0.01) / (1 - pc)


# experiment 2: ----
# what are the factors that influence input sample size?

# select which tests to run
# all tests
# tests <- c('base', 'CV', 'PU', 'C', 'SU', 'nSU250', 'nSU500')
# selected tests
tests <- c('base')

# full run
if(isTRUE(full_run)){
  # test expansion weighting, selex, & pop'n structure
  if('base' %in% tests){
    runtime_base <- test_base(X, d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters, cov_strc = c('iid', '1DAR1'))
  } else{tictoc::tic();runtime_base <- tictoc::toc()}
  
  # test pop'n unit structure (spread around mean category)
  if('CV' %in% tests){
    runtime_CV <- test_CV(X, d, pu, pc, su_num, su_samp, p_su_samp, iters, cov_strc = c('iid', '1DAR1'))
  } else{tictoc::tic();runtime_CV <- tictoc::toc()}

  # test number of categories (i.e., longevity, growth)
  if('C' %in% tests){
    runtime_C <- test_C(X, d, pu, pu_cv, su_num, su_samp, p_su_samp, iters, cov_strc = c('iid', '1DAR1'))
  } else{tictoc::tic();runtime_C <- tictoc::toc()}
  
  # test number of sampling units (i.e., number of hauls)
  if('SU' %in% tests){
    runtime_SU <- test_SU(X, d, pu, pc, pu_cv, su_samp, p_su_samp, iters, cov_strc = c('iid', '1DAR1'))
  } else{tictoc::tic();runtime_SU <- tictoc::toc()}
  
  # test sample size within sampling units (for 250 sampling units)
  if('nSU250' %in% tests){
    runtime_nSU250 <- test_nSU(X, d, pu, pc, pu_cv, su_num, iters, 'S250', cov_strc = c('iid', '1DAR1'))
  } else{tictoc::tic();runtime_nSU250 <- tictoc::toc()}

} else{ # test runs
  # test expansion weighting, selex, & pop'n structure
  if('base' %in% tests){
    runtime_base <- test_base(sim_reps, d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters, cov_strc = c('iid', '1DAR1'))
  } else{tictoc::tic();runtime_base <- tictoc::toc()}
  
  # test pop'n unit structure (spread around mean category)
  if('CV' %in% tests){
    runtime_CV <- test_CV(sim_reps, d, pu, pc, su_num, su_samp, p_su_samp, iters, cov_strc = c('iid', '1DAR1'))
  } else{tictoc::tic();runtime_CV <- tictoc::toc()}

  # test number of categories (i.e., longevity, growth)
  if('C' %in% tests){
    runtime_C <- test_C(sim_reps, d, pu, pu_cv, su_num, su_samp, p_su_samp, iters, cov_strc = c('iid', '1DAR1'))
  } else{tictoc::tic();runtime_C <- tictoc::toc()}
  
  # test number of sampling units (i.e., number of hauls)
  if('SU' %in% tests){
    runtime_SU <- test_SU(sim_reps, d, pu, pc, pu_cv, su_samp, p_su_samp, iters, cov_strc = c('iid', '1DAR1'))
  } else{tictoc::tic();runtime_SU <- tictoc::toc()}
  
  # test sample size within sampling units (for 250 sampling units)
  if('nSU250' %in% tests){
    runtime_nSU250 <- test_nSU(sim_reps, d, pu, pc, pu_cv, su_num, iters, 'S250', cov_strc = c('iid', '1DAR1'))
  } else{tictoc::tic();runtime_nSU250 <- tictoc::toc()}

}


# calc runtime for X replicates
if(isTRUE(full_run)){
  cat("Run took", round(((runtime_base$toc - runtime_base$tic) +
                           (runtime_CV$toc - runtime_CV$tic) +
                           (runtime_C$toc - runtime_C$tic) +
                           (runtime_SU$toc - runtime_SU$tic) +
                           (runtime_nSU250$toc - runtime_nSU250$tic)) / 60 / 60, digits = 1), "hrs")

} else{
    cat("Run is estimated to take", round(((runtime_base$toc - runtime_base$tic) +
                                            (runtime_CV$toc - runtime_CV$tic) +
                                            (runtime_C$toc - runtime_C$tic) +
                                            (runtime_SU$toc - runtime_SU$tic) +
                                            (runtime_nSU250$toc - runtime_nSU250$tic)) / 60 / 60  * X / sim_reps, digits = 1), "hrs")
}

