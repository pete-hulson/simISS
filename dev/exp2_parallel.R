
# load/source stuff ----
library(tidyverse)
library(future)
source(here::here('R', 'functions.R'))


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
pu <- 25

# number of population categories (e.g., ages or length bins)
pc <- 15

# CV around mean within pop'n unit (e.g., CV in mean age of school)
pu_cv <- 0.25

# pop'n exponential decay (e.g., lnM, tied to inverse of number of pop'n categories so that pop'n = 0.01 at largest category)
d <- log(0.01) / (1 - pc)


# experiment 2: ----
# what are the factors that influence input sample size?

## test expansion weighting & pop'n structure ----
test1 <- function(sim_reps, d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters){
  #start timer
  tictoc::tic()
  # run simulation
  rr_exp <- purrr::map(1:sim_reps, ~rep_sim(d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters))
  # save & plot results
  plot_exp(rr_exp)
  # end timer
  runtime <- tictoc::toc()
  return(runtime)
}

## test pop'n unit structure (spread around mean category) ----
test2 <- function(sim_reps, d, pu, pc, su_num, su_samp, p_su_samp, iters){
  #start timer
  tictoc::tic()
  # set levels of cv
  pu_cv_test <- c(0.1, 0.25, 1, 100)
  # run simulation
  rr_cv <- purrr::map(1:sim_reps, ~purrr::map(1:length(pu_cv_test), ~rep_sim(d, pu, pc, pu_cv = pu_cv_test[.], su_num, su_samp, p_su_samp, iters)))
  # save & plot results
  plot_sim(rr = rr_cv, 
           plot_name = 'cv', 
           test_vec = pu_cv_test, 
           test_name = "CV",
           test_lab = 'Population unit CV around mean category',
           plot_nss = TRUE,
           fact_perc = TRUE)
  # end timer
  runtime <- tictoc::toc()
  return(runtime)
}

## test number of pop'n units ----
test3 <- function(sim_reps, d, pc, pu_cv, su_num, su_samp, p_su_samp, iters){
  #start timer
  tictoc::tic()
  # set numbers of pop'n units
  npu_test <- c(25, 100, 250, 1000)
  # run simulation
  rr_npu <- purrr::map(1:sim_reps, ~purrr::map(1:length(npu_test), ~rep_sim(d, pu = npu_test[.], pc, pu_cv, su_num, su_samp, p_su_samp, iters)))
  # save & plot results
  plot_sim(rr = rr_npu, 
           plot_name = 'Npu', 
           test_vec = npu_test, 
           test_name = "N",
           test_lab = 'Number of population units')
  # end timer
  runtime <- tictoc::toc()
  return(runtime)
}

## test number of categories (i.e., longevity, growth) ----
test4 <- function(sim_reps, d, pu, pu_cv, su_num, su_samp, p_su_samp, iters){
  #start timer
  tictoc::tic()
  # set numbers of categories
  cat_test <- c(10, 15, 25, 50)
  # run simulation
  rr_cat <- purrr::map(1:sim_reps, ~purrr::map(1:length(cat_test), ~rep_sim(d, pu, pc = cat_test[.], pu_cv, su_num, su_samp, p_su_samp, iters)))
  # save & plot results
  plot_sim(rr = rr_cat, 
           plot_name = 'Ncat', 
           test_vec = cat_test,
           test_name = "N",
           test_lab = 'Number of categories within the population')
  # end timer
  runtime <- tictoc::toc()
  return(runtime)
}

## test number of sampling units (i.e., number of hauls) ----
test5 <- function(sim_reps, d, pu, pc, pu_cv, su_samp, p_su_samp, iters){
  #start timer
  tictoc::tic()
  # set number of sampling units
  su_test <- c(100, 250, 500, 1000)
  # run simulation
  rr_su <- purrr::map(1:sim_reps, ~purrr::map(1:length(su_test), ~rep_sim(d, pu, pc, pu_cv, su_num = su_test[.], su_samp, p_su_samp, iters)))
  # save & plot results
  plot_sim(rr = rr_su, 
           plot_name = 'Nsu', 
           test_vec = su_test, 
           test_name = "N",
           test_lab = 'Number of sampling units')
  # end timer
  runtime <- tictoc::toc()
  return(runtime)
}

## test sample size within sampling units ----
test6 <- function(sim_reps, d, pu, pc, pu_cv, su_num, iters, plot_name){
  #start timer
  tictoc::tic()
  # set sample size within sampling units
  samp_test <- c(20, 50, 100, 250)
  # run simulation
  rr_samp <- purrr::map(1:sim_reps, ~purrr::map(1:length(samp_test), ~rep_sim(d, pu, pc, pu_cv, su_num, su_samp = c(samp_test[.], 10), p_su_samp = c(1, 0), iters)))
  # save & plot results
  plot_sim(rr = rr_samp, 
           plot_name = plot_name, 
           test_vec = samp_test, 
           test_name = "n",
           test_lab = 'Number of samples within a sampling unit')
  # end timer
  runtime <- tictoc::toc()
  return(runtime)
}

## run tests in parallel ----

# number of simulation replicates for testing iss axes of influence
sim_reps <- 10

# set up parallel fcn
run_tests <- function(sim_reps, d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters){
  # runtime_exp %<-% test1(sim_reps, d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters) %seed% TRUE
  runtime_cv %<-% test2(sim_reps, d, pu, pc, su_num, su_samp, p_su_samp, iters) %seed% TRUE
  runtime_npu %<-% test3(sim_reps, d, pc, pu_cv, su_num, su_samp, p_su_samp, iters) %seed% TRUE
  runtime_cat %<-% test4(sim_reps, d, pu, pu_cv, su_num, su_samp, p_su_samp, iters) %seed% TRUE
  runtime_su %<-% test5(sim_reps, d, pu, pc, pu_cv, su_samp, p_su_samp, iters) %seed% TRUE
  runtime_samp_h250 %<-% test6(sim_reps, d, pu, pc, pu_cv, su_num, iters, 'samp_h250') %seed% TRUE
  runtime_samp_h500 %<-% test6(sim_reps, d, pu, pc, pu_cv, su_num = 500, iters, 'samp_h500') %seed% TRUE
  runtimes <- c((runtime_exp$toc - runtime_exp$tic),
                (runtime_cv$toc - runtime_cv$tic),
                (runtime_npu$toc - runtime_npu$tic),
                (runtime_cat$toc - runtime_cat$tic),
                (runtime_su$toc - runtime_su$tic),
                (runtime_samp_h250$toc - runtime_samp_h250$tic),
                (runtime_samp_h500$toc - runtime_samp_h500$tic))
  return(runtimes)
}

# run in parallel
tictoc::tic()
future::plan(multisession, workers = 6)
run_tests(sim_reps, d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters)
future::plan(sequential)
runtime_test <- tictoc::toc()

# calc runtime for 1000 simulations
(runtime_test$toc - runtime_test$tic) / (60 * sim_reps) * 1000 / 60







