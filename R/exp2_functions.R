#' function to  run experiment 2 tests in parallel
#' 
#' @param sim_reps number of simulation replicates desired
#' @param d population exponential decay parameters
#' @param pu number of population units (e.g., number of schools)
#' @param pc number of population categories (e.g., ages or lengths)
#' @param pu_cv CV in mean category within a population unit (e.g., spread in ages around mean age within a given school) 
#' @param su_num total number of sampling units
#' @param su_samp vector of sampling unit sample sizes
#' @param p_su_samp vector of probabilities for sampling unit sample sizes
#' @param iters number of iterations that sample pop'n
#' @param cov_strc logistic-normal covariance structure options ("iid" and/or "1DAR1")
#' 
#' @return runtime for experiement 2 tests
#' 
#' @export
#' 
run_exp2_tests <- function(sim_reps, d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters, cov_strc){
  require(future)
  
  runtime_base %<-% test_base(sim_reps, d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters, cov_strc) %seed% TRUE
  runtime_CV %<-% test_CV(sim_reps, d, pu, pc, su_num, su_samp, p_su_samp, iters, cov_strc) %seed% TRUE
  runtime_C %<-% test_C(sim_reps, d, pu, pu_cv, su_num, su_samp, p_su_samp, iters, cov_strc) %seed% TRUE
  runtime_SU %<-% test_SU(sim_reps, d, pu, pc, pu_cv, su_samp, p_su_samp, iters, cov_strc) %seed% TRUE
  runtime_nSU %<-% test_nSU(sim_reps, d, pu, pc, pu_cv, su_num, iters, cov_strc) %seed% TRUE
  runtimes <- c((runtime_base$toc - runtime_base$tic),
                (runtime_CV$toc - runtime_CV$tic),
                (runtime_C$toc - runtime_C$tic),
                (runtime_SU$toc - runtime_SU$tic),
                (runtime_nSU$toc - runtime_nSU$tic))
}

#' function to test expansion weighting, selectivity, & pop'n structure
#' 
#' @param sim_reps number of simulation replicates desired
#' @param d population exponential decay parameters
#' @param pu number of population units (e.g., number of schools)
#' @param pc number of population categories (e.g., ages or lengths)
#' @param pu_cv CV in mean category within a population unit (e.g., spread in ages around mean age within a given school) 
#' @param su_num total number of sampling units
#' @param su_samp vector of sampling unit sample sizes
#' @param p_su_samp vector of probabilities for sampling unit sample sizes
#' @param iters number of iterations that sample pop'n
#' @param cov_strc logistic-normal covariance structure options ("iid" and/or "1DAR1")
#' 
#' @return runtime for test
#' 
#' @export
#' 
test_base <- function(sim_reps, d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters, cov_strc){
  #start timer
  tictoc::tic()
  # run simulation
  rr_base <- purrr::map(1:sim_reps, ~rep_sim(d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters, cov_strc))
  # save & plot results
  plot_base(rr_base)
  # end timer
  runtime <- tictoc::toc()
}

#' function to test pop'n unit structure (spread around mean category determined by CV)
#' 
#' @param sim_reps number of simulation replicates desired
#' @param d population exponential decay parameters
#' @param pu number of population units (e.g., number of schools)
#' @param pc number of population categories (e.g., ages or lengths)
#' @param su_num total number of sampling units
#' @param su_samp vector of sampling unit sample sizes
#' @param p_su_samp vector of probabilities for sampling unit sample sizes
#' @param iters number of iterations that sample pop'n
#' @param cov_strc logistic-normal covariance structure options ("iid" and/or "1DAR1")
#' 
#' @return runtime for test
#' 
#' @export
#' 
test_CV <- function(sim_reps, d, pu, pc, su_num, su_samp, p_su_samp, iters, cov_strc){
  #start timer
  tictoc::tic()
  # set levels of cv
  pu_cv_test <- c(0.1, 0.25, 1, 100)
  # run simulation
  rr_cv <- purrr::map(1:sim_reps, ~purrr::map(1:length(pu_cv_test), ~rep_sim(d, pu, pc, pu_cv = pu_cv_test[.], su_num, su_samp, p_su_samp, iters, cov_strc)))
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

#' function to test number of categories (i.e., longevity, growth)
#' 
#' @param sim_reps number of simulation replicates desired
#' @param d population exponential decay parameters
#' @param pu number of population units (e.g., number of schools)
#' @param pu_cv CV in mean category within a population unit (e.g., spread in ages around mean age within a given school) 
#' @param su_num total number of sampling units
#' @param su_samp vector of sampling unit sample sizes
#' @param p_su_samp vector of probabilities for sampling unit sample sizes
#' @param iters number of iterations that sample pop'n
#' @param cov_strc logistic-normal covariance structure options ("iid" and/or "1DAR1")
#' 
#' @return runtime for test
#' 
#' @export
#' 
test_C <- function(sim_reps, d, pu, pu_cv, su_num, su_samp, p_su_samp, iters, cov_strc){
  #start timer
  tictoc::tic()
  # set numbers of categories
  cat_test <- c(10, 15, 25, 50)
  # run simulation
  rr_cat <- purrr::map(1:sim_reps, ~purrr::map(1:length(cat_test), ~rep_sim(d, pu, pc = cat_test[.], pu_cv, su_num, su_samp, p_su_samp, iters, cov_strc)))
  # save & plot results
  plot_sim(rr = rr_cat, 
           plot_name = 'Ncat', 
           test_vec = cat_test,
           test_name = "C",
           test_lab = 'Number of categories within the population')
  # end timer
  runtime <- tictoc::toc()
  
  return(runtime)
}

#' function to  test number of sampling units (i.e., number of hauls)
#' 
#' @param sim_reps number of simulation replicates desired
#' @param d population exponential decay parameters
#' @param pu number of population units (e.g., number of schools)
#' @param pc number of population categories (e.g., ages or lengths)
#' @param pu_cv CV in mean category within a population unit (e.g., spread in ages around mean age within a given school) 
#' @param su_samp vector of sampling unit sample sizes
#' @param p_su_samp vector of probabilities for sampling unit sample sizes
#' @param iters number of iterations that sample pop'n
#' @param cov_strc logistic-normal covariance structure options ("iid" and/or "1DAR1")
#' 
#' @return runtime for test
#' 
#' @export
#' 
test_SU <- function(sim_reps, d, pu, pc, pu_cv, su_samp, p_su_samp, iters, cov_strc){
  #start timer
  tictoc::tic()
  # set number of sampling units
  su_test <- c(100, 250, 500, 1000)
  # run simulation
  rr_su <- purrr::map(1:sim_reps, ~purrr::map(1:length(su_test), ~rep_sim(d, pu, pc, pu_cv, su_num = su_test[.], su_samp, p_su_samp, iters, cov_strc)))
  # save & plot results
  plot_sim(rr = rr_su, 
           plot_name = 'Nsu', 
           test_vec = su_test, 
           test_name = "S",
           test_lab = 'Number of sampling units')
  # end timer
  runtime <- tictoc::toc()
  
  return(runtime)
}

#' function to  test sample size within sampling units
#' 
#' @param sim_reps number of simulation replicates desired
#' @param d population exponential decay parameters
#' @param pu number of population units (e.g., number of schools)
#' @param pc number of population categories (e.g., ages or lengths)
#' @param pu_cv CV in mean category within a population unit (e.g., spread in ages around mean age within a given school) 
#' @param su_num total number of sampling units
#' @param iters number of iterations that sample pop'n
#' @param cov_strc logistic-normal covariance structure options ("iid" and/or "1DAR1")
#' 
#' @return runtime for test
#' 
#' @export
#' 
test_nSU <- function(sim_reps, d, pu, pc, pu_cv, su_num, iters, cov_strc){
  #start timer
  tictoc::tic()
  # set sample size within sampling units
  samp_test <- c(10, 20, 50, 100)
  # run simulation
  rr_samp <- purrr::map(1:sim_reps, ~purrr::map(1:length(samp_test), ~rep_sim(d, pu, pc, pu_cv, su_num, su_samp = c(samp_test[.], 10), p_su_samp = c(1, 0), iters, cov_strc)))
  # save & plot results
  plot_sim(rr = rr_samp, 
           plot_name = 'nsu_250', 
           test_vec = samp_test, 
           test_name = "n",
           test_lab = 'Number of samples within a sampling unit')
  # end timer
  runtime <- tictoc::toc()
  
  return(runtime)
}

#' function to replicate simulation of sampling population comprised of subunits with different compositions
#'
#' @param d population exponential decay parameters
#' @param pu number of population units (e.g., number of schools)
#' @param pc number of population categories (e.g., ages or lengths)
#' @param pu_cv CV in mean category within a population unit (e.g., spread in ages around mean age within a given school) 
#' @param su_num total number of sampling units
#' @param su_samp vector of sampling unit sample sizes
#' @param p_su_samp vector of probabilities for sampling unit sample sizes
#' @param iters number of iterations that sample pop'n
#' @param cov_strc logistic-normal covariance structure options ("iid" and/or "1DAR1")
#' 
#' @return list of input sample size (by expansion complexity and selectivity form) and nominal sample size (nss)
#' 
#' @export
#' 
rep_sim <- function(d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters, cov_strc){
  
  # get simulated pop'n
  sim_popn <- get_popn(d, pu, pc, pu_cv, plot = FALSE)
  
  # run sim loop
  rr_sim <- purrr::map(1:iters, ~sim_comp(su_num, sim_popn, su_samp, p_su_samp))
  
  # estimate statistics & join with other results
  stats <- est_stats(rr_sim, sim_popn, cov_strc) %>% 
    tidytable::left_join(do.call(mapply, c(list, rr_sim, SIMPLIFY = FALSE))$nss %>% 
                           tidytable::map_df(., ~as.data.frame(.x), .id = "sim") %>% 
                           tidytable::summarise(mean_nss = mean(nss),
                                                .by = selex_type)) %>% 
    tidytable::left_join(do.call(mapply, c(list, rr_sim, SIMPLIFY = FALSE))$popn_strctr[[1]]) 

  # return results
  list(stats = stats)
  
}

#' function to save & plot results for expansion weighting, selectivity, & pop'n structure
#' 
#' @param rr list of results from simulation
#' 
#' @return saves dataframe results to 'output' folder and plots to 'figs' folder
#' 
#' @export
#' 
plot_base <- function(rr){
  
  # unlist results
  do.call(mapply, c(list, rr, SIMPLIFY = FALSE))$stats %>% 
    tidytable::map_df(., ~as.data.frame(.x), .id = "rep") -> res
  
  # save results
  saveRDS(res,
          file = here::here('output', 'exp2_base.rds'))

  # plot results by population structure
  plot_dat <- res %>% 
    tidytable::select(-c(sigma_1DAR1, rho_1DAR1, ess_DM)) %>% 
    tidytable::pivot_longer(cols = c(iss, sigma_iid, theta), names_to = 'param', values_to = 'stat') %>% 
    tidytable::mutate(param = case_when(param == 'iss' ~ 'Mult(ISS)',
                                       param == 'sigma_iid' ~ 'LogisticN(\u03C3)',
                                       param == 'theta' ~ 'DM(\u03B8)')) %>% 
    tidytable::mutate(popn_strctr = factor(popn_strctr, levels = c('recruitment pulse', 'multimodal', 'unimodal')),
                      param = factor(param, levels = c('Mult(ISS)', 'LogisticN(\u03C3)', 'DM(\u03B8)')))

  plot <- ggplot(data = plot_dat, aes(x = selex_type, y = stat, fill = comp_type)) +
    geom_boxplot(aes(fill = comp_type), position = position_dodge(0.4), width = 0.5, alpha = 0.7, outliers = FALSE) +
    facet_grid(param ~ popn_strctr, scales = 'free_y') +
    scale_fill_manual(values = c(scico::scico(3, palette = 'roma')[1], scico::scico(3, palette = 'roma')[2])) +
    theme_bw() +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    xlab('Selectivity shape') +
    ylab('Composition pdf statistic') +
    labs(fill = 'Expansion type') +
    theme(legend.position = "top")

  ggsave(filename = "exp2_base.png",
         plot = plot,
         path = here::here("figs"),
         width = 6.5,
         height = 5,
         units = "in")
  
}

#' function to save & plot results for simulations of tested effects
#' 
#' @param rr list of results
#' @param plot_name character string of name to save plot as
#' @param test_vec vector of tested values for simulation
#' @param test_name character string to name the test being conducted in plots
#' @param test_lab character string to name the x-axis in plot
#' @param plot_nss boolean, to plot the sample size as a horizontal line (default = FALSE)
#' @param fact_perc boolean, whether the facet factor for plotting should be converted to percent or not (default = FALSE)
#' 
#' @return saves dataframe results to 'output' folder and plots to 'figs' folder
#' 
#' @export
#' 
plot_sim <- function(rr, plot_name, test_vec, test_name, test_lab, plot_nss = FALSE, fact_perc = FALSE){

  # unlist results
  if(isTRUE(fact_perc)){
    res <- purrr::map(1:length(rr), ~(do.call(mapply, c(list, rr[[.]], SIMPLIFY = FALSE))$stats %>% 
                                        tidytable::map_df(., ~as.data.frame(.x), .id = "test") %>% 
                                        tidytable::mutate(test = dplyr::case_when(test == 1 ~ scales::percent(test_vec[1]),
                                                                                  test == 2 ~ scales::percent(test_vec[2]),
                                                                                  test == 3 ~ scales::percent(test_vec[3]),
                                                                                  test == 4 ~ scales::percent(test_vec[4]))))) %>% 
      tidytable::map_df(., ~as.data.frame(.x), .id = "rep") %>% 
      tidytable::mutate(facet = factor(paste0(test_name, " = ", test), levels = paste0(test_name, " = ", scales::percent(test_vec))))
  } else{
    res <- purrr::map(1:length(rr), ~(do.call(mapply, c(list, rr[[.]], SIMPLIFY = FALSE))$stats %>% 
                                        tidytable::map_df(., ~as.data.frame(.x), .id = "test") %>% 
                                        tidytable::mutate(test = dplyr::case_when(test == 1 ~ test_vec[1],
                                                                                  test == 2 ~ test_vec[2],
                                                                                  test == 3 ~ test_vec[3],
                                                                                  test == 4 ~ test_vec[4])))) %>% 
      tidytable::map_df(., ~as.data.frame(.x), .id = "rep") %>% 
      tidytable::mutate(facet = factor(paste0(test_name, " = ", test), levels = paste0(test_name, " = ", test_vec)))
  }
  
  # save results
  saveRDS(res,
          file = here::here('output', paste0('exp2_', plot_name, '.rds')))

  # plot source simulated results by population structure
  plot_dat <- res %>% 
    tidytable::select(-c(sigma_1DAR1, rho_1DAR1, ess_DM)) %>% 
    tidytable::pivot_longer(cols = c(iss, sigma_iid, theta), names_to = 'param', values_to = 'stat') %>% 
    tidytable::mutate(param = case_when(param == 'iss' ~ 'Mult(ISS)',
                                        param == 'sigma_iid' ~ 'LogisticN(\u03C3)',
                                        param == 'theta' ~ 'DM(\u03B8)')) %>% 
    tidytable::left_join(res %>% 
                           tidytable::summarise(nss = mean(mean_nss)) %>% 
                           tidytable::mutate(param = 'Mult(ISS)')) %>% 
    tidytable::mutate(popn_strctr = factor(popn_strctr, levels = c('recruitment pulse', 'multimodal', 'unimodal')),
                      param = factor(param, levels = c('Mult(ISS)', 'LogisticN(\u03C3)', 'DM(\u03B8)')))

 if(isTRUE(plot_nss)){
   suppressWarnings(plot <- ggplot(data = plot_dat, aes(x = facet, y = stat, fill = comp_type)) +
                      geom_boxplot(aes(fill = comp_type), position = position_dodge(0.4), width = 0.5, alpha = 0.7, outliers = FALSE) +
                      geom_hline(data = plot_dat,
                                 aes(yintercept = nss),
                                 colour = scico::scico(3, palette = 'roma')[3]) +
                      facet_grid(param ~ popn_strctr, scales = 'free_y') +
                      scale_fill_manual(values = c(scico::scico(3, palette = 'roma')[1], scico::scico(3, palette = 'roma')[2])) +
                      theme_bw() +
                      xlab(test_lab) +
                      scale_x_discrete(guide = guide_axis(angle = 45)) +
                      ylab('Composition pdf statistic') +
                      labs(fill = 'Expansion type') +
                      theme(legend.position = "top"))
  } else{
    plot <- ggplot(data = plot_dat, aes(x = facet, y = stat, fill = comp_type)) +
      geom_boxplot(aes(fill = comp_type), position = position_dodge(0.4), width = 0.5, alpha = 0.7, outliers = FALSE) +
      facet_grid(param ~ popn_strctr, scales = 'free_y') +
      scale_fill_manual(values = c(scico::scico(3, palette = 'roma')[1], scico::scico(3, palette = 'roma')[2])) +
      theme_bw() +
      xlab(test_lab) +
      scale_x_discrete(guide = guide_axis(angle = 45)) +
      ylab('Composition pdf statistic') +
      labs(fill = 'Expansion type') +
      theme(legend.position = "top")
    
  }
  
  ggsave(filename = paste0('exp2_', plot_name, '.png'),
         plot = plot,
         path = here::here("figs"),
         width = 6.5,
         height = 5,
         units = "in")
}

