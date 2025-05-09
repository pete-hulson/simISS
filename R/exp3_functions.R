#' function to run bootstrap tests in parallel
#' 
#' @param bs_iters number of iterations to run bootstrap test over
#' @param d population exponential decay parameters
#' @param pu number of population units (e.g., number of schools)
#' @param pc number of population categories (e.g., ages or lengths)
#' @param pu_cv CV in mean category within a population unit (e.g., spread in ages around mean age within a given school) 
#' @param su_num total number of sampling units
#' @param su_samp vector of sampling unit sample sizes
#' @param p_su_samp vector of probabilities for sampling unit sample sizes
#' @param iters number of bootstrap iterations within a sampling event
#' @param numCore number of cores available for parallel processing (function developed for 10 or 7 cores)
#' 
#' @return runtimes for bootstrap tests
#' 
#' @export
#' 
run_bs_test <- function(bs_iters, d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters, numCore){
  
  # define parallel function runs
  if(numCore > 10){ # using 10 cores
    run1 %<-% test_bs(bs_iters, d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters) %seed% TRUE
    run2 %<-% test_bs(bs_iters, d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters) %seed% TRUE
    run3 %<-% test_bs(bs_iters, d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters) %seed% TRUE
    run4 %<-% test_bs(bs_iters, d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters) %seed% TRUE
    run5 %<-% test_bs(bs_iters, d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters) %seed% TRUE
    run6 %<-% test_bs(bs_iters, d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters) %seed% TRUE
    run7 %<-% test_bs(bs_iters, d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters) %seed% TRUE
    run8 %<-% test_bs(bs_iters, d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters) %seed% TRUE
    run9 %<-% test_bs(bs_iters, d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters) %seed% TRUE
    run10 %<-% test_bs(bs_iters, d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters) %seed% TRUE
  } else{ # using 7 cores
    run1 %<-% test_bs(bs_iters, d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters) %seed% TRUE
    run2 %<-% test_bs(bs_iters, d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters) %seed% TRUE
    run3 %<-% test_bs(bs_iters, d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters) %seed% TRUE
    run4 %<-% test_bs(bs_iters, d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters) %seed% TRUE
    run5 %<-% test_bs(bs_iters, d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters) %seed% TRUE
    run6 %<-% test_bs(bs_iters, d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters) %seed% TRUE
    run7 %<-% test_bs(bs_iters, d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters) %seed% TRUE
  }
  
  # combine results
  if(numCore > 10){ # using 10 cores
    true_stats <- run1$true_stats %>% 
      tidytable::mutate(run = 1) %>% 
      tidytable::bind_rows(run2$true_stats %>% 
                             tidytable::mutate(run = 2)) %>% 
      tidytable::bind_rows(run3$true_stats %>% 
                             tidytable::mutate(run = 3)) %>% 
      tidytable::bind_rows(run4$true_stats %>% 
                             tidytable::mutate(run = 4)) %>% 
      tidytable::bind_rows(run5$true_stats %>% 
                             tidytable::mutate(run = 5)) %>% 
      tidytable::bind_rows(run6$true_stats %>% 
                             tidytable::mutate(run = 6)) %>% 
      tidytable::bind_rows(run7$true_stats %>% 
                             tidytable::mutate(run = 7)) %>% 
      tidytable::bind_rows(run8$true_stats %>% 
                             tidytable::mutate(run = 8)) %>% 
      tidytable::bind_rows(run9$true_stats %>% 
                             tidytable::mutate(run = 9)) %>% 
      tidytable::bind_rows(run10$true_stats %>% 
                             tidytable::mutate(run = 10))
    
    bs_stats <- run1$bs_stats %>% 
      tidytable::mutate(run = 1) %>% 
      tidytable::bind_rows(run2$bs_stats %>% 
                             tidytable::mutate(run = 2)) %>% 
      tidytable::bind_rows(run3$bs_stats %>% 
                             tidytable::mutate(run = 3)) %>% 
      tidytable::bind_rows(run4$bs_stats %>% 
                             tidytable::mutate(run = 4)) %>% 
      tidytable::bind_rows(run5$bs_stats %>% 
                             tidytable::mutate(run = 5)) %>% 
      tidytable::bind_rows(run6$bs_stats %>% 
                             tidytable::mutate(run = 6)) %>% 
      tidytable::bind_rows(run7$bs_stats %>% 
                             tidytable::mutate(run = 7)) %>% 
      tidytable::bind_rows(run8$bs_stats %>% 
                             tidytable::mutate(run = 8)) %>% 
      tidytable::bind_rows(run9$bs_stats %>% 
                             tidytable::mutate(run = 9)) %>% 
      tidytable::bind_rows(run10$bs_stats %>% 
                             tidytable::mutate(run = 10))
    
    popn_strctr <- run1$popn_strctr %>% 
      tidytable::mutate(run = 1) %>% 
      tidytable::bind_rows(run2$popn_strctr %>% 
                             tidytable::mutate(run = 2)) %>% 
      tidytable::bind_rows(run3$popn_strctr %>% 
                             tidytable::mutate(run = 3)) %>% 
      tidytable::bind_rows(run4$popn_strctr %>% 
                             tidytable::mutate(run = 4)) %>% 
      tidytable::bind_rows(run5$popn_strctr %>% 
                             tidytable::mutate(run = 5)) %>% 
      tidytable::bind_rows(run6$popn_strctr %>% 
                             tidytable::mutate(run = 6)) %>% 
      tidytable::bind_rows(run7$popn_strctr %>% 
                             tidytable::mutate(run = 7)) %>% 
      tidytable::bind_rows(run8$popn_strctr %>% 
                             tidytable::mutate(run = 8)) %>% 
      tidytable::bind_rows(run9$popn_strctr %>% 
                             tidytable::mutate(run = 9)) %>% 
      tidytable::bind_rows(run10$popn_strctr %>% 
                             tidytable::mutate(run = 10))
    
    res <- list(true_stats = true_stats, bs_stats = bs_stats, popn_strctr = popn_strctr)
    
    # save results
    saveRDS(res,
            file = here::here('output', paste0('exp3_bs.rds')))
    
  } else{ # using 7 cores
    true_stats <- run1$true_stats %>% 
      tidytable::mutate(run = 1) %>% 
      tidytable::bind_rows(run2$true_stats %>% 
                             tidytable::mutate(run = 2)) %>% 
      tidytable::bind_rows(run3$true_stats %>% 
                             tidytable::mutate(run = 3)) %>% 
      tidytable::bind_rows(run4$true_stats %>% 
                             tidytable::mutate(run = 4)) %>% 
      tidytable::bind_rows(run5$true_stats %>% 
                             tidytable::mutate(run = 5)) %>% 
      tidytable::bind_rows(run6$true_stats %>% 
                             tidytable::mutate(run = 6)) %>% 
      tidytable::bind_rows(run7$true_stats %>% 
                             tidytable::mutate(run = 7))
    
    bs_stats <- run1$bs_stats %>% 
      tidytable::mutate(run = 1) %>% 
      tidytable::bind_rows(run2$bs_stats %>% 
                             tidytable::mutate(run = 2)) %>% 
      tidytable::bind_rows(run3$bs_stats %>% 
                             tidytable::mutate(run = 3)) %>% 
      tidytable::bind_rows(run4$bs_stats %>% 
                             tidytable::mutate(run = 4)) %>% 
      tidytable::bind_rows(run5$bs_stats %>% 
                             tidytable::mutate(run = 5)) %>% 
      tidytable::bind_rows(run6$bs_stats %>% 
                             tidytable::mutate(run = 6)) %>% 
      tidytable::bind_rows(run7$bs_stats %>% 
                             tidytable::mutate(run = 7))
    
    popn_strctr <- run1$popn_strctr %>% 
      tidytable::mutate(run = 1) %>% 
      tidytable::bind_rows(run2$popn_strctr %>% 
                             tidytable::mutate(run = 2)) %>% 
      tidytable::bind_rows(run3$popn_strctr %>% 
                             tidytable::mutate(run = 3)) %>% 
      tidytable::bind_rows(run4$popn_strctr %>% 
                             tidytable::mutate(run = 4)) %>% 
      tidytable::bind_rows(run5$popn_strctr %>% 
                             tidytable::mutate(run = 5)) %>% 
      tidytable::bind_rows(run6$popn_strctr %>% 
                             tidytable::mutate(run = 6)) %>% 
      tidytable::bind_rows(run7$popn_strctr %>% 
                             tidytable::mutate(run = 7))
    
    res <- list(true_stats = true_stats, bs_stats = bs_stats, popn_strctr = popn_strctr)
    
    # save results
    saveRDS(res,
            file = here::here('output', paste0('exp3_bs.rds')))
    
  }
  
  # plot results
  # box vs point plot
  dat_true <- res$true_stats %>%
    tidytable::left_join(res$popn_strctr)  %>%
    tidytable::summarise(iss = psych::harmonic.mean(rss, zero = FALSE),
                         sigma_iid = psych::harmonic.mean(sigma_iid, zero = FALSE),
                         theta = psych::harmonic.mean(theta, zero = FALSE),
                         .by = c(selex_type, comp_type, popn_strctr)) %>%
    tidytable::pivot_longer(cols = c(iss, sigma_iid, theta), names_to = 'param', values_to = 'stat') %>% 
    tidytable::mutate(param = case_when(param == 'iss' ~ 'Mult(ISS)',
                                        param == 'sigma_iid' ~ 'LogisticN(\u03C3)',
                                        param == 'theta' ~ 'DM(\u03B8)')) %>% 
    tidytable::mutate(popn_strctr = factor(popn_strctr, levels = c('recruitment pulse', 'multimodal', 'unimodal')),
                      param = factor(param, levels = c('Mult(ISS)', 'LogisticN(\u03C3)', 'DM(\u03B8)')))
  
  dat_bs <- res$bs_stats %>%
    tidytable::left_join(res$popn_strctr) %>%
    tidytable::pivot_longer(cols = c(iss_bs, sigma_iid_bs, theta_bs), names_to = 'param', values_to = 'stat') %>% 
    tidytable::mutate(param = case_when(param == 'iss_bs' ~ 'Mult(ISS)',
                                        param == 'sigma_iid_bs' ~ 'LogisticN(\u03C3)',
                                        param == 'theta_bs' ~ 'DM(\u03B8)')) %>% 
    tidytable::mutate(popn_strctr = factor(popn_strctr, levels = c('recruitment pulse', 'multimodal', 'unimodal')),
                      param = factor(param, levels = c('Mult(ISS)', 'LogisticN(\u03C3)', 'DM(\u03B8)')))
  
  bs_plot1 <- ggplot(data = dat_bs, aes(x = selex_type, y = stat, fill = comp_type)) +
    geom_boxplot(position = position_dodge(0.4), width = 0.5, alpha = 0.5) +
    geom_point(data = dat_true, position = position_dodge(0.4), shape = 24, size = 2) +
    facet_grid(param ~ popn_strctr, scales = 'free_y') +
    scale_fill_manual(values = c(scico::scico(3, palette = 'roma')[1], scico::scico(3, palette = 'roma')[2])) +
    theme_bw() +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    xlab('Selectivity shape') +
    ylab('Composition pdf statistic') +
    labs(fill = 'Expansion type') +
    theme(legend.position = "top")
  
  ggsave(filename = "exp3_bs1.png",
         plot = bs_plot1,
         path = here::here("figs"),
         width = 6.5,
         height = 5,
         units = "in")
  
  
  # point plot of means
  dat_all <- res$true_stats %>%
    tidytable::left_join(res$popn_strctr) %>%
    tidytable::pivot_longer(cols = c(rss, sigma_iid, theta), names_to = 'param', values_to = 'stat_true') %>% 
    tidytable::mutate(param = case_when(param == 'rss' ~ 'Mult(ISS)',
                                        param == 'sigma_iid' ~ 'LogisticN(\u03C3)',
                                        param == 'theta' ~ 'DM(\u03B8)')) %>% 
    tidytable::mutate(popn_strctr = factor(popn_strctr, levels = c('recruitment pulse', 'multimodal', 'unimodal')),
                      param = factor(param, levels = c('Mult(ISS)', 'LogisticN(\u03C3)', 'DM(\u03B8)'))) %>% 
    tidytable::summarise(stat_true = psych::harmonic.mean(stat_true, zero = FALSE),
                         .by = c(selex_type, comp_type, popn_strctr, param)) %>% 
    tidytable::left_join(res$bs_stats %>%
                           tidytable::left_join(res$popn_strctr) %>%
                           tidytable::pivot_longer(cols = c(iss_bs, sigma_iid_bs, theta_bs), names_to = 'param', values_to = 'stat_bs') %>% 
                           tidytable::mutate(param = case_when(param == 'iss_bs' ~ 'Mult(ISS)',
                                                               param == 'sigma_iid_bs' ~ 'LogisticN(\u03C3)',
                                                               param == 'theta_bs' ~ 'DM(\u03B8)')) %>% 
                           tidytable::mutate(popn_strctr = factor(popn_strctr, levels = c('recruitment pulse', 'multimodal', 'unimodal')),
                                             param = factor(param, levels = c('Mult(ISS)', 'LogisticN(\u03C3)', 'DM(\u03B8)'))) %>% 
                           tidytable::summarise(stat_bs = psych::harmonic.mean(stat_bs, zero = FALSE),
                                                .by = c(selex_type, comp_type, popn_strctr, param)))
  
  
  bs_plot2 <- ggplot(data = dat_all, aes(x = stat_true, y = stat_bs, color = comp_type, shape = popn_strctr)) +
    geom_point(size = 2) +
    facet_wrap( ~ param, scales = 'free', strip.position = 'right', ncol = 1) +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
    theme_bw() +
    scale_color_manual(values = c(scico::scico(3, palette = 'roma')[1], scico::scico(3, palette = 'roma')[2])) +
    xlab('True statistic') +
    ylab('Bootstrap estimated statistic') +
    labs(color = 'Expansion type',
         shape = 'Population structure')
  
  ggsave(filename = "exp3_bs2.png",
         plot = bs_plot2,
         path = here::here("figs"),
         width = 6.5,
         height = 5,
         units = "in")
  
  
  
}

#' function that tests bootstrap method in experiment 3
#' 
#' @param bs_iters number of iterations to run bootstrap test over
#' @param d population exponential decay parameters
#' @param pu number of population units (e.g., number of schools)
#' @param pc number of population categories (e.g., ages or lengths)
#' @param pu_cv CV in mean category within a population unit (e.g., spread in ages around mean age within a given school) 
#' @param su_num total number of sampling units
#' @param su_samp vector of sampling unit sample sizes
#' @param p_su_samp vector of probabilities for sampling unit sample sizes
#' @param iters number of bootstrap iterations within a sampling event
#' 
#' @return list of results including realized sample size (RSS) for generated pop'n, input sample size (ISS) for bootstrap replicates, and pop'n structure (popn_strctr)
#' 
#' @export
#' 
test_bs <- function(bs_iters, d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters){
  
  #start timer
  tictoc::tic()
  # run simulation
  rr_bs <- purrr::map(1:bs_iters, ~bs_sim(d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters))
  # end timer
  runtime <- tictoc::toc()
  
  # unlist results
  do.call(mapply, c(list, rr_bs, SIMPLIFY = FALSE))$true_stats %>% 
    tidytable::map_df(., ~as.data.frame(.x), .id = "sim") -> true_stats
  do.call(mapply, c(list, rr_bs, SIMPLIFY = FALSE))$bs_stats %>% 
    tidytable::map_df(., ~as.data.frame(.x), .id = "sim") -> bs_stats
  do.call(mapply, c(list, rr_bs, SIMPLIFY = FALSE))$popn_strctr %>% 
    tidytable::map_df(., ~as.data.frame(.x), .id = "sim") -> popn_strctr
  
  res <- list(true_stats = true_stats, bs_stats = bs_stats, popn_strctr = popn_strctr)
  
  res
}

#' function to bootstrap a single sampling event realization
#' 
#' @param samp_ev list of results from running the sim_comp() function (after running the get_popn() fcn)
#' 
#' @return list of bootstrap realized sample size and comp
#' 
#' @export
#' 
bs_samp_event <- function(samp_ev){
  
  # bootstrap the sample realization
  
  # get number of sampling units
  su_num <- length(unique(samp_ev$N_su$samp_event))
  
  # bootsrap realized sampling units (e.g., hauls)
  bs_su <- data.frame(samp_event = sample(1:su_num, su_num, replace = TRUE)) %>% 
    tidytable::mutate(id = .I)
  
  # bootstrap samples within sampling unit (e.g., ages/lengths)
  bs_n_su <- tidytable::expand_grid(bs_samp_event = 1:su_num,
                                    cat = 1:pc,
                                    selex_type = samp_ev$nss$selex_type) %>% 
    tidytable::left_join(bs_su %>% 
                           tidytable::left_join(samp_ev$n_su) %>% 
                           tidytable::select(bs_samp_event = id, selex_type, cat, samp) %>% 
                           tidytable::uncount(samp) %>% 
                           tidytable::mutate(cat = sample(cat, .N, replace = TRUE), .by = c(bs_samp_event, selex_type)) %>% 
                           tidytable::summarize(bs_samp = .N, .by = c(bs_samp_event, selex_type, cat))) %>% 
    tidytable::mutate(bs_samp = case_when(is.na(bs_samp) ~ 0, .default = bs_samp))
  
  # compute proportion of sample size across sampling units
  bs_n_su %>% 
    tidytable::summarise(bs_samp_se = sum(bs_samp), .by = c(bs_samp_event, selex_type)) %>% 
    tidytable::mutate(bs_p_samp_se = bs_samp_se / sum(bs_samp_se), .by = selex_type) -> bs_p_samp_su
  
  # calculate abundance proportions
  bs_su %>% 
    tidytable::left_join(samp_ev$N_su) %>% 
    tidytable::select(bs_samp_event = id, bs_N_su = N_su, selex_type) %>% 
    tidytable::mutate(bs_prop_N = bs_N_su / sum(bs_N_su), .by = selex_type) -> bs_N_su
  
  # compute bootstrap composition
  bs_n_su %>% 
    # join sampling event abundance
    tidytable::left_join(bs_N_su) %>% 
    # join sampling event sample size
    tidytable::left_join(bs_p_samp_su) %>% 
    tidytable::summarise(bs_samp_wtd = sum(bs_samp * bs_p_samp_se * bs_prop_N), # weight samples
                         bs_samp_unwtd = sum(bs_samp), # do not weight samples
                         .by = c(cat, selex_type)) %>% 
    # compute compositions
    tidytable::mutate(bs_samp_p_wtd = bs_samp_wtd / sum(bs_samp_wtd),
                      bs_samp_p_unwtd = bs_samp_unwtd / sum(bs_samp_unwtd), .by = selex_type) %>% 
    tidytable::select(cat, selex_type, bs_samp_p_wtd, bs_samp_p_unwtd) -> bs_comp
  
  # calculate the bootstrap realized sample size
  rss_bs <- bs_comp %>% 
    tidytable::left_join(samp_ev$comp) %>% 
    tidytable::summarise(bs_rss_wtd = sum(samp_p_wtd * (1- samp_p_wtd)) / sum((bs_samp_p_wtd - samp_p_wtd) ^ 2),
                         bs_rss_unwtd = sum(samp_p_unwtd * (1- samp_p_unwtd)) / sum((bs_samp_p_unwtd - samp_p_unwtd) ^ 2),
                         .by = selex_type)
  
  # get sample size
  nss_bs <- bs_n_su %>% 
    tidytable::summarise(nss = sum(bs_samp), .by = c(selex_type))
  
  
  list(rss_bs = rss_bs, bs_comp = bs_comp, nss_bs = nss_bs)
}

#' function to replicate bootstrap replication of sampling events
#'
#' @param d population exponential decay parameters
#' @param pu number of population units (e.g., number of schools)
#' @param pc number of population categories (e.g., ages or lengths)
#' @param pu_cv CV in mean category within a population unit (e.g., spread in ages around mean age within a given school) 
#' @param su_num total number of sampling units
#' @param su_samp vector of sampling unit sample sizes
#' @param p_su_samp vector of probabilities for sampling unit sample sizes
#' @param iters number of bootstrap iterations
#'
#' @return list of simulated pop'n realized sample size and input sample size from bootstrap of realized sampling evennt
#' 
#' @export
#' 
bs_sim <- function(d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters){
  
  # get a sample realization ----
  
  # generate the pop'n
  sim_popn <- get_popn(d, pu, pc, pu_cv, plot = FALSE)
  
  # generate the sampling event
  samp_ev <- sim_comp(su_num, sim_popn, su_samp, p_su_samp)
  
  # calculate the 'true' statistic for the sampling event
  
  # set up data list
  data <- list(exp = sim_popn$p_true %>% 
                 tidytable::select(-N_c), 
               obs = samp_ev$comp %>% 
                 tidytable::pivot_longer(cols = c('samp_p_wtd', 'samp_p_unwtd'), names_to = 'comp_type', values_to = 'p_obs') %>% 
                 tidytable::mutate(comp_type = case_when(comp_type == 'samp_p_wtd' ~ 'wtd',
                                                         .default = 'unwtd'),
                                   sim = 1),
               N = samp_ev$nss)
  
  # combinations of selectivity/composition expansion types tested
  combs <- tidytable::expand_grid(selex = unique(data$obs$selex_type), 
                                  comp = unique(data$obs$comp_type))
  
  
  # multinomial statistic (rss)
  mult_sim <- data$obs %>% 
    tidytable::left_join(data$exp) %>% 
    tidytable::summarise(rss = sum(p_true * (1- p_true)) / sum((p_obs - p_true) ^ 2),
                         .by = c(selex_type, comp_type))
  
  # logistic-normal statistics (sigma & rho)
  # estimate parameters
  rr_logistN <- purrr::map(1:dim(combs)[1],
                           ~est_logistic_normal(data = data, 
                                                selex_t = combs$selex[.],
                                                comp_t = combs$comp[.]))
  # get results
  logistN_sim <- do.call(mapply, c(list, rr_logistN, SIMPLIFY = FALSE))$res %>% 
    tidytable::map_df(., ~as.data.frame(.x), .id = "comb") %>% 
    tidytable::select(-comb)
  
  # dirichlet-multinomial statistic (theta)
  # estimate parameter
  rr_DM <- suppressWarnings(purrr::map(1:dim(combs)[1],
                                       ~est_dirmult(data, 
                                                    selex_t = combs$selex[.],
                                                    comp_t = combs$comp[.])))
  
  # unlist results
  DM_sim <- do.call(mapply, c(list, rr_DM, SIMPLIFY = FALSE))$res %>% 
    tidytable::map_df(., ~as.data.frame(.x), .id = "comb") %>% 
    tidytable::select(-comb)
  
  # put 'true' stats together
  true_stats <- mult_sim %>% 
    tidytable::left_join(logistN_sim) %>% 
    tidytable::left_join(DM_sim) %>% 
    tidytable::left_join(data$N)
  
  
  # run bootstrap of sampling event ----
  rr_bs <- purrr::map(1:iters, ~bs_samp_event(samp_ev))
  
  # set up data list
  data <- list(exp = samp_ev$comp %>% 
                 tidytable::pivot_longer(cols = c('samp_p_wtd', 'samp_p_unwtd'), names_to = 'comp_type', values_to = 'p_true') %>% 
                 tidytable::mutate(comp_type = case_when(comp_type == 'samp_p_wtd' ~ 'wtd',
                                                         .default = 'unwtd')), 
               obs = do.call(mapply, c(list, rr_bs, SIMPLIFY = FALSE))$bs_comp %>% 
                 tidytable::map_df(., ~as.data.frame(.x), .id = "sim") %>% 
                 tidytable::pivot_longer(cols = c('bs_samp_p_wtd', 'bs_samp_p_unwtd'), names_to = 'comp_type', values_to = 'p_obs') %>% 
                 tidytable::mutate(comp_type = case_when(comp_type == 'bs_samp_p_wtd' ~ 'wtd',
                                                         .default = 'unwtd')),
               N = do.call(mapply, c(list, rr_bs, SIMPLIFY = FALSE))$nss_bs %>% 
                 tidytable::map_df(., ~as.data.frame(.x), .id = "sim"))
  
  # get bootstrap statistics
  
  # multinomial statistic (iss)
  mult_bs <- do.call(mapply, c(list, rr_bs, SIMPLIFY = FALSE))$rss_bs %>% 
    tidytable::map_df(., ~as.data.frame(.x), .id = "sim") %>%     
    tidytable::summarise(bs_iss_wtd = psych::harmonic.mean(bs_rss_wtd, zero = FALSE),
                         bs_iss_unwtd = psych::harmonic.mean(bs_rss_unwtd, zero = FALSE),
                         .by = selex_type) %>% 
    tidytable::pivot_longer(cols = c('bs_iss_wtd', 'bs_iss_unwtd'), names_to = 'comp_type', values_to = 'iss_bs') %>% 
    tidytable::mutate(comp_type = case_when(comp_type == 'bs_iss_wtd' ~ 'wtd',
                                            .default = 'unwtd'))
  
  # logistic-normal statistics (sigma & rho)
  # estimate parameters
  rr_logistN <- purrr::map(1:dim(combs)[1],
                       ~est_logistic_normal(data = list(exp = data$exp[comp_type == combs$comp[.]],
                                                        obs = data$obs,
                                                        N = data$N), 
                                            selex_t = combs$selex[.],
                                            comp_t = combs$comp[.]))
  # get results
  logistN_bs <- do.call(mapply, c(list, rr_logistN, SIMPLIFY = FALSE))$res %>% 
    tidytable::map_df(., ~as.data.frame(.x), .id = "comb") %>% 
    tidytable::select(-comb)
    tidytable::rename(sigma_iid_bs = sigma_iid, sigma_1DAR1_bs = sigma_1DAR1, rho_1DAR1_bs = rho_1DAR1)
  
  # dirichlet-multinomial statistic (theta)
  # estimate parameter
  rr_DM <- suppressWarnings(purrr::map(1:dim(combs)[1],
                                       ~est_dirmult(data = list(exp = data$exp[comp_type == combs$comp[.]],
                                                                obs = data$obs,
                                                                N = data$N), 
                                                    selex_t = combs$selex[.],
                                                    comp_t = combs$comp[.])))
  # unlist results
  DM_bs <- do.call(mapply, c(list, rr_DM, SIMPLIFY = FALSE))$res %>% 
    tidytable::map_df(., ~as.data.frame(.x), .id = "comb") %>% 
    tidytable::select(selex_type, comp_type, theta_bs = theta, ess_DM_bs = ess_DM)
  
  # put bootstrap stats together
  bs_stats <- mult_bs %>% 
    tidytable::left_join(logistN_bs) %>% 
    tidytable::left_join(DM_bs)
  
  
  # output
  list(true_stats = true_stats, bs_stats = bs_stats, popn_strctr = samp_ev$popn_strctr)
}

