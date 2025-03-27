#' function to set up a population comprised of subunits with different compositions
#' 
#' @param d population exponential decay parameters
#' @param pu number of population units (e.g., number of schools)
#' @param pc number of population categories (e.g., ages or lengths)
#' @param pu_cv CV in mean category within a population unit (e.g., spread in ages around mean age within a given school)
#' @param replace_pc boolean, whether to sample mean pop'n category for pop'n unit with replacement or not (default = TRUE)
#' @param plot boolean, whether to output plot of generated pop'n (default = TRUE)
#' @param plot_name if plot desired, add name to write plot as (default = 'gen_popn')
#' 
#' @return List of population unit compositions (p_pu), relative abundance of population units (p_popn), 
#' and the combined population composition (p_true)
#' 
#' @export
#' 
get_popn <- function(d, pu, pc, pu_cv, replace_pc = TRUE, plot = TRUE, plot_name = 'gen_popn'){
  
  # step 1: set up pop'n with exponential decay adjusted by selex ----
  
  # uniform selex
  popn_u <- data.frame(cat = 1:pc, popn = exp(-(1:pc - 1) * d), selex_type = 'uniform', selex = 1)
  # asympotic selex
  # generate category at 50% selex (restricted to first 25% of categories)
  c_50 <- stats::runif(1, 1, 0.25 * pc)
  # generate slope (so that not knife edged)
  delta <- runif(1, 0.5, 2)
  # compute selectivity
  selex_a <- data.frame(cat = 1:pc, selex = 1 / (1 + exp(-(1:pc - c_50) / delta)))
  # adjust pop'n
  popn_a <- popn_u %>% 
    tidytable::select(-selex) %>%  
    tidytable::left_join(selex_a) %>% 
    tidytable::mutate(popn = popn * selex,
                      selex_type = 'asymptotic')
  # dome-shaped selex
  # generate category at max selex (restricted to 2nd quantile of categories)
  c_max <- stats::runif(1, 0.25 * pc, 0.5 * pc)
  # generate slope (so that not knife edged)
  delta <- runif(1, 2, 5)
  # calculate power parameter
  p <- 0.5 * (sqrt(c_max ^ 2 + (4 * delta ^ 2)) - c_max)
  # compute selectivity
  selex_d <- data.frame(cat = 1:pc, selex = (1:pc / c_max)^(c_max / p) * exp( (c_max - 1:pc) / p))
  # adjust pop'n
  popn_d <- popn_u %>% 
    tidytable::select(-selex) %>% 
    tidytable::left_join(selex_d) %>% 
    tidytable::mutate(popn = popn * selex,
                      selex_type = 'dome-shaped')
  # put all together
  popn <- popn_u %>% 
    tidytable::bind_rows(popn_a) %>% 
    tidytable::bind_rows(popn_d)
  
  
  # step 2: set up pop'n units ----
  
  # get mean category for pop'n unit (e.g., mean age of school)
  mu_cat <- sample(1:pc, pu, replace = replace_pc)
  
  # set up pop'n unit composition with normal distribution
  ptwid_pu <- purrr::map(1:pu, ~data.frame(cat = 1:pc, p_dist = dnorm(1:pc, mean = mu_cat[.], sd = mu_cat[.] * pu_cv) / sum(dnorm(1:pc, mean = mu_cat[.], sd = mu_cat[.] * pu_cv)))) %>% 
    tidytable::map_df(., ~as.data.frame(.x), .id = "popn_unit")
  
  # set up pop'n units taking into account seelx and pop'n size
  p_cpu <- ptwid_pu %>% 
    tidytable::left_join(popn) %>% 
    tidytable::mutate(N_cpu = p_dist * popn,
                      p_cpu = N_cpu / sum(N_cpu), .by = c(popn_unit, selex_type)) %>% 
    tidytable::select(popn_unit, cat, N_cpu, p_cpu, selex_type)
  
  # get pop'n unit relative proportions to use for resampling pop'n units
  p_pu <- ptwid_pu %>% 
    tidytable::left_join(popn) %>% 
    tidytable::summarise(N_pu = sum(p_dist * popn), .by = c(popn_unit, selex_type)) %>% 
    tidytable::mutate(p_pu = N_pu / sum(N_pu), .by = selex_type)
  
  
  # step 3: get 'true' pop'n composition ----
  
  p_true <- ptwid_pu %>% 
    tidytable::left_join(popn) %>% 
    tidytable::summarise(N_c = sum(p_dist * popn), .by = c(cat, selex_type)) %>% 
    tidytable::mutate(p_true = N_c / sum(N_c), .by = selex_type)
  
  
  # classify population unit structure as 'recruitment pulse', 'multimodal', or 'unimodal'
  p_true %>% 
    tidytable::select(-N_c) %>% 
    tidytable::left_join(p_true %>% 
                           tidytable::select(-N_c) %>% 
                           tidytable::mutate(cat = cat + 1) %>% 
                           tidytable::rename(p_true_1 = p_true)) %>% 
    tidytable::drop_na() %>% 
    tidytable::mutate(test = p_true - p_true_1) %>% 
    tidytable::select(cat, test, selex_type) -> .test
  
  p_true %>% 
    tidytable::filter(cat <= 2) %>% 
    tidytable::summarise(rec_test = sum(p_true), .by = selex_type) %>% 
    tidytable::left_join(.test %>% 
                           tidytable::left_join(.test %>% 
                                                  tidytable::mutate(cat = cat + 1) %>% 
                                                  tidytable::rename(test_1 = test)) %>% 
                           tidytable::drop_na() %>% 
                           tidytable::filter(test < 0 & test_1 > 0) %>% 
                           tidytable::summarise(mode_test = .N, .by = selex_type)) %>% 
    tidytable::mutate(mode_test = tidytable::replace_na(mode_test, 0)) -> test
  
  popn_strctr <- test %>% 
    tidytable::mutate(popn_strctr = case_when(rec_test > 0.45 ~ 'recruitment pulse',
                                              rec_test <= 0.45 ~ case_when(mode_test > 1 ~ 'multimodal', 
                                                                           .default = 'unimodal'))) %>% 
    tidytable::select(selex_type, popn_strctr)
  
  # if desired, plot generated pop'n
  if(isTRUE(plot)){
    
    p_cpu %>% 
      tidytable::rename(N = N_cpu) %>% 
      tidytable::select(-p_cpu) %>% 
      tidytable::bind_rows(p_true %>% 
                             tidytable::mutate(popn_unit = 'combined') %>% 
                             tidytable::select(popn_unit, cat, N = N_c, selex_type)) -> plot_dat
    
    popn_plot <- ggplot(data = plot_dat, aes(x = as.factor(cat), y = N, fill = popn_unit)) +
      geom_bar(stat = 'identity') +
      facet_grid(popn_unit ~ selex_type, scale = 'free_y') +
      theme_bw() +
      labs(fill = "Population unit") +
      xlab('Category') +
      ylab('Population unit abundance') +
      scico::scale_fill_scico_d(palette = 'roma') +
      scale_x_discrete(guide = guide_axis(angle = 90))
    
    ggsave(filename = paste0(plot_name, ".png"),
           plot = popn_plot,
           path = here::here("figs"),
           width = 6.5,
           height = 5,
           units = "in")
  }
  
  # list(p_pu = p_pu, p_popn = p_popn, p_true = p_true, popn_strctr = popn_strctr)
  list(p_cpu = p_cpu, p_pu = p_pu, p_true = p_true, popn_strctr = popn_strctr)
}

#' function to sample simulation population and compute composition
#' 
#' @param su_num total number of sampling units
#' @param sim_popn simulated population (obtained from get_popn() function)
#' @param su_samp vector of sampling unit sample sizes
#' @param p_su_samp vector of probabilities for sampling unit sample sizes
#' 
#' @return list of population level composition after sampling (comp) and total sample size (nss)
#' 
#' @export
#' 
sim_comp <- function(su_num, sim_popn, su_samp, p_su_samp){
  
  # step 1: generate sampling unit samples ----
  # determine which population unit is being sampled for each selex type (based on relative abundance)
  selex_type <- sim_popn$popn_strctr$selex_type
  num_selex <- length(selex_type)
  samp_pu <- NULL
  for(i in 1:num_selex){
    rand <- data.frame(samp_event = 1:su_num,
                       popn_unit = sample(x = unique(sim_popn$p_pu$popn_unit), 
                                          size = su_num, 
                                          replace = TRUE, 
                                          prob = sim_popn$p_pu$p_pu[sim_popn$p_pu$selex_type == selex_type[i]])) %>% 
      tidytable::mutate(selex_type = selex_type[i])
    
    samp_pu <- samp_pu %>% 
      tidytable::bind_rows(rand)
  }
  
  # define the sampling events and join sampled pop'n units
  tidytable::expand_grid(samp_event = 1:su_num, selex_type = sim_popn$popn_strctr$selex_type) %>% 
    tidytable::left_join(samp_pu) %>% 
    # join population unit composition
    tidytable::left_join(sim_popn$p_cpu) %>% 
    # generate samples within categories based on population unit composition across sampling events
    tidytable::mutate(samp = stats::rmultinom(1, 
                                              # generate sampling unit number of samples
                                              sample(x = su_samp, size = 1, replace = TRUE, prob = p_su_samp),
                                              p_cpu),
                      .by = c(samp_event, selex_type)) -> n_su
  
  # generate log-normal sampling unit abundance and calculate proportions
  samp_pu %>% 
    tidytable::mutate(N_su = exp(stats::rnorm(1, 0, 1)), .by = samp_event) %>% 
    tidytable::select(samp_event, N_su, selex_type) -> samp_N
  
  
  # step 2: compute weightings ----
  # compute proportion of sample size across sampling units
  n_su %>% 
    tidytable::summarise(samp_se = sum(samp), .by = c(samp_event, selex_type)) %>% 
    tidytable::mutate(p_samp_se = samp_se / sum(samp_se), .by = selex_type) -> p_samp_su
  
  # calculate proportions of sampling unit abundance
  samp_N %>% 
    tidytable::mutate(prop_N = N_su / sum(N_su), .by = selex_type) -> N_su 
  
  # compute total sample size
  n_su %>% 
    tidytable::summarise(nss = sum(samp), .by = selex_type) -> nss
  
  # step 3: compute generated comps ----
  # get generated composition results
  n_su %>% 
    # join sampling event abundance
    tidytable::left_join(N_su) %>% 
    # join sampling event sample size
    tidytable::left_join(p_samp_su) %>% 
    tidytable::summarise(samp_wtd = sum(samp * p_samp_se * prop_N), # weight samples
                         samp_unwtd = sum(samp), # do not weight samples
                         .by = c(cat, selex_type)) %>% 
    # compute compositions
    tidytable::mutate(samp_p_wtd = samp_wtd / sum(samp_wtd),
                      samp_p_unwtd = samp_unwtd / sum(samp_unwtd), .by = selex_type) %>% 
    tidytable::select(cat, samp_p_wtd, samp_p_unwtd, selex_type) -> comp
  
  
  # output
  list(comp = comp, nss = nss, popn_strctr = sim_popn$popn_strctr, n_su = n_su, N_su = N_su)
  
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
#' 
#' @return list of input sample size (by expansion complexity and selectivity form) and nominal sample size (nss)
#' 
#' @export
#' 
rep_sim <- function(d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters){
  
  # get simulated pop'n
  sim_popn <- get_popn(d, pu, pc, pu_cv, plot = FALSE)
  
  # run sim loop
  rr_sim <- purrr::map(1:iters, ~sim_comp(su_num, sim_popn, su_samp, p_su_samp))
  
  # unlist results
  do.call(mapply, c(list, rr_sim, SIMPLIFY = FALSE))$comp %>% 
    tidytable::map_df(., ~as.data.frame(.x), .id = "sim") -> res_sim
  
  # compute realized sample size
  res_sim %>% 
    tidytable::left_join(sim_popn$p_true) %>% 
    tidytable::summarise(rss_wtd = sum(p_true * (1- p_true)) / sum((samp_p_wtd - p_true) ^ 2),
                         rss_unwtd = sum(p_true * (1- p_true)) / sum((samp_p_unwtd - p_true) ^ 2),
                         .by = c(sim, selex_type)) -> rss_sim
  
  # compute input sample size & nominal sample size
  rss_sim %>% 
    tidytable::left_join(do.call(mapply, c(list, rr_sim, SIMPLIFY = FALSE))$nss %>% 
                           tidytable::map_df(., ~as.data.frame(.x), .id = "sim")) %>% 
    tidytable::summarise(iss_wtd = psych::harmonic.mean(rss_wtd, zero = FALSE),
                         iss_unwtd = psych::harmonic.mean(rss_unwtd, zero = FALSE),
                         mean_nss = mean(nss),
                         .by = selex_type) %>% 
    tidytable::left_join(do.call(mapply, c(list, rr_sim, SIMPLIFY = FALSE))$popn_strctr[[1]]) -> iss_sim
  
  list(iss_sim = iss_sim)
  
}

#' function to save & plot results for expansion and pop'n structure effects
#' 
#' @param rr list of results from simulation
#' 
#' @return saves dataframe results to 'output' folder and plots to 'figs' folder
#' 
#' @export
#' 
plot_exp <- function(rr){
  
  # unlist results
  do.call(mapply, c(list, rr, SIMPLIFY = FALSE))$iss_sim %>% 
    tidytable::map_df(., ~as.data.frame(.x), .id = "rep") -> res
  
  # save results
  saveRDS(res,
          file = here::here('output', 'exp2_test1.rds'))
  
  # plot results by population structure
  plot_dat <- res %>% 
    tidytable::rename(Wtd = iss_wtd, Unwtd = iss_unwtd) %>% 
    tidytable::pivot_longer(cols = c(Wtd, Unwtd, mean_nss)) %>% 
    tidytable::mutate(popn_strctr = factor(popn_strctr, levels = c('recruitment pulse', 'multimodal', 'unimodal', 'combined')))
  
  plot <- ggplot(data = plot_dat %>% 
                   tidytable::filter(name %in% c('Wtd', 'Unwtd')), aes(x = selex_type, y = value, fill = name)) +
    geom_boxplot(alpha = 0.7) +
    geom_hline(yintercept = as.numeric(plot_dat %>% 
                                         tidytable::filter(name %in% c('mean_nss')) %>% 
                                         tidytable::summarise(nss = mean(value))),
               linewidth = 1,
               colour = scico::scico(3, palette = 'roma')[3]) +
    facet_grid(name ~ popn_strctr) +
    scale_fill_manual(values = c(scico::scico(3, palette = 'roma')[1], scico::scico(3, palette = 'roma')[2])) +
    theme_bw() +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    xlab('Selectivity shape') +
    ylab('Input Sample Size (ISS)') +
    guides(fill = 'none')
  
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
#' @param plot_nss boolean, to plot the sample size as a horizontal line
#' @param fact_perc boolean, whether the facet factor for plotting should be converted to percent or not
#' 
#' @return saves dataframe results to 'output' folder and plots to 'figs' folder
#' 
#' @export
#' 
plot_sim <- function(rr, plot_name, test_vec, test_name, test_lab, plot_nss = FALSE, fact_perc = FALSE){
  
  # unlist results
  if(isTRUE(fact_perc)){
    res <- purrr::map(1:length(rr), ~(do.call(mapply, c(list, rr[[.]], SIMPLIFY = FALSE))$iss_sim %>% 
                                        tidytable::map_df(., ~as.data.frame(.x), .id = "test") %>% 
                                        tidytable::mutate(test = dplyr::case_when(test == 1 ~ scales::percent(test_vec[1]),
                                                                                  test == 2 ~ scales::percent(test_vec[2]),
                                                                                  test == 3 ~ scales::percent(test_vec[3]),
                                                                                  test == 4 ~ scales::percent(test_vec[4]))))) %>% 
      tidytable::map_df(., ~as.data.frame(.x), .id = "rep") %>% 
      tidytable::mutate(facet = factor(paste0(test_name, " = ", test), levels = paste0(test_name, " = ", scales::percent(test_vec)))) %>% 
      tidytable::pivot_longer(cols = c(iss_wtd, iss_unwtd, mean_nss))
  } else{
    res <- purrr::map(1:length(rr), ~(do.call(mapply, c(list, rr[[.]], SIMPLIFY = FALSE))$iss_sim %>% 
                                        tidytable::map_df(., ~as.data.frame(.x), .id = "test") %>% 
                                        tidytable::mutate(test = dplyr::case_when(test == 1 ~ test_vec[1],
                                                                                  test == 2 ~ test_vec[2],
                                                                                  test == 3 ~ test_vec[3],
                                                                                  test == 4 ~ test_vec[4])))) %>% 
      tidytable::map_df(., ~as.data.frame(.x), .id = "rep") %>% 
      tidytable::mutate(facet = factor(paste0(test_name, " = ", test), levels = paste0(test_name, " = ", test_vec))) %>% 
      tidytable::pivot_longer(cols = c(iss_wtd, iss_unwtd, mean_nss))
  }
  
  # save results
  saveRDS(res,
          file = here::here('output', paste0('exp2_', plot_name, '.rds')))
  
  # plot source simulated results by population structure
  .plot_dat <- res %>% 
    tidytable::mutate(popn_strctr = factor(popn_strctr, levels = c('recruitment pulse', 'multimodal', 'unimodal')))
  
  .plot_dat %>% 
    tidytable::filter(name != 'mean_nss') %>% 
    tidytable::left_join(.plot_dat %>% 
                           tidytable::filter(name == 'mean_nss') %>% 
                           tidytable::summarise(nss = mean(value), .by = facet)) %>% 
    tidytable::mutate(name = case_when(name == 'iss_wtd' ~ 'Wtd',
                                       name == 'iss_unwtd' ~ 'Unwtd')) -> plot_dat
  if(isTRUE(plot_nss)){
    plot <- ggplot(data = plot_dat, aes(x = facet, y = value, fill = name)) +
      geom_boxplot(alpha = 0.7) +
      geom_hline(data = plot_dat,
                 aes(yintercept = nss),
                 linewidth = 1,
                 colour = scico::scico(3, palette = 'roma')[3]) +
      facet_grid(name ~ popn_strctr, scales = 'free_y') +
      scale_fill_manual(values = c(scico::scico(3, palette = 'roma')[1], scico::scico(3, palette = 'roma')[2])) +
      theme_bw() +
      xlab(test_lab) +
      scale_x_discrete(guide = guide_axis(angle = 45)) +
      ylab('Input Sample Size (ISS)') +
      guides(fill = 'none')
  } else{
    plot <- ggplot(data = plot_dat, aes(x = facet, y = value, fill = name)) +
      geom_boxplot(alpha = 0.7) +
      facet_grid(name ~ popn_strctr, scales = 'free_y') +
      scale_fill_manual(values = c(scico::scico(3, palette = 'roma')[1], scico::scico(3, palette = 'roma')[2])) +
      theme_bw() +
      xlab(test_lab) +
      scale_x_discrete(guide = guide_axis(angle = 45)) +
      ylab('Input Sample Size (ISS)') +
      guides(fill = 'none')
    
  }
  
  ggsave(filename = paste0('exp2_', plot_name, '.png'),
         plot = plot,
         path = here::here("figs"),
         width = 6.5,
         height = 5,
         units = "in")
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
  
  list(rss_bs = rss_bs, bs_comp = bs_comp)
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
  
  # get a sample realization
  
  # generate the pop'n
  sim_popn <- get_popn(d, pu, pc, pu_cv, plot = FALSE)
  
  # generate the sampling event
  samp_ev <- sim_comp(su_num, sim_popn, su_samp, p_su_samp)
  
  # calculate the realized sample size for the sampling event
  samp_ev$comp %>% 
    tidytable::left_join(sim_popn$p_true) %>% 
    tidytable::summarise(rss_wtd = sum(p_true * (1- p_true)) / sum((samp_p_wtd - p_true) ^ 2),
                         rss_unwtd = sum(p_true * (1- p_true)) / sum((samp_p_unwtd - p_true) ^ 2),
                         .by = selex_type) -> rss_se
  
  # run bootstrap of sampling event
  rr_bs <- purrr::map(1:iters, ~bs_samp_event(samp_ev))
  
  # unlist results
  do.call(mapply, c(list, rr_bs, SIMPLIFY = FALSE))$rss_bs %>% 
    tidytable::map_df(., ~as.data.frame(.x), .id = "sim") -> res_bs
  
  # calculate bootstrap iss
  iss_bs <- res_bs %>%     
    tidytable::summarise(bs_iss_wtd = psych::harmonic.mean(bs_rss_wtd, zero = FALSE),
                         bs_iss_unwtd = psych::harmonic.mean(bs_rss_unwtd, zero = FALSE),
                         .by = selex_type)
  
  list(rss_se = rss_se, iss_bs = iss_bs, popn_strctr = samp_ev$popn_strctr)
}

#' function to test expansion weighting, selectivity, & pop'n structure
#' 
#' @export
#' 
test_base <- function(sim_reps, d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters){
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

#' function to test pop'n unit structure (spread around mean category determined by CV)
#' 
#' @export
#' 
test_CV <- function(sim_reps, d, pu, pc, su_num, su_samp, p_su_samp, iters){
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

#' function to test number of pop'n units
#' 
#' @export
#' 
test_PU <- function(sim_reps, d, pc, pu_cv, su_num, su_samp, p_su_samp, iters){
  #start timer
  tictoc::tic()
  # set numbers of pop'n units
  npu_test <- c(25, 100, 250, 500)
  # run simulation
  rr_npu <- purrr::map(1:sim_reps, ~purrr::map(1:length(npu_test), ~rep_sim(d, pu = npu_test[.], pc, pu_cv, su_num, su_samp, p_su_samp, iters)))
  # save & plot results
  plot_sim(rr = rr_npu, 
           plot_name = 'Npu', 
           test_vec = npu_test, 
           test_name = "PU",
           test_lab = 'Number of population units')
  # end timer
  runtime <- tictoc::toc()
  return(runtime)
}

#' function to test number of categories (i.e., longevity, growth)
#' 
#' @export
#' 
test_C <- function(sim_reps, d, pu, pu_cv, su_num, su_samp, p_su_samp, iters){
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
           test_name = "C",
           test_lab = 'Number of categories within the population')
  # end timer
  runtime <- tictoc::toc()
  return(runtime)
}

#' function to  test number of sampling units (i.e., number of hauls)
#' 
#' @export
#' 
test_SU <- function(sim_reps, d, pu, pc, pu_cv, su_samp, p_su_samp, iters){
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
           test_name = "S",
           test_lab = 'Number of sampling units')
  # end timer
  runtime <- tictoc::toc()
  return(runtime)
}

#' function to  test sample size within sampling units
#' 
#' @export
#' 
test_nSU <- function(sim_reps, d, pu, pc, pu_cv, su_num, iters, plot_name){
  #start timer
  tictoc::tic()
  # set sample size within sampling units
  samp_test <- c(10, 20, 50, 100)
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

#' function to  run experiment 2 tests in parallel
#' 
#' @export
#' 
run_exp2_tests <- function(sim_reps, d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters){
  runtime_base %<-% test_base(sim_reps, d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters) %seed% TRUE
  runtime_CV %<-% test_CV(sim_reps, d, pu, pc, su_num, su_samp, p_su_samp, iters) %seed% TRUE
  runtime_PU %<-% test_PU(sim_reps, d, pc, pu_cv, su_num, su_samp, p_su_samp, iters) %seed% TRUE
  runtime_C %<-% test_C(sim_reps, d, pu, pu_cv, su_num, su_samp, p_su_samp, iters) %seed% TRUE
  runtime_SU %<-% test_SU(sim_reps, d, pu, pc, pu_cv, su_samp, p_su_samp, iters) %seed% TRUE
  runtime_nSU_250 %<-% test6(sim_reps, d, pu, pc, pu_cv, su_num, iters, 'S250') %seed% TRUE
  runtime_nSU_500 %<-% test6(sim_reps, d, pu, pc, pu_cv, su_num = 500, iters, 'S500') %seed% TRUE
}

#' function that tests bootstrap method in experiment 3
#' 
#' @export
#' 
test_bs <- function(bs_iters, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters){
  
  #start timer
  tictoc::tic()
  # run simulation
  rr_bs <- purrr::map(1:bs_iters, ~bs_sim(d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters))
  # end timer
  runtime <- tictoc::toc()
  
  # unlist results
  do.call(mapply, c(list, rr_bs, SIMPLIFY = FALSE))$rss_se %>% 
    tidytable::map_df(., ~as.data.frame(.x), .id = "sim") -> rss_se
  do.call(mapply, c(list, rr_bs, SIMPLIFY = FALSE))$iss_bs %>% 
    tidytable::map_df(., ~as.data.frame(.x), .id = "sim") -> iss_bs
  do.call(mapply, c(list, rr_bs, SIMPLIFY = FALSE))$popn_strctr %>% 
    tidytable::map_df(., ~as.data.frame(.x), .id = "sim") -> popn_strctr
  
  res <- list(rss_se = rss_se, iss_bs = iss_bs, popn_strctr = popn_strctr)
  
  res
}


#' function to run bootstrap tests in parallel
#' 
#' @export
#' 
run_bs_test <- function(bs_iters, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters){
  
  # define parallel function runs
  run1 %<-% test_bs(bs_iters, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters) %seed% TRUE
  run2 %<-% test_bs(bs_iters, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters) %seed% TRUE
  run3 %<-% test_bs(bs_iters, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters) %seed% TRUE
  run4 %<-% test_bs(bs_iters, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters) %seed% TRUE
  run5 %<-% test_bs(bs_iters, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters) %seed% TRUE
  run6 %<-% test_bs(bs_iters, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters) %seed% TRUE
  run7 %<-% test_bs(bs_iters, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters) %seed% TRUE
  run8 %<-% test_bs(bs_iters, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters) %seed% TRUE
  run9 %<-% test_bs(bs_iters, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters) %seed% TRUE
  run10 %<-% test_bs(bs_iters, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters) %seed% TRUE
  
  # combine results
  rss_se <- run1$rss_se %>% 
    tidytable::bind_rows(run2$rss_se %>% 
                           tidytable::mutate(sim = sim + bs_iters)) %>% 
    tidytable::bind_rows(run3$rss_se %>% 
                           tidytable::mutate(sim = sim + 2 * bs_iters)) %>% 
    tidytable::bind_rows(run4$rss_se %>% 
                           tidytable::mutate(sim = sim + 3 * bs_iters)) %>% 
    tidytable::bind_rows(run5$rss_se %>% 
                           tidytable::mutate(sim = sim + 4 * bs_iters)) %>% 
    tidytable::bind_rows(run6$rss_se %>% 
                           tidytable::mutate(sim = sim + 5 * bs_iters)) %>% 
    tidytable::bind_rows(run7$rss_se %>% 
                           tidytable::mutate(sim = sim + 6 * bs_iters)) %>% 
    tidytable::bind_rows(run8$rss_se %>% 
                           tidytable::mutate(sim = sim + 7 * bs_iters)) %>% 
    tidytable::bind_rows(run9$rss_se %>% 
                           tidytable::mutate(sim = sim + 8 * bs_iters)) %>% 
    tidytable::bind_rows(run10$rss_se %>% 
                           tidytable::mutate(sim = sim + 9 * bs_iters))
  
  iss_bs <- run1$iss_bs %>% 
    tidytable::bind_rows(run2$iss_bs %>% 
                           tidytable::mutate(sim = sim + bs_iters)) %>% 
    tidytable::bind_rows(run3$iss_bs %>% 
                           tidytable::mutate(sim = sim + 2 * bs_iters)) %>% 
    tidytable::bind_rows(run4$iss_bs %>% 
                           tidytable::mutate(sim = sim + 3 * bs_iters)) %>% 
    tidytable::bind_rows(run5$iss_bs %>% 
                           tidytable::mutate(sim = sim + 4 * bs_iters)) %>% 
    tidytable::bind_rows(run6$iss_bs %>% 
                           tidytable::mutate(sim = sim + 5 * bs_iters)) %>% 
    tidytable::bind_rows(run7$iss_bs %>% 
                           tidytable::mutate(sim = sim + 6 * bs_iters)) %>% 
    tidytable::bind_rows(run8$iss_bs %>% 
                           tidytable::mutate(sim = sim + 7 * bs_iters)) %>% 
    tidytable::bind_rows(run9$iss_bs %>% 
                           tidytable::mutate(sim = sim + 8 * bs_iters)) %>% 
    tidytable::bind_rows(run10$iss_bs %>% 
                           tidytable::mutate(sim = sim + 9 * bs_iters))
  
  popn_strctr <- run1$popn_strctr %>% 
    tidytable::bind_rows(run2$popn_strctr %>% 
                           tidytable::mutate(sim = sim + bs_iters)) %>% 
    tidytable::bind_rows(run3$popn_strctr %>% 
                           tidytable::mutate(sim = sim + 2 * bs_iters)) %>% 
    tidytable::bind_rows(run4$popn_strctr %>% 
                           tidytable::mutate(sim = sim + 3 * bs_iters)) %>% 
    tidytable::bind_rows(run5$popn_strctr %>% 
                           tidytable::mutate(sim = sim + 4 * bs_iters)) %>% 
    tidytable::bind_rows(run6$popn_strctr %>% 
                           tidytable::mutate(sim = sim + 5 * bs_iters)) %>% 
    tidytable::bind_rows(run7$popn_strctr %>% 
                           tidytable::mutate(sim = sim + 6 * bs_iters)) %>% 
    tidytable::bind_rows(run8$popn_strctr %>% 
                           tidytable::mutate(sim = sim + 7 * bs_iters)) %>% 
    tidytable::bind_rows(run9$popn_strctr %>% 
                           tidytable::mutate(sim = sim + 8 * bs_iters)) %>% 
    tidytable::bind_rows(run10$popn_strctr %>% 
                           tidytable::mutate(sim = sim + 9 * bs_iters))
  
  
  res <- list(rss_se = rss_se, iss_bs = iss_bs, popn_strctr = popn_strctr)
  
  # save results
  saveRDS(res,
          file = here::here('output', paste0('exp3_bs.rds')))
  
  # plot results
  true_iss <- res_bs$rss_se %>%  
    tidytable::left_join(res_bs$popn_strctr) %>% 
    tidytable::summarise(Wtd = psych::harmonic.mean(rss_wtd, zero = FALSE),
                         Unwtd = psych::harmonic.mean(rss_unwtd, zero = FALSE),
                         .by = c(selex_type, popn_strctr)) %>% 
    tidytable::pivot_longer(cols = c(Wtd, Unwtd)) %>% 
    tidytable::rename(type = name, iss = value) %>% 
    tidytable::mutate(popn_strctr = factor(popn_strctr, levels = c('recruitment pulse', 'multimodal', 'unimodal')))
  
  bs_iss <- res_bs$iss_bs %>%  
    tidytable::left_join(res_bs$popn_strctr) %>%
    tidytable::rename(Wtd = bs_iss_wtd, Unwtd = bs_iss_unwtd) %>% 
    tidytable::pivot_longer(cols = c(Wtd, Unwtd)) %>% 
    tidytable::rename(type = name, iss = value) %>% 
    tidytable::mutate(popn_strctr = factor(popn_strctr, levels = c('recruitment pulse', 'multimodal', 'unimodal')))
  
  bs_plot_popn <- ggplot(data = bs_iss, aes(x = type, y = iss, fill = type)) + 
    geom_boxplot(alpha = 0.7) +
    geom_point(data = true_iss, shape = 24, size = 2, fill = 'white', aes(x = type, y = iss)) +
    facet_grid(selex_type ~ popn_strctr, scales = 'free_y') +
    scico::scale_color_scico_d(palette = 'roma') +
    scico::scale_fill_scico_d(palette = 'roma') +
    theme_bw() +
    guides(fill = 'none') +
    ylab("Input sample size (ISS)") +
    xlab("Composition expansion type")
  
  ggsave(filename = "exp3_bs.png",
         plot = bs_plot_popn,
         path = here::here("figs"),
         width = 6.5,
         height = 5,
         units = "in")
  
}

