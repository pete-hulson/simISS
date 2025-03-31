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
