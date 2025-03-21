#' function to set up a population comprised of subunits with different compositions
#' 
#' @param d population exponential decay parameters
#' @param pu number of population units (e.g., number of schools)
#' @param pc number of population categories (e.g., ages or lengths)
#' @param pu_cv CV in mean category within a population unit (e.g., spread in ages around mean age within a given school)
#' @param plot boolean, whether to output plot of generated pop'n (default = TRUE)
#' 
#' @return List of population unit compositions (p_pu), relative abundance of population units (p_popn), 
#' and the combined population composition (p_true)
#' 
#' @export
#' 
get_popn <- function(d, pu, pc, pu_cv, plot = TRUE){
  
  # get mean category for pop'n unit (e.g., mean age of school)
  mu_cat <- sample(1:pc, pu)
  
  # set up pop'n units
  p_pu <- purrr::map(1:pu, ~data.frame(cat = 1:pc, p_pu = dnorm(1:pc, mean = mu_cat[.], sd = mu_cat[.] * pu_cv) / sum(dnorm(1:pc, mean = mu_cat[.], sd = mu_cat[.] * pu_cv)))) %>% 
    tidytable::map_df(., ~as.data.frame(.x), .id = "popn_unit")
  
  # set up pop'n with decay
  popn <- data.frame(cat = 1:pc, popn = exp(-(1:pc - 1) * d))
  
  # get pop'n unit relative proportions
  p_popn <- p_pu %>% 
    tidytable::left_join(popn) %>% 
    tidytable::summarise(rel_popn = sum(p_pu * popn), .by = popn_unit) %>% 
    tidytable::mutate(p_popn = rel_popn / sum(rel_popn)) %>% 
    tidytable::select(popn_unit, rel_popn, p_popn)
  
  # get 'true' pop'n composition
  p_true <- p_pu %>% 
    tidytable::left_join(p_popn) %>% 
    tidytable::summarise(popn_tot = sum(p_pu * rel_popn), .by = c(cat)) %>% 
    tidytable::mutate(p_true = popn_tot / sum(popn_tot)) %>% 
    tidytable::select(cat, p_true)
  
  # if desired, plot generated pop'n
  if(isTRUE(plot)){
    p_pu %>% 
    tidytable::bind_rows(p_true %>% 
                           tidytable::mutate(popn_unit = 'combined') %>% 
                           tidytable::select(popn_unit, cat, p_pu = p_true)) %>% 
      tidytable::rename(comp = p_pu) -> plot_dat
    
    popn_plot <- ggplot(data = plot_dat, aes(x = as.factor(cat), y = comp, fill = popn_unit)) +
      geom_bar(stat = 'identity') +
      facet_wrap(~popn_unit, ncol = 1) +
      theme_bw() +
      xlab('category') +
      scico::scale_fill_scico_d(palette = 'roma')
    
    ggsave(filename = "gen_popn.png",
           plot = popn_plot,
           path = here::here("figs"),
           width = 6.5,
           height = 5,
           units = "in")
  }
  
  list(p_pu = p_pu, p_popn = p_popn, p_true = p_true)
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

  # generate sampling unit category samples
  # define the sampling event
  data.frame(samp_event = 1:su_num) %>% 
    # determine which population unit is being sampled (based on relative abundance)
    data.frame(popn_unit = sample(x = sim_popn$p_popn$popn_unit, 
                                  size = su_num, 
                                  replace = TRUE, 
                                  prob = sim_popn$p_popn$p_popn)) %>% 
    # join population unit composition
    tidytable::left_join(sim_popn$p_pu) %>% 
    # generate samples within categories based on population unit composition across sampling events
    tidytable::mutate(samp = stats::rmultinom(1, 
                                              # generate sampling unit number of samples
                                              sample(x = su_samp, size = 1, replace = TRUE, prob = p_su_samp),
                                              p_pu),
                      .by = samp_event) -> n_su
  
  # compute proportion of sample size across sampling units
  n_su %>% 
    tidytable::summarise(samp_se = sum(samp), .by = samp_event) %>% 
    tidytable::mutate(p_samp_se = samp_se / sum(samp_se)) -> p_samp_su

  # generate log-normal sampling unit abundance (related to population unit abundance) and calculate proportions
  data.frame(samp_event = 1:su_num, N_su = exp(stats::rnorm(su_num, 0, 1))) %>% 
    tidytable::mutate(prop_N = N_su / sum(N_su)) -> N_su 

  # compute total sample size
  n_su %>% 
    tidytable::summarise(nss = sum(samp)) -> nss
  
  # get generated composition results
  n_su %>% 
    # join sampling event abundance
    tidytable::left_join(N_su) %>% 
    # join sampling event sample size
    tidytable::left_join(p_samp_su) %>% 
    tidytable::summarise(samp_wtd = sum(samp * p_samp_se * prop_N), # weight samples
                         samp_unwtd = sum(samp), # do not weight samples
                         .by = cat) %>% 
    # compute compositions
    tidytable::mutate(samp_p_wtd = samp_wtd / sum(samp_wtd),
                      samp_p_unwtd = samp_unwtd / sum(samp_unwtd)) %>% 
    tidytable::select(cat, samp_p_wtd, samp_p_unwtd) -> comp


  # output
  list(comp = comp, nss = nss)
  
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
#' @return list of input sample size (by expansion complexity) and nominal sample size (nss)
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
                         .by = sim) -> rss_sim
  
  # compute input sample size & nominal sample size
  rss_sim %>% 
    tidytable::left_join(do.call(mapply, c(list, rr_sim, SIMPLIFY = FALSE))$nss %>% 
                           tidytable::map_df(., ~as.data.frame(.x), .id = "sim")) %>% 
    tidytable::summarise(iss_wtd = psych::harmonic.mean(rss_wtd, zero = FALSE),
                         iss_unwtd = psych::harmonic.mean(rss_unwtd, zero = FALSE),
                         mean_nss = mean(nss)) -> iss_sim
  
  list(iss_sim = iss_sim)
  
}