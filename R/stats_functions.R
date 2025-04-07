#' Function to estimate composition statistics
#'
#' @param rr list of replicated results from sim_comp() function
#' @param sim_popn list of results for simulated population from get_popn() function

#' 
#' @return estimates of input sample size and logistic-normal parameters
#' 
#' @export
#'
est_stats <- function(rr, sim_popn){
  
  # unlist results
  res_sim <- do.call(mapply, c(list, rr, SIMPLIFY = FALSE))$comp %>% 
    tidytable::map_df(., ~as.data.frame(.x), .id = "sim") %>% 
    tidytable::pivot_longer(cols = c('samp_p_wtd', 'samp_p_unwtd')) %>% 
    tidytable::rename(comp_type = name, p_obs = value) %>% 
    tidytable::mutate(comp_type = case_when(comp_type == 'samp_p_wtd' ~ 'wtd',
                                            .default = 'unwtd'))
  
  # input sample size statistic ----
  iss_sim <- res_sim %>% 
    tidytable::left_join(sim_popn$p_true) %>% 
    tidytable::summarise(rss = sum(p_true * (1- p_true)) / sum((p_obs - p_true) ^ 2),
                         .by = c(sim, selex_type, comp_type)) %>% 
    tidytable::summarise(iss = psych::harmonic.mean(rss, zero = FALSE),
                         .by = c(selex_type, comp_type))
  
  
  # logistic-normal statistics ----
  
  # set up data for logistic-normal
  data <- list(exp = sim_popn$p_true %>% 
                 tidytable::select(-N_c), 
               obs = res_sim,
               iss = iss_sim)
  
  # remove 0's
  if(any(data$obs == 0) || any(data$exp == 0)) {
    # small constant
    eps <- 1e-6
    
    # expected
    data$exp <- data$exp %>% 
      # add constant in
      tidytable::mutate(p_true = case_when(p_true == 0 ~ eps,
                                           .default = p_true)) %>% 
      # renormalize
      tidytable::mutate(p_true = p_true / sum(p_true), 
                        .by = c(selex_type))
    
    # observed
    data$obs <- data$obs %>% 
      # add constant in
      tidytable::mutate(p_obs = case_when(p_obs == 0 ~ eps,
                                          .default = p_obs)) %>% 
      # renormalize
      tidytable::mutate(p_obs = p_obs / sum(p_obs), 
                        .by = c(sim, selex_type, comp_type))
  }
  
  # estimate stats
  
  # dataframe of selectivity/composition expansion types tested
  combs <- tidytable::expand_grid(selex = unique(data$iss$selex_type), 
                                  comp = unique(data$iss$comp_type))
  
  # run for iid
  rr_iid <- purrr::map(1:dim(combs)[1],
                       ~est_logistic_normal(start_sigma = log(10),
                                            cov_strc = 'iid',
                                            data = data, 
                                            selex_t = combs$selex[.],
                                            comp_t = combs$comp[.]))
  
  # run for 1DAR1
  rr_1DAR1 <- purrr::map(1:dim(combs)[1],
                         ~est_logistic_normal(start_sigma = log(10),
                                              start_rho = 0.1,
                                              cov_strc = '1DAR1',
                                              data = data, 
                                              selex_t = combs$selex[.],
                                              comp_t = combs$comp[.]))
  
  # unlist results
  # iid
  do.call(mapply, c(list, rr_iid, SIMPLIFY = FALSE))$res %>% 
    tidytable::map_df(., ~as.data.frame(.x), .id = "comb") %>% 
    tidytable::select(selex_type, comp_type, sigma_iid = sigma) %>% 
    tidytable::left_join(
      # 1DAR1
      do.call(mapply, c(list, rr_1DAR1, SIMPLIFY = FALSE))$res %>% 
        tidytable::map_df(., ~as.data.frame(.x), .id = "comb") %>% 
        tidytable::select(selex_type, comp_type, sigma_1DAR1 = sigma, rho_1DAR1 = rho)) -> logistN_sim
  
  
  
  # put results together
  iss_sim %>% 
    tidytable::left_join(logistN_sim) -> res
  
  res
}


#' Compute negative log likelihood for a logistic-normal
#'
#' @param data data list of expected ('true'), observed (simulated), and realized sample size
#' @param sigma variance parameter (default = NULL)
#' @param rho correlation parameter (default = NULL)
#' @param cov_strc covariance structure options ("iid", "1DAR1")
#' 
#' @return negative log-likelihood value
#' 
#' @export
#'
nll_logistic_normal <- function(data,
                                sigma = NULL,
                                rho = NULL,
                                cov_strc = NULL) {
  
  # Get observed and predicted data here
  obs <- data$obs # observed
  exp <- data$exp # predicted
  
  # Get dimensions
  Ncat <- length(unique(obs$cat)) # number of categories
  
  # calculate sigma
  sigma2 <- data$iss %>% 
    # tidytable::mutate(sigma2 = sigma^2 / iss) # matt's way where divide by iss
    tidytable::mutate(sigma2 = sigma^2)
  
  # set up covariance matrix
  # for iid covariance structure
  if(cov_strc == "iid"){
    covmat <- diag(rep(as.numeric(sigma2 %>% 
                                    tidytable::select(sigma2)), Ncat))
  }
  # for 1DAR1 covariance structure
  if(cov_strc == "1DAR1"){
    corrMatrix <- matrix(0, nrow = Ncat, ncol = Ncat)
    for (i in 1:Ncat) {
      for (j in 1:Ncat) {
        # Calculate the correlation based on the lag distance
        corrMatrix[i, j] <- rho^(abs(i - j))
      } # end i
    } # end j
    covmat <-  as.numeric(sigma2 %>% 
                            tidytable::select(sigma2)) * corrMatrix
  }
  
  # do logistic transformation on observed values
  tmp_Obs <- obs %>% 
    tidytable::filter(cat != Ncat) %>% 
    tidytable::left_join(obs %>% 
                           tidytable::filter(cat == Ncat) %>% 
                           tidytable::select(sim, selex_type, comp_type, p_Ncat = p_obs)) %>% 
    tidytable::mutate(tmp_Obs = log(p_obs) - log(p_Ncat)) %>% 
    tidytable::select(sim, cat, selex_type, comp_type, tmp_Obs)
  
  # do logistic transformation on expected values
  mu <- exp %>% 
    tidytable::filter(cat != Ncat) %>% 
    tidytable::mutate(mu = log(p_true) - 
                        log(as.numeric(exp %>% 
                                         tidytable::filter(cat == Ncat) %>% 
                                         tidytable::select(p_true)))) %>% 
    tidytable::select(cat, selex_type, mu)
  
  # compute negative log-likelihood
  negloglik <- as.numeric(tmp_Obs %>% 
                            tidytable::left_join(mu) %>% 
                            tidytable::summarise(negloglik = -1 * RTMB::dmvnorm(tmp_Obs, mu, covmat[-nrow(covmat), -ncol(covmat)], log = TRUE),
                                                 .by = sim) %>% 
                            tidytable::summarise(negloglik = sum(negloglik)))
  
  return(negloglik)
}

#' Function to estimate logistic-normal parameters
#'
#' @param start_sigma starting value for sigma (default = NULL)
#' @param start_rho Starting value for correlation (default = NULL, bound between -1 and 1)
#' @param cov_strc covariance structure options ("iid", "1DAR1")
#' @param data data list of expected ('true'), observed (simulated), and realized sample size
#' @param selex_t selectivity option (default = NULL)
#' @param comp_t composition expansion option (default = NULL)
#' 
#' @return estimated parameters for logistic-normal distribution
#' 
#' @export
#'
est_logistic_normal <- function(start_sigma = NULL,
                                start_rho = NULL,
                                cov_strc = NULL,
                                data = NULL, 
                                selex_t = NULL,
                                comp_t = NULL){
  
  # for 'iid' covariance structure
  if(cov_strc == 'iid'){
    # Define new likelihood function to minimize
    negloglik <- function(pars, data){
      sigma <- exp(pars[1]) # gets squared later 
      nLL <- nll_logistic_normal(data = data, 
                                 sigma = sigma,
                                 cov_strc = cov_strc)
      return(nLL)
    } # end likelihood function
    
    # minimize
    fit <- nlminb(start = start_sigma, 
                  objective = negloglik, 
                  data = list(obs = data$obs %>% 
                                tidytable::filter(selex_type == selex_t,
                                                  comp_type == comp_t),
                              exp = data$exp %>% 
                                tidytable::filter(selex_type == selex_t),
                              iss =  data$iss %>% 
                                tidytable::filter(selex_type == selex_t,
                                                  comp_type == comp_t)), 
                  hessian = TRUE, 
                  control = list(iter.max = 1e7,
                                 eval.max = 1e7, 
                                 rel.tol = 1e-10))
    
    # list results
    res <- list(res = data.frame(sigma = exp(fit$par), selex_type = selex_t, comp_type = comp_t))
  }
  
  # for 1DAR1 covariance structure
  if(cov_strc == '1DAR1'){
    # Define new likelihood function to minimize
    rho_trans <- function(x) 2/(1+ exp(-2 * x)) - 1
    negloglik <- function(pars, data){
      sigma <- exp(pars[1]) # gets squared later 
      rho <- rho_trans(pars[2]) # constrain estimation to -1 and 1
      nLL <- nll_logistic_normal(data = data, 
                                 sigma = sigma,
                                 rho = rho,
                                 cov_strc = cov_strc)
      return(nLL)
    }
    
    # minimize
    fit <- nlminb(start = c(start_sigma, start_rho), 
                  objective = negloglik, 
                  data = list(obs = data$obs %>% 
                                tidytable::filter(selex_type == selex_t,
                                                  comp_type == comp_t),
                              exp = data$exp %>% 
                                tidytable::filter(selex_type == selex_t),
                              iss =  data$iss %>% 
                                tidytable::filter(selex_type == selex_t,
                                                  comp_type == comp_t)), 
                  hessian = TRUE, 
                  control = list(iter.max = 1e7,
                                 eval.max = 1e7, 
                                 rel.tol = 1e-10))
    
    # list results
    res <- list(res = data.frame(sigma = exp(fit$par[1]), rho = fit$par[2], selex_type = selex_t, comp_type = comp_t))
    
  }
  
  # return results
  res
}
