#' Function to estimate composition statistics
#'
#' @param rr_sim list of replicated results from sim_comp() function
#' @param sim_popn list of results for simulated population from get_popn() function
#' @param cov_strc logistic-normal covariance structure options ("iid" and/or "1DAR1")
#' 
#' @return estimates of input sample size and logistic-normal parameters
#' 
#' @export
#'
est_stats <- function(rr_sim, sim_popn, cov_strc){
  
  # unlist results
  res_sim <- do.call(mapply, c(list, rr_sim, SIMPLIFY = FALSE))$comp %>% 
    tidytable::map_df(., ~as.data.frame(.x), .id = "sim") %>% 
    tidytable::pivot_longer(cols = c('samp_p_wtd', 'samp_p_unwtd'), names_to = 'comp_type', values_to = 'p_obs') %>% 
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
  
  # set up data list
  data <- list(exp = sim_popn$p_true %>% 
                 tidytable::select(-N_c), 
               obs = res_sim,
               N = do.call(mapply, c(list, rr_sim, SIMPLIFY = FALSE))$nss %>% 
                 tidytable::map_df(., ~as.data.frame(.x), .id = "sim"))

  
  # combinations of selectivity/composition expansion types tested
  combs <- tidytable::expand_grid(selex = unique(data$obs$selex_type), 
                                  comp = unique(data$obs$comp_type))
  
  # run for iid
  if('iid' %in% cov_strc){
    rr_iid <- purrr::map(1:dim(combs)[1],
                         ~est_logistic_normal(cov_strc = 'iid',
                                              data = data, 
                                              selex_t = combs$selex[.],
                                              comp_t = combs$comp[.]))
  }
  
  # run for 1DAR1
  if('1DAR1' %in% cov_strc){
    rr_1DAR1 <- purrr::map(1:dim(combs)[1],
                           ~est_logistic_normal(cov_strc = '1DAR1',
                                                data = data, 
                                                selex_t = combs$selex[.],
                                                comp_t = combs$comp[.]))
  }
  
  # unlist results
  if('iid' %in% cov_strc & '1DAR1' %in% cov_strc){
    # iid
    do.call(mapply, c(list, rr_iid, SIMPLIFY = FALSE))$res %>% 
      tidytable::map_df(., ~as.data.frame(.x), .id = "comb") %>% 
      tidytable::select(selex_type, comp_type, sigma_iid) %>% 
      tidytable::left_join(
        # 1DAR1
        do.call(mapply, c(list, rr_1DAR1, SIMPLIFY = FALSE))$res %>% 
          tidytable::map_df(., ~as.data.frame(.x), .id = "comb") %>% 
          tidytable::select(-comb)) -> logistN_sim
  }
  if(length(cov_strc) == 1 & 'iid' %in% cov_strc){
    # iid
    do.call(mapply, c(list, rr_iid, SIMPLIFY = FALSE))$res %>% 
      tidytable::map_df(., ~as.data.frame(.x), .id = "comb") %>% 
      tidytable::select(selex_type, comp_type, sigma_iid) -> logistN_sim
  }
  if(length(cov_strc) == 1 & '1DAR1' %in% cov_strc){
    # 1DAR1
    do.call(mapply, c(list, rr_1DAR1, SIMPLIFY = FALSE))$res %>% 
      tidytable::map_df(., ~as.data.frame(.x), .id = "comb") %>% 
      tidytable::select(selex_type, comp_type, sigma_1DAR1, rho_1DAR1) -> logistN_sim
  }
  
  # dirichlet-multinomial statistic ----
  # run model
  rr_DM <- purrr::map(1:dim(combs)[1],
                      ~est_dirmult(rr_sim, 
                                   sim_popn, 
                                   selex_t = combs$selex[.],
                                   comp_t = combs$comp[.]))
  
  # unlist results
  DM_sim <- do.call(mapply, c(list, rr_DM, SIMPLIFY = FALSE))$res %>% 
    tidytable::map_df(., ~as.data.frame(.x), .id = "comb") %>% 
    tidytable::select(-comb)
  
  # put results together ----
  iss_sim %>% 
    tidytable::left_join(logistN_sim) %>% 
    tidytable::left_join(DM_sim) -> res
  
  res
}

#' Function to estimate logistic-normal parameters
#'
#' @param cov_strc covariance structure options ("iid", "1DAR1")
#' @param data data list of expected ('true') and observed (simulated) composition, and total sample size
#' @param selex_t selectivity option (default = NULL)
#' @param comp_t composition expansion option (default = NULL)
#' 
#' @return estimated parameters for logistic-normal distribution
#' 
#' @export
#'
est_logistic_normal <- function(cov_strc = NULL,
                                data = NULL, 
                                selex_t = NULL,
                                comp_t = NULL){
  
  # set up data
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
  
  # observed
  o <- data$obs %>%
    tidytable::filter(selex_type == selex_t,
                      comp_type == comp_t) %>% 
    tidytable::pivot_wider(names_from = cat, values_from = p_obs) %>% 
    tidytable::select(-c(sim, selex_type, comp_type)) %>% 
    as.matrix(.)
  
  # true/expected
  e <- (data$exp %>% 
          tidytable::filter(selex_type == selex_t))$p_true
  
  # run RTMB models
  
  # define iid likelihood functions
  dlogistN_iid <- function(pars){
    library(RTMB)
    "c" <- RTMB::ADoverload("c")
    "[<-" <- RTMB::ADoverload("[<-")
    # load starting values and data
    RTMB::getAll(data, pars)
    # Get dimensions
    n_c = length(e) # number of categories
    # set up covariance matrix
    covmat = diag(rep(sigma^2, n_c))
    # do logistic transformation on observed values
    tmp_Obs = log(o[, -n_c]) - log(o[, n_c])
    # do logistic transformation on expected values
    mu = log(e[-n_c]) - log(e[n_c])
    # compute negative log-likelihood:
    nll = sum(-1 * RTMB::dmvnorm(tmp_Obs, mu, Sigma = covmat[-nrow(covmat), -ncol(covmat)], log = TRUE))
    return(nll)
  }
  # define 1DAR1 likelihood functions
  dlogistN_1DAR1 <- function(pars){
    library(RTMB)
    "c" <- RTMB::ADoverload("c")
    "[<-" <- RTMB::ADoverload("[<-")
    # load starting values and data
    RTMB::getAll(data, pars)
    # Get dimensions
    n_c = length(e) # number of categories
    # set up covariance matrix
    corrMatrix <- matrix(0, nrow = n_c, ncol = n_c)
    # constrain rho estimation to -1 and 1
    rho_trans <- function(x) 2/(1+ exp(-2 * x)) - 1
    rho_t <- rho_trans(rho)
    for (i in 1:n_c) {
      for (j in 1:n_c) {
        # Calculate the correlation based on the lag distance
        corrMatrix[i, j] <- rho_t^(abs(i - j))
      } # end i
    } # end j
    covmat <- sigma^2 * corrMatrix
    # do logistic transformation on observed values
    tmp_Obs = log(o[, -n_c]) - log(o[, n_c])
    # do logistic transformation on expected values
    mu = log(e[-n_c]) - log(e[n_c])
    # compute negative log-likelihood:
    nll = sum(-1 * RTMB::dmvnorm(tmp_Obs, mu, Sigma = covmat[-nrow(covmat), -ncol(covmat)], log = TRUE))
    return(nll)
  }

  # for 'iid' covariance structure
  if(cov_strc == 'iid'){
    # set up data/parameters
    data <- list(o = o, e = e)
    pars <- list(sigma = log(10))
    # make model
    obj <- RTMB::MakeADFun(dlogistN_iid, pars)
    # run model
    invisible(capture.output(opt <- stats::nlminb(obj$par, obj$fn, obj$gr,
                                                  control = list(iter.max = 1e5, eval.max = 1e5, rel.tol = 1e-15))))
    # get output
    res <- list(res = data.frame(sigma_iid = as.numeric(exp(opt$par)), selex_type = selex_t, comp_type = comp_t))
  }
  
  # for 1DAR1 covariance structure
  if(cov_strc == '1DAR1'){
    # set up data/parameters
    data <- list(o = o, e = e)
    pars <- list(sigma = log(10),
                 rho = 0.1)
    # make model
    obj <- RTMB::MakeADFun(dlogistN_1DAR1, pars)
    # run model
    invisible(capture.output(opt <- stats::nlminb(obj$par, obj$fn, obj$gr,
                                                  control = list(iter.max = 1e5, eval.max = 1e5, rel.tol = 1e-15))))
    # get output
    res <- list(res = data.frame(sigma_1DAR1 = as.numeric(exp(opt$par[1])),
                                 rho_1DAR1 = as.numeric(2/(1+ exp(-2 * opt$par[2])) - 1), 
                                 selex_type = selex_t, comp_type = comp_t))
  }
  
  # return results
  res
}

#' Function to estimate Dirichlet-Multinomial parameter
#'
#' @param data data list of expected ('true') and observed (simulated) composition, and total sample size
#' @param selex_t selectivity option (default = NULL)
#' @param comp_t composition expansion option (default = NULL)
#' 
#' @return estimated parameters for D-M distribution
#' 
#' @export
#'
est_dirmult <- function(data,
                        selex_t = NULL,
                        comp_t = NULL){
  
  # set up data
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
  
  # observed
  o <- data$obs %>% 
    tidytable::filter(selex_type == selex_t,
                      comp_type == comp_t) %>% 
    tidytable::pivot_wider(names_from = cat, values_from = p_obs) %>% 
    tidytable::select(-c(sim, selex_type, comp_type)) %>% 
    as.matrix(.)

  # total sample size
  N <- data$N %>% 
    tidytable::filter(selex_type == selex_t) %>% 
    tidytable::pull(nss) %>% 
    as.matrix(.)

  # true/expected
  e <- (data$exp %>% 
          tidytable::filter(selex_type == selex_t))$p_true

  # run model
  # define model
  ddirmult_test <- function(pars) {
    library(RTMB)
    "c" <- RTMB::ADoverload("c")
    "[<-" <- RTMB::ADoverload("[<-")
    # load starting values and data
    RTMB::getAll(data, pars)
    # Set up likelihood variables
    n_c = length(e) # number of categories
    dirichlet_Parm = exp(ln_theta) * N # Dirichlet alpha parameters  
    theta = exp(ln_theta)
    ess = (1 + theta * N) / (1 + theta)
    # set up pdf
    logres = lgamma(N + 1)
    for(c in 1:n_c) logres = logres - lgamma(N * o[, c] + 1) # integration constant
    logres = logres + lgamma(dirichlet_Parm) - lgamma(N + dirichlet_Parm) # 2nd term in formula
    # Summation in 3rd term in formula
    for(c in 1:n_c) {
      logres = logres + lgamma(N * o[, c] + dirichlet_Parm * e[c])
      logres = logres - lgamma(dirichlet_Parm * e[c])
    } # end c
    nll = sum(-1 * logres)
    RTMB::REPORT(ess)
    return(nll)
  }
  # set up data/parameters
  data <- list(o = o, N = N, e = e)
  pars <- list(ln_theta = log(3))
  # make model
  obj <- RTMB::MakeADFun(ddirmult_test, pars)
  # optimize
  invisible(capture.output(opt <- stats::nlminb(obj$par, obj$fn, obj$gr,
                                                control = list(iter.max = 1e5, eval.max = 1e5, rel.tol = 1e-15))))
  # get output
  res <- list(res = data.frame(theta = as.numeric(exp(opt$par)), 
                               ess_DM = mean(obj$report()$ess), 
                               selex_type = selex_t, comp_type = comp_t))
  res
  
}
