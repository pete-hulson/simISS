

# load/source stuff ----
library(tidyverse)
source(here::here('R', 'base_functions.R'))


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
pu <- 5

# number of population categories (e.g., ages or length bins)
pc <- 15

# CV around mean within pop'n unit (e.g., CV in mean age of school)
pu_cv <- 0.25

# pop'n exponential decay (e.g., lnM, tied to inverse of number of pop'n categories so that pop'n = 0.01 at largest category)
d <- log(0.01) / (1 - pc)


# experiment 1: ----
# do we get unbiased estimates of pop'n composition sampling different schools?
# what is effect of expansion complexity?

## get simulated pop'n ----
# note: sample pop'n category sampling without replacement for example plots
sim_popn <- get_popn(d, pu, pc, pu_cv, replace_pc = FALSE, plot_name = 'exp1_popn')

## run sim loop ----
tictoc::tic()
rr_exp1 <- purrr::map(1:iters, ~sim_comp(su_num, sim_popn, su_samp, p_su_samp))
runtime_exp1 <- tictoc::toc()

# unlist results
do.call(mapply, c(list, rr_exp1, SIMPLIFY = FALSE))$comp %>% 
  tidytable::map_df(., ~as.data.frame(.x), .id = "sim") -> res_sim

exp1_res <- list(sim_popn = sim_popn, res_sim = res_sim)


# observed
get_comps <- function(comps) {
  comps %>%
    filter(selex_type == 'uniform') %>%
    pull(samp_p_wtd)
}

get_iss <- function(iss) {
  iss %>%
    filter(selex_type == 'uniform') %>%
    pull(nss)
}

# observed
o <- do.call(rbind, lapply(rr_exp1, function(x) get_comps(x$comp)))
iss <- do.call(rbind, lapply(rr_exp1, function(x) get_iss(x$nss)))


# true
e <- (exp1_res$sim_popn$p_true %>% filter(selex_type == 'uniform'))$p_true

ddirmult = function(obs, pred, Ntotal, ln_theta, give_log = TRUE) {
  # Set up function variables
  n_c = length(obs) # number of categories
  p_exp = pred # expected values container
  p_obs = obs # observed values container
  dirichlet_Parm = exp(ln_theta) * Ntotal # Dirichlet alpha parameters
  
  # set up pdf
  logres = lgamma(Ntotal + 1)
  for(c in 1:n_c) logres = logres - lgamma(Ntotal*p_obs[c]+1) # integration constant
  logres = logres + lgamma(dirichlet_Parm) - lgamma(Ntotal+dirichlet_Parm) # 2nd term in formula
  
  # Summation in 3rd term in formula
  for(c in 1:n_c) {
    logres = logres + lgamma(Ntotal*p_obs[c] + dirichlet_Parm*p_exp[c])
    logres = logres - lgamma(dirichlet_Parm * p_exp[c])
  } # end c
  
  if(give_log == TRUE) return(logres)
  else return(exp(logres))
} # end function


ddirmult_all <- function(pars) {
  RTMB::getAll(data, pars)
  n_yrs = nrow(data$o)
  theta = exp(ln_theta)
  nll = 0
  for(i in 1:n_yrs) nll = nll - ddirmult(o[i,], e, iss[i,], ln_theta, TRUE)
  ess = (1 + theta * iss) / (1 + theta)
  RTMB::REPORT(ess)
  return(nll)
}

data <- list(o = o, iss = iss, e = e)
pars <- list(ln_theta = log(3))

obj <- RTMB::MakeADFun(ddirmult_all, pars)

opt <- stats::nlminb(obj$par, obj$fn, obj$gr,
                     control = list(iter.max = 1e5, eval.max = 1e5, rel.tol = 1e-15))


try_improve <- tryCatch(expr =
                          for(i in 1:3) {
                            g = as.numeric(obj$gr(opt$par))
                            h = optimHess(opt$par, fn = obj$fn, gr = obj$gr)
                            opt$par = opt$par - solve(h,g)
                            opt$objective = obj$fn(opt$par)
                          }
                        , error = function(e){e}, warning = function(w){w})

sd <- RTMB::sdreport(obj)
sd

rep <- obj$report(obj$env$par)
rep
