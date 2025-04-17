library(tidyverse)
library(future)
source(here::here('R', 'base_functions.R'))
source(here::here('R', 'stats_functions.R'))
source(here::here('R', 'exp2_functions.R'))
source(here::here('R', 'exp3_functions.R'))


# test logisticN ----

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
pu <- 25

# number of population categories (e.g., ages or length bins)
pc <- 15

# CV around mean within pop'n unit (e.g., CV in mean age of school)
pu_cv <- 0.25

# pop'n exponential decay (e.g., lnM, tied to inverse of number of pop'n categories so that pop'n = 0.01 at largest category)
d <- log(0.01) / (1 - pc)


# get simulated pop'n

# trying to get a dome-shaped unimodal popn...

sim_popn <- get_popn(d, pu, pc, pu_cv, plot = FALSE)

sim_popn$popn_strctr

# run sim loop
rr_sim <- purrr::map(1:iters, ~sim_comp(su_num, sim_popn, su_samp, p_su_samp))



# unlist results
res_sim <- do.call(mapply, c(list, rr_sim, SIMPLIFY = FALSE))$comp %>% 
  tidytable::map_df(., ~as.data.frame(.x), .id = "sim") %>% 
  tidytable::pivot_longer(cols = c('samp_p_wtd', 'samp_p_unwtd'), names_to = 'comp_type', values_to = 'p_obs') %>% 
  tidytable::mutate(comp_type = case_when(comp_type == 'samp_p_wtd' ~ 'wtd',
                                          .default = 'unwtd'))


# remove sims with 0's
sim_rm <- res_sim %>% 
  tidytable::filter(p_obs == 0) %>% 
  tidytable::distinct(sim, selex_type) %>% 
  tidytable::mutate(rem = 1)
  
res_sim <- res_sim %>% 
  tidytable::left_join(sim_rm) %>% 
  tidytable::filter(is.na(rem)) %>% 
  tidytable::select(-rem)


# logistic-normal statistics ----

# set up data list
data <- list(exp = sim_popn$p_true %>% 
               tidytable::select(-N_c), 
             obs = res_sim,
             N = do.call(mapply, c(list, rr_sim, SIMPLIFY = FALSE))$nss %>% 
               tidytable::map_df(., ~as.data.frame(.x), .id = "sim") %>% 
               tidytable::left_join(sim_rm) %>% 
               tidytable::filter(is.na(rem)) %>% 
               tidytable::select(-rem))


# combinations of selectivity/composition expansion types tested
combs <- tidytable::expand_grid(selex = unique(data$obs$selex_type), 
                                comp = unique(data$obs$comp_type))

# run for iid
rr_iid <- purrr::map(1:dim(combs)[1],
                     ~est_logistic_normal(cov_strc = 'iid',
                                          data = data, 
                                          selex_t = combs$selex[.],
                                          comp_t = combs$comp[.]))

# unlist results
do.call(mapply, c(list, rr_iid, SIMPLIFY = FALSE))$res %>% 
  tidytable::map_df(., ~as.data.frame(.x), .id = "comb") %>% 
  tidytable::select(selex_type, comp_type, sigma_iid) -> logistN_sim

logistN_sim



ggplot(data = sim_popn$p_true, aes(x = cat, y = p_true)) +
  geom_bar(stat = 'identity', aes(col = selex_type), alpha = 0) +
  geom_line(data = res_sim %>% 
              tidytable::summarise(p_obs = mean(p_obs), .by = c(cat, selex_type, comp_type)),
            aes(x = cat, y = p_obs, linetype = comp_type)) +
  facet_wrap(~selex_type, ncol = 1, scales = 'free_y')

res_sim %>% 
  tidytable::summarise(p_obs = mean(p_obs), .by = c(cat, selex_type, comp_type))


sim_popn$p_true %>% 
  select(-N_c) %>% 
  tidytable::pivot_wider(names_from = selex_type, values_from = p_true)

res_sim %>% 
  tidytable::filter(selex_type == 'dome-shaped',
                    comp_type == 'wtd') %>% 
  select(-c(selex_type, comp_type)) %>% 
  tidytable::pivot_wider(names_from = sim, values_from = p_obs)


res_sim %>% 
  tidytable::filter(selex_type == 'dome-shaped',
                    comp_type == 'wtd') %>% 
  select(-c(selex_type, comp_type)) %>% 
  filter(cat == 15,
         p_obs == 0)





selex_t = 'dome-shaped'
comp_t = 'wtd'
# set up data
# remove 0's
if(any(data$obs == 0) || any(data$exp == 0)) {
  # small constant
  eps <- 1e-4
  
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

sort(log(o[, n_c]))

# for 'iid' covariance structure
# set up data/parameters
data <- list(o = o, e = e)
pars <- list(sigma = log(10))
# make model
obj <- RTMB::MakeADFun(dlogistN_iid, pars)
# run model
invisible(capture.output(opt <- stats::nlminb(obj$par, obj$fn, obj$gr,
                                              control = list(iter.max = 1e5, eval.max = 1e5, rel.tol = 1e-15))))
# get output
list(res = data.frame(sigma_iid = as.numeric(exp(opt$par)), selex_type = selex_t, comp_type = comp_t))







# dirichlet-multinomial statistic ----
# run model
rr_DM <- purrr::map(1:dim(combs)[1],
                    ~est_dirmult(data,
                                 selex_t = combs$selex[.],
                                 comp_t = combs$comp[.]))

# unlist results
DM_sim <- do.call(mapply, c(list, rr_DM, SIMPLIFY = FALSE))$res %>% 
  tidytable::map_df(., ~as.data.frame(.x), .id = "comb") %>% 
  tidytable::select(-comb)






