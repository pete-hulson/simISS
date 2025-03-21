
library(tidyverse)
source(here::here('R', 'functions.R'))



# set up experiment parameters ----

## simulation/bootstrap parameters ----

# number of simulation/bootstrap iterations
iters <- 1000

# number of bootstrap replicates (for testing whether bootstrap provides unbiased estimate of 'true' iss)
bs_reps <- 10

## sampling parameters ----

# number of sampling units (e.g., hauls)
su_num <- 500

# number of samples per sampling unit (e.g., ages/lengths)
# vector to test differences in number of samples across sampling unit
su_samp <- c(20, 10)

# vector pf probabilities to select options of number of samples in sampling unit 
p_su_samp <- c(0.9, 0.1)

## population unit parameters ----

# number of population units sampled (e.g., number of schools/year classes)
pu <- 5

# number of population categories (e.g., ages or length bins)
pc <- 20

# CV around mean within pop'n unit (e.g., CV in mean age of school)
pu_cv <- 0.2

# pop'n exponential decay (e.g., lnM, tied to inverse of number of pop'n categories so that pop'n = 0.01 at largest category)
d <- log(0.01) / (1 - pc)

# experiment 1: ----
# do we get unbiased estimates of pop'n composition sampling different schools?
# what is effect of expansion complexity?


## get simulated pop'n ----

sim_popn <- get_popn(d, pu, pc, pu_cv)


## run sim loop ----
rr_sim <- purrr::map(1:iters, ~sim_comp(su_num, sim_popn, su_samp, p_su_samp))

do.call(mapply, c(list, rr_sim, SIMPLIFY = FALSE))$comp %>% 
  tidytable::map_df(., ~as.data.frame(.x), .id = "sim") -> res_sim

## plot results ----
res_sim %>% 
  tidytable::pivot_longer(cols = c(samp_p_wtd, samp_p_unwtd)) %>% 
  tidytable::select(cat, comp_type = name, comp = value) %>% 
  tidytable::bind_rows(res_sim %>% 
                         tidytable::filter(sim == sample(1:iters, 1)) %>% 
                         tidytable::select(-sim) %>% 
                         tidytable::rename(rand_wtd = samp_p_wtd,
                                           rand_unwtd = samp_p_unwtd) %>% 
                         tidytable::pivot_longer(cols = c(rand_wtd, rand_unwtd)) %>% 
                         tidytable::rename(comp_type = name, comp = value)) %>% 
  tidytable::bind_rows(sim_popn$p_true %>% 
                         tidytable::mutate(comp_type = 'true') %>% 
                         tidytable::rename(comp = p_true)) -> plot_dat

sim_plot <- ggplot(data = plot_dat, aes(x = as.factor(cat), y = comp)) +
  geom_bar(data = plot_dat %>% tidytable::filter(comp_type %in% c('true')), stat = 'identity', 
           aes(fill = comp_type), alpha = 0.5) +
  geom_boxplot(data = plot_dat %>% tidytable::filter(comp_type %in% c('samp_p_wtd', 'samp_p_unwtd')), 
               aes(fill = comp_type)) +
  geom_line(data = plot_dat %>% tidytable::filter(comp_type %in% c('rand_wtd', 'rand_unwtd')) %>% tidytable::left_join(data.frame(comp_type = c('rand_wtd', 'rand_unwtd'), col = c(scico::scico(3, palette = 'roma')[1], scico::scico(3, palette = 'roma')[2]))), 
            aes(x = cat, col = comp_type), linetype = 'dashed') +
  scico::scale_fill_scico_d(palette = 'roma')  +
  scale_color_manual(values = c(scico::scico(3, palette = 'roma')[1], scico::scico(3, palette = 'roma')[2])) +
  theme_bw() +
  xlab('category')

ggsave(filename = "sim_popn.png",
       plot = sim_plot,
       path = here::here("figs"),
       width = 6.5,
       height = 5,
       units = "in")

# notes:
# so long as the schools are represented across the haul samples in proportion to their relative abundance, 
# the mean across iterations in unbiased, whether weighted or not, so,
# you don't need to weight by haul's catch even if the number of samples within a haul aren't proportional to catch
# however, for any given iteration, the composition resulting from weighted or unweighted comps can be different.
# the uncertainty around weighted comps is larger than unweighted comps

# experiment 2: ----
# what are the factors that influence input sample size?

## test expansion weighting ----

# number of simulation replicates (for testing iss axes of influence)
sim_reps <- 1000

# run simulation
tictoc::tic()
rr_exp <- purrr::map(1:sim_reps, ~rep_sim(d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters))
tictoc::toc()

# unlist results
do.call(mapply, c(list, rr_exp, SIMPLIFY = FALSE))$iss_sim %>% 
  tidytable::map_df(., ~as.data.frame(.x), .id = "rep") -> res_exp

# plot results
res_exp %>% 
  tidytable::pivot_longer(cols = c(iss_wtd, iss_unwtd, mean_nss)) -> plot_dat

exp_plot <- ggplot(data = plot_dat, aes(x = name, y = value, fill = name)) +
  geom_boxplot(data = plot_dat %>% 
                 tidytable::filter(name %in% c('iss_wtd', 'iss_unwtd'))) +
  geom_hline(yintercept = as.numeric(plot_dat %>% 
              tidytable::filter(name %in% c('mean_nss')) %>% 
              tidytable::summarise(nss = mean(value))),
             linewidth = 1,
             colour = scico::scico(3, palette = 'roma')[3]) +
  scale_fill_manual(values = c(scico::scico(3, palette = 'roma')[1], scico::scico(3, palette = 'roma')[2])) +
  theme_bw() +
  xlab(NULL)

ggsave(filename = "exp_sim.png",
       plot = exp_plot,
       path = here::here("figs"),
       width = 6.5,
       height = 5,
       units = "in")


## other things to test: ----
# comp structure
# school structure (i.e., cv)
# number of pop'n units
# number of categories (i.e., longevity, growth)




# weighting reduces the iss, both result in a sim that results in reduction from nss



# experiment 2: can a given realization of a sample give the 'true' iss back? ----

# number of bootstrap replicates
bs_iters <- 10

## get a sample realization ----

# generate log-normal catch
c_realize <- data.frame(haul = 1:haul_num, c_h = exp(stats::rnorm(haul_num, 0, 1))) 

# determine which school is being sampled (based on relative abundance)
sc_smp <- colSums(rmultinom(haul_num, 1, p_sc) * 1:3)

# generate a realization of samples across hauls
rr <- purrr::map(1:haul_num, ~gen_haul_samp(sc_smp[.]))

do.call(mapply, c(list, rr, SIMPLIFY = FALSE))$samp_h %>% 
  tidytable::map_df(., ~as.data.frame(.x), .id = "haul") %>% 
  tidytable::rename(cat = V1,
                    samp_h = V2) -> haul_realize

# compute the realization composition
haul_realize %>% 
  tidytable::left_join(haul_realize %>% 
                         tidytable::summarise(haul_samp = sum(samp_h), .by = haul) %>% 
                         tidytable::mutate(prop_s = haul_samp / sum(haul_samp))) %>% 
  tidytable::left_join(c_realize %>% 
                         tidytable::mutate(prop_c = c_h / sum(c_h))) %>% 
  tidytable::summarise(samp_wtd = sum(samp_h * prop_s * prop_c), # weight samples by sample and catch proportion
                       samp_unwtd = sum(samp_h), # do not weight samples
                       .by = cat) %>% 
  tidytable::mutate(real_p_wtd = samp_wtd / sum(samp_wtd),
                    real_p_unwtd = samp_unwtd / sum(samp_unwtd)) %>% 
  tidytable::select(cat, real_p_wtd, real_p_unwtd) -> p_realize


## define bootstrap functions ----

bs_haul_samp <- function(haul_num,
                         haul_realize,
                         c_realize){
  
  # bootstrap hauls
  bs_h <- data.frame(haul = sample(1:haul_num, haul_num, replace = TRUE)) %>% 
    tidytable::mutate(id = .I)
  
  # get bootstrap haul sample size results
  bs_h %>% 
    tidytable::left_join(haul_realize) %>% 
    tidytable::summarise(haul_samp = sum(samp_h), .by = id) %>% 
    tidytable::mutate(prop_s = haul_samp / sum(haul_samp)) %>% 
    tidytable::rename(haul = id) -> bs_prop_s
  
  bs_prop_s %>% 
    tidytable::summarise(nss = sum(haul_samp)) -> bs_nss
  
  # get catch proportions
  bs_h %>% 
    tidytable::left_join(c_realize) %>% 
    tidytable::mutate(prop_c = c_h / sum(c_h)) %>% 
    tidytable::select(haul = id, c_h, prop_c) -> bs_prop_c
  
  # get bootstrap comp results
  bs_h %>% 
    tidytable::left_join(haul_realize) %>% 
    tidytable::select(haul = id, cat, samp_h) %>% 
    tidytable::left_join(bs_prop_s) %>% 
    tidytable::left_join(bs_prop_c) %>% 
    tidytable::summarise(samp_wtd = sum(samp_h * prop_s * prop_c), # weight samples by sample and catch proportion
                         samp_unwtd = sum(samp_h), # do not weight samples
                         .by = cat) %>% 
    tidytable::mutate(samp_p_wtd = samp_wtd / sum(samp_wtd),
                      samp_p_unwtd = samp_unwtd / sum(samp_unwtd)) %>% 
    tidytable::select(cat, samp_p_wtd, samp_p_unwtd) -> bs_comp_it
  
  list(bs_comp_it = bs_comp_it, bs_nss = bs_nss)
}

sim_bs <- function(haul_num, iters){
  
  # get a sample realization
  
  # generate log-normal catch
  c_realize <- data.frame(haul = 1:haul_num, c_h = exp(stats::rnorm(haul_num, 0, 1))) 
  
  # determine which school is being sampled (based on relative abundance)
  sc_smp <- colSums(rmultinom(haul_num, 1, p_sc) * 1:3)
  
  # generate a realization of samples across hauls
  rr <- purrr::map(1:haul_num, ~gen_haul_samp(sc_smp[.]))
  
  do.call(mapply, c(list, rr, SIMPLIFY = FALSE))$samp_h %>% 
    tidytable::map_df(., ~as.data.frame(.x), .id = "haul") %>% 
    tidytable::rename(cat = V1,
                      samp_h = V2) -> haul_realize
  
  # compute the realization composition
  haul_realize %>% 
    tidytable::left_join(haul_realize %>% 
                           tidytable::summarise(haul_samp = sum(samp_h), .by = haul) %>% 
                           tidytable::mutate(prop_s = haul_samp / sum(haul_samp))) %>% 
    tidytable::left_join(c_realize %>% 
                           tidytable::mutate(prop_c = c_h / sum(c_h))) %>% 
    tidytable::summarise(samp_wtd = sum(samp_h * prop_s * prop_c), # weight samples by sample and catch proportion
                         samp_unwtd = sum(samp_h), # do not weight samples
                         .by = cat) %>% 
    tidytable::mutate(real_p_wtd = samp_wtd / sum(samp_wtd),
                      real_p_unwtd = samp_unwtd / sum(samp_unwtd)) %>% 
    tidytable::select(cat, real_p_wtd, real_p_unwtd) -> p_realize
  
  # bootstrap the realization
  rr <- purrr::map(1:iters, ~bs_haul_samp(haul_num, haul_realize, c_realize))
  
  
  do.call(mapply, c(list, rr, SIMPLIFY = FALSE))$bs_comp_it %>% 
    tidytable::map_df(., ~as.data.frame(.x), .id = "sim") -> res_bs
  
  res_bs %>% 
    tidytable::left_join(p_realize) %>% 
    tidytable::summarise(rss_wtd = sum(real_p_wtd * (1- real_p_wtd)) / sum((samp_p_wtd - real_p_wtd) ^ 2),
                         rss_unwtd = sum(real_p_unwtd * (1- real_p_unwtd)) / sum((samp_p_unwtd - real_p_unwtd) ^ 2),
                         .by = sim) -> rss_res_bs
  
  rss_res_bs %>% 
    tidytable::left_join(do.call(mapply, c(list, rr, SIMPLIFY = FALSE))$bs_nss %>% 
                           tidytable::map_df(., ~as.data.frame(.x), .id = "sim")) %>% 
    tidytable::summarise(iss_wtd = psych::harmonic.mean(rss_wtd, zero = FALSE),
                         iss_unwtd = psych::harmonic.mean(rss_unwtd, zero = FALSE),
                         mean_nss = mean(nss)) -> iss_res_bs
  
  list(iss_res_bs = iss_res_bs)
}

## run bootstrap ----

purrr::map(1:bs_iters, ~sim_bs(haul_num, iters))


















# base way to do exp 1 ----

samp_it_wtd <- NULL
samp_it_unwtd_s <- NULL
samp_it_unwtd_c <- NULL
samp_it_unwtd <- NULL

for(r in 1:iters){
  
  # determine which school is being sampled (based on relative abundance)
  sc_smp <- colSums(rmultinom(haul_num, 1, p_sc) * 1:3)
  
  samp <- NULL
  c <- NULL
  
  for(h in 1:haul_num){
    if(sc_smp[h] == 1){
      # generate haul sample size (90% prob that 20, otherwise 10)
      haul_samp <- if(rbinom(1, 1, 0.9) == 1){20}else{10}
      # generate multinomial sample based on school composition
      samp_h <- rmultinom(1, haul_samp, p_s1)
      # generate random catch based on relative school size in pop'n
      c_h <- pop_size[1] * runif(1, c_lci, c_uci)
    }
    if(sc_smp[h] == 2){
      haul_samp <- if(rbinom(1, 1, 0.9) == 1){20}else{10}
      samp_h <- rmultinom(1, haul_samp, p_s2)
      c_h <- pop_size[2] * runif(1, c_lci, c_uci)
    }
    if(sc_smp[h] == 3){
      haul_samp <- if(rbinom(1, 1, 0.9) == 1){20}else{10}
      samp_h <- rmultinom(1, haul_samp, p_s3)
      c_h <- pop_size[3] * runif(1, c_lci, c_uci)
    }
    samp <- cbind(samp, samp_h)
    c <- c(c, c_h)
  }
  
  prop_c <- c / sum(c)
  prop_s <- colSums(samp) / sum(samp)
  
  samp_p_wtd <- rowSums(samp * prop_s * prop_c) / sum(samp * prop_s * prop_c)
  samp_p_unwtd_s <- rowSums(samp * prop_c) / sum(samp * prop_c)
  samp_p_unwtd_c <- rowSums(samp * prop_s) / sum(samp * prop_s)
  samp_p_unwtd <- rowSums(samp) / sum(samp)
  
  samp_it_wtd <- rbind(samp_it_wtd, samp_p_wtd)
  samp_it_unwtd_s <- rbind(samp_it_unwtd_s, samp_p_unwtd_s)
  samp_it_unwtd_c <- rbind(samp_it_unwtd_c, samp_p_unwtd_c)
  samp_it_unwtd <- rbind(samp_it_unwtd, samp_p_unwtd)
  
}


# plot results
plot(p_true, ylim = c(0, 0.777), ylab = 'Composition')
lines(colMeans(samp_it_wtd), col = 'darkgreen')
lines(colMeans(samp_it_unwtd_c), col = 'darkblue')

rand_row <- sample(1:iters, 1)
lines(samp_it_wtd[rand_row,], col = 'green', lty = 2)
lines(samp_it_unwtd[rand_row,], col = 'blue', lty = 2)

