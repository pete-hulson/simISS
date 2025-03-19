
library(tidyverse)



# experiment 1: pop'n comprised of composit schools ----

## set up sim -----
# number of bootstrap iterations
iters <- 1000

# number of hauls sampled
haul_num <- 500

# catch limits on random uniform scalar
c_lci = 0.5
c_uci = 2

# school relative pop'n size
pop_size <- c(240, 60, 15)
p_sc <- pop_size / sum(pop_size)

# school compositions
p_s1 <- c(0.1, 0.8, 0.1, 0, 0, 0, 0, 0)
p_s2 <- c(0, 0.05, 0.1, 0.7, 0.1, 0.05, 0, 0)
p_s3 <- c(0, 0, 0, 0.05, 0.1, 0.7, 0.1, 0.05)

# compute true pop'n composition
p_true <- colSums(rbind(p_s1, p_s2, p_s3) * pop_size) / sum(pop_size)


## define sim functions ----

# function to generate samples within a haul
gen_haul_samp <- function(sc_smp){
  if(sc_smp == 1){
    # generate haul sample size (90% prob that 20, otherwise 10)
    haul_samp <- if(stats::rbinom(1, 1, 0.9) == 1){20}else{10}
    # generate multinomial sample based on school composition
    samp_h <- cbind(1:length(p_s1), stats::rmultinom(1, haul_samp, p_s1))
  }
  if(sc_smp == 2){
    haul_samp <- if(stats::rbinom(1, 1, 0.9) == 1){20}else{10}
    samp_h <- cbind(1:length(p_s2), stats::rmultinom(1, haul_samp, p_s2))
  }
  if(sc_smp == 3){
    haul_samp <- if(stats::rbinom(1, 1, 0.9) == 1){20}else{10}
    samp_h <- cbind(1:length(p_s3), stats::rmultinom(1, haul_samp, p_s3))
  }
  list(samp_h = samp_h)
}


# funtion to generate and compute composition
sim_comp <- function(haul_num, p_sc){
  
  # generate log-normal catch
  data.frame(haul = 1:haul_num, c_h = exp(stats::rnorm(haul_num, 0, 1))) %>% 
    tidytable::mutate(prop_c = c_h / sum(c_h)) -> prop_c  
  
  # determine which school is being sampled (based on relative abundance)
  sc_smp <- colSums(rmultinom(haul_num, 1, p_sc) * 1:3)
  
  # generate samples across hauls
  rr <- purrr::map(1:haul_num, ~gen_haul_samp(sc_smp[.]))
  
  # get haul sample size results
  do.call(mapply, c(list, rr, SIMPLIFY = FALSE))$samp_h %>% 
    tidytable::map_df(., ~as.data.frame(.x), .id = "haul") %>% 
    tidytable::summarise(haul_samp = sum(V2), .by = haul) %>% 
    tidytable::mutate(prop_s = haul_samp / sum(haul_samp)) -> prop_s
  
  prop_s %>% 
    tidytable::summarise(nss = sum(haul_samp)) -> nss
  
  # get generated comp results
  do.call(mapply, c(list, rr, SIMPLIFY = FALSE))$samp_h %>% 
    tidytable::map_df(., ~as.data.frame(.x), .id = "haul") %>% 
    tidytable::rename(cat = V1,
                      samp_h = V2) %>% 
    tidytable::left_join(prop_s) %>% 
    tidytable::left_join(prop_c) %>% 
    tidytable::summarise(samp_wtd = sum(samp_h * prop_s * prop_c), # weight samples by sample and catch proportion
                         samp_unwtd = sum(samp_h), # do not weight samples
                         .by = cat) %>% 
    tidytable::mutate(samp_p_wtd = samp_wtd / sum(samp_wtd),
                      samp_p_unwtd = samp_unwtd / sum(samp_unwtd)) %>% 
    tidytable::select(cat, samp_p_wtd, samp_p_unwtd) -> comp_it
  
  list(comp_it = comp_it, nss = nss)
  
}


## run sim loop ----
rr <- purrr::map(1:iters, ~sim_comp(haul_num, p_sc))


do.call(mapply, c(list, rr, SIMPLIFY = FALSE))$comp_it %>% 
  tidytable::map_df(., ~as.data.frame(.x), .id = "sim") -> res


## plot results ----

# set up plot data
rand_row <- sample(1:iters, 1)

data.frame(cat = 1:length(p_true), name = 'true', value = p_true) %>% 
  tidytable::bind_rows(res %>% 
                         tidytable::summarise(mean_wtd = mean(samp_p_wtd),
                                              mean_unwtd = mean(samp_p_unwtd), .by = cat) %>% 
                         tidytable::pivot_longer(cols = c(mean_wtd, mean_unwtd)) %>% 
                         tidytable::bind_rows(res %>% 
                                                tidytable::filter(sim == rand_row) %>% 
                                                tidytable::select(-sim) %>% 
                                                tidytable::rename(rand_wtd = samp_p_wtd,
                                                                  rand_unwtd = samp_p_unwtd) %>% 
                                                tidytable::pivot_longer(cols = c(rand_wtd, rand_unwtd)))) %>% 
  tidytable::rename(comp_type = name) %>% 
  tidytable::mutate(lty = tidytable::case_when(comp_type == 'mean_unwtd' ~ 'solid',
                                               comp_type == 'mean_wtd' ~ 'dashed',
                                               comp_type == 'rand_unwtd' ~ 'solid',
                                               comp_type == 'rand_wtd' ~ 'dashed')) -> plot_dat

# plot
ggplot(data = plot_dat %>% tidytable::filter(comp_type != 'true'), 
       aes(x = cat, y = value, col = comp_type)) +
  geom_line(linewidth = 0.75, aes(linetype = lty)) +
  geom_point(data = plot_dat %>% tidytable::filter(comp_type == 'true'), size = 2.5) +
  theme_bw() +
  scico::scale_color_scico_d(palette = 'roma') +
  guides(linetype = 'none')

# so long as the schools are represented across the haul samples in proportion to their relative abundance, 
# the mean across iterations in unbiased, whether weighted or not, so,
# you don't need to weight by hauls catch even if the number of samples within a haul aren't proportional to catch
# however, for any given iteration, the composition resulting from weighted or unweighted comps can be different

## compute the input sample size ----


res %>% 
  tidytable::left_join(data.frame(cat = 1:length(p_true), p_true = p_true)) %>% 
  tidytable::summarise(rss_wtd = sum(p_true * (1- p_true)) / sum((samp_p_wtd - p_true) ^ 2),
                       rss_unwtd = sum(p_true * (1- p_true)) / sum((samp_p_unwtd - p_true) ^ 2),
                       .by = sim) -> rss_res

rss_res %>% 
  tidytable::left_join(do.call(mapply, c(list, rr, SIMPLIFY = FALSE))$nss %>% 
                         tidytable::map_df(., ~as.data.frame(.x), .id = "sim")) %>% 
  tidytable::summarise(iss_wtd = psych::harmonic.mean(rss_wtd, zero = FALSE),
                       iss_unwtd = psych::harmonic.mean(rss_unwtd, zero = FALSE),
                       mean_nss = mean(nss)) -> iss_res

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

