
# load/source stuff ----
library(tidyverse)
library(future)
source(here::here('R', 'functions.R'))


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
pu <- 25

# number of population categories (e.g., ages or length bins)
pc <- 15

# CV around mean within pop'n unit (e.g., CV in mean age of school)
pu_cv <- 0.25

# pop'n exponential decay (e.g., lnM, tied to inverse of number of pop'n categories so that pop'n = 0.01 at largest category)
d <- log(0.01) / (1 - pc)


# experiment 3: can a given realization of a sample give the 'true' iss back? ----

test_bs <- function(bs_iters, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters){
  #start timer
  tictoc::tic()
  ## run simulation ----
  rr_bs <- purrr::map(1:bs_iters, ~bs_sim(d, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters))
  # end timer
  runtime <- tictoc::toc()
  
  ## unlist & save results ----
  do.call(mapply, c(list, rr_bs, SIMPLIFY = FALSE))$rss_se %>% 
    tidytable::map_df(., ~as.data.frame(.x), .id = "sim") -> rss_se
  do.call(mapply, c(list, rr_bs, SIMPLIFY = FALSE))$iss_bs %>% 
    tidytable::map_df(., ~as.data.frame(.x), .id = "sim") -> iss_bs
  do.call(mapply, c(list, rr_bs, SIMPLIFY = FALSE))$popn_strctr %>% 
    tidytable::map_df(., ~as.data.frame(.x), .id = "sim") -> popn_strctr
  
  res <- list(rss_se = rss_se, iss_bs = iss_bs, popn_strctr = popn_strctr)
  
  res
}


## run tests in parallel ----

# number of bootstrap replicates
bs_iters <- 10

# set up parallel fcn
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
  
  res
}

# run in parallel
tictoc::tic()
future::plan(multisession, workers = 10)
res_bs <- run_bs_test(bs_iters, pu, pc, pu_cv, su_num, su_samp, p_su_samp, iters)
future::plan(sequential)
runtime_test <- tictoc::toc()

# calc runtime test
(runtime_test$toc - runtime_test$tic) / (60 * bs_iters) * 500 / 60

# save results
saveRDS(res_bs,
        file = here::here('output', paste0('res_bs.rds')))

## plot results ----

true_iss <- res_bs$rss_se %>%  
  tidytable::left_join(res_bs$popn_strctr) %>% 
  tidytable::summarise(Wtd = psych::harmonic.mean(rss_wtd, zero = FALSE),
                       Unwtd = psych::harmonic.mean(rss_unwtd, zero = FALSE),
                       .by = c(selex_type, popn_strctr)) %>% 
  tidytable::bind_rows(res_bs$rss_se %>%  
                         tidytable::left_join(res_bs$popn_strctr) %>% 
                         tidytable::mutate(popn_strctr = 'combined') %>% 
                         tidytable::summarise(Wtd = psych::harmonic.mean(rss_wtd, zero = FALSE),
                                              Unwtd = psych::harmonic.mean(rss_unwtd, zero = FALSE),
                                              .by = c(selex_type, popn_strctr))) %>% 
  tidytable::pivot_longer(cols = c(Wtd, Unwtd)) %>% 
  tidytable::rename(type = name, iss = value) %>% 
  tidytable::mutate(popn_strctr = factor(popn_strctr, levels = c('recruitment pulse', 'multimodal', 'unimodal', 'combined')))

bs_iss <- res_bs$iss_bs %>%  
  tidytable::left_join(res_bs$popn_strctr) %>%
  tidytable::rename(Wtd = bs_iss_wtd, Unwtd = bs_iss_unwtd) %>% 
  tidytable::bind_rows(res_bs$iss_bs %>%  
                         tidytable::left_join(res_bs$popn_strctr) %>%
                         tidytable::mutate(popn_strctr = 'combined') %>% 
                         tidytable::rename(Wtd = bs_iss_wtd, Unwtd = bs_iss_unwtd)) %>% 
  tidytable::pivot_longer(cols = c(Wtd, Unwtd)) %>% 
  tidytable::rename(type = name, iss = value) %>% 
  tidytable::mutate(popn_strctr = factor(popn_strctr, levels = c('recruitment pulse', 'multimodal', 'unimodal', 'combined')))

bs_plot_popn <- ggplot(data = bs_iss, aes(x = type, y = iss, fill = type)) + 
  geom_boxplot() +
  geom_point(data = true_iss, shape = 24, size = 2, fill = 'white', aes(x = type, y = iss)) +
  facet_grid(selex_type ~ popn_strctr, scales = 'free_y') +
  scico::scale_color_scico_d(palette = 'roma') +
  scico::scale_fill_scico_d(palette = 'roma') +
  theme_bw() +
  guides(fill = 'none') +
  ylab("ISS") +
  xlab("Composition expansion type")

ggsave(filename = "bs_test.png",
       plot = bs_plot_popn,
       path = here::here("figs"),
       width = 6.5,
       height = 5,
       units = "in")
