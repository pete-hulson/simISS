#' Title Dirichlet Mutlinomial Likelihood
#' From https://github.com/James-Thorson/CCSRA/blob/main/inst/executables/CCSRA_v9.cpp
#' @param obs Vector of observed values in proportions
#' @param pred Vector or predicted values in proportions
#' @param Ntotal Input sample size scalar
#' @param ln_theta Weighting parameter in log space
#' @param give_log Whether or not likelihood is in log space
#'
#' @returns
#' @export
#'
#' @examples
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


#' Title Estimate parameters for a logistic normal with different correlation structures
#'
#' @param data data list (obs, exp, sex_labels (S1, S2), a vector of wts (input sample size))
#' @param start_sigma Starting value for sigma (gets squared in the estimation process)
#' @param start_rho_age Starting value for correlation ages (bound between -1 and 1)
#' @param start_rho_sex Starting value for correlation sexes (bound between -1 and 1)
#' @param cov_strc covariance structure options ("iid", "1DAR1_Across", "2DAR1_Across)
#' @param const_percent constant percentage of minimum value in dataset != 0 to use when there are 0s in data set
#' @param iter_limit Iteration limit for nlm
#' @param gradient_tol Gradient tolerance for when to stop for nlm
#'
#' @return
#' @export
#'
#' @examples
Est_Logistic_Normal <- function(data, 
                                start_sigma = 0.1, 
                                start_rho_age = 0.1,
                                start_rho_sex = 0.1,
                                const_percent = 0.1,
                                cov_strc,
                                iter_limit = 1e7,
                                gradient_tol = 1e-10) {
  
  rho_trans <- function(x) 2/(1+ exp(-2 * x)) - 1
  
  # if the proportions do not sum to 1, renormalize
  if(!sum(rowSums(data$obs)) != ncol(data$obs) ||
     !sum(rowSums(data$exp)) != ncol(data$exp)) {
    data$obs <- t(apply(data$obs, 1, function(x) x / sum(x))) # need to transpose because apply switches it
    data$exp <- t(apply(data$exp, 1, function(x) x / sum(x)))
    print("Normalizing datasets, proportions do not sum to 1")
  } # end if
  
  if(any(data$obs == 0) || any(data$exp == 0) ) {
    eps <- 1e-6
    data$obs[data$obs == 0] <- eps # Add constant in
    data$exp[data$exp == 0] <- eps # Add constant in
    data$obs <- t(apply(data$obs, 1, function(x) x / sum(x))) # renormalize
    data$exp <- t(apply(data$exp, 1, function(x) x / sum(x))) # renormalize
    print("The Logistic-normal cannot accomodate 0s. Adding a small constant based on x percent of the smallest number.")
  } else{
    eps <- 0 # no constant added
  } # end if else
  
  # Set up dimensions
  Nsexes <- length(data$sex_labels)
  
  if(cov_strc == "iid") {
    # Define new likelihood function to minimize
    negloglik <- function(pars, data){
      sigma <- exp(pars[1]) # gets squared later 
      nLL <- nLL_Logistic_Normal(data = data, sigma = sigma, 
                                 rho_s = 0, rho_a = 0,
                                 cov_strc = cov_strc)
      return(nLL)
    } # end likelihood function
    
    # Minimize here
    fit <- nlminb(start = start_sigma, objective = negloglik, 
                  data = data, hessian = TRUE, 
                  control = list(iter.max = iter_limit, eval.max = iter_limit, rel.tol = 1e-10))
    
    # output results
    est_sigma <- exp(fit$par) # get estimated sigma
    mvn_cov <- diag(rep(est_sigma, ncol(data$obs)))  # Compute MVN covariance matrix
    
    out <- list(sigma = est_sigma, MVN_Covariance = mvn_cov, fit = fit, 
                data = data, AIC = get_AIC(fit), eps = eps, cov_strc = cov_strc) # return outputs
  } # covariance = iid
  
  if(cov_strc == "1DAR1_Across") {
    # Define new likelihood function to minimize
    negloglik <- function(pars, data){
      sigma <- exp(pars[1]) # gets squared later 
      rho_age <- rho_trans(pars[2]) # constrain estimation to -1 and 1
      nLL <- nLL_Logistic_Normal(data = data, sigma = sigma, 
                                 rho_a = rho_age, rho_s = 0,
                                 cov_strc = cov_strc)
      return(nLL)
    } # end likelihood function
    
    # Minimize here
    fit <- nlminb(start = c(start_sigma, start_rho_age), objective = negloglik, 
                  data = data, hessian = TRUE, 
                  control = list(iter.max = iter_limit, eval.max = iter_limit, rel.tol = 1e-10))
    
    # Output results
    est_sigma <- exp(fit$par[1]) # get estimated sigma
    est_rho_age <- rho_trans(fit$par[2]) # get estimated correlation
    mvn_cov <- est_sigma^2 * get_AR1_CorrMat(n = ncol(data$obs), rho = est_rho_age) # MVN covariance
    out <- list(sigma = est_sigma, rho_a = est_rho_age, MVN_Covariance = mvn_cov, 
                fit = fit, data = data, AIC = get_AIC(fit), eps = eps, cov_strc = cov_strc)  # return outputs
  } # covariance = 1DAR1 Across
  
  if(cov_strc == "1DAR1_Within") {
    
    # Define new likelihood function to minimize
    negloglik <- function(pars, data){
      sigma <- exp(pars[1]) # gets squared later 
      rho_sex <- rho_trans(0) # constrain estimation to -1 and 1
      rho_age <- rho_trans(pars[2]) # constrain estimation to -1 and 1
      nLL <- nLL_Logistic_Normal(data = data, sigma = sigma, 
                                 rho_s = rho_sex, rho_a = rho_age, 
                                 cov_strc = cov_strc)
      return(nLL)
    } # end likelihood function
    
    # Minimize here
    fit <- nlminb(start = c(start_sigma, start_rho_age), objective = negloglik, 
                  data = data, hessian = TRUE, 
                  control = list(iter.max = iter_limit, eval.max = iter_limit, rel.tol = 1e-10))
    
    # Output results
    est_sigma <- exp(fit$par[1]) # get estimated sigma
    est_rho_age <- rho_trans(fit$par[2]) # get estimated correlation
    mvn_cov <- est_sigma^2 * get_2d_correlation(data = data, rho_a = est_rho_age, rho_s = 0)
    
    # return outputs
    out <- list(sigma = est_sigma, rho_s = 0, rho_a = est_rho_age, 
                MVN_Covariance = mvn_cov, fit = fit, 
                data = data, AIC = get_AIC(fit), eps = eps, cov_strc = cov_strc) 
    
  } # end 1DAR1 Within
  
  if(cov_strc == "2DAR1_Across") {
    # Define new likelihood function to minimize
    negloglik <- function(pars, data){
      sigma <- exp(pars[1]) # gets squared later 
      rho_sex <- rho_trans(pars[2]) # constrain estimation to -1 and 1
      rho_age <- rho_trans(pars[3]) # constrain estimation to -1 and 1
      nLL <- nLL_Logistic_Normal(data = data, sigma = sigma, 
                                 rho_s = rho_sex, rho_a = rho_age, 
                                 cov_strc = cov_strc)
      return(nLL)
    } # end likelihood function
    
    # Minimize here
    fit <- nlminb(start = c(start_sigma, start_rho_age, start_rho_sex), objective = negloglik, 
                  data = data, hessian = TRUE, 
                  control = list(iter.max = iter_limit, eval.max = iter_limit, rel.tol = 1e-10))
    
    # Output results
    est_sigma <- exp(fit$par[1]) # get estimated sigma
    est_rho_sex <- rho_trans(fit$par[2]) # get estimated correlation
    est_rho_age <- rho_trans(fit$par[3]) # get estimated correlation
    mvn_cov <- est_sigma^2 * get_2d_correlation(data = data, rho_a = est_rho_age, rho_s = est_rho_sex)
    
    # return outputs
    out <- list(sigma = est_sigma, rho_s = est_rho_sex, rho_a = est_rho_age, 
                MVN_Covariance = mvn_cov, fit = fit, 
                data = data, AIC = get_AIC(fit), eps = eps, cov_strc = cov_strc) 
  } # covariance = 2DAR1
  
  return(out)
}

#' Title negative log likelihood for a logistic-normal
#'
#' @param data data list (requires a matrix of obs, exp, a vector of sex_labels (S1, S2), a vector of wts (input sample size))
#' @param sigma variance parameter
#' @param rho correlation parameters (rho_s, rho_a)
#' @param cov_strc covariance structure options ("iid", "1DAR1_Across", "2DAR1_Across)
#'
#' @return
#' @export
#'
#' @examples
nLL_Logistic_Normal <- function(data, sigma, rho_s, rho_a, cov_strc) {
  
  # Get observed and predicted data here
  obsmat <- data$obs # observed
  expmat <- data$exp # predicted
  
  # Get dimensions
  Nsexes <- length(data$sex_labels)
  Nbins <- ncol(obsmat)
  Nyears <- nrow(obsmat)
  
  # Inter-annual weighting
  wts <- data$wts
  negloglik <- 0
  
  for(i in 1:Nyears) {
    
    sigma2_y <- sigma^2 / wts[i] # calculate sigma here
    
    # Construct Covariance Matrices here
    if(cov_strc == "iid") covmat <- diag(rep(sigma2_y,Nbins)) # iid matrix
    if(cov_strc %in% c("1DAR1_Across", "1DAR1_Within")) covmat <- sigma2_y * get_AR1_CorrMat(n = Nbins, rho = rho_a) # AR1 within or across
    if(cov_strc == "2DAR1_Across") {
      # Define correlation matrices
      # corr_total <- kronecker(mat_s, ar1_a) # kronecker of the two matrices
      corr_total <- get_2d_correlation(data = data, rho_a = rho_a, rho_s = rho_s ) # get 2d correlation
      covmat <- sigma^2 * corr_total # covariance (multiply by the 2dar1 correlation) 
    } # if 2DAR1
    
    # do logistic transformation on observed values
    tmp_Obs = log(obsmat[i,-ncol(obsmat)])
    tmp_Obs = tmp_Obs - log(obsmat[i,ncol(obsmat)])
    # do logistic transformation on expected values
    mu = log(expmat[i,-ncol(expmat)]) # remove last bin since it's known
    mu = mu - log(expmat[i,ncol(expmat)]) # calculate log ratio
    negloglik <- negloglik - RTMB::dmvnorm(tmp_Obs, mu, covmat[-nrow(covmat), -ncol(covmat)], log = TRUE)
    
  } # end i loop
  
  
  # Now, create K matrix for the MVN covariance
  # Kmat <- cbind(diag(Nbins-1),-1) # diagonal 1s stacked with -1 on the end
  # Vmat <- Kmat %*% (covmat %*% t(Kmat)) # get V matrix for use in nLL
  # Vinv <- solve(Vmat) # inverse of V matrix
  # Vinvdiag <- diag(Vinv) # diagonal of inverse of V matrix
  
  # Create w matrix (log (Ob / OB)) - log( (Eb / EB))
  # ww <- log(sweep(obsmat[,-Nbins,drop=F],1,obsmat[,Nbins],'/')) -
  # log(sweep(expmat[,-Nbins,drop=F],1,expmat[,Nbins],'/'))
  
  # Calculate nLL
  # negloglik <- 0.5*Nyears*(Nbins-1)*log(2*pi) + sum(log(obsmat)) + 0.5*Nyears*log(det(Vmat))
  # for(i in 1:Nyears) negloglik <- negloglik + (0.5/(wts[i]^2))*(ww[i,] %*% Vinv) %*% ww[i,]
  return(negloglik)
}


#' Title Create a correlation matrix that induces two dimensional correlation across ages and sexes
#'
#' @param data data list (requires at least sex_labels, bin width, and column labels (i.e., sex age combination labels))
#' @param rho_a correlation by age
#' @param rho_s correlation by sex
#'
#' @return
#' @export
#'
#' @examples
get_2d_correlation <- function(data, 
                               rho_a, 
                               rho_s) {
  require(Matrix)
  require(readr)
  
  # Get labels
  sex_labels <- data$sex_labels
  bin_width <- data$bin_width
  labels <- data$Labels
  n_sexes <- length(sex_labels) # get number of sexes
  
  # Set up storage objects 
  label_list <- list() # labels
  corr_list <- list() # correlation matrices
  cross_corr_list <- list() # cross correlation matrices
  
  # Loop through to get labels and construct matrices
  for(s in 1:n_sexes) {
    # get labels
    label_list[[s]] <- labels[grep(sex_labels[s], labels)] 
    # Create correlation matrices
    corr_list[[s]] <- matrix(data = 0, 
                             nrow = length(label_list[[s]]), 
                             ncol = length(label_list[[s]]),
                             dimnames = list(label_list[[s]], label_list[[s]]))
  } # end s
  
  names(label_list) <- sex_labels # name label list
  names(corr_list) <- sex_labels # name correlation list
  
  # Get all combinations of sexes to create cross correlation matrices
  comb_indices <- combn(sex_labels, 2)
  cross_comb_labels <- vector() # create labels 
  for(i in 1:ncol(comb_indices)) {
    cross_comb_names <- comb_indices[,i] # extract unique combinations
    cross_comb_labels[i] <- paste(cross_comb_names[1], cross_comb_names[2], sep = "_") # create labels
    # Create cross correlation matrix
    cross_corr_list[[i]] <- matrix(data = 0, 
                                   nrow = length(unlist(label_list[names(label_list) == cross_comb_names[1]])), # Number of rows
                                   ncol = length(unlist(label_list[names(label_list) == cross_comb_names[2]])), # Number of columns
                                   dimnames = list(unlist(label_list[names(label_list) == cross_comb_names[1]]), # row names
                                                   unlist(label_list[names(label_list) == cross_comb_names[2]])) # column names
    )
  } # end i
  
  names(cross_corr_list) <- cross_comb_labels # name cross correlation matrix
  
  
  # Construct correlation via paths -----------------------------------------
  # Construct paths here first
  path_list <- list()
  for(l in 1:length(label_list)) {
    path_vec <- vector() # empty vector to store paths
    sex_vec_tmp <- label_list[[l]] # Extract list
    for(i in 1:length(sex_vec_tmp) - 1) {
      path_vec[i] <- paste(paste(sex_vec_tmp[i], "->", sex_vec_tmp[i+1]), 
                           paste(", " ,  0, sep = ""), sep = "")
    } # end i
    path_list[[l]] <- path_vec
  } # end label list
  
  # Loop through specified paths to construct correlation matrices
  for(s in 1:n_sexes) {
    # Get model path via scan
    model <- scan(text = path_list[[s]], what = list(path = "", offset = 1), 
                  sep = ",", strip.white = TRUE, quiet = TRUE,
                  comment.char = "#", fill = TRUE)
    
    # Clean up model paths
    model$path <- gsub("\\t", " ", model$path) # remove back ticks
    model$par[model$par == ""] <- NA
    model <- cbind(path = model$path, offset = model$offset) # get structure to loop through
    
    for(i in seq_len(nrow(model))) {
      paths <- model[i, 1] # get path
      offset <- as.numeric(model[i, 2]) # get lag from, to (i.e., 39 -> 41, 1 = lag1 of rho^2 between these)
      
      # munge to get path patterns
      vars <- gsub(pattern = " ", replacement = "", x = paths)
      vars <- sub("-*>", "->", sub("<-*", "<-", vars))
      vars <- sub("<->|<-", "->", vars)
      vars <- strsplit(vars, "->")[[1]]    
      
      # Convert variable names to index (should be square matrix)
      from_idx <- which(vars[1] == colnames(corr_list[[s]])) # get from index for correlation matrix
      to_idx <- which(vars[2] == colnames(corr_list[[s]])) # get to index for correlation matrix
      
      # Apply the AR1 structure
      corr_list[[s]][from_idx, to_idx:ncol(corr_list[[s]])] <- rho_a^abs(from_idx - to_idx:ncol(corr_list[[s]]) - offset) # from, to
      corr_list[[s]][to_idx:ncol(corr_list[[s]]), from_idx] <- rho_a^abs(min(to_idx):nrow(corr_list[[s]]) - from_idx - offset) # to, from 
    } # end i loop
    diag(corr_list[[s]]) <- 1 # Set diagonals to 1
  } # end s
  
  
  # Construct cross correlation matrices ------------------------------------
  # Construct cross paths here first
  cross_path_list <- list()
  for(l in 1:ncol(comb_indices)) {
    path_vec <- vector() # empty vector to store paths
    sex_label_list_tmp <- label_list[comb_indices[,l]] # extract labels based on correlation combinations
    # find difference from first bin of a given sex with the first bin of the next sex
    diff <- abs(parse_number(str_remove_all(sex_label_list_tmp[[1]][1], paste(sex_labels, collapse = "|"))) - 
                  parse_number(str_remove_all(sex_label_list_tmp[[2]][1], paste(sex_labels, collapse = "|")))) 
    offset <- diff / bin_width # get offset to lag by
    
    # Get the shorter vector to sequence through
    if(length(sex_label_list_tmp[[1]]) < length(sex_label_list_tmp[[2]])) {
      short_vec_seq <- sex_label_list_tmp[[1]]
    }else short_vec_seq <- sex_label_list_tmp[[2]]
    
    for(i in 1:length(short_vec_seq)) {
      path_vec[i] <- paste(paste(sex_label_list_tmp[[1]][i], "->", sex_label_list_tmp[[2]][i]), 
                           paste(", " ,  offset, sep = ""), sep = "")
    } # end i
    cross_path_list[[l]] <- path_vec
  } # end label list
  
  for(i in 1:length(cross_corr_list)) {
    
    # Get model path via scan
    model_cross <- scan(text = cross_path_list[[i]], 
                        what = list(path = "", offset = 1), 
                        sep = ",", strip.white = TRUE, quiet = TRUE,
                        comment.char = "#", fill = TRUE)
    
    # Clean up model paths
    model_cross$path <- gsub("\\t", " ", model_cross$path) # remove back ticks
    model_cross$par[model_cross$par == ""] <- NA
    model_cross <- cbind(path = model_cross$path, offset = model_cross$offset) # get structure to loop through
    
    for(j in seq_len(nrow(model_cross))) {
      paths <- model_cross[j, 1] # get path
      offset <- as.numeric(model_cross[j, 2]) # get lag from, to (i.e., 39 -> 41, 1 = lag1 of rho^2 between these)
      
      # munge to get path patterns
      vars <- gsub(pattern = " ", replacement = "", x = paths)
      vars <- sub("-*>", "->", sub("<-*", "<-", vars))
      vars <- sub("<->|<-", "->", vars)
      vars <- strsplit(vars, "->")[[1]]    
      
      # Convert variable names to index (should be square matrix)
      from_idx <- which(vars[1] == rownames(cross_corr_list[[i]])) # get from index for correlation matrix
      to_idx <- which(vars[2] == colnames(cross_corr_list[[i]])) # get to index for correlation matrix
      # Because there is potential for a non-square matrix, we need to account for this
      cross_corr_list[[i]][from_idx, to_idx:ncol(cross_corr_list[[i]])] <- rho_s * rho_a^abs(from_idx - to_idx:ncol(cross_corr_list[[i]]) - offset)
      cross_corr_list[[i]][to_idx:nrow(cross_corr_list[[i]]), from_idx] <- rho_s * rho_a^abs(min(to_idx):nrow(cross_corr_list[[i]]) - from_idx - offset)
      
    } # end i loop
  } # end i
  
  # Create joint correlation here
  if(n_sexes == 2) {
    joint_corr <- rbind(cbind(corr_list[[1]], cross_corr_list[[1]]),
                        cbind(t(cross_corr_list[[1]]), corr_list[[2]]))
  } # 2 sexes (collapses to a kronecker when nages are equal across sexes)
  
  if(n_sexes == 3) {
    # Construct entire correlation here
    joint_corr <- rbind(cbind(corr_list[[1]], cross_corr_list[[1]], cross_corr_list[[2]]),
                        cbind(t(cross_corr_list[[1]]), corr_list[[2]], cross_corr_list[[3]]),
                        cbind(t(cross_corr_list[[2]]), t(cross_corr_list[[3]]), corr_list[[3]]))
  } # 3 sexes(collapses to a kronecker when nages are equal across sexes)
  return(joint_corr)
} # end function