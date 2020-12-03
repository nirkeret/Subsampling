library(survival)
library(Rcpp)
library(Epi)
library(multipleNCC)
library(msm)
sourceCpp("coxsub.cpp", verbose=TRUE)

information_matrix <- function(beta,weights=NULL,times,status,covariates) {
  if(is.null(weights)) {weights = rep(1,n)}
  expb = exp(beta %*% t(covariates))
  tmp = -rcpp_information_matrix(beta,weights,expb,times,status,covariates)
  tmp = tmp + t(tmp)
  diag(tmp) = diag(tmp)/2
  return(tmp)
}

score_residual_matrix <- function(beta,weights=NULL,times,status,covariates) {
  n <- length(times)
  if(is.null(weights)) {weights = rep(1,n)}
  expb = exp(beta %*% t(covariates))
  c = n - sum(status)
  return(rcpp_score_residual(beta,weights,expb,times,status,covariates,c))
}


set.seed(1)

n = 15000
# q=1000
n_controls = 2
r <- 6
N_samples <- 100

PL_res <- matrix(NA,nrow = N_samples,ncol = r)
q_vec = rep(NA,N_samples)
unif_res <- matrix(NA,nrow = N_samples,ncol = r)
opt_A_res <- matrix(NA,nrow = N_samples,ncol = r)
opt_score_cens_res <- matrix(NA,nrow = N_samples,ncol = r)
ncc_naive_res <- matrix(NA,nrow = N_samples,ncol = r)
ncc_samu_res <- matrix(NA,nrow = N_samples,ncol = r)


b <- c(3,-5,1,-1,1,-3,rep_len(c(-1,1/2),r-6))/10
upper_unif = c(4,4,4,4,4,4,rep_len(c(1,6),r-6))
# upper_unif = c(1,6,2,2,1,6,rep_len(c(1,6),r-6))

for(m in 1:N_samples)
{
  c = rexp(n,0.2)
  X_vec = runif(n*r,0,rep(upper_unif,each = n))
  X = matrix(X_vec,nrow = n,ncol = r)
  #creating correlation:
  # X[,4] = 0.5*X[,2] + 0.5*X[,1] +rnorm(n,0,0.1)
  # X[,5] = X[,1] + rnorm(n,0,1)
  # X[,6] = X[,1] + rnorm(n,1,1.5)
  
  
  linear_comb = X %*% b
  rates = c(0.001,0.05)
  knots = c(0,6)
  obs_rates = exp(linear_comb) %*% rates
  y = apply(obs_rates,1,rpexp,n=1,t=knots)
  v = pmin(c,y)
  delta = v == y
  n_events = sum(delta)
  n_cens <- n - n_events
  q = round(n_events*n_controls)

  q_vec[m] = q
  
  #sorting the times
  ord = order(v)
  v = v[ord]
  delta = delta[ord]
  X = X[ord,]
  
  ## full sample
  fit = coxph(Surv(v,delta) ~ X)
  cens_ind = which(delta == 0)
  PL_res[m,] = coef(fit)
  
  
  #pilot estimate - uniform probabilities
  unif_cont_ind = sample(cens_ind,q,replace = T)
  samp_ind_unif = c(unif_cont_ind,setdiff(1:n,cens_ind))
  cens_weights_unif = length(cens_ind)/ q
  X_samp_unif = X[samp_ind_unif,]
  v_samp_unif = v[samp_ind_unif]
  delta_samp_unif = delta[samp_ind_unif]
  weights_unif = ifelse(delta_samp_unif,1,cens_weights_unif)
  fit_samp_unif = coxph(Surv(v_samp_unif,delta_samp_unif) ~ X_samp_unif,weights = weights_unif,robust = F)
  
  unif_res[m,] <- coef(fit_samp_unif)

  #  NCC

  ncc_samp = ccwc(exit=v,fail = delta,controls = n_controls,silent = T)
  fit_naive_ncc = clogit(ncc_samp$Fail~X[ncc_samp$Map,] + strata(ncc_samp$Set))
  ncc_naive_res[m,] <- coef(fit_naive_ncc)
  sampstat = rep(0,n)
  sampstat[ncc_samp$Map[ncc_samp$Fail==0]] = 1
  sampstat[delta] = 2
  fit_samu_ncc = wpl(Surv(v,delta) ~ .,data.frame(X),m=n_controls,samplestat = sampstat)
  ncc_samu_res[m,] = coef(fit_samu_ncc)

  
  # L optimality probabilities calculation + estimation
  Score_u = score_residual_matrix(unif_res[m,],times = v,status = delta,covariates = X)
  score_resids_sums <- sqrt(rowSums(Score_u^2))
  samp_prob_opt_score_cens <- score_resids_sums/sum(score_resids_sums)
  
  # random sampling with score-optimal weights among censored
  samp_ind_cens = sample(1:n_cens,q,replace = T,prob = samp_prob_opt_score_cens)
  samp_ind_opt = c(cens_ind[samp_ind_cens],setdiff(1:n,cens_ind))
  cens_weights_opt = (1/(samp_prob_opt_score_cens * q))[samp_ind_cens]
  X_samp_opt = X[samp_ind_opt,]
  v_samp_opt = v[samp_ind_opt]
  delta_samp_opt = delta[samp_ind_opt]
  weights_opt = c(cens_weights_opt,rep(1,n_events))
  fit_samp_opt = coxph(Surv(v_samp_opt,delta_samp_opt) ~ X_samp_opt,weights = weights_opt,robust = F)

  opt_score_cens_res[m,] <- coef(fit_samp_opt)
  
  #A optimality  full-data-based / subsample-based
  I_u = information_matrix(unif_res[m,],times = v,status = delta,covariates = X)
  # I_u = information_matrix(unif_res[m,],times = v[samp_ind_unif],status = delta[samp_ind_unif],covariates = X[samp_ind_unif,],weights = weights_unif)

  I_inv <- solve(I_u)
  resids_A <- I_inv %*% t(Score_u)
  resids_A_sums <-  as.vector(sqrt(colSums(resids_A^2)))
  samp_prob_opt_A <- resids_A_sums/sum(resids_A_sums)

  # random sampling with A-optimal probabilities among censored
  samp_ind_cens = sample(1:n_cens,q,replace = T,prob = samp_prob_opt_A)
  samp_ind_opt = c(cens_ind[samp_ind_cens],setdiff(1:n,cens_ind))
  cens_weights_opt = (1/(samp_prob_opt_A * q))[samp_ind_cens]
  X_samp_opt = X[samp_ind_opt,]
  v_samp_opt = v[samp_ind_opt]
  delta_samp_opt = delta[samp_ind_opt]
  weights_opt = c(cens_weights_opt,rep(1,n_events))
  fit_samp_opt = coxph(Surv(v_samp_opt,delta_samp_opt) ~ X_samp_opt,weights = weights_opt,robust = F)

  opt_A_res[m,] <- coef(fit_samp_opt)
  
}

bmat <- matrix(b,nrow = N_samples,ncol = r,byrow = T)
RMSE_real_b = c(mean(sqrt(rowSums((PL_res - bmat)^2))),
               mean(sqrt(rowSums((opt_score_cens_res - bmat)^2))),
               mean(sqrt(rowSums((opt_A_res - bmat)^2))),
               mean(sqrt(rowSums((unif_res - bmat)^2))),
               mean(sqrt(rowSums((ncc_naive_res - bmat)^2))),
              mean(sqrt(rowSums((ncc_samu_res - bmat)^2))))

RMSE_PL_b = c(0,mean(sqrt(rowSums((opt_score_cens_res - PL_res)^2))),
             mean(sqrt(rowSums((opt_A_res - PL_res)^2))),
             mean(sqrt(rowSums((unif_res - PL_res)^2))),
             mean(sqrt(rowSums((ncc_naive_res - PL_res)^2))),
             mean(sqrt(rowSums((ncc_samu_res - PL_res)^2))))

sd(q_vec)
