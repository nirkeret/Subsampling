library(survival)
library(Rcpp)
library(Epi)
library(multipleNCC)
library(msm)
sourceCpp("coxsub-LT.cpp", verbose=TRUE)

information_matrix <- function(beta,weights=NULL,times,truncs,status,covariates) {
  if(is.null(weights)) {weights = rep(1,n)}
  R  = c(rep(0,length(times)),rep(1,length(truncs)))
  expb = exp(beta %*% t(covariates))
  weights = c(weights,weights)
  times = c(times,truncs)
  status = c(status,status)
  covariates = as.matrix(rbind(covariates,covariates))
  expb = c(expb,expb)
  ord = order(times)
  weights = weights[ord]
  times = times[ord]
  status = status[ord]
  covariates = covariates[ord,]
  R = R[ord]
  expb = expb[ord]
  tmp = -rcpp_information_matrix(beta,weights,expb,times,R,status,covariates)
  tmp = tmp + t(tmp)
  diag(tmp) = diag(tmp)/2
  return(tmp)
}

score_residual_matrix <- function(beta,weights=NULL,times,truncs,status,covariates) {
  n <- length(times)
  if(is.null(weights)) {weights = rep(1,n)}
  expb = exp(beta %*% t(covariates))
  c = n - sum(status)
  R  = c(rep(0,length(times)),rep(1,length(truncs)))
  weights = c(weights,weights)
  times = c(times,truncs)
  status = c(status,status)
  covariates = as.matrix(rbind(covariates,covariates))
  expb = c(expb,expb)
  ord = order(times)
  rn = rank(times)
  
  weights = weights[ord]
  times = times[ord]
  status = status[ord]
  covariates = covariates[ord,]
  R = R[ord]
  expb = expb[ord]
  RL = rn[(n+1):(2*n)]
  RL = c(RL,rep(0,n))
  RL = RL[ord]
  return(-rcpp_score_residual(beta,weights,expb,times,R,RL - 1,status,covariates,c))
}


set.seed(1)

n = 15000
# q=1000
n_controls = 2
r <- 6
N_samples <- 100
results_table <- NULL

q_vec = rep(NA,N_samples)
PL_res <- matrix(NA,nrow = N_samples,ncol = r)
unif_res <- matrix(NA,nrow = N_samples,ncol = r)
opt_A_res <- matrix(NA,nrow = N_samples,ncol = r)
ncc_naive_res <- matrix(NA,nrow = N_samples,ncol = r)
ncc_samu_res <- matrix(NA,nrow = N_samples,ncol = r)
opt_score_cens_res <- matrix(NA,nrow = N_samples,ncol = r)

b <- c(3,-5,1,-1,1,-3,rep_len(c(-1,1/2),r-6))/10
# upper_unif = c(1,6,2,2,1,6,rep_len(c(1,6),r-6))
upper_unif = c(4,4,4,4,4,4,rep_len(c(1,6),r-6))


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
  R = runif(n,0,0.9*v)
  n_events = sum(delta)
  n_cens <- n - n_events
  q = round(n_events*n_controls)
  
  q_vec[m] = q
  
  #sorting the times
  ord = order(v)
  v = v[ord]
  delta = delta[ord]
  X = X[ord,]
  R = R[ord]
  
  ## full sample
  fit = coxph(Surv(R,v,delta,type = "counting") ~ X)
  PL_time[m] = tmp$toc - tmp$tic
  cens_ind = which(delta == 0)
  PL_res[m,] = coef(fit)
 
  #pilot estimate - uniform probabilities
  
  unif_cont_ind = sample(cens_ind,q,replace = T)
  samp_ind_unif = c(unif_cont_ind,setdiff(1:n,cens_ind))
  cens_weights_unif = length(cens_ind)/ q
  X_samp_unif = X[samp_ind_unif,]
  v_samp_unif = v[samp_ind_unif]
  R_samp_unif = R[samp_ind_unif]
  delta_samp_unif = delta[samp_ind_unif]
  weights_unif = ifelse(delta_samp_unif,1,cens_weights_unif)
  fit_samp_unif = coxph(Surv(R_samp_unif,v_samp_unif,delta_samp_unif,type = "counting") ~ X_samp_unif,weights = weights_unif,robust = F)
  
  unif_res[m,] <- coef(fit_samp_unif)

  #  NCC
  
  ncc_samp = ccwc(entry =R ,exit=v,fail = delta,controls = n_controls,silent = T)

  sampstat = rep(0,n)
  sampstat[ncc_samp$Map[ncc_samp$Fail==0]] = 1
  sampstat[delta] = 2
  fit_samu_ncc = wpl(Surv(R,v,delta,type = "counting") ~ .,data.frame(X),m=n_controls,samplestat = sampstat)
  ncc_samu_res[m,] = coef(fit_samu_ncc)
  
  # L optimality probabilities calculation + estimation
  Score_u = score_residual_matrix(unif_res[m,],times = v,truncs = R, status = delta,covariates = X)
  score_mat_time[m] = tmp$toc - tmp$tic
  
  score_resids_sums <- sqrt(rowSums(Score_u^2))
  samp_prob_opt_score_cens <- score_resids_sums/sum(score_resids_sums)
  
  # random sampling with score-optimal probabilities among censored
  samp_ind_cens = sample(1:n_cens,q,replace = T,prob = samp_prob_opt_score_cens)
  samp_ind_opt = c(cens_ind[samp_ind_cens],setdiff(1:n,cens_ind))
  cens_weights_opt = (1/(samp_prob_opt_score_cens * q))[samp_ind_cens]
  X_samp_opt = X[samp_ind_opt,]
  v_samp_opt = v[samp_ind_opt]
  R_samp_opt = R[samp_ind_opt]
  delta_samp_opt = delta[samp_ind_opt]
  weights_opt = c(cens_weights_opt,rep(1,n_events))
  fit_samp_opt = coxph(Surv(R_samp_opt,v_samp_opt,delta_samp_opt,type = "counting") ~ X_samp_opt,weights = weights_opt,robust = F,init = unif_res[m,])
  
  opt_score_cens_res[m,] <- coef(fit_samp_opt)
  
  
  #A optimalityf  ull-data-based / subsample-based

  # I_u = information_matrix(unif_res[m,],times = v,truncs = R,status = delta,covariates = X)
  I_u = information_matrix(unif_res[m,],times = v_samp_unif,truncs =R_samp_unif, status = delta_samp_unif,covariates = X_samp_unif,weights = weights_unif)

  I_inv <- solve(I_u)
  resids_A <- I_inv %*% t(Score_u)
  resids_A_sums <-  as.vector(sqrt(colSums(resids_A^2)))
  samp_prob_opt_A <- resids_A_sums/sum(resids_A_sums)

  # random sampling with A-optimal weights among censored
  samp_ind_cens = sample(1:n_cens,q,replace = T,prob = samp_prob_opt_A)
  samp_ind_opt = c(cens_ind[samp_ind_cens],setdiff(1:n,cens_ind))
  cens_weights_opt = (1/(samp_prob_opt_A * q))[samp_ind_cens]
  X_samp_opt = X[samp_ind_opt,]
  v_samp_opt = v[samp_ind_opt]
  R_samp_opt = R[samp_ind_opt]
  delta_samp_opt = delta[samp_ind_opt]
  weights_opt = c(cens_weights_opt,rep(1,n_events))
  fit_samp_opt = coxph(Surv(R_samp_opt,v_samp_opt,delta_samp_opt,type = "counting") ~ X_samp_opt,weights = weights_opt,robust = F,init = unif_res[m,])

  opt_A_res[m,] <- coef(fit_samp_opt)
  
}


bmat <- matrix(b,nrow = N_samples,ncol = r,byrow = T)
RMSE_real_b = c(mean(sqrt(rowSums((PL_res - bmat)^2))),
               mean(sqrt(rowSums((opt_score_cens_res - bmat)^2))),
               mean(sqrt(rowSums((opt_A_res - bmat)^2))),
               mean(sqrt(rowSums((unif_res - bmat)^2))),
              mean(sqrt(rowSums((ncc_samu_res - bmat)^2))))

RMSE_PL_b = c(0,mean(sqrt(rowSums((opt_score_cens_res - PL_res)^2))),
             mean(sqrt(rowSums((opt_A_res - PL_res)^2))),
             mean(sqrt(rowSums((unif_res - PL_res)^2))),
             mean(sqrt(rowSums((ncc_samu_res - PL_res)^2))))
             
mean(q_vec)
sd(q_vec)

