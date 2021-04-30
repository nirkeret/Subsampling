library(survival)
library(MASS)
library(Rcpp)
library(Epi)
library(multipleNCC)
library(msm)
sourceCpp("cox-subsampling.cpp", verbose=TRUE)

information_score_matrix <- function(beta,weights=NULL,times,truncs=NULL,status,covariates,samp_prob=NULL,information_mat=T,score_res_mat=T) {
  if(is.null(weights)) {weights = c(-1)}
  if(is.null(truncs)) {truncs = c(-1)}
  if(is.null(samp_prob))
  {
    ret = rcpp_score_wrapper(beta,weights,times,truncs,status,covariates,1*information_mat,1*score_res_mat)
    if(information_mat)
    {
      tmp = ret$hess  
      ret_hess = tmp + t(tmp)
      diag(ret_hess) = diag(tmp)
      ret[["hess"]] = -ret_hess 
    }
    return(ret)
  }
  if(samp_prob == "A")
  {
    ret = rcpp_score_wrapper(beta,weights,times,truncs,status,covariates,1,1)
    tmp = ret$hess  
    ret_hess = tmp + t(tmp)
    diag(ret_hess) = diag(tmp)
    ret[["hess"]] = -ret_hess 
    ret[["samp_prob"]] = rcpp_A_OPT(ret[["residual"]],solve(ret[["hess"]]))
    return(ret)
  }
  if(samp_prob == "L")
  {
    if(information_mat)
    {
      ret = rcpp_score_wrapper(beta,weights,times,truncs,status,covariates,1,1)
      tmp = ret$hess  
      ret_hess = tmp + t(tmp)
      diag(ret_hess) = diag(tmp)
      ret[["hess"]] = -ret_hess 
    } else
    {
      ret = rcpp_score_wrapper(beta,weights,times,truncs,status,covariates,0,1)
    }
    ret[["samp_prob"]] = rcpp_L_OPT(ret[["residual"]])
    return(ret)
  }
}

set.seed(1)

n = 15000
r <- 6
N_samples <- 100
n_controls = 4


q_vec = rep(NA,N_samples)
PL_res <- matrix(NA,nrow = N_samples,ncol = r)
unif_res <- matrix(NA,nrow = N_samples,ncol = r)
opt_A_res <- matrix(NA,nrow = N_samples,ncol = r)
ncc_naive_res <- matrix(NA,nrow = N_samples,ncol = r)
ncc_samu_res <- matrix(NA,nrow = N_samples,ncol = r)
opt_L_res <- matrix(NA,nrow = N_samples,ncol = r)

var_opt_L = vector("list",length = N_samples)
var_opt_A = vector("list",length = N_samples)
var_unif = vector("list",length = N_samples)
var_NCC_c = vector("list",length = N_samples)
var_NCC_S = vector("list",length = N_samples)
var_PL = vector("list",length = N_samples)

cens_rate <- rep(NA,N_samples)

b <- c(3,-5,1,-1,1,-3,rep_len(c(-1,1/2),r-6))/10
upper_unif = c(4,4,4,4,4,4,rep_len(c(1,1),r-6))
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

  R = runif(n,0,0.9*v)
  
  q_vec[m] = q
  
  cens_rate[m] <- 1-mean(delta)
  
  ## full sample
  # fit = coxph(Surv(v,delta) ~ X)  # no truncation
  fit = coxph(Surv(R,v,delta,type = "counting") ~ X)  #with truncation
  PL_res[m,] = coef(fit)
  var_PL[[m]] = fit$var
  
  cens_ind = which(delta == 0)
  #pilot estimate - uniform weights
  unif_cont_ind = sample(cens_ind,q,replace = T)
  samp_ind_unif = c(unif_cont_ind,setdiff(1:n,cens_ind))
  cens_weights_unif = length(cens_ind)/ q
  weights_unif = ifelse(delta,1,cens_weights_unif)
  # fit_samp_unif = coxph(Surv(time=v[samp_ind_unif],event=delta[samp_ind_unif]) ~ X[samp_ind_unif,],weights = weights_unif[samp_ind_unif],robust = F)  #no truncation
  fit_samp_unif = coxph(Surv(R[samp_ind_unif],v[samp_ind_unif],delta[samp_ind_unif],type = "counting") ~ X[samp_ind_unif,],weights = weights_unif[samp_ind_unif],robust = F)  #with truncation
  unif_res[m,] <- coef(fit_samp_unif)
  
  # tmpU = information_score_matrix(unif_res[m,],weights = weights_unif[samp_ind_unif],times = v[samp_ind_unif],status = delta[samp_ind_unif],covariates = X[samp_ind_unif,])  #no truncation
  tmpU = information_score_matrix(unif_res[m,],weights = weights_unif[samp_ind_unif],times = v[samp_ind_unif],truncs = R[samp_ind_unif],status = delta[samp_ind_unif],covariates = X[samp_ind_unif,]) #with truncation
  
  Score_U = tmpU$residual * n_cens
  phi_mat_U = cov(Score_U)
  I_inv_U = solve(tmpU$hess)
  var_unif[[m]] = I_inv_U + I_inv_U %*% phi_mat_U %*% I_inv_U/q


  #  NCC
  # ncc_samp = ccwc(exit=v,fail = delta,controls = n_controls,silent = T) #no truncation
  ncc_samp = ccwc(entry = R,exit=v,fail = delta,controls = n_controls,silent = T) #with truncation
  
  # fit_naive_ncc = coxph(Surv(time=v[ncc_samp$Map],event=ncc_samp$Fail)~X[ncc_samp$Map,] + strata(ncc_samp$Set))  #no truncation
  fit_naive_ncc = coxph(Surv(R[ncc_samp$Map],v[ncc_samp$Map],ncc_samp$Fail,type = "counting")~X[ncc_samp$Map,] + strata(ncc_samp$Set))  #with truncation

  ncc_naive_res[m,] <- coef(fit_naive_ncc)
  var_NCC_c[[m]] = fit_naive_ncc$var

  sampstat = rep(0,n)
  sampstat[ncc_samp$Map[ncc_samp$Fail==0]] = 1
  sampstat[delta] = 2
  
  # fit_samu_ncc = wpl(Surv(time=v,event=delta) ~ .,data.frame(X),m=n_controls,samplestat = sampstat)  #no truncation
  fit_samu_ncc = wpl(Surv(R,v,delta,type = "counting") ~ .,data.frame(X),m=n_controls,samplestat = sampstat)  #with truncation

  ncc_samu_res[m,] = coef(fit_samu_ncc)
  var_NCC_S[[m]] = fit_samu_ncc$var

  # L optimality weights calculation + estimation

  # tmp_L1 = information_score_matrix(unif_res[m,],times = v,status = delta,covariates = X, samp_prob = "L",information_mat = F) #no truncation
  tmp_L1 = information_score_matrix(unif_res[m,],times = v,truncs = R,status = delta,covariates = X, samp_prob = "L",information_mat = F)  #with truncation
  delta_L = delta[tmp_L1$ord]
  cens_ind_ord = which(!delta_L)
  
  # random sampling with score-optimal weights among censored
  samp_ind_cens = sample(1:n_cens,q,replace = T,prob = tmp_L1$samp_prob)
  samp_ind_opt = c(cens_ind_ord[samp_ind_cens],which(delta_L))
  cens_weights_opt = (1/(tmp_L1$samp_prob * q))[samp_ind_cens]
  weights_opt = c(cens_weights_opt,rep(1,n_events))
  
  samp_ord = tmp_L1$ord[samp_ind_opt]
  
  # fit_samp_opt_L = coxph(Surv(time=v[samp_ord],event=delta[samp_ord]) ~ X[samp_ord,],weights = weights_opt,robust = F,init = unif_res[m,])  #no truncation
  fit_samp_opt_L = coxph(Surv(R[samp_ord],v[samp_ord],delta[samp_ord],type = "counting") ~ X[samp_ord,],weights = weights_opt,robust = F,init = unif_res[m,])  #with truncation
  
  opt_L_res[m,] <- coef(fit_samp_opt_L)
  
  # tmp_L2 = information_score_matrix(opt_L_res[m,],weights = weights_opt,times = v[samp_ord],status = delta[samp_ord],covariates = X[samp_ord,]) #no truncation
  tmp_L2 = information_score_matrix(opt_L_res[m,],weights = weights_opt,truncs = R[samp_ord],times = v[samp_ord],status = delta[samp_ord],covariates = X[samp_ord,]) #with truncation
  ind_rm = (1:n_events)+q
  order_rm = tmp_L2$ord[!(tmp_L2$ord %in% ind_rm)]
  Score_L = tmp_L2$residual / tmp_L1$samp_prob[samp_ind_cens][order_rm]
  phi_mat_L = cov(Score_L)
  I_inv_L = solve(tmp_L2$hess)
  var_opt_L[[m]] = I_inv_L + I_inv_L %*% phi_mat_L %*% I_inv_L/q
  
  #A optimality

  # tmp_A1 = information_score_matrix(unif_res[m,],times = v,status = delta,covariates = X, samp_prob = "A")  #no truncation
  tmp_A1 = information_score_matrix(unif_res[m,],times = v,truncs = R,status = delta,covariates = X, samp_prob = "A")  #with truncation

  delta_A = delta[tmp_A1$ord]
  cens_ind_ord = which(!delta_A)

  
  # random sampling with A-optimal weights among censored
  samp_ind_cens = sample(1:n_cens,q,replace = T,prob = tmp_A1$samp_prob)
  samp_ind_opt = c(cens_ind_ord[samp_ind_cens],which(delta_A))
  cens_weights_opt = (1/(tmp_A1$samp_prob * q))[samp_ind_cens]
  weights_opt = c(cens_weights_opt,rep(1,n_events))
  
  samp_ord = tmp_A1$ord[samp_ind_opt]
  
  # fit_samp_opt_A = coxph(Surv(time=v[samp_ord],event=delta[samp_ord]) ~ X[samp_ord,],weights = weights_opt,robust = F,init = unif_res[m,])  #no truncation
  fit_samp_opt_A = coxph(Surv(R[samp_ord],v[samp_ord],delta[samp_ord],type = "counting") ~ X[samp_ord,],weights = weights_opt,robust = F,init = unif_res[m,])  #with truncation
  
  opt_A_res[m,] <- coef(fit_samp_opt_A)
  
  # tmp_A2 = information_score_matrix(opt_A_res[m,],weights = weights_opt,times = v[samp_ord],status = delta[samp_ord],covariates = X[samp_ord,])   #no truncation
  tmp_A2 = information_score_matrix(opt_A_res[m,],weights = weights_opt,truncs = R[samp_ord],times = v[samp_ord],status = delta[samp_ord],covariates = X[samp_ord,]) #with truncation
  ind_rm = (1:n_events)+q
  order_rm = tmp_A2$ord[!(tmp_A2$ord %in% ind_rm)]
  Score_A = tmp_A2$residual / tmp_A1$samp_prob[samp_ind_cens][order_rm]
  phi_mat_A = cov(Score_A)
  I_inv_A = solve(tmp_A2$hess)
  var_opt_A[[m]] = I_inv_A + I_inv_A %*% phi_mat_A %*% I_inv_A/q

  if(m%%10 == 0) print(m)
}

#############
bmat <- matrix(b,nrow = N_samples,ncol = r,byrow = T)
RMSE_real_b = c(mean(sqrt(rowSums((PL_res - bmat)^2))),
               mean(sqrt(rowSums((opt_L_res - bmat)^2))),
               mean(sqrt(rowSums((opt_A_res - bmat)^2))),
               mean(sqrt(rowSums((unif_res - bmat)^2))),
               mean(sqrt(rowSums((ncc_naive_res - bmat)^2))),
              mean(sqrt(rowSums((ncc_samu_res - bmat)^2))))

RMSE_PL_b = c(0,mean(sqrt(rowSums((opt_L_res - PL_res)^2))),
             mean(sqrt(rowSums((opt_A_res - PL_res)^2))),
             mean(sqrt(rowSums((unif_res - PL_res)^2))),
             mean(sqrt(rowSums((ncc_naive_res - PL_res)^2))),
             mean(sqrt(rowSums((ncc_samu_res - PL_res)^2))))
             
mean(q_vec)
sd(q_vec)

D = data.frame("A"=NA,"B"=NA,"C"=NA,"Method"=c("PL","L-opt","A-opt","Uniform","NCC-c","NCC-S"),"RMSE beta_0" = RMSE_real_b,"RE"=RMSE_real_b/RMSE_real_b[1],"RMSE_PL" = RMSE_PL_b)

D

#comparing empirical variance to estimated variance
L_emp = cov(opt_L_res)
L_sub = Reduce("+", var_opt_L) / N_samples

A_emp = cov(opt_A_res)
A_sub = Reduce("+", var_opt_A) / N_samples

U_emp = cov(unif_res)
U_sub = Reduce("+", var_unif) / N_samples

NCC_c_emp = cov(ncc_naive_res)
NCC_c_est = Reduce("+", var_NCC_c) / N_samples

NCC_S_emp = cov(ncc_samu_res)
NCC_S_est = Reduce("+", var_NCC_S) / N_samples

PL_emp = cov(PL_res)
PL_est = Reduce("+", var_PL) / N_samples


round(L_emp,5)
round(L_sub,5)
round(A_emp,5)
round(A_sub,5)
round(U_emp,5)
round(U_sub,5)
round(NCC_c_emp,5)
round(NCC_c_est,5)
round(NCC_S_emp,5)
round(NCC_S_est,5)
round(PL_emp,5)
round(PL_est,5)


sub = function(x,y){x-y}
mean(sapply(mapply(sub,var_unif,var_PL,SIMPLIFY = F),norm,"f"))
mean(sapply(mapply(sub,var_opt_L,var_PL,SIMPLIFY = F),norm,"f"))
mean(sapply(mapply(sub,var_opt_A,var_PL,SIMPLIFY = F),norm,"f"))

