library(survival)
library(MASS)
library(Rcpp)
library(Epi)
library(multipleNCC)
library(msm)
sourceCpp("cox-subsampling.cpp", verbose=TRUE)

information_score_matrix = function(beta,weights=NULL,times,truncs=NULL,status,covariates,samp_prob=NULL,information_mat=T,score_res_mat=T) {
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


r <- 6 #number of covariates
N_samples <- 500
n_controls = 3 #number of sampled censored observations per failure time


PL_time = score_mat_time  = unif_time1 = unif_time2 = opt_A_time = opt_L_time = ccwc_time =  ncc_naive_time = ncc_samu_time = rep(NA,N_samples)
q_vec = rep(NA,N_samples)
PL_res = unif_res = opt_A_res = opt_L_res= ncc_naive_res = ncc_samu_res = matrix(NA,nrow = N_samples,ncol = r)
coverage = matrix(NA,nrow = N_samples,ncol = 6)

var_PL = var_opt_L = var_opt_A = var_unif = var_NCC_c = var_NCC_S = vector("list",length = N_samples)

cens_rate <- rep(NA,N_samples)

b <- c(3,-5,1,-1,1,-3,rep_len(c(-1,1/2),r-6))/10 #regression coefficients
upper_unif = c(4,4,4,4,4,4,rep_len(c(1,1),r-6))  #covariates upper bound
# upper_unif = c(1,6,2,2,1,6,rep_len(c(1,6),r-6))  #covariates upper bound for setting "B"

for(m in 1:N_samples)
{
  n = 35000
  c = rexp(n,0.2)  #sampling censoring times
  X_vec = runif(n*r,0,rep(upper_unif,each = n))  #sampling covariates
  X = matrix(X_vec,nrow = n,ncol = r)
  #creating correlation:  for setting "C"
  # X[,4] = 0.5*X[,2] + 0.5*X[,1] +rnorm(n,0,0.1)
  # X[,5] = X[,1] + rnorm(n,0,1)
  # X[,6] = X[,1] + rnorm(n,1,1.5)
  
  linear_comb = X %*% b
  rates = c(0.001,0.025)  #baseline hazard rates
  knots = c(0,6)  #baseline hazard change points 
  obs_rates = exp(linear_comb) %*% rates
  y = apply(obs_rates,1,rpexp,n=1,t=knots)  #sampling event times from piecewise exponential
  v = pmin(c,y)
  delta = v == y

  R = runif(n,0,quantile(v,0.75))  #sampling recruitment (truncation) times
  n = 15000
  obs = sample(which(v > R),n)  #sampling n observations out of those we manage to observe (due to truncation)
  v = v[obs]
  R = R[obs]
  X = X[obs,]
  delta = delta[obs]
  
  
  n_events = sum(delta)
  n_cens <- n - n_events
  q = n_events*n_controls  #number of sampled censored observations per failure time
  
  q_vec[m] = q
  
  ## full sample - full partial-likelihood
  tic()
  # fit = coxph(Surv(v,delta) ~ X)   #no truncation
  fit = coxph(Surv(R,v,delta,type = "counting") ~ X) #with truncation
  tmp = toc(quiet = T)
  PL_time[m] = tmp$toc - tmp$tic
  PL_res[m,] = coef(fit)
  var_PL[[m]] = fit$var
  
  cover1 = b <= PL_res[m,] + qnorm(0.975)*sqrt(diag(var_PL[[m]]))
  cover2 = b >= PL_res[m,] - qnorm(0.975)*sqrt(diag(var_PL[[m]]))
  
  coverage[m,1] = sum(cover1 & cover2)
  
  cens_ind = which(delta == 0)
  #pilot estimate - uniform sampling
  tic("Uniform1")
  unif_cont_ind = sample(cens_ind,q,replace = T)  #sampling q censored observations
  samp_ind_unif = c(unif_cont_ind,setdiff(1:n,cens_ind)) #joining sampled censored with all events
  cens_weights_unif = length(cens_ind)/ q
  weights_unif = ifelse(delta,1,cens_weights_unif)
  # fit_samp_unif = coxph(Surv(time=v[samp_ind_unif],event=delta[samp_ind_unif]) ~ X[samp_ind_unif,],weights = weights_unif[samp_ind_unif],robust = F)  #no truncation
  fit_samp_unif = coxph(Surv(R[samp_ind_unif],v[samp_ind_unif],delta[samp_ind_unif],type = "counting") ~ X[samp_ind_unif,],weights = weights_unif[samp_ind_unif],robust = F)  #with truncation
  
  tmp=toc(quiet = T)
  unif_time1[m] = tmp$toc - tmp$tic
  unif_res[m,] = coef(fit_samp_unif)
  
  tic("Uniform2")   #estimating the variance
  # tmpU = information_score_matrix(unif_res[m,],weights = weights_unif[samp_ind_unif],times = v[samp_ind_unif],status = delta[samp_ind_unif],covariates = X[samp_ind_unif,])  # no trucation
  tmpU = information_score_matrix(unif_res[m,],weights = weights_unif[samp_ind_unif],times = v[samp_ind_unif],truncs = R[samp_ind_unif],status = delta[samp_ind_unif],covariates = X[samp_ind_unif,]) #with truncation
  
  Score_U = tmpU$residual * n_cens
  phi_mat_U = cov(Score_U)
  I_inv_U = solve(tmpU$hess)
  var_unif[[m]] = I_inv_U + I_inv_U %*% phi_mat_U %*% I_inv_U/q
  
  tmp=toc(quiet = T)
  unif_time2[m] = tmp$toc - tmp$tic

  cover1 = b <= unif_res[m,] + qnorm(0.975)*sqrt(diag(var_unif[[m]]))
  cover2 = b >= unif_res[m,] - qnorm(0.975)*sqrt(diag(var_unif[[m]]))
  
  coverage[m,4] = sum(cover1 & cover2)

  #  NCC - naive NCC and Samuelsen's NCC

  tic()
  # ncc_samp = ccwc(exit=v,fail = delta,controls = n_controls,silent = T) #no truncation
  ncc_samp = ccwc(entry = R,exit=v,fail = delta,controls = n_controls,silent = T)  #with truncation

  tmp = toc(quiet = T)
  ccwc_time[m] = tmp$toc - tmp$tic
  tic()
  # fit_naive_ncc = coxph(Surv(time=v[ncc_samp$Map],event=ncc_samp$Fail)~X[ncc_samp$Map,] + strata(ncc_samp$Set))  # no truncation
  fit_naive_ncc = coxph(Surv(R[ncc_samp$Map],v[ncc_samp$Map],ncc_samp$Fail,type = "counting")~X[ncc_samp$Map,] + strata(ncc_samp$Set))  #with truncation

  tmp = toc(quiet = T)
  ncc_naive_time[m] = tmp$toc - tmp$tic
  ncc_naive_res[m,] <- coef(fit_naive_ncc)
  var_NCC_c[[m]] = fit_naive_ncc$var
  
  cover1 = b <= ncc_naive_res[m,] + qnorm(0.975)*sqrt(diag(var_NCC_c[[m]]))
  cover2 = b >= ncc_naive_res[m,] - qnorm(0.975)*sqrt(diag(var_NCC_c[[m]]))
  
  coverage[m,5] = sum(cover1 & cover2)
  
  tic()
  sampstat = rep(0,n)
  sampstat[ncc_samp$Map[ncc_samp$Fail==0]] = 1
  sampstat[delta] = 2
  # fit_samu_ncc = wpl(Surv(time=v,event=delta) ~ .,data.frame(X),m=n_controls,samplestat = sampstat)  # no truncation
  fit_samu_ncc = wpl(Surv(R,v,delta,type = "counting") ~ .,data.frame(X),m=n_controls,samplestat = sampstat)  # with truncation
  tmp = toc(quiet = T)
  ncc_samu_time[m] = tmp$toc - tmp$tic
  ncc_samu_res[m,] = coef(fit_samu_ncc)
  var_NCC_S[[m]] = fit_samu_ncc$var
  
  cover1 = b <= ncc_samu_res[m,] + qnorm(0.975)*sqrt(diag(var_NCC_S[[m]]))
  cover2 = b >= ncc_samu_res[m,] - qnorm(0.975)*sqrt(diag(var_NCC_S[[m]]))
  
  coverage[m,6] = sum(cover1 & cover2)

  # L optimality sampling
  tic()

  # tmp_L1 = information_score_matrix(unif_res[m,],times = v,status = delta,covariates = X, samp_prob = "L",information_mat = F)  # no truncation
  tmp_L1 = information_score_matrix(unif_res[m,],times = v,truncs = R,status = delta,covariates = X, samp_prob = "L",information_mat = F)  # with truncation
  delta_L = delta[tmp_L1$ord]
  cens_ind_ord = which(!delta_L)

  # random sampling with L-optimal probabilities from censored
  samp_ind_cens = sample(1:n_cens,q,replace = T,prob = tmp_L1$samp_prob)  #sampling q censored observations
  samp_ind_opt = c(cens_ind_ord[samp_ind_cens],which(delta_L)) #joining the sampled censored with the failure times
  cens_weights_opt = (1/(tmp_L1$samp_prob * q))[samp_ind_cens]
  weights_opt = c(cens_weights_opt,rep(1,n_events))
  
  samp_ord = tmp_L1$ord[samp_ind_opt]
  
  # fit_samp_opt_L = coxph(Surv(time=v[samp_ord],event=delta[samp_ord]) ~ X[samp_ord,],weights = weights_opt,robust = F,init = unif_res[m,])  #no truncation
  fit_samp_opt_L = coxph(Surv(R[samp_ord],v[samp_ord],delta[samp_ord],type = "counting") ~ X[samp_ord,],weights = weights_opt,robust = F,init = unif_res[m,])  #with truncation
  
  opt_L_res[m,] = coef(fit_samp_opt_L)
  
  # estimating the variance
  # tmp_L2 = information_score_matrix(opt_L_res[m,],weights = weights_opt,times = v[samp_ord],status = delta[samp_ord],covariates = X[samp_ord,])  #no truncation
  tmp_L2 = information_score_matrix(opt_L_res[m,],weights = weights_opt,truncs = R[samp_ord],times = v[samp_ord],status = delta[samp_ord],covariates = X[samp_ord,]) #with truncation
  ind_rm = (1:n_events)+q
  order_rm = tmp_L2$ord[!(tmp_L2$ord %in% ind_rm)]
  Score_L = tmp_L2$residual / tmp_L1$samp_prob[samp_ind_cens][order_rm]
  phi_mat_L = cov(Score_L)
  I_inv_L = solve(tmp_L2$hess)
  var_opt_L[[m]] = I_inv_L + I_inv_L %*% phi_mat_L %*% I_inv_L/q
  
  tmp = toc(quiet = T)
  opt_L_time[m] = tmp$toc - tmp$tic
  
  cover1 = b <= opt_L_res[m,] + qnorm(0.975)*sqrt(diag(var_opt_L[[m]]))
  cover2 = b >= opt_L_res[m,] - qnorm(0.975)*sqrt(diag(var_opt_L[[m]]))
  
  coverage[m,2] = sum(cover1 & cover2)
  
  
  #A optimality sampling
  tic()
  # tmp_A1 = information_score_matrix(unif_res[m,],times = v,status = delta,covariates = X, samp_prob = "A")  # no truncation
  tmp_A1 = information_score_matrix(unif_res[m,],times = v,truncs = R,status = delta,covariates = X, samp_prob = "A") # with truncation

  delta_A = delta[tmp_A1$ord]
  cens_ind_ord = which(!delta_A)

  
  # random sampling with A-optimal probabilities from censored
  samp_ind_cens = sample(1:n_cens,q,replace = T,prob = tmp_A1$samp_prob) # sampling q censored observations
  samp_ind_opt = c(cens_ind_ord[samp_ind_cens],which(delta_A)) # joining censored with failure times
  cens_weights_opt = (1/(tmp_A1$samp_prob * q))[samp_ind_cens]
  weights_opt = c(cens_weights_opt,rep(1,n_events))
  
  samp_ord = tmp_A1$ord[samp_ind_opt]
  
  # fit_samp_opt_A = coxph(Surv(time=v[samp_ord],event=delta[samp_ord]) ~ X[samp_ord,],weights = weights_opt,robust = F,init = unif_res[m,])  # no truncation
  fit_samp_opt_A = coxph(Surv(R[samp_ord],v[samp_ord],delta[samp_ord],type = "counting") ~ X[samp_ord,],weights = weights_opt,robust = F,init = unif_res[m,])  # with truncation
  
  opt_A_res[m,] = coef(fit_samp_opt_A)
  
  #estimating the variance
  # tmp_A2 = information_score_matrix(opt_A_res[m,],weights = weights_opt,times = v[samp_ord],status = delta[samp_ord],covariates = X[samp_ord,]) #no truncation
  tmp_A2 = information_score_matrix(opt_A_res[m,],weights = weights_opt,truncs = R[samp_ord],times = v[samp_ord],status = delta[samp_ord],covariates = X[samp_ord,]) #with truncation
  ind_rm = (1:n_events)+q
  order_rm = tmp_A2$ord[!(tmp_A2$ord %in% ind_rm)]
  Score_A = tmp_A2$residual / tmp_A1$samp_prob[samp_ind_cens][order_rm]
  phi_mat_A = cov(Score_A)
  I_inv_A = solve(tmp_A2$hess)
  var_opt_A[[m]] = I_inv_A + I_inv_A %*% phi_mat_A %*% I_inv_A/q

  tmp = toc(quiet = T)
  opt_A_time[m] = tmp$toc - tmp$tic
  
  cover1 = b <= opt_A_res[m,] + qnorm(0.975)*sqrt(diag(var_opt_A[[m]]))
  cover2 = b >= opt_A_res[m,] - qnorm(0.975)*sqrt(diag(var_opt_A[[m]]))
  
  coverage[m,3] = sum(cover1 & cover2)

}


bmat <- matrix(b,nrow = N_samples,ncol = r,byrow = T)
RMSE_real_b = c(mean(sqrt(rowSums((PL_res - bmat)^2))),
               mean(sqrt(rowSums((opt_L_res - bmat)^2))),
               mean(sqrt(rowSums((opt_A_res - bmat)^2))),
               mean(sqrt(rowSums((unif_res - bmat)^2))),
               mean(sqrt(rowSums((ncc_naive_res - bmat)^2))),
              mean(sqrt(rowSums((ncc_samu_res - bmat)^2))))

RMSE_real_b_sd = c(sd(sqrt(rowSums((PL_res - bmat)^2))),
                sd(sqrt(rowSums((opt_L_res - bmat)^2))),
                sd(sqrt(rowSums((opt_A_res - bmat)^2))),
                sd(sqrt(rowSums((unif_res - bmat)^2))),
                sd(sqrt(rowSums((ncc_naive_res - bmat)^2))),
                sd(sqrt(rowSums((ncc_samu_res - bmat)^2))))

rmse_PL = sqrt(rowSums((PL_res - bmat)^2))
RR_real = c(mean((sqrt(rowSums((opt_L_res - bmat)^2)))/rmse_PL),
            mean((sqrt(rowSums((opt_A_res - bmat)^2)))/rmse_PL),
            mean((sqrt(rowSums((unif_res - bmat)^2)))/rmse_PL),
            mean((sqrt(rowSums((ncc_naive_res - bmat)^2)))/rmse_PL),
            mean((sqrt(rowSums((ncc_samu_res - bmat)^2)))/rmse_PL))

RMSE_PL_b = c(0,mean(sqrt(rowSums((opt_L_res - PL_res)^2))),
             mean(sqrt(rowSums((opt_A_res - PL_res)^2))),
             mean(sqrt(rowSums((unif_res - PL_res)^2))),
             mean(sqrt(rowSums((ncc_naive_res - PL_res)^2))),
             mean(sqrt(rowSums((ncc_samu_res - PL_res)^2))))

RMSE_PL_b_sd = c(0,sd(sqrt(rowSums((opt_L_res - PL_res)^2))),
              sd(sqrt(rowSums((opt_A_res - PL_res)^2))),
              sd(sqrt(rowSums((unif_res - PL_res)^2))),
              sd(sqrt(rowSums((ncc_naive_res - PL_res)^2))),
              sd(sqrt(rowSums((ncc_samu_res - PL_res)^2))))
             
time_vec = c(mean(PL_time),
             mean(unif_time1 + opt_L_time),
             mean(unif_time1 + opt_A_time),
             mean(unif_time1 + unif_time2),
             mean(ccwc_time + ncc_naive_time),
             mean(ccwc_time + ncc_samu_time))

time_vec_sd = c(sd(PL_time),
             sd(unif_time1 + opt_L_time),
             sd(unif_time1 + opt_A_time),
             sd(unif_time1 + unif_time2),
             sd(ccwc_time + ncc_naive_time),
             sd(ccwc_time + ncc_samu_time))

mean(q_vec)
sd(q_vec)

D = data.frame("A"=NA,"B"=c(paste0(round(mean(q_vec))," (",round(sd(q_vec)),")"),rep(NA,5)),"Method"=c("PL","L-opt","A-opt","Uniform","NCC-c","NCC-S"),
               "RMSE beta_0" = paste0(format(round(RMSE_real_b,3),nsmall = 3)," (",format(round(RMSE_real_b_sd,3),nsmall=3),")"),
               "RR"=format(round(RMSE_real_b/RMSE_real_b[1],3),nsmall = 3),"RMSE_PL" = paste0(format(round(RMSE_PL_b,3),nsmall=3)," (",format(round(RMSE_PL_b_sd,3),nsmall = 3),")"),
               "CR" = format(round(colSums(coverage) / (r*N_samples),3),nsmall = 3),"Running Time in sec." = paste0(format(round(time_vec,3),nsmall = 3)," (",format(round(time_vec_sd,3),nsmall = 3),")"))
D
print(xtable(D,digits = 3), include.rownames=FALSE)
