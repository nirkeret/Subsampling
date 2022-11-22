library(survival)
library(Rcpp)

sourceCpp("cox-subsampling.cpp", verbose=TRUE)


# function that runs internally in the main function
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

uniform_sampling = function(n, q0, cens_ind) 
  {
  unif_cont_ind = sample(cens_ind,q0,replace = T)  #sampling q0 censored observations
  samp_ind_unif = c(unif_cont_ind,setdiff(1:n,cens_ind)) #joining sampled censored with all events
  cens_weights_unif = length(cens_ind)/ q0
  weights_unif = ifelse(D,1,cens_weights_unif)
  return(list("samp" = samp_ind_unif, "weights" = weights_unif))
}

calc_uniform_variance = function(tmpU, n_cens, q0) 
{
  Score_U = tmpU$residual * n_cens
  phi_mat_U = cov(Score_U)
  I_inv_U = solve(tmpU$hess)
  var_unif = I_inv_U + I_inv_U %*% phi_mat_U %*% I_inv_U/q0
  return(var_unif)
}

random_sampling = function(n_cens, n_events, q, tmp1, cens_ind_ord, D_ord)
{
  samp_ind_cens = sample(1:n_cens,q,replace = T,prob = tmp1$samp_prob)  #sampling q censored observations
  samp_ind_opt = c(cens_ind_ord[samp_ind_cens],which(D_ord)) #joining the sampled censored with the failure times
  cens_weights_opt = (1/(tmp1$samp_prob * q))[samp_ind_cens]
  weights_opt = c(cens_weights_opt,rep(1,n_events))
  samp_ord = tmp1$ord[samp_ind_opt]
  return(list("samp_ord" = samp_ord, "samp" = samp_ind_cens, "weights" = weights_opt))
}

get_sampling = function(method, U_coef, V, R, D, X) 
{
  get_info_matrix = method == "A"
  if(is.null(R))
  {
    tmp1 = information_score_matrix(U_coef,times = V,status = D,covariates = X, samp_prob = method,information_mat = get_info_matrix)  
  }else
  {
    tmp1 = information_score_matrix(U_coef,times = V,truncs = R,status = D,covariates = X, samp_prob = method,information_mat = get_info_matrix)  
  }
  D_ord = D[tmp1$ord]
  cens_ind_ord = which(!D_ord)
  
  # random sampling with L/A-optimal probabilities from censored
  rand_sampling = random_sampling(n_cens, n_events, q, tmp1, cens_ind_ord, D_ord)
  samp_ord = rand_sampling$samp_ord
  
  if(is.null(R))
  {
    fit_samp_opt = coxph(Surv(time=V[samp_ord],event=D[samp_ord]) ~ X[samp_ord,],weights = rand_sampling$weights,robust = F,init = U_coef)  
    tmp2 = information_score_matrix(coef(fit_samp_opt),weights = rand_sampling$weights,times = V[samp_ord],status = D[samp_ord],covariates = X[samp_ord,]) 
  }else
  {
    fit_samp_opt = coxph(Surv(R[samp_ord],V[samp_ord],D[samp_ord],type = "counting") ~ X[samp_ord,],weights = rand_sampling$weights,robust = F,init = U_coef)  
    tmp2 = information_score_matrix(coef(fit_samp_opt),weights = rand_sampling$weights,truncs = R[samp_ord],times = V[samp_ord],status = D[samp_ord],covariates = X[samp_ord,]) 
  }
  opt_coef = coef(fit_samp_opt)
  names(opt_coef) = colnames(X)
  ind_rm = (1:n_events)+q
  order_rm = tmp2$ord[!(tmp2$ord %in% ind_rm)]
  Score = tmp2$residual / tmp1$samp_prob[rand_sampling$samp][order_rm]
  phi_mat = cov(Score)
  I_inv = solve(tmp2$hess)
  var_opt = I_inv + I_inv %*% phi_mat %*% I_inv/q
  ret = list("coef" = opt_coef, "var" = var_opt)
  return(ret)
}

# V = observed times
# D = status (T/F)
# X = covariate matrix
# R = recruitment (left truncation) times
# q0 = number of subsampled censored observations for the uniform estimator (being the pilot estimator for L and A). Defaults to be the same as q.
# q = number of subsampled censored observations in the second round for the L and A 
# method: one of "U"/"L"/"A" (standing for uniform, L-optimal, A-optimal)

subsampling_cox = function(V, D, X, R=NULL, q, q0 = q, method)
{
  #input validation
  if(!(method %in% c("U","L","A"))) 
  {
    stop("method has to be one of U, L or A")
  }
  if((method == "U") & (q != q0))
  {
    stop("when method set to U, q0 must be equal to q")
  }
  
  n = length(V)
  cens_ind = which(D == 0)
  n_cens = length(cens_ind)
  n_events = n - n_cens
  
  #uniform sampling
  uni_samp = uniform_sampling(n, q0, cens_ind)
  
  if(is.null(R))
  {
    fit_samp_unif = coxph(Surv(time = V[uni_samp$samp], event = D[uni_samp$samp]) ~ X[uni_samp$samp,],weights = uni_samp$weights[uni_samp$samp],robust = F)  
  }else
  {
    fit_samp_unif = coxph(Surv(R[uni_samp$samp],v[uni_samp$samp],delta[uni_samp$samp],type = "counting") ~ X[uni_samp$samp,],weights = uni_samp$weights[uni_samp$samp],robust = F)  
  }
  U_coef = coef(fit_samp_unif)
  names(U_coef) = colnames(X)
  
  if(method == "U")
  {
    if(is.null(R))
    {
      tmpU = information_score_matrix(U_coef,weights = weights_unif[samp_ind_unif],times = v[samp_ind_unif],status = delta[samp_ind_unif],covariates = X[samp_ind_unif,])  
    }else
    {
      tmpU = information_score_matrix(U_coef,weights = weights_unif[samp_ind_unif],times = v[samp_ind_unif],truncs = R[samp_ind_unif],status = delta[samp_ind_unif],covariates = X[samp_ind_unif,]) 
    }
    
    var_unif = calc_uniform_variance(tmpU, n_cens, q0)
    ret = list("coef" = U_coef, "var" = var_unif)
    return(ret)
  }
  else
  {
    sampling = get_sampling(method, U_coef, V, R, D, X)
    return(sampling)
  }
}

