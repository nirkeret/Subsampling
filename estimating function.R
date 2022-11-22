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

# V = observed times
# D = status (T/F)
# X = covariate matrix
# R = recruitment (left truncation) times
# q0 = number of subsampled censored observations for the uniform estimator (being the pilot estimator for L and A). Defaults to be the same as q.
# q = number of subsampled censored observations in the second round for the L and A 
# method: one of "U"/"L"/"A" (standing for uniform, L-optimal, A-optimal)

subsampling_cox = function(V, D, X, R=NULL, q, q0 = q, method)
{
  if(!(method %in% c("U","L","A"))) 
  {
    stop("method has to be one of U, L or A")
  }
  ## add check if method == U then q == q0
  n = length(V)
  cens_ind = which(D == 0)
  n_cens = length(cens_ind)
  n_events = n - n_cens
  
  #uniform sampling
  unif_cont_ind = sample(cens_ind,q0,replace = T)  #sampling q0 censored observations
  samp_ind_unif = c(unif_cont_ind,setdiff(1:n,cens_ind)) #joining sampled censored with all events
  cens_weights_unif = length(cens_ind)/ q0
  weights_unif = ifelse(D,1,cens_weights_unif)
  
  if(is.null(R))
  {
    fit_samp_unif = coxph(Surv(time = V[samp_ind_unif], event = D[samp_ind_unif]) ~ X[samp_ind_unif,],weights = weights_unif[samp_ind_unif],robust = F)  
  }else
  {
    fit_samp_unif = coxph(Surv(R[samp_ind_unif],v[samp_ind_unif],delta[samp_ind_unif],type = "counting") ~ X[samp_ind_unif,],weights = weights_unif[samp_ind_unif],robust = F)  
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
    
    ## variance calc
    Score_U = tmpU$residual * n_cens
    phi_mat_U = cov(Score_U)
    I_inv_U = solve(tmpU$hess)
    var_unif = I_inv_U + I_inv_U %*% phi_mat_U %*% I_inv_U/q0
    ret = list("coef" = U_coef, "var" = var_unif)
    return(ret)
  }
  if(method == "L")
  {
    if(is.null(R))
    {
      tmp_L1 = information_score_matrix(U_coef,times = V,status = D,covariates = X, samp_prob = "L",information_mat = F)  
    }else
    {
      tmp_L1 = information_score_matrix(U_coef,times = V,truncs = R,status = D,covariates = X, samp_prob = "L",information_mat = F)  
    }
    D_L = D[tmp_L1$ord]
    cens_ind_ord = which(!D_L)
    
    # random sampling with L-optimal probabilities from censored
    samp_ind_cens = sample(1:n_cens,q,replace = T,prob = tmp_L1$samp_prob)  #sampling q censored observations
    samp_ind_opt = c(cens_ind_ord[samp_ind_cens],which(D_L)) #joining the sampled censored with the failure times
    cens_weights_opt = (1/(tmp_L1$samp_prob * q))[samp_ind_cens]
    weights_opt = c(cens_weights_opt,rep(1,n_events))
    
    samp_ord = tmp_L1$ord[samp_ind_opt]
    
    if(is.null(R))
    {
      fit_samp_opt_L = coxph(Surv(time=V[samp_ord],event=D[samp_ord]) ~ X[samp_ord,],weights = weights_opt,robust = F,init = U_coef)  
      tmp_L2 = information_score_matrix(coef(fit_samp_opt_L),weights = weights_opt,times = V[samp_ord],status = D[samp_ord],covariates = X[samp_ord,]) 
    }else
    {
      fit_samp_opt_L = coxph(Surv(R[samp_ord],V[samp_ord],D[samp_ord],type = "counting") ~ X[samp_ord,],weights = weights_opt,robust = F,init = U_coef)  
      tmp_L2 = information_score_matrix(coef(fit_samp_opt_L),weights = weights_opt,truncs = R[samp_ord],times = V[samp_ord],status = D[samp_ord],covariates = X[samp_ord,]) 
    }
    L_coef = coef(fit_samp_opt_L)
    names(L_coef) = colnames(X)
    ind_rm = (1:n_events)+q
    order_rm = tmp_L2$ord[!(tmp_L2$ord %in% ind_rm)]
    Score_L = tmp_L2$residual / tmp_L1$samp_prob[samp_ind_cens][order_rm]
    phi_mat_L = cov(Score_L)
    I_inv_L = solve(tmp_L2$hess)
    var_opt_L = I_inv_L + I_inv_L %*% phi_mat_L %*% I_inv_L/q
    ret = list("coef" = L_coef, "var" = var_opt_L)
    return(ret)
  }
  if(method == "A")
  {
    if(is.null(R))
    {
      tmp_A1 = information_score_matrix(U_coef,times = V,status = D,covariates = X, samp_prob = "A")  
    }else
    {
      tmp_A1 = information_score_matrix(U_coef,times = V,truncs = R,status = D,covariates = X, samp_prob = "A")  
    }
    D_A = D[tmp_A1$ord]
    cens_ind_ord = which(!D_A) 
    
    # random sampling with A-optimal probabilities from censored
    samp_ind_cens = sample(1:n_cens,q,replace = T,prob = tmp_A1$samp_prob)  #sampling q censored observations
    samp_ind_opt = c(cens_ind_ord[samp_ind_cens],which(D_A)) #joining the sampled censored with the failure times
    cens_weights_opt = (1/(tmp_A1$samp_prob * q))[samp_ind_cens]
    weights_opt = c(cens_weights_opt,rep(1,n_events))
    
    samp_ord = tmp_A1$ord[samp_ind_opt]
    if(is.null(R))
    {
      fit_samp_opt_A = coxph(Surv(time=V[samp_ord],event=D[samp_ord]) ~ X[samp_ord,],weights = weights_opt,robust = F,init = U_coef)  
      tmp_A2 = information_score_matrix(coef(fit_samp_opt_A),weights = weights_opt,times = V[samp_ord],status = D[samp_ord],covariates = X[samp_ord,]) 
    }else
    {
      fit_samp_opt_A = coxph(Surv(R[samp_ord],V[samp_ord],D[samp_ord],type = "counting") ~ X[samp_ord,],weights = weights_opt,robust = F,init = U_coef)  
      tmp_A2 = information_score_matrix(coef(fit_samp_opt_A),weights = weights_opt,truncs = R[samp_ord],times = V[samp_ord],status = D[samp_ord],covariates = X[samp_ord,]) 
    }
    A_coef = coef(fit_samp_opt_A)
    names(A_coef) = colnames(X)
    ind_rm = (1:n_events)+q
    order_rm = tmp_A2$ord[!(tmp_A2$ord %in% ind_rm)]
    Score_A = tmp_A2$residual / tmp_A1$samp_prob[samp_ind_cens][order_rm]
    phi_mat_A = cov(Score_A)
    I_inv_A = solve(tmp_A2$hess)
    var_opt_A = I_inv_A + I_inv_A %*% phi_mat_A %*% I_inv_A/q
    ret = list("coef" = A_coef, "var" = var_opt_A)
    return(ret)
  }
}

