#install.packages(c("fda", "MASS"))
library(fda)
library(MASS)

# response            : a matrix containing response curves (center)
# response_l          : a matrix containing lower limits of the response curves
# response_u          : a matrix containing upper limits of the response curves
# response_r          : a matrix containing half-ranges of the response curves
# predictor           : a list whose elements consist of predictors (center)
# predictor_l         : a list whose elements consist of lower limits of the predictors
# predictor_u         : a list whose elements consist of upper limits of the predictors
# predictor_r         : a list whose elements consist of half-ranges of the predictors
# predictor_test      : a list whose elements consist of predictors (center) for the test sample
# predictor_test_l    : a list whose elements consist of lower limits of the predictors for the test sample
# predictor_test_u    : a list whose elements consist of lower upper of the predictors for the test sample
# predictor_test_r    : a list whose elements consist of half-ranges of the predictors for the test sample
# nbf_response        : number of basis function for approximating response variable
# nbf_vec_predictors  : a vector of number of basis function corresponding to predictors
# B                   : number of Monte Carlo simulations to be used in the MCM method

iv_fof = function(response, response_l, response_u, response_r, predictor, predictor_l, predictor_u, predictor_r,
                 predictor_test, predictor_test_l, predictor_test_u, predictor_test_r,
                 nbf_response, nbf_vec_predictors, B)
{
  # Number of predictors
  np = length(predictor_test)
  
  # Discrete time points
  dtp = seq(0, 1, length=ncol(response))
  
  # B-spline for response
  B_spline_basis_y = create.bspline.basis(c(0,1), nbasis = nbf_response)
  B_spline_basis_fun_y = eval.basis(dtp, B_spline_basis_y)
  
  # B-spline for predictors
  B_spline_basis_x = vector("list",)
  B_spline_basis_funs_x = vector("list",)
  
  for(i in 1:np){
    B_spline_basis_x[[i]] = create.bspline.basis(c(0,1), nbasis = nbf_vec_predictors[i])
    B_spline_basis_funs_x[[i]] = eval.basis(dtp, B_spline_basis_x[[i]])
  }
  
  #Inner products
  Inner_prod = vector("list",)
  
  for(i in 1:np)
    Inner_prod[[i]] = inprod(B_spline_basis_x[[i]], B_spline_basis_x[[i]])
  
  # Weight argument
  w_arg = matrix(dtp, nrow = nrow(response), ncol = ncol(response), byrow=T)
  
  # Weight matrices of the response and predictors
  W_y = t(smooth.basis(argvals=t(w_arg), y=t(response), fdParobj=B_spline_basis_y)$fd$coefs)
  W_yl = t(smooth.basis(argvals=t(w_arg), y=t(response_l), fdParobj=B_spline_basis_y)$fd$coefs)
  W_yu= t(smooth.basis(argvals=t(w_arg), y=t(response_u), fdParobj=B_spline_basis_y)$fd$coefs)
  W_yr= t(smooth.basis(argvals=t(w_arg), y=t(response_r), fdParobj=B_spline_basis_y)$fd$coefs)
  
  W_yl_out = W_yl
  W_yu_out = W_yu
  
  W_x = vector("list",)
  W_xl = vector("list",)
  W_xu = vector("list",)
  W_xr = vector("list",)
  
  W_x_test = vector("list",)
  W_xl_test = vector("list",)
  W_xu_test = vector("list",)
  W_xr_test = vector("list",)
  
  for(i in 1:np){
    W_x[[i]] = t(smooth.basis(argvals=t(w_arg), y=t(predictor[[i]]), fdParobj=B_spline_basis_x[[i]])$fd$coefs)
    W_xl[[i]] = t(smooth.basis(argvals=t(w_arg), y=t(predictor_l[[i]]), fdParobj=B_spline_basis_x[[i]])$fd$coefs)
    W_xu[[i]] = t(smooth.basis(argvals=t(w_arg), y=t(predictor_u[[i]]), fdParobj=B_spline_basis_x[[i]])$fd$coefs)
    W_xr[[i]] = t(smooth.basis(argvals=t(w_arg), y=t(predictor_r[[i]]), fdParobj=B_spline_basis_x[[i]])$fd$coefs)
    
    W_x_test[[i]] = t(smooth.basis(argvals=t(w_arg), y=t(predictor_test[[i]]), fdParobj=B_spline_basis_x[[i]])$fd$coefs)
    W_xl_test[[i]] = t(smooth.basis(argvals=t(w_arg), y=t(predictor_test_l[[i]]), fdParobj=B_spline_basis_x[[i]])$fd$coefs)
    W_xu_test[[i]] = t(smooth.basis(argvals=t(w_arg), y=t(predictor_test_u[[i]]), fdParobj=B_spline_basis_x[[i]])$fd$coefs)
    W_xr_test[[i]] = t(smooth.basis(argvals=t(w_arg), y=t(predictor_test_r[[i]]), fdParobj=B_spline_basis_x[[i]])$fd$coefs)
  }
  
  # Mean of variables
  mean_y = apply(W_y, 2, mean)
  mean_yl = apply(W_yl, 2, mean)
  mean_yu = apply(W_yu, 2, mean)
  mean_yr = apply(W_yr, 2, mean)
  
  mean_x = vector("list",)
  mean_xl = vector("list",)
  mean_xu = vector("list",)
  mean_xr = vector("list",)
  
  mean_x_test = vector("list",)
  mean_xl_test = vector("list",)
  mean_xu_test = vector("list",)
  mean_xr_test = vector("list",)
  
  for(i in 1:np){
    mean_x[[i]] = apply(W_x[[i]], 2, mean)
    mean_xl[[i]] = apply(W_xl[[i]], 2, mean)
    mean_xu[[i]] = apply(W_xu[[i]], 2, mean)
    mean_xr[[i]] = apply(W_xr[[i]], 2, mean)
    
    mean_x_test[[i]] = apply(W_x_test[[i]], 2, mean)
    mean_xl_test[[i]] = apply(W_xl_test[[i]], 2, mean)
    mean_xu_test[[i]] = apply(W_xu_test[[i]], 2, mean)
    mean_xr_test[[i]] = apply(W_xr_test[[i]], 2, mean)
  }
  
  # Centered variables
  for(i in 1:nrow(response)){
  W_y[i,] = W_y[i,] - mean_y
  W_yl[i,] = W_yl[i,] - mean_yl
  W_yu[i,] = W_yu[i,] - mean_yu
  W_yr[i,] = W_yr[i,] - mean_yr
  }
  
  for(i in 1:np){
    for(j in 1:nrow(response)){
    W_x[[i]][j,] = W_x[[i]][j,] - mean_x[[i]]
    W_xl[[i]][j,] = W_xl[[i]][j,] - mean_xl[[i]]
    W_xu[[i]][j,] = W_xu[[i]][j,] - mean_xu[[i]]
    W_xr[[i]][j,] = W_xr[[i]][j,] - mean_xr[[i]]
    
    W_x_test[[i]][j,] = W_x_test[[i]][j,] - mean_x_test[[i]]
    W_xl_test[[i]][j,] = W_xl_test[[i]][j,] - mean_xl_test[[i]]
    W_xu_test[[i]][j,] = W_xu_test[[i]][j,] - mean_xu_test[[i]]
    W_xr_test[[i]][j,] = W_xr_test[[i]][j,] - mean_xr_test[[i]]
    }
  }

  
  # Matrices for the regressions
  Reg_mat = vector("list",)
  Reg_mat_l = vector("list",)
  Reg_mat_u = vector("list",)
  Reg_mat_r = vector("list",)
  
  Reg_mat_test = vector("list",)
  Reg_mat_test_l = vector("list",)
  Reg_mat_test_u = vector("list",)
  Reg_mat_test_r = vector("list",)
  
  for(i in 1:np){
    Reg_mat[[i]] = W_x[[i]] %*% Inner_prod[[i]]
    Reg_mat_l[[i]] = W_xl[[i]] %*% Inner_prod[[i]]
    Reg_mat_u[[i]] = W_xu[[i]] %*% Inner_prod[[i]]
    Reg_mat_r[[i]] = W_xr[[i]] %*% Inner_prod[[i]]
    
    Reg_mat_test[[i]] = W_x_test[[i]] %*% Inner_prod[[i]]
    Reg_mat_test_l[[i]] = W_xl_test[[i]] %*% Inner_prod[[i]]
    Reg_mat_test_u[[i]] = W_xu_test[[i]] %*% Inner_prod[[i]]
    Reg_mat_test_r[[i]] = W_xr_test[[i]] %*% Inner_prod[[i]]
  }
  
  Reg_mat = cbind(1, do.call(cbind, Reg_mat))
  Reg_mat_l = cbind(1, do.call(cbind, Reg_mat_l))
  Reg_mat_u = cbind(1, do.call(cbind, Reg_mat_u))
  Reg_mat_r = cbind(1, do.call(cbind, Reg_mat_r))
  Reg_mat_bcrm = cbind(1, cbind(Reg_mat, Reg_mat_r))
  
  Reg_mat_test = cbind(1, do.call(cbind, Reg_mat_test))
  Reg_mat_test_l = cbind(1, do.call(cbind, Reg_mat_test_l))
  Reg_mat_test_u = cbind(1, do.call(cbind, Reg_mat_test_u))
  Reg_mat_test_r = cbind(1, do.call(cbind, Reg_mat_test_r))
  Reg_mat_test_bcrm = cbind(1, cbind(Reg_mat_test, Reg_mat_test_r))
  
  ###################################### MCM Start ######################################
  Bhat_mcm = array(dim=c((sum(nbf_vec_predictors)+1), nbf_response, B))
  for(boot in 1:B){
    ### Generating MCM data
    Ymcm = matrix(0, nrow=nrow(response), ncol=ncol(response))
    for(i in 1:nrow(response)){
      for(j in 1:ncol(response)){
        Ymcm[i,j] = runif(1, min=response_l[i,j], max=response_u[i,j])
      }
    }
    
    Xmcm = vector("list",)
    for(k in 1:np){
      Xmcm1 = matrix(0, nrow=nrow(predictor[[k]]), ncol=ncol(predictor[[k]]))
      for(i in 1:nrow(predictor[[k]])){
        for(j in 1:ncol(predictor[[k]])){
          Xmcm1[i,j] <- runif(1, min=predictor_l[[k]][i,j], max=predictor_u[[k]][i,j])
        }
      }
      Xmcm[[k]] = Xmcm1
    }
    
    W_y_mcm = t(smooth.basis(argvals=t(w_arg), y=t(Ymcm), fdParobj=B_spline_basis_y)$fd$coefs)
    
    W_x_mcm = vector("list",)
    for(i in 1:np)
      W_x_mcm[[i]] = t(smooth.basis(argvals=t(w_arg), y=t(Xmcm[[i]]), fdParobj=B_spline_basis_x[[i]])$fd$coefs)
    
    mean_y_mcm = apply(W_y_mcm, 2, mean)
    
    mean_x_mcm = vector("list",)
    
    for(i in 1:np)
      mean_x_mcm[[i]] = apply(W_x_mcm[[i]], 2, mean)
    
    for(i in 1:nrow(response))
      W_y_mcm[i,] = W_y_mcm[i,] - mean_y_mcm
    
    for(i in 1:np){
      for(j in 1:nrow(response))
        W_x_mcm[[i]][j,] = W_x_mcm[[i]][j,] - mean_x_mcm[[i]]

    }
    
    Reg_mat_mcm = vector("list",)
    for(i in 1:np)
      Reg_mat_mcm[[i]] = W_x_mcm[[i]] %*% Inner_prod[[i]]
    
    Reg_mat_mcm = cbind(1, do.call(cbind, Reg_mat_mcm))
    
    Bhat_mcm[,,boot] = ginv(t(Reg_mat_mcm)%*%Reg_mat_mcm) %*% t(Reg_mat_mcm)%*%W_y_mcm
  }
  
  Bhat_MCM1 = Bhat_mcm[,,1]
  for(i in 2:B){
    Bhat_MCM1 = Bhat_MCM1 + Bhat_mcm[,,i]
  }
  Bhat_MCM = Bhat_MCM1/B
  ###################################### MCM End ###################################### 
  
  # Model estimation
  coeff_c = ginv(t(Reg_mat)%*%Reg_mat) %*% t(Reg_mat)%*%W_y
  coeff_l = ginv(t(Reg_mat_l)%*%Reg_mat_l) %*% t(Reg_mat)%*%W_yl
  coeff_u = ginv(t(Reg_mat_u)%*%Reg_mat_u) %*% t(Reg_mat_u)%*%W_yu
  coeff_r = ginv(t(Reg_mat_r)%*%Reg_mat_r) %*% t(Reg_mat_r)%*%W_yr
  coeff_bcrm_c = ginv(t(Reg_mat_bcrm)%*%Reg_mat_bcrm) %*% t(Reg_mat_bcrm)%*%W_y
  coeff_bcrm_r = ginv(t(Reg_mat_bcrm)%*%Reg_mat_bcrm) %*% t(Reg_mat_bcrm)%*%W_yr
  
  pred_l1 = Reg_mat_test_l %*% coeff_l
  pred_u1 = Reg_mat_test_u %*% coeff_u
  
  pred_cm_l1 = Reg_mat_test_l %*% coeff_c
  pred_cm_u1 = Reg_mat_test_u %*% coeff_c
  
  pred_bcrm_c = Reg_mat_test_bcrm %*% coeff_bcrm_c
  pred_bcrm_r = Reg_mat_test_bcrm %*% coeff_bcrm_r
  
  pred_mcm_l1 = Reg_mat_test_l %*% Bhat_MCM
  pred_mcm_u1 = Reg_mat_test_u %*% Bhat_MCM
  
  for(i in 1:nrow(predictor_test[[1]])){
    pred_l1[i,] = pred_l1[i,] + mean_yl
    pred_u1[i,] = pred_u1[i,] + mean_yu
    
    pred_cm_l1[i,] = pred_cm_l1[i,] + mean_yl
    pred_cm_u1[i,] = pred_cm_u1[i,] + mean_yu
    
    pred_bcrm_c[i,] = pred_bcrm_c[i,] + mean_y
    pred_bcrm_r[i,] = pred_bcrm_r[i,] + mean_yr
    
    pred_mcm_l1[i,] = pred_mcm_l1[i,] + mean_yl
    pred_mcm_u1[i,] = pred_mcm_u1[i,] + mean_yu
  }
  
  pred_l = pred_l1 %*% t(B_spline_basis_fun_y)
  pred_u = pred_u1 %*% t(B_spline_basis_fun_y)
  
  pred_cm_l = pred_cm_l1 %*% t(B_spline_basis_fun_y)
  pred_cm_u = pred_cm_u1 %*% t(B_spline_basis_fun_y)
  
  pred_bcrm_l = pred_bcrm_c %*% t(B_spline_basis_fun_y) -
    pred_bcrm_r %*% t(B_spline_basis_fun_y)
  pred_bcrm_u = pred_bcrm_c %*% t(B_spline_basis_fun_y) +
    pred_bcrm_r %*% t(B_spline_basis_fun_y)
  
  pred_mcm_l = pred_mcm_l1 %*% t(B_spline_basis_fun_y)
  pred_mcm_u = pred_mcm_u1 %*% t(B_spline_basis_fun_y)
  
  return(list("flm_l" = pred_l, "flm_u" = pred_u, "cm_l" = pred_cm_l, "cm_u" = pred_cm_u,
              "bcrm_l" = pred_bcrm_l, "bcrm_u" = pred_bcrm_u, "mcm_l" = pred_mcm_l,
              "mcm_u" = pred_mcm_u, "Bhat_mcm" = Bhat_mcm, "mean_yl" = mean_yl,
              "mean_yu" = mean_yu, "B_spline_basis_fun_y" = B_spline_basis_fun_y, "W_yl" = W_yl_out,
              "W_yu" = W_yu_out, "Reg_mat_l" = Reg_mat_l, "Reg_mat_u" = Reg_mat_u, "Bhat_MCM" = Bhat_MCM))
  
}
