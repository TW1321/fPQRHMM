
#' Quantile regression  using fPQR+HMM
#'
#'  This function is used to fit quantile regression models using the fPQR and HMM algorithms and is based
#'  off of code written by Civieta et al. in 2021.
#' @param X_0  predictor matrix
#' @param Y_0 response matrix
#' @param q  quantile you are trying to predict
#'
#' @return A list containing the following
#' \itemize{
#'   \item  rbind(intercept,beta): estimated regression coefficients for X_0
#'   \item  final_mod: Final fitted model for quantile regression of Y_0 on score matrix
#'   \item decomp_time: Amount of time for fPQR to decompose X
#'   \item select_time: Amount of time for QREM_vs (aka HMM) to run
#'   \item fit_time: Amount of time for QREM to run
#'   \item convert_time: Amount of time to convert coefficients for score matrix into coefficients for X_0
#'  }
#'
#' @export
#'
#' @examples X_0= MASS::mvrnorm(100, mu=rep(0,100), Sigma=diag(100))
#'  beta=rep(1,100)
#'  Y_0= X_0%*%beta+rnorm(100)
#'  fPQR_mod= fPQR_fit(X_0, Y_0, 0.5)
fPQR_fit <- function(X_0,Y_0,q) {

  n_iter=0
  X_init=scale(X_0,scale=FALSE)
  Y_init=scale(Y_0, scale=FALSE)
  n=nrow(X_0)
  p=ncol(X_0)
  r=ncol(Y_0)
  n_components=p

  x_weights=matrix(0, nrow=p, ncol=n_components)
  y_weights=matrix(0, nrow=r, ncol=n_components)
  x_scores=matrix(0,  nrow=n, ncol=n_components)
  y_scores=matrix(0,  nrow=n, ncol=n_components)
  x_loadings=matrix(0, nrow=p, ncol=n_components)
  y_loadings=matrix(0, nrow=r, ncol=n_components)

  start=proc.time()[3]
  while(sum(X_init^2)>10^-10 & sum(Y_init^2)>10^-10) {
    qcov_xy=matrix(0,nrow=p, ncol=r)
    for (i in 1:r){
      qcov_xy[,i]=qcov_li(X_init,Y_init[,i],q)
    }
    svd_cov=svd(qcov_xy)
    x_weight=svd_cov$u[,1]
    y_weight=svd_cov$v[1,]
    flip_weights=signflip_svd(x_weight, y_weight)
    x_weight=flip_weights[[1]]
    y_weight=flip_weights[[2]]

    x_score=X_init%*%x_weight
    y_score=Y_init%*%y_weight

    x_loading=(t(X_init)%*%x_score)/(t(x_score)%*%x_score)[1]
    y_loading=(t(Y_init)%*%x_score)/(t(x_score)%*%x_score)[1]

    X_init=X_init- x_score%*%t(x_loading)
    Y_init=Y_init- x_score%*%t(y_loading)

    x_weights[,n_iter+1]=x_weight
    y_weights[,n_iter+1]=y_weight

    x_scores[,n_iter+1]=x_score
    y_scores[,n_iter+1]=y_score

    x_loadings[,n_iter+1]=x_loading
    y_loadings[,n_iter+1]=y_loading

    n_iter=n_iter+1

  }
  decomp_time=proc.time()[3]-start

  x_weights=x_weights[,1:n_iter]
  y_weights=y_weights[,1:n_iter]

  x_scores=x_scores[,1:n_iter]
  y_scores=y_scores[,1:n_iter]

  x_loadings=x_loadings[,1:n_iter]
  y_loadings=y_loadings[,1:n_iter]

  x_rot=x_weights%*%solve(t(x_loadings)%*%x_weights)
  y_rot=y_weights%*%solve(t(y_loadings)%*%y_weights)

  intercept=matrix(0, nrow=1, ncol=r)
  beta=matrix(0, nrow=p, ncol=r)

  select_time=0
  fit_time=0
  convert_time=0
  for (i in 1:r){
    Y_col=Y_0[,i]
    dat_train=data.frame(Y_col,x_scores)
    start=proc.time()[3]
    select=QREM::QREM_vs(dat_train,1,2:ncol(dat_train),qn=q)
    select_time=select_time+proc.time()[3]-start
    var_selec=select$fittedSEMMS$gam.out$nn

    if(length(var_selec)==0){
      start=proc.time()[3]
      final_mod=QREM::QREM(stats::lm, Y_col~1, dat_train, qn=q)
      fit_time=fit_time+proc.time()[3]-start
      intercept[i]=final_mod$coef$beta[1]
    }else{
      dfsemms=dat_train[,c(1,1+var_selec)]
      start=proc.time()[3]
      final_mod=QREM::QREM(stats::lm, Y_col~., dfsemms, qn=q)
      fit_time=fit_time+proc.time()[3]-start
      intercept[i]=final_mod$coef$beta[1]
      start=proc.time()[3]
      beta[,i]=as.matrix(x_rot[,var_selec])%*%as.matrix(final_mod$coef$beta[-1])
      convert_time=convert_time+proc.time()[3]-start
    }
  }

  fPQR_fit=list(rbind(intercept,beta),final_mod,decomp_time,select_time, fit_time, convert_time)
}


#' Quantile Covariance
#'
#' This function calculates quantile covariance between two variables
#' using the notion of quantile covariance developed by Li et al. in 2015
#' @param X covariate matrix
#' @param Y response matrix
#' @param q  quantile you are predicting
#'
#' @return qcov_li: quantile covariance of X and Y
#' @export
#'
#' @examples X_0=MASS::mvrnorm(100, mu=rep(0,100), Sigma=diag(100))
#'  beta=rep(1,100)
#'  Y_0= X_0%*%beta+rnorm(100)
#'  qcov= qcov_li(X_0, Y_0, 0.5)
qcov_li<-function(X,Y, q){

  n=nrow(X)
  psi=q-((Y-stats::quantile(Y,q))<0)
  X_center=scale(X, scale=FALSE)
  qcov_li=(1/n)*t((t(psi)%*%X_center))
}


#' Sign flip (based on fPQR code from Civieta et al.)
#'
#' This function flips the signs of u and v to ensure consistency among different solutions for singular value decomposition, and is based off of code written by Civieta et al. in 2021.
#'
#' @param u first column of U from singular value decomposition UDV^t
#' @param v first row of V from singular value decomposition UDV^t
#'
#' @return List of two items
#' \itemize{
#'   \item  u and v if largest value of u is positive.
#'   \item  -u and -v if largest value of u is negative .
#'  }
#' @export
#'
#' @examples X=MASS::mvrnorm(100, mu=rep(0,100), Sigma=diag(100))
#' svd_X=svd(X)
#' col_weight=svd_X$u[,1]
#' row_weight=svd_X$v[1,]
#' flip_weights=signflip_svd(col_weight, row_weight)
signflip_svd=function(u,v){

  biggest_value_index=which(abs(u)==max(abs(u)))
  largest_sign=sign(u[biggest_value_index])
  signflip_svd=list(u*largest_sign, v*largest_sign)
}
#' Check loss
#'
#' This function calculates check loss error for a fitted quantile regression model.
#' @param beta_pred coefficients for model
#' @param X_test test set for predictor
#' @param Y_test test set for response
#' @param q quantile you are predicting
#'
#' @return check_lossTE= check loss error of model given by beta_pred on test set given by X_test and Y_test
#' @export
#'
#' @examples  X_0=MASS::mvrnorm(100, mu=rep(0,100), Sigma=diag(100))
#'  beta=rep(1,100)
#'  Y_0= X_0%*%beta+rnorm(100)
#'  fPQR_mod= fPQR_fit(X_0, Y_0, 0.5)
#'  est_coef=fPQR_mod[[1]]
#'  check_losserror=check_lossTE(est_coef, X_0, Y_0,0.5)
check_lossTE=function(beta_pred, X_test, Y_test,q){
  test_pred=beta_pred[1]+as.matrix(X_test)%*%beta_pred[-1]
  u=Y_test-test_pred
  check_lossTE=mean(u*(q-(u<0)))
}
