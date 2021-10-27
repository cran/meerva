#======================================================================================================
#======= meerva.sim.nrm_210503.R                                                               ========
#======================================================================================================
#==========  normal data generating function  =======================================================
#====================================================================================================
#' Simulate Linear Regression Data with Measurement Errors in Outcome and Predictors 
#' 
#' @description 
#' The meerva package is designed to analyze data with measurement error when there is a 
#' validation subsample.  The merva.sim.nrm function generates a simulated data set for the 
#' linear regression 
#' setting demonstrating the data form expected for input to the meervad.fit function.
#' This simulation function first generates 4 reference predictors based upon 
#' a multivariate normal distribution, with variance-covariance specified by the 
#' user.  The first two predictors are dichotomized to have probabilites specified 
#' by the user.  This results in two class and two quantitative reference predictor 
#' variables.   
#' The response variable may have a surrogate with differential measurement error.  
#' There is one yes/no surrogate predictor variable involving error in place of 
#' one of the yes/no reference predictors, and one quantitative surrogate predictor 
#' variable involving error in place of one of the quantitative reference predictors.  
#' The simulated data are not necessarily realistic, but their analysis shows how 
#' even with rather strong measurement error the method yields reasonable solutions.   
#' The method is able to handle different types of measurement error without the user 
#' having to specify any relationship between the reference variables measured without 
#' error and the surrogate variables measured with error. 
#'    
#' @param n The full dataset size. 
#' @param m The validation subsample size (m < n). 
#' @param beta A vector of length 5 for the true regression parameter for the linear 
#'   regression model with 5 predictors including the intercept.  
#' @param alpha1 a vector of length four determining the measurement error for the outcome.
#'   if x1==1 then the error has mean alpha1[1] and variance alpha1[2].
#'   if x1==0 then the error has mean alpha1[3] and variance alpha1[4].
#' @param alpha2 A vector describing the correct classification probabilities for the surrogate for x1.
#'   If the outcome variable has positive error, then alpha2[1] and alpha2[2] 
#'   are the probabilities of correct classification when x1 is 1 or 0.  
#'   If the outcome variable has negative error, then alpha2[3] and alpha2[4] 
#'   are the probabilities of correct classification when x1 is 1 or 0.  
#' @param bx3s1 A vector of length 5 determining the relation between the reference variable x3 
#'   and the mean and SD of the surrogate x3s1.  
#'   Roughly, bx3s1[1] determines a minimal measurement error SD,   
#'   conditional on x3 bx3s1[2] determines a rate of increase in SD for values of x3 greater than bx3s1[3], 
#'   bx3s1[4] is a value above which the relation between x3 and the mean of x3s is determined by the power bx3s1[5].
#'   The mean values for x3s1 are rescaled to have mean 0 and variance 1. 
#' @param bx3s2 A vector of length 3 determining scale in x3s and potentially x3s2, 
#'   a second surrogate for xs.  
#'   Roughly, bx3s2[1] takes the previously determined mean for x3s1 
#'   using bx3s1 and multiples by bx3s2[1].
#'   Conditional on x3, x3s2 has mean  bx3s2[2] * x3 and variance bx3s2[3]. 
#' @param sd The sd of outcome y
#' @param fewer When set to 1 x3s1 and x4 will be collapsed to one 
#'   variable in the surrogate set.  This demonstrates how the method
#'   works when there are fewer surrogate variables than reference 
#'   variables.  If bx3s2 is specified such that there are 
#'   duplicate surrogate variables for the reference variable x3 
#'   then the number of surrogate predictors will not be reduced.
#' @param bx12 Bernoulli probabilities for reference variables x1 and x2.  
#'   A vector of length 2, default is c(0.25, 0.15). If mncor (see below)
#'   is positive the correlations between these Bernoulli and continuous
#'   predictors remains positive.  
#' @param mncor Correlation of the columns in the x matrix before 
#'   x1 and x2 are dichotomized to Bernoulli random variables. 
#'   Default is 0.
#' @param sigma A 4x4 varaince-covarniance matrix for the 
#'   multivarite normal dsitribution used to derive the 4 
#'   reference predictor variables.  
#'   
#' @return meerva.sim.nrm returns a list containing vectors and matrices 
#'   which can be used as example input to the meerva.fit function.  
#'      
#' @export
#' 
#' @importFrom stats rbinom rnorm var 
#'
#' @seealso
#'   \code{\link{meerva.sim.block}} , \code{\link{meerva.sim.brn}} , \code{\link{meerva.sim.cox}} , \code{\link{meerva.fit}} 
#'  
#' @examples 
#' # Simulate linear regression data with mesurement errors
#' simd = meerva.sim.nrm(beta=c(-0.5, 0.5, 0.2, 1, 0.5), 
#'     alpha1=c(-0.05, 0.1, 0.05, 0.1), alpha2=c(0.95, 0.91, 0.9, 0.9),
#'     bx3s1=c(0.05, 0, 0, NA, NA), bx3s2 = c(1.1, 0.9, 0.05) )
#'                
#' simd = meerva.sim.nrm(beta=c(-0.5, 0.5, 0.2, 1, 0.5), 
#'     alpha1=c(-0.05, 0.1, 0.05, 0.1), alpha2=c(0.95, 0.91, 0.9, 0.9),
#'     bx3s1=c(0.05, 0, 0, NA, NA), bx3s2 = c(1.1, NA, NA), fewer=1 )
#' # Copy the data vectors and matrices to input to meerva.fit
#' x_val  = simd$x_val
#' y_val  = simd$y_val
#' xs_val = simd$xs_val
#' ys_val = simd$ys_val
#' xs_non = simd$xs_non
#' ys_non = simd$ys_non
#' 
#' # Analyze the data and display results  
#' nrmout = meerva.fit(x_val, y_val, xs_val, ys_val, xs_non, ys_non )
#' summary(nrmout)
#'
meerva.sim.nrm = function(n=4000, m=400, beta=c(-0.50, 0.5, 0.2, 1, 0.5),
                          alpha1 = c(0, 0, 0, 0 ), alpha2 = c(1, 1, 1, 1),
                          bx3s1  = c( NA, NA,  NA, NA, NA), bx3s2 = c(NA, NA, NA),
                          sd=1, fewer=0 , bx12=c(0.25, 0.15), mncor=0, sigma=NULL ) {

  if (!is.null(sigma)) {
    if (matrixcalc::is.positive.definite(sigma) == 0) {
      print("Sigma is not positive definite and will be ignored.")
      sgima = NULL 
    }
  }
  if (is.null(sigma)) {
    if (is.na(mncor)) { mncor = 0 }
    mncor = min(1,max(mncor,0))
    sigma = matrix(mncor, nrow=4, ncol=4)
    diag(sigma) = c(1,1,1,1)
  } 
  x = mvtnorm::rmvnorm(n, mean=c(0,0,0,0), sigma=sigma) 
  x[,1] = 1*(x[,1] >= qnorm(1-bx12[1]) )
  x[,2] = 1*(x[,2] >= qnorm(1-bx12[2]) )
  x1 = x[,1] ; x2 = x[,2] ; x3 = x[,3] ; x4 = x[,4]
  xb = cbind(1,x) %*% beta
  
  y = as.vector( xb + sd * rnorm(n) )

  if ( is.na(alpha1[1])    ) { alpha1[1] = 0 }
  if ( is.na(alpha1[2])    ) { alpha1[2] = 0 }
  if ( alpha1[2] < 0       ) { alpha1[2] = 0 }
  if ( is.na(alpha1[3])    ) { alpha1[3] = -alpha1[1] }
  if ( is.na(alpha1[4])    ) { alpha1[4] =  alpha1[2] }
  ye_m  = alpha1[1] * (x1==1)  + alpha1[3] * (x1==0)
  ye_sd = alpha1[2] * (x1==1)  + alpha1[4] * (x1==0)
  ys    = y + ye_m + ye_sd * rnorm(n)

  pr_x1s = alpha2[1] * ((x1==1) & (y>=xb)) + (1-alpha2[2] ) * ((x1==0) & (y>=xb)) +
           alpha2[3] * ((x1==1) & (y <xb)) + (1-alpha2[4] ) * ((x1==0) & (y <xb)) ;
  x1s = rbinom(n,1,pr_x1s)

  x3s1 = x3
  x3s2 = rep(0,length=n) ;
  if ( is.na(bx3s1[1] ) == 0) {
    if ( is.na(bx3s1[4]) == 0 ) { x3 = (x3 <=  bx3s1[4]) * bx3s1[4]  + (x3 >  bx3s1[4]) * x3 }
    x3s1_s = bx3s1[1] + bx3s1[2] * (x3 >= bx3s1[3]) * (x3 - bx3s1[3]) ;  ## plot(x3, x3s1_s)
    if ( is.na(bx3s1[4]) == 1 ) { x3s1_m =  x3 }
    if ( is.na(bx3s1[4]) == 0 ) { x3s1_m = (x3 - bx3s1[4])^bx3s1[5] ;              ## take power
    x3s1_m = (x3s1_m - mean(x3s1_m))/sqrt(var(x3s1_m)) ; }               ## standardize # plot(x3s1, x3)

    if ( is.na(bx3s2[1]) == 0 ) { x3s1_m = bx3s2[1] * x3s1_m }
    x3s1   = x3s1_m + x3s1_s * rnorm(n) ;                                ## plot(x3, x3s1)
    if ( (is.na(bx3s2[2]) == 0) & (is.na(bx3s2[3]) == 0) ) { x3s2 = bx3s2[2] * x3 + bx3s2[3] * rnorm(n) }
  }

  temp = sample(n, m) ; 
  id_val = temp[order(temp)]  ;  # id_val[1:10] ; 
  id_non = c(1:n)[-id_val]    ;  # id_non[1:20]
  
  x = cbind(x1, x2, x3, x4) ; 
  simd = cbind(x1, x2, x3, x4, y, x1s, x3s1, x3s2, ys)  
  
  temp = 1 ; 
  if ( (is.na(bx3s2[2]) == 0) & (is.na(bx3s2[3]) == 0) ) { x3s = simd[,7:8] ; temp = 0 ; 
  }  else                                                { x3s = simd[,7]   } 
  
  if ((fewer==1) & (temp==1)) {    
    xs = cbind(x1s=simd[,6], x2=simd[,2], x3s = (simd[,7] + simd[,4])/2 )  
    }
  else { xs = cbind(x1s=simd[,6], x2=simd[,2], x3s, x4=simd[,4]) }

  x_val = x [id_val,]           
  y_val = y [id_val ]            
  xs_val= xs[id_val,]          
  ys_val= ys[id_val ]           
  xs_non= xs[-id_val,]     
  ys_non= ys[-id_val ]    
  
  rlist = list(simd=simd, id_val=id_val, id_non=id_non, x_val=x_val, y_val=y_val, xs_val=xs_val, ys_val=ys_val, xs_non=xs_non, ys_non=ys_non) 
  
  return( rlist )
}

## ========================================================================================== ##
## ========================================================================================== ##
