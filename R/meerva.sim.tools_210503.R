#======================================================================================================
#======= meerva.sim.tools_210503.R                                                               ========
#======================================================================================================
#' Simulation of meerva used to Analyze Data with Measurement Error
#' 
#' @description
#'  The meerva package is designed to analyze data with measurement error when 
#'  there is a validation subsample randomly selected from the full sample.  
#'  The method assumes surrogate variables measured with error are available 
#'  for the full sample, and reference variables measured with little or no 
#'  error are available for this randomly chosen subsample of the full sample.  
#'  Measurement errors may be differential or non differential, in any or all 
#'  predictors (simultaneously) as well as outcome.   
#' 
#'  The meerva.sim.block lets the user specify a model with measurement error, 
#'  and then simulate and analyze many datasets to 
#'  examine the model fits and judge how the method works.
#'  Data sets are generated according to 3 functions for simulating
#'  Cox PH, linear and logistic regression models.  These functions generate 
#'  data sets with 4 reference predictor variables and from 3 to 5 surrogate 
#'  predictor variables.  The user can 
#'  consider, program and simulate data sets of greater complexity
#'  but these examples provided with the package should serve as a 
#'  reasonable introduction to the robustness of the method.  
#' 
#' @param simfam The family for the underlying regression model 
#' to be simulated, amongst "binomial", "gaussian" and "Cox".
#' @param nsims Number of datasets to be simulated 
#' @param seed A seed for the R random number generator.  The default is 0 in which case the 
#'   program random selects and records the seed so one ca replicate simulation studies.
#' @param n The full dataset size. 
#' @param m The validation subsample size (m < n). 
#' @param beta A vector of length 5 for the true regression parameter for the linear 
#'   regression model with 5 predictors including the intercept. For the Cox model
#'   beta[0] is not estimated but determines a basal event rate.   
#' @param alpha1 A vector of length four determining the measurement error or 
#'   misclassification probabilities for the outcome surrogate ys.  
#'   Usage is slightly different 
#'   for the different  simfam  values "gaussian", "binomial" and "Cox".  See the
#'   help pages for meerva.sim.brn, meerva.sim.cox and meerva.sim.nrm 
#'   for clarification.
#' @param alpha2 A vector describing the correct classification probabilities for x1s, 
#'   the surrogate for x1.  
#'   Usage is slighlty different 
#'   for the different  simfam  values "gaussian", "binomial" and "Cox".  See the
#'   help pages for meerva.sim.brn, meerva.sim.cox and meerva.sim.nrm 
#'   for clarification.
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
#' @param bx12 Bernoulli probabilities for reference variables x1 and x2.  
#'   A vector of length 2, default is c(0.25, 0.15).  If mncor (see below)
#'   is positive the correlations between these Bernoulli and continuous
#'   predictors remains positive.   
#' @param sd In case of simfam == "gaussain" for linear regression, the sd of outcome y. 
#'   In case of simfam == "Cox" for Cox PH regression, the multiplicative error 
#'   term for ys, the surrogate for the time to event y 
#'   (ys = log(sd * a (random variable) * y).
#' @param fewer When set to 1 x3s1 and x4 will be collapsed to one 
#'   variable in the surrogate set.  This demonstrates how the method
#'   works when there are fewer surrogate variables than reference 
#'   variables.  If bx3s2 is specified such that there are 
#'   duplicate surrogate variables for the reference variable x3 
#'   then the number of surrogate predictors will not be reduced.
#' @param mncor Correlation of the columns in the x matrix before 
#'   x1 and x2 are dichotomized to Bernoulli random variables. 
#'   Default is 0.
#' @param sigma A 4x4 varaince-covarniance matrix for the 
#'   multivarite normal dsitribution used to derive the 4 
#'   reference predictor variables.
#' @param vmethod Method for robust estimation of variance covariance matrices needed 
#'   for calculation of the augmented estimates (beta aug).
#'   0 for JK or jackknife (slowest but more accurate), 
#'   1 for IJK or the infinitesimal JK using the R default dfbeta's
#'   2 for IJK using an alternate formula for the dfbeta, and 
#'   3 for all three of these methods to be used
#'   NA to let the program choose a stronger, faster method.
#' @param jksize leave out number for grouped jackknife used for non validation data 
#'   The default is 0 where the program chooses jksize so that the number of leave out 
#'   groups is about validation subsample size.
#' @param compare 1 to compare gamma_val with gamma_ful (default) or 0 with gamma_non.
#' @param diffam inidcates a cutoff if for a "guassian" family in surrogate a "binomial" 
#'   famliy is to be similated for the refernce model.  For example, the
#'   surrogate outcome could be an estimated probit (or logit) based upon
#'   a convolutional neural network. Normal data are simulated and
#'   y_val is repalced by 1*(y_val >= diffam).  Default is NA and
#'   the surrogate and reference have the same model form.  Only 
#'   for use with vmethod of 0 or 1.   
#'
#'@author Walter Kremers (kremers.walter@mayo.edu)
#'
#' @seealso
#'    \code{\link{meerva.fit}} , \code{\link{meerva.sim.brn}} , \code{\link{meerva.sim.cox}} , \code{\link{meerva.sim.nrm}} 
#'
#' @return
#' meerva.sim.block returns a list object of class meerva.sim.  
#' The list will contain summary information used to simulate the data, and
#' for each data set simulated with measurement error,
#' the augmented estimates based upon the full data set accounting for measurement error,
#' estimates based upon reference variables from the validation subsample,
#' estimates based upon the surrogate variables from the whole sample,
#' along with estimated variances for these estimates.
#' These can be inspected by the user directly or by as shown in the example.  
#' 
#' @export
#' 
#' @importFrom stats runif 
#' @importFrom utils timestamp 
#'
#' @examples
#' # Simulation study for logistic reg data with 
#' # differential misclassification in outcome 
#' # and a predictor and measurement error in 
#' # another predictor.  nsims=10 is as an 
#' # example only.  Try running nsims=100 or 
#' # 1000, but be prepared to wait a little while. 
#' sim.binomial = meerva.sim.block(simfam="binomial", 
#'     nsims=10, seed=0, n=4000, m=400, 
#'     beta = c(-0.5, 0.5, 0.2, 1, 0.5) , 
#'     alpha1 = c(0.95, 0.90, 0.90, 0.95), 
#'     alpha2 = c(0.98,0.94,0.95,0.95), 
#'     bx3s1=c(0.05, 0, 0, NA, NA) , 
#'     bx3s2 = c(NA,NA,NA) , 
#'     vmethod=2, jksize=0, compare=1) 
#'     
#' plot(sim.binomial) 
#' summary(sim.binomial, 1) 
#'
#' # Simulation study for linear reg data.   
#' # For this example there are more surrogate 
#' # predictors than reference predictors.  
#' # nsims=10 is as an example only.  Try 
#' # running nsims=100 or 1000, but be 
#' # prepared to wait a little while.   
#' sim.gaussianm = meerva.sim.block(simfam="gaussian", 
#'     nsims=10, seed=0, n=4000, m=400, 
#'     beta = c(-0.5, 0.5, 0.2, 1, 0.5) , 
#'     alpha1 = c(-0.05, 0.1, 0.05, 0.1) , 
#'     alpha2 = c(0.98,0.94,0.95,0.95) , 
#'     bx3s1=c(0.05, 0, 0, NA, NA) , 
#'     bx3s2  = c(1.1,0.9,0.05) ,  
#'     sd=1, fewer=0, 
#'     vmethod=1, jksize=0, compare=1) 
#'
#' plot(sim.gaussianm)
#' summary(sim.gaussianm)
#'  
#' # Simulation study for Cox PH data.  
#' # For this example there are fewer surrogates 
#' # than reference variables yet they provide 
#' # information to decrease the variance in the 
#' # augmented estimate.  nsims=10 is as an 
#' # example only.  Try running nsims=100 or 
#' # 1000, but be prepared to wait a little 
#' # while.   
#' sim.coxphf = meerva.sim.block(simfam="Cox", 
#'     nsims=10, seed=0, n=4000, m=400, 
#'     beta   = c(-0.5, 0.5, 0.2, 1, 0.5) , 
#'     alpha1 = c(0.95,0.90,0.90,0.95)  , 
#'     alpha2 = c(0.98,0.94,0.94,0.98) , 
#'     bx3s1  = c(0.05,0,0,NA,NA) , 
#'     bx3s2  = c(1.1, NA, NA) , 
#'     sd=0.1, fewer=1, 
#'     vmethod=1, jksize=0, compare=1 ) 
#'
#' plot(sim.coxphf)
#' summary(sim.coxphf)
#'  
meerva.sim.block = function(simfam="gaussian", nsims=100, seed=0, n=4000, m=400, 
      beta = c(-0.5, 0.5, 0.2, 1, 0.5) , 
      alpha1 = c(-0.05, 0.1, 0.05, 0.1) , 
      alpha2 = c(0.98,0.94,0.95,0.95) , 
      bx3s1 =c(0.05, 0, 0, NA, NA) , 
      bx3s2 = c(NA,NA,NA) , 
      bx12=c(0.25, 0.15) , 
      sd=1 , fewer=0, mncor=0, sigma=NULL , vmethod=NA , jksize=0 , compare=1, diffam=NA ) {

  if (seed == 0) { seed = round(runif(1)*1000000000) }
  set.seed(seed) ; 
  
  if (nsims < 2) { nsims = 2 }
  
  if ((vmethod==2) & (simfam == "Cox")) {vmethod=1 ; cat("vmethod cannot be 2 for Cox model.  vmethed set to 1\n") }
  
  if (is.na(vmethod)) {
    if (simfam=="binomial") { vmethod = 2
    } else { vmethod = 1 }
  }
  
  if ( is.na(simfam) ) { simfam = "gaussian" }  
  if (simfam %in% c("gausian" , "gause")) { simfam = "gaussian" } 
  
  if (is.na(mncor)) { mncor = 0 }
  mncor = min(1,max(mncor,0))
  
  for (nsim in 1:nsims) {
    if (nsim == 1)  { cat(" before ", nsim, " of ", nsims, " iterations ", timestamp(), "\n") }  
    
    if ( simfam == "binomial" ) { simd = meerva.sim.brn(n, m, beta, alpha1, alpha2, bx3s1, bx3s2, bx12=bx12, fewer=fewer, mncor=mncor, sigma=sigma) }  
    if ( simfam == "gaussian" ) { simd = meerva.sim.nrm(n, m, beta, alpha1, alpha2, bx3s1, bx3s2, bx12=bx12, sd=sd, fewer=fewer, mncor=mncor, sigma=sigma) } 
    if ( simfam == "Cox"      ) { simd = meerva.sim.cox(n, m, beta, alpha1, alpha2, bx3s1, bx3s2, bx12=bx12, sd=sd, fewer=fewer, mncor=mncor, sigma=sigma) }  
#    if ( simfam == "test1"    ) { simd = meerva.sim.test(n, m, beta, alpha1, alpha2, opt=1) } 
#    if ( simfam == "test2"    ) { simd = meerva.sim.test(n, m, beta, alpha1, alpha2, opt=2) }     
    
    x_val  = simd$x_val
    y_val  = simd$y_val
    xs_val = simd$xs_val
    ys_val = simd$ys_val
    xs_non = simd$xs_non
    ys_non = simd$ys_non
    if ((simfam == "gaussian") & (!is.na(diffam[1]))) { 
      ys_val = 1 * (ys_val >= diffam)
      ys_non = 1 * (ys_non >= diffam) 
      if (nsim==1) { tmp = round(mean(ys_non),3) ; cat(' mean yx_non about ', tmp, "\n") } 
    }
     
    if ( simfam == "Cox" ) {
      e_val  = simd$e_val
      es_val = simd$es_val
      es_non = simd$es_non    
    }
    
    if (vmethod %in% c(0, 3)) {
      if ( simfam == "Cox" ) {
        meerva.0 <- meerva.fit(x_val, y_val, xs_val, ys_val, xs_non, ys_non, 
                                     familyr="Cox", e_val, es_val, es_non, vmethod=0, jksize=jksize, compare=compare )  
      } else if ((simfam == "gaussian") & (!is.na(diffam[2])) & (diffam[2]==0) ) {
        meerva.0 <- meerva.fit(x_val, y_val, xs_val, ys_val, xs_non, ys_non, familyr="gaussian", familys="gaussian", 
                               vmethod=0, jksize=jksize, compare=compare) 
      } else {meerva.0 <- meerva.fit(x_val, y_val, xs_val, ys_val, xs_non, ys_non, vmethod=0, jksize=jksize, compare=compare) 
      }
    } 
    
    if (vmethod %in% c(1, 3)) {
      if (simfam == "Cox") {
        meerva.1 <- meerva.fit(x_val, y_val, xs_val, ys_val, xs_non, ys_non,  
                               familyr="Cox", e_val, es_val, es_non, vmethod=1, jksize=jksize, compare=compare )  
      } else if ((simfam == "gaussian") & (!is.na(diffam[2])) & (diffam[2]==0)) {
        meerva.1 <- meerva.fit(x_val, y_val, xs_val, ys_val, xs_non, ys_non, familyr="gaussian", familys="gaussian", 
                              jksize=jksize, compare=compare )
      } else { meerva.1 = meerva.fit(x_val, y_val, xs_val, ys_val, xs_non, ys_non, vmethod=1, jksize=jksize, compare=compare) }
    }
    
    if (vmethod %in% c(2, 3)) {
     meerva.2 = meerva.fit(x_val, y_val, xs_val, ys_val, xs_non, ys_non, vmethod=2, jksize=jksize, compare=compare) 
    }
    
    if ( (nsim %in% c(1, 10, 50)) | (((nsim %% 100) == 0) & (nsim<=1000)) | (((nsim %% 1000) == 0) & (nsim<=10000)) | (nsim == nsims) ) { 
      cat(" after ", nsim, " of ", nsims, " iterations ", timestamp(), "\n")  
    } 
  
    if        (vmethod == 2) { meerva.0 = meerva.2 
    } else if (vmethod == 1) { meerva.0 = meerva.1 
    }
    
    if (nsim==1) { 
      beta_augs   =  meerva.0$coef_beta[1,]  ;  beta_aug_vars   = meerva.0$var_beta[1,]  ;  
      beta_vals   =  meerva.0$coef_beta[2,]  ;  beta_val_vars   = meerva.0$var_beta[2,]  ;  
      gamma_bigs  =  meerva.0$coef_gamma[1,] ;  gamma_big_varjs = meerva.0$var_gamma[1,] ; 
                                                gamma_big_vars  = meerva.0$var_gamma[2,] ; 
      gamma_vals  =  meerva.0$coef_gamma[2,] ;  gamma_val_vars  = meerva.0$var_gamma[3,] ; 
      if  (vmethod==3)                      { beta_aug_1s = meerva.1$coef_beta[1,]  ;  beta_aug_1_vars = meerva.1$var_beta[1,]  } 
      if ((vmethod==3) & (simfam != "Cox")) { beta_aug_2s = meerva.2$coef_beta[1,]  ;  beta_aug_2_vars = meerva.2$var_beta[1,]  } 
    }
    
    if (nsim!=1) {
      beta_augs  = rbind( beta_augs , meerva.0$coef_beta[1,] ) ; beta_aug_vars   = rbind( beta_aug_vars  , meerva.0$var_beta[1,] ) ;  
      beta_vals  = rbind( beta_vals , meerva.0$coef_beta[2,] ) ; beta_val_vars   = rbind( beta_val_vars  , meerva.0$var_beta[2,] ) ;   
      gamma_bigs = rbind( gamma_bigs, meerva.0$coef_gamma[1,]) ; gamma_big_varjs = rbind( gamma_big_varjs, meerva.0$var_gamma[1,] ) ; 
                                                                 gamma_big_vars  = rbind( gamma_big_vars , meerva.0$var_gamma[2,] ) ;
      gamma_vals = rbind( gamma_vals, meerva.0$coef_gamma[2,]) ; gamma_val_vars  = rbind( gamma_val_vars , meerva.0$var_gamma[3,] ) ;
      if  (vmethod==3) { 
        beta_aug_1s = rbind( beta_aug_1s, meerva.1$coef_beta[1,] ) ; beta_aug_1_vars = rbind(beta_aug_1_vars, meerva.1$var_beta[1,] ) }
      if ((vmethod==3) & (simfam != "Cox")) {
        beta_aug_2s = rbind( beta_aug_2s, meerva.2$coef_beta[1,] ) ; beta_aug_2_vars = rbind(beta_aug_2_vars, meerva.2$var_beta[1,] ) } 
    } 
  }
  
  listname=list(simfam=simfam, mncor=mncor, vmethod=vmethod, jksize=jksize, compare=compare, nsims=nsims, seed=seed, nm=c(n,m),   
                beta=beta, alpha1=alpha1, alpha2=alpha2, bx3s1=bx3s1, bx3s2=bx3s2, bx12=bx12, sd=sd, mncor=mncor, sigma=sigma,
      beta_augs       = beta_augs      , beta_vals      = beta_vals ,
      beta_aug_vars   = beta_aug_vars  , beta_val_vars  = beta_val_vars, 
      gamma_bigs      = gamma_bigs     , gamma_vals     = gamma_vals,
      gamma_big_varjs = gamma_big_varjs, gamma_val_vars = gamma_val_vars, gamma_big_vars = gamma_big_vars )  

  if (vmethod==3) { 
    listtemp = list( beta_aug_1s = beta_aug_1s, beta_aug_1_vars = beta_aug_1_vars ) ;   
    listname=c( listname, listtemp ) 
  }
  if ((vmethod==3) & (simfam != "Cox")) { 
    listtemp = list( beta_aug_2s = beta_aug_2s, beta_aug_2_vars = beta_aug_2_vars ) ;   
    listname=c( listname, listtemp ) 
  }
  class(listname) <- c("meerva.sim")  
  return(listname)
} 

##############################################################################################################################
##############################################################################################################################
#=========== summarize simulation results =========================================================================;

#' Summarize Information for a meerva.sim Simulation Study Output Object
#'
#' @param object Output object from the simulations study program meerva.sim.block
#' @param short 0 to produce extensive output summary, 1 to produce only a table of biases and MSEs
#' @param round number of decimal places to round to in some tables, NA for R default
#' @param ... further arguments 
#'
#' @return A summary print 
#' 
#' @export
#' 
#' @seealso 
#'   \code{\link{meerva.sim.block}} 
#'
summary.meerva.sim <- function(object, short=0, round=NA, ...) {
  
  listnamex = deparse(substitute(object))                        
  simfam  = object$simfam   
  mncor   = object$mncor 
  vmethod = object$vmethod    
  jksize = object$jksize 
  compare= object$compare 
  nsims  = object$nsims 
  seed   = object$seed 
  nm     = object$nm
  beta   = object$beta
  alpha1 = object$alpha1 
  alpha2 = object$alpha2 
  bx3s1  = object$bx3s1
  bx3s2  = object$bx3s2
  bx12   = object$bx12 
  sd     = object$sd
  mncor  = object$mncor
  sigma  = object$sigma
  beta_augs   = object$beta_augs       
  beta_vals   = object$beta_vals 
  gamma_vals  = object$gamma_vals 
  beta_aug_vars   = object$beta_aug_vars 
  beta_val_vars   = object$beta_val_vars 
  gamma_val_vars  = object$gamma_val_vars 
  gamma_bigs      = object$gamma_bigs
  gamma_big_varjs = object$gamma_big_varjs 
  gamma_big_vars  = object$gamma_big_vars   
  if (vmethod==3) { beta_aug_1s     = object$beta_aug_1s ; 
                    beta_aug_1_vars = object$beta_aug_1_vars }
  if ((vmethod==3)&(simfam!="Cox")) { beta_aug_2s     = object$beta_aug_2s ; 
                                    beta_aug_2_vars = object$beta_aug_2_vars }

  cat( "\n") ;   
  cat( "======================= Simulation parameters ===================================================\n\n" ) ; 
  
#  names(listnamex)  = c("list name =") ; print(listnamex) ; cat("\n") ; 
  cat("list name = " , listnamex, "\n\n" ) ; 
  # names(simfam) = c("glm family for analysis") ; print( simfam ) ; cat( "\n" ) ;  
  cat(paste0("R glm family = ", simfam, "\n\n"))
  cat(paste0("mncor = ", mncor, "\n\n"))
  cat(paste0("VCOV method vmethod =  ", vmethod))  
  if (vmethod==0) { cat(paste0(" , JK")) 
  } else if (vmethod<=1) { cat(paste0(" , dfbeta"))
  } else if (vmethod==2) { cat(paste0(" , IJK")) 
  } else if (vmethod==3) { cat(paste(", (compare multiple VCOV methods)"))
  }   
  cat(paste0("\n\n"))
  if (vmethod %in% c(0,3) ) { cat(paste0("jksize =  ", jksize, "\n\n"))  }
  cat(paste0("Comparison group =  ", compare)) 
  if (compare==1) { cat(paste0(" , ful")) } else if (compare==0) { cat(paste0(" , non"))  } ; cat(paste0("\n\n"))
  cat(paste0("Number of simulations = ", nsims, "\n\n")) ;
  cat(paste0("seed = ", seed, "\n\n"))  
  cat(paste0("Full and val sample sizes of ", nm[1], " and ", nm[2], "\n\n"))
  print(rbind(beta)) ; cat( "\n") ; 
  print(rbind(alpha1)) ; cat( "\n") ; 
  print(rbind(alpha2)) ; cat( "\n") ; 
  print(rbind(bx3s1)) ;  cat( "\n") ; 
  print(rbind(bx3s2)) ; cat( "\n") ; 
  if ("bx12" %in% names(object))  {print(rbind(bx12)) ; cat( "\n") ; }
#  if ("sd"   %in% names(object))  
  cat(paste0(" SD = ", sd, "\n\n")) 
  cat(paste0(" mncor = ", mncor, "\n\n")) 
  if (!is.null(sigma)) {print(sigma) ; cat( "\n") ; } 
  
  #===================================================================================================
  #====================== Estimate Averages  =========================================================

#  beta
  if (simfam == "Cox") { beta = beta[2:length(beta)]  }
#  beta
  
  betadim = dim(object$beta_augs)[2]
  gammadim = dim(object$gamma_bigs)[2]
  
  if (gammadim < betadim) { betag = c( beta[1:(betadim-2)] , (beta[(betadim-1)] + beta[betadim]) ) 
  } else if (gammadim > betadim) { betag = c( beta[1:(betadim-2)], beta[(betadim-1)]/2, beta[(betadim-1)]/2, beta[betadim] ) 
  } else { betag = beta[1:gammadim]
  }
  
  betadim = dim(beta_vals)[2]
  gammadim = dim(gamma_bigs)[2]
  
  beta_aug_true   = subvec(beta_augs, beta)    ;  beta_aug_mse = colMeans(beta_aug_true^2  ) ; 
  beta_val_true   = subvec(beta_vals, beta)    ;  beta_val_mse = colMeans(beta_val_true^2  ) 
  if (vmethod==3) { beta_aug_1_true = subvec(beta_aug_1s, beta) ;  beta_aug_1_mse = colMeans(beta_aug_1_true^2  ) }
  if ((vmethod==3) & (simfam!="Cox")) { 
                  beta_aug_2_true = subvec(beta_aug_2s, beta) ;  beta_aug_2_mse = colMeans(beta_aug_2_true^2  ) }  

  gamma_big_true  = subvec(gamma_bigs, betag)  ;  gamma_big_mse = colMeans(gamma_big_true^2  )
  gamma_val_true  = subvec(gamma_vals, betag)  ;  gamma_val_mse = colMeans(gamma_val_true^2  ) 
  
  if (short == 0) {
    cat( "================================================================================================\n" ) ; 
    cat( "======================= Estimate averages ======================================================\n\n" ) ; 
    
    aves = rbind(beta=beta, beta_val=colMeans(beta_vals), beta_aug=colMeans(beta_augs))
    if  (vmethod==3) { aves = rbind(aves, beta_aug_1=colMeans(beta_aug_1s))  }
    if ((vmethod==3) & (simfam!="Cox")) { aves = rbind(aves, beta_aug_2=colMeans(beta_aug_2s))  }
    print(aves) ;     cat( "\n") ; 
    
    if (compare==1) { print( rbind(gamma_ful=colMeans(gamma_bigs ), gamma_val=colMeans(gamma_vals) ) ) 
    }   else        { print( rbind(gamma_non=colMeans(gamma_bigs ), gamma_val=colMeans(gamma_vals) ) ) }
    cat( "\n") ;     
  }
#======================= Bias ===================================================================; 
  
  biasbeta = rbind(beta_val_bias=colMeans(beta_val_true), beta_aug_bias=colMeans(beta_aug_true)) 
  if (vmethod==3) { biasbeta = rbind(biasbeta, beta_aug_1_bias=colMeans(beta_aug_1_true)) }
  if ((vmethod==3)&(simfam!="Cox")) { biasbeta = rbind(biasbeta, beta_aug_2_bias=colMeans(beta_aug_2_true)) }  

  if (compare==1) { biasgamma = rbind(gamma_ful_bias=colMeans(gamma_big_true ), gamma_val_bias=colMeans(gamma_val_true) ) 
  }   else        { biasgamma = rbind(gamma_non_bias=colMeans(gamma_big_true ), gamma_val_bias=colMeans(gamma_val_true) ) }
  
  if (short == 0) {
    cat( "======================= Bias ===================================================================\n\n" ) ; 
    print( biasbeta ) ;  cat( "\n") ; 
    print( biasgamma) ;  cat( "\n") ; 
  }
#======================= Standard Deviations ===================================================\n\n" ) ; 
 
  sdbeta = rbind( beta_val_sd =apply(beta_vals ,2,stats::sd), beta_aug_sd =apply(beta_augs, 2,stats::sd) )  
  if (vmethod==3) { sdbeta = rbind( sdbeta, beta_aug_1_sd =apply(beta_aug_1s, 2,stats::sd) )  } 
  if ((vmethod==3)&(simfam!="Cox")) { sdbeta = rbind( sdbeta, beta_aug_2_sd =apply(beta_aug_2s, 2,stats::sd) )  }  
  
  if (compare==1) { sdgamma = rbind( gamma_ful_sd=apply(gamma_bigs,2,stats::sd), gamma_val_sd=apply(gamma_vals,2,stats::sd) ) 
  }  else         { sdgamma = rbind( gamma_non_sd=apply(gamma_bigs,2,stats::sd), gamma_val_sd=apply(gamma_vals,2,stats::sd) ) }
  
  if (short == 0) {
    cat( "======================= Standard Deviations ===================================================\n\n" ) ;   
    print( sdbeta ) ;  cat( "\n") ; 
    print( sdgamma) ;  cat( "\n") ; 
  }
#======================= Average Standard Errors  ===============================================# 
  
  sebeta = rbind( beta_val_sd=colMeans(sqrt(beta_val_vars)), beta_aug_sd =colMeans(sqrt(beta_aug_vars)) ) 
  if (vmethod==3) { sebeta = rbind( sebeta, beta_aug_1_sd =colMeans(sqrt(beta_aug_1_vars)) ) }   
  if ((vmethod==3)&(simfam!="Cox")) { sebeta = rbind( sebeta, beta_aug_2_sd =colMeans(sqrt(beta_aug_2_vars)) ) }     
  
  if (compare==1) { segamma = sqrt( rbind( gamma_ful_sdj=colMeans(gamma_big_varjs), gamma_val_sd=colMeans(gamma_val_vars), gamma_ful_sd=colMeans(gamma_big_vars) ) ) 
  }   else        { segamma = sqrt( rbind( gamma_non_sdj=colMeans(gamma_big_varjs), gamma_val_sd=colMeans(gamma_val_vars), gamma_non_sd=colMeans(gamma_big_vars) ) ) }
  
  if (short == 0) {
    cat( "======================= Average Standard Errors  ===============================================\n\n" ) ;   
    print( sebeta  ) ; cat( "\n" )
    print( segamma ) ; cat( "\n" )
  }
  
  if (short == -1) {
    cat( "================================================================================================\n" ) ; 
    cat( "============= Check consistency of MSE, bias and var of estimates ==============================\n\n" ) ; 
  
    cat( "  beta_aug \n" ) ; 
    print( compmse( beta_aug_true ,  ivals = beta_augs ) ) ; cat( "\n") ;    ## is.vector(beta_val_true) ; is.matrix(beta_val_true)
    if  (vmethod==3) {  cat( "  beta_aug_1 \n" ) ; print( compmse( beta_aug_1_true ,  ivals = beta_aug_1s ) ) ; cat( "\n") }
    if ((vmethod==3)&(simfam!="Cox")) { cat("  beta_aug_2 \n" ) ; print( compmse( beta_aug_2_true ,  ivals = beta_aug_2s ) ) ; cat( "\n") }
    cat("  beta_val\n")  ;  print( compmse( beta_val_true ,  ivals = beta_vals ) ) ; cat( "\n") ;    
    if (compare==1) { cat("  gamma_ful\n") } else { cat("  gamma_non\n") }
    print( compmse( gamma_big_true ,  ivals = gamma_bigs ) ) ; cat( "\n") ;  
    cat( "  gamma_val\n" ) ; 
    print( compmse( gamma_val_true ,  ivals = gamma_vals ) ) ; cat( "\n") ;  
  }
  
#================================================================================================ 
#======================= Square Root MSEs =======================================================

  rmsebeta = rbind(beta_val_rmse=beta_val_mse, beta_aug_rmse=beta_aug_mse)  
  if  (vmethod==3) { rmsebeta = rbind(rmsebeta, beta_aug_1_rmse=beta_aug_1_mse) }  
  if ((vmethod==3)&(simfam!="Cox")) { rmsebeta = rbind(rmsebeta, beta_aug_2_rmse=beta_aug_2_mse) } 
  rmsebeta= sqrt(rmsebeta) ; 
  
  if (compare==1) { rmsegamma = sqrt( rbind(gamma_ful_rmse=gamma_big_mse, gamma_val_rmse=gamma_val_mse) ) 
  }  else         { rmsegamma = sqrt( rbind(gamma_non_rmse=gamma_big_mse, gamma_val_rmse=gamma_val_mse) )  }
  
  if (short == 0) {
    cat( "================================================================================================\n" ) 
    cat( "======================= Square Root MSEs =======================================================\n\n" ) ;  
    print( rmsebeta ) ;  cat( "\n") ; 
    print( rmsegamma) ;  cat( "\n") ; 
  }

  if (short == 0) {
    cat( "================================================================================================\n" ) 
    cat( "======================= Compare means and squared errors =======================================\n\n" ) ;
  
    beta_aug_val_mse   = beta_aug_true^2     - beta_val_true^2
    if (vmethod==3) { beta_aug_aug_1_mse = beta_aug_true^2 - beta_aug_1_true^2 } 
    if ((vmethod==3)&(simfam!="Cox")) { beta_aug_aug_2_mse = beta_aug_true^2 - beta_aug_2_true^2 } 
  
    cat( " Compare MSE beta_aug minus MSE bata_val (Paired t-test) \n" ) ;
    print( myttest( beta_aug_val_mse ) ) ; cat( "\n") ;
  
    if (vmethod==3) { cat( " Compare MSE aug minus MSE aug_1 (Paired t-test) \n" ) ;      
                    print( myttest( beta_aug_aug_1_mse ) ) ; cat( "\n") ;             } 
    if ((vmethod==3)&(simfam!="Cox")) { cat( " Compare MSE aug minus MSE aug_2 (Paired t-test) \n" ) ;      
                    print( myttest( beta_aug_aug_2_mse ) ) ; cat( "\n") ;             }
    
    cat( " Bias in beta_val \n" ) ;    print( myttest(beta_val_true) ) ; cat( "\n") 

    cat( " Bias in beta_aug \n" ) ;    print( myttest(beta_aug_true) ) ; cat( "\n") 
    
    if (vmethod==3) {
    cat( " Bias in beta_aug_1 \n" ) ;  print( myttest(beta_aug_1_true) ) ; cat( "\n") } 
    
    if ((vmethod==3) & (simfam!="Cox")) {
    cat( " Bias in beta_aug_2 \n" ) ;  print( myttest(beta_aug_2_true) ) ; cat( "\n") }
    
    if (compare==1) { cat( " Bias in gamma_ful\n") } else { cat( " Bias in gamma_non\n") }
    print( myttest(gamma_big_true) ) ; cat( "\n") ;
  }
  if (short == 0) {
    cat( "================================================================================================\n" ) 
    cat( "======== beta and gamma 95% Confidence Interval coverage probabilities =========================\n\n" ) ; 

    cat( " Coverage for beta_aug 95% CI \n" ) ;
    print( coverage(beta_augs, beta_aug_vars, beta) ) ; cat( "\n") ;## estimates=beta_augs ; vars = beta_aug_vars ; 

    if (vmethod==3) {   cat( " Coverage for beta_aug_1 95% CI \n" ) ;
      print( coverage(beta_aug_1s, beta_aug_1_vars, beta) ) ; cat( "\n") ; }
    
    if ((vmethod==3) & (simfam!="Cox")) { cat( " Coverage for beta_aug_2 95% CI \n" ) ;
      print( coverage(beta_aug_2s, beta_aug_2_vars, beta) ) ; cat( "\n") ; }
  
    cat( " Coverage for beta_val 95% CI \n" ) ;  
    print( coverage(beta_vals, beta_val_vars, beta) ) ; cat( "\n") ;
    
    if (compare==1) { cat( " Coverage for gamma_ful 95% CI \n" ) } else { cat( " Coverage for gamma_non 95% CI \n" ) }
    print( coverage(gamma_bigs, gamma_big_varjs, betag) ) ; cat( "\n") ;  

    if (vmethod==3) { cat( " Comparison of 95% coverages between beta_aug (1) and beta_aug_1 (2) \n" ) ;  
      print( coverage2(beta_augs, beta_aug_vars, beta_aug_1s, beta_aug_1_vars, beta) ) ; cat( "\n") ; }  

    if ((vmethod==3)&(simfam!="Cox")) { cat( " Comparison of 95% coverages between beta_aug (1) and beta_aug_2 (2) \n" ) ;  
      print( coverage2(beta_augs, beta_aug_vars, beta_aug_2s, beta_aug_2_vars, beta) ) ; cat( "\n") ; }  
  }

  if ( (short==1)  &  (dim(biasbeta)[2] == dim(biasgamma)[2]) ) { 
    cat( "=========== Bias, SD and MSE for estimators ===================================================\n\n" )  
    if (!is.na(round)) { tabel = rbind(round(biasbeta,round), round(biasgamma,round), round(sdbeta,round),
                                       round(sdgamma,round), round(rmsebeta,round), round(rmsegamma,round))  
    } else { tabel = rbind(biasbeta, biasgamma, sdbeta, sdgamma, rmsebeta, rmsegamma) } 
    print(tabel) ; cat("\n") ; 
  } else if  (short==1)  { 
    cat( "=========== Bias, SD and MSE for estimators ===================================================\n\n" )  
    if (!is.na(round)) { tabel1 = rbind(round(biasbeta,round), round(sdbeta,round), round(rmsebeta,round)) 
    } else { tabel1 = rbind(biasbeta, sdbeta, rmsebeta) }
    if (!is.na(round)) { tabel2 = rbind(round(biasgamma,round), round(sdgamma,round), round(rmsegamma,round)) 
    } else { tabel2 = rbind(biasgamma, sdgamma, rmsegamma) } 
    print(tabel1) ; print(tabel2) ; cat("\n") ; 
  }
  
  cat( "================================================================================================\n\n" )  
} 

##############################################################################################################################
##############################################################################################################################

#' Calculate relative differences
#'
#' @param a Object 1 for comparison 
#' @param b Object 2 for comparison 
#'
#' @return A object with relative differences
#'
reldif = function(a,b) { 2*(a-b)/(a+b) }

#=====================================================================================================#
#=========Define function to subtract vector from each row of a matrix ===============================#

#' Subtract a vector from each row of a matrix 
#'
#' @param x A matrix
#' @param v A vector of length dim[x](2)
#'
#' @return A matrix 
#'
subvec = function(x,v) { return( t(t(x) - v) ) }

#######################################################################################################
#========== check that MSE = Bios^2 + var ============================================================#

#' Comapre bias, var and MSE 
#'
#' @param ibias Matrix of bias's
#' @param ivals Matrix of var's
#'
#' @return A matrix 
#'
compmse = function( ibias, ivals ) {
  bias2  = colMeans(ibias)^2                           ## bias square is.vector(colMeans(ibias)^2) 
  var_   = apply(ivals,2, stats::var) * (dim(ivals)[1] - 1)/dim(ivals)[1]   ## variance   is.vector(varg) 
  mse_   = colMeans(ibias^2)                           ## MSE   -  is.vector(mseg) 
  diff   = mse_ - bias2 - var_ 
  lengths = cbind( length(bias2) , length(var_) , length(mse_) , length(diff ) ) ; 
  summary =  rbind( bias2, var_, mse_ , diff)  
  return( summary )              
} 

#########################################################################################################
#========== Define function to do 1-sample t-test on each column of a matrix ===========================#

#' A simple summary description 
#'
#' @param x A data matrix 
#' @param beta0 A vector for H0
#'
#' @return A matrix
#' 
#' @importFrom stats qt pt 
#'
myttest = function(x, beta0=NULL) {
  if (length(beta0) == 0) { beta0=c(rep(0,dim(x)[2])) }  
  ns = dim(x)[1] ;                    
  ns = colSums(!is.na(x))
  means = colMeans(x) 
  vars  = (colMeans (x^2) - colMeans(x) ^2 ) * ns / (ns-1)
  sds   = sqrt(vars) 
  ses   = sds / sqrt(ns) 
  ttests = (means - beta0) / ses ; 
  lcls  = means - qnorm(0.975) * ses 
  ucls  = means + qnorm(0.975) * ses 
  lcls  = means - qt(0.975, ns-1) * ses 
  ucls  = means + qt(0.975, ns-1) * ses   
  ps = 2*pt(-abs(ttests), ns-1)              # 2 * pt(-2,100)
  #  return(list(n=ns, mean=means, sd=sds, lcl=lcls, ucl=ucls, t=ttests,p=ps))
  return(cbind(n=ns, mean=means, beta0=beta0, sd=sds, lcl=lcls, ucl=ucls, t=ttests, p=ps))
} 

##########################################################################################################
#========== calcualte CIs for CI coverage probablities from simulation study ============================#

#' Calculate coverage probabilties 
#'
#' @param estimates A matrix of estimates 
#' @param vars A matrix of varainces
#' @param beta The beta under H0
#' @param round The decimal places for rounding
#'
#' @return 95% CI coverage probabilities 
#' 
#' @importFrom stats pbinom 
#'
coverage <- function(estimates, vars, beta, round=3) {
  if (length(beta) == 0) { beta=c(rep(0,dim(estimates)[2])) }
  lcls = estimates - qnorm(0.975) * sqrt(vars)
  ucls = estimates + qnorm(0.975) * sqrt(vars)
  include = ( (subvec(lcls, beta) <= 0) & ( subvec(ucls, beta) >= 0) ) * 1 
  coverage     = colMeans(include) ; 
  n_ = colSums( !is.na(include) )   ; 
  #  coverage_lcl = coverage - qnorm(0.975) * sqrt( (coverage * colMeans(1-include))  / dim(include)[1] )  
  #  coverage_ucl = coverage + qnorm(0.975) * sqrt( (coverage * colMeans(1-include))  / dim(include)[1] )  
  coverage_lcl = coverage - qnorm(0.975) * sqrt( (coverage * (1 - coverage))  / n_ )  
  coverage_ucl = coverage + qnorm(0.975) * sqrt( (coverage * (1 - coverage))  / n_ )  
  coverage_lcl =  coverage_lcl * (coverage_ucl >= 0)  +  0 * (coverage_ucl < 0)    
  coverage_ucl =  coverage_ucl * (coverage_ucl <= 1)  +  1 * (coverage_ucl > 1)
  coverage     =  round(coverage, round)
  coverage_lcl =  round(coverage_lcl, round) 
  coverage_ucl =  round(coverage_ucl, round)    
  return( t( rbind(n_, coverage, coverage_lcl,  coverage_ucl) ) ) 
}

###############################################################################################
#========== compare coverage probabilities ===================================================#

#' Compare coverage probabilities between two estimators on matched simulations 
#'
#' @param estimates1 Matrix of estimates for 1 
#' @param vars1 Matrix of vars for 1 
#' @param estimates2 Matrix of estimates for 2 
#' @param vars2 Matrix of vars for 2 
#' @param beta The beta under H0
#' @param round The decimal places for rounding
#' 
#' @return Comparison of coverage probabilities between two estimators.
#'
coverage2 <- function(estimates1, vars1, estimates2, vars2, beta, round=3) {
  if (length(beta) == 0) { beta=c(rep(0,dim(estimates1)[2])) }
  lcls1 = estimates1 - qnorm(0.975) * sqrt(vars1)
  ucls1 = estimates1 + qnorm(0.975) * sqrt(vars1)
  include1 = ( (subvec(lcls1, beta) <= 0) & ( subvec(ucls1, beta) >= 0) ) * 1 
  coverage1     = colMeans(include1) ; 
  n1 = colSums( !is.na(include1) )   ; 
  
  lcls2 = estimates2 - qnorm(0.975) * sqrt(vars2)
  ucls2 = estimates2 + qnorm(0.975) * sqrt(vars2)
  include2 = ( (subvec(lcls2, beta) <= 0) & ( subvec(ucls2, beta) >= 0) ) * 1 
  coverage2     = colMeans(include2) ;                          # coverage1[1:10,]
  n2 = colSums( !is.na(include2) )   ; 
  
  sum1b2 = colSums(include1 > include2) ;
  sum2b1 = colSums(include2 > include1) ;  
  ndif = sum1b2 + sum2b1 ; 
  pval = 2*pbinom((sum1b2 <sum2b1)*sum1b2 + (sum1b2 >= sum2b1)*sum2b1 ,ndif,0.5) ; ## pval = pbinom(0, 5, 0.5)   
  return( rbind(sum1b2=sum1b2, sum2b1=sum2b1, ndif=ndif, n=n1, pval=pval) )
}

################  BEGIN SIMULATION PLOT FUNCTION     ############################################################################

#' Plot results for meerva.sim output object generated by meerva.sim.block function
#'
#' @param x A meerva.sim class object 
#' @param violin 1 to produce a violin plot instead of a boxplot
#' @param ... further arguments
#'
#' @return This displays a plot 
#' 
#' @export
#'
plot.meerva.sim = function(x, violin=0, ...) {
  
  beta = x$beta 
  if (x$simfam == "Cox") { beta = beta[2:length(beta)]  }
  
  betadim = dim(x$beta_augs)[2]
  gammadim = dim(x$gamma_bigs)[2]
  
  if (gammadim < betadim) { betag = c( beta[1:(betadim-2)] , (beta[(betadim-1)] + beta[betadim]) ) 
  } else if (gammadim > betadim) { betag = c( beta[1:(betadim-2)], beta[(betadim-1)]/2, beta[(betadim-1)]/2, beta[betadim] ) 
  } else { betag = beta[1:gammadim]
  }
  
  fun01 = function(x, item) { 
    if        (item==1) { 
      beta_augs = t(t(x$beta_augs) - beta) 
      betas = as.data.frame(beta_augs) 
    } else if (item==2) { 
      beta_vals = t(t(x$beta_vals) - beta)
      betas = as.data.frame(beta_vals) 
    } else if (item==3) { 
      gamma_fuls = t(t(x$gamma_bigs) - betag)
      betas = as.data.frame(gamma_fuls) 
    } else if (item==4) { 
      gamma_nons = t(t(x$gamma_bigs) - betag)
      betas = as.data.frame(gamma_nons) 
    }
    
    betas <- tidyr::pivot_longer(betas, names_to="Group", values_to="Value", cols= names(betas))
    
    betas$Group <- ifelse(betas$Group == "(Intercept)", "B0",
                   ifelse(betas$Group %in% c("x_valx1", "xs_nonx1s", "xs_fulx1s" ), "B1",
                   ifelse(betas$Group %in% c("x_valx2", "xs_nonx2"  , "xs_nonx2s", "xs_fulx2" , "xs_fulx2s"), "B2",
                   ifelse(betas$Group %in% c("xs_fulx3s2", "xs_nonx3s2"), "B3d",                        
                   ifelse(betas$Group %in% c("x_valx3", "xs_nonx3s", "xs_nonx3s1", "xs_fulx3s", "xs_fulx3s1"), "B3",
                   ifelse(betas$Group %in% c("x_valx4", "xs_nonx4"  , "xs_nonx4s", "xs_fulx4"  , "xs_fulx4s"), "B4" , "B5"))))))

    if        (item==1) { betas$var = "beta_aug"  
    } else if (item==2) { betas$var = "beta_val"  
    } else if (item==3) { betas$var = "gamma_ful" 
    } else if (item==4) { betas$var = "gamma_non" 
    }  
    return(betas) 
  }  
  
  beta_aug  = fun01(x,1)
  beta_val  = fun01(x,2)
  if (x$compare==1) { 
    gamma_ful = fun01(x,3)
    plotdata = rbind(beta_aug, beta_val, gamma_ful) 
  } else { 
    gamma_non = fun01(x,4)
    plotdata = rbind(beta_aug, beta_val, gamma_non) 
  } 
  
  Group = plotdata$Group
  Value = plotdata$Value
  var   = plotdata$var 
  
  tsize = 12 
  if (violin==0) {
#     ggplot2::geom_boxplot(notch = TRUE) +
#      ggplot2::ggplot(data = plotdata, ggplot2::aes(x = Group, y = Value, fill = var)) +      
      ggplot2::ggplot(mapping=ggplot2::aes(x = Group, y = Value, fill = var) ) +      
      ggplot2::geom_hline(yintercept = 0, col = "lightgrey") +
      ggplot2::geom_boxplot() +
      ggplot2::theme_classic() +
      ggplot2::theme(panel.border = ggplot2::element_rect(color = "black", fill = NA), 
            axis.text = ggplot2::element_text(size  = tsize, color = "black"), 
            axis.title = ggplot2::element_text(size = tsize),
            plot.title = ggplot2::element_text(size = tsize)) +
      ggplot2::scale_fill_manual(values = c("blue", "darkgreen", "red"), name = "Model") +
      ggplot2::guides(fill = ggplot2::guide_legend(label.theme = ggplot2::element_text(size = tsize))) +
      ggplot2::xlab("Coefficient") +
      ggplot2::ylab("Estimate - True") +
      ggplot2::ggtitle("Simultion results for augmented estimates")
  } else {
    ggplot2::ggplot(data = plotdata, ggplot2::aes(x = Group, y = Value, fill = var)) +
      ggplot2::geom_hline(yintercept = 0, col = "lightgrey") +
      ggplot2::geom_violin() +
      ggplot2::theme_classic() +
      ggplot2::theme(panel.border = ggplot2::element_rect(color = "black", fill = NA), 
                     axis.text = ggplot2::element_text(size  = tsize, color = "black"), 
                     axis.title = ggplot2::element_text(size = tsize),
                     plot.title = ggplot2::element_text(size = tsize)) +
      ggplot2::scale_fill_manual(values = c("blue", "darkgreen", "red"), name = "Model") +
      ggplot2::guides(fill = ggplot2::guide_legend(label.theme = ggplot2::element_text(size = tsize))) +
      ggplot2::xlab("Coefficient") +
      ggplot2::ylab("Estimate - True") +
      ggplot2::ggtitle("Simultion results for augmented estimates")
  }
  #  return()
} 

##############################################################################################################################
#=================================  END  ====================================================================================#
##############################################################################################################################

