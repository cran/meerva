#==========================================================================================================#
#===== meerva_210412.R                                                                               ======#
#===== Analysis of MEasuremnt ERror data with VAlidation subsample                                   ======#
#==========================================================================================================#
#' Analysis of Data with Measurement Error Using a Validation Subsample
#' 
#' @description
#' The meerva package is designed to analyze data with measurement error when there is a
#' validation subsample randomly selected from the full sample.  The method assumes
#' surrogate variables measured with error are available for the full sample,
#' and reference variables measured with little or no error are available for this randomly
#' chosen subsample of the full sample.  Measurement errors may be differential or
#' non differential, in any or all predictors (simultaneously) as well as outcome.
#' The "augmented" estimates derived by meerva are based upon the multivariate correlation between regression
#' models based upon the reference variables and the surrogate variables in the validation subset.  Because the
#' validation subsample is chosen at random whatever biases are imposed by measurement error, non-differential
#' or differential, are reflected in this correlation and can be used to derives estimates for the reference
#' variables using data from the whole sample.
#'
#' Intuitively one expects there to be at least one surrogate for each reference variable but the method is based
#' upon multivariate correlations and therefore also works if there are more or fewer surrogate than reference variables.
#' The package fits linear, logistic or Cox regression models individually to the reference variales and to the
#' surrogate varibles, then combines the results to descibe a model in terms of the reference variables based
#' upon the entire dataset.
#'
#' @param x_val A matrix object including reference predictor variables (and predictors "without" error) in validation subsample.
#'   This and other x_ matrices must not include any missing values (NA).
#'   All data vectors and matrices must be numerical.  For categorical variables one
#'   should first construct corresponding numerical variables to represent these categories.
#' @param y_val A vector object for the reference outcome variable in validation subsample.
#'   This and other y_ vectors must not include any missing values (NA).
#' @param xs_val A matrix object including surrogate predictors (and predictors "without" error) in validation subsample
#' @param ys_val A vector object for the surrogate outcome variable in validation sample.
#' @param xs_non A matrix object including surrogate predictors (and predictors "without" error) in NON validation data
#' @param ys_non A vector object for the surrogate outcome variable in the NON validation sample.
#' @param e_val A vector object for the survival data reference event outcome variable in validation subsample.
#'   This and the other e_ vectors are optional.
#'   The e_ vectors are required when analyzing survival data based upon an underlying Cox regression model (survival package).
#'   This and other e_ vectors must not include any missing values (NA).
#' @param es_val A vector object for the survival data surrogate event outcome variable in validation subsample.
#' @param es_non A vector object for the survival data surrogate event outcome variable in NON validation data.
#' @param id_val A vector object identifying clusters in case of multiple records per subject in the validation subsample.
#'   This and id_non are optional. They must not include any missing values (NA).
#'   No subjects should be included in both the validation subsample and the NON validation data.
#' @param id_non A vector object identifying clusters in case of multiple records per subject in the NON validation data.
#' @param weights_val A vector object with weights used in model fit of the validation subsample.
#'   This can be used, for example, to describe inverse sampling probability weights.
#'   Note, when fitting the "binomial" or logistic model, weights for weights_val and weights_non
#'   must be integer.  This is a restriction of the glm.fit routine called from meerva.  The user may rescale or round the
#'   weights to achieve integers.  By using robust variance estimates meerva provides correct variance estimates.
#' @param weights_non A vector object with weights used in model fit of the NON validation subsample.
#'   This and weights_val, can be used, for example, to down weight records from patients with multiple records.
#' @param familyr The family for the underlying regression model amongst "binomial", "gaussian" and "Cox".
#'   If not specified the program chooses among these three based upon a simple data inspection. 
#'   In principle, though not (yet) implemented here, the regression model for the reference variables may 
#'   be of a different type than for the surrogate variables.  For example the reference outcome could be yes/no 
#'   in nature while the surrogate outcome could be a numeric, and the method would continue to work.
#' @param vmethod Method for robust estimation of variance covariance matrices needed for calculation of the augmented estimates (beta aug).
#'   0, 1 or 2 determines JK (slow), IJK using dfbeta of glm or coxph, or IJK using an alternate formula for dfbeta.
#'   Recommendations:  For "gaussian" use 1, for "Cox" use 1 for speed and 0 for accuracy,
#'   and for "binomial" use 2 for speed, and 0 for accuracy.
#' @param jksize Number of elements to leave out number in each cycle of the grouped jackknife
#'   for non validation data.  The default is 0 where the program chooses jksize
#'   so that the number of leave out groups is about validation subsample size.
#'   For the grouped jackknife the program randomly sorts the non validation subsample.  To get
#'   the exact same results twice one can set the seed for the random generator with the
#'   statement set.seed(seed) for some value of seed, and to get a "random" seed one can
#'   first run the statement seed = round(runif(1)*1000000000) .
#' @param compare 1 to compare gamma_val with gamma_ful (default) or 0 with gamma_non.
#'   Comparisons of gamma_val with gamma_ful is consistent with the principle of the
#'   validation set being a subsample of the entire dataset.  This assures the 
#'   correlations between gamma_val and beta_val are representative of what
#'   should be the case based upon the whole dataset.  If there were an 
#'   external validation sample where one could be reasonably certain 
#'   that the correlation between gamma_val and beta_val would be
#'   representative then one could also use this method.
#'   
#' @details
#' As currently implemented the package requires the data to be input as
#' vectors and matrices with no missing values (NA).
#' All data vectors and matrices must be numerical.  For categorical variables one
#' should first construct corresponding numerical variables to represent these categories.
#' Note, variables thought of as measured without error should be included in both the reference variable set and
#' the surrogate variable set.  Such variables may be thought of as perfect surrogates.  This applies for both
#' outcome variables and predictor variables.  For the Cox model both the time to event and the event indicator
#' may be measured with error.
#'
#' The length of the vectors for the validation subsample must all be the same, and be the same as the number of rows
#' in the predictor matrices for the validation subsample.  Data for sample elements not included in the validation
#' subsample are referred to as NON validation data and are to be included in separate vectors and matrix.  Here, too, the
#' length of all vectors must be the same as number of rows in the predictor matrix.  The columns in the data matrix for the
#' validation subsample surrogates must be logically the same as the columns in the data matrix for the
#' NON validation surrogates.
#'
#' The data for analysis may include weights, for example to account for non identical
#' sampling probabilities when selecting the subsample, by taking weights as the inverse of these probabilities.
#' The data may include cluster identifiers in case of multiple observations on study participants.
#' Weights may also be used to lessen the influence of individuals with multiple observations.
#'
#' Internally the analysis uses robust variance estimation which accounts for deviations from the usual regression
#' model assumptions for the surrogate variables, and accounts for multiple observations per patient.
#'
#' This package came out of our work analyzing electronic health records data, where different sources, e.g diagnosis codes
#' and natural language processing, may provide different surrogate variables.  Reference variables were obtained by
#' manual chart review.  For our datasets to date with tens of thousands of patients the analyses take a few seconds
#' when run on a PC.
#'
#' In the examples we generate simulated data of the form expected for input, call the main program,
#' and summarize the output.
#'
#' @return
#' meerva.fit returns an object of class meerva which contains
#' the augmented estimates based upon the full data set accounting for measurement error,
#' estimates based upon reference variables from the validation subsample,
#' estimates based upon the surrogate variables from the whole sample,
#' along with robust variance-covariances matrix estimates for these estimates.
#' This meerva class list contains the following objects.
#'
#' Call     --- The call used to invoke meerva.fit.
#'
#' FitInput --- A list with
#'
#' --- familyr --- The type of regression model fit to the data.
#'
#' --- compare --- The input parameter compare, 1 to compare the validation
#'            data with the whole dataset, or 0 to compare with the NON validation data.
#'
#' --- comparec --- A short text interpretation of compare.
#'
#' --- vmethod --- The method used to estimate the variance-covariance
#'            matrices needed for calculation of the estimates.
#'
#' --- vmethodc --- A short text description of vmethod.
#'
#' --- n_val    --- The number of observations in the validation subsample.
#'
#' --- n_ful    --- The number of observations in the whole dataset.
#'
#' --- n_val_id --- The number of clusters identified by id_val in the validation subsample.
#'
#' --- n_ful_id --- The number of clusters identified by id_val and id_non in the whole dataset.
#'
#' --- dim_beta --- The number of parameters in the regression model
#'           for reference variables including a possible intercept.
#'
#' --- dim_gamma --- The number of parameters in the regression
#'           model for surrogate variables including a possible intercept.
#'
#' names_x  ---  The reference variable predictors used in analysis.
#'
#' names_xs ---  The surrogate variable predictors used in analysis.
#'
#' names_y  ---  The reference outcome variable used in analysis.
#'
#' names_ys ---  The surrogate outcome variable used in analysis.
#'
#' coef_beta --- The regression parameter estimates for the reference
#'    variables including both beta_val based upon the reference variables alone
#'    (available only in the validation subsample) and beta_aug, the augmented 
#'    estimates based upon the reference variables in the validation subsample
#'    augmented by the surrogate variables in the whole dataset.
#'
#' coef_gamma --- The regression parameter estimates for the surrogate
#'    variables for both gamma_val derived using dataset elements included
#'    in the validation subsample, and
#'    either gamma_ful or gamma_non, derived using either the whole 
#'    sample or the NON validation data.
#'
#' var_beta --- Robust variance estimates for coef_beta, which are
#'    also included in vcov_beta and vcov_beta_val.
#'
#' var_gamma --- Robust variance estimates for coef_gamma, which are
#'    also included in vcov_gamma.
#'
#' vcov_beta_aug --- Robust variance-covariance estimates for beta_aug of coef_beta.
#'
#' vcov_beta_val ---  Robust variance-covariance estimates for beta_val of coef_beta.
#'
#' vcov_beta_val_naive --- Naive variance-covariance estimates for beta_val
#'    of coef_beta obtained without any consideration of clustering optionally
#'    described by input parameters id_val and id_non.
#'
#' vcov_gamma_ful --- Robust variance-covariance estimates for gamma_ful of coef_gamma.
#'
#' or vcov_gamma_non --- Robust variance-covariance estimates for gamma_non of coef_gamma.
#'
#' vcov_gamma_ful_naive --- Naive variance-covariance estimates for
#'    gamma_ful of coef_gamma obtained without any consideration of
#'    clustering optionally described by input parameters id_val and id_non.
#'
#' or vcov_gamma_non_naive --- Like vcov_gamma_ful_naive but for gamma_non.
#'
#' omega --- The robust covariance estimate between beta_val and either
#'    gamma_ful or gamma_non, which is integral for derivation of beta_aug.
#'
#' omega_cor --- The robust correlation estimate between beta_val and
#'    either gamma_ful or gamma_non, which reflects the relative amount
#'    of information on reference variable estimates contained in the
#'    surrogate variables.
#'
#' kappa --- The robust variance covariance estimate of either
#'    (gamma_val - gamma_ful) or (gamma_val - gamma_non), which is
#'    integral for derivation of beta_aug.
#'
#' @author Walter Kremers (kremers.walter@mayo.edu)
#'
#' @seealso
#'   \code{\link{meerva.sim.block}} , \code{\link{meerva.sim.nrm}} , \code{\link{meerva.sim.brn}} , \code{\link{meerva.sim.cox}}
#'
#' @references
#'  Chen Y-H, Chen H.  A Unified Approach to Regression Analysis under Double-Sampling Designs.
#'  Journal of the Royal Statistical Society. Series B (Statistical Methodology) , 2000 (62) 449-460.
#'
#'  Chen Y-H. Cox regression in cohort studies with validation sampling.
#'  Journal of the Royal Statistical Society. Series B (Statistical Methodology), 2002 64, 51-62.
#'
#'  Wang X, Wang QH. Semiparametric linear transformation model with differential measurement error
#'  and validation sampling. J Multivariate Anal. 2015;141:67-80.
#'
#'  Tong JY, Huang J, Chubak J, et al. An augmented estimation procedure for EHR-based association
#'  studies accounting for differential misclassification. J Am Med Inform Assn. 2020;27(2):244-253.
#'
#' @export
#' 
#' @importFrom stats glm glm.fit glm.control binomial gaussian vcov dfbeta cov cor residuals 
#'
#' @examples
#' #======================================================
#' 
#' # Simulate logistic regression data with measurement error
#' simd = meerva.sim.brn(n=4000, m=400,
#'      beta = c(-0.5, 0.5, 0.2, 1, 0.5) , 
#'      alpha1 = c(0.95, 0.90, 0.90, 0.95) , 
#'      alpha2 = c(0.98,0.94,0.95,0.95) , 
#'      bx3s1 = c(0.05, 0, 0, NA, NA) , 
#'      bx3s2 = c(NA,NA,NA) )
#' 
#' # Read the simulated data to input data format
#' x_val  = simd$x_val
#' y_val  = simd$y_val
#' xs_val = simd$xs_val
#' ys_val = simd$ys_val
#' xs_non = simd$xs_non
#' ys_non = simd$ys_non
#' 
#' # Analyze the data
#' brn.me = meerva.fit(x_val, y_val, xs_val, ys_val, xs_non, ys_non)
#' summary(brn.me)
#' 
#' #======================================================
#' 
#' # Simulate linear regression data with measurement error
#' simd = meerva.sim.nrm(n=4000, m=400,
#'      beta=c(-0.5,0.5,0.2,1,0.5),
#'      alpha1=c(-0.05,0.1,0.05,0.1), 
#'      alpha2=c(0.95,0.91,0.9,0.9),
#'      bx3s1= c(0.05, 0, 0, NA, NA), 
#'      bx3s2=c(1.1,0.9,0.05),
#'      sd=5)
#'
#' # Read the simulated data to input data format
#' x_val  = simd$x_val
#' y_val  = simd$y_val
#' xs_val = simd$xs_val
#' ys_val = simd$ys_val
#' xs_non = simd$xs_non
#' ys_non = simd$ys_non
#'
#' # Analyze the data
#' nrm.me = meerva.fit(x_val, y_val, xs_val, ys_val, xs_non, ys_non)
#' summary(nrm.me)
#'
#' #======================================================
#' # Simulate Cox regression data with measurement error
#' simd = meerva.sim.cox(n=4000, m=400,
#'      beta   = c(-0.5, 0.5, 0.2, 1, 0.5) ,
#'      alpha1 = c(0.95,0.90,0.90,0.95)  ,
#'      alpha2 = c(0.98,0.94,0.94,0.98) ,
#'      bx3s1  = c(0.05,0,0,NA,NA) ,
#'      bx3s2  = c(1.1, NA, NA) ,
#'      sd=0.1)
#'
#' x_val  = simd$x_val
#' y_val  = simd$y_val
#' xs_val = simd$xs_val
#' ys_val = simd$ys_val
#' xs_non = simd$xs_non
#' ys_non = simd$ys_non
#' e_val  = simd$e_val
#' es_val = simd$es_val
#' es_non = simd$es_non
#'
#' cox.me = meerva.fit(x_val, y_val, xs_val, ys_val, xs_non, ys_non,
#'                     e_val, es_val, es_non)
#' summary(cox.me)
#'
#' #======================================================
#'
meerva.fit = function(x_val, y_val, xs_val, ys_val, xs_non, ys_non,
                  e_val=NULL, es_val=NULL, es_non=NULL, id_val=NULL, id_non=NULL,
                  weights_val=NULL, weights_non=NULL, familyr=NULL, vmethod=NULL, jksize=0, compare=1) {

  if ( is.null(familyr)) {familyr = NA}
  if ( is.null(vmethod)) { vmethod = NA }
  if ( is.na(familyr) & (!is.null(e_val)) ) { familyr = "Cox" }
  if ( !is.null(e_val) ) { familyr = "Cox" }
  if ( is.na(familyr) & ( length(table(ys_non)) == 2 ) ) { familyr = "binomial" }
  if ( is.na(familyr) & ( length(table(ys_non)) >  2 ) ) { familyr = "gaussian" }
  if ( familyr == "poison" ) { familyr = "poisson" }
  if ( familyr %in% c("gausian", "gause", "normal") ) { familyr = "gaussian" }
  if ( familyr %in% c("COX", "cox", "coxph") ) { familyr = "Cox" }
  if (is.null(vmethod) | (is.na(vmethod) )) {
    if (familyr == "binomial" ) { vmethod = 2
    } else                      { vmethod = 1
    }
  }
  if (compare %in% c("ful", "full")) { compare = 1 } else if (compare=="non") { compare = 0 }
  if ( ((familyr=="Cox")  &  (vmethod==2)) | ((familyr=="Cox")  &  (vmethod=="ijk")) ) {vmethod=1}
  if (vmethod %in% c("dfbeta")) {vmethod=1} else if (vmethod %in% c("ijk")) {vmethod=2} else if (vmethod %in% c("jk")) {vmethod=0}
  
  if ((vmethod==2) & (!is.null(weights_val))) { 
    vmethod=1 
    print("Warning, when assigning weights vmethod=2 is changed to vmethod=1.") 
  }

  ## get number of predictors for models, inclding intercepts
  if (familyr == "Cox") {
    if (is.vector(x_val))  { dim_mod  = 1
    } else {                 dim_mod  = dim(x_val )[2] }        ## number of predictors for reference set
    if (is.vector(xs_val)) { dim_mods = 1
    } else {                 dim_mods = dim(xs_val)[2] }        ## number of predictors for surrogate set
  } else {
    if (is.vector(x_val))  { dim_mod  = 2
    } else {                 dim_mod  = dim(x_val )[2] + 1 }
    if (is.vector(xs_val)) { dim_mods = 2
    } else {                 dim_mods = dim(xs_val)[2] + 1 }
  }
  dim_modt = dim_mod + dim_mods

  ## sort matrices, vectors by id_val and id_non to be able to easly leave out clusters for IJK and JK =================
  if ( (!is.null(id_val)) != (!is.null(id_non)) ) {print("Warning, only one of id_val and id_non is assinged!")}
  if ( (!is.null(id_val)) &  (!is.null(id_non)) )  {
    order_val = order(id_val)
    order_non = order(id_non)
    id_val = id_val[order_val]
    id_non = id_non[order_non]
    x_val  = x_val [order_val,]
    y_val  = y_val [order_val]
    if (!is.null(weights_val)) { weights_val = weights_val[ order_val ] }
    xs_val = xs_val[order_val,]
    ys_val = ys_val[order_val]
    xs_non = xs_non[order_non,]
    ys_non = ys_non[order_non]
    if (!is.null(weights_non)) { weights_non = weights_non[ order_non ] }
    if (!is.null(e_val)) {
      e_val  = e_val[order_val]
      es_val = es_val[order_val]
      es_non = es_non[order_non] ; # length(es_non)
    }
  }

  n_val = length(y_val)             ## val sample size
  n_non = length(ys_non)            ## non val sample size
  n_ful = n_val + n_non             ## full sample size

  ## set up ful datasets ==================================
  if (is.vector(xs_val) | is.vector(xs_non)) { xs_ful = c( as.vector(xs_val), as.vector(xs_non) )
  } else                                     { xs_ful = rbind(xs_val, xs_non)  }
  ys_ful = c(as.vector(ys_val), as.vector(ys_non) )
  if (!is.null(id_val))  {id_ful = c(as.vector(id_val), as.vector(id_non) ) }  ; # table(table(id_val))
  if (familyr== "Cox") { es_ful  = c(as.vector(es_val), as.vector(es_non) ) }
  if (!is.null(weights_non)) { weights_ful = c( as.vector(weights_val), as.vector(weights_non) )
  } else { weights_ful=NULL }

  ## set up vactor id_ful_u as consecutive, possibly repeating, intergers as id, e.g. 1, 1, 2, 3, 3, 3, 4,... n for ful dataset =======
  if ( (!is.null(id_val)) & (!is.null(id_non)) )  {
    n_val_id = as.integer( table(table(id_val))[1] )    ## number of levels in id_val # as.vector(n_val_id) ; as.vector class(n_val_id[1])
    n_non_id = as.integer( table(table(id_non))[1] )    ## number of levels in id_non

    id_val_u = c(rep(0,n_val))    ## working id_val with consecutive, possibly repeating, intergers
    for (i_ in 1:n_val) {
      if (i_ == 1) { count = 1
      } else if (id_val[i_] != last) { count = count + 1 }
      id_val_u[i_] = count
      last = id_val[i_]
    }
    id_non_u = c(rep(0,n_non))    ## working id_non with consecutive, possibly repeating, intergers
    for (i_ in 1:n_non) {
      if (i_ == 1) { count = 1
      } else if (id_non[i_] != last) { count = count + 1 }
      id_non_u[i_] = count
      last = id_non[i_]
    }
  } else {
    n_val_id = n_val
    n_non_id = n_non
    id_val_u = c(1:n_val)
    id_non_u = c(1:n_non)
  }
  n_ful_id = n_val_id + n_non_id
  id_non_ = id_non_u + n_val_id
  id_ful_u = c(id_val_u, id_non_ )

  # fit models for reference, surrogates, validaiton and non-validaiton datasets
  # model REFERENCE measurements from validation (set only)
  if (familyr == "Cox") {
    ref_val = survival::coxph(survival::Surv(y_val , e_val ) ~ x_val , weights=weights_val )
    sur_val = survival::coxph(survival::Surv(ys_val, es_val) ~ xs_val, weights=weights_val )
    if (compare==1) { sur_ful = survival::coxph(survival::Surv(ys_ful, es_ful) ~ xs_ful, weights=weights_ful )
    }   else        { sur_non = survival::coxph(survival::Surv(ys_non, es_non) ~ xs_non, weights=weights_non ) }
  } else {
    ref_val = glm(y_val  ~ x_val , weights=weights_val, family=familyr)         ######## your model statement here #######
    sur_val = glm(ys_val ~ xs_val, weights=weights_val, family=familyr)         ###### and in other calls like this ######
    if (compare==1) { sur_ful = glm(ys_ful ~ xs_ful, weights=weights_ful, family=familyr)
    } else          { sur_non = glm(ys_non ~ xs_non, weights=weights_non, family=familyr) }
  }

  # stats for reference from validation dataset
  beta_val      = ref_val$coefficients
  vcov_beta_val_naive = vcov(ref_val)
  var_beta_val_naive  = diag(vcov_beta_val_naive)

  # stats for surrogates from validation dataset only
  gamma_val     = sur_val$coefficients
  vcov_gamma_val_naive = vcov( sur_val )
  var_gamma_val_naive  = diag(vcov_gamma_val_naive)

  # stats for surrogates from ful or non dataset
  if (compare==1) {
    gamma_ful     = sur_ful$coefficients
    vcov_gamma_ful_naive = vcov( sur_ful )
    var_gamma_ful_naive = diag( vcov_gamma_ful_naive )
    gamma_big = gamma_ful
  } else  {
    gamma_non     = sur_non$coefficients
    vcov_gamma_non_naive = vcov( sur_non )
    var_gamma_non_naive = diag( vcov_gamma_non_naive )
    gamma_big = gamma_non
  }

  #====== estimates for sigma, omega, kappa and vcov(gamma_ful)  =============================================#
  #------ i.e. vcov  beta_val, (gamma_val - gamma_ful) and gamma_ful -----------------------------------------#

  if  ( (familyr %in% c("binomial", "gaussian", "Cox"))  & (vmethod == 1) )  {
    #====== using fdbeta (beta influence) ====================================================================#
    #    print("Entering loop vmethod=1")
    if (familyr == "Cox") {
      dfbeta1 = residuals(ref_val,  type="dfbeta")
      dfbeta2 = residuals(sur_val,  type="dfbeta")
      if (compare==1) { dfbeta3 = residuals(sur_ful,  type="dfbeta") } else {dfbeta3 = residuals(sur_non,  type="dfbeta")  }
    } else {
      dfbeta1 = dfbeta(ref_val)
      dfbeta2 = dfbeta(sur_val)
      if (compare==1) { dfbeta3 = dfbeta(sur_ful) } else { dfbeta3 = dfbeta(sur_non) }
    }
    if (!is.null(id_val)) {
      dfbeta1 = dfbetac(id_val_u, dfbeta1)
      dfbeta2 = dfbetac(id_val_u, dfbeta2)
      if (compare==1) {dfbeta3 = dfbetac(id_ful_u, dfbeta3 )              # dim(dfbeta3)
      } else          {dfbeta3 = dfbetac(id_non_u, dfbeta3) }
    }
    if (compare==1) {
      dfbeta3val   = dfbeta3[1:dim(dfbeta2)[1], ]                         # dim(dfbeta2) ; dim(dfbeta3val)
      dfbeta3non   = dfbeta3[(dim(dfbeta2)[1]+1):dim(dfbeta3)[1], ]       # dim(dfbeta3non) ; dim(dfbeta3)
      gammadiffval = dfbeta2 - dfbeta3val                                 # dim(dfbata2) ; dim(dfbeta3val ) # dim( gammadiffval )
      gammadiff    = rbind(gammadiffval, -dfbeta3non)                     #  t(a - t(b))     ; # dim(gammadiff)
      omega  = t(dfbeta1)   %*% gammadiffval - as.matrix(colMeans(dfbeta1  )) %*% colMeans(gammadiffval)
      kappa  = t(gammadiff) %*% gammadiff    - as.matrix(colMeans(gammadiff)) %*% colMeans(gammadiff   )
      kappa_val = t(gammadiffval) %*% gammadiffval - as.matrix(colMeans(gammadiffval)) %*% colMeans(gammadiffval)
    } else {
      omega  = t(dfbeta1) %*% dfbeta2 - as.matrix(colMeans(dfbeta1)) %*% colMeans(dfbeta2)
      kappa_val = t(dfbeta2) %*% dfbeta2 - as.matrix(colMeans(dfbeta2)) %*% colMeans(dfbeta2)
      kappa_non = t(dfbeta3) %*% dfbeta3 - as.matrix(colMeans(dfbeta3)) %*% colMeans(dfbeta3)
      kappa  = kappa_non + kappa_val
    }
    sigma     = t(dfbeta1) %*% dfbeta1 - as.matrix(colMeans(dfbeta1)) %*% colMeans(dfbeta1)
    sigma_gam = t(dfbeta3) %*% dfbeta3 - as.matrix(colMeans(dfbeta3)) %*% colMeans(dfbeta3)
    sdb = sqrt(diag(sigma))  ;  sdd = sqrt(diag(kappa_val)) ;  omega_cor = t(t( omega / sdb ) / sdd )
    if (compare==1) { vcov_gamma_ful = sigma_gam } else { vcov_gamma_non = sigma_gam }

  } else if ( (familyr %in% c("binomial", "gaussian")) & (vmethod==2)) {
    #======= using formula to estimate beta influence =============================================================#
    #    print("Entering loop vmethod=2")
    x_val1  = cbind(1,x_val)
    xs_val1 = cbind(1,xs_val)
#    xs_non1 = cbind(1,xs_non)
    z_val  = x_val1  %*% beta_val
    zs_val = xs_val1 %*% gamma_val
#    zs_non = xs_non1 %*% gamma_non
    if (compare==1) {
      xs_ful1 = cbind(1,xs_ful)
      zs_ful = xs_ful1 %*% gamma_ful
    }  else        {
      xs_non1 = cbind(1,xs_non)
      zs_non = xs_non1 %*% gamma_non
    }

    sigma  = matrix(0, nrow=dim_mod , ncol=dim_mod )
    omega  = matrix(0, nrow=dim_mod , ncol=dim_mods)
    kappa  = matrix(0, nrow=dim_mods, ncol=dim_mods)
    kappa_val = matrix(0, nrow=dim_mods, ncol=dim_mods)
    kappa_non = matrix(0, nrow=dim_mods, ncol=dim_mods)
    sigma_gam = matrix(0, nrow=dim_mods, ncol=dim_mods)
    dfbeta1 = matrix(0,nrow=dim(x_val)[1] ,ncol=dim(x_val1)[2])
    dfbeta2 = matrix(0,nrow=dim(xs_val)[1],ncol=dim(xs_val1)[2])         # dim(dfbeta2)
    dfbeta3val = matrix(0,nrow=dim(xs_val)[1],ncol=dim(xs_val1)[2])      # dim(dfbeta1m)
    dfbeta3non = matrix(0,nrow=n_non, ncol=dim(xs_val1)[2])

    if (familyr == "binomial") {
      for (i in 1:n_val) {
        dfbeta1[i,] = matrix(vcov_beta_val_naive %*% as.vector((y_val[i] - (1/(1+exp(-z_val[i] )))) %*% as.numeric(x_val1[i,])), nrow=1)
        dfbeta2[i,] =        vcov_gamma_val_naive %*%  as.vector(  (ys_val[i] - (1/(1+exp(-zs_val[i])))) %*% as.numeric(xs_val1[i,]) )
        if (compare==1) {
          dfbeta3val[i,] = vcov_gamma_ful_naive %*%  as.vector(  (ys_ful[i] - (1/(1+exp(-zs_ful[i])))) %*% as.numeric(xs_ful1[i,]) )
          #          dfbeta2m[i,] = dfbeta2m[i,] - dfbeta4m[i,]
        }
      }
      for (i in 1:n_non) {
        if (compare==1) {
          i_ = n_val + i
          dfbeta3non[i,]  = vcov_gamma_ful_naive %*% as.vector( (ys_ful[i_] - (1/(1+exp(-zs_ful[i_])))) %*% as.numeric(xs_ful1[i_,]) )
#          dfbeta3non[i,]  = vcov_gamma_ful_naive %*%  as.vector( (ys_non[i_] - (1/(1+exp(-zs_non[i_])))) %*% as.numeric(xs_non1[i_,]) )
        } else {
          dfbeta3non[i,]  = vcov_gamma_non_naive %*% as.vector( (ys_non[i ] - (1/(1+exp(-zs_non[i ])))) %*% as.numeric(xs_non1[i ,]) )
        }
      }
    } else if (familyr == "gaussian") {
      vcovref_val = vcov_beta_val_naive  / var(ref_val$residuals)     ## update
      vcovsur_val = vcov_gamma_val_naive / var(sur_val$residuals)
      if (compare==1) { vcovsur_ful = vcov_gamma_ful_naive / var(sur_ful$residuals)
      } else { vcovsur_non = vcov_gamma_non_naive / var(sur_non$residuals) }
      for (i in 1:n_val) {
        dfbeta1[i,] = vcovref_val %*%  as.vector(  (y_val[i]  - z_val[i] ) %*% as.numeric(x_val1[i,] ) )
        dfbeta2[i,] = vcovsur_val %*%  as.vector(  (ys_val[i] - zs_val[i]) %*% as.numeric(xs_val1[i,]) )
        if (compare==1) { dfbeta3val[i,] = vcovsur_ful %*%  as.vector(  (ys_ful[i] - zs_ful[i]) %*% as.numeric(xs_ful1[i,]) ) }
      }
      for (i in 1:n_non) {
        if (compare==1) {
#          i_ = n_val + 1
          i_ = n_val + i
          dfbeta3non[i,] = vcovsur_ful %*%  as.vector(  (ys_ful[i_] - zs_ful[i_]) %*% as.numeric(xs_ful1[i_,]) )
        } else {
          dfbeta3non[i,] = vcovsur_non %*%  as.vector(  (ys_non[i ] - zs_non[i ]) %*% as.numeric(xs_non1[i ,]) )
        }
      }
    }

    if (compare==1) { dfbeta3 = rbind(dfbeta3val, dfbeta3non)
    } else { dfbeta3 = dfbeta3non }

    if (!is.null(id_val)) {
      dfbeta1 = dfbetac(id_val_u, dfbeta1)
      dfbeta2 = dfbetac(id_val_u, dfbeta2)
      dfbeta3val = dfbetac(id_val_u, dfbeta3val)
      dfbeta3non = dfbetac(id_non_u, dfbeta3non)
      if (compare==1) { dfbeta3 = dfbetac(id_ful_u, dfbeta3)
      } else { dfbeta3 = dfbetac(id_non_u, dfbeta3)  }
    }

    sigma  = t(dfbeta1) %*% dfbeta1  -  as.matrix(colMeans(dfbeta1)) %*% colMeans(dfbeta1)
    if (compare==1) {
      gammadiffval = dfbeta2 - dfbeta3val
      gammadiff    = rbind(gammadiffval, -dfbeta3non)
      omega  = t(dfbeta1)   %*% gammadiffval - as.matrix(colMeans(dfbeta1  )) %*% colMeans(gammadiffval)
      kappa     = t(gammadiff   ) %*% gammadiff    - as.matrix(colMeans(gammadiff   )) %*% colMeans(gammadiff   )
      kappa_val = t(gammadiffval) %*% gammadiffval - as.matrix(colMeans(gammadiffval)) %*% colMeans(gammadiffval)
      vcov_gamma_ful = t(dfbeta3) %*% dfbeta3 - as.matrix(colMeans(dfbeta3)) %*% colMeans(dfbeta3)
    }  else {
      omega  = t(dfbeta1) %*% dfbeta2 - as.matrix(colMeans(dfbeta1)) %*% colMeans(dfbeta2)
      kappa_val = t(dfbeta2) %*% dfbeta2 - as.matrix(colMeans(dfbeta2)) %*% colMeans(dfbeta2)
      kappa_non = t(dfbeta3non) %*% dfbeta3non - as.matrix(colMeans(dfbeta3non)) %*% colMeans(dfbeta3non)
      kappa  = kappa_non + kappa_val
      vcov_gamma_non = kappa_non
    }
    sigma_gam = t(dfbeta3) %*% dfbeta3 - as.matrix(colMeans(dfbeta3)) %*% colMeans(dfbeta3)
    sdb = sqrt(diag(sigma))  ;  sdd = sqrt(diag(kappa_val)) ;  omega_cor = t(t( omega / sdb ) / sdd )
  } else {
    #====== using JK -- slowest but most generally applicable to different models =================================#
    if (familyr== "Cox") {coxcontrol = survival::coxph.control() }
    glmcontrol = glm.control()

    #------ Jackknife val data ------------------------------------------------------------------------------------#

    for (i_ in 1:n_val_id) {
      index_LO = c(1:n_val)[ id_val_u == i_  ]      # index to leave out ## CHECK and repalce code below if faster
      x_val_j  = x_val [ -index_LO ,]         # x_val_j[1:8,]     # validation set of x, leave out obs i dim(x_val_j)
      y_val_j  = y_val [ -index_LO ]          # validation set of y,
      xs_val_j = xs_val[ -index_LO ,]         # validation set of xs
      ys_val_j = ys_val[ -index_LO ]          # validation set of ys
      xs_ful_j = rbind(xs_val_j, xs_non)
      ys_ful_j = c(ys_val_j, ys_non )
      if( !is.null(weights_val) ) { weights_val_j = weights_val[ -index_LO ]
      } else { weights_val_j = NULL }
      if( !is.null(weights_non) ) { weights_ful_j = c( weights_val_j , weights_non )
      } else { weights_ful_j = NULL }

      if (familyr == "Cox") {
        e_val_j   = e_val [ -index_LO ]
        es_val_j  = es_val[ -index_LO ]
        es_ful_j  = c( es_val_j, es_non )      # length(es_ful_j)
        ref_val_j = survival::coxph.fit(x_val_j ,cbind( y_val_j,  e_val_j) , weights=weights_val_j, strata=NULL, init=beta_val ,
                              control=coxcontrol, method="efron", rownames=NULL)
        sur_val_j = survival::coxph.fit(xs_val_j,cbind(ys_val_j, es_val_j) , weights=weights_val_j, strata=NULL, init=gamma_val,
                              control=coxcontrol, method="efron", rownames=NULL)
        if (compare==1) { sur_ful_j = survival::coxph.fit(xs_ful_j,cbind(ys_ful_j, es_ful_j) , weights=weights_ful_j, strata=NULL,
                                                init=gamma_ful, control=coxcontrol, method="efron", rownames=NULL) }
      } else if ( (familyr == "gaussian") & (vmethod == 0) ) {
        x_val_j1   = cbind(1, x_val_j )
        xs_val_j1  = cbind(1, xs_val_j)
        xs_ful_j1  = cbind(1, xs_ful_j)
        if (is.null(weights_val)) {
          ref_val_j = list( coefficients = solve(t(x_val_j1 ) %*% x_val_j1 ) %*% t(x_val_j1 ) %*% y_val_j  )
          sur_val_j = list( coefficients = solve(t(xs_val_j1) %*% xs_val_j1) %*% t(xs_val_j1) %*% ys_val_j )
          if (compare==1) { sur_ful_j = list( coefficients = solve(t(xs_ful_j1) %*% xs_ful_j1) %*% t(xs_ful_j1) %*% ys_ful_j ) }
        } else {
          ref_val_j = glm.fit(x_val_j1 , y_val_j , weights=weights_val_j, family = gaussian(), control = glmcontrol)  ## why not family=familyr ??
          sur_val_j = glm.fit(xs_val_j1, ys_val_j, weights=weights_val_j, family = gaussian(), control = glmcontrol)
          if (compare==1) { sur_ful_j = glm.fit(xs_ful_j1, ys_ful_j,weights=weights_ful_j,family=gaussian(),control=glmcontrol) }
        }
      } else if ( familyr == "binomial" ) {
        x_val_j1   = cbind(1, x_val_j )
        xs_val_j1  = cbind(1, xs_val_j)
        xs_ful_j1  = cbind(1, xs_ful_j)
        ref_val_j = glm.fit(x_val_j1 , y_val_j , start = beta_val , weights=weights_val_j, family = binomial(link = "logit"),
                            control = glmcontrol)
        sur_val_j = glm.fit(xs_val_j1, ys_val_j, start = gamma_val, weights=weights_val_j, family = binomial(link = "logit"),
                            control = glmcontrol)
        if (compare==1) { sur_ful_j = glm.fit(xs_ful_j1, ys_ful_j, weights=weights_ful_j, start = gamma_ful,
                                               family = binomial(link = "logit"), control = glmcontrol) }
      }  else {
        ref_val_j = glm(y_val_j  ~ x_val_j , weights=weights_val_j, start = beta_val , family=familyr)
        sur_val_j = glm(ys_val_j ~ xs_val_j, weights=weights_val_j, start = gamma_val, family=familyr)
        if (compare==1) { sur_ful_j = glm(ys_ful_j ~ xs_ful_j, weights=weights_ful_j, start = gamma_ful, family=familyr) }
      }

      beta_val_j  = ref_val_j$coefficients
      gamma_val_j = sur_val_j$coefficients
      if (compare==1) { gamma_ful_j = sur_ful_j$coefficients }

      if (compare==1) {
        if (i_==1) { jklist =                cbind( t(beta_val_j), t( gamma_val_j - gamma_ful_j ) , t(gamma_ful_j) )   }
        if (i_!=1) { jklist = rbind( jklist, cbind( t(beta_val_j), t( gamma_val_j - gamma_ful_j ) , t(gamma_ful_j) ) ) }
      } else {
        if (i_==1) { jklist =                cbind( t(beta_val_j), t( gamma_val_j - gamma_non   ) , t(gamma_non  ) )   }
        if (i_!=1) { jklist = rbind( jklist, cbind( t(beta_val_j), t( gamma_val_j - gamma_non   ) , t(gamma_non  ) ) ) }
      }
    }

    #------ Jackknife NON val data --------------------------------------------------------------------------------#

    beta_val_jk = colMeans(jklist[,1:dim_mod])                               ## get mean for jk betas in val data

    # set default jksize so same num of iterations for non and val jk estimates
    if (jksize==0) { jksize = ceiling(n_non_id/n_val_id) }                      ## default leave out size
    jkgrpsn = ceiling( n_non_id / jksize )                                      ## number of leave out groups
    if (jksize==1)  { index_LO = c(rep(1:jkgrpsn)) }                            ## index Leave Out list
    if (jksize!=1)  { index_LO = sample( rep( 1:jkgrpsn , jksize )[ 1:n_non_id ] , n_non_id ) } ;  ##  Leave Out groups ;
    # table(index_LO) ; lo_list[1:50]

    #---------------------------------------------------------

    for (i_ in 1:jkgrpsn) {

      if ( !is.null(id_non) ) {
        index_LO_i = c(1:n_non_id)[ (index_LO==i_) ]               # id_non_u Leave Outs for iteration i_
        drop    = c(1:n_non)[id_non_u %in% index_LO_i]             # drop list
        ys_non_j = ys_non[-drop]
        xs_non_j = xs_non[-drop,]
        if (familyr=="Cox") { es_non_j = es_non[-drop] }
        if( !is.null(weights_non) ) { weights_non_j = weights_non[-drop] }
      } else {                                ## CHECK to define id_non_u when is.null(id_non)==TRUE
        ys_non_j = ys_non[ (index_LO != i_) ]     ; length(ys_non_j)
        xs_non_j = xs_non[ (index_LO != i_), ]    ; dim(xs_non_j)
        if (familyr == "Cox") { es_non_j = es_non[(index_LO != i_)] }
        if( !is.null(weights_non) ) { weights_non_j = weights_non[ (index_LO != i_) ] }
      }

      if (compare==1) {                        # Compare Bval with Gval-Gful
        xs_big_j = rbind(xs_val, xs_non_j)
        ys_big_j = c(as.vector(ys_val), as.vector(ys_non_j) )
        if (familyr == "Cox") { es_big_j = c(as.vector(es_val), as.vector(es_non_j) ) }
        if( !is.null(weights_non) ) { weights_big_j = c( as.vector(weights_val), as.vector(weights_non_j) )
        } else { weights_big_j = NULL }
      } else {                                # Compare Bval with Gval-Gnon
        xs_big_j = xs_non_j
        ys_big_j = ys_non_j
        if (familyr == "Cox") { es_big_j = es_non_j }
        if( !is.null(weights_non) ) { weights_big_j = as.vector(weights_non_j)
        } else { weights_big_j = NULL }
      }
      ## dim(xs_val) ; dim(xs_non_j) ; dim(xs_big_j)

      if (familyr == "Cox") {
        if (vmethod == 0)  {sur_big_j = survival::coxph.fit(xs_big_j,cbind(ys_big_j, es_big_j) , weights=weights_big_j, strata=NULL,
                                               init=gamma_big, control=coxcontrol, method="efron", rownames=NULL)
        } else { sur_big_j = survival::coxph(survival::Surv(ys_big_j, es_big_j) ~ xs_big_j , weights=weights_big_j, init=gamma_big )
        }
      } else if ((familyr == "gaussian") & (vmethod==0)) {
        xs_big_j1 = cbind(1, xs_big_j)    ; dim(xs_big_j) ; xs_big_j1[1:25,]
        if (is.null(weights_non)) {
          sur_big_j = list( coefficients = solve(t(xs_big_j1) %*% xs_big_j1) %*% t(xs_big_j1) %*% ys_big_j )
        } else {
          glm.fit(xs_big_j1, ys_big_j, weights=weights_big_j, start = gamma_big, family = familyr, control = glmcontrol)
        }
      } else if ((familyr == "binomial") & (vmethod==0)) {
        xs_big_j1 = cbind(1, xs_big_j)
        sur_big_j = glm.fit(xs_big_j1, ys_big_j, weights=weights_big_j, start = gamma_big, family = binomial(link = "logit"),
                            control = glmcontrol)
      } else { sur_big_j = glm(ys_big_j ~ xs_big_j, weights=weights_big_j, start = gamma_big, family=familyr)
      }

      gamma_big_j = sur_big_j$coefficients

      jklist = rbind( jklist, cbind( t(beta_val_jk), t( gamma_val - gamma_big_j ) , t(gamma_big_j) ) )
    }

    #---------- final calculations for jackknife estimates for for sigma, omega, kappa and vcov(gamma_ful/non) ------------#

    jklistrows= dim(jklist)[1]
    jklistval = jklist[1:n_val_id, 1:dim_modt]
    sigma     = (n_val_id - 1) * cov(jklistval[1:n_val_id, 1:dim_mod])
    omega     = (n_val_id - 1) * cov(jklistval)[1:dim_mod, (dim_mod+1):dim_modt]     ## cov between beta_val and (gamma_val - gamma_ful/non)
    omega_cor =                  cor(jklistval)[1:dim_mod, (dim_mod+1):dim_modt]     ## cov between beta_val and (gamma_val - gamma_ful/non)
    if (compare==1) {
      jklistdiff = jklist[, (dim_mod +1):dim_modt ]
      jklistgam  = jklist[, (dim_modt+1):(dim_modt+dim_mods)]
      kappa  = (jklistrows - 1) * cov( jklistdiff )
      sigma_gam = (jklistrows - 1) * cov( jklistgam )
      vcov_gamma_ful = sigma_gam
    } else {
      jklistdiffval = jklist[1:n_val              , (dim_mod+1):dim_modt ]
      jklistdiffnon = jklist[(n_val+1):jklistrows , (dim_mod+1):dim_modt ]
      jklistgam     = jklist[(n_val+1):jklistrows , (dim_modt+1):(dim_modt+dim_mods)]
      kappa = (n_val - 1) * cov( jklistdiffval ) + (jklistrows - n_val - 1) * cov( jklistdiffnon )
      sigma_gam = (jklistrows - n_val - 1) * cov(jklistgam)
      vcov_gamma_non = sigma_gam
    }
  }

  #-----------  End of JK loop -----------

  ##### get augmented estimates ##################################################################################

  gamma_dif = (gamma_val - gamma_big)

  beta_aug       = beta_val - t( omega  %*% solve( kappa ) %*% gamma_dif )
  vcov_beta_aug  = sigma - omega %*% solve( kappa ) %*% t(omega)

  #====== combined estimates ===============================================================#
  beta_aug = as.vector(beta_aug) ;         ## beta_aug_a = as.vector(beta_aug_a) ; beta_aug_b = as.vector(beta_aug_b)
  coef_beta  = rbind( beta_aug  , beta_val  )         ## depending on models gamma may be longer or shorter than beta
  if (compare==1) { coef_gamma = rbind( gamma_ful , gamma_val ) } else { coef_gamma = rbind( gamma_non , gamma_val ) }

  #====== combine variances ================================================================#
  var_beta_aug   = diag( vcov_beta_aug  )
  var_beta_val  = diag(sigma)
  var  = rbind(var_beta_aug, var_beta_val=var_beta_val, var_beta_val_naive)  #  is.vector(var_gamma_big) ; is.vector(var_gamma_big_naive)
  if (compare==1) {
    var_gamma_ful = diag(vcov_gamma_ful)
    varg = rbind( var_gamma_ful, var_gamma_ful_naive, var_gamma_val_naive)
  } else {
    var_gamma_non = diag(vcov_gamma_non)
    varg = rbind( var_gamma_non, var_gamma_non_naive, var_gamma_val_naive)
  }

  #=============  RETURN Information  ==============================================================================

  Call <- match.call()
  indx <- match(c("x_val","y_val","xs_val","ys_val","xs_non","ys_non","jksize","vmethod","familyr",
                  "e_val","es_val","es_non","id_val","id_non"),
                  names(Call), nomatch=0)

  if (compare == 1) { comparec = "ful" } else { comparec = "non" }
  if (vmethod <= 0) { vmethodc = "jk"} else if (vmethod==1) {vmethodc="ijk (dfbeta)"} else {vmethodc="ijk (alt)"}
  if (!is.null(id_val)) { ns = cbind(n_val=n_val, n_ful=n_ful, n_val_id, n_ful_id)
  } else { ns = cbind(n_val=n_val, n_ful=n_ful) }

  if (vmethod <= 0) { inputs = cbind( familyr=familyr, compare=compare, comparec=comparec, vmethod=vmethod, vmethodc=vmethodc,
                                   jksize=jksize, ns, dim_beta=dim_mod, dim_gamma=dim_mods )
  }  else      { inputs = cbind( familyr=familyr, compare=compare, comparec=comparec, vmethod=vmethod, vmethodc=vmethodc,
                                 ns, dim_beta=dim_mod, dim_gamma=dim_mods )
  }

  names_x  = colnames(x_val)  ; if (familyr != "Cox") { names_x =c("(Intercept)", names_x ) }
  names_xs = colnames(xs_val) ; if (familyr != "Cox") { names_xs=c("(Intercept)", names_xs) }
  names_y  = colnames(y_val)  ;
  names_ys = colnames(ys_val) ;

  if (compare==1) {
    returnlist = list(Call=Call, FitInput=inputs, names_x=names_x, names_xs=names_xs, names_y=names_y, names_ys=names_ys,
                      coef_beta=coef_beta, coef_gamma=coef_gamma, var_beta=var, var_gamma=varg,
                      vcov_beta_aug=vcov_beta_aug, vcov_beta_val=sigma, vcov_beta_val_naive=vcov_beta_val_naive,
                      vcov_gamma_ful=vcov_gamma_ful, vcov_gamma_ful_naive=vcov_gamma_ful_naive,
                      omega=omega, omega_cor=omega_cor, kappa=kappa )
  } else {
    returnlist = list(Call=Call, FitInput=inputs, names_x=names_x, names_xs=names_xs, names_y=names_y, names_ys=names_ys,
                      coef_beta=coef_beta, coef_gamma=coef_gamma, var_beta=var, var_gamma=varg,
                      vcov_beta_aug=vcov_beta_aug, vcov_beta_val=sigma, vcov_beta_val_naive=vcov_beta_val_naive,
                      vcov_gamma_non=vcov_gamma_non, vcov_gamma_non_naive=vcov_gamma_non_naive,
                      omega=omega, omega_cor=omega_cor, kappa=kappa )
  }
  class(returnlist) <- c("meerva")
  return(returnlist)
}

################  END MAIN PROGRAM       #####################################################################################

################  BEGIN ANALYSIS SUMMARY TOOL     ############################################################################

#' Summarize Information for a meerva Output Object
#'
#' @param object A meerva class object for summary.
#' @param ... further arguments 
#'
#' @return Summarize output
#' 
#' @export
#' 
#' @importFrom stats qnorm pnorm 
#'
summary.meerva = function(object, ...) {
  compare= object$FitInput[2]

  cat(paste0("\n"))
  print(object$Call) ; cat(paste0("\n"))
  print(object$FitInput) ; cat(paste0("\n"))
  #  cat(paste0(object$FitInput)) ; cat(paste0("\n"))

  ztest = function(estimate, var, names, alpha=0.05) { # estimate = object$coef_beta[1,] ; var=object$var_beta[1,] ; names = object$names_x ;
    estimate = as.vector(estimate)
    se = as.vector(sqrt(var))
    lcl = estimate + qnorm(alpha/2) * se
    ucl = estimate - qnorm(alpha/2) * se
    z = estimate/se
    p = 2*pnorm(-abs(z))
    summary = cbind(estimate, se, lcl, ucl, z, p)
    rownames(summary) = names
    print(summary)
    cat(paste0("\n"))
  }

  cat(paste0(" Estimates for Beta using beta_aug (references augmented with surrogates)\n"))
  ztest(object$coef_beta[1,], object$var_beta[1,], object$names_x )

  cat(paste0(" Estimates for Beta using beta_val (references alone\n"))
  ztest( object$coef_beta[2,], object$var_beta[2,], object$names_x )

  cat(paste0(" Effective multiplicative increase in sample size by using augmented estimates\n"))
  Increase = object$var_beta[2,]/ object$var_beta[1,]
  print(Increase)
  cat("\n")

  if (compare==1) { cat(paste0(" Estimates for Gamma using gamma_ful (surrogates alone)\n not for direct comparisons \n"))
  }   else        { cat(paste0(" Estimates for Gamma using gamma_non (surrogates alone)\n not for direct comparisons \n"))
  }
  ztest( object$coef_gamma[1,], object$var_gamma[1,], object$names_xs )

  if (compare==1) { cat(paste0(" Correlations between beta_val and (gamma_val - gamma_ful)  \n"))
  }   else        { cat(paste0(" Correlations between beta_val and gamma_val \n"))
  }
  print( object$omega_cor )
  cat(paste0("\n"))
}

################  END analysis summary tool     ##############################################################################

################  BEGIN ANALYSIS PRINT TOOL     ############################################################################

#' Print Minimal Summary Information for a meerva Output Object
#'
#' @param x A meerva class object for printing 
#' @param ... further arguments 
#'
#' @return Print output
#' 
#' @export
#' 
#' @importFrom stats qnorm pnorm 
#'
print.meerva = function(x, ...) {
  compare= x$FitInput[2]

  cat(paste0("\n"))
  print(x$Call) ; cat(paste0("\n"))
  print(x$FitInput) ; cat(paste0("\n"))
  #  cat(paste0(x$FitInput)) ; cat(paste0("\n"))

  ztest = function(estimate, var, names, alpha=0.05) { # estimate = x$coef_beta[1,] ; var=x$var_beta[1,] ; names = x$names_x ;
    estimate = as.vector(estimate)
    se = as.vector(sqrt(var))
    lcl = estimate + qnorm(alpha/2) * se
    ucl = estimate - qnorm(alpha/2) * se
    z = estimate/se
    p = 2*pnorm(-abs(z))
    summary = cbind(estimate, se, lcl, ucl, z, p)
    rownames(summary) = names
    print(summary)
    cat(paste0("\n"))
  }
  cat(paste0(" Estimates for Beta using beta_aug \n"))
  ztest(x$coef_beta[1,], x$var_beta[1,], x$names_x )
}

##############################################################################################################################

#' Sum dfbeta s According to id_vector Clusters
#'
#' @param id_vector Input id vector
#' @param dfbeta    dfbeta s for sandwich
#'
#' @return dfbeta by id_vector clusters
#'
#'
dfbetac = function(id_vector, dfbeta) {
  n_indeces = sum( table(table(id_vector)) )
  n_ = length(id_vector)

  dfbetac = matrix(0, nrow=n_indeces, ncol=dim(dfbeta)[2])
  for (i_ in 1:n_) {
    if (i_ == 1) {
      mdfbeta = dfbeta[i_,]
      last = id_vector[i_]
    } else {
      if (id_vector[i_] == last) { mdfbeta = rbind(mdfbeta, dfbeta[i_,])
      } else {
        if (is.vector(mdfbeta)) { dfbetac[id_vector[i_-1],] = mdfbeta   # dfbetac[1:5,]
        } else { dfbetac[id_vector[i_-1],] = colSums(mdfbeta) }
        mdfbeta = dfbeta[i_,]
        last = id_vector[i_]
      }
    }
    if (i_ == n_) {
      if (is.vector(mdfbeta)) { dfbetac[id_vector[i_],] = mdfbeta
      } else { dfbetac[id_vector[i_],] = colSums(mdfbeta)
      }
    }
  }
  rownames(dfbetac) = NULL  ;
  return(dfbetac)
}

#===================================================================================

#' Analysis of Data with Measurement Error Using a Validation Subsample
#'
#' @description
#' The meerva package performs regression analyses on data
#' with measurement error when there is a
#' validation subsample.  The functional .fit program is meerva.fit. The meerva function
#' is intended for future development and use as a wrapper for meerva.fit.
#' Please help(meerva.fit).
#'
#' @author Walter Kremers (kremers.walter@mayo.edu)
#'
#' @return Describes future development of the meerva package. 
#'
#' @export
#'
#' @seealso
#'  meerva.fit , meerva.sim.block , meerva.sim.nrm , meerva.sim.brn , meerva.sim.cox
#'
meerva = function() {
cat("meerva is intended for future development and use as a wrapper \n") ;
cat("for meerva.fit, which analyzes regression data with measurement \n")
cat("error when there is a validation subsample. \n\n")
cat("Try help(meerva.fit) \n\n")
}

