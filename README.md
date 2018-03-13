# maxadjAUCpen
R code to maximize adjusted AUC with optional penalization

This code is related to methods described in <http://biostats.bepress.com/uwbiostat/paper421/>. 

See below for documentation of this function.

This code is related to that in the R package "maxaAUC", but the code provided here allows for penalization. 

---------------------------------------------------------

__Imports:__ survival, aucm, Rsolnp, Hmisc, foreach

__Description__

Often there is interest in combining several predictors or biomarkers into a linear combination for diagnosis, prognosis or screening. This can be done by targeting measures of predictive capacity. In the presence of a discrete covariate, such as batch or study center, an appropriate summary of discriminatory performance is the covariate-adjusted area under the receiver operating characteristic curve (AUC), or aAUC. This function implements the SaAUC method: it estimates a linear combination of predictors by maximizing a smooth approximation to the aAUC. In addition, when the discrete covariate is a nuisance variable, there may be interest in fitting a combination with good performance overall but lower variability in covariate-specific performance. This function allows the user to penalize variability in the covariate-specific AUCs when fitting the linear combination.

__Usage__

`maxadjAUC(outcome, predictors, covariate, lambda=0, lambdaseq=10^seq(from=log10(0.1), to=log10(200), length=20), initialval="rGLM", approxh = 1/3, conditional=FALSE, tolval = 1e-6, stepsz = 1e-5, warnfileTR="warnTR.txt", warnfileCV="warnCV.txt")`

__Arguments__

_outcome_ A vector of outcome (disease) indicators for each observation (1 for diseased, 0 for non-diseased). Missing values are not allowed. predictors A numeric matrix with one row for each observation and one column for each candidate predictor. Missing values are not allowed. The columns of the matrix will be (re)named "V1", "V2", ....

_covariate_ A numeric vector of covariate values for each observation. The covariate should have a limited number of values (i.e., it should be a discrete covariate). Missing values are not allowed.

_lambda_ The value of the penalty parameter (see 'Details'). Must be non-negative; default value is 0, i.e., no penalization.

_lambdaseq_ A vector of penalty parameter values to be evaluated. Default is 10^seq(from=log10(0.1), to=log10(200), length=20). If there is interest in only one value of the penalty parameter, use lambda to identify that value, and set lambdaseq=NULL.

_initialval_ Starting values of the predictor combination for the SaAUC algorithm. Default value is "rGLM", which means that estimates from robust logistic regression, specifically the method of Bianco and Yohai (implemented via the aucm package), are used as starting values. If any other value of initialval is given, or if robust logistic regression fails to converge, estimates from standard logistic regression are used as starting values. For both robust and standard logistic regression, the covariate will be included as a stratifying variable.

_approxh_ The tuning parameter for the smooth approximation to the covariate-specific AUC is the ratio of the standard deviation of the linear combination (based on the starting values) to n_c^approxh, where n_c is the number of observations with covariate value c. In particular, larger values of approxh will provide a better approximation to the AUC, though estimation may become unstable if approxh is too large. Default 1/3.

_conditional_ A logical value indicating whether standard logistic regression should be conditional if TRUE (i.e., survival::clogit) or unconditional if FALSE (stats::glm). Note that if the number of covariate strata is large (the situation to which conditional logistic regression is commonly applied) and lambdaseq is not NULL, the cross-validation procedure described below will be time-consuming. Default is FALSE.

_tolval_ Controls the tolerance on feasibility and optimality for the optimization procedure (performed by solnp in the Rsolnp package). Default 1e-6.

_stepsz_ Controls the step size for the optimization procedure (performed by solnp in the Rsolnp package). Default 1e-5.

_warnfileTR_ The name (including path) of the file where warnings from the analysis in the training data for lambdaseq will be sent. This includes warnings resulting from (i) the application of the SaAUC method for each value in lambdaseq and (ii) estimating the covariate-specific AUCs in the training data for the SaAUC fitted combination. That is, if lambdaseq is not NULL, the SaAUC method will be applied for each value of the penalty parameter in the vector lambdaseq and the covariate-specific AUCs for the resulting combination will be calculated in the training data; any warnings produced will be sent to warnfileTR. As a result, any warnings will not be printed in the console. If the function is run multiple times without changing warnfileTR, the file will be appended. Default value is "warnTR.txt", which will be created in the working directory. 

_warnfileCV_ Same as warnfileTR, but this includes any warnings resulting from the crossvalidation procedure (described in ‘Details’). This procedure includes (i) fitting standard and robust logistic regression in each of the cross-validation training datasets, and applying the SaAUC method for each value in lambdaseq to each of the cross-validation training datasets and (ii) estimating the covariate-specific AUCs for the SaAUC combination fitted to each of the cross-validation training datasets; any resulting warnings will be sent to warnfileCV. Default value is "warnCV.txt".

__Details__

The function seeks to optimize a smooth approximation to the (penalized) covariate-adjusted AUC:
SaAUC - \lambda \sum_{c=1}^m w_c(SAUC_c - SaAUC)^2 
where SAUC_c is the smooth approximation to the covariate-specific AUC, w_c are covariate-specific weights, \lambda is the penalty parameter, and SaAUC is the smooth approximation to the covariate-adjusted AUC (SaAUC = \sum_{c=1}^m w_c SAUC_c), for a covariate with m values in the data.

If a vector of penalty parameter values is provided via lambdaseq, the function output includes two plots: one contains the results for each penalty parameter value in the training data, and the other presents the results for a cross-validation procedure where for each value of the penalty parameter, observations in all but one of the m covariate strata are used as training data and the observations in the remaining covariate stratum are used as test data (termed "leave one covariate out crossvalidation", or LOCOCV).

In both plots, the left y-axis corresponds to the AUC and the right x-axis corresponds to the variability in covariate-specific AUCs around the aAUC. The gray lines correspond to the covariate-specific AUCs for the SaAUC method for each value of the penalty parameter (plotted on the x-axis). In the cross-validation plot these covariate-specific AUCs are estimated in the test stratum. The black solid line on each plot corresponds to the aAUC for the SaAUC method for each value of the penalty parameter. The black dashed and dot-dashed lines correspond to the aAUCs for standard and robust logistic regression. The red solid lines correspond to the variability in the covariate-specific AUCs around (i) the aAUC in the training data, and (ii) the aAUC based on the training strata in LOCOCV. The blue solid line in the LOCOCV plot corresponds to the variability in the covariatespecific AUCs around the aAUC based on the test strata. In the training data plot, the red dashed and dot-dashed lines correspond to the variability in the covariate-specific AUCs around the aAUC in the training data for standard and robust logistic regression.

The function is set up to allow for parallel processing via foreach::foreach. To incorporate parallel processing, users must have (e.g., to set up 4 clusters)

`library(doParallel)`

`cl = makeCluster(4)`

`registerDoParallel(cl)`

`maxadjAUC(...)`

`stopCluster(cl)`

If users do not register clusters, a warning will be issued:

`Warning message: executing %dopar% sequentially: no parallel backend registered`

__Value__

A list will be returned with the following components:

_NumCov_ The number of covariate strata used after removing concordant strata.

_FittedCombs_ A list containing four fitted combinations: InitialVal (either robust logistic regression, if initialval="rGLM", standard unconditional logistic regression, if initialval is not "rGLM" and conditional=FALSE, or standard conditional logistic regression, if initialval is not "rGLM" and conditional=TRUE), NormGLM (standard unconditional or conditional logistic regression, depending on conditional), NormrGLM (robust logistic regression), MaxSaAUCSupplied (SaAUC approach when lambda is used as the penalty parameter). All fitted combination vectors are normalized. If robust logistic regression fails to converge, standard unconditional or conditional logistic regression (depending on conditional) is used instead, and a warning is given.

_aAUCTR_ A vector of the aAUC in the training data for the four fitted combinations.

_varTR_ A vector of the variability in the covariate-specific AUCs around the aAUC in the training data for the four fitted combinations. 

Furthermore, if lambdaseq is not NULL the list will additionally include the following components:

_TRrslts_ A matrix of training data results including a column with the value of the penalty parameter ("lambda"), the aAUC in the training data ("aAUCTR"), the variability in the covariate-specific AUCs around the aAUC in the training data ("varTR"), m columns with the covariate-specific AUCs in the training data, and convergence results for the SaAUC approach, where 0 indicates convergence ("TRconv").

_CVrslts_ A matrix of cross-validation results including a column with the value of the penalty parameter ("lambda"), the cross-validated aAUC ("aAUCCV"), the variability in the covariate-specific AUCs in the cross-validation test strata around the aAUC in the cross-validation training strata ("varCVTR"), the variability in the covariate-specific AUCs in the cross-validation test strata around the aAUC in the cross-validation test strata ("varCVTE"), and m columns with the covariate-specific AUCs in each of the m cross-validation test strata.

_CVconv_ A matrix of convergence results for the SaAUC approach and robust logistic regression for each value of the penalty parameter. A value of 0 indicates convergence. Note that if robust logistic regression fails to converge in any of the m cross-validation "training" datasets, it will fail for all values of lambda since robust logistic regression is not related to lambda.

_plotOut_ A plot of training data results and LOCOCV results, as described above.

__Note__

The function automatically removes any covariate strata that are concordant on the outcome (i.e., all 0 or all 1).

Warnings are issued if the SaAUC algorithm does not converge for lambda in the training data or if the robust logistic regression does not converge in the training data. Furthermore, if a sequence of penalty parameter values is provided via lambdaseq, warnings are issued if the SaAUC algorithm fails to converge for some value(s) of the penalty parameter in the training data or in cross-validation, or if the robust logistic regression does converge for some subset of strata in the cross-validation procedure. Finally, warnings are issued if warnfileTR or warnfileCV exist and are non-empty. If such a warning is given, the relevant file (warnfileTR and/or warnfileCV) should be examined. The standard unconditional or conditional logistic regression algorithm may not converge, producing a warning. If such a convergence failure occurs, the "GLM" results will be affected, as will the "rGLM" results if the robust logistic model also fails to converge, and the "SaAUC" results if initialval is not "rGLM" or if the robust logistic model fails to converge. Thus, users should be alert to any convergence failures.

__References__

Bianco, A.M. and Yohai, V.J. (1996) Robust estimation in the logistic regression model. In Robust statistics, data analysis, and computer intensive methods (ed H. Rieder), pp 17-34. Springer.

Janes, H. and Pepe, M.S. (2009) Adjusting for covariate effects on classification accuracy using the covariate-adjusted receiver operating characteristic curve. Biometrika, pages 1-12.

Meisner, A., Parikh, C.R., and Kerr, K.F. (2017). Developing biomarker combinations in multicenter studies via direct maximization and penalization. UW Biostatistics Working Paper Series, Working Paper 421.

__See Also__

`rlogit, solnp, foreach, %dopar%`

__Examples__

`## takes a few minutes to run`

`expit <- function(x){`

`exp(x)/(1+exp(x))`

`}`

`set.seed(1)`

`covar <- rep(c(1:5),each=200)`

`x1 <- rnorm(1000,0,rep(runif(5,0.8,1.2),each=200))`

`x2 <- rnorm(1000,0,rep(runif(5,0.8,1.2),each=200))`

`x3 <- rnorm(1000,0,rep(runif(5,0.8,1.2),each=200))`

`x4 <- rnorm(1000,0,rep(runif(5,0.8,1.2),each=200))`

`covint <- rep(runif(5,-1.5,1.5), each=200)`

`y <- rbinom(1000,1,expit(covint + 1*x1 - 1*x2 + 1*x3 - 1*x4))`

`X <- cbind(x1,x2,x3,x4)`

`output <- maxadjAUC(outcome=y, predictors=X, covariate=covar, lambda=0, lambdaseq=10^seq(from=log10(0.1), to=log10(200), length=10), initialval="rGLM", approxh = 1/3, conditional=FALSE, tolval = 1e-6, stepsz = 1e-5, warnfileTR="warnTR.txt", warnfileCV="warnCV.txt")`

`plot.new()`

`output[["plotOut"]]`

