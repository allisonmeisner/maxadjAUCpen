##########
### Code to maximize adjusted AUC with optional penalization
##########

############ Define helper functions
globalVariables('ind')

saveAllWarnings <- function(expr, logFile="warning_log.R") {
  withCallingHandlers(expr,
                      warning=function(w) {
                        cat(conditionMessage(w), "\n\n",  file=logFile, append=TRUE)
                        invokeRestart("muffleWarning")
                      })
}

saauc <- function(betavec, Xmat, yvec, hc, lambdapar){
  wts <- prop.table(table(yvec, Xmat[,1]),1)[2,]
  costval <- rep(NA,length(wts))

  betadiffs <- list();
  for(i in 1:length(wts)){
    XDc <- Xmat[(Xmat[,1]==i & yvec==1),]
    XnDc <- Xmat[(Xmat[,1]==i & yvec==0),]
    diffs <- matrix(NA, nrow=nrow(XDc)*nrow(XnDc), ncol=ncol(Xmat)-1)
    for(j in 1:ncol(diffs)){
      diffs[,j] <- as.vector(outer(XDc[,j+1],XnDc[,j+1],"-"))
    }
    diffsB <- diffs %*% betavec
    costval[i] = -sum(pnorm(diffsB/hc[i]))/nrow(diffs)
  }

  totcost = (costval %*% wts) + lambdapar*(wts %*% (costval - (costval %*% wts))^2)
  as.numeric(totcost)
}

"%ni%" <- Negate("%in%")

constr <- function(betavec, Xmat, yvec, hc, lambdapar){
  norm(matrix(betavec), type="F")
}
####################

maxadjAUC <- function(outcome, predictors, covariate, lambda=0,
                      lambdaseq=10^seq(from=log10(0.1), to=log10(200), length=20),
                      initialval="rGLM", approxh = 1/3, conditional=FALSE, tolval = 1e-6, stepsz = 1e-5,
                      warnfileTR="warnTR.txt", warnfileCV="warnCV.txt"){

  ###### Give error if any observations missing
  if(max(is.na(outcome))==1 | max(is.na(predictors))==1 | max(is.na(covariate))==1){
    stop("missing values not allowed")
  }
  ###### Give error if wrong variable types
  if(!is.matrix(predictors)){
    stop("predictors must be a matrix")
  }
  if(!is.numeric(covariate) | !is.numeric(outcome) | !is.numeric(predictors)){
    stop("outcome, covariate, and predictors must be numeric")
  }
  ###### Give error if outcome is not 0/1
  outcomevalue <- unique(outcome)
  if((length(outcomevalue) != 2) | (min(outcomevalue) != 0) | (max(outcomevalue) != 1)){
    stop("outcome must be 0/1 variable")
  }
  ###### Give error if any inputs have different lengths
  if(length(outcome) != length(covariate) | length(outcome) != nrow(predictors)){
    stop("outcome, covariate, and predictors are not the same length")
  }

  ### Re-level covariate so it's 1,...,M
  covval <- as.factor(covariate)
  levels(covval) <- rank(as.numeric(levels(covval)))
  cord <- as.numeric(covval)

  ###### Remove covariate values that are concordant
  numcases <- table(outcome, cord)[2,]
  numobs <- table(cord)
  dropcovar <- which(pmax(numcases,numobs-numcases)==numobs)
  onew <- outcome[which(cord %ni% dropcovar)]
  pnew <- predictors[which(cord %ni% dropcovar), ]
  cnew <- cord[which(cord %ni% dropcovar)]
  numcov <- length(table(cnew))

  ### Re-level covariate so it's 1,...,numcov
  covvalnew <- as.factor(cnew)
  levels(covvalnew) <- rank(as.numeric(levels(covvalnew)))
  cnew <- as.numeric(covvalnew)

  ###### Get initial estimates
  data <- data.frame(onew, cnew, pnew)
  names(data) <- c("outcome","covariate",paste("V",c(1:ncol(pnew)),sep=""))
  varnames <- names(data)
  prednames <- varnames[-c(1:2)]
  ### Logistic regression
  if(conditional==FALSE){
    glmcoef <- glm(as.formula(paste(varnames[1], " ~ ", paste(prednames, collapse=" + "),
                                    " + factor(",varnames[2],")", sep="")),
                   family="binomial", data=data)$coef[2:(length(prednames)+1)]
  }else{
    glmcoef <- survival::clogit(as.formula(paste(varnames[1], " ~ ", paste(prednames, collapse=" + "),
                    " + strata(",varnames[2],")", sep="")), data=data)$coef[1:length(prednames)]
  }
  ### Robust logistic regression
  robustmod <- aucm::rlogit(as.formula(paste(varnames[1], " ~ ", paste(prednames, collapse=" + "),
                                             " + factor(",varnames[2],")", sep="")), dat=data)
  if(robustmod$convergence==TRUE){
    rglmcoef <- robustmod$coef[2:(length(prednames)+1)]
  }else{
    rglmcoef <- glmcoef
  }

  if(initialval=="rGLM"){
    beta0 = rglmcoef/norm(matrix(rglmcoef),type="F")
  }else{
    beta0 = glmcoef/norm(matrix(glmcoef),type="F")
  }
  normglm = glmcoef/norm(matrix(glmcoef),type="F")
  normrglm = rglmcoef/norm(matrix(rglmcoef),type="F")
  hvalc <- rep(NA,numcov)
  for(i in 1:numcov){
    hvalc[i] = sd(as.matrix(pnew[which(cnew==i),]) %*% matrix(beta0,ncol=1))/((table(cnew)[i])^approxh)
  }

  #### Results for supplied value of lambdapar
  rsltsupplied <- Rsolnp::solnp(beta0, saauc, eqfun=constr, eqB=1, Xmat=data[,-1], yvec=data[,1], hc=hvalc,
                        lambdapar=lambda,control=list("outer.iter"=10^3, "inner.iter"=10^4, "delta"=stepsz,
                                                      "tol"=tolval,"trace"=0))

  #### Performance in training data
  AUCcTRrglm <- rep(NA, numcov)
  AUCcTRglm <- rep(NA, numcov)
  AUCcTRsuppl <- rep(NA, numcov)
  wtvals <-  prop.table(table(onew, cnew),1)[2,]

  for(i in 1:numcov){
    AUCcTRrglm[i] = Hmisc::somers2(pnew[which(cnew==i),] %*% normrglm, onew[which(cnew==i)])[1]
    AUCcTRglm[i] = Hmisc::somers2(pnew[which(cnew==i),] %*% normglm, onew[which(cnew==i)])[1]
    AUCcTRsuppl[i] = Hmisc::somers2(pnew[which(cnew==i),] %*% rsltsupplied$pars, onew[which(cnew==i)])[1]
  }
  aAUCTRrglm <- AUCcTRrglm %*% wtvals
  varTRrglm <-(wtvals %*% (AUCcTRrglm - AUCcTRrglm %*% wtvals)^2)

  aAUCTRsuppl <- AUCcTRsuppl %*% wtvals
  varTRsuppl <- (wtvals %*% (AUCcTRsuppl - AUCcTRsuppl %*% wtvals)^2)

  aAUCTRglm <- AUCcTRglm %*% wtvals
  varTRglm <- (wtvals %*% (AUCcTRglm - AUCcTRglm %*% wtvals)^2)

  ##### Looking at multiple lambda values
  if(length(lambdaseq)>0){
    lambdaseqTR <- function(lambdaval){
      estpar <- Rsolnp::solnp(beta0, saauc, eqfun=constr, eqB=1, Xmat=data[,-1], yvec=data[,1], hc=hvalc,
            lambdapar=lambdaval,control=list("outer.iter"=10^3, "inner.iter"=10^4, "delta"=stepsz,
                                             "tol"=tolval,"trace"=0))

      costLambda <- rep(NA, numcov)
      for(i in 1:numcov){
        costLambda[i] = Hmisc::somers2(pnew[which(cnew==i),] %*% estpar$pars, onew[which(cnew==i)])[1]
      }
      c(estpar$convergence, costLambda %*% wtvals, sqrt(wtvals %*% ((costLambda - (costLambda %*% wtvals))^2)),
        costLambda)
    }
    rsltsTR <- foreach(j=1:length(lambdaseq), .inorder=FALSE, .combine="rbind",
                       .packages=c("Hmisc","Rsolnp"),
                       .export=c("saveAllWarnings","saauc","%ni%","constr",
                                 "beta0","data","hvalc","stepsz","tolval","numcov",
                                 "pnew","cnew","onew","wtvals")) %dopar% {
                         cbind(lambdaseq[j],
                               matrix(saveAllWarnings(lambdaseqTR(lambdaval=lambdaseq[j]),logFile=warnfileTR),
                                      ncol=numcov+3))
                       }

    lambdaseqCV <- function(lambdaval){
      convnow=0; rglmfailCV=0
      rsltsCV <- matrix(NA, nrow=numcov, ncol=2)
      for(i in 1:numcov){
        pCV <- pnew[which(cnew!=i),]
        oCV <- onew[which(cnew!=i)]
        cCV <- cnew[which(cnew!=i)]

        ### Re-level covariate so it's 1,...,numcov-1
        covvalCV <- as.factor(cCV)
        levels(covvalCV) <- rank(as.numeric(levels(covvalCV)))
        cCV <- as.numeric(covvalCV)
        ### Create dataset, estimate combinations
        dataCV <- data.frame(oCV, cCV, pCV)
        names(dataCV) <- c("outcome","covariate",paste("V",c(1:ncol(pCV)),sep=""))
        varnamesCV <- names(dataCV)
        prednamesCV <- varnamesCV[-c(1:2)]
        ## Logistic regression
        if(conditional==FALSE){
          glmcoefCV <- glm(as.formula(paste(varnamesCV[1], " ~ ", paste(prednamesCV, collapse=" + "),
                                          " + factor(",varnamesCV[2],")", sep="")),
                         family="binomial", data=dataCV)$coef[2:(length(prednamesCV)+1)]
        }else{
          glmcoefCV <- survival::clogit(as.formula(paste(varnamesCV[1], " ~ ", paste(prednamesCV, collapse=" + "),
                                                       " + strata(",varnamesCV[2],")", sep="")),
                                        data=dataCV)$coef[1:length(prednamesCV)]
        }
        ## Robust logistic regression
        robustmodCV <- aucm::rlogit(as.formula(paste(varnamesCV[1], " ~ ", paste(prednamesCV, collapse=" + "),
                                                   " + factor(",varnamesCV[2],")", sep="")), dat=dataCV)
        if(robustmodCV$convergence==TRUE){
          rglmcoefCV <- robustmodCV$coef[2:(length(prednamesCV)+1)]
        }else{
          rglmcoefCV <- glmcoefCV
          rglmfailCV <- 1
        }

        if(initialval=="rGLM"){
          beta0CV = rglmcoefCV/norm(matrix(rglmcoefCV),type="F")
        }else{
          beta0CV = glmcoefCV/norm(matrix(glmcoefCV),type="F")
        }
        normglmCV = glmcoefCV/norm(matrix(glmcoefCV),type="F")
        normrglmCV = rglmcoefCV/norm(matrix(rglmcoefCV),type="F")
        hvalcCV <- rep(NA,numcov-1)
        for(k in 1:(numcov-1)){
          hvalcCV[k] = sd(as.matrix(pCV[which(cCV==k),]) %*% matrix(beta0CV,ncol=1))/((table(cCV)[k])^approxh)
        }

        estCV <- Rsolnp::solnp(beta0CV, saauc, eqfun=constr, eqB=1, Xmat=dataCV[,-1], yvec=dataCV[,1], hc=hvalcCV,
                              lambdapar=lambdaval,control=list("outer.iter"=10^3, "inner.iter"=10^4, "delta"=stepsz,
                                                               "tol"=tolval,"trace"=0))

        #####
        wtvalsCV <- prop.table(table(oCV,cCV),1)[2,]

        convnow <- max(convnow, estCV$convergence)
        perfc <- rep(NA, numcov-1)
        for(l in 1:(numcov-1)){
          perfc[l]=Hmisc::somers2(pCV[which(cCV==l),] %*% estCV$pars, oCV[which(cCV==l)])[1]
        }
        perfnew <- Hmisc::somers2(pnew[which(cnew==i),] %*% estCV$pars, onew[which(cnew==i)])[1]
        rsltsCV[i,] = c(perfnew, (perfnew - (wtvalsCV %*% perfc))^2)
      }
      as.numeric(c(rsltsCV,convnow,rglmfailCV))
    }
    doCV <- foreach(ind=1:length(lambdaseq), .inorder=FALSE, .combine="rbind",
                    .packages=c("Hmisc","Rsolnp","aucm", "survival"),
                    .export=c("saveAllWarnings","saauc","%ni%","constr",
                              "numcov","pnew","onew","cnew","conditional",
                              "initialval","approxh","stepsz","tolval")) %dopar% {
                      cbind(lambdaseq[ind],
                            matrix(saveAllWarnings(lambdaseqCV(lambdaval=lambdaseq[ind]),logFile=warnfileCV),
                                   ncol=2*numcov+2))
                    }

    aAUCCV <- apply(doCV, 1, function(x) x[2:(numcov+1)] %*% wtvals)
    varCVTR <- apply(doCV, 1, function(x) wtvals %*% x[(numcov+2):(2*numcov+1)])
    varCVTE <- apply(doCV, 1, function(x) wtvals %*% (x[2:(numcov+1)] - (x[2:(numcov+1)] %*% wtvals))^2)

    #### Make plots
    layout(matrix(c(1,0,2,0,3),ncol=5), widths=c(3.5,0.5,3.5,0.2,2.5))

    par(xpd=TRUE)
    plot(log10(rsltsTR[,1]), rsltsTR[,3],lty=1,lwd=2,type="l",ylim=c(min(rsltsTR[,c(3,5:(numcov+4))],
                                                                       aAUCTRglm,aAUCTRrglm),
                                                                   max(rsltsTR[,c(3,5:(numcov+4))],
                                                                       aAUCTRglm,aAUCTRrglm)),
         main="Training Data",xlab="log10(lambda)",ylab="AUC")
    Hmisc::minor.tick(ny=5)
    for(j in 1:numcov){
      lines(log10(rsltsTR[,1]), rsltsTR[,j+4],lty=1,lwd=1,type="l",col="gray")
    }
    lines(log10(rsltsTR[,1]), rsltsTR[,3],lty=1,lwd=2,type="l")
    abline(h=aAUCTRglm, lwd=1.5, lty=2,xpd=FALSE)
    abline(h=aAUCTRrglm, lwd=1.5, lty=4,xpd=FALSE)

    par(new = T)
    plot(log10(rsltsTR[,1]), 100*(rsltsTR[,4]), col="red",type="l",lwd=2, axes=F, xlab=NA, ylab=NA,
         ylim=c(0,max(100*(rsltsTR[,4]),sqrt(varTRglm)*100,sqrt(varTRrglm)*100)))
    axis(side = 4)
    mtext(side = 4, line = 3, 'SD (x 100)',cex=0.7)
    abline(h=sqrt(varTRglm)*100, lwd=1.5, lty=2, col="red",xpd=FALSE)
    abline(h=sqrt(varTRrglm)*100, lwd=1.5, lty=4, col="red",xpd=FALSE)

    plot(log10(doCV[,1]), aAUCCV, lty=1,lwd=2,type="l",ylim=c(min(doCV[,c(2:(numcov+1))]),
                                                               max(doCV[,c(2:(numcov+1))])),
         main="LOCOCV",xlab="log10(lambda)",ylab="AUC")
    Hmisc::minor.tick(ny=5)
    for(j in 1:numcov){
      lines(log10(doCV[,1]), doCV[,j+1],lty=1,lwd=1,type="l",col="gray")
    }
    lines(log10(doCV[,1]), aAUCCV,lty=1,lwd=2,type="l")

    par(new = T)
    plot(log10(doCV[,1]), 100*sqrt(varCVTR), col="red",type="l",lwd=2, axes=F, xlab=NA, ylab=NA,
         ylim=c(0,max(100*sqrt(varCVTR),100*sqrt(varCVTE))))
    lines(log10(doCV[,1]), 100*sqrt(varCVTE), col="blue",type="l",lwd=2)
    axis(side = 4)
    mtext(side = 4, line = 3, 'SD (x 100)',cex=0.7)

    plot(1, type = "n", axes=FALSE, xlab="", ylab="")
    legend("center",c("Adjusted","Covariate-specific","SD (relative to TR)","SD (relative to TE)","GLM","rGLM"), lty=c(1,1,1,1,2,4), lwd=c(2,1,2,2,1,1),
           col=c(1,"gray","red","blue","black","black"), cex=0.85)

    plotOut <- recordPlot()

    if(max(rsltsTR[,2])>0){
      warning("SaAUC algorithm failed to converge for some value(s) of lambda in training data")
    }
    if(max(doCV[,(2*numcov)+2])>0){
      warning("SaAUC algorithm failed to converge for some value(s) of lambda for some strata in CV")
    }
    if(max(doCV[,(2*numcov)+3])>0){
      warning("rGLM failed to converge for some strata in CV; GLM used instead")
    }
    if(file.exists(warnfileTR)){
      if(file.info(warnfileTR)$size > 0){
        warning("warnings produced in training data (see warnfileTR)")
      }
    }
    if(file.exists(warnfileCV)){
      if(file.info(warnfileCV)$size > 0){
        warning("warnings produced in CV (see warnfileCV)")
      }
    }
  }
  if(rsltsupplied$convergence > 0){
    warning("SaAUC algorithm failed to converge for supplied value of lambda in training data")
  }
  if(!(robustmod$convergence)){
    warning("rGLM algorithm failed to converge in training data; GLM used instead")
  }

  if(length(lambdaseq)>0){
    return(list(NumCov=numcov, FittedCombs=list(InitialVal=beta0, NormGLM=normglm, NormrGLM=normrglm,
                                           MaxSaAUCSupplied=rsltsupplied$pars),
                aAUCTR=c("aAUCTRrGLM"=aAUCTRrglm, "aAUCTRGLM"=aAUCTRglm, "aAUCTRsuppl"=aAUCTRsuppl),
                varTR=c("varTRrGLM"=varTRrglm,"varTRGLM"=varTRglm,"varTRsuppl"=varTRsuppl),
                TRrslts=cbind("lambda"=rsltsTR[,1], "aAUCTR"=rsltsTR[,3], "varTR"=rsltsTR[,4]^2,
                              "AUCcTR"=rsltsTR[,5:(numcov+4)], "TRconv"=as.numeric(rsltsTR[,2]>0)),
                CVrslts=cbind("lambda"=doCV[,1], "aAUCCV"=aAUCCV, "varCVTR"=varCVTR, "varCVTE"=varCVTE, "AUCcCVTE"=doCV[,2:(numcov+1)]),
                CVconv=cbind("lambda"=doCV[,1], "SaAUC"=as.numeric(doCV[,(2*numcov)+2]>0), "rGLM"=doCV[,(2*numcov)+3]),
                plotOut=plotOut))
  }else{
    return(list(NumCov=numcov, FittedCombs=list(InitialVal=beta0, NormGLM=normglm, NormrGLM=normrglm,
                                                MaxSaAUCSupplied=rsltsupplied$pars),
                aAUCTR=c("aAUCTRrGLM"=aAUCTRrglm, "aAUCTRGLM"=aAUCTRglm, "aAUCTRsuppl"=aAUCTRsuppl),
                varTR=c("varTRrGLM"=varTRrglm,"varTRGLM"=varTRglm,"varTRsuppl"=varTRsuppl)))
  }
}
