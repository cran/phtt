

## Methods ========================================================================================

Eup <- function(formula,
    additive.effects = c("none", "individual", "time", "twoways"),
    dim.criterion    = c("PC1", "PC2", "PC3", "IC1", "IC2" , "IC3",
"IPC1", "IPC2", "IPC3" , "ED"),
    d.max            = NULL,
    sig2.hat         = NULL,
    factor.dim       = NULL,
    double.iteration = TRUE,
    start.beta       = NULL,
    max.iteration    = 500,
    convergence      = 1e-6,
    restrict.mode    = c("restrict.factors", "restrict.loadings"),
    ...){
  UseMethod("Eup")
}

print.Eup <- function(x,...){
  cat("Call:\n")
  print(x$call)
  
  cat("\nCoeff(s) of the Observed Regressor(s) :\n\n")
  slope.para <- x$slope.para
  if(x$is.intercept){
    inter <- matrix(x$Intercept, 1, 1)
    colnames(inter) <- ""
    rownames(inter) <- "(Intercept)"
    slope.para <- rbind(signif(inter,digits=3), signif(slope.para,digits=3))
    slope.para <- signif(slope.para, 3)
  }
  print(t(slope.para))
  cat("\nAdditive Effects Type: ", as.name(x$additive.effects)," \n")
  cat("\nDimension of the Unobserved Factors:", x$used.dim," \n")
  cat("\nNumber of iterations:", x$Nbr.iteration,"\n")
}



coef.Eup <- function(object,...){
    if(object$is.intercept)
    Intercept <- object$Intercept
    else Intercept <- NULL
    
    Slope.Coef <- object$slope.para
    
    if(object$additive.effects== "individual"| object$additive.effects== "twoways")
    Add.Ind.Eff <- object$Add.Ind.Eff
    else Add.Ind.Eff <- NULL

    if(object$additive.effects== "time"| object$additive.effects== "twoways")
    Add.Tim.Eff <- object$Add.Tim.Eff
    else Add.Tim.Eff <- NULL

    Common.factors <- object$unob.factors
    
    Ind.loadings.param <- object$ind.loadings
    
    Time.varying.ind.eff <- object$unob.fact.stru
    
    Factor.Dim <- object$used.dim

    Var.shares.of.loadings.param      <- numeric(Factor.Dim)
    Total.var.loadings.param          <- sum(apply(Ind.loadings.param,2,var))
    for(i in 1:Factor.Dim){
      Var.shares.of.loadings.param[i] <- round(var(c(Ind.loadings.param[,i]))/Total.var.loadings.param,
                    digits=4)*100
    }
    coef.list <- list(
        Intercept                    = Intercept,
        Slope.Coef                   = Slope.Coef,
        Add.Ind.Eff                  = Add.Ind.Eff, 
        Add.Tim.Eff                  = Add.Tim.Eff, 
        Common.factors               = Common.factors, 
        Ind.loadings.param           = Ind.loadings.param,
        Var.shares.of.loadings.param = Var.shares.of.loadings.param,
        Time.varying.ind.eff         = Time.varying.ind.eff,
        Factor.Dim                   = Factor.Dim)  
        
    return(coef.list)
}

residuals.Eup <- resid.Eup <- function(object,...){
  Residual.mat <- object$residuals
  return(Residual.mat)
}

summary.Eup <- function(object,...){
  ## Residuals:
  Res.outpt <- signif((summary(as.vector(object$residuals))), digits=3)[-4]
  names(Res.outpt) <- c("Min", "1Q", "Median", "3Q", "Max")
  yy <- sum(diag(crossprod(object$orig.Y - mean(object$orig.Y))))
  ee <- sum(diag(crossprod(object$residuals)))
  R2 <- signif(1 - ee/yy, 4)
  
  ## Add-Effect-Type:
  eff              <- matrix(object$additive.effects)
  colnames(eff)    <- ""
  rownames(eff)    <- ""
  
  ## Coefficients:
  TAB  <-  Eup.inference(Eup.Obj=object)$inf.result
  TAB <- signif(TAB, 3)
  
  ## Result:
  result        <- list(Res.outpt    = Res.outpt,
                        coefficients = TAB,
                        R2 = R2,
                        Eup.obj      = object)                
  class(result) <- "summary.Eup"
  result
}

print.summary.Eup <- function(x, ...){

  ## Call
  cat("Call:\n")
  print(x$Eup.obj$call)
  ## Residuals:
  cat("\nResiduals:\n")
  print(x$Res.outpt)
  cat("\n")
  ## Beta-Coeffs
  cat("\n Slope-Coefficients:\n")
  printCoefmat(x$coefficients)
  cat("\nAdditive Effects Type: ", as.name(x$Eup.obj$additive.effects)," \n")
  cat("\nDimension of the Unobserved Factors:", x$Eup.obj$used.dim," \n")
  #cat("\nOptimized Factor Dimension:         ", x$Eup.obj$optimal.dim," \n")
  cat("\nResidual standard error:", signif(x$Eup.obj$sig2.hat, 4), "on",
            x$Eup.obj$degrees.of.freedom, "degrees of freedom, ", "\nR-squared:", x$R2,"\n")
}

## print.summary.Eup <- function(x, ...){
##   ## Call
##   cat("Call:\n")
##   print(x$Eup.obj$call)
##   ## Residuals:
##   cat("\nResiduals:\n")
##   print(x$Res.outpt)
##   cat("\n")
##   ## Beta-Coeffs
##   cat("\n Slope-Coefficients:\n")
##   printCoefmat(x$coefficients)
  
##   cat("\nAdditive Effects Type: ",                   as.name(x$Eup.obj$additive.effects)," \n")
##   cat("\nUsed Dimension of the Unobserved Factors:", x$Eup.obj$used.dim)
## #  cat("\nOptimized Factor Dimension:              ", x$Eup.obj$optimal.dim," \n") 
##   cat("\nResidual standard error:",             signif(x$Eup.obj$sig2.hat, digits=3), "on", 
##                                                 x$Eup.obj$degrees.of.freedom, "degrees of freedom \n")
##   cat("Multiple R-squared:",                    signif(x$R2,digits=3),"\n")
## }

## print.summary.Eup <- function(x, ...){
##   ## Call
##   cat("Call:\n")
##   print(x$Eup.obj$call)
##   ## Residuals:
##   cat("\nResiduals:\n")
##   print(x$Res.outpt)
##   cat("\n")
##   ## Beta-Coeffs
##   cat("\n Slope-Coefficients:\n")
##   printCoefmat(x$coefficients)
  
##   cat("\nAdditive Effects Type: ", as.name(x$Eup.obj$additive.effects)," \n")
##   cat("\nDimension of the Unobserved Factors:", x$Eup.obj$used.dim," \n")
##   #cat("\nOptimized Factor Dimension:         ", x$Eup.obj$optimal.dim," \n")
  
##   cat("\nResidual standard error:", signif(x$Eup.obj$sig2.hat,digits=3), "on", 
##             x$Eup.obj$degrees.of.freedom, "degrees of freedom, ", "\nR-squared:", signif(x$R2,digits=3),"\n")
## }

plot.summary.Eup <- function(x,...){
  if(is.null(x$Eup.obj$unob.factors) & x$Eup.obj$additive.effects=="none"){
    stop("Neither an estimated factor structure nor additive effects to plot.")
  }
  if(!is.null(x$Eup.obj$unob.factors)){
    if(x$Eup.obj$additive.effects=="none"){
      par(mfrow=c(1,2))
      matplot(x$Eup.obj$unob.factors,
            main=paste("Estimated Factors\n(Used Dimension = ",x$Eup.obj$used.dim,")",sep=""),
            xlab="Time",ylab="", type="l",...)
      matplot(x$Eup.obj$unob.fact.stru,
            main=paste("Estimated Factor-Structure"),
            xlab="Time",ylab="", type="l",...)
    par(mfrow=c(1,1))
    }
    if(x$Eup.obj$additive.effects=="time"){
      par(mfrow=c(1,3))
      plot.ts(x$Eup.obj$Add.Tim.Eff, main="Additive Time Effect", ylab="",xlab="Time",...)
      matplot(x$Eup.obj$unob.factors,
            main=paste("Estimated Factors\n(Used Dimension = ",x$Eup.obj$used.dim,")",sep=""),
            xlab="Time",ylab="", type="l",...)
      matplot(x$Eup.obj$unob.fact.stru,
            main=paste("Estimated Factor-Structure"),
            xlab="Time",ylab="", type="l",...)
    par(mfrow=c(1,1))
    }
    if(x$Eup.obj$additive.effects=="twoways"){
      par(mfrow=c(1,4))
      plot.ts(x$Eup.obj$Add.Tim.Eff, main="Additive Time Effect", ylab="",xlab="Time",...)
      matplot(matrix(rep(x$Eup.obj$Add.Ind.Eff,each=x$Eup.obj$dat.dim[1]),
                     nrow=x$Eup.obj$dat.dim[1],ncol=x$Eup.obj$dat.dim[2]),
                     main="Additive Individual Effects", ylab="",xlab="Time", type="l", ...)
      matplot(x$Eup.obj$unob.factors,
            main=paste("Estimated Factors\n(Used Dimension = ",x$Eup.obj$used.dim,")",sep=""),
            xlab="Time",ylab="", type="l",...)
      matplot(x$Eup.obj$unob.fact.stru,
            main=paste("Estimated Factor-Structure"),
            xlab="Time",ylab="", type="l",...)
    par(mfrow=c(1,1))
    }
    if(x$Eup.obj$additive.effects=="individual"){
      par(mfrow=c(1,3))
      matplot(matrix(rep(x$Eup.obj$Add.Ind.Eff,each=x$Eup.obj$dat.dim[1]),
                     nrow=x$Eup.obj$dat.dim[1],ncol=x$Eup.obj$dat.dim[2]),
                     main="Additive Individual Effects", ylab="",xlab="Time", type="l", ...)
      matplot(x$Eup.obj$unob.factors,
            main=paste("Estimated Factors\n(Used Dimension = ",x$Eup.obj$used.dim,")",sep=""),
            xlab="Time",ylab="", type="l",...)
      matplot(x$Eup.obj$unob.fact.stru,
            main=paste("Estimated Factor-Structure"),
            xlab="Time",ylab="", type="l",...)
    par(mfrow=c(1,1))
    }
  }else{
    if(x$Eup.obj$additive.effects=="time"){
      par(mfrow=c(1,1))
      plot.ts(x$Eup.obj$Add.Tim.Eff, main="Additive Time Effect", ylab="",xlab="Time",...)
      par(mfrow=c(1,1))
    }
    if(x$Eup.obj$additive.effects=="twoways"){
      par(mfrow=c(1,2))
      plot.ts(x$Eup.obj$Add.Tim.Eff, main="Additive Time Effect", ylab="",xlab="Time",...)
      matplot(matrix(rep(x$Eup.obj$Add.Ind.Eff,each=x$Eup.obj$dat.dim[1]),
                     nrow=x$Eup.obj$dat.dim[1],ncol=x$Eup.obj$dat.dim[2]),
              main="Additive Individual Effects", ylab="",xlab="Time", type="l", ...)
      par(mfrow=c(1,1))
    }
    if(x$Eup.obj$additive.effects=="individual"){
      par(mfrow=c(1,1))
      matplot(matrix(rep(x$Eup.obj$Add.Ind.Eff,each=x$Eup.obj$dat.dim[1]),
                     nrow=x$Eup.obj$dat.dim[1],ncol=x$Eup.obj$dat.dim[2]),
                     main="Additive Individual Effects", ylab="",xlab="Time", type="l", ...)
    par(mfrow=c(1,1))
    }
  }
}

