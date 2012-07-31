################################## Eup slope inference ####################
## Input:
#     dat.matrix	= the data in matrix form where the first colomn
#				  contain NT vector of Y second one the NT vector 
#				  of the first regressor X
# 		dat.dim	= the dimension of the data N and T
# 		used.d	= used dimension d in the interative procedure
#	 	beta.Eup	= the estimated slope estimator for given d
# 		factors	= common factors after scaling according the used 
#				  restriction
# 		loadings	= individual loading parameters after scaling
#				 (according restriction)
#		residuals	= the residual terms
## Output:
#		Eup slope Estimate
#		std
#		Pr(>|z|)
##########################################################################

Eup.inference <- function(Eup.Obj){

## collect informations from the Eup.Obj

	y 	<- Eup.Obj$dat.matrix[, 1, drop = FALSE]
	x 	<- Eup.Obj$dat.matrix[,-1, drop = FALSE]
	nr  <- Eup.Obj$dat.dim[1]
	nc	<- Eup.Obj$dat.dim[2]
	P	 <- Eup.Obj$dat.dim[3]
  	d   <- Eup.Obj$used.dim
	slope <- Eup.Obj$slope.para
  Intercept <- matrix(Eup.Obj$Intercept, 1, 1)
  is.intercept <- Eup.Obj$is.intercept
  OvMeans <- Eup.Obj$OvMeans
	F	 <- Eup.Obj$unob.factors
	A	 <- Eup.Obj$ind.loadings
	sig2.hat <- Eup.Obj$sig2.hat
  



## Projection matrix of the factors F
  
  I.TxT <- diag(1, ncol= nr, nrow= nr)
  if(d==0) M.F   <- I.TxT
  else M.F   <- I.TxT - tcrossprod(F)/nr

## Projection matrix of the loadings A

  I.NxN <- diag(1, ncol= nc, nrow= nc)
  if(d==0) M.A   <- I.NxN
  else{
    S  <- diag(diag(crossprod(A))^{-0.5}, d)
    AS <- tcrossprod(A,S)
    P.A   <- tcrossprod(AS)
    M.A   <- I.NxN - P.A
  }
    

## write the x matrices in a list: each regressor is written in a 
## list component
	X.mat.list <- NULL
	for(p in 1:P) X.mat.list[[p]] <- matrix(x[,p], nr, nc)

  # Z_i = M.F * X_i - sum{M.F * X_k*a_ik}/n
	Z.list	<- sapply (X.mat.list, function(X) M.F %*% X %*% M.A 
				, simplify = FALSE)

  # construct the matrix D= sum Z_i'Z_i/NT
	Z		<- sapply (Z.list, function(Z) c(Z), simplify = TRUE)
	ZZ     <- crossprod(Z)/(nr*nc)
	inv.ZZ <- solve(ZZ)
	asy.var <- (inv.ZZ * sig2.hat)/(nr*nc)
  mpp <- sqrt(diag(asy.var))  
  
  ## Add intercept if it exists in the formula
 if(is.intercept){
  rownames(Intercept) <- "(Intercept)"
  colnames(Intercept) <- " "
  slope <- rbind(Intercept, slope)
  asy.var.inter <- sig2.hat/(nc*nr) + 
          matrix(OvMeans[-1], 1, P, byrow = TRUE)%*%asy.var%*%matrix(OvMeans[-1], P, 1)
  mpp11 <-  sqrt(c(asy.var.inter))
  mpp <- c(mpp11, mpp)
  }
  
	test <- slope/mpp 
	p.value <- (1 - pnorm(abs(test)))*2
	inf.result <- cbind(slope,
                            mpp,
                            test,
                            p.value)
	colnames(inf.result) <- c("Estimate", "Std.Err", "Z value", "Pr(>z)")
	result <- list(ZZ = ZZ, inv.ZZ = inv.ZZ, inf.result=inf.result, sig2.hat=sig2.hat) 	
  return(result)
}
