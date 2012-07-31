####################### Entire Updated Estimation #########################
# This function is constructed to estimate the iterated paramters according 
# to Bai(2009), and Bada and Kneip (2010). 
# Input:
#	  1. dat.matrix	=  is a matrix where the first colomn containing 
#		 the(NTx1) vector of Ys and the remaining colomns containing 
#		the (NTxP) vector of Xs
# 	  2. dat.dim         dat.dim[1] is nr (nbr of rows) and dat.dim[2] 
#		 isnc (nbr of colomns) of the panel matrix Y.
#	  3. dim.criterion	= c("PC1", "PC2", "PC3", "IC1", "IC2" , "IC3"
#		, "IPC1", "IPC2", "ED", "ER", "GR")
#	  4. factor.dim	= the number of factors if it is known (standard is 
#		NULL)	 
#	  5. d.max 		= the maximum number of factos (an argument needed 
#		for the dimension selction)
#	  6. sig2.hat = is an argument needed for the dimension selction
#	  9. double.iteration =  logical argument, if TRUE the function will
#		estimate the optimal dimension in an outer iteration 
#		the convergence of all the paramters is obtained in a double 
#		iteration. If FALSE the dimension will be estimated parallelly 
#		with beta, lambda and F in order to reduce the number of 
#		computation (desadvantege: convergence to a local a local 
#		optimum).
#		This argument will be neglected if factor.dim is specified.
# Output:
#	  1. $PCA	= svd.pca object calculated at the optimal dimension 
#			(or given factor.dim)
#	  2. $beta	= optimal dimension according to the dimensionalty 
#			criterion
#	  3. $opt.d = the slope estimator of the observed regressors 
#			(beta.eup)
#	  4. $nbr.iterations = number of iteration
##########################################################################


FUN.Eup <- function(dat.matrix, dat.dim
		, double.iteration = double.iteration
		, dim.criterion, factor.dim, d.max, sig2.hat
		, start.beta, max.iteration, convergence){


#### data 
	y 	<- dat.matrix[, 1, drop = FALSE]
	x 	<- dat.matrix[,-1, drop = FALSE]
	nr 	<- dat.dim[1]
	nc	<- dat.dim[2]
	P	<- dat.dim[3]# or ncol(x)

#### if  d.max not given then d.max will be setted to sqrt(min(N, T))
	if(is.null(d.max)) d.max <- round(sqrt(min(nr, nc)))

#### start value #####
	beta.0 <- start.beta
	if(is.null(beta.0)) 
		{
		ymats <- matrix(y, nr, nc)
		xmats <- matrix(x, nr, (nc*P))
		zmats <- cbind(ymats, xmats)
		trzma <- zmats - svd.pca(zmats, given.d = d.max )$Q.fit
		tryxm <- matrix(trzma, (nr*nc), (P+1))
		try   <- tryxm[, 1, drop = FALSE]
		trx   <- tryxm[, -1, drop = FALSE]
		beta.0 <- coef(lm(try ~ -1 + trx ))
		}

#### calculate the inverse once in order to reduce computations of the 
#### iterated slope estimator
	FUN.ols.beta <- function(updated.y, x, inv.xx.x){
		beta <- tcrossprod(inv.xx.x, t(updated.y))
		}

	xx 	   <- crossprod(x)
	inv.xx   <- solve(xx)
	inv.xx.x  <- inv.xx%*%t(x)


#### define given.d = factor.dim and set factor.dim= d.max if it is null
	given.d <- factor.dim
	if(is.null(factor.dim)) factor.dim <- d.max

#### the inner Iteration function
	inner.iteration <- function(y, x, inv.xx.x =inv.xx.x 
				  , beta.0 = beta.0
				  , factor.dim = factor.dim, d.max=d.max
				  , sig2.hat= sig2.hat
				  , double.iteration = double.iteration
				  , max.iteration = max.iteration, convergence = convergence
				  , past.iterations = past.iterations, i=1){
  # Iteration (0): initial cumputations
	# w.0 = y-x%*%beta.0 
		w.0 <- y - tcrossprod(x, t(beta.0))

  	# W.0: write w.0 in a matrix form 
		W.0 <- matrix(w.0, nr, nc)

 	# PCA.0: PCA.0 computation, OptDim.0 and y.fitted
		PCA.0    <- svd.pca(W.0, given.d=factor.dim)
		if(!double.iteration && is.null(given.d)){
			OptDim.0   <- EstDim(PCA.0
					   , dim.criterion = dim.criterion 
					   , d.max = d.max, sig2.hat= sig2.hat)
			opt.dim.0  <- OptDim.0[,2]
			factor.dim <- opt.dim.0
			y.fitted.0 <- tcrossprod(PCA.0$L[, 0:opt.dim.0
					, drop = FALSE])%*%W.0
		}

		else y.fitted.0 = PCA.0$Q.fit

   # Iteration (+1)
  	# y.updated.0: updat y.updated.0 = y - fs.0
		y.updated.1 <- y -  c(y.fitted.0)

  	# beta.1: OLS.1 computation for the computed fs.0 in interation 0 
		beta.1 <- FUN.ols.beta(y.updated.1, x, inv.xx.x) 

  	# convergence condition
		if(all( abs((beta.0 - beta.1)) < convergence)| (i + past.iterations) == max.iteration){
			if((i + past.iterations) == max.iteration) {
			warning(paste("The maximal number of iterations is achieved ", max.iteration), call. = FALSE)
				if(double.iteration && is.null(given.d)){
				OptDim.0   <- EstDim(PCA.0
					   , dim.criterion = dim.criterion 
					   , d.max = d.max, sig2.hat= sig2.hat)
				opt.dim.0  <- OptDim.0[,2]
				factor.dim <- opt.dim.0
				y.fitted.0 <- tcrossprod(PCA.0$L[, 0:opt.dim.0
					, drop = FALSE])%*%W.0
				}
			
			}
			Result <- list(PCA=PCA.0, beta=beta.0
			,factor.dim=factor.dim, Nbr.Iterations = i)
			Result
			}

		else inner.iteration(y=y, x=x, inv.xx.x =inv.xx.x 
			,beta.0 = beta.1,factor.dim = factor.dim 
			,d.max=d.max,sig2.hat= sig2.hat 
			,double.iteration = double.iteration
			,max.iteration = max.iteration, convergence = convergence
			,past.iterations = past.iterations, (i+1))
	}

####  integrating the inner iterateion in the outer iteration function

	entire.iteration <- function(y=y, x=x, inv.xx.x =inv.xx.x 
				    , beta.0 = beta.0, factor.dim = factor.dim 
				    , d.max=d.max, sig2.hat= sig2.hat 
				    , double.iteration = double.iteration
				    , max.iteration = max.iteration
				    , convergence = convergence
				    , past.iterations = 0, l =1){
	# first inner iteration 
		In.Iter.0 <- inner.iteration(y=y, x=x, inv.xx.x =inv.xx.x 
				, beta.0 = beta.0, factor.dim = factor.dim
				, d.max=d.max, sig2.hat= sig2.hat
				, double.iteration = double.iteration
				, max.iteration = max.iteration, convergence = convergence
				, past.iterations = past.iterations
				, i =1)
		pca.d0 	  <- In.Iter.0$PCA
		beta.d0 	  <- In.Iter.0$beta 
		opt.d0	  <- In.Iter.0$factor.dim
		nbr.iteration <- In.Iter.0$Nbr.Iterations + past.iterations

		# if double.iteration is TRUE select new opt.d iteratively
		if(double.iteration && is.null(given.d)){	
		  # 1 new optimal dimension
		opt.dim1 <- EstDim(pca.d0, dim.criterion = dim.criterion 
				, d.max = d.max, sig2.hat= sig2.hat)
		opt.d1  <- opt.dim1[,2]
		  # convergence condition
			if(opt.d1==opt.d0| nbr.iteration >= max.iteration){
				#if(l >= d.max) warning(expression(" The maximal
				#	number of outer iterations is 
				#	achieved (d.max)"))
				Result  <- list(y=y, x=x, dat.dim= dat.dim
						, PCA = pca.d0, beta= beta.d0
						, opt.d =opt.d1
						, nbr.iterations= nbr.iteration)
				Result
				}
			else entire.iteration(y=y, x=x, inv.xx.x =inv.xx.x 
					, beta.0 = beta.d0, factor.dim = opt.d1
					, d.max=d.max, sig2.hat= sig2.hat
					, double.iteration = double.iteration
					, past.iterations = nbr.iteration
					, max.iteration = max.iteration, convergence = convergence
					, l = (l+1))
		}
		else {
		Result <- list(y=y, x=x, dat.dim= dat.dim
			, PCA = pca.d0, beta= beta.d0, opt.d =opt.d0
			, nbr.iterations=nbr.iteration)
		Result
		}
	}

###### entire iteration result 

	Result	 <- entire.iteration(y=y, x=x, inv.xx.x =inv.xx.x 
			    , beta.0 = beta.0 , factor.dim = factor.dim
			    , d.max=d.max, sig2.hat= sig2.hat 
			    , double.iteration = double.iteration
			    , max.iteration = max.iteration, convergence = convergence
			    , past.iterations = 0, l =1)

	Result
  }





############################### Eup.default ###############################
# The rapper function of FUN.Eup(). Here, the data can be given in a formula
# the user has the possibility to transform the data before iterating the estimator
###########################################################################
Eup.default <- function(formula,
		additive.effects = c("none", "individual", "time", "twoways"),
		dim.criterion	 = c("PC1", "PC2", "PC3", "IC1", "IC2" , "IC3", "IPC1", "IPC2", "IPC3" , "ED"),
		d.max            = NULL,
		sig2.hat         = NULL,
		factor.dim       = NULL,
		double.iteration = TRUE,
		start.beta       = NULL,
		max.iteration    = 500,
		convergence      = 1e-6,
		restrict.mode    = c("restrict.factors", "restrict.loadings"),
                ...){

### substruct data from fomrmula and perfome a transformation according 
### to additive.effects 

  # check fomula

	if(!class(formula)=="formula"){
		stop("\n Argument >>formula<< needs a formula-object like 
					y~x1+... where the elements are matrices")
	}
  
  # names 
  
  names  <- names(model.frame(formula))

  # define additive.effects
	
	additive.effects <- match.arg(additive.effects)
	sig2.hat.dim <- sig2.hat

  ## Transform the response variable as well as the 'P' regressors and give them in a list where ech 
    ## componente contains also a list with:
    #       1- "Tr" Name of the transformation
    #       2- "I" Logical variable if ther is intercept or no
    #       3- "ODM" Original Data matrix
    #       4- "TDM" Transformed Data in a matrix
    #       5- "TDV" Transformed Data in a NT x 1 Vector
    #       6- "TRm" Sublist with
    #           a- "OVc" Overall Constant
    #           b- "InC" time constant individual effects
    #           c- "TiVC" additive time varying effects

	PF.obj	<- FUN.Pformula(formula = formula
						, effect = additive.effects)    
	nc 		<- ncol(PF.obj[[1]]$ODM)#N
	nr 		<- nrow(PF.obj[[1]]$ODM)#T
	P  		<- length(PF.obj) - 1
  	intercept <- PF.obj[[1]]$I

  # prepare for FUN.default (give the data in a large matrix: (y, x1, ...xP) where y, xp are  NT x1 Vectors  )

	dat.dim 	  <- c(nr, nc, P)
  dat.matrix	  <- sapply(1:(P+1), function(i) PF.obj[[i]]$TDV)

	dim.criterion <- match.arg(dim.criterion)
	
  # Estimation results

	tr.model.est <- FUN.Eup(dat.matrix		= dat.matrix
					, dat.dim		= dat.dim
					, double.iteration= double.iteration
					, start.beta	= start.beta
					, dim.criterion 	= dim.criterion
					, factor.dim	= factor.dim
					, d.max		= d.max
					, max.iteration   = max.iteration
					, convergence     = convergence
					, sig2.hat = sig2.hat)

    # Eup beta and Nbr.iteration

	Nbr.iteration	<- tr.model.est$nbr.iterations
	beta.Eup		<- as.matrix(tr.model.est$beta)
  	colnames(beta.Eup) <- ""
  	rownames(beta.Eup) <- names[-1]

    # additive effects

	# calculate the constant
	OvConst <- sapply(1:(P+1),  function(i) PF.obj[[i]]$TRm$OVc)
	ConsCoef <- OvConst[1] - OvConst[-1]%*%beta.Eup 
  
	# calculate the additive indivudal effects

	ind.means <- sapply(1:(P+1),  function(i) PF.obj[[i]]$TRm$InC)
	Ind.Eff <- ind.means[,1] - ind.means[,-1]%*%beta.Eup - c(ConsCoef)

	# calculate the additive time effecs

	tim.means <- sapply(1:(P+1),  function(i) PF.obj[[i]]$TRm$TiVC)
	Tim.Eff <- tim.means[,1] - tim.means[,-1]%*%beta.Eup - c(ConsCoef)


    # factor dimension 
	used.dim <- tr.model.est$opt.d
	proposed.dim <- ifelse(is.null(factor.dim), "Not Indicated"
					, factor.dim)
	if(is.null(factor.dim)){
		optimal.dim <- tr.model.est$opt.d
		}
	else{
		optimal.dim <- EstDim(tr.model.est$PCA
						, dim.criterion = dim.criterion 
						, d.max = d.max, sig2.hat= sig2.hat)[,2]
		}

    # factor structure and resuduals 

	restrict.fs.a.resid	<- restrict.pca(tr.model.est$PCA
							, restrict.mode=restrict.mode)
	fs.and.resid		<- restrict.fs.a.resid$orig.values
	factors			<- restrict.fs.a.resid$factors
	loadings			<- restrict.fs.a.resid$loadings
	unob.fact.stru		<- restrict.fs.a.resid$fitted.values
	residuals			<- fs.and.resid - unob.fact.stru
  
    # fitted values
    
  orig.Y <- PF.obj[[1]]$ODM
  fitted.values <- orig.Y - residuals
  degrees.of.freedom <- (nr*nc - (nr+nc)*used.dim - P - 
                    intercept -
                    nc*(additive.effects == "individual"| additive.effects == "twoways") - 
                    nr*(additive.effects == "time"| additive.effects == "twoways"))
  sig2.hat <- sum(diag(crossprod(residuals)))/degrees.of.freedom

## Results
				
	final.result <- list(
    dat.matrix = dat.matrix
    , formula = formula
    , dat.dim = dat.dim
    , slope.para = beta.Eup
    , names = names
    , is.intercept = intercept
    , additive.effects = additive.effects
    , Intercept = c(ConsCoef)
    , Add.Ind.Eff = c(Ind.Eff)
    , Add.Tim.Eff = c(Tim.Eff)
    , unob.factors = factors
    , ind.loadings = loadings
    , unob.fact.stru = unob.fact.stru
    , used.dim= used.dim
    , proposed.dim= proposed.dim
    , optimal.dim = optimal.dim
    , factor.dim = factor.dim
    , d.max = d.max
    , dim.criterion = dim.criterion
    , OvMeans = OvConst
    , ColMean = ind.means
    , RowMean = tim.means
    , max.iteration = max.iteration
    , convergence = convergence
    , start.beta = start.beta
    , Nbr.iteration= Nbr.iteration    		
    , fitted.values = fitted.values
    , orig.Y = orig.Y
    , residuals = residuals
    , sig2.hat.dim = sig2.hat.dim
    , sig2.hat = sig2.hat
    , degrees.of.freedom = degrees.of.freedom
    , call = match.call())
  class(final.result) <- "Eup"
  return(final.result)  
      }
