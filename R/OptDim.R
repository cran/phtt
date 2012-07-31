##################################### Bai ##################################
bai.dim.opt <- function(Obj, d.max = NULL, sig2.hat = NULL)
	{
	nr    <- Obj$nr
	nc    <- Obj$nc
	w 	<- Obj$V.d/(nr*nc)
	d.seq <- Obj$d.seq
	C     <- min(nr, nc)
	max.rk <- length(w)
	d.max <- ifelse(is.null(d.max),min(round(sqrt(C)), max.rk ), d.max)
	if(d.max > max.rk)
		{
		warning(expression("The given maximun dimension 'd.max' is larger than the postive eingen values of the corresponding covariance matrix. We use maximum rank as 'd.max' "))
		d.max <- max.rk
		}

	ref   <- ifelse(is.null(sig2.hat),w[(d.max+1)], sig2.hat)
	rmw   <- w - ref
	dmax  <- length(rmw[rmw>= 0]) - 1

	PC1 = w + ref * d.seq * ((nr+nc)/(nr*nc)) * log((nr*nc)/(nr + nc))
	d.opt.PC1 <- order(PC1)[1]-1	
	PC2 = w + ref *d.seq * ((nr+nc)/(nr*nc)) * log(C)
	d.opt.PC2 <- order(PC2)[1]-1	
	PC3 = w + ref * d.seq * (log(C)/C)
	d.opt.PC3 <- order(PC3)[1]-1	
        ## number of eigenvalues greater than zero: 
        d.max.p1.gzero <- length(w[1:(d.max+1)][w[1:(d.max+1)]>0])
        ##
	IC1 = log(w[1:d.max.p1.gzero]) + d.seq[1:d.max.p1.gzero] * ((nr+nc)/(nr*nc)) * log((nr*nc)/(nr + nc))
	d.opt.IC1 <- order(IC1)[1]-1	
	IC2 = log(w[1:d.max.p1.gzero]) + d.seq[1:d.max.p1.gzero] * ((nr+nc)/(nr*nc)) * log(C)
	d.opt.IC2 <- order(IC2)[1]-1	
	IC3 = log(w[1:d.max.p1.gzero]) + d.seq[1:d.max.p1.gzero] * (log(C)/C)
	d.opt.IC3 <- order(IC3)[1]-1	

	IPC1 = w + ref * d.seq * (nr/(4*log(log(nr)))) *((nr+nc)/(nr*nc)) * log((nr*nc)/(nr + nc))
	d.opt.IPC1 <- order(IPC1)[1]-1	
	IPC2 = w + ref *d.seq * (nr/(4*log(log(nr)))) *((nr+nc)/(nr*nc)) * log(C)
	d.opt.IPC2 <- order(IPC2)[1]-1
	IPC3 = w + ref * d.seq * (nr/(4*log(log(nr)))) *((nr+nc - ref)/(nr*nc)) *(log(nr*nc))
	d.opt.IPC3 <- order(IPC3)[1]-1

	result <- matrix(c(d.opt.PC1, d.opt.PC2, d.opt.PC3, d.opt.IC1
				, d.opt.IC2, d.opt.IC3, d.opt.IPC1, d.opt.IPC2
				, d.opt.IPC3, w[d.opt.PC1+1], w[d.opt.PC2+1]
				, w[d.opt.PC3+1], w[d.opt.IC1 + 1], w[d.opt.IC2+1]
				, w[d.opt.IC3+1], w[d.opt.IPC1+1], w[d.opt.IPC2+1]
				, w[d.opt.IPC3+1]), 9,2)
	all.crit <- c("PC1", "PC2",  "PC3", "IC1"
				, "IC2", "IC3", "IPC1", "IPC2", "IPC3")
	Result <- data.frame(I(all.crit), result
				, rep.int(ref, 9), rep(min(dmax, d.max), 9))
	colnames(Result) <- c("Criterion ", "Optimal Dimension", "sd2"
				, "sd2.hat.ref", "d.ref.max")
	return(Result)
	}


B.OptDim <- function(Obj, criteria = c("PC1","PC2","PC3","IC1","IC2","IC3"
			, "IPC1","IPC2","IPC3") , d.max = NULL, sig2.hat = NULL){
	# what is Obj?
	if(class(Obj)=="svd.pca"|class(Obj)=="fsvd.pca") obj <- Obj
	else{
		if(class(Obj)=="pca.fit"|class(Obj)=="fpca.fit"){
			nr  <- Obj$data.dim[1]
			nc  <- Obj$data.dim[2]
			V.d <- Obj$Sd2*(nr*nc)
			d.seq = seq.int(0, (length(V.d)-1))
			obj <- list(V.d = V.d, nr = nr, nc = nc, d.seq = d.seq)
			}	
		else{
			if(is.regular.panel(Obj)) obj <- svd.pca(Obj)
			else{
				if(!is.vector(Obj[[1]])|!is.numeric(Obj[[1]])
				  |!is.numeric(Obj[[2]])|length(Obj[[2]])!=2)
				  stop(c("'Obj' does not have the correct form."))
				else{# the function can deal with a list containing a vector of RSS for each d in the first listcomponent and the dimension as 2 dimensional vector in the second component 
					nr <- Obj[[2]][1]
					nc <- Obj[[2]][2]
					V.d <- Obj[[1]]
					d.seq = seq.int(0, (length(V.d)-1))
					obj <- list(V.d =V.d, nr = nr, nc = nc
						, d.seq = d.seq)
					}
				}
			}
		}
	result <- bai.dim.opt(obj, d.max = d.max, sig2.hat = sig2.hat)
	criteria <- match.arg(criteria, several.ok = TRUE)
	return(result[result[,1] %in% criteria, ])
	}

#####################################################################################################################

onatski.dim.opt <- function(svd.pca.obj, d.max = NULL)
	{
	nr      <- svd.pca.obj$nr
	nc      <- svd.pca.obj$nc
	w 	  <- svd.pca.obj$V.d/(nr*nc)
	ev      <- svd.pca.obj$E/nr
	max.rk  <- length(ev)

	exa.ev  <- c(sum(ev), ev)
	d.max   <- ifelse(is.null(d.max), max.rk-5, min(d.max, (max.rk-5)))	
	j       <- max(1, (d.max - 1))

	c.reg   <- as.matrix(seq((j - 1), (j + 3))^{(2/3)}, 4, 1)
	delta   <- 2* abs(coef(lm.fit(c.reg, ev[j:(j+4)])))
	dist.ev <- exa.ev[1:(max.rk-1)] -  exa.ev[2:max.rk] - delta

	FUN.o.dim <- function(o.dist.ev, d)
		{
		if(o.dist.ev[d] < 0) (d - 2)
		else {
			if (d == (max.rk-1)){
				(d - 2)
				warning(expression("ED faild to find an appropriate dimension")) 
				}
			else FUN.o.dim(o.dist.ev, (d+1))
			}
		}
	
	d.opt.i <- FUN.o.dim(dist.ev, 1)

	if(d.opt.i==d.max) {
			result <- data.frame(I("ED")
					, matrix(c(d.opt.i, w[d.opt.i+1]), 1, 2))
			colnames(result) <- c("Criterion", "Optimal Dimension", "sd2")
			return(result)
			}
	else {
		if (d.opt.i == (max.rk-1)) {
			result <- data.frame(I("ED")
					, matrix(c(d.opt.i, w[d.opt.i+1]), 1, 2))
			colnames(result) <- c("Criterion", "Optimal Dimension", "sd2")
			return(result)
			warning(expression("ED faild to find an appropriate dimension")) 
			}
		else onatski.dim.opt(svd.pca.obj, d.opt.i)
		}
	}

O.OptDim <- function(Obj, d.max = NULL){
	# what is Obj?
	if(class(Obj)=="svd.pca"|class(Obj)=="fsvd.pca") obj <- Obj
	else{
		if(class(Obj)=="pca.fit"|class(Obj)=="fpca.fit"){
			nr  <- Obj$data.dim[1]
			nc  <- Obj$data.dim[2]
			V.d <- Obj$Sd2*(nr*nc)
			E   <- Obj$eigen.values*(nr*nc)
			obj <- list(V.d = V.d, nr = nr, nc = nc, E = E)
			}	
		else{
			if(is.regular.panel(Obj)) obj <- svd.pca(Obj)
			else{
				if(!is.vector(Obj[[1]])|!is.numeric(Obj[[1]])
				  |!is.numeric(Obj[[2]])|length(Obj[[2]])!=2)
				  stop(c("'Obj' does not have the correct form."))
				else{# the function can deal with a list containing a vector of RSS for each d in the first listcomponent and the dimension as 2 d-vector 'c(nr, nc)' in the second component 
					nr  <- Obj[[2]][1]
					nc  <- Obj[[2]][2]
					V.d <- Obj[[1]][-length(Obj[[1]])]
					E   <- -diff(Obj[[1]]-Obj[[1]][1])
					obj <- list(V.d = V.d, nr = nr, nc=nc, E = E)
					}
				}
			}
		}
	result <- onatski.dim.opt(obj, d.max = d.max)
	return(result)
	}

#####################################################################################################################



RH.dim.opt <- function(svd.pca.obj, d.max = NULL)
	{
	nr      <- svd.pca.obj$nr
	nc      <- svd.pca.obj$nc
	w 	  <- svd.pca.obj$V.d/(nr*nc)
	ev      <- svd.pca.obj$E/nr
	d.seq   <- svd.pca.obj$d.seq
	max.rk  <- length(w)

	exa.ev  <- c(sum(ev), ev)	
	C       <- min(nr, nc)
	d.max   <- ifelse(is.null(d.max),min(round(sqrt(C)), max.rk), d.max)
	if(d.max > max.rk)
		{
		warning(expression("The given maximun dimension 'd.max' is larger than the postive eingen values of the corresponding covariance matrix. We use maximum rank as 'd.max' "))
		d.max <- max.rk
		}

        ## number of eigenvalues greater than zero: 
        d.max.p1.gzero <- length(w[1:(d.max+1)][w[1:(d.max+1)]>0])
        d.max.p2.gzero <- length(w[1:(d.max+2)][w[1:(d.max+2)]>0])
        d.max.p3.gzero <- length(w[1:(d.max+3)][w[1:(d.max+3)]>0])
        ##
	ER = exa.ev[1:(d.max+1)]/exa.ev[2:(d.max+2)]
	d.opt.ER <- order(ER, na.last = FALSE)[d.max]
	GR = (log(w[1:d.max.p1.gzero]) - log(w[2:d.max.p2.gzero]))/(log(w[2:d.max.p2.gzero]) - log(w[3:d.max.p3.gzero]))
	d.opt.GR <- order(GR, na.last = FALSE)[d.max]

	result <- matrix(c(d.opt.ER, d.opt.GR, w[d.opt.ER+1], w[d.opt.GR+1]), 2, 2)
	Result <- data.frame(I(c("ER", "GR")),  result, rep(d.max, 2))
	colnames(Result) <- c("Criterion", "Optimal Dimension", "sd2", "d.ref.max")
	return(Result)
	}


RH.OptDim <- function(Obj, criteria = c("ER", "GR"), d.max = NULL){
	# what is Obj?
	if(class(Obj)=="svd.pca"|class(Obj)=="fsvd.pca") obj <- Obj
	else{
		if(class(Obj)=="pca.fit"|class(Obj)=="fpca.fit"){
			nr    <- Obj$data.dim[1]
			nc    <- Obj$data.dim[2]
			V.d   <- Obj$Sd2*(nr*nc)
			E     <- Obj$eigen.values*(nr*nc)
			d.seq <- seq.int(0, (length(V.d)-1))
			obj <- list(V.d = V.d, nr = nr, nc = nc, E = E)
			}	
		else{
			if(is.matrix(Obj)) obj <- svd.pca(Obj)
			else{
				if(!is.vector(Obj[[1]])|!is.numeric(Obj[[1]])
				  |!is.numeric(Obj[[2]])|length(Obj[[2]])!=2)
				  stop(c("'Obj' does not have the correct form."))
				else{# the function can deal with a list containing a vector of RSS for each d in the first listcomponent and the dimension as 2 d-vector 'c(nr, nc)' in the second component 
					nr  <- Obj[[2]][1]
					nc  <- Obj[[2]][2]
					V.d <- Obj[[1]][-length(Obj[[1]])]
					E   <- -diff(Obj[[1]]-Obj[[1]][1])
					d.seq = seq.int(0, (length(Obj[[1]])-1))
					obj <- list(V.d = V.d, nr = nr, nc = nc
						, E = E, d.seq = d.seq)
					}
				}
			}
		}

	result <- RH.dim.opt(obj, d.max = d.max)
	criteria <- match.arg(criteria, several.ok = TRUE)
	return(result[result[,1] %in% criteria, ])
	}

#####################################################################################################################
KSS.dim.opt <- function(obj, sig2.hat = NULL, alpha=0.01, factor.dim = NULL, d.max = NULL){
  # kann direkt mit fsvd.pca-objekten arbeiten
  nr       <- obj$nr
  nc       <- obj$nc
  spar.low <- obj$spar.low
  dat      <- obj$Q.orig
  dat.smth <- obj$Q.orig.smth
  w        <- obj$V.d/(nr*nc)  
  evec     <- obj$L
  Eval     <- obj$V.d
  max.rk   <- length(Eval)

 if(spar.low ==0) dat.smth <- NULL
  
### calculate traces
  
  I.smth1 <- smooth.Pspline(x=seq.int(0,1, length.out = nr), y = diag(rep(1,nr)), spar = spar.low, method = 1)$ysmth
  I.smth2 <- smooth.Pspline(x=seq.int(0,1, length.out = nr), y = I.smth1,         spar = spar.low, method = 1)$ysmth
  tr.dim.zero     <- sum(diag(I.smth2))
  tr.dim.zero.sqr <- sum(diag(I.smth2)^2)

  P             <- diag(1, nr) - tcrossprod(evec)
  pca.fit.p.obj <- eigen(P)
  W             <- pca.fit.p.obj[[2]]
  ## P.E is left out in further computations, since P.E[1:T]==rep(1,T)
  P.E           <- pca.fit.p.obj[[1]]

  W.smth  <- smooth.Pspline(x=seq.int(0,1, length.out = nr), y = W,               spar = spar.low, method = 1)$ysmth
  I.smth1 <- smooth.Pspline(x=seq.int(0,1, length.out = nr), y = diag(rep(1,nr)), spar = spar.low, method = 1)$ysmth
  I.smth2 <- smooth.Pspline(x=seq.int(0,1, length.out = nr), y = I.smth1,         spar = spar.low, method = 1)$ysmth
  tr.dim.zero     <- sum(diag(I.smth2))
  tr.dim.zero.sqr <- sum(diag(I.smth2)^2)

  diag.Wsmt <- diag(crossprod(W.smth)) # t(x) %*% y

  tr1 <- c(tr.dim.zero,     (sum(diag.Wsmt)   - cumsum(diag.Wsmt)))
  tr2 <- c(tr.dim.zero.sqr, (sum(diag.Wsmt^2) - cumsum(diag.Wsmt^2)))


### determine / calculate sig2.hat
  d.max0 = d.max
  sig2.hat0 = sig2.hat

  if(is.null(sig2.hat)){
	# estimation of sig2.hat: Classical, if one wants to use the KSS-Criterion for non-smoothed dat
	if(is.null(dat.smth)| !is.null(d.max)){
  		if(is.null(d.max)) d.max <- round(sqrt(min(nr, nc)))
		sig2.hat <- w[d.max+1]*(nc*nr)/(nr*nc - (nr + nc)*d.max - 1)
	}
	# estimation of sig2.hat: Variance-Estimator see Section 3.4 (KSS) 
	else{
		tr		<- (nr + sum(diag(I.smth2)) - 2*sum(diag(I.smth1)))
		sig2.hat	<- sum((dat-dat.smth)^2)/((nc-1)*tr)
                
#	rice np.variance estimator
#		sig2.hat	<- sum(diff(dat)^2)/(2*(nr - 1)*nc)	
	}

  }

## calculate the criteria of all dimensions:
  delta       <- (Eval - (nc-1) * sig2.hat * tr1[1:max.rk])/(sig2.hat * sqrt(2*nc*tr2[1:max.rk]))
  thres1      <- qnorm(1-alpha)
  thres2      <- sqrt(2*log(min(nr, nc)))# default alpha = NULL / falls alpha != 0 dann, werden beide beide berechnet
  level2      <- 1 - pnorm(thres2)
  crit1       <- delta - thres1
  crit2       <- delta - thres2
  d.opt.KSS1  <- length(crit1[crit1 > 0])# minus 1, weil start bei dim = 0
                                         # plus  1, weil nur die dim, die das crit nicht erfüllen.
  d.opt.KSS2  <- length(crit2[crit2 > 0])
  if(is.null(dat.smth)& is.null(d.max0) & is.null(sig2.hat0)){
	sig2.hat <- w[d.opt.KSS1+1]*(nc*nr)/(nr*nc - (nr + nc)*d.opt.KSS1 - 1)
  	delta       <- (Eval - (nc-1) * sig2.hat * tr1[1:max.rk])/(sig2.hat * sqrt(2*nc*tr2[1:max.rk]))
 	thres1      <- qnorm(1-alpha)
 	thres2      <- sqrt(2*log(min(nr, nc)))# default alpha = NULL / falls alpha != 0 dann, werden beide beide berechnet
  	level2      <- 1 - pnorm(thres2)
  	crit1       <- delta - thres1
 	crit2       <- delta - thres2
  	d.opt.KSS1  <- length(crit1[crit1 > 0])
  }

  if(!is.null(factor.dim)){
    used.dim.C1 <- factor.dim
    used.dim.C2 <- factor.dim
  }else{
    used.dim.C1 <- d.opt.KSS1
    used.dim.C2 <- d.opt.KSS2
  }
  result1     <- c(d.opt.KSS1, used.dim.C1, w[(d.opt.KSS1+1)], sig2.hat, alpha )
  result2     <- c(d.opt.KSS2, used.dim.C2, w[(d.opt.KSS2+1)], sig2.hat, level2)# level2: p.value of thres2 (thres2: alternativ crit.value)

  result      <- rbind(result1, result2)
  Result      <- vector("list", 2)
  Result[[1]] <- data.frame(I(c("KSS.C", "KSS.C2")), result)
  colnames(Result[[1]]) <- c("Criterion", "Optimal Dimension", "Used Dimension", "sd2.rest", "sd2.hat", "level")
  rownames(Result[[1]]) <- c("KSS.1", "KSS.2")

  Result[[2]] <- list(Test.Stat  = round(delta[1],2), p.value = round(1-pnorm(delta[1]), 2),
                      crit.value = round(thres1, 2),  sig.level = round(alpha, 2))
return(Result)
}

## KSS.OptDim() =====================================================================================

KSS.OptDim <- function(Obj,
                       criteria    = c("KSS.C", "KSS.C2"),
                       sig2.hat    = NULL,
                       alpha       = 0.01, 
                       d.max       = NULL,
                       factor.dim  = NULL, 
			     spar = NULL){
  ## what is Obj?
  if(class(Obj)=="svd.pca"|class(Obj)=="fsvd.pca"){
    if(class(Obj)=="fsvd.pca") obj <- Obj
    else{
      ## Liste um spar.low und Q.orig.smth erweitern:
      nr          <- Obj$nr
      nc          <- Obj$nc
      spar.low    <- ifelse(is.null(spar), 0, spar)          
      Q.orig      <- Obj$Q.orig
      Q.orig.smth <- NULL
      L           <- Obj$L
      V.d         <- Obj$V.d
      
      obj <- list(nr = nr, nc = nc, spar.low = spar.low, Q.orig = Q.orig,
                  Q.orig.smth = Q.orig.smth, L = L, V.d = V.d)
    }    
  }
  if(class(Obj)=="pca.fit"|class(Obj)=="fpca.fit"){
    if(class(Obj)=="fpca.fit"){
      ## umbenennungen zu fsvd.pca-Elementen
      nr          <- Obj$data.dim[1]
      nc          <- Obj$data.dim[2]
      spar.low    <- Obj$spar.low
      Q.orig      <- Obj$orig.values
      Q.orig.smth <- Obj$orig.values.smth
      L           <- Obj$L
      V.d         <- Obj$Sd2*(nr*nc)
      
      obj <- list(nr = nr, nc = nc, spar.low = spar.low, Q.orig = Q.orig,
                  Q.orig.smth = Q.orig.smth, L = L, V.d = V.d)
    }else{
      nr          <- Obj$data.dim[1]
      nc          <- Obj$data.dim[2]
      spar.low    <- ifelse(is.null(spar), 0, spar)   
      Q.orig      <- Obj$orig.values
      Q.orig.smth <- NULL
      L           <- Obj$L
      V.d         <- Obj$Sd2*(nr*nc)
      
      obj <- list(nr = nr, nc = nc, spar.low = spar.low, Q.orig = Q.orig,
                  Q.orig.smth = Q.orig.smth, L = L, V.d = V.d)
    } 
  }else{
    if(is.matrix(Obj)){
      ## fPCA
      Obj        <- fpca.fit(dat = Obj, given.d = factor.dim, spar = spar)
      ## umbenennungen zu fsvd.pca-Elementen
      nr          <- Obj$data.dim[1]
      nc          <- Obj$data.dim[2]
      spar.low    <- Obj$spar.low
      Q.orig      <- Obj$orig.values
      Q.orig.smth <- Obj$orig.values.smth
      L           <- Obj$L
      V.d         <- Obj$Sd2*(nr*nc)
      
      obj <- list(nr = nr, nc = nc, spar.low = spar.low, Q.orig = Q.orig,
                  Q.orig.smth = Q.orig.smth, L = L, V.d = V.d)
      
    }
  }
  result        <- KSS.dim.opt(obj, sig2.hat = sig2.hat, alpha = alpha, factor.dim = factor.dim, d.max = d.max)
  criteria      <- match.arg(criteria, several.ok = TRUE)
  Result        <- vector("list", 2)
  Result[[1]]   <- result[[1]][result[[1]][,1] %in% criteria, ]
  Result[[2]]   <- result[[2]]
  return(Result)
}
##===========================================================================================================



############## Est dim intern  ######################################################################################


EstDim <- function(Obj, 
                   dim.criterion = c("PC1", "PC2", "PC3",
                     "IC1", "IC2", "IC3",
                     "IPC1","IPC2", "IPC3",
                     "KSS.C",
                     "ED",  "ER",  "GR"),
                   d.max,
                   factor.dim,
                   sig2.hat,
			 spar = NULL,
                   level = 0.01
                   )
  {
  ## missing parameters
       if(missing(factor.dim)) factor.dim  <- NULL
  	 if(missing(d.max))      d.max       <- NULL
 	 if(missing(sig2.hat))   sig2.hat    <- NULL
	 if(missing(spar))       spar        <- NULL

  ## estimation 
	 criteria <- match.arg(dim.criterion)
         
 	 est.dim       <- switch(criteria,
                          PC1  = try(B.OptDim(Obj, criteria = c("PC1")
                            , d.max = d.max, sig2.hat = sig2.hat)),
                          PC2  = try(B.OptDim(Obj, criteria = c("PC2")
                            , d.max = d.max, sig2.hat = sig2.hat)),
                          PC3  = try(B.OptDim(Obj, criteria = c("PC3")
                            , d.max = d.max, sig2.hat = sig2.hat)),
                          IC1  = try(B.OptDim(Obj, criteria = c("IC1")
                            , d.max = d.max, sig2.hat = sig2.hat)),
                          IC2  = try(B.OptDim(Obj, criteria = c("IC2")
                            , d.max = d.max, sig2.hat = sig2.hat)),
                          IC3  = try(B.OptDim(Obj, criteria = c("IC3")
                            , d.max = d.max, sig2.hat = sig2.hat)),
                          IPC1 = try(B.OptDim(Obj, criteria = c("IPC1")
                            , d.max = d.max, sig2.hat = sig2.hat)),
                          IPC2 = try(B.OptDim(Obj, criteria = c("IPC2")
                            , d.max = d.max, sig2.hat = sig2.hat)),
                          IPC3 = try(B.OptDim(Obj, criteria = c("IPC3")
                            , d.max = d.max, sig2.hat = sig2.hat)),
                          ED   = try(O.OptDim(Obj, d.max = d.max)),				
                          ER   = try(RH.OptDim(Obj, criteria = c("ER")
                            , d.max = d.max)),
                          GR   = try(RH.OptDim(Obj, criteria = c("GR")
                            , d.max = d.max)),
                          KSS.C  = try(KSS.OptDim(Obj, criteria = c("KSS.C")
				    , sig2.hat = sig2.hat, alpha=level, spar = spar
                                    , factor.dim=factor.dim, d.max=d.max)[[1]])
                          )
	est.dim
       }

############## Est dim compair  ######################################################################################



OptDim.default <- function(Obj, 
                   criteria = c("PC1", "PC2", "PC3",
                     "IC1", "IC2", "IC3",
                     "IPC1","IPC2", "IPC3",
                     "KSS.C",
                     "ED",  "ER",  "GR"),
			standardize = FALSE,
			d.max,
			sig2.hat,
			spar,
			level = 0.01
                   )
  {
  ## if Eup or KSS classes
	if(class(Obj)=="KSS"| class(Obj)== "Eup") Obj <- Obj$unob.fact.stru + Obj$residuals
	else is.regular.panel(Obj, stopper = TRUE)

  ## missing parameters
  	 if(missing(d.max))      d.max       <- NULL
 	 if(missing(sig2.hat))   sig2.hat    <- NULL
	 if(missing(spar))       spar        <- NULL

  ## standardize
	 if(standardize) Obj <- standardize(Obj)

  ## estimation 
	 criteria <- match.arg(criteria, several.ok = TRUE)
         
 	 est.dim       <- function(criteria, d.max, sig2.hat, level){
				switch(criteria,
                          PC1  = try(B.OptDim(Obj, criteria = c("PC1")
                            , d.max = d.max, sig2.hat = sig2.hat)),
                          PC2  = try(B.OptDim(Obj, criteria = c("PC2")
                            , d.max = d.max, sig2.hat = sig2.hat)),
                          PC3  = try(B.OptDim(Obj, criteria = c("PC3")
                            , d.max = d.max, sig2.hat = sig2.hat)),
                          IC1  = try(B.OptDim(Obj, criteria = c("IC1")
                            , d.max = d.max, sig2.hat = sig2.hat)),
                          IC2  = try(B.OptDim(Obj, criteria = c("IC2")
                            , d.max = d.max, sig2.hat = sig2.hat)),
                          IC3  = try(B.OptDim(Obj, criteria = c("IC3")
                            , d.max = d.max, sig2.hat = sig2.hat)),
                          IPC1 = try(B.OptDim(Obj, criteria = c("IPC1")
                            , d.max = d.max, sig2.hat = sig2.hat)),
                          IPC2 = try(B.OptDim(Obj, criteria = c("IPC2")
                            , d.max = d.max, sig2.hat = sig2.hat)),
                          IPC3 = try(B.OptDim(Obj, criteria = c("IPC3")
                            , d.max = d.max, sig2.hat = sig2.hat)),
                          ED   = try(O.OptDim(Obj, d.max = d.max)),				
                          ER   = try(RH.OptDim(Obj, criteria = c("ER")
                            , d.max = d.max)),
                          GR   = try(RH.OptDim(Obj, criteria = c("GR")
                            , d.max = d.max)),
                          KSS.C  = try(KSS.OptDim(Obj, criteria = c("KSS.C")
				    , sig2.hat = sig2.hat, alpha=level
                             , factor.dim= NULL, d.max=d.max, spar = spar)[[1]])
                          )
	}

	Result <- sapply(criteria, est.dim, d.max = d.max, sig2.hat = sig2.hat, level = level, simplify = FALSE)
	summary <- sapply(1:length(criteria), function(l) Result[[l]][,1:2])
	Summary <- matrix(as.numeric(summary[2,]), ncol(summary))
	colnames(Summary) <- " "
	rownames(Summary) <- summary[1, ]
	Summary  <- t(Summary )
	Result$summary <- Summary 
	Result$criteria <- criteria
	Result$BaiNgC <- criteria%in%c("PC1", "PC2", "PC3", "IC1", "IC2", "IC3")
	Result$BaiC   <- criteria%in%c("IPC1", "IPC2", "IPC3")
	Result$KSSC    <- criteria%in%c("KSS.C")
	Result$OnatC  <- criteria%in%c("ED")
	Result$RHC    <- criteria%in%c("ER", "GR")
	Result$cl     <- match.call()
	Result$obj    <- Obj
	structure(Result, class = "OptDim")
}


# ####################### Methods ########################
	OptDim <- function(Obj, 
                   criteria = c("PC1", "PC2", "PC3",
                     "IC1", "IC2", "IC3",
                     "IPC1","IPC2", "IPC3",
                     "KSS.C",
                     "ED",  "ER",  "GR"),
			standardize = FALSE,
			d.max,
			sig2.hat,
			spar,
			level = 0.01){ UseMethod("OptDim")}


## Print
print.OptDim <- function(x,...){
   	cat("Call: ")
	cl <- x$cl
   	print(cl)

	ll <- length(x$criteria)

	if(sum(x$KSSC) > 0){
	cat("\n---------")
	cat("\nCriterion of Kneip, Sickles, and Song (2009):\n\n")
		KSS2009 <- numeric(0)
		for(l in 1:ll) if(x$KSSC[l]) KSS2009 <- rbind(KSS2009, x[[l]][,1:2])
		dimKSS2009 <- matrix(KSS2009[, 2], 1, length(KSS2009[, 2]))
		colnames(dimKSS2009) <- "KSS.C"
		rownames(dimKSS2009) <- " "
		print(dimKSS2009)
	}
	
	if(sum(x$BaiNgC) > 0){
	cat("\n---------")
	if(sum(x$BaiNgC) > 1) cat("\nCriteria of Bai and Ng (2002):\n\n")
	else  cat("\nCriterion of Bai and Ng (2002):\n\n")
		BN2002 <- numeric(0)
		for(l in 1:ll) if(x$BaiNgC[l]) BN2002 <- rbind(BN2002, x[[l]][,1:2])
		dimBN2002 <- matrix(BN2002[, 2], 1, length(BN2002[, 2]))
		colnames(dimBN2002) <- BN2002[, 1]
		rownames(dimBN2002) <- " "
		print(dimBN2002)
	}
	if(sum(x$RHC) > 0){
	cat("\n--------")
	if(sum(x$RHC) > 1) cat("\nCriteria of Ahn and Horenstein (2008):\n\n")
	else cat("\nCriterion of Ahn and Horenstein (2008):\n\n") 
		AH2008 <- numeric(0)
		for(l in 1:ll) if(x$RHC[l]) AH2008 <- rbind(AH2008, x[[l]][,1:2])
		dimAH2008 <- matrix(AH2008[, 2], 1, length(AH2008[, 2]))
		colnames(dimAH2008) <- AH2008[, 1]
		rownames(dimAH2008) <- " "
		print(dimAH2008)
	}
	if(sum(x$BaiC) > 0){
	cat("\n---------")
	if(sum(x$BaiC) > 1) cat("\nCriteria of Bai (2004):\n\n") 
	else cat("\nCriterion of Bai (2004):\n\n")
		B2004 <- numeric(0)
		for(l in 1:ll) if(x$BaiC[l]) B2004 <- rbind(B2004, x[[l]][,1:2])
		dimB2004 <- matrix(B2004[, 2], 1, length(B2004[, 2]))
		colnames(dimB2004) <- B2004[, 1]
		rownames(dimB2004) <- " "
		print(dimB2004)
	}
	if(sum(x$OnatC) > 0){
	cat("\n---------")
	cat("\nCriterion of Onatski (2009):\n\n")
		O2009 <- numeric(0)
		for(l in 1:ll) if(x$OnatC[l]) O2009 <- rbind(O2009, x[[l]][,1:2])
		dimO2009 <- matrix(O2009[, 2], 1, length(O2009[, 2]))
		colnames(dimO2009) <- O2009[, 1]
		rownames(dimO2009) <- " "
		print(dimO2009)
	}

 }
	
## Print
summary.OptDim <- function(object,...){
   	cat("Call: ")
	cl <- object$cl
   	print(cl)

	ll <- length(object$criteria)

	if(sum(object$KSSC) > 0){
	cat("\n---------")
	cat("\nSequential testing of Kneip, Sickles, and Song (2009):\n\n")
		KSS2009 <- numeric(0)
		for(l in 1:ll) if(object$KSSC[l]) KSS2009 <- rbind(KSS2009, object[[l]])
		dimKSS2009 <- KSS2009[, c(2, 6)]
		colnames(dimKSS2009) <- c("Estimate", "Level")
		rownames(dimKSS2009) <- " KSS.C"
		print(round(dimKSS2009,3))
		cat("\n\nUsed Std Err: ")
		cat(round(KSS2009[1, 5], 3))
	}
	
	if(sum(object$BaiNgC) > 0){
	cat("\n\n---------")
	if(sum(object$BaiNgC) > 1) cat("\nCriteria of Bai and Ng (2002):\n\n")
	else  cat("\nCriterion of Bai and Ng (2002):\n\n")
		BN2002 <- numeric(0)
		for(l in 1:ll) if(object$BaiNgC[l]) BN2002 <- rbind(BN2002, object[[l]])
		dimBN2002 <- BN2002[, 2:3]
		colnames(dimBN2002) <- c("Estimate", "Std Err")
		rownames(dimBN2002) <- BN2002[, 1]
		print(round(dimBN2002, 3))
		cat("\n\nUsed d.max:")
		cat(BN2002[1, 5])
	}
	if(sum(object$RHC) > 0){
	cat("\n\n--------")
	if(sum(object$RHC) > 1) cat("\nCriteria of Ahn and Horenstein (2008):\n\n")
	else cat("\nCriterion of Ahn and Horenstein (2008):\n\n") 
		AH2008 <- numeric(0)
		for(l in 1:ll) if(object$RHC[l]) AH2008 <- rbind(AH2008, object[[l]])
		dimAH2008 <- AH2008[, 2:3]
		colnames(dimAH2008) <-  c("Estimate", "Std Err")
		rownames(dimAH2008) <- AH2008[,1]
		print(round(dimAH2008, 3))
		cat("\n\nUsed d.max:")
		cat(AH2008[1, 4])	
	}

	if(sum(object$BaiC) > 0){
	cat("\n\n---------")
	if(sum(object$BaiC) > 1) cat("\nCriteria of Bai (2004):\n\n") 
	else cat("\n\nCriterion of Bai (2004):\n\n")
		B2004 <- numeric(0)
		for(l in 1:ll) if(object$BaiC[l]) B2004 <- rbind(B2004, object[[l]])
		dimB2004 <- B2004[, 2:3]
		colnames(dimB2004) <- c("Estimate", "Std Err")
		rownames(dimB2004) <- B2004[, 1]
		print(round(dimB2004, 3))
		cat("\n\nUsed d.max:")
		cat(B2004[1, 5])
	}

	if(sum(object$OnatC) > 0){
	cat("\n\n---------")
	cat("\nCriterion of Onatski (2009):\n\n")
		O2009 <- numeric(0)
		for(l in 1:ll) if(object$OnatC[l]) O2009 <- rbind(O2009, object[[l]])
		dimO2009 <- O2009[, 2:3]
		colnames(dimO2009) <- c("Estimate", "Std Err")
		rownames(dimO2009) <- O2009[1, 1]
		print(round(dimO2009, 3))
	}


 }

plot.OptDim <- function(x, main, border, col, ...){
	Resultcrit <- x$summary

	dims <- as.numeric(levels(as.factor(Resultcrit)))
	dcol <- rainbow(length(dims))

	Obj <- x$obj
	nr  <- nrow(Obj)
	nc  <- ncol(Obj)
	dual = FALSE
	if(nr > nc) dual = TRUE
	if(dual)  Q <- t(Obj)%*%Obj
	else Q <- Obj%*%t(Obj)
	eig <- eigen(Q, only.values = TRUE)[[1]]
	d.max = min(nc, nr)
	perceig <- eig/sum(eig)
	perccum <- cumsum(perceig)
	ycords <- c(((3/2)*perceig[1]-0.5*perceig[2]), perceig[-d.max])
	yycoreds <- ycords 
	yycoreds[1] <- (ycords[1] - perceig[1])/2 + perceig[1]


	critperdim <- sapply(dims, function(x) paste(colnames(Resultcrit)[Resultcrit[1,]== x], collapse = ",    "))


	if(missing(main)) main = "Screeplot"
	if(missing(border)) border = FALSE
	if(missing(col)) col = "lightslategrey"
	if(dims[1]==0 && perceig[2] > 0) Ylim <- (3/2)*perceig[1]-0.5*perceig[2]
	else Ylim <- perceig[1]

	xaxis <- barplot(perceig,  main = main, border = border 
		,axes = FALSE, col = col, ylim = c(0, Ylim)
		, ylab = "Proportion of variance"
		, xlab = "Ordered eigenvalues")
	
	axis(1, xaxis[1:max(dims)], 1:max(dims), tick = FALSE)
	axis(2, perceig[1:max(dims)], paste((round(perceig[1:max(dims)], 3)*100), "%"))
	axis(4, yycoreds[(dims+1)], dims)

	for(r in 1:length(dims)){
	if(dims[r]!=0) abline(h = yycoreds[(dims[r]+1)], lty = 3, col = "gray")
	}
	text(rep(d.max, length(dims)), (yycoreds[(dims+1)]+ 0.01) , critperdim, adj = 1, cex = 1.1)
	}



