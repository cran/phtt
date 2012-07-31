\name{KSS}
\alias{KSS}
\alias{KSS.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{KSS-Routine}
\description{Estimation of Panel Data Models with Heterogeneous Time Trends}
\usage{
KSS(formula,
    additive.effects = c("none", "individual", "time", "twoways"),
    consult.dim.crit = FALSE,
    d.max            = NULL,
    sig2.hat         = NULL,
    factor.dim       = NULL,
    level            = 0.01,
    spar             = NULL,
    CV               = FALSE,
    convergence      = 1e-6,
    restrict.mode    = c("restrict.factors","restrict.loadings"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{An object of class 'formula'.}
  \item{additive.effects}{
    Type of Data Transformations:
    \itemize{
      \item \code{"none"}: for no transformation
      \item \code{"individual"}: for within transformation
      \item \code{"time"}: for between transformation
      \item \code{"twoways"}: for twoways transformation
  }}
  \item{consult.dim.crit}{logical.
    \itemize{
      \item{If \code{consult.dim.crit} is \code{FALSE} (default) and
	\code{factor.dim} is \code{NULL}:
	Only the dimensionality criterion of Kneip, Sickles & Song 2012
	is used.}
      \item{If \code{consult.dim.crit} is \code{TRUE} and
	\code{factor.dim} is \code{NULL}:
	All implemented dimensionality criteria as implemented in the
	function \code{OptDim()} are computed and the user
	has to select one proposed dimension via a GUI.}
    }}
  \item{d.max}{A maximal dimension needed for some
    dimensionality-criteria that are implemented in the function
    \code{OptDim()}. The default (\code{d.max=NULL}) yields to an internal
    selection of \code{d.max}.}
  \item{sig2.hat}{Standard deviation of the error-term. The default
    (\code{sig2.hat=NULL}) yields to an internal
    estimation of \code{sig2.hat}.}
  \item{factor.dim}{Dimension of Factor-Structure. The default
    (\code{factor.dim=NULL}) yields to an internal
    estimation of \code{factor.dim}.}
  \item{level}{Significance-level for Dimensionality-Criterion of Kneip,
    Sickles & Song 2012.}
  \item{spar}{Smoothing parameter for spline smoothing of the
    residuals. If (\code{spar=NULL}) (default) and \code{CV=FALSE}
    \code{spar} is determined via generalized cross validation (GCV).}
  \item{CV}{logical. Selects the procedure for the determination of the
    smoothing parameter \code{spar}.
    \itemize{
      \item{If \code{CV=FALSE} (default) and \code{spar=NULL}:
	The smoothing parameter \code{spar} is determined by GCV.}
      \item{If \code{CV=TRUE} and \code{spar=NULL}:
      The smoothing parameter \code{spar} is determined by Leave-one-out
      cross validation (CV).}
    }}
  \item{convergence}{Convergence criterion for the CV-optimization of
    the smoothing parameter \code{spar}. Default is
    \code{convergence=1e-6}.}
  \item{restrict.mode}{
    Type of Restriction on the Factor-Structure:
    \itemize{
      \item \code{"restrict.factors"}: Factors are restricted to have an
      euclidean norm of 1.
      \item \code{"restrict.loadings"}: Factor-Loadings are restricted
      to have an euclidean norm of 1.
    }}
  \item{...}{Additional arguments to be passed to the low level functions.}
}
\details{
  'KSS' is a function to estimate panel data models with unobserved
  heterogeneous time trends v_i(t). The considered model in Kneip,
  Sickles & Song (2012) is given by
  \eqn{$Y_{it}=\theta_{t}+\sum_{j=1}^P\beta_{j}
    X_{itj}+v_i(t)+\epsilon_{it}\quad i=1,\dots,n; t=1,\dots,T.$}{
        Y_{it}=\theta_{t}+\sum_{j=1}^P\beta_{j}
	X_{itj}+v_i(t)+\epsilon_{it}, i=1,...,n; t=1,...,T.}
      Where the individual time trends, v_i(t), are assumed to come from
      a finite dimensional factor model
      \eqn{$v_i(t)=\sum_{l=1}^d\lambda_{il}f_l(t),\quad\lambda_{il}\in
	R,\quad f_l\in L^2[0,T].$}{v_i(t)=\sum_{l=1}^d\lambda_{il}f_l(t),
	 \lambda_{il}\in R,  f_l\in L^2[0,T].}
      The unobserved functions v_i(t) can be interpreted as smooth
      functions of a continuous argument t, as well as stochastic
      processes for discrete argument t.
      \itemize{
	\item{\code{formula} }{Usual 'formula'-object. If you wish to
	  estimate a model without an intercept use '-1' in the
	  formula-specification. Each Variable has to be given as a
	  TxN-matrix. Missing values are not allowed.}
	\item{\code{additive.effects}}{
	  \itemize{
	    \item \code{"none"}: The data is not transformed, except for an
	    eventually subtraction of the overall mean; if the model is
	    estimated with an intercept. The assumed model can be written as
	    \eqn{$Y_{it}=\mu+\sum_{j=1}^P\beta_{j}
	      X_{itj}+v_i(t)+\epsilon_{it}\quad i=1,\dots,n; t=1,\dots,T.$}{
              Y_{it}=\mu+\sum_{j=1}^P\beta_{j}
	      X_{itj}+v_i(t)+\epsilon_{it}, i=1,...,n; t=1,...,T.}
	    The parameter 'mu' is set to zero if '-1' is used in \code{formula}.
	    \item \code{"individual"}: This is the "within"-model, which 
	    assumes that there are time-constant individual effects,
	    tau_i, besides the individual time trends v_i(t). The
	    model can be written as
	    \eqn{$Y_{it}=\mu+\sum_{j=1}^P\beta_{j}
	      X_{itj}+v_i(t)+\alpha_{i}+\epsilon_{it}\quad i=1,\dots,n; t=1,\dots,T.$}{
              Y_{it}=\mu+\sum_{j=1}^P\beta_{j}
	      X_{itj}+v_i(t)+\alpha_{i}+\epsilon_{it}, i=1,...,n; t=1,...,T.}
	    The parameter 'mu' is set to zero if '-1' is used in \code{formula}.
	    \item \code{"time"}: This is the "between"-model, which assumes
	    that there is a common (for all individuals) time trend, beta_0(t). The
	    model can be written as
	    \eqn{$Y_{it}=\mu+\theta_{t}+\sum_{j=1}^P\beta_{j}
	      X_{itj}+v_i(t)+\epsilon_{it}\quad i=1,\dots,n; t=1,\dots,T.$}{
              Y_{it}=\mu+\theta_{t}+\sum_{j=1}^P\beta_{j}
	      X_{itj}+v_i(t)+\epsilon_{it}, i=1,...,n; t=1,...,T.}
	    The parameter 'mu' is set to zero if '-1' is used in \code{formula}.
	    \item \code{"twoways"}: This is the "twoways"-model ("within" &
	    "between"), which assumes that there are time-constant
	    individual effects, tau_i, and a common time trend,
	    beta_0(t). The model can be written as
	    \eqn{$Y_{it}=\mu+\theta_{t}+\sum_{j=1}^P\beta_{j}
	      X_{itj}+\alpha_i+v_i(t)+\epsilon_{it}\quad i=1,\dots,n; t=1,\dots,T.$}{
              Y_{it}=\mu+\theta_{t}+\sum_{j=1}^P\beta_{j}
	      X_{itj}+\alpha_i+v_i(t)+\epsilon_{it}, i=1,...,n; t=1,...,T.}
	    The parameter 'mu' is set to zero if '-1' is used in \code{formula}.
      }}}
    }
    \value{
      'KSS' returns an object of 'class' '"KSS"'.

      An object of class '"KSS"' is a list containing at least the
      following components:
      \itemize{
	\item{\code{dat.matrix}: }{Whole data set stored within a
	  (N*T)x(p+1)-Matrix, where P is the number of independent
	  variables without the intercept.}
	\item{\code{dat.dim}: }{Vector of length 3: c(T,N,p)}
	\item{\code{slope.para}: }{Beta-parameters}
	\item{\code{beta.V}: }{Covariance matrix of the beta-parameters.}
	\item{\code{names}: }{Names of the dependent and independent variables.}
	\item{\code{is.intercept}: }{Used an intercept in the formula?: TRUE or FALSE}
	\item{\code{additive.effects}: }{Additive effect type. One of: "none","individual","time", "twoways".}
	\item{\code{Intercept}: }{Intercept-parameter}
	\item{\code{Add.Ind.Eff}: }{Estimated values of additive individual effects.}
	\item{\code{Add.Tim.Eff}: }{Estimated values of additive time effects.}
	\item{\code{unob.factors}: }{Txd-matrix of estimated unobserved common
	  factors, where 'd' is the number of used factors.} 
	\item{\code{ind.loadings}: }{Nxd-matrix of loadings parameters.}
	\item{\code{unob.fact.stru}: }{TxN-matrix of the estimated factor
	  structure. Each column represents an estimated individual unobserved time trend.}
	\item{\code{used.dim}: }{Used dimensionality of the factor structure.}
	\item{\code{optimal.dim}: }{List of proposed dimensionalities.}
	\item{\code{fitted.values}: }{Fitted values.}
	\item{\code{orig.Y}: }{Original values of the dependent variable.}       
	\item{\code{residuals}: }{Residuals}    
	\item{\code{sig2.hat}: }{Estimated variance of the error term.}     
	\item{\code{degrees.of.freedom}: }{Degrees of freedom of the residuals.}
	\item{\code{call}}
      }
}
\references{
  \itemize{
    \item Kneip, A., Sickles, R. C., Song, W., 2012 \dQuote{A New Panel Data
      Treatment for Heterogeneity in Time Trends}, \emph{Econometric Theory}
}}
\author{
  Dominik Liebl
}
% \note{
% %%  ~~further notes~~
% }
% %% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
Eup
}
\examples{
## See the example in 'help(Cigar)' in order to take a look at the
## data set Cigar

##########
## DATA ##
##########

data(Cigar)
## Panel-Dimensions:
N <- 46
T <- 30
## Dependent variable:
  ## Cigarette-Sales per Capita
  l.Consumption    <- log(matrix(Cigar$sales, T,N))
## Independent variables:
  ## Consumer Price Index
  cpi        <- matrix(Cigar$cpi, T,N)
  ## Real Price per Pack of Cigarettes 
  l.Price  <- log(matrix(Cigar$price, T,N)/cpi)
  ## Real Disposable Income per Capita  
  l.Income    <- log(matrix(Cigar$ndi,   T,N)/cpi)

## Estimation:
KSS.fit      <- KSS(l.Consumption~l.Price+l.Income, CV=TRUE)
(KSS.fit.sum <- summary(KSS.fit))
plot(KSS.fit.sum)
}

