\name{Eup}
\alias{Eup}
\alias{Eup.default}
\title{Eup-Routine}
\description{Estimation of Panel Data Models with Interactive Fixed Effects.}
\usage{
Eup(formula,
    additive.effects = c("none", "individual", "time", "twoways"),
    dim.criterion    = c("PC1", "PC2", "PC3", "BIC3","IC1", "IC2" , "IC3",
                         "IPC1", "IPC2", "IPC3"),
    d.max            = NULL,
    sig2.hat         = NULL,
    factor.dim       = NULL,
    double.iteration = TRUE,
    start.beta       = NULL,
    max.iteration    = 500,
    convergence      = 1e-6,
    restrict.mode    = c("restrict.factors", "restrict.loadings"),
    ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{An object of class 'formula' where the arguments are matrices. The number of rows has to be equal to the temporal dimension and the number of columns has to be equal to the number of individuals. The details of model specification are given under 'Details'.}
  \item{additive.effects}{
    Type of Data Transformations:
    \itemize{
      \item \code{"none"}: for no transformation
      \item \code{"individual"}: for within transformation
      \item \code{"time"}: for between transformation
      \item \code{"twoways"}: for twoways transformation
  }}
  \item{dim.criterion}{The dimensionality criterion to be used if \code{factor.dim} is left unspecified. The default criterion is \code{"PC1"}}.
  \item{d.max}{Maximal dimension used in the dimensionality-criteria of Bai
    (2009). The default (\code{d.max=NULL}) yields to an internal selection of
    d.max.}
  \item{sig2.hat}{The squared standard deviation of the error-term required for the computation of some dimensionality criteria. The user can specify it in instead of \code{d.max}. The default (\code{sig2.hat=NULL}) yields to an internal estimation.}
  \item{factor.dim}{Dimension of Factor-Structure, pre-specified
    by the user. The default (\code{factor.dim=NULL}) yields to an internal estimation.}
  \item{double.iteration}{logical. If \code{FALSE} the update of the factor dimension \code{d} will be done simultaneously with remaining model parameters without alternating between inner and outer iteration. This may speed up computations, but the convergence is less stable than in the default setting.}
  \item{start.beta}{allows the user to give a vector of starting values for the slope parameters.}
  \item{max.iteration}{controls the maximum number of iterations. The default is '500'.}
  \item{convergence}{Convergence condition of the estimators. The default is '1e-6'.}
  \item{restrict.mode}{
    Type of Restriction on the Factor-Structure:
    \itemize{
      \item \code{"restrict.factors"}: Factors are restricted to have an
      euclidean norm of 1.
      \item \code{"restrict.loadings"}: Factor-Loadings are restricted
      to have an euclidean norm of 1.
    }
  }
  \item{...}{Additional arguments to be passed to the low level functions.}
}
\details{
'Eup' is a function to estimate equidistant panel data models with unobserved
  multiple time varying individual effects. The considered model is given by
  \eqn{$Y_{it}= \sum_{j=1}^P\beta_{j}
    X_{itj}+ v_{it} + \epsilon_{it}\quad i=1,\dots,n; t=1,\dots,T.$}{
        Y_{it}= \sum_{j=1}^P\beta_{j}X_{itj}+v_{it}+ \epsilon_{it}, i=1,...,n; t=1,...,T.}
  Where the individual time trends, v_{it}, are assumed to come from a finite dimensional factor model: 
   \eqn{$v_{it}=\sum_{l=1}^d\lambda_{il}f_{lt},\quad \lambda_{il}, f_{lt} \in	R.$}{v_{it}=\sum_{l=1}^d\lambda_{il}f_{lt}.}
         \itemize{
	\item{\code{formula} }{Usual 'formula'-object. If you wish to
	  estimate a model without an intercept use '-1' in the
	  formula-specification. Each Variable has to be given as a
	  TxN-matrix. Missing values are not allowed.}
	\item{\code{additive.effects}}{
	  \itemize{
	    \item \code{"none"}: The data is not transformed, except for a subtraction of the overall mean, if the model is estimated with an intercept. The assumed model can be written as
	    \eqn{$Y_{it}= \mu + \sum_{j=1}^P\beta_{j}
	      X_{itj}+v_{it}+\epsilon_{it}\quad i=1,\dots,n; t=1,\dots,T.$}{
              Y_{it}=\mu+\sum_{j=1}^P\beta_{j}
	      X_{itj}+v_{it}+\epsilon_{it}, i=1,...,n; t=1,...,T.}
	    The parameter 'mu' is set to zero if '-1' is used in \code{formula}.
	    \item \code{"individual"}: This is the "within"-model, which 
	    assumes that there are time-constant individual effects,
	    alpha_i, besides the individual time trends v_{it}. The
	    model can be written as
	    \eqn{$Y_{it}= \mu +\alpha_{i} + \sum_{j=1}^P\beta_{j}
	      X_{itj}+v _{it} + \epsilon_{it} \quad i=1,\dots,n; t=1,\dots,T.$}{
              Y_{it}=\mu + \alpha_{i}+ \sum_{j=1}^P\beta_{j}
	      X_{itj} + v_{it} + \epsilon_{it}, i=1,...,n; t=1,...,T.}
	    The parameter 'mu' is set to zero if '-1' is used in \code{formula}.
	    \item \code{"time"}: This is the "between"-model, which assumes
	    that there is a common time trend (for all individuals), theta_t. The
	    model can be written as
	    \eqn{$Y_{it}=\mu+ \theta_t + \sum_{j=1}^P\beta_{j}
	      X_{itj}+v_i(t)+\epsilon_{it}\quad i=1,\dots,n; t=1,\dots,T.$}{
              Y_{it}=\mu+ \theta_t + \sum_{j=1}^P\beta_{j}
	      X_{itj}+v_i(t)+\epsilon_{it}, i=1,...,n; t=1,...,T.}
	    The parameter 'mu' is set to zero if '-1' is used in \code{formula}.
	    \item \code{"twoways"}: This is the "twoways"-model ("within" &
	    "between"), which assumes that there are time-constant
	    individual effects, alpha_i, and a common time trend,
	    theta_t. The model can be written as
	    \eqn{$Y_{it}=\mu+ \alpha_{i} + \theta_t +\sum_{j=1}^P\beta_{j}
	      X_{itj}+\tau_i+v_i(t)+\epsilon_{it}\quad i=1,\dots,n; t=1,\dots,T.$}{
              Y_{it}=\mu+ \alpha_{i} + \theta_t +\sum_{j=1}^P\beta_{j}
	      X_{itj}+\tau_i+v_i(t)+\epsilon_{it}, i=1,...,n; t=1,...,T.}
	    The parameter 'mu' is set to zero if '-1' is used in \code{formula}.
	  }}}
	  Inferences about the slope parameters can be obtained by using the method \code{summary()}. The type of correlation and heteroskedasticity in the idiosyncratic errors can be specified by
choosing the corresponding number for the argument \code{error.type = c(1, 2, 3, 4, 5, 6, 7, 8)} in \code{summary()}, where 
	\itemize{
	\item \code{1}: indicates the presence of i.i.d. errors, 
	\item \code{2}: indicates the presence of cross-section heteroskedasticity with  \eqn{$n/T \to 0$}{n/T \to 0}, 
	\item \code{3}: indicates the presence of  cross-section correlation and heteroskedasticity with \eqn{$n/T \to 0$}{n/T \to 0},
	\item \code{4}: indicates the presence of  heteroskedasticity in the time dimension with \eqn{$T/n \to 0$}{T/n \to 0},
	\item \code{5}: indicates the presence of  correlation and heteroskedasticity in the time dimension with \eqn{$T/n \to 0$}{T/n \to 0},
	\item \code{6}: indicates the presence of  both time and cross-section dimensions with \eqn{$T/n^2 \to$}{T/n^2 \to} and \eqn{$n/T^2 \to 0$}{n/T^2 \to 0},
	\item \code{7}: indicates the presence of  both time and cross-section dimensions with \eqn{$n/T \to c > 0$}{n/T \to c > 0}, and
	\item \code{8}: indicates the presence of  correlation and heteroskedasticity in both time and cross-section dimensions with \eqn{$n/T \to c > 0$}{n/T \to c > 0}.
	}
	The default is \code{1}.  In presence of serial correlations
	(cases 5 and 8), the kernel weights required for estimating the
	long-run covariance  can be externally specified by given a
	vector of weights in the argument \code{kernel.weights}. By
	default, the function uses internally the linearly decreasing
	weights of Newey and West (1987) and a truncation at the lower
	integer part of \eqn{$min(\sqrt{n},\sqrt{T})$}{\min(\sqrt{n},\sqrt{T})}. If case 7 or 8 are chosen, the method \code{summary()} calculates the realization of the bias corrected estimators and gives appropriate inference. The bias corrected coefficients can be called by using the method \code{coef()} to the object produced by \code{summary()}. 
}
\value{
      'Eup' returns an object of 'class' '"Eup"' containing the
      following components:
      \itemize{
    \item{dat.matrix: }{Whole data set stored within a
	  (N*T)x(p+1)-Matrix, where P is the number of independent
	  variables without the intercept.}
    \item{formula: }{returns the used formula object.}
    \item{dat.dim: }{Vector of length 3: c(T,N,p)}
    \item{slope.para: }{Beta-parameters}
    \item{names: }{Names of the dependent and independent variables.}
    \item{is.intercept: }{logical.Used an intercept in the formula?: TRUE or FALSE}
    \item{additive.effects: }{Additive effect type. One of: "none","individual","time", "twoways".}
    \item{Intercept: }{Intercept-parameter. Tacks the value 0 if it is not specified in the model.}
    \item{Add.Ind.Eff: } {Estimated values of additive individual effects. If additive individual effects are not specified in the model, the function returns a vector of zeros.}
    \item{Add.Tim.Eff: } {Estimated values of additive time effects. If this effects are not specified in the model, the function returns a vector of zeros.}
    \item{unob.factors: } {Txd-matrix of estimated unobserved common
	  factors, where 'd' is the number of used factors.}
    \item{ind.loadings: } {Nxd-matrix of loadings parameters.}
    \item{unob.fact.stru: } {TxN-matrix of the estimated factor
	  structure. Each column represents an estimated individual unobserved time trend.}
    \item{used.dim: }{Used dimension 'd' to calculate the factor structure.}
    \item{proposed.dim: }{Indicates whether the user has specified the factor dimension or not.}
    \item{optimal.dim: }{The optimal dimension calculated internally.}
    \item{factor.dim: } {The user-specified factor dimension. Default is \code{NULL}}
    \item{d.max: } {The maximum number of factors used to estimate the optimal dimension.}
    \item{dim.criterion: } {The used dimensionality criterion.}
    \item{OvMeans: }{A vector that contains the overall means of the observed variables (Y and X).}
    \item{ColMean: } {A matrix that contains the column means of the observed variables (Y and X).}
    \item{RowMean: } {A matrix that contains the row means of the observed variables (Y and X).}
    \item{max.iteration: } {The maximum number of iterations. The default is '500'.}
    \item{convergence: }{The convergence condition. The default is '1e-6'.}
    \item{start.beta: A vector of user-specified starting values for the estimation of the beta-parameters. Default is \code{NULL}.} 
    \item{Nbr.iteration: } {Number of iterations required for the computation.} 		
    \item{fitted.values: } {Fitted values.}
    \item{orig.Y: } {Original values of the dependent variable.} 
    \item{residuals: } {Original values of the dependent variable.}   
    \item{sig2.hat.dim: }{user-specified variance estimator of the errors. Default is \code{NULL}.}
    \item{sig2.hat: } {Estimated variance of the error term.}     
    \item{degrees.of.freedom: }{Degrees of freedom of the residuals.}
    \item{call}
      }
}
\references{
  \itemize{
    \item Bai, J., 2009 \dQuote{Panel data models with interactive fixed
      effects}, \emph{Econometrica}
    \item Bada, O. and Kneip, A., 2014  \dQuote{Parameter Cascading for Panel Models with Unknown Number of Unobserved Factors: An Application to the Credit Spread Puzzle}, \emph{Computational Statistics \& Data Analysis (forthcoming)}
}}
\author{
  Oualid Bada
}
% \note{
% %%  ~~further notes~~
% }
% %% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
KSS, OptDim
}
\examples{
## See the example in 'help(Cigar)' in order to take a look at the
## data set 'Cigar'

##########
## DATA ##
##########

data(Cigar)
## Panel-Dimensions:
N <- 46
T <- 30
## Dependent variable:
  ## Cigarette-Sales per Capita
  d.l.Consumption    <- diff(log(matrix(Cigar$sales, T,N)))
## Independent variables:
  ## Consumer Price Index
  cpi        <- matrix(Cigar$cpi, T,N)
  ## Real Price per Pack of Cigarettes 
  d.l.Price  <- diff(log(matrix(Cigar$price, T,N)/cpi))
  ## Real Disposable Income per Capita  
  d.l.Income    <- diff(log(matrix(Cigar$ndi,   T,N)/cpi))
  
## Estimation:
Eup.fit      <- Eup(d.l.Consumption~d.l.Price+d.l.Income)
(Eup.fit.sum <- summary(Eup.fit))

## Plot the components of the estimated individual effects

plot(Eup.fit.sum)
}
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
