\name{checkSpecif}
\alias{checkSpecif}
\alias{checkSpecification}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Tests to check the model specifications
}
\description{
Tests to check the model specifications
}
\usage{
checkSpecif(obj1, obj2, level = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj1}{
    \itemize{
      \item{If \code{obj2} is left unspecified the specification test
      proposed by Kneip, Sickles, and Song (2012) is computed. In this
      case \code{obj1} can be an object of class '\code{KSS}' or
      '\code{Eup}'. The given \code{KSS}- or \code{Eup}-model needs unspecified
      factor dimensions (\code{factor.dim=NULL}). See also the Details.}
      \item{If \code{obj2} is specified by an object of class
      '\code{Eup}' the Hausman-type test proposed by Bai (2009) is
      computed, which compares the model in \code{obj1} with the model in
      \code{obj2}. The Hausman-type test of Bai applies to
      '\code{Eup}'-objects only. See also the Details.}
    }}
  \item{obj2}{An object of class '\code{Eup}'. If \code{obj2} is left
    unspecified the test proposed by Kneip, Sickles, and Song (2012) is
    computed.
  }
  \item{level}{The significance level.}
}
\details{This function is equipped with two types of specification-tests:
  \itemize{
    \item{The first specification-tests is the Hausman-type test
      proposed by Bai (2009), which is computed if \code{obj1} as well
      as \code{obj2} are specified by '\code{Eup}'-objects. In this case
      the model given to the arguments \code{obj1} and \code{obj2} are
      compared with each other. Note that this test assumes that the
      (unobserved) true factor dimension is lower or equal to 2; see
      Section 9 in Bai (2009) for more details.
      
      Given the assumption that there is only one unobserved common factor:
      \itemize{
	\item{Null-Hypothesis: 'The unobserved common factor is a
	  classical \code{individual} or \code{time} effect'.}
	\item{Alternative-Hypothesis: 'The unobserved common factor is
	  an arbitrary process'.}
      }
      Given the assumption that there are two unobserved common factors:
      \itemize{
	\item{Null-Hypothesis: 'The two unobserved common factors are classical
	  \code{twoways} effects'.}
	\item{Alternative-Hypothesis: 'The two unobserved common factors are
	  arbitrary processes'.}
    }} 
    \item{The second specification-test tests the existence of an additional
      factor structure beyond a classical additive effects model; as suggested
      in Kneip, Sickles, and Song (2012), which is applied if only
      \code{obj1} is specified and \code{obj2} is left unspecified. This
      test can be used for '\code{Eup}'-objects as well as for
      '\code{KSS}'-objects.
      \itemize{
	\item{Null-Hypothesis: 'There are no unobserved common factors
      beyond the classical \code{individual}, \code{time}, or
      \code{twoways} effects'.}
	\item{Alternative-Hypothesis: 'There are additional unobserved
	  common factors'.}
      }
}}}
%\value{}
\references{
  \itemize{
    \item Bai, J., 2009 \dQuote{Panel data models with interactive fixed
      effects}, \emph{Econometrica}
    \item Kneip, A., Sickles, R. C., Song, W., 2012 \dQuote{A New Panel Data
      Treatment for Heterogeneity in Time Trends}, \emph{Econometric Theory}
}}
\author{
  Oualid Bada, Dominik Liebl
}
% \note{
% %%  ~~further notes~~
% }
% %% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
KSS, Eup, OptDim
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
  l.Consumption      <- log(matrix(Cigar$sales, T,N))
  d.l.Consumption    <- diff(l.Consumption)
## Independent variables:
  ## Consumer Price Index
  cpi           <- matrix(Cigar$cpi, T,N)
  ## Real Price per Pack of Cigarettes 
  l.Price       <- log(matrix(Cigar$price, T,N)/cpi)
  d.l.Price     <- diff(l.Price)
  ## Real Disposable Income per Capita  
  l.Income      <- log(matrix(Cigar$ndi,   T,N)/cpi)
  d.l.Income    <- diff(l.Income)

#####################################################################
## Testing the Sufficiency of a classical 'twoways' effects model: ##
## Hausman-type Test of Bai (2009)                                 ##
#####################################################################

## Model under the null Hypothesis:
twoways.obj <- Eup(d.l.Consumption ~ -1 + d.l.Price + d.l.Income,
	factor.dim = 0, additive.effects = "twoways")

## Model under the alternative Hypothesis:
not.twoways.obj <- Eup(d.l.Consumption ~ -1 + d.l.Price + d.l.Income,
	factor.dim = 2, additive.effects = "none")

###########
## Test: ##
###########

## (This test returns an error message, since the (unobserved) true
## factor dimension is probably greater than 2.)
\dontrun{
checkSpecif(obj1 = twoways.obj, obj2 = not.twoways.obj, level = 0.01)
}

#####################################################################
## Testing the Existence of additional (unobserved) common Factors ##
## Specification Test of Kneip, Sickles, and Song (2012)           ##
#####################################################################

## For the model of Bai (2009):
Eup.obj <- Eup(d.l.Consumption ~ -1 + d.l.Price + d.l.Income,
	   additive.effects = "twoways")

## Test:
checkSpecif(Eup.obj, level = 0.01)

## For the model of Kneip, Sickles, and Song (2012):
KSS.obj <- KSS(l.Consumption ~ -1 + l.Price + l.Income,
	   additive.effects = "twoways")

## Test:
checkSpecif(KSS.obj, level = 0.01)}
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line