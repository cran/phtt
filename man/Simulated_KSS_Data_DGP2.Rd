\name{Simulated Data for the KSS-Model: DGP2}
\docType{data}
\alias{DGP2}
\title{Simulated Panel-Data Set with Polynomial Factor Structure and \emph{endogenous} regressors.}
\description{
  A Panel-Data Sets with:
  
  \emph{time-index} : t=1,\ldots,T=30

  \emph{individual-index} : i=1,\ldots,N=60

  This panel-data set has a polynomial factor structure (3 common
  factors) and \emph{endogenous} regressors.
}

\usage{
data(DGP2)
}

\format{A data frame containing :
  \describe{
    \item{Y}{dependent variable as N*T-vector}
    \item{X1}{first regressor as N*T-vector}
    \item{X2}{second regressor as N*T-vector}
    \item{CF.1}{first (unobserved) common factor: \eqn{$CF.1(t)=1$}{CF.1(t)=1}}
    \item{CF.2}{second (unobserved) common factor: \eqn{$CF.2(t)=\frac{t}{T}$}{CF.2(t)=(t/T)}}
    \item{CF.3}{thrid (unobserved) common factor: \eqn{$CF.3(t)=\left(\frac{t}{T}\right)^2$}{CF.3(t)=(t/T)^2}}    
    \item{}{Remark: The time-index t is running "faster" than the
    individual-index i such that e.g. Y_{it} is ordered as:
    \eqn{$Y_{11},Y_{12},\ldots,Y_{1T},Y_{21},Y_{22},\ldots$}{Y_{11},Y_{12},\ldots,Y_{1T},Y_{21},Y_{22},\ldots}}
  }
}

\details{
  The panel-data set DPG2 is simulated according to the simulation-study
  in Kneip, Sickles & Song (2012):
  \eqn{$Y_{it}=\beta_{1}X_{it1}+\beta_{2}X_{it2}+v_i(t)+\epsilon_{it}\quad
    i=1,\dots,n;\quad
    t=1,\dots,T$}{Y_{it}=\beta_{1}X_{it1}+\beta_{2}X_{it2}+v_i(t)+\epsilon_{it},
    i=1,\dots,n, t=1,\dots,T}
  -Slope parameters:
  \eqn{$\beta_{1}=\beta_{2}=0.5$}{beta_{1}=beta_{2}=0.5}

  -Time varying individual effects being second order polynomials:
  \eqn{$v_i(t)=\theta_{i0}+\theta_{i1}\frac{t}{T}+\theta_{i2}\left(\frac{t}{T}\right)^2$}{v_i(t)=theta_{i0}+theta_{i1}*frac{t}{T}+theta_{i2}*(frac{t}{T})^2}
  Where theta_{i1}, theta_{i1}, and theta_{i1} are iid as N(0,4)

  The Regressors X_{it}=(X_{it1},X_{it2})' are simulated from a bivariate VAR model:
  \eqn{$X_{it}=R X_{i,t-1}+\eta_{it}\quad\textrm{with}\quad
  R=\left(\begin{array}{cc}0.4&0.05\\0.05&0.4\end{array}\right)\quad\textrm{and}\quad
  \eta_{it}\sim N(0,I_2)$}{X_{it}=R X_{i,t-1}+eta_{it} with
  R=matrix(c(0.4,0.05,0.05,0.4),2,2) and eta_{it}~N(0,I_2)}

  After this simulation, the N regressor-series
  \eqn{$(X_{1i1},X_{2i1})',\dots,(X_{1iT},X_{2iT})'$}{(X_{1i1},X_{2i1})',\dots,(X_{1iT},X_{2iT})'}
  are additionally shifted such that there are three different
  mean-value-clusters. Such that every third of the N regressor-series
  fluctuates around on of the following mean-values
  \eqn{$\mu_1=(5,5)',\;\mu_2=(7.5,7.5)',\textrm{ and }\;\mu_3=(10,10)'$}{mu_1=(5,5)', mu_2=(7.5,7.5)', and mu_3=(10,10)'}
  
  In this Panel-Data Set the regressor X_{it2} is made endogenous by the
  re-definition:
  \eqn{$X_{it2}:=X_{it2}+0.5 v_i(t)$}{X_{it2}:=X_{it2}+0.5*v_i(t)}

   See Kneip, Sickles & Song (2012) for more details.
}

\references{
  \itemize{
    \item Kneip, A., Sickles, R. C., Song, W., 2012 \dQuote{A New Panel Data
      Treatment for Heterogneity in Time Trends}, \emph{Econometric Theory}
  }
}
\author{
  Dominik Liebl
}
\keyword{datasets}
\examples{
data(DGP2)

## Dimensions
N    <- 60
T    <- 30

## Observed Variables
Y    <- matrix(DGP2$Y,  nrow=T,ncol=N)
X1   <- matrix(DGP2$X1, nrow=T,ncol=N)
X2   <- matrix(DGP2$X2, nrow=T,ncol=N)

## Unobserved common factors
CF.1 <- DGP2$CF.1[1:T]
CF.2 <- DGP2$CF.2[1:T]
CF.3 <- DGP2$CF.3[1:T]

## Take a look at the simulated data set DGP2:
par(mfrow=c(2,2))
matplot(Y,  type="l", xlab="Time", ylab="", main="Depend Variable")
matplot(X1, type="l", xlab="Time", ylab="", main="First Regressor")
matplot(X2, type="l", xlab="Time", ylab="", main="Second Regressor")
## Usually unobserved common factors:
matplot(matrix(c(CF.1,
                 CF.2,
                 CF.3), nrow=T,ncol=3),
        type="l", xlab="Time", ylab="", main="Unobserved Common Factors")
par(mfrow=c(1,1))

## Esimation
KSS.fit      <- KSS(Y~-1+X1+X2)
(KSS.fit.sum <- summary(KSS.fit))

plot(KSS.fit.sum)
}