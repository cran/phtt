\name{OptDim}
\alias{OptDim}
\alias{OptDim.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimation of the Factor Dimension}
\description{Functions for the Estimation of the Factor Dimension}
\usage{
OptDim(Obj, 
       criteria = c("PC1", "PC2", "PC3", "BIC3",
                     "IC1", "IC2", "IC3",
                     "IPC1","IPC2", "IPC3",
                     "ABC.IC1", "ABC.IC2",
                     "KSS.C",
                     "ED",  "ER",  "GR"),
       standardize = FALSE,
       d.max,
       sig2.hat,
       spar,
       level = 0.01, 
       c.grid = seq(0, 5, length.out = 128),
       T.seq, n.seq)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Obj}{The function requires either a Txn matrix or an object with class "`Eup"' or "`KSS"'.}
  \item{criteria}{A character vector that contains the desired criteria to be used. If it is left unspecified, the function returns the result of all 16 criteria. }
  \item{standardize}{logical. If \code{TRUE} the input variable will be standardized. Default is \code{FALSE}.}
  \item{d.max}{Maximal dimension used in the dimensionality-criteria of Bai
    (2009). The default (\code{d.max=NULL}) yields to an internal selection of
    d.max.}
  \item{sig2.hat}{The squared standard deviation of the error-term required for the computation of some dimensionality criteria. The user can specify it in instead of \code{d.max}. The default (\code{sig2.hat=NULL}) yields to an internal estimation.}
  \item{spar}{Smoothing parameter used to calculate the criterion of Kneip, Sickles, and Song (2012). The default is \code{NULL}, which leads to internal computation.}
  \item{level}{The significance level used for the criterion of Kneip, Sickles, and Song (2012). The default is 0.01.}
  \item{c.grid}{Required only for computing \code{"ABC.IC1"} and \code{"ABC.IC2"}. It specifies the grid interval in which the scaling parameter of the penalty terms in \code{"ABC.IC1"} and \code{"ABC.IC2"} are calibrated. Default is \code{c.grid  =seq(0, 5, length.out = 128)}.}
  \item{T.seq}{Required only for computing \code{"ABC.IC1"} and \code{"ABC.IC2"}. It can be a vector containing different dimensions for T or an integer indicating the length of the sequence to be considered in calibrating  \code{"ABC.IC1"} and \code{"ABC.IC2"}. If it is left unspecified, the function determines internally a sequence of the form \code{seq((T-C), T)}, where \code{C} is the square root of \code{min{T,900}}. }
  \item{n.seq}{Required only for computing \code{"ABC.IC1"} and \code{"ABC.IC2"}. It can be a vector containing different dimensions for n or an integer indicating the length of the sequence to be considered in calibrating  \code{"ABC.IC1"} and \code{"ABC.IC2"}. If it is left unspecified, the function determines internally a sequence of the form \code{seq((n-D), n)}, where \code{D} is the square root of \code{min{n,900}}.}
}
\details{The function 'OptDim' allows for a comparison of the optimal factor dimensions obtained from different panel criteria (in total 13). This criteria are adjusted for panel data with diverging T and N.}
\value{
'OptDim' returns an object of 'class' '"OptDim"' containing a list with the
      following components:
\item{criteria:}{The name of the criteria specified by the user.}
\item{PC1:}{If specified in \code{criteria} a table is returned with the optimal dimension, the empirical standard deviation of the residuals, and some other informations required internally by the criterion, such as \code{d.max} and/ or \code{sig2.hat}.}
\item{PC2:}{If specified in \code{criteria} a table is returned with the optimal dimension, the empirical standard deviation of the residuals, and some other informations required internally by the criterion, such as \code{d.max} and/ or \code{sig2.hat}.}
\item{PC3:}{If specified in \code{criteria} a table is returned with the optimal dimension, the empirical standard deviation of the residuals, and some other informations required internally by the criterion, such as \code{d.max} and/ or \code{sig2.hat}.}
\item{IC1:}{If specified in \code{criteria} a table is returned with the optimal dimension, the empirical standard deviation of the residuals, and some other informations required internally by the criterion, such as \code{d.max}.}
\item{IC2:}{If specified in \code{criteria} a table is returned with the optimal dimension, the empirical standard deviation of the residuals, and some other informations required internally by the criterion, such as \code{d.max}.}
\item{IC3:}{If specified in \code{criteria} a table is returned with the optimal dimension, the empirical standard deviation of the residuals, and some other informations required internally by the criterion, such as \code{d.max}.}
\item{IPC1:}{If specified in \code{criteria} a table is returned with the optimal dimension, the empirical standard deviation of the residuals, and some other informations required internally by the criterion, such as \code{d.max} and/ or \code{sig2.hat}.}
\item{IPC2:}{If specified in \code{criteria} a table is returned with the optimal dimension, the empirical standard deviation of the residuals, and some other informations required internally by the criterion, such as \code{d.max} and/ or \code{sig2.hat}.}
\item{IPC3:}{If specified in \code{criteria} a table is returned with the optimal dimension, the empirical standard deviation of the residuals, and some other informations required internally by the criterion, such as \code{d.max} and/ or \code{sig2.hat}.}
\item{KSS.C:}{If specified in \code{criteria} a table is returned with the optimal dimension, the empirical standard deviation of the residuals, and some other informations required internally by the criterion, such as \code{d.max} and/ or \code{sig2.hat}.}
\item{ED:}{If specified in \code{criteria} a table is returned with the optimal dimension, the empirical standard deviation of the residuals, and some other informations required internally by the criterion, such as \code{d.max}.}
\item{ER:}{If specified in \code{criteria} a table is returned with the optimal dimension, the empirical standard deviation of the residuals, and some other informations required internally by the criterion, such as \code{d.max}.}
\item{GR:}{If specified in \code{criteria} a table is returned with the optimal dimension, the empirical standard deviation of the residuals, and some other informations required internally by the criterion, such as \code{d.max}.}
\item{summary:}{A table (in a matrix form) containing all the estimated dimensions obtained by the specified criteria.}
\item{BaiNgC:}{A logical vector required for further internal computations.}
\item{BaiC:}{A logical vector required for further internal computations.}
\item{KSSC:}{A logical vector required for further internal computations.}
\item{OnatC:}{A logical vector required for further internal computations.}
\item{RHC:}{A logical vector required for further internal computations.}
\item{obj:}{The argument '\code{Obj}' given to the function \code{OptDim()}.}
\item{cl:}{Object of mode "call".}
}
\references{
  \itemize{
      \item  Ahn, S. C., Horenstein, A. R. 2009 \dQuote{Eigenvalue ratio test for the number of factors}, \emph{Working Paper}
    \item Bai, J., 2009 \dQuote{Panel data models with interactive fixed
      effects}, \emph{Econometrica}
    \item Bai, J. 2004 \dQuote{Estimating cross-section common stochastic trends in nonstationary data},  \emph{Journal of Econometrics}
    \item Bai, J., Ng, S. 2009 \dQuote{Determining the number of factors in approximated factor models}, \emph{Econometrica}
    \item Kneip, A., Sickles, R. C., Song, W., 2012 \dQuote{A New Panel Data
      Treatment for Heterogeneity in Time Trends}, \emph{Econometric Theory}
    \item Onatski, A. 2010 \dQuote{Determining the number of factors from empirical distribution of eigenvalues}, \emph{The Review of Economics and Statistics}
}}
\author{
  Oualid Bada
}
% \note{
% %%  ~~further notes~~
% }

% %% ~Make other sections like Warning with \section{Warning }{....} ~

% \seealso{plot.OptDim, 
% }
\examples{
## See the example in 'help(Cigar)' in order to take a look at the
## data set 'Cigar'

##########
## DATA ##
##########

data(Cigar)
N <- 46
T <- 30


## Data: Cigarette-Sales per Capita
l.Consumption    <- log(matrix(Cigar$sales, T,N))

## Calculation is based on the covariance matrix of l.Consumption

OptDim(l.Consumption)

## Calculation is based on the correlation matrix of l.Consumption

OptDim(l.Consumption, standardize = TRUE)

## Display the magnitude of the eigenvalues in percentage of the total variance

plot(OptDim(l.Consumption))
}
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
