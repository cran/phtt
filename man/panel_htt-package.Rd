\name{phtt-package}
\Rdversion{1.1}
\alias{phtt-package}
\alias{phtt}
\docType{package}
\title{Panel Data Analysis with Heterogeneous Time Trends}
\description{
  The R-package \code{phtt} provides estimation procedures for panel data
  with large dimensions n, T, and general forms of unobservable
  heterogeneous effects. Particularly, the estimation procedures are
  those of Bai (2009) and Kneip, Sickles, and Song (2012), which
  complement one another very well: both models assume the unobservable
  heterogeneous effects to have a factor structure. The method of Bai
  (2009) assumes that the factors are stationary, whereas the method of
  Kneip et al. (2012) allows the factors to be
  non-stationary. Additionally, the \code{phtt} package provides a wide range
  of dimensionality criteria in order to estimate the number of the
  unobserved factors simultaneously with the remaining model
  parameters.}
\details{
  \tabular{ll}{
    Package: \tab phtt\cr
    Type: \tab Package\cr
    Version: \tab 3.1\cr
    Date: \tab 2013-09-26\cr
    License: \tab GPL-2\cr
    LazyLoad: \tab yes \cr
  }
}
\author{
  Oualid Bada, 
  Dominik Liebl
}
\references{
  \itemize{
    \item Kneip, A., Sickles, R. C., Song, W., 2012 \dQuote{A new panel data
      treatment for heterogeneity in time trends}, \emph{Econometric Theory}
    \item Bai, J., 2009 \dQuote{Panel data models with interactive fixed
      effects}, \emph{Econometrica}
}}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
%~~ the R documentation directory ~~
% \keyword{package}
%, panel data, factor models}
%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
% \examples{

% }
