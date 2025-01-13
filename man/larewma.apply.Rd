\name{larewma.apply}
\alias{larewma.apply}
\title{Apply a LAR-EWMA control chart}
\description{The function can be used to compute and plot the LAR-EWMA
  control statistic.
It also provides the diagnostic information described in the reference below.}
\usage{
larewma.apply(chart, y, mu, plot = TRUE )
}
\arguments{
  \item{chart}{representation of a LAR-EWMA control chart as returned by  \code{larewma}.}
  \item{y}{vector or matrix of data to be charted}
  \item{mu}{in-control mean vector. If missing, \code{mu=0} is assumed.}
  \item{plot}{logical. If  \code{TRUE} (default), the LAR-EWMA control
    statistic and its decision interval are plotted.}
}
\value{A list showing for each instant of time for which an alarm is signalled:
\item{t}{time of the alarm}
\item{W}{LAR-EWMA control statistic computed at time \code{t}.}
\item{stat}{diagnostic table containing the standardized control
  statistics \eqn{S_{t,k}} and the corresponding estimates of
  beta coefficients.}
}
\references{
  Capizzi, Giovanna and Masarotto, Guido (2010),
  Least Angle Regression Control Chart for Multidimensional Data,
  \emph{Technometrics}, forthcoming.
}

\author{Giovanna Capizzi and Guido Masarotto}
\seealso{\code{\link{larewma}}}
\examples{
data(drie)
x <- seq(-2.5,2.5,by=0.5)
mu <- 1.55+0.62*(x^2-2.5)
Sigma <- diag(rep(0.4^2,11))
F <- cbind(1,x,x^2-2.5)
h <- 4.995
lambda <- 0.2
chart <- larewma(F,Sigma,lambda,h=h)
larewma.apply(chart,drie,mu)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% __ONLY ONE__ keyword per line
\keyword{htest}
\keyword{multivariate}
\keyword{statistical process control}

