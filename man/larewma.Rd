\name{larewma}
\alias{larewma}

\title{LAR-EWMA design}
\description{The function is used to create an object representing a LAR-EWMA control chart.}
\usage{
larewma(F, S, lambda, arl0, h, only.mean = FALSE, horiz = 20000,
A = 3, B = 0.6, csi = 5, burnin = 1000, iter = 50000, seed = 11226457)
}
\arguments{
  \item{F}{matrix of known constants defining an approximate relation
    between mean shifts and regression coefficients.}
  \item{S}{in-control variance and covariance matrix.}
  \item{lambda}{smoothing parameter \code{(0 < lambda <= 1)}.}
  \item{arl0}{in-control ARL; if not missing, control limit \code{h} leading to the specified
    \code{arl0} is computed via stochastic approximation (see details).}
  \item{h}{control limit. It can be missing if \code{arl0} is given.}
  \item{only.mean}{logical. If \code{TRUE}, only the process mean is
    tested for a change.}
  \item{horiz}{length of the in-control realization used to estimate
    means and standard deviations, say \eqn{a_k} and \eqn{b_k}, of the \eqn{S_{t,k}}
    control statistics. See **** and **** (2010).}
  \item{A, B, csi, burnin, iter}{suitable tuning constants of stochastic
    approximation. See details.}
  \item{seed}{the seed to be used for random sampling: see \code{\link{.Random.seed}}.
The current value of .Random.seed will be preserved.}
}
\details{When argument \code{arl0} is given,  the control limit \code{h} is computed
using the Polyak-Ruppert stochastic approximation algorithm and following
the implementation details given in Capizzi and Masarotto (2008).
In particular, \code{burnin+iter} values of \code{h} are generated
using the following recursion
  \deqn{h_{i+1}=\max\{0,h_i+g_i \min[(RL_i-\code{arl0})/\code{arl0},\code{csi}/g_i)\}}{h_{i+1}=max(0,h_i+g_i min((RL_i-arl0)/arl0,csi/g_i))}
where \eqn{RL_i} denotes a simulated in-control run length and \eqn{g_i=A/i^B}.
The control limit estimate is equal to the mean of the last
\code{iter} \eqn{h_i} values.
}
\value{A list of all the function arguments, with \code{h} eventually computed
as described before, plus
\item{ab}{matrix containing the values of \eqn{a_k} and \eqn{b_k}.}
\item{ipar,par,stat0}{internal representation of the control chart used to communicate
with the C functions.}
}
\references{
  Capizzi, Giovanna and Masarotto, Guido (2010),
  Least Angle Regression Control Chart for Multidimensional Data,
  \emph{Technometrics}, forthcoming.

  Capizzi, G. and Masarotto, G. (2008), Practical Design of Generalized Likelihood
  Ratio Control Charts for Autocorrelated Data, \emph{Technometrics}, \bold{50}, pp.357-370.

  Polyak, B. T. (1990), New Method of Stochastic Approximation Type, \emph{Automatic Remote Control},
  \bold{51}, pp. 937--946.
}



\author{Giovanna Capizzi and Guido Masarotto}

\seealso{\code{\link{larewma.simrl}}, \code{\link{larewma.apply}}}
\examples{
# Linear profile
# l = 0.2
# Desired IC ARL= 200
x <- seq(-1,1,by=0.5)
F <- cbind(1,x)
S <- 3*diag(rep(1,length(x)))
# Only a small number of stochastic approximation iteration
chart <- larewma(F,S,0.2,200,iter=5000)
chart$h
chart$ab
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% __ONLY ONE__ keyword per line
\keyword{htest}
\keyword{multivariate}
\keyword{statistical process control}
