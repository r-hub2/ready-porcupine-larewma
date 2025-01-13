\name{larewma.simrl}
\alias{larewma.simrl}
\title{Run-length simulation}
\description{This function simulates \code{nsim} run-lengths of a LAR-EWMA control chart}
\usage{larewma.simrl(chart, nsim, tau = 0, delta = 0, omega = chart$S)}
\arguments{
  \item{chart}{representation of a LAR-EWMA control chart as returned by  \code{larewma}.}
  \item{nsim}{number of run lengths to simulate.}
  \item{tau}{time of the shift in process parameters.}
  \item{delta}{vector specifying the shifts in the mean vector.
    The \code{delta} argument can also be a scalar;
    in this case, a mean shift of \code{(delta,...,delta)}
    is assumed.}
  \item{omega}{out-of-control variance and covariance matrix.}
}
\value{A integer vector of \code{nsim} run lengths.}
\references{
  Capizzi, Giovanna and Masarotto, Guido (2010),
  Least Angle Regression Control Chart for Multidimensional Data,
  \emph{Technometrics}, forthcoming.
}

\author{Giovanna Capizzi and Guido Masarotto}
\seealso{\code{\link{larewma}}}
\examples{
# DRIE example
# Simulation of the run-length distribution 
# for a shift occurring at tau=51
x <- seq(-2.5,2.5,by=0.5)
F <- cbind(1,x,x^2-2.5)
S <- diag(rep(0.4^2,length(x)))
lambda <- 0.2
h <- 4.995
chart <- larewma(F,S,lambda,h=h)
nsim <- 2000
tau <- 51

# Shift in the intercept
rl <- larewma.simrl(chart,nsim,tau=tau,delta=0.05)
cat("\nPercent of false alarms before tau =",100*mean(rl<tau),
  "\nE(RL|RL>=tau) =",mean((rl-tau+1)[rl>=tau]),
  "\nsd(RL|RL>=tau) =",sd((rl-tau+1)[rl>=tau]),"\n")

# Shift in beta3
rl <- larewma.simrl(chart,nsim,tau=tau,delta=0.05*(x^2-2.5))
cat("\nPercent of false alarms before tau =",100*mean(rl<tau),
  "\nE(RL|RL>=tau) =",mean((rl-tau+1)[rl>=tau]),
  "\nsd(RL|RL>=tau) =",sd((rl-tau+1)[rl>=tau]),"\n")

# Shift in beta3 and omega
rl <- larewma.simrl(chart,nsim,tau=tau,delta=0.05*(x^2-2.5),omega=2*S)
cat("\nPercent of false alarms before tau =",100*mean(rl<tau),
  "\nE(RL|RL>=tau) =",mean((rl-tau+1)[rl>=tau]),
  "\nsd(RL|RL>=tau) =",sd((rl-tau+1)[rl>=tau]),"\n")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% __ONLY ONE__ keyword per line
\keyword{htest}
\keyword{multivariate}
\keyword{statistical process control}