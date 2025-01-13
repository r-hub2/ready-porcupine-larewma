\name{drie}
\alias{drie}
\docType{data}
\title{Simulated DRIE profiles}
\description{
This data set contains five in-control profiles followed by nine out-control
profiles given by Zou \emph{et al.} (see references below).  
}
\usage{data(drie)}
\format{
  A \eqn{14 \times 11}{14x11}  matrix. Each row contains a profile. 
}
\references{
Zou, C., Tsung, F., Wang, Z. (2007), Monitoring general linear profiles
using multivariate Exponentially Weighted Moving Average schemes,
\emph{Technometris}, \bold{49}, 395-408.
}
\examples{
data(drie)
# Plot with different line types for the
# in-control and out-of-control profiles
x <- seq(-2.5,2.5,by=0.5)
plot(ts(t(drie),start=-2.5,deltat=0.5),plot.type="s",xlab="x",
ylab="drie profiles",
lty=c(rep("solid",5),rep("dashed",9))) 
legend(1,-0.1,c("in-control","out-control"),lty=c("solid","dashed"),bty="n")
}
\keyword{datasets}
