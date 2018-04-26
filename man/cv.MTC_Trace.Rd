\name{cv.MTC_Trace}
\alias{cv.MTC_Trace}
\title{
Cross-validation of multi-task classification with low-rank structure
}

\description{
  k-fold cross-validation for MTC_Trace generates a cv plot and estimates one 
  parameter. The default values of the arguments are indicated in the
  usage section. 
}

\usage{
cv.MTC_Trace(X, Y, lam2 = 0, opts = list(init = 0, tol = 10^-3, 
maxIter = 1000), stratify = FALSE, nfolds = 5, lam1 = 10^seq(2, -4, -1))
}

\arguments{
  \item{X}{
    a set of feature matrixes
}
  \item{Y}{
    a set of binary responses \eqn{\in \{-1,1\}}
}
  \item{lam2}{
    a positive constant \eqn{\lambda_2} to improve the generalization performance
}
  \item{stratify}{
   \code{stratify==TRUE} is used for stratified cross-validation
}
  \item{opts}{
    options of solver
}
  \item{nfolds}{
    number of folds
}
  \item{lam1}{
    a positive sequence of \eqn{\lambda_1} to control the low-rank structure
}
}

\details{
  The function first computes \code{nfolds} solution paths, then calculates
  the mean error for each possible choice of \eqn{\lambda_1}, and finally selects the solution
  with the minimum error. Please note, \eqn{\lambda_2} is pre-defined by
  users to avoid over-fitting, and not selected by cross-validation. 
}

\value{
  \item{lam1 }{\eqn{\lambda_1} sequence}
  \item{lam2 }{\eqn{\lambda_2} value}
  \item{lam1.min }{the selected \eqn{\lambda_1} with the lowest
    cross-validation error}  
  \item{cvm }{cross-validation error for each possible choice of \eqn{\lambda_1} }
}

\author{han.cao@zi-mannheim.de}



\seealso{
\code{\link{MTC_Trace}}
\code{\link{MTR_Trace}}
}

\examples{
#load the data
load('./RMTL/data/Simulated_Classification_Trace.rda')

#specify parameters
opts=list(init=0,  tol=10^-6, maxIter=10000)
lam1=10^seq(2,-3, -0.05)

#perform cv
cv <- cv.MTC_Trace(X, Y, lam1=lam1, opts=opts)

#best parameter
cv$lam1.min

#plot results
plot(cv)
}

\keyword{ low-rank }
\keyword{ classification }
