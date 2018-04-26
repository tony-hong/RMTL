\name{cv.MTC_CMTL}
\alias{cv.MTC_CMTL}

\title{
Cross-validation of multi-task classification with clustered structure
}

\description{
  k-fold cross-validation for MTC_CMTL generates a cv plot and estimates two 
  parameters. The default values of the arguments are indicated in the
  usage section. 
}

\usage{
cv.MTC_CMTL(X, Y, k = 2, opts = list(init = 0, tol = 10^-3,
maxIter = 1000), stratify = FALSE, nfolds = 5, lam1 = 10^seq(1, -4, -1),
 lam2 = 10^seq(1, -4, -1))
}

\arguments{
  \item{X}{
    a set of feature matrixes
}
  \item{Y}{
    a set of binary responses \eqn{\in \{-1,1\}}
}
  \item{k}{
    a positive number to modulate the stucture of clusters 
}
  \item{opts}{
    options of solver
}
  \item{stratify}{
   \code{stratify==TRUE} is used for stratified cross-validation
}  \item{nfolds}{
    number of folds
}
  \item{lam1}{
    a positive sequence of \eqn{\lambda_1} to balance the clustering and data fitting effect
}
  \item{lam2}{
    a positive constant \eqn{\lambda_2} to improve the generalization performance
}
}

\details{
  The function first computes \code{nfolds} solution paths, then calculates
  the mean error for each pair of \eqn{\lambda_1} and
  \eqn{\lambda_2}, and finally selects the pair with the minimum error. 
}

\value{
  \item{lam1 }{\eqn{\lambda_1} sequence}
  \item{lam2 }{\eqn{\lambda_2} sequence}
  \item{lam1.min }{the selected \eqn{\lambda_1} with the lowest
    cross-validation error}
  \item{lam2.min }{the selected \eqn{\lambda_2} with the lowest
    cross-validation error}
  \item{cvm }{cross-validation error for each possible choices of \eqn{\lambda_1} and \eqn{\lambda_2} }
}

\author{han.cao@zi-mannheim.de}

\seealso{
\code{\link{MTC_CMTL}}
\code{\link{MTR_CMTL}}
}

\examples{
#load the data
load('./RMTL/data/Simulated_Classification_CMTL.rda')

#cross-validation
cv <- cv.MTC_CMTL(X, Y, k=2)

#selected parameters
cv$lam1.min
cv$lam2.min

#plot the cross-validation error 
plot(cv)
}

\keyword{ clustering }
\keyword{ classification }
