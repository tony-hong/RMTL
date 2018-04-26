\name{cv.MTR_Graph}
\alias{cv.MTR_Graph}
\title{
Cross-validation of multi-task regression with network structure
}

\description{
  k-fold cross-validation for MTR_Graph generates a cv plot and estimates one 
  parameter. The default values of the arguments are indicated in the
  usage section. 
}

\usage{
cv.MTR_Graph(X, Y, G, opts = list(init = 0, tol = 10^-3, 
maxIter = 1000), nfolds = 5, lam1 = 10^seq(3, -2, -1), lam2 = 0)
}

\arguments{ 
  \item{X}{
    a set of feature matrixes
}
  \item{Y}{
    a set of contiuous responses
}
  \item{G}{
    a matrix to encode graph information. For more detail, please refer to \code{\link{MTR_Graph}}
}
  \item{lam2}{ 
    a positive constant \eqn{\lambda_2} to improve the generalization performance
}
  \item{opts}{
    options of solver
}
  \item{nfolds}{
    number of folds
}
  \item{lam1}{
    a positive sequence of \eqn{\lambda_1} to control the network constraint
}
}


\details{
  The function first trains a model using each possible choice of \eqn{\lambda_1} in each fold, then calculates
  the mean error across folds, and finally selects the solution
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
\code{\link{MTR_Graph}}
\code{\link{MTC_Graph}}
}
\examples{
#load the data
load('./RMTL/data/Simulated_Regression_Graph.rda')

#speficy options
lam1=10^seq(2,-5, -0.05)

#cross-validation
cv <- cv.MTR_Graph(X, Y, G=G, lam1=lam1)

#selected parameter
cv$lam1.min

#plot the cross-validation error 
plot(cv)
}

\keyword{ graph }
\keyword{ regression }
