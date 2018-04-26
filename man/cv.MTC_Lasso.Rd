\name{cv.MTC_Lasso}
\alias{cv.MTC_Lasso}

\title{
Cross-validation of multi-task classification with L1 regularizer
}

\description{
  k-fold cross-validation for MTC_Lasso generates a cv plot and estimates one 
  parameter. The default values of the arguments are indicated in the
  usage section. 
}

\usage{
cv.MTC_Lasso(X, Y, lam2 = 0, opts = list(init = 0, tol = 10^-3, 
maxIter = 1000), stratify = FALSE, nfolds = 5, lam1 = 10^seq(1, -5, -1))
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
  \item{opts}{
    options of solver
}
  \item{stratify}{
   \code{stratify==TRUE} is used for stratified cross-validation
    
}
  \item{nfolds}{
    number of folds
}
  \item{lam1}{
    a positive sequence of \eqn{\lambda_1} to control the sparse constraint
}
}

\details{
  The function will compute \code{nfolds} solution paths, and calculate
  the mean error for each candidate \eqn{\lambda_1}, then selects the solution
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
\code{\link{MTC_Lasso}}
\code{\link{MTR_Lasso}}
}

\examples{
#load the data
load('./RMTL/data/Simulated_Classification_Lasso.rda')

#specify the parameters
opts=list(init=0,  tol=10^-6, maxIter=10000)
lam1=10^seq(1,-6, -0.1)

#perform cv
cv <- cv.MTC_Lasso(X, Y, lam1=lam1, opts=opts)

#best parameter
cv$lam1.min

#plot cv errors
plot(cv)
}

\keyword{ sparse}
\keyword{ classification }
