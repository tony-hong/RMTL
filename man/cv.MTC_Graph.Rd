\name{cv.MTC_Graph}
\alias{cv.MTC_Graph}
\title{
Cross-validation of multi-task classification with network structure
}

\description{
  k-fold cross-validation for MTC_Graph generates a cv plot and estimates one 
  parameter. The default values of the arguments are indicated in the
  usage section. 
}

\usage{
cv.MTC_Graph(X, Y, G, opts = list(init = 0, tol = 10^-3, 
maxIter = 1000), stratify = FALSE, nfolds = 5, lam1 = 10^seq(1, -4, -1),
lam2=0)
}

\arguments{
  \item{X}{
    a set of feature matrixes
}
  \item{Y}{
    a set of binary responses \eqn{\in \{-1,1\}}
}
  \item{G}{
    a matrix to encode graph information. For more detail, please refer to \code{\link{MTC_Graph}}
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
\code{\link{MTC_Graph}}
\code{\link{MTR_Graph}}
}

\examples{
#load the data
load('./RMTL/data/Simulated_Classification_Graph.rda')

#specify the parameters
lam1=10^seq(2,-5, -0.05)

#cross-validation
cv <- cv.MTC_Graph(X, Y, G=G, lam1=lam1)

#best parameters
cv$lam1.min

#plot the cv error 
plot(cv)
}

\keyword{ graph }
\keyword{ classification }
