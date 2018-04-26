\name{MTR_Lasso}
\alias{MTR_Lasso}

\title{
Multi-task regression with low-rank structure
}

\description{
  This formulation extends the lasso to the multi-task scenario such that all
  models are penalized according to the same \eqn{L_1} strength, and all
  unimportant coefficients are shrunken to 0, to achieve the global optimun.
}

\usage{
MTR_Lasso(X, Y, ...)
}

\arguments{
  \item{X}{
    a set of feature matrixes
}
  \item{Y}{
    a set of continuous responses
}
  \item{Lam1}{
    a parameter \eqn{\lambda_1} to control the
    sparsity. The default is 1 
}
  \item{lam2}{
    a parameter \eqn{\lambda_2} to improve the
    generalization performance. The default is 0
}
  \item{opts}{
    options of the solver. The default is
    \code{list(init = 0, tol = 10^-3, maxIter = 1000)}
}
}

\details{
  \deqn{\min\limits_{W} \sum_{i}^{t}{\frac{1}{n_i}||Y_i- X_iW_i||^2 +
  \lambda_1||W||_1 + \lambda_2||W||_F^2 }}

  \eqn{X} and \eqn{Y} are the sets of feature matrixes and continuous responses respectively, 
  \eqn{W} is the coefficient matrix, and
  \eqn{t} is the number of tasks. Accordingly, \eqn{Y_i}, \eqn{X_i}, \eqn{W_i}
  and \eqn{n_i} refer to the  data, model parameter set and
  the number of subjects for task \eqn{i}. Note \eqn{W_i}
  is the \eqn{i}th column of \eqn{W}, and \eqn{C_i} is the \eqn{i}th
  term of \eqn{C}. \eqn{||\circ||_F} is the
  Frobenius norm. \eqn{||\circ||_1} is the L1 norm.
}

\value{
  The function will return a trained MTR_Lasso model
  \item{W}{a matrix of features' coefficients}
  \item{Obj}{historical record of objective values} 
  \item{fitted.values}{predictive value of the training data.}
  \item{residuals}{the residuals of the training data. For each subject
  \eqn{i}, the residual is \eqn{y_i-\hat{y_i}}}
  \item{lam1}{\eqn{\lambda_1} value}
  \item{lam2}{\eqn{\lambda_2} value}
  \item{opts}{options of the solver}
  \item{dim}{size of feature matrix of each task}
  \item{features}{feature names}
}

\references{
Beck, A., & Teboulle, M. (2009). A fast iterative shrinkage-thresholding algorithm for linear inverse problems. SIAM journal on imaging  sciences, 2(1), 183-202.

Tibshirani, R. (1996). Regression shrinkage and selection 
lasso. Journal of the Royal Statistical Society. Series B (Methodological), 267-288.
}




\seealso{
 \code{\link{MTC_Lasso}}
 \code{\link{cv.MTR_Lasso}}
}

\examples{
#load the data
####
#X, Y: training data
#tX, tY: test data
#W: ground truth
####
load('./RMTL/data/Simulated_Regression_Lasso.rda')

#speficy options
opts=list(init=0,  tol=10^-6, maxIter=10000)
lam1=10^seq(1,-6, -0.1)

#cross-validation
cv <- cv.MTR_Lasso(X, Y, lam1=lam1, opts=opts)

#training with the warm start
opt <- opts
for (i in 1: length(lam1)){
    r <- MTR_Lasso(X, Y, lam1=lam1[i], opts=opt) 
    opt$init=1;
    opt$W0=r$W;
    opt$C0=r$C;
    if (lam1[i]==cv$lam1.min) break    
}

#predict on new dataset
predict(r, tX)

#show results
par(mfrow=c(1,2))
library(fields)
image.plot(t(W!=0), xlab='tasks', ylab='features', main="ground truth")
image.plot(t(opt$W0!=0), xlab='tasks', ylab='features', main="Sparse")

#extract more information about the model
print(r)
plotObj(r)
r$fitted.values
r$residuals
r$opts
r$features
r$dim
}

\keyword{ sparse }
\keyword{ regression }
