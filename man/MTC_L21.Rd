\name{MTC_L21}
\alias{MTC_L21}

\title{
Multi-task classification with joint feature selection
}

\description{
  This formulation constraints all models to select or reject the same set of
  features simultaneously. Therefore, the solution only contains
  features which are consistently important to all tasks.
}

\usage{
MTC_L21(X, Y, ...)
}

\arguments{
  \item{X}{
    a set of feature matrixes
}
  \item{Y}{
    a set of binary responses \eqn{\in \{-1,1\}}
}
  \item{Lam1}{
    a positive constant \eqn{\lambda_1} to control the group
    sparsity. The default is 0.01  
}
  \item{lam2}{
    a positive constant \eqn{\lambda_2} to improve the
    generalization performance. The default is 0
}
  \item{opts}{
    options of the solver. The default is
    \code{list(init = 0, tol = 10^-3, maxIter = 1000)}
}

}

\details{
  \deqn{\min\limits_{W,C}
  \sum_{i}^{t}{\frac{1}{n_i}sum(log(1+e^{<-Y_i^T, X_iW_i+C_i>}))} +
  \lambda_1||W||_{2,1} + \lambda_2||W||_F^2 }

  \eqn{X} and \eqn{Y} are the sets of feature matrixes and binary responses respectively,  
  \eqn{W} is the coefficient matrix, and
  \eqn{t} is the number of tasks. Accordingly, \eqn{Y_i}, \eqn{X_i}, \eqn{W_i}
  and \eqn{n_i} refer to the  data, model parameter set and
  the number of subjects for task \eqn{i}. Note \eqn{W_i}
  is the \eqn{i}th column of \eqn{W}, and \eqn{C_i} is the \eqn{i}th
  term of \eqn{C}. \eqn{||\circ||_F} is the Frobenius norm.  

  \eqn{||W||_{2,1}=\sum{_i}{||W[i,]||_2}} aims to create the group
  sparse structure in the feature space. In the multi-task learning
  scenario, the same feature of all tasks form a group, such that features
  in the same group are equally penalized, while results in the group-wise sparsity.

}


\value{
  The function will return a trained MTC_L21 model
  \item{W}{a matrix of features' coefficients}
  \item{C}{a constant vector(intercept) of all models}
  \item{Obj}{historical record of objective values} 
  \item{fitted.values}{predictive scores(probability) of the training data.}
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

\author{han.cao@zi-mannheim.de}

\seealso{
 \code{\link{MTR_L21}}
 \code{\link{cv.MTR_L21}}
}

\examples{
#load the data
####
#X, Y: training data
#tX, tY: test data
#W: ground truth
####
load('./RMTL/data/Simulated_Classification_L21.rda')

#specify the parameters
opts=list(init=0,  tol=10^-6, maxIter=10000)
lam1=10^seq(1,-6, -0.1)

#cross-validation
cv <- cv.MTC_L21(X, Y, lam1=lam1, opts=opts)

#training with warm start
opt <- opts
for (i in 1: length(lam1)){
    r <- MTC_L21(X, Y, lam1=lam1[i], opts=opt) 
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
image.plot(t(opt$W0!=0), xlab='tasks', ylab='features', main="L21")

#extract more information about the model
print(r)
plotObj(r)
r$fitted.values
r$residuals
r$opts
r$features
r$dim
}

\keyword{ feature selection }
\keyword{ classification }
