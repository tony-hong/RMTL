\name{RMTL-package}
\alias{RMTL-package}
\alias{RMTL}
\docType{package}

\title{
\packageTitle{RMTL}
}

\description{
\packageDescription{RMTL}
}

\details{
  This package provides 10 multi-task learning 
  algorithms (5 classification and 5 regression), which incorporate five
  regularization strategies for knowledge transfering among tasks. All
  algorithms share the same framework:
  
  \deqn{\min\limits_{W}
    \sum_{i}^{t}{\frac{1}{n_i}L(W_i|X_i, Y_i)} + \Omega(W)}

  where \eqn{L(\circ)} is the loss function (logistic loss
  for classification or least square loss for regression). \eqn{X} and \eqn{Y} are sets of 
  feature matrixes and responses respectively, \eqn{W} is the coefficient matrix, and \eqn{t} is the
  number of tasks. Therefore \eqn{Y_i}, \eqn{X_i}, \eqn{W_i} and \eqn{n_i}
  refer to the data, model parameter vector and the number of subject for
  task \eqn{i}. Note \eqn{W_i} is the \eqn{i}th column of \eqn{W}.

  Knowledge exchange is achieved via the convex term \eqn{\Omega(W)}, which 
  jointly modulates models according their specific functionalities. In this package, 
  5 common regularization methods are implemented to incorperate different priors, i.e. 
  sparse structure (\eqn{\Omega(W)=||W||_1}), joint feature selection
  (\eqn{\Omega(W)=||W||_{2,1}}), low-rank structure
  (\eqn{\Omega(W)=||W||_*}), network-based relatedness across tasks
  (\eqn{\Omega(W)=||WG||_F^2}) and task clustering
  (\eqn{\Omega(W)=\lambda_1\eta(1+\eta)tr(W(\eta I+M)^{-1}W^T)}).
  
  For all algorithms, we implemented an solver based on the accelerated
  gradient descent method, which takes advantage of information from the
  previous two iterations to calculate the current gradient and then
  achieves an improved convergent rate. To solve the non-smooth and convex
  regularizer, the proximal operator is applied. Moreover, backward
  line search is used to determine the appropriate step-size in each
  iteration. Overall, the solver achieves a complexity of \eqn{O(\frac{1}{k^2})} 
   and is optimal among first-order gradient descent methods.
  
  \code{opts = list(init = 0, tol = 10^-3, maxIter = 1000)}

  These options are used to control the optimization procedure and can
  be customized by users. \code{opts$init} specifies the starting point
  of the gradient descent algorithm, \code{opts$tol} controls tolerance of
  the acceptable precision of solution to terminate the algorithm, and
  \code{opts$maxIter} is the maximum number of iterations.
  
  For \code{opts$init}, two options are provided. \code{opts$init==0}
  refers to 0 matrix. And \code{opts$init==1} refers to the
  user-specific starting point. If specified, the algorithm will attempt
  to access \code{opts$W0} and \code{opts$C0} as the starting point, which have to
  be given by users in advance. Otherwise, errors are reported. Particularly, the setting
  \code{opts$init==1} is key to warm-start technique for sparse model training.
  

  Only two core functions are necessary for application of each algorithm. 
  The naming of functions follows the structure:
  
  1) Training procedure: MTR_<name_of_prior> or MTC_<name_of_prior>

  i.e. MTR_L21 or MTC_L21
  
  2) cross-validation procedure: cv.MTR_<name_of_prior> or
  cv.MTC_<name_of_prior>
  
  i.e. cv.MTR_L21 or cv.MTC_L21

  Here, "MTR" is short for multi-task regression, and "MTC" is short for
  multi-task classification. Available <name_of_prior> include L21, Lasso,
  Trace, Graph and CMTL.
}

\author{
\packageAuthor{RMTL}
Maintainer: \packageMaintainer{RMTL}
}

\references{
Beck, A., & Teboulle, M. (2009). A fast iterative shrinkage-thresholding
algorithm for linear inverse problems. SIAM journal on imaging
sciences, 2(1), 183-202.

Parikh, N., & Boyd, S. (2014). Proximal algorithms. Foundations and
Trends® in Optimization, 1(3), 127-239.  

Nesterov, Y. (2007). Gradient methods for minimizing composite objective
function.
}

\keyword{ package }

\examples{

#classification example
#load the data
####
#X, Y: training data
#tX, tY: test data
#W: ground truth
####
load('./RMTL/data/Simulated_Classification_L21.rda')

#specify the options
myopt<-list(init = 0, tol = 10^-3, maxIter = 1000)
lam1=10^seq(1,-6, -0.1)

#perform cross-validation and training using the user-specific options
cv <- cv.MTC_L21(X, Y, opts=myopt, lam1=lam1)
m <- MTC_L21(X, Y, opts=myopt, lam1=cv$lam1.min)

#predict on test data
predict(m, tX)

#--------------------------------------------------------------
# regression example
#load the data
####
#X, Y: training data
#tX, tY: test data
#W: ground truth
####
load('./RMTL/data/Simulated_Regression_L21.rda')

#perform cv and train using the default options
cv <- cv.MTC_L21(X, Y)
m <- MTC_L21(X, Y, lam1=cv$lam1.min)

#predict on test data
predict(m, tX)
}
