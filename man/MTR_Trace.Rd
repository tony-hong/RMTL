\name{MTR_Trace}
\alias{MTR_Trace}

\title{
Multi-task regression with trace norm regularizer
}

\description{
  This formulation constraints all models to a low-rank subspace. With
  increasing  penalty(\eqn{\lambda_1}), the correlation between models increases.
}

\usage{
MTR_Trace(X, Y, ...)
}

\arguments{
  \item{X}{
    a set of feature matrixes
}
  \item{Y}{
    a set of continuous responses
}
  \item{Lam1}{
    a positive number (\eqn{\lambda_1}) to constraint the rank
    of \eqn{W}. The default is 1 
}
  \item{lam2}{
    a positive number (\eqn{\lambda_2}) to improve the
    generalization performance. The default is 0
}
  \item{opts}{
    options of the solver. The default is
    \code{list(init = 0, tol = 10^-3, maxIter = 1000)}
}
}
 
\details{
  \deqn{\min\limits_{W} \sum_{i}^{t}{\frac{1}{n_i}||Y_i- X_iW_i||^2 +
  \lambda_1||W||_* + \lambda_2||W||_F^2 }}

  \eqn{X} and \eqn{Y} are the sets of feature matrixes and continuous responses respectively, 
  \eqn{W} is the coefficient matrix, and
  \eqn{t} is the number of tasks. Accordingly, \eqn{Y_i}, \eqn{X_i}, \eqn{W_i}
  and \eqn{n_i} refer to the  data, model parameter set and
  the number of subjects for task \eqn{i}. Note \eqn{W_i}
  is the \eqn{i}th column of \eqn{W}, and \eqn{C_i} is the \eqn{i}th
  term of \eqn{C}. \eqn{||\circ||_F} is the Frobenius norm. \eqn{||W||_*}
  is the trace norm of \eqn{W}.
}

\value{
  The function will return a trained MTR_Trace model
  \item{W}{a matrix of features' coefficients}
  \item{Obj}{historical record of objective values} 
  \item{fitted.values}{predictive values of the training data.}
  \item{residuals}{the residuals of the training data. For each subject
  \eqn{i}, the residual is \eqn{y_i-\hat{y_i}}}
  \item{lam1}{\eqn{\lambda_1} value}
  \item{lam2}{\eqn{\lambda_2} value}
  \item{opts}{options of the solver}
  \item{dim}{size of feature matrix of each task}
  \item{features}{feature names}
}

\references{
Ji, S., & Ye, J. (2009, June). An accelerated gradient method for trace norm minimization. In Proceedings of the 26th annual international conference on machine learning (pp. 457-464). ACM.

Pong, T. K., Tseng, P., Ji, S., & Ye, J. (2010). Trace norm
regularization: Reformulations, algorithms, and multi-task learning. SIAM Journal on Optimization, 20(6), 3465-3489.
}

\author{han.cao@zi-mannheim.de}

\seealso{
 \code{\link{MTC_Trace}}
 \code{\link{cv.MTR_Trace}}
}

\examples{
#load the data
####
#X, Y: training data
#tX, tY: test data
#W: ground truth
####
load('./RMTL/data/Simulated_Regression_Trace.rda')

#speficy options
opts=list(init=0,  tol=10^-6, maxIter=10000)
lam1=10^seq(2,-3, -0.05)


#cross-validation
cv <- cv.MTR_Trace(X, Y, lam1=lam1, opts=opts)

#training
opt <- opts
for (i in 1: length(lam1)){
    r <- MTR_Trace(X, Y, lam1=lam1[i], opts=opt) 
    opt$init=1;
    opt$W0=r$W;
    opt$C0=r$C;
    if (lam1[i]==cv$lam1.min) break    
}

#predict on new dataset
predict(r, tX)

#compare the learnt model with the groud truth
par(mfrow=c(1,2))
library(fields)
image.plot(cor(W), xlab='tasks', ylab='tasks', main="ground truth")
image.plot(cor(opt$W0), xlab='tasks', ylab='tasks', main="Trace")

#extract more information about the model
print(r)
plotObj(r)
r$fitted.values
r$residuals
r$opts
r$features
r$dim
}

\keyword{ regression }
\keyword{ low-rank }
