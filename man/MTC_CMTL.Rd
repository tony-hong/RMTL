\name{MTC_CMTL}
\alias{MTC_CMTL}

\title{
Multi-task classification with clustered structure
}

\description{
  The formulation combines the objective of logistic loss and convex form of k-means clustering, and therefore 
  is able to detect a clustered structure among tasks. The parameter (\eqn{\lambda_1}) is used to
  balance the clustering and data fitting effects, and can be learned via cross-validation.
}
\usage{
MTC_CMTL(X, Y, ...)
}

\arguments{
  \item{X}{
    a set of feature matrixes
}
  \item{Y}{
    a set of binary responses \eqn{\in \{-1,1\}}
}
  \item{k}{
    a positive number to modulate the stucture of clusters with the default of 2
}
  \item{Lam1}{
    a positive constant \eqn{\lambda_1} to control the clustering effect with default of 0.5  
}
  \item{lam2}{
    a positive constant \eqn{\lambda_2} to improve the
    generalization performance with default of 0.5
}

  \item{opts}{
    options of the solver. The default is
    \code{list(init = 0, tol = 10^-3, maxIter = 1000)}
}
}

\details{
  \deqn{\min\limits_{W,C}
  \sum_{i}^{t}{\frac{1}{n_i}sum(log(1+e^{<-Y_i^T, X_iW_i+C_i>}))} +
  \lambda_1\eta(1+\eta)tr(W(\eta I+M)^{-1}W^T)}
  \deqn{s.t. \quad tr(M)=k, M\preceq I, M\in S^t_+,
    \eta=\frac{\lambda_2}{\lambda_1}}

  \eqn{X} and \eqn{Y} are the sets of feature matrixes and binary responses respectively, 
  \eqn{W} is the coefficient matrix, and
  \eqn{t} is the number of tasks. Accordingly, \eqn{Y_i}, \eqn{X_i}, \eqn{W_i}
  and \eqn{n_i} refer to the  data, model parameter set and
  the number of subjects for task \eqn{i}. Note \eqn{W_i} is the \eqn{i}th
  column of \eqn{W}. \eqn{k} is the number of clusters. 

  Before training, \eqn{k}, \eqn{\lambda_1}, \eqn{\lambda_2} need to be specified by
  users in advance. 
}

\value{
  The function will return a trained MTC_CMTL model
  \item{W}{a matrix of features' coefficients}
  \item{C}{a constant vector(intercept) of all models}
  \item{Obj}{historical record of objective values}
  \item{fitted.values}{predictive scores(probability) of the training data.}
  \item{residuals}{the residuals of the training data. For each subject
  \eqn{i}, the residual is \eqn{y_i-\hat{y_i}}}
  \item{lam1}{\eqn{\lambda_1} value}
  \item{lam2}{\eqn{\lambda_2} value}
  \item{k}{the number of clusters}
  \item{opts}{options of the solver}
  \item{dim}{size of feature matrix of each task}
  \item{features}{feature names}
}

\references{
Zhou, J., Chen, J., & Ye, J. (2011). Clustered multi-task learning via
alternating structure optimization. In Advances in neural information
processing systems (pp. 702-710).

Jacob, L., Vert, J. P., & Bach, F. R. (2009). Clustered multi-task
learning: A convex formulation. In Advances in neural information
processing systems (pp. 745-752).
}

\author{han.cao@zi-mannheim.de}

\seealso{
 \code{\link{MTR_CMTL}}
 \code{\link{cv.MTC_CMTL}}
}

\examples{
#load the data
####
#X, Y: training data
#tX, tY: test data
#W: ground truth
#k: a possitive number to modulate the clustered structure
####
load('./RMTL/data/Simulated_Classification_CMTL.rda')

#cross-validation
cv <- cv.MTC_CMTL(X, Y, k=2)

#training
m <- MTC_CMTL(X, Y, lam1=cv$lam1.min, lam2=cv$lam2.min, k=2)

#predict on new dataset
predict(m, tX)

#compare the learnt model with the groud truth
library(fields)
par(mfrow=c(1,2))
image.plot(cor(W), xlab='tasks', ylab='features', main="ground truth")
image.plot(cor(m$W), xlab='tasks', ylab='features', main="CMTL")

#extract more information about the model
print(m)
plotObj(m)
m$fitted.values
m$residuals
m$opts
m$features
m$dim
}

\keyword{ clustering }
\keyword{ classification }
