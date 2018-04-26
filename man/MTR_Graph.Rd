\name{MTR_Graph}
\alias{MTR_Graph}

\title{
Multi-task regression with network structure
}

\description{
  This formulation constraints the models' relatedness according to the
  pre-defined graph \code{G}. If the penalty is heavy enough, the difference of
  connected tasks is0.
}

\usage{
MTR_Graph(X, Y, G, ...)
}

\arguments{
  \item{X}{
    a set of feature matrixes
}
  \item{Y}{
    a set of continuous responses
}
  \item{Lam1}{
    a positive number \eqn{\lambda_1} to control the network constraint with default of 1 
}
  \item{lam2}{
    a positive number \eqn{\lambda_2} to improve the
    generalization performance with default of 0
    
}
  \item{G}{
    a matrix to encode the network information. For more
  information, please refer to the \code{\link{Details}} section
}
  \item{opts}{
    options of the solver. The default is
    \code{list(init = 0, tol = 10^-3, maxIter = 1000)}
}
}

\details{
  \deqn{\min\limits_{W} \sum_{i}^{t}{\frac{1}{n_i}||Y_i- X_iW_i||^2 +
  \lambda_2||WG||_F^2 + \lambda_1||W||_F^2 }}

  \eqn{X} and \eqn{Y} are the sets of feature matrixes and continuous responses respectively,   
  \eqn{W} is the coefficient matrix, and
  \eqn{t} is the number of tasks. Accordingly, \eqn{Y_i}, \eqn{X_i}, \eqn{W_i}
  and \eqn{n_i} refer to the  data, model parameter set and
  the number of subjects for task \eqn{i}. Note \eqn{W_i} is the \eqn{i}th column of \eqn{W}, and
  \eqn{C_i} is the \eqn{i}th term of \eqn{C}. \eqn{||\circ||_F} is the
  Frobenius norm.  
   

  
  \eqn{G} represents the task-task similarity and needs to be specified by
  users. Intuitively, \eqn{||WG||_F^2} equals to an accumulation of
  differences bewteen pairwise tasks, i.e.
  \eqn{||WG||_F^2=\sum||W_{\mathcal{x}}-W_{\mathcal{y}}||_2^2}, where \eqn{\mathcal{x}} and 
  \eqn{\mathcal{y}} are connected tasks. Theoretically, 
  \eqn{||WG||_F^2} = \eqn{tr(WLW^T)}, where \eqn{L} is the graph Laplacian. Therefore,
  penalizing such term improves the task relatedness. 
 
  However, it is nontrivial to design \eqn{G}. Here, we give three common examples to demonstrate the construction of \eqn{G}.
    
  1) Assume your tasks are subject to orders i.e. temporal or spatial order. Such order forces the order-oriented
  smoothness across tasks such that the adjacent models(tasks) are similar, then your penalty can be designed as:
  \deqn{||WG||_F^2=\sum_i^{t-1}||W_{i}-W_{i+1}||_2^2}
  where the size of \eqn{G} is \eqn{(t+1) \times t}
  \deqn{G_{ij}=
    \left\{
      \begin{array}{ll} 
        1 & i=j\\
        -1 & i=j+1\\
        0 & otherwise
      \end{array}
    \right.
  }
  2) The so-called mean-regularized multi-task learning. In this formulation, each model is forced to approximate the mean of all models,
  \deqn{||WG||_F^2=\sum_i^{t}||W_{i}-\frac{1}{t}\sum_j^t{W_j}||_2^2}
  where the size of \eqn{G} is \eqn{t \times t}
  \deqn{G_{ij}=
    \left\{
      \begin{array}{ll} 
        \frac{t-1}{t} & i=j\\
        -\frac{1}{t} & others
      \end{array}
    \right.
  }
  3) Assume your tasks are related according to a given graph \eqn{g=(N,E)}
  where \eqn{E} is the edge set, and \eqn{N} is the node set. Then the
  penalty is
  \deqn{||WG||_F^2=\sum_i^{||E||}||W_{\alpha^i}-W_{\beta^i}||_2^2}
  where \eqn{\alpha \in N} and \eqn{\beta \in N} are connected tasks in
  the graph. 
  The size of \eqn{G} is \eqn{t \times ||E||}
  For each column i of G:
  \deqn{G_{ji} =
    \left\{
      \begin{array}{ll} 
        1 & j=\alpha^i\\
        -1 & j=\beta^i\\
	0 & otherwise
      \end{array}
    \right.
  }
  For more examples, the users are refered to the rich literatures of graph Laplacian.
}

\value{
  The function will return a trained MTR_Graph model
  \item{W}{a matrix of features' coefficients}
  \item{Obj}{historical record of objective values} 
  \item{fitted.values}{predictive values of the training data.}
  \item{residuals}{the residuals of the training data. For each subject
  \eqn{i}, the residual is \eqn{y_i-\hat{y_i}}}
  \item{lam1}{\eqn{\lambda_1} value}
  \item{lam2}{\eqn{\lambda_2} value}
  \item{G}{graph information}
  \item{opts}{options of the solver}
  \item{dim}{size of feature matrix of each task}
  \item{features}{feature names}
}

\references{
Evgeniou, T., & Pontil, M. (2004, August). Regularized multi-task
learning. In Proceedings of the tenth ACM SIGKDD international
conference on Knowledge discovery and data mining (pp. 109-117). ACM.

Widmer, C., Kloft, M., Görnitz, N., & Rätsch, G. (2012). Efficient
training of graph-regularized Multitask SVMs. Machine Learning and
Knowledge Discovery in Databases, 633-647.
}

\author{han.cao@zi-mannheim.de}

\seealso{
 \code{\link{MTC_Graph}}
 \code{\link{cv.MTR_Graph}}
}

\examples{
#load the data
####
#X, Y: training data
#tX, tY: test data
#W: ground truth
#G: the network information
####
load('./RMTL/data/Simulated_Regression_Graph.rda')

#speficy options
lam1=10^seq(2,-5, -0.05)

#cross-validation
cv <- cv.MTR_Graph(X, Y, G=G, lam1=lam1)

#training
r <- MTR_Graph(X, Y, G=G, lam1=cv$lam1.min) 

#predict on new dataset
predict(r, tX)

#compare the learnt model with the groud truth
par(mfrow=c(1,2))
library(fields)
image.plot(cor(W), xlab='tasks', ylab='tasks', main="ground truth")
image.plot(cor(r$W), xlab='tasks', ylab='tasks', main="Graph")

#extract more information about the model
print(r)
plotObj(r)
r$fitted.values
r$residuals
r$opts
r$features
r$dim
}

\keyword{ graph }
\keyword{ regression }
