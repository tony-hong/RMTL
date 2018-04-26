LS_CMTL <- function (X, Y, lam1, lam2, k, opts){

#------------------------------------------------------------
# private functions
bsa_ihb <- function(a,b,r,u){
# initilization
break_flag <- 0;
t_l <- a/b; t_u <- (a - u)/b;
T <- c(t_l, t_u)
t_L <- -Inf; t_U <- Inf;
g_tL <- 0; g_tU <- 0;

iter = 0
while (length(T)!=0){
    iter <- iter + 1
    g_t <- 0
    t_hat <- median(T)
    
    U <- t_hat < t_u
    M <- (t_u <= t_hat) & (t_hat <= t_l)

    if (sum(U)){
       g_t <- g_t + t(b[U]) %*% u[U] 
    }

    if (sum(M)){
        g_t <- g_t + sum(b[M] %*% (a[M] - t_hat * b[M]))
    }
    
    if (g_t > r){
        t_L <- t_hat
        T <- T[T > t_hat]
        g_tL = g_t
    } else if (g_t < r){
        t_U <- t_hat
        T <- T[T < t_hat]
        g_tU <- g_t
    } else{
        t_star <- t_hat
        break_flag <- 1
        break
    }
}
if (!break_flag){
     t_star <- t_L - (g_tL -r) * (t_U - t_L)/(g_tU - g_tL)     
}
temp <- sapply(a - t_star*b, function(x)max(0,x))
x_star <- sapply(temp, function(x) min(x,u))
return(list(x_star,t_star,iter))
}

singular_projection <- function(Msp, k){
requireNamespace('MASS')
eig <- eigen(Msp, symmetric=FALSE)
EVector <- eig$vector
EValue  <- eig$values
Pz <- Re(EVector)
diag_EValue <- Re(EValue)
DiagSigz <- bsa_ihb(diag_EValue, rep(1,length(diag_EValue)),
                    k, rep(1,length(diag_EValue)))
DiagSigz <- DiagSigz[[1]]
Mzp <- Pz %*% diag(DiagSigz) %*% t(Pz)
Mzp_Pz <- Pz
Mzp_DiagSigz <- DiagSigz
return(list(Mzp, Mzp_Pz, Mzp_DiagSigz))
}

gradVal_eval <- function (W, M){
requireNamespace('MASS')
requireNamespace('psych')
IM = (eta * diag(task_num) + M)
invIM <- MASS::ginv(IM)
invEtaMWt = invIM %*% t(W)

grad_W <- sapply(c(1:task_num),
    function(x) t(X[[x]]) %*% (X[[x]] %*% W[,x]-Y[[x]]) / nrow(X[[x]]))+
    2 * c * t(invEtaMWt)
grad_M = - c * (t(W) %*% W %*% invIM %*% invIM )      #M component

funcVal = sum(sapply(c(1:task_num),
    function(x)0.5 * mean((Y[[x]] - X[[x]] %*% W[,x])^2)))+
    c * psych::tr( W %*% invEtaMWt)
return(list(grad_W, grad_M, funcVal))
}

funVal_eval <- function (W, M_Pz, M_DiagSigz){
requireNamespace('psych')
invIM = M_Pz %*% (diag( 1/(eta + M_DiagSigz))) %*% t(M_Pz)
invEtaMWt = invIM %*% t(W);
return(sum(sapply(c(1:task_num),
    function(x)0.5 * mean((Y[[x]] - X[[x]] %*% W[,x])^2))) +
    c * psych::tr( W %*% invEtaMWt))
}
    
#------------------------------------------------------------

task_num <- length (X);
dimension = dim(X[[1]])[2];
Obj <- vector(); 

#initialize a starting point
if(opts$init==0){
   W0 <- matrix(0, nrow=dimension, ncol=task_num);
   M0 <- diag (task_num) * k / task_num;
}else if(opts$init==1){
   W0 <- opts$W0
   M0 <- opts$M0
}    

#precomputation
eta <- lam2 / lam1;
c <- lam1 * eta * (1 + eta);


bFlag <- 0; 
Wz <- W0;
Mz <- M0;
Wz_old <- W0;
Mz_old <- M0;

t <- 1;
t_old <- 0;
iter <- 0;
gamma <- 1;
gamma_inc <- 2;

while (iter < opts$maxIter){
    alpha <- (t_old - 1) /t;
    
    Ws <- (1 + alpha) * Wz - alpha * Wz_old;
    Ms <- (1 + alpha) * Mz - alpha * Mz_old;

    # compute function value and gradients of the search point
    r <- gradVal_eval(Ws, Ms);
    gWs <- r[[1]]
    gMs <- r[[2]]
    Fs <- r[[3]]
    
    # the Armijo Goldstein line search scheme
    while (TRUE){
        Wzp = Ws - gWs/gamma;
        r <- singular_projection (Ms - gMs/gamma, k);
        Mzp <- r[[1]]
        Mzp_Pz <- r[[2]]
        Mzp_DiagSigz <- r[[3]]
        Fzp <- funVal_eval  (Wzp, Mzp_Pz, Mzp_DiagSigz);
        
        delta_Wzp <- Wzp - Ws;
        delta_Mzp <- Mzp - Ms;
        r_sum <- (norm(delta_Wzp, 'f')^2 + norm(delta_Mzp, 'f')^2)/2;
        
        Fzp_gamma = Fs + sum(delta_Wzp* gWs) + 
            +sum(delta_Mzp * gMs)+
            +gamma * r_sum;
        
        if (r_sum <=1e-20){
            bFlag=1; 
            break;
        }
        
        if (Fzp <= Fzp_gamma) break else {gamma = gamma * gamma_inc}
  
    }
    
    Mz_old <- Mz;
    Wz_old = Wz;
    Wz = Wzp;
    Mz <- Mzp;
    
    Obj = c(Obj, Fzp);
    
    #test stop condition.
    if (bFlag) break;
    if (iter>=2){
        if (abs( Obj[length(Obj)] - Obj[length(Obj)-1] ) <= opts$tol)
            break;
    }
    
    iter = iter + 1;
    t_old = t;
    t = 0.5 * (1 + (1+ 4 * t^2)^0.5);

}

W = Wzp;
M = Mzp;

return(list(W=W, Obj=Obj, M=M))
}

MTR_CMTL <- function(X, Y, ...) UseMethod("MTR_CMTL")

MTR_CMTL.default <- function(X, Y, lam1=0.5, lam2=0.5,k=2,
             opts=list(init=0, tol=10^-3, maxIter=1000)){
    task_num <- length(X)
    features <- colnames(X[[1]])
    X <- lapply(X, function(x) as.matrix(x))
    Y <- lapply(Y, function(x) as.numeric(x))

    r <- LS_CMTL(X, Y, lam1, lam2, k, opts)
    r$fitted.values <- lapply(c(1:task_num),
        function(x) X[[x]] %*% r$W[,x])
    r$residuals <-lapply(c(1:task_num),
        function(x) Y[[x]] - r$fitted.values[[x]])
    r$call <- match.call()
    r$lam1 <- lam1
    r$lam2 <- lam2
    r$opts <- opts
    r$dim <- sapply(X, function(x)dim(x))
    r$features <- features
    r$k <- k
    class(r) <- "MTR_CMTL"
    return(r)
}

print.MTR_CMTL <- function(x)
{
    cat("\nHead Coefficients:\n")
    print(head(x$W))
    cat("Call:\n")
    print(x$call)
    cat("Formulation:\n")
    print('SUM_i Loss_i(W) + lam1*eta(1+eta)tr(W(eta*I+M)^{-1}W^T) \n s.t. tr(M)=k, M<=I, M \\in S{^t}{_+}, eta=lam2/lam1')
    
}

predict.MTR_CMTL <- function(m, newdata=NULL)
{
    if(is.null(newdata))
        Y <- m$fitted.values
    else{
        task_num <- length(newdata)
        y <- lapply(c(1:task_num),
        function(x) newdata[[x]] %*% m$W[,x])
    }
    return(y)
}

cv.MTR_CMTL <- function(X, Y,k=2,
             opts=list(init=0, tol=10^-3, maxIter=1000),
             nfolds=5,lam1=10^seq(3,-2, -1), lam2=10^seq(3,-2, -1)){
task_num <- length(X)
X <- lapply(X, function(x) as.matrix(x))
Y <- lapply(Y, function(x) as.numeric(x))

cvPar <- getCVPartition(X, Y, nfolds, FALSE)
cvm <- matrix(0, length(lam1), length(lam2))

#cv
for (i in 1:nfolds){
    cv_Xtr <- cvPar[[i]][[1]];
    cv_Ytr <- cvPar[[i]][[2]];
    cv_Xte <- cvPar[[i]][[3]];
    cv_Yte <- cvPar[[i]][[4]];

    
    for (i2 in 1: length(lam2)){
        cv_opt <- opts;
        for (i1 in 1: length(lam1)){
            r <- LS_CMTL(cv_Xtr, cv_Ytr, lam1[i1], lam2[i2],
                          k, cv_opt)
            cv_opt$init=1;
            cv_opt$W0=r$W;
            cv_opt$M0=r$M;
            cvm[i1, i2] = cvm[i1, i2]+mean(sapply(c(1:task_num),
                function(x)mean((cv_Xte[[x]] %*% r$W[,x]-cv_Yte[[x]])^2)))
        }
    }
}

cvm = cvm/nfolds
best_idx <- which(cvm==min(cvm), arr.ind=T)
cv <- list(lam1=lam1, lam2=lam2, lam1.min=lam1[best_idx[1]],
           lam2.min=lam2[best_idx[2]], cvm=cvm)
class(cv) <- "cv.MTR_CMTL"
return(cv)
}


plot.cv.MTR_CMTL <- function(x){
requireNamespace('fields')
par(oma=c( 0,0,0,5))
image(x$cvm, col=rev(heat.colors(128)), xlab='lambda1',
      ylab='lambda2', axes=FALSE, cex.lab = 1)
axis(1, at=seq(from=0,to=1,length=nrow(x$cvm)),
     labels=x$lam1,cex.axis=1)
axis(2, at=seq(from=0,to=1,length=ncol(x$cvm)),
     labels=x$lam2,cex.axis=1)
par(oma=c( 0,0,0,1))
fields::image.plot(x$cvm, col=rev(heat.colors(128)), legend.only=TRUE)
}
