LR_Lasso <- function (X, Y, lam1, lam2, opts){
#------------------------------------------------------------
# private functions
l1_projection <- function (W, lambda ){
p <- abs(W) - lambda/2
p[p<0] <- 0
Wp <- sign(W) * p
return(Wp)
}

gradVal_eval <- function (W, C){
r <- lapply(c(1:task_num),
    function(x)unit_grad_eval( W[, x], C[x], X[[x]], Y[[x]]))
grad_W <- sapply(r, function(x)x[[1]]) + 2* lam2 * W
grad_C <- sapply(r, function(x)x[[2]])
funcVal = sum(sapply(r, function(x)x[[3]])) + lam2 * norm(W, 'f')^2
return(list(grad_W, grad_C, funcVal))
}    

funVal_eval <- function (W, C){
return(sum(sapply(c(1:task_num),
    function(x)unit_funcVal_eval(W[, x], C[x], X[[x]], Y[[x]]))) +
    lam2 * norm(W, 'f')^2)
}
    
nonsmooth_eval <- function (W, lam1){
return(lam1*sum(abs(W)))
}

unit_grad_eval <- function( w, c, x, y){
#gradient and logistic evaluation for each task
weight <- 1/length(y)
l <- -y*(x %*% w + c)
lp <- l
lp[lp<0] <- 0
funcVal <- sum(weight * ( log( exp(-lp) +  exp(l-lp) ) + lp ))
b <- (-weight*y)*(1 - 1/ (1+exp(l)))
grad_c <- sum(b)
grad_w <- t(x) %*% b
return(list(grad_w, grad_c, funcVal))
}

unit_funcVal_eval <- function ( w, c, x, y){
#function value evaluation for each task 
weight <- 1/length(y)
l <- -y*(x %*% w + c)
lp <- l
lp[lp<0] <- 0
return(sum(weight * ( log( exp(-lp) +  exp(l-lp) ) + lp )))
}
    
#------------------------------------------------------------
    
task_num <- length (X);
dimension = dim(X[[1]])[2];
subjects <- dim(X[[1]])[1];
Obj <- vector(); 

#initialize a starting point
if(opts$init==0){
   W0 <- matrix(0, nrow=dimension, ncol=task_num);
   C0 <- rep(0, task_num);
}else if(opts$init==1){
   W0 <- opts$W0
   C0 <- opts$C0
}

bFlag <- 0; 
Wz <- W0;
Cz <- C0;
Wz_old <- W0;
Cz_old <- C0;

t <- 1;
t_old <- 0;
iter <- 0;
gamma <- 1;
gamma_inc <- 2;

while (iter < opts$maxIter){
    alpha <- (t_old - 1) /t;
    
    Ws <- (1 + alpha) * Wz - alpha * Wz_old;
    Cs <- (1 + alpha) * Cz - alpha * Cz_old;
    
    # compute function value and gradients of the search point
    r <- gradVal_eval(Ws, Cs);
    gWs <- r[[1]]
    gCs <- r[[2]]
    Fs <- r[[3]]


    # the Armijo Goldstein line search scheme
    while (TRUE){
        Wzp <- l1_projection(Ws - gWs/gamma, 2 * lam1 / gamma);
        Czp <- Cs - gCs/gamma;
        Fzp <- funVal_eval  (Wzp, Czp);
        
        delta_Wzp <- Wzp - Ws;
        delta_Czp <- Czp - Cs;
        nrm_delta_Wzp <- norm(delta_Wzp, 'f')^2;
        nrm_delta_Czp <- sum(delta_Czp * delta_Czp);
        r_sum <- (nrm_delta_Wzp+nrm_delta_Czp)/2;
        
        Fzp_gamma = Fs + sum(delta_Wzp* gWs) + 
            sum(delta_Czp * gCs) + gamma/2 * nrm_delta_Wzp +
            gamma/2 * nrm_delta_Czp;
        
        if (r_sum <=1e-20){
            bFlag=1; 
            break;
        }
        if (Fzp <= Fzp_gamma) break else {gamma = gamma * gamma_inc}
  
    }
    
    Wz_old = Wz;
    Cz_old = Cz;
    Wz = Wzp;
    Cz = Czp;
    Obj = c(Obj, Fzp + nonsmooth_eval(Wz, lam1));
    
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
C = Czp;

return(list(W=W, C=C, Obj=Obj))
}


MTC_Lasso <- function(X, Y, ...) UseMethod("MTC_Lasso")


MTC_Lasso.default <- function(X, Y, lam1=0.01, lam2=0,
             opts=list(init=0,  tol=10^-3, maxIter=1000)){
    task_num <- length(X)
    X <- lapply(X, function(x) as.matrix(x))
    Y <- lapply(Y, function(x) as.numeric(x))

    r <- LR_Lasso(X, Y, lam1, lam2, opts)
    fitted <- lapply(c(1:task_num),
        function(x) exp(X[[x]] %*% r$W[,x] + r$C[x]))
    r$fitted.values <- lapply(fitted, function(x) x/(1+x))
    r$residuals <-lapply(c(1:task_num),
        function(x) Y[[x]]-(round(r$fitted.values[[x]])-0.5)*2)
    r$call <- match.call()
    r$lam1 <- lam1
    r$lam2 <- lam2
    r$opts <- opts
    r$dim <- sapply(X, function(x)dim(x))
    r$features <- colnames(X[[1]])
    class(r) <- "MTC_Lasso"
    return(r)
}

print.MTC_Lasso <- function(x)
{

    cat("\nHead Coefficients:\n")
    print(head(x$W))
    cat("Call:\n")
    print(x$call)
    cat("Formulation:\n")
    print('SUM_i Loss_i(W) + lam1*||W||_1 + lam2*||W||{_2}{^2}')
}

predict.MTC_Lasso <- function(m, newX=NULL)
{
    if(is.null(newX))
        y <- m$fitted.values
    else{
        task_num <- length(newX)
        y <- lapply(c(1:task_num),
            function(x) exp(newX[[x]] %*% m$W[,x] + m$C[x]))
        y <- lapply(y, function(x) x/(1+x))
    }
    return(y)
}

cv.MTC_Lasso <- function(X, Y, lam2=0,
             opts=list(init=0,  tol=10^-3, maxIter=1000),
             stratify=FALSE, nfolds=5,lam1=10^seq(1,-5, -1)){
task_num <- length(X)
X <- lapply(X, function(x) as.matrix(x))
Y <- lapply(Y, function(x) as.numeric(x))

cvPar <- getCVPartition(X, Y, nfolds, stratify)
cvm <- rep(0, length(lam1));

#cv
for (i in 1:nfolds){
    cv_Xtr <- cvPar[[i]][[1]];
    cv_Ytr <- cvPar[[i]][[2]];
    cv_Xte <- cvPar[[i]][[3]];
    cv_Yte <- cvPar[[i]][[4]];

    cv_opt <- opts;
    for (p_idx in 1: length(lam1)){
        r <- LR_Lasso(cv_Xtr, cv_Ytr, lam1[p_idx], lam2, cv_opt)
        cv_opt$init=1;
        cv_opt$W0=r$W;
        cv_opt$C0=r$C;
        cv_fitted <- lapply(c(1:task_num),
            function(x) exp(cv_Xte[[x]] %*% r$W[,x] + r$C[x]))
        cv_fitted <- lapply(cv_fitted, function(x) x/(1+x))
        cv_err<-mean(sapply(c(1:task_num),
            function(x) {mean(cv_Yte[[x]]!=(round(cv_fitted[[x]])-0.5)*2)}))
        cvm[p_idx] = cvm[p_idx]+cv_err
    }
}
cvm = cvm/nfolds
best_idx <- which(cvm==min(cvm))[1]
cv <- list(lam1=lam1, lam1.min=lam1[best_idx], lam2=lam2, cvm=cvm)
class(cv) <- "cv.MTC_Lasso"
return(cv)
}


plot.cv.MTC_Lasso <- function(x){
plot(log10(x$lam1), x$cvm, xlab="log10(lambda1)",
     ylab="error")
}
