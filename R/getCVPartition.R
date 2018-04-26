getCVPartition <- function(X, Y, cv_fold, stratify){
task_num = length(X);
cvPar = {};
for (cv_idx in 1: cv_fold){
    # buid cross validation data splittings for each task.
    cv_Xtr = {};
    cv_Ytr = {};
    cv_Xte = {};
    cv_Yte = {};

    if(stratify)
    #stratified cross validation
    for (t in 1: task_num){
        task_sample_size <- length(Y[[t]]);
    
        ct <- which(Y[[t]]<1);
        cs <- which(Y[[t]]>=1);
        ct_idx <- seq(cv_idx, length(ct), cv_fold);
        cs_idx <- seq(cv_idx, length(cs), cv_fold);
        te_idx <- c(ct[ct_idx], cs[cs_idx]);
        tr_idx <- seq(1,task_sample_size)[
            !is.element(1:task_sample_size, te_idx)];

        cv_Xtr[[t]] = X[[t]][tr_idx,]
        cv_Ytr[[t]] = Y[[t]][tr_idx]
        cv_Xte[[t]] = X[[t]][te_idx,]
        cv_Yte[[t]] = Y[[t]][te_idx]
    }
   else{
    for (t in 1: task_num){
        task_sample_size <- length(Y[[t]]);
    
        te_idx <- seq(cv_idx, task_sample_size, by=cv_fold)
        tr_idx <- seq(1,task_sample_size)[
            !is.element(1:task_sample_size, te_idx)];

        cv_Xtr[[t]] = X[[t]][tr_idx,]
        cv_Ytr[[t]] = Y[[t]][tr_idx]
        cv_Xte[[t]] = X[[t]][te_idx,]
        cv_Yte[[t]] = Y[[t]][te_idx]
    }
   }
    
    cvPar[[cv_idx]]=list(cv_Xtr, cv_Ytr, cv_Xte, cv_Yte);
}
return(cvPar)
}

