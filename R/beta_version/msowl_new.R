# Multistate outcome weighted learning
# Bakoyannis Biometrika 2023 

# libraries 
library(survival)
library(kernlab)
library(pbapply)  # parallel
library(parallel)
library(matrixcalc)
#setwd("H:/msowl_hiv/programs")
#setwd("/Users/giorgosbakoyannis/Documents/MSOWL_HIV")
#dat <- read.csv("example_data.csv")

#Weighted SVM solver
#Code adapted from DTRlearn2 package (source code: https://github.com/ychen178/DTRlearn2/blob/master/R/wsvm_solve.R)
wsvm_solve <-function(X, A, wR, kernel='linear', sigma=0.05, lambda=1, e=1e-7) {
  if (kernel=='linear') {
    K = X %*% t(X)
    if (is.vector(X)) K = t(X) %*% X
  }
  else if (kernel=='rbf'){
    rbf = rbfdot(sigma = sigma)
    K = kernelMatrix(rbf, X)
  }
  y = A * sign(wR)
  H = y %*% t(y) * K
  n = length(A)
  C = 1/(2*n*lambda)

  # solve linear program
  solution <-
  tryCatch(
    ipop(
       c = rep(-1, n), 
       H = H, 
       A = t(y), 
       b = 0, 
       l = numeric(n), 
       u = C*abs(wR), 
       r = 0
      ), 
  error=function(er) er
  )
  
  if ("error" %in% class(solution)) {
    model <- NA
    class(model) <- "ipop function error"
  } else {
    alpha = primal(solution)
    alpha1 = alpha * y

    if (kernel=='linear'){
      w = t(X) %*% alpha1
      fitted = X %*% w
    } else if (kernel=='rbf'){
      fitted = K %*% alpha1
    }
    rm = y - fitted
    Imid = (alpha < C-e) & (alpha > e)
    rmid = rm[Imid==1]
    if (sum(Imid)>0){
      bias = mean(rmid)
    } else {
      Iup = ((alpha<e)&(A==-sign(wR)))|((alpha>C-e)&(A==sign(wR)))
      Ilow = ((alpha<e)&(A==sign(wR)))|((alpha>C-e)&(A==-sign(wR)))
      rup = rm[Iup]
      rlow = rm[Ilow]
      bias = (min(rup)+max(rlow))/2
    }

    if (kernel=='linear') {
      model = list(beta0=bias, beta=w, alpha1=alpha1)
      class(model)<-'linear'
    } else if (kernel=='rbf') {
      model = list(beta0=bias, sigma=sigma, Z=X, alpha1=alpha1)
      class(model) = 'rbf'
    }
  }
  return (model)
}

zeta_t <- function(t, data, w=c(0, 1, 0), fit, S, T){
  ##############################################################################
  # Utility function for the optimal treatment rule estimation function
  ##############################################################################
  # Calc the bracketed expr from the Value / Empirical Risk function 
  # \int_0^{\tau} { \frac{Y_w(t)I(C >= min(T,T))}{exp(-\Lambda_0(min(\tilde{T},t))} dm(t) } 
  # which calcs the Mean Patient-Weighted Multistate Survival Time, 
  # Inflated by Inverse Probability of Censoring Survival Times 
  ##############################################################################
  # t = tt[4] # just for debugging / playing around
  # baseline hazard function 
  hazard <- basehaz(fit, centered=FALSE)
  
  #Get the value that corresponds to the maximum t that is <= "time"
  element <- function(t){
    max(1,max(1:length(hazard$time)*(hazard$time<=t)))
  }
  
  # patient death indicator per state 
  events <- list()
  for(j in 1:length(S)){
    if(!(S[j] %in% T)){
      # not terminal state 
      events[[j]] <- (data$s1==S[j] & data$t1<=t & data$t2>t)
    } else {
      # terminal state 
      
      # boolean for each event: other events occur before particular event-time t 
      events[[j]] <- (data$s2==S[j] & data$t2<=t)
    }
  }
  # weight the events by patient-preference of state
  event <- do.call(cbind, events)%*%w

  Tt <- (1*!(data$s2 %in% T))*t + (1*(data$s2 %in% T))*mapply(min, data$t2, t) #T^t

  elements <- sapply(Tt, element)
  G_t <- exp(-hazard$hazard[elements])

  # inverse probability of censoring weighting ?? 
  res <- event/G_t
  
  return(res)
}

ms <- function(t1, t2, s1, s2, a){
  #if (missing(id))
  #  stop("Missing id argument")
  if (missing(t1))
    stop("Missing values in `t1'")
  if (missing(t2))
    stop("Missing values in `t2'")
  if (missing(s1))
    stop("Missing values in `s1'")
  if (missing(s2))
    stop("Missing values in `s2'")
  if (max(t1 >= t2)==1)
    stop("t1 >= t2")
  if (max(!(unique(a) %in% c(-1, 1)))==1)
    stop("Treatment var `a' is not equal to -1 or 1")
  msdat <- cbind(t1 = t1, t2 = t2, s1 = s1, s2 = s2, a = a)
  class(msdat) <- "msdat"
  msdat
}


calculate_reward = function(
    formula, 
    id, 
    w = c(0, 1, 0), 
    tau = 3, 
    data,
    owl_method = "MSOWL"
) {
  
  # print(data)
  
  # recover the structure of the data 
  msout <- model.frame(formula, data)[[1]]
  covs <- model.matrix(formula, data)[, -1]
  
  data <- cbind.data.frame(
    id = data$id,
    as.data.frame(cbind(msout, covs))
  )
  
  feat <- colnames(data)[(ncol(msout) + 2):ncol(data)]
  rm(msout, covs)
  data <- data[order(data[,id], data[,"t1"]),]
  
  # print("calculate_reward")
  # print("data")
  # print(data)
  
  
  
  ## state space
  S <- with(data, sort(unique(c(s1, s2))))
  
  ## absorbing state subspace
  T <- S[!(S %in% sort(unique(data$s1)))]
  
  # validate state patient-preference weights 
  if(length(S)!=length(w)){
    stop("length(w) not equal to the total number of states")
  }
  
  if(min(w)<0 | max(w)>1 | sum(w)>=length(w)){
    stop("min(w)<0 or max(w)>1 or sum(w)>=length(w)")
  }
  
  # Calculate reward
  # Estimate censoring distribution
  
  # censor status (variable often named "delta")
  c.delta <- 1*(data$s1==data$s2 & data$t2<tau)
  
  # propensity score 
  dat <- data[data$t1==0,]
  dat <- dat[order(dat[, id]),]
  pi_n <- mean(dat$a==1)
  
  # MSOWL INTEGRATES THE SURVIVAL OF EACH PATIENT OVER TIME UNTIL CENSOR TIME 
  # REQUIRES MORE COMPUTATIONALLY COSTLY ANALYSIS OF HAZARD 
  # USES FUNCTION zeta_t() 
  if(owl_method == "MSOWL") {
    
    # Censoring Distribution. 
    # nelson-aalen estimator for baseline hazard (via coxph model) 
    form <- as.formula("Surv(time=t1, time2=t2, event=c.delta, type='counting') ~ 1")
    fit <- coxph(form, data=data)
    
    # combine 
    # (1) imputed survival times from both censored and uncensored data
    # (2) true un-censored survival times 
    tt <- sort(
      unique(
        c(
          survfit(fit, se.fit = FALSE)$time, # 
          data[data$s1!=data$s2, "t2"]       # state transition times 
        )
      )
    )
    tt <- c(0,tt[tt<=tau])
    
    # zeta_t: compute hazard function (all pairwise combinations )
    event <- matrix(
      sapply(
          tt, 
          zeta_t, 
          data=data, 
          w=w, 
          fit=fit, 
          S=S, 
          T=T
          # cl=num_cores
        ),
        byrow=TRUE, 
        nrow=length(tt)
      )
    
    # per-individual integration 
    # hazard(t) = events(t) * dm(t)
    dm <- diff(c(tt, tau), lag=1) # dm 
    xi <- colSums(event*dm)  
    
    # this just transposes from row to column ... 
    xi <- matrix(tapply(X=xi, INDEX = data$id, FUN=sum))
    
  } else if (owl_method == "ICO") {
    # event_status = !c.delta
    event_status <- 1*(data$s1!=data$s2 & data$t2<tau)
    
    form <- as.formula("Surv(time=t1, time2=t2, event=c.delta, type='counting') ~ 1")
    fit <- coxph(form, data=data)
    hazard <- basehaz(fit, centered=FALSE) # cumulative hazard 
    
    print(length(event_status))
    print(length(hazard$hazard))
    print(length(data$t1))
    xi <- (event_status * (data$t2 - data$t1)) / exp(-hazard$hazard)
    
  } else if (owl_method == "DR") {
    # DR: W_i = R(Y_i, ind_cens_i, \hat{S_C},\hat{E_T})
    # R(Y,\delta,S_C,E_\tilde{E}) = \frac{\Delta Y}{S_C(Y|A,X)}-\int{E_\tilde{T}}
    
    # model survival function S_T(T|A,X) 
    # model E_T(T_tilde | T > t, A, X)
    # model S_C(t|A,X)  censoring distribution
    # W_i = R(Y_i, \Delta_i, S_C, E_T)
  }
  
  Wt <- xi/(pi_n*(dat$a==1) + (1 - pi_n)*(dat$a==-1))
  
  complete_data = data.frame(
    X <- as.matrix(dat[,feat]),
    A <- dat$a, 
    Y = dat$t2 - dat$t1, 
    c.delta = c.delta,
    Wt
  )
  
  retlist = list(
    dat, 
    complete_data,
    feat
  )
  
  return(retlist)
}


#Estimator of the value function of an ITR
#function(formula, id, w = c(0, 1, 0), tau = 3, data, lambda = 1){
msowl.val <- function(
    formula, 
    id, 
    w = c(0, 1, 0), 
    tau = 3, 
    data,
    rule, 
    fixed.rules = FALSE, 
    SE = FALSE, 
    nboot = 100,
    trim = NULL,
    reward = NULL
){
  # print(rule)
  # print(class(rule))
  # if( !(class(rule) %in%  c("msowl", "data.frame") ) ) 
    # print(paste(typeof(rule), class(rule)))
    # stop("`rule' should either be output of msowl() run, or a data.frame of {id, treatments}")
  if(!(fixed.rules %in% c(TRUE, FALSE)))
    stop("`fixed.rules' should be TRUE or FALSE")
  if(!(SE %in% c(TRUE, FALSE)))
    stop("`SE' should be TRUE or FALSE")
  if(!is.null(trim)){
    if(trim < 0 | trim > 0.5){
      stop("`trim' should be between >=0 and <= 0.5")
    }
  }
  V_d <- function(data){
  
    if (is.null(reward)){
      # print("calculating reward")
      cres = calculate_reward(
        formula, 
        id, 
        w = w, 
        tau = tau, 
        data = data,
        owl_method = "MSOWL"
      )
    } else {
      # print("using precalculated reward")
      cres = reward
    }
    
    # observed data
    dat = cres[[1]] # input data
    complete_data = cres[[2]] # weighted-outcome
    feat = cres[[3]] # features (covariates)
    # print("msowl366:dat")
    # print(dat)
    
    X <- as.matrix(dat[,feat]) # covariates 
    A <- dat$a                 # observed treatment 
    Wt <- complete_data$Wt     # outcome weights 
    
    # compute individualised treatments from rule
    if (class(rule) == "msowl"){
      if (rule$kernel == "linear"){
        form <- as.formula(paste("~", strsplit(as.character(rule$call), "~", fixed = TRUE))[[3]])
        # print("msowl377:form")
        # print(form)
        # print("dat")
        # print(dat)
        Z <- model.matrix(form, dat)
        f_n <- Z%*%rule$beta_opt
      } else if (rule$kernel == "rbf"){
        rbf <- rbfdot(sigma = rule$sigma)
        # kernel, new_covariates, training_covariates
        K <- kernelMatrix(rbf, as.matrix(dat[,feat]), rule$fit$Z)
        f_n <- rule$fit$beta0 + K%*%rule$fit$alpha1
      }
    } else if(class(rule) == "data.frame"){
      f_n <- rule$A_hat # data.frame [id, A_hat]
    } else{
      # print(paste("owl_method", owl_method))
      print("ERROR ITR NOT OF KNOWN CLASS")
      print(class(rule))
      return(NA)
    }
    
    # Value: compute the empirical expectation (mean)
    # of inverse-probability-of-censoring-weighted outcome 
    # for patients whose observed treatment matches 
    # the prescribed treatment from the itr estimator.
    V_dn <- mean(Wt*(A*f_n>=0))
    
    # Fixed One-Size-Fits-All Rules
    if(fixed.rules){
      f <- 1
      V_1 <- mean(Wt*(A*f>=0))
      f <- -1
      V_m1 <- mean(Wt*(A*f>=0))
      res <- c(V_dn, V_1, V_m1)
      names(res) <- c("V(dn)", "V(1)", "V(-1)")
    } else {
      res <- V_dn
      names(res) <- "V(dn)"
    }
    return(res)
  }
  
  # Compute Value of Treatment Rule on Population Sample
  est <- V_d(data = data)
  
  # Compute Standard Error of Value Estim. via Bootstrap procedure 
  if(SE){
    pb <- txtProgressBar(
      min = 0,      # Minimum value of the progress bar
      max = nboot,  # Maximum value of the progress bar
      style = 3,    # Progress bar style (also available style = 1 and style = 2)
      width = 50,   # Progress bar width. Defaults to getOption("width")
      char = "="    # Character used to create the bar
    )
    
    if(!is.character(data$id)){
      data$id <- as.character(data$id)
    }
    clusters <- sort(unique(data[,id]))
    bres <- matrix(NA, nrow=nboot, ncol=(fixed.rules*3 + !fixed.rules*1))
    for(b in 1:nboot){
      index <- sample(1:length(clusters), length(clusters), replace=TRUE)
      cl <- clusters[index]
      bdat <- NULL
      for(j in 1:length(cl)){
        aa <- data[data[,id] %in% cl[j],]
        reps <- table(cl[1:j])[names(table(cl[1:j]))==cl[j]] - 1
        if(reps > 0){
          aa[,id] <- paste(aa[1,id], 1*reps, sep = ".")
        }
        bdat <- rbind(bdat, aa)
      }
      bVal <- try(V_d(data = bdat), silent = TRUE)
      if(class(bVal) != "try-error"){
        bres[b,] <- bVal
      }
      setTxtProgressBar(pb, b)
    }
    close(pb) # Close the connection

    if(fixed.rules){
      nfail <- sum(is.na(bres[,1]))
      SE <- apply(bres, 2, sd, na.rm = TRUE)

      res <- cbind(est, SE, est - qnorm(0.975)*SE, est + qnorm(0.975)*SE)
      colnames(res) <- c("estimate", "SE", "ll", "ul")
      res <- round(res, 3)

      C <- rbind(c(1, -1, 0),
                 c(1, 0, -1))
      est <- C%*%est
      SE <- sqrt(diag(C%*%var(bres, na.rm = TRUE)%*%t(C)))
      res2 <- cbind(est, est - qnorm(0.975)*SE, est + qnorm(0.975)*SE)
      colnames(res2) <- c("estimate", "ll", "ul")
      rownames(res2) <- c("V(dn) - V(1)", "V(dn) - V(-1)")
      res2 <- round(res2, 3)

      res <- list(res, res2, n.failed = nfail)
    } else {
      nfail <- sum(is.na(bres[,1]))
      SE <- sd(bres, na.rm = TRUE)
      res <- c(est, SE, est - qnorm(0.975)*SE, est + qnorm(0.975)*SE)
      names(res) <- c("estimate", "SE", "ll", "ul")
      res <- list(round(res, 3), n.failed = nfail)
    }
  } else {
    res <- est
  }

  return(res)
}

msowl <- function(
    formula, 
    id, 
    w = c(0, 1, 0), 
    tau = 3, 
    data, 
    lambda = 1,
    kernel="linear",
    sigma=5,
    jackknife = FALSE, 
    trim = NULL,
    owl_method = "MSOWL",
    reward = NULL,
    debug=T
){
  
  if (is.null(reward)){
    cres = calculate_reward(
      formula, 
      id, 
      w = w, 
      tau = tau, 
      data = data,
      owl_method = owl_method
    )
  } else {
    # print("using precalculated reward")
    cres = reward
  }
  
  dat = cres[[1]] 
  complete_data = cres[[2]]
  feat = cres[[3]]
  
  X <- as.matrix(dat[,feat])
  A <- dat$a
  Wt <- complete_data$Wt
  
  # obtain \hat{f}(x) by minimizing 
  # sum_i=1^n W_i \frac{\phi(A_i f(X_i))}{\pi(A_i;X_i)} + \lambda_n ||f||^2
  fit <- wsvm_solve(
    X=X, 
    A=A, 
    wR=Wt, 
    kernel=kernel, 
    sigma=sigma, 
    lambda=lambda
  )

  # exit if error 
  if(unlist(gregexpr('error', class(fit)))[1] !=-1){
      res <- NA
      class(res) <- "error"
      if (debug){
        print("error!")
        print(fit) 
      }
      return(res)
  }
  
  if (kernel == "linear"){
    
    # format linear itr parameters
    beta_opt <- c(fit$beta0, fit$beta)
    names(beta_opt) <- c("(constant)", feat)
    
    # compute treatments for the given population from learned itr
    f_n <- fit$beta0 + X%*%fit$beta
    
    # value (empirical expectation of patients whose 
    # observed treatment corresponds with that assigned by learned ITR
    V_dn <- mean(Wt*(A*f_n>=0))
    
    # format output 
    res <- list(
      beta_opt = beta_opt, 
      Value = V_dn, 
      fit = fit, 
      call = formula,
      kernel=kernel
    )
    
  } else if (kernel=="rbf"){
    # compute treatments for given population from learned itr 
    rbf <- rbfdot(sigma = sigma)
    K <- kernelMatrix(rbf, as.matrix(dat[,feat]), fit$Z)
    f_n <- fit$beta0 + K%*%fit$alpha1
    V_dn <- mean(Wt*(A*f_n>=0))
    
    # format output 
    res <- list(
      # beta_opt = beta_opt, 
      Value = V_dn, 
      fit = fit, 
      call = formula,
      kernel=kernel,
      sigma=sigma
    )
  }

  # Variance Estimator using JackKnife method 
  # Only for Linear Rules ... 
  if( kernel == "linear" & isTRUE(jackknife)){
    ids <- sort(unique(dat$id))
    f <- rep(NA, times=nrow(dat))
    pb <- txtProgressBar(
      min = 0,      # Minimum value of the progress bar
      max = length(ids),  # Maximum value of the progress bar
      style = 3,    # Progress bar style (also available style = 1 and style = 2)
      width = 50,   # Progress bar width. Defaults to getOption("width")
      char = "=")   # Character used to create the bar
    
    # leave-one-out 
    for(i in ids){
      fit1 <- wsvm_solve(
        X=X[dat$id!=i,], 
        A=A[dat$id!=i], 
        wR=Wt[dat$id!=i],
        kernel='linear', 
        lambda=lambda
      )
      
      if(unlist(gregexpr('error', class(fit1)))[1]==-1){
        f[dat$id==i] <- fit1$beta0 + X[dat$id==i,]%*%fit1$beta
      }
      setTxtProgressBar(pb, i)
    }
    close(pb) # Close the connection
    
    # V_jack <- mean(Wt*(A*f>=0))
    ### Sean code added. 
    # Issue with wsvm_solve returning NAs for some jacknifed samples (soln: rm NA from output)
    vecVal=Wt*(A*f>=0)
    V_jack <- mean(vecVal,na.rm=TRUE)
    SD_jack=sd(vecVal,na.rm=TRUE)
    
    nMiss_jack=sum(is.na(vecVal))
    
    res <- c(
      res, 
      list(
        Value_jack = V_jack,
        SD_jack=SD_jack,
        nMiss_jack=nMiss_jack
      )
    )
  }
  
  class(res) <- "msowl"
  
  return(res)
}

# Select the optimal lambda through leave-one-out cross-validation
select.lambda <- function(
    formula, 
    id, 
    w, 
    tau, 
    data, 
    lambdas = c(0.01, 0.1, 0.5, 1, 5, 10, 50, 100)
){
  V_jack <- rep(NA,length(lambdas))
  SD_jack <- rep(NA,length(lambdas))
  nMiss_jack <- rep(NA,length(lambdas))

  for(i in 1:length(lambdas)){
    print(paste("lambda", lambdas[i], sep = " = "))
  
    fit <- msowl(
      formula, 
      id = id, 
      w = w, 
      tau = tau, 
      data = data,
      lambda = lambdas[i], 
      jackknife = TRUE
    )
    
    if(unlist(gregexpr('error', class(fit)))[1]==-1){
      print(fit$Value_jack, fit$SD_jack, fit$nMiss_jack)
      V_jack[i] <- fit$Value_jack
      SD_jack[i] <- fit$SD_jack
      nMiss_jack <- fit$nMiss_jack
    } else{
      print("error ... :( ")
    }
  }
  #Should pick lambda better with standard deviation as well
  lambda <- min(lambdas[V_jack==max(V_jack, na.rm = TRUE)], na.rm = TRUE)
  res <- list(
    lambda = lambda, 
    details = cbind(lambdas, V_jack,SD_jack,nMiss_jack)
  )
  return(res)
}


#### Examples ####

#select.lambda(ms(t1 = t1, t2 = t2, s1 = s1, s2 = s2, a = A) ~ Z1 + Z2,
#                 id = "id", w = c(0, 1, 0), tau = 3, data = dat,
#                 lambdas = c(0.01, 0.1, 0.5, 1, 5, 10, 50, 100))

#fit <- msowl(ms(t1 = t1, t2 = t2, s1 = s1, s2 = s2, a = A) ~ Z1 + Z2,
#             id = "id", w = c(0, 1, 0), tau = 3, data = dat,
#             lambda = 1)

#Value <- msowl.val(ms(t1 = t1, t2 = t2, s1 = s1, s2 = s2, a = A) ~ Z1 + Z2,
#                   id = "id", w = c(0, 1, 0), tau = 3, rule = fit, fixed.rules = TRUE, data = dat)


# ----------------------------------------------
# Utility Functions for MSOWL model
# ----------------------------------------------
# Select the optimal lambda through leave-one-out cross-validation
select.lambda.new <- function(
    formula,
    id,
    w,
    tau,
    data,
    lambdas = c(0.01, 0.1, 0.5, 1, 5, 10, 50, 100),
    debug=F
){
  if (debug) print("select.lambda.new")
  V <- rep(NA,length(lambdas))
  # SD_jack <- rep(NA,length(lambdas))
  # nMiss_jack <- rep(NA,length(lambdas))
  
  for(i in 1:length(lambdas)){
    if (debug) print(paste("lambda", lambdas[i], sep = " = "))
    
    fit <- msowl(
      formula,
      id = id,
      w = w,
      tau = tau,
      data = data,
      lambda = lambdas[i],
      jackknife = F
    )
    if (debug) print(fit$Value)
    
    if(unlist(gregexpr('error', class(fit)))[1]==-1){
      
      V[i] <- fit$Value
      # V_jack[i] <- fit$Value_jack
      # SD_jack[i] <- fit$SD_jack
      # nMiss_jack <- fit$nMiss_jack
    }
  }
  #Should pick lambda better with standard deviation as well
  # lambda <- min(lambdas[V_jack==max(V_jack, na.rm = TRUE)], na.rm = TRUE)
  lambda <- min(lambdas[V==max(V,na.rm=T)])
  res <- list(
    lambda = lambda,
    details = cbind(lambdas, V)
    # details = cbind(lambdas, V_jack,SD_jack,nMiss_jack)
  )
  return(res)
}
# ----------------------------------------------
# tune lambda parameter
# ----------------------------------------------
# lambda_error = 0
# tuned_lambda <-
#   tryCatch(
#     # select.lambda.new(
#     select.lambda.new(
#       my_formula,
#       id = "id",
#       w = c(1,0),
#       tau = tau,
#       data = data.train,
#       lambdas = c(0.01, 0.05, 0.1, 0.5, 1, 5, 10, 50, 100)
#       # lambdas = c(0.01, 0.1, 1, 10)
#   ),
#   error=function(er) er
# )
# return(tuned_lambda)
# 
# if ("error" %in% class(tuned_lambda)) {
#   print("*--------------------------------*")
#   print("ERROR IN LAMBDA TUNING!")
#   print(my_input)
#   print("*--------------------------------*")
#  tuned_lambda = 0.1  
#  tuned_lambda$lambda = 0.1
#  lambda_error = 1
# }
# if(debug)print(tuned_lambda)



precalculate_reward = function(
    data,
    n,
    my_formula,
    tau,
    n.batch,
    owl_method="MSOWL",
    debug=F
) {
  # print("# -----------------------------------")
  if (debug) print("pre-calculate reward in small batches")
  
  n_groups = n / n.batch
  
  data.batches = split(
    as.data.frame(data),
    sample.int(
      n=n_groups,
      size=nrow(data),
      replace=TRUE
    )
  )
  
  cres = pbapply::pbsapply(
    # cres = sapply(
    X=data.batches,
    FUN=calculate_reward,
    formula = my_formula,
    id = "id",
    w = c(1,0),
    tau = tau,
    owl_method=owl_method
  )
  # print(cres)
  
  res = list(
    dat = bind_rows(cres[1,]),
    complete_data = bind_rows(cres[2,]),
    feat = cres[3,][[1]]
  )
  
  return(res)
}

# ----------------------------------
# Divide & Conquer Value Function
# ----------------------------------
# n_groups = n.test / 200
# 
# data.test.batches = split(
#   as.data.frame(data.test),
#   sample.int(
#     n=n_groups,
#     size=nrow(data.test),
#     replace=TRUE
#   )
# )
# 
# msowl_test = pbapply::pbsapply(
#   data.test.batches,
#   msowl.val,
#   formula=my_formula,
#   id="id",
#   w = c(1,0),
#   tau=tau,
#   rule=msowl_train,
#   fixed.rules=T
# )
# 

# ----------------------------------------------
# DIVIDE & CONQUER for COX Regression Model
# ----------------------------------------------
# cox_dq_helper = function(
#     data,
#     cox_fit,
#     formula=my_formula,
#     id="id",
#     w = c(1,0),
#     tau=tau
# ){
#   cox.test.predictions <- coxph_itr_predict(
#     cox_fit = cox_fit,
#     data = data,
#     tau=tau
#   )
#   
#   cval = msowl.val(
#     data=data,
#     formula=my_formula,
#     id="id",
#     w = c(1,0),
#     tau=tau,
#     rule=cox.test.predictions,
#     fixed.rules=F
#   )
# }
# 
# cox.test.val = pbapply::pbsapply(
#   X=data.test.batches,
#   FUN=cox_dq_helper,
#   cox_fit = cox.train,
#   formula=my_formula,
#   id="id",
#   w = c(1,0),
#   tau=tau
# )
# ############################################