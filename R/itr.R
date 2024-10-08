# Estimation of an optimal individualized treatment rule
library(survival)
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# print(getwd())
# source("./zeta_t.R")
# source("./wsvm_solve.R")
# source("./V_d.R")

itr <- function(
    data, 
    feat, 
    w=c(0, 1, 0), 
    tau=3, 
    kernel='linear', 
    lambda=1, 
    sigma=1, 
    SE=TRUE
){
  
  # TODO generalise to other kernels ...
  if(!(kernel %in% c('linear', 'rbf'))){
    stop("Only 'linear' and 'rbf' (radial basis function or Gaussian) kernels are supported")
  }
  ## state space
  S <- with(data, sort(unique(c(s1, s2))))
  
  ## absorbing state subspace
  T <- S[!(S %in% sort(unique(data$s1)))]
  
  # ill-formed state-weights
  if(length(S)!=length(w)){
    stop("length(w) not equal to the total number of states")
  }
  
  # invalid 
  if(max(data$t1>data$t2)==1){
    stop("t1 > t2")
  }
  
  # ill-formed state-weights
  if(min(w)<0 | max(w)>1 | sum(w)>=length(w)){
    stop("min(w)<0 or max(w)>1 or sum(w)>=length(w)")
  }
  
  #Calculate reward
  
  #Estimate censoring distribution 
  c.delta <- 1*(data$s1==data$s2 & data$t2<tau)
  fit <- coxph(Surv(time=t1, time2=t2, event=c.delta,type="counting") ~ 1, data=data)
  
  # predicted survival curve for cox model
  Haz <- basehaz(fit, centered=FALSE) 
  
  # pseudo-observations = merge(predicted_from_censored_times, uncensored_event_times)
  # print(Haz$time)
  # print(data[data$s1!=data$s2,"t2"])
  tt <- sort(unique(c(Haz$time, data[data$s1!=data$s2,"t2"]))) # uncensored 
  tt <- c(0,tt[tt<=tau]) # all pseudo-observed event-times, ordered
  
  
  # zeta_t
  event <- matrix(
    pbapply::pbsapply(
      tt, 
      zeta_t, 
      data=data, 
      w=w, 
      hazard=Haz, 
      S=S, 
      T=T
    ), 
    byrow=TRUE, 
    nrow=length(tt)
  )
  
  tt <- sort(unique(c(survfit(fit, se.fit = FALSE)$time,data[data$s1!=data$s2, "t2"])))
  tt <- c(0,tt[tt<=tau])
  event <- matrix(
    sapply(tt, zeta_t, data = data, w = w, hazard= Haz, S = S, T = T),
    byrow = TRUE, 
    nrow = length(tt)
  )
  
  dm <- diff(c(tt, tau), lag=1)
  xi <- colSums(event*dm)
  xi <- matrix(tapply(X=xi, INDEX = data$id, FUN=sum))
  
  # propensity score 
  dat <- data[data$t1==0,]
  dat <- dat[order(dat$id),]
  pi_n <- mean(dat$A==1)
  
  # Outcome Weights (IPW)
  Wt <- xi/(pi_n*(dat$A==1) + (1 - pi_n)*(dat$A==-1))
  
  # Covariates 
  X <- as.matrix(dat[,feat])
  
  # treatment
  A <- dat$A
  
  # Value Estimation
  if(kernel=='linear'){
    fit <- wsvm_solve(X=X, A=A, wR=Wt, kernel='linear', lambda=lambda)
    # print(summary(fit))
    
    if(!is.na(fit$beta0)){
      beta_opt <- c(fit$beta0, fit$beta)
      V_opt <- V_d(data=data, w=w, tau=tau, dec.fun=fit, feat=feat, SE=SE)
      res <- list(
        beta_opt=beta_opt,
        V_opt=V_opt$V_n,
        se_V_opt=V_opt$se,
        fit=fit
      )
      # class(res) = "linear"
    } else {
      res <- NA
      class(res) <- "error"
    }
  } else {
    fit <- wsvm_solve(X=X, A=A, wR=Wt, kernel='rbf', sigma=sigma, lambda=lambda)
    if(!is.na(fit$beta0)){
      V_opt <- V_d(data=data, w=w, tau=tau, dec.fun=fit, feat=feat, SE=FALSE)
      res <- list(V_opt=V_opt, fit=fit)
      class(res) = "rbf"
    } else {
      res <- NA
      class(res) <- "error"
    }
  }
  return(res)
}


