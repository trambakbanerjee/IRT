assess.func <- function(zv,q,ns=10000,delta=0.1,ucb=TRUE){
  m <- length(zv)
  pi <- epsest.func(zv,0,1)
  dens.func <- function(zv,data,h) mean(dnorm((zv-data)/h))/h
  h.zv <- density(zv,from = min(zv)-10,to=max(zv)+10,n=m)$bw
  dens.est <- sapply(zv, dens.func, data=zv, h=h.zv)
  lfdr.est <- (1-pi)*dnorm(zv)/dens.est
  
  zv.null <- rnorm(ns)
  dens.null <- sapply(zv.null, dens.func, data=zv, h=h.zv)
  grid.null <- (1-pi)*dnorm(zv.null)/dens.null
  lfdr.cdf <- sapply(lfdr.est, function(x) mean(grid.null<=x))
  
  lfdr.sort <- sort(lfdr.est)
  if (ucb){
    ucb.bern <- function(p) max(sqrt(p*(1-p)*(1-pi)*m*log(1/delta)),1)
    fdp.est <- sapply(lfdr.sort, 
                      function(x) (ucb.bern(mean(grid.null<=x))+sum(lfdr.cdf>=1-mean(grid.null<=x)))/max(1,sum(lfdr.est<=x)))
  }else{
    fdp.est <- sapply(lfdr.sort, function(x) (1+sum(lfdr.cdf>=1-mean(grid.null<=x)))/max(1,sum(lfdr.est<=x)))
  }
  thr <- lfdr.sort[max(which(fdp.est<=q))]
  de <- rep(0,m); de[which(lfdr.est<=thr)]<-1
  estat <- m*de/(1+sum(lfdr.cdf>=1-mean(grid.null<=thr)))
  
  return(list(lfdr=lfdr.est, cdf=lfdr.cdf, cdf.func=grid.null, thr=thr, de=de, ev=estat))
}

assess.BH.func <- function(pv,q,mirror=TRUE){
  m <- length(pv)
  pv.sort <- sort(pv)
  pvi <- pv.sort/1:m
  thr <- pv.sort[max(which(pvi<=q/m))]
  de <- rep(0,m); de[which(pv<=thr)]<-1
  if (mirror){
    estat <- m*de/(1+sum(pv>=1-thr))
  }else{
    estat <- de/thr
  }
  
  return(list(thr=thr, de=de, ev=estat))
}

assess.Lfdr.func <- function(zv,q){
  m <- length(zv)
  pi <- epsest.func(zv,0,1)
  dens.func <- function(zv,data,h) mean(dnorm((zv-data)/h))/h
  h.zv <- density(zv,from = min(zv)-10,to=max(zv)+10,n=m)$bw
  dens.est <- sapply(zv, dens.func, data=zv, h=h.zv)
  lfdr.est <- (1-pi)*dnorm(zv)/dens.est
  lfdr.sort <- sort(lfdr.est)
  fdp.est <- sapply(lfdr.sort, function(x) (sum(lfdr.est*(lfdr.est<=x)))/max(1,sum(lfdr.est<=x)))
  thr <- lfdr.sort[max(which(fdp.est<=q))]
  de <- rep(0,m); de[which(lfdr.est<=thr)]<-1
  estat <- m*de/(sum(lfdr.est*(lfdr.est<=thr)))
  
  return(list(lfdr=lfdr.est, thr=thr, de=de, ev=estat))
}

assess.t.func <- function(zv,df,q,ns=10000,delta=0.1,ucb=TRUE){
  m <- length(zv)
  pi <- epsest.func(zv,0,1)
  dens.func <- function(zv,data,h) mean(dnorm((zv-data)/h))/h
  h.zv <- density(zv,from = min(zv)-10,to=max(zv)+10,n=m)$bw
  dens.est <- sapply(zv, dens.func, data=zv, h=h.zv)
  lfdr.est <- (1-pi)*dt(zv,df)/dens.est
  
  zv.null <- rt(ns,df)
  dens.null <- sapply(zv.null, dens.func, data=zv, h=h.zv)
  grid.null <- (1-pi)*dt(zv.null,df)/dens.null
  lfdr.cdf <- sapply(lfdr.est, function(x) mean(grid.null<=x))
  
  lfdr.sort <- sort(lfdr.est)
  if (ucb){
    ucb.bern <- function(p) max(sqrt(p*(1-p)*(1-pi)*m*log(1/delta)),1)
    fdp.est <- sapply(lfdr.sort, 
                      function(x) (ucb.bern(mean(grid.null<=x))+sum(lfdr.cdf>=1-mean(grid.null<=x)))/max(1,sum(lfdr.est<=x)))
  }else{
    fdp.est <- sapply(lfdr.sort, function(x) (1+sum(lfdr.cdf>=1-mean(grid.null<=x)))/max(1,sum(lfdr.est<=x)))
  }
  thr <- lfdr.sort[max(which(fdp.est<=q))]
  de <- rep(0,m); de[which(lfdr.est<=thr)]<-1
  estat <- m*de/(1+sum(lfdr.cdf>=1-mean(grid.null<=thr)))
  
  return(list(lfdr=lfdr.est, cdf=lfdr.cdf, cdf.func=grid.null, thr=thr, de=de, ev=estat))
}

ebh.func<-function(ev, q)
{ 
  # the input: 
  # ev: the e-values
  # q: the control FDR level
  # the output :
  # thr: the e-value threshold
  # de: the decision rule

  m=length(ev)
  st.ev<-sort(ev,decreasing = TRUE)   
  evi<-st.ev*1:m
  hps<-rep(0,m)
  k<-max(which(evi>=(m/q)))
  ek<-st.ev[k]
  hps[which(ev>=ek)]<-1

  return (list(thr=ek, de=hps))
}

uebh.func<-function(ev, q,s)
{ 
  # the input: 
  # ev: the e-values
  # q: the control FDR level
  # s: seed for replicability
  # the output :
  # thr: the e-value threshold
  # de: the decision rule
  
  set.seed(s)
  u = runif(1)
  m=length(ev)
  st.ev<-sort(ev,decreasing = TRUE)   
  evi<-st.ev*1:m
  hps<-rep(0,m)
  k<-max(which(evi>=(u*m/q)))
  ek<-st.ev[k]
  hps[which(ev>=ek)]<-1
  
  return (list(thr=ek, de=hps))
}

bh.func<-function(pv, q)
{ 
  # the input 
  # pv: the p-values
  # q: the FDR level
  # the output 
  # nr: the number of hypothesis to be rejected
  # th: the p-value threshold
  # re: the index of rejected hypotheses
  # ac: the index of accepted hypotheses
  # de: the decision rule
  
  m=length(pv)
  st.pv<-sort(pv)   
  #print(length(st.pv))
  #print(m)
  pvi<-st.pv/1:m
  hps<-rep(0, m)
  if (max(pvi<=(q/m))==0)
  {
    k<-0
    pk<-1
    reject<-NULL
    accept<-1:m
  }
  else
  {
    k<-max(which(pvi<=(q/m)))
    pk<-st.pv[k]
    reject<-which(pv<=pk)
    accept<-which(pv>pk)
    hps[reject]<-1
  }
  y<-list(nr=k, th=pk, re=reject, ac=accept, de=hps)
  return (y)
}

cauchy.comb<-function(pv,w){
  
  test_stat = sum(w*tan((0.5-pv)*pi))
  cauchy_pval = pcauchy(test_stat, location = 0, scale = 1, lower.tail = FALSE, log.p = FALSE)
  return(cauchy_pval)
}

AR1_chol<-function(p,rho=0.5)
{
  R = matrix(0,p,p);                # allocate p x p matrix
  R[1,] = rho^(0:(p-1));        # formula for 1st row
  cc = sqrt(1 - rho^2); # scaling factor: c^2 + rho^2 = 1
  R2 = cc*R[1,];        # formula for 2nd row
  for(j in 2 : p){    # shift elements in 2nd row for remaining rows
    R[j, j:p] = R2[1:(p-j+1)]
  }
  return(R)
}

GenPMat_corr <- function(M = 1000,Sigma_chol.1,Sigma_chol.2,
                         n = 2,
                         r = 2,
                         all.zero.frac = 0.8,
                         alternative.frac = 0.1,
                         mu = c(3,-3,4,-4,6,-6),
                         one.sided = F) 
{
  combs <- as.matrix(expand.grid(data.frame(matrix(rep(c(0, 1), 
                                                       n), nrow = 2))))
  category <- rowSums(combs)
  n.cases <- sum(category >= r)
  weights <- rep((1 - alternative.frac - all.zero.frac)/(2^n - n.cases - 1),
                 nrow(combs))
  weights[category == 0] <- all.zero.frac
  weights[category >= r] <- rep(alternative.frac/n.cases, n.cases)
  
  mean.matrix <- combs[sample(1:nrow(combs), M, replace = T, prob = weights), ]
  
  
  
  truth.pc <- rowSums(mean.matrix) >= r
  
  
  ## This is the signal matrix
  mean.matrix <- apply(mean.matrix, 1, function(v) {
    
    v[v == 1] <- sample(mu, sum(v == 1), replace = T)
    return(v)
  })
  mean.matrix <- t(mean.matrix)
  
  
  noise.matrix = matrix(rnorm(M*n),M,n)
  
  ## generate dependency across screening tests and agents
  U = Sigma_chol.1#cholesky(rho1*matrix(1,M,M)+(1-rho1)*diag(M))
  V = Sigma_chol.2#chol(rho2*matrix(1,n,n)+(1-rho2)*diag(n))
  
  zmat = mean.matrix+t(U)%*%noise.matrix%*%V
  
  signs <- 2 * (zmat >= 0) - 1
  raw.pvalues <- pnorm(abs(zmat), lower.tail = F)
  
  if (one.sided) {
    pvalue.matrix <- signs * raw.pvalues + (1 - signs)/2
  } else {
    pvalue.matrix <- 2 * raw.pvalues 
  }
  
  return(list(truth.pc = truth.pc, 
              pvalue.mat = pvalue.matrix,
              zmat = zmat))
}



##### From Nikos's Github

divide_Ps_by_Ws <- function(Ps, Ws) {
  pmin(ifelse(Ws > 0, Ps / Ws, 1), 1)
}



#' The tau-weighted BH multiple testing procedure
#'
#' @param Ps   Numeric vector of unadjusted p-values.
#' @param Ws   Numeric vector of multiple testing weights
#' @param tau  Numeric (default = 1), the level at which tau-censoring is applied. Forced to be <= tau_pi0 if Storey=TRUE.
#' @param tau_pi0  Numeric (default = 0.5), the threshold at which Storey is applied
#' @param Storey  Bool (default: FALSE): is the procedure pi0 adaptive or not?
#'
#' @return Vector of adjusted p-values
tau_weighted_bh <-
  function(Ps,
           Ws,
           tau = 1,
           tau_pi0 = 0.5,
           Storey = FALSE) {
    if (length(Ws) == 1) {
      Ws <- rep(Ws, length(Ps))
    }
    if (Storey) {
      pi0_hat <-  weighted_storey_pi0(Ps, Ws, tau = tau_pi0)
      Ws <- Ws / pi0_hat
      tau <- min(tau, tau_pi0)
    }
    weighted_pvals <- divide_Ps_by_Ws(Ps, Ws)
    weighted_pvals[Ps > tau] <- Inf
    adj_p <- stats::p.adjust(weighted_pvals, method = "BH")
    adj_p
  }

#' ep-BH
#'
#' @param Ps   Numeric vector of unadjusted p-values.
#' @param Es   Numeric vector of e-values
#' @param alpha   Significance level at which to apply method
#' @param Storey  Bool (default: FALSE): is the procedure pi0 adaptive or not?
#'
#' @return Binary vector of rejected/non-rejected hypotheses.
#' @export
ep_BH <- function(Ps, Es, alpha, Storey = FALSE) {
  if (Storey) {
    pi0_hat <- weighted_storey_pi0(Ps, 1, tau = 0.5)
    tau_cens <- 0.5
  } else {
    pi0_hat <- 1
    tau_cens <- 1
  }
  Es <- Es / pi0_hat
  adj_p <- tau_weighted_bh(Ps, Es, tau = tau_cens, Storey = FALSE)
  adj_p <= alpha
}




