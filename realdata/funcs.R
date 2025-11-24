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




