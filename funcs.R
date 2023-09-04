
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




