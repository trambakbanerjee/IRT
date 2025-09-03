
library(Rfast)
library(poolr)
library(harmonicmeanp)
library(heavytailcombtest)

source('funcs.R')


# 5 dependent agents (equi-correlated with positive correlation) and
# conservative p-values
nd = 5
m = 1000
reps = 500
alph = 0.02
null_prop = 0.8
rho1 = 0
rho2 = 0.5
df = c(3,5,10,15,20)

fisher_fdp = imputed_fdp = cauchy_fdp = IRT_fdp = trunct_fdp = hm_fdp = matrix(0,reps,length(df))
fisher_etp = imputed_etp = cauchy_etp = IRT_etp = trunct_etp = hm_etp = matrix(0,reps,length(df))

alphd = rep(0.01,nd)
U = chol(rho1*matrix(1,m,m)+(1-rho1)*diag(m))
V = chol(rho2*matrix(1,nd,nd)+(1-rho2)*diag(nd))

for(rr in 1:length(df)){
  
  for(r in 1:reps){
    set.seed(r)
    y = matrix(rnorm(m*nd),m,nd)
    p = runif(m)
    mu = sapply(1:m,function(i) 0*(p[i]<=null_prop)+rnorm(1,3,1)*(p[i]>null_prop))
    MU = t(sapply(1:m,function(i) rep(mu[i],nd)))
    theta = 1*(p>null_prop)
    x = MU+t(U)%*%y%*%V
    pvals = sapply(1:nd,function(i) pt(x[,i],df[rr],ncp=0,lower.tail = FALSE))
    
    
    ### IRT
    bh_decisions =  sapply(1:nd,function(i) bh.func(pvals[,i],alphd[i])$de)
    
    estat = m*sapply(1:nd,function(i) bh_decisions[,i]/(alphd[i]*max(sum(bh_decisions[,i],1))))
    
    fedeval_mean = rowmeans(estat)
    ebh = ebh.func(fedeval_mean,alph)
    IRT_fdp[r,rr] = sum((1-theta)*ebh$de)/max(sum(ebh$de),1)
    IRT_etp[r,rr] = sum(theta*ebh$de)/max(sum(theta),1)
    
     ##### Fisher ################################
    pooled_pvals = sapply(1:m,function(i) fisher(pvals[i,])$p)
    bh_pooled = bh.func(pooled_pvals,alph)$de
    fisher_fdp[r,rr] = sum((1-theta)*bh_pooled)/max(sum(bh_pooled),1)
    fisher_etp[r,rr] = sum(theta*bh_pooled)/max(sum(theta),1)
    
    ### Cauchy combination p-values (Liu and Xie, 2020)
    cauchy_pvals = sapply(1:m,function(i) cauchy.comb(pvals[i,],
                                                      rep(1/nd,nd)))
    bh_cauchy = bh.func(cauchy_pvals,alph)$de
    cauchy_fdp[r,rr] = sum((1-theta)*bh_cauchy)/max(sum(bh_cauchy),1)
    cauchy_etp[r,rr] = sum(theta*bh_cauchy)/max(sum(theta),1)
    
    imputed = matrix(0,m,nd)
    for(j in 1:nd){
      
      tj = (alphd[j]/m)*sum(bh_decisions[,j])
      imputed[,j] = sapply(1:m,function(i) (tj)*(bh_decisions[i,j]==1)+
                             1*(bh_decisions[i,j]==0))
    }
    cauchy_pvals = sapply(1:m,function(i) cauchy.comb(imputed[i,],
                                                      rep(1/nd,nd)))
    bh_meanimpute = bh.func(cauchy_pvals,alph)$de
    imputed_fdp[r,rr] = sum((1-theta)*bh_meanimpute)/max(sum(bh_meanimpute),1)
    imputed_etp[r,rr] = sum(theta*bh_meanimpute)/max(sum(theta),1)
    
    ### Truncated p-values (Gui et. al. 2023)
    trunct_pvals = sapply(1:m,function(i) combination.test(pvals[i,], 
                                                           method = 't', tail.idx = 1, truncate = T, truncate.threshold = 0.9))
    bh_trunct = bh.func(trunct_pvals,alph)$de
    trunct_fdp[r,rr] = sum((1-theta)*bh_trunct)/max(sum(bh_trunct),1)
    trunct_etp[r,rr] = sum(theta*bh_trunct)/max(sum(theta),1)
    
    ### Harmonic mean p-values (Wilson 2019)
    hm_pvals = sapply(1:m,function(i) p.hmp(pvals[i,],L=nd,
                                            multilevel = FALSE))
    bh_hm = bh.func(hm_pvals,alph)$de
    hm_fdp[r,rr] = sum((1-theta)*bh_hm)/max(sum(bh_hm),1)
    hm_etp[r,rr] = sum(theta*bh_hm)/max(sum(theta),1)
    
    ##########################################################
    print(r)
  }
}

library(ggplot2)
library(ggpubr)
library(latex2exp)

plotdata1<- as.data.frame(c(colmeans(IRT_fdp),colmeans(fisher_fdp),
                            colmeans(imputed_fdp),colmeans(cauchy_fdp)))
names(plotdata1)<-"FDR"
plotdata1$type<-as.factor(rep(c('IRT','Fisher','Cauchy + Imputed','Cauchy'),each=length(df)))
plotdata1$d<- rep(df,4)
g1<-ggplot()+geom_line(data=plotdata1,aes(x=d,y=FDR,color=type),size=1)+
  geom_point(data=plotdata1,aes(x=d,y=FDR,color=type,shape=type),size=4,fill=NA)+
  geom_hline(yintercept = alph,linetype='dotted', size=1,col = 'black')+
  scale_y_continuous(limits = c(0,2*alph))+
  ylab('FDP')+theme_bw()+xlab(TeX('$\\nu$'))+
  theme(legend.position='top',legend.title=element_blank(),
        legend.background = element_rect(fill="white",
                                         size=1, linetype="solid", colour ="black"),
        legend.text=element_text(size=12),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

plotdata2<- as.data.frame(c(colmeans(IRT_etp),colmeans(fisher_etp),
                            colmeans(imputed_etp),colmeans(cauchy_etp)))
names(plotdata2)<-"etp"
plotdata2$type<-as.factor(rep(c('IRT','Fisher','Cauchy + Imputed','Cauchy'),each=length(df)))
plotdata2$d<- rep(df,4)
g2<-ggplot()+geom_line(data=plotdata2,aes(x=d,y=etp,color=type),size=1)+
  geom_point(data=plotdata2,aes(x=d,y=etp,color=type,shape=type),size=4,fill=NA)+
  ylab('ETP')+theme_bw()+xlab(TeX('$\\nu$'))+
  theme(legend.position='top',legend.title=element_blank(),
        legend.background = element_rect(fill="white",
                                         size=1, linetype="solid", colour ="black"),
        legend.text=element_text(size=12),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggarrange(g1,g2,ncol=2,nrow=1,common.legend = TRUE)

save.image('scenario_7_b.Rdata')
