
library(Rfast)
library(poolr)
library(harmonicmeanp)
library(heavytailcombtest)
library(Matrix)

source('funcs.R')

# d1 studies are equi-correlated but give p-values
# d2 studies are independent but only give 0/1 decisions
# d1 studies independent of d2 studies.

nd = 5
d1 = 2
d2 = nd-d1
m = 1000
reps = 2000
alph = 0.01#0.005
null_prop = 0.95
rho1 = 0
rho2_vec = c(0,0.1,0.3,0.5,0.7,0.9)


onestudy_fdp = array(0,c(reps,length(rho2_vec),nd))
IRT_fdp =  cauchy_fdp = cauchy_fdp_d1 = meanimputed_fdp = matrix(0,reps,length(rho2_vec))
onestudy_etp = array(0,c(reps,length(rho2_vec),nd))
IRT_etp = cauchy_etp = cauchy_etp_d1 = meanimputed_etp = matrix(0,reps,length(rho2_vec))

alphd = rep(0.01,nd)


for(rr in 1:length(rho2_vec)){
  
  rho2=rho2_vec[rr]
  U = chol(rho1*matrix(1,m,m)+(1-rho1)*diag(m))
  S1 = rho2*matrix(1,d1,d1)+(1-rho2)*diag(d1)
  S2 = diag(d2)
  S = as.matrix(bdiag(S1,S2))
  V = chol(S)
  
  for(r in 1:reps){
    set.seed(r)
    y = matrix(rnorm(m*nd),m,nd)
    p = runif(m)
    mu = sapply(1:m,function(i) 0*(p[i]<=null_prop)+rnorm(1,3,1)*(p[i]>null_prop & p[i]<=0.975)+
                  rnorm(1,-3,1)*(p[i]>0.975))
    MU = t(sapply(1:m,function(i) rep(mu[i],nd)))
    theta = 1*(abs(mu)>0)
    x = MU+t(U)%*%y%*%V
    pvals = sapply(1:nd,function(i) 2*pnorm(-abs(x[,i]),0,1))
    
    bh_decisions =  sapply(1:nd,function(i) bh.func(pvals[,i],alphd[i])$de)
    
    ### Best single study
    bh_onestudy = sapply(1:nd,function(i) bh.func(pvals[,i],alphd[i])$de)
    onestudy_fdp[r,rr,] = sapply(1:nd,function(i) sum((1-theta)*bh_onestudy[,i])/max(sum(bh_onestudy[,i]),1))
    onestudy_etp[r,rr,] = sapply(1:nd,function(i) sum(theta*bh_onestudy[,i])/max(sum(theta),1))
    
    
    ### Cauchy (d1) + IRT (d2)
    #1. Cauchy combination p-values (Liu and Xie, 2020) on the first d1 studies
    cauchy_pvals_d1 = sapply(1:m,function(i) cauchy.comb(pvals[i,1:d1],
                                                         rep(1/d1,d1)))
    
    #2. IRT on the remaining studies
    estat = m*sapply((d1+1):nd,function(i) bh_decisions[,i]/(alphd[i]*max(sum(bh_decisions[,i],1))))
    fedeval_mean = rowmeans(estat)
    
    #3. Use ep-BH on Cauchy p-vals and IRT evals
  
    ebh = ep_BH(cauchy_pvals_d1,fedeval_mean,alph)
    IRT_fdp[r,rr] = sum((1-theta)*ebh)/max(sum(ebh),1)
    IRT_etp[r,rr] = sum(theta*ebh)/max(sum(theta),1)
    
    ### All Cauchy pvals
    cauchy_pvals = sapply(1:m,function(i) cauchy.comb(pvals[i,],
                                                      rep(1/nd,nd)))
    bh_cauchy = bh.func(cauchy_pvals,alph)$de
    cauchy_fdp[r,rr] = sum((1-theta)*bh_cauchy)/max(sum(bh_cauchy),1)
    cauchy_etp[r,rr] = sum(theta*bh_cauchy)/max(sum(theta),1)
    
    ### Cauchy pvals using only d1 studies
    bh_cauchy_d1 = bh.func(cauchy_pvals_d1,alph)$de
    cauchy_fdp_d1[r,rr] = sum((1-theta)*bh_cauchy_d1)/max(sum(bh_cauchy_d1),1)
    cauchy_etp_d1[r,rr] = sum(theta*bh_cauchy_d1)/max(sum(theta),1)
    
    ### Cauchy pvals using only d1 studies+Mean imputation pvals for d2 studies
    imputed = matrix(0,m,d2)
    for(j in 1:d2){
      
      tj = (alphd[d1+j]/m)*sum(bh_decisions[,d1+j])
      imputed[,j] = sapply(1:m,function(i) (tj)*(bh_decisions[i,d1+j]==1)+
                             1*(bh_decisions[i,d1+j]==0))
    }
    cauchy_pvals = sapply(1:m,function(i) cauchy.comb(c(pvals[i,1:d1],imputed[i,]),
                                                      rep(1/nd,nd)))
    bh_meanimpute = bh.func(cauchy_pvals,alph)$de
    meanimputed_fdp[r,rr] = sum((1-theta)*bh_meanimpute)/max(sum(bh_meanimpute),1)
    meanimputed_etp[r,rr] = sum(theta*bh_meanimpute)/max(sum(theta),1)
    
    ##########################################################
    print(r)
  }
}

library(ggplot2)
library(ggpubr)
library(latex2exp)

plotdata1<- as.data.frame(c(rep(mean(sapply(1:nd,function(i) colmeans(onestudy_fdp[,,i]))),
                                length(rho2_vec)),colmeans(meanimputed_fdp),colmeans(IRT_fdp),colmeans(cauchy_fdp),colmeans(cauchy_fdp_d1)))
names(plotdata1)<-"FDR"
plotdata1$type<-as.factor(rep(c('Single Study','Cauchy + Imputed','Cauchy + IRT','Cauchy','Cauchy d1'),each=length(rho2_vec)))
plotdata1$d<- rep(rho2_vec,5)
g1<-ggplot()+geom_line(data=plotdata1,aes(x=d,y=FDR,color=type),size=1)+
  geom_point(data=plotdata1,aes(x=d,y=FDR,color=type,shape=type),size=4,fill=NA)+
  geom_hline(yintercept = alph,linetype='dotted', size=1,col = 'black')+
  scale_y_continuous(limits = c(0,2*alph))+
  ylab('FDP')+theme_bw()+xlab(TeX('$\\rho_1$'))+scale_x_continuous(limits = c(0, 0.9),
                                                                 breaks = rho2_vec) +
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

plotdata2<- as.data.frame(c(rep(max(sapply(1:nd,function(i) colmeans(onestudy_etp[,,i]))),
                                length(rho2_vec)),colmeans(meanimputed_etp),
                            colmeans(IRT_etp),colmeans(cauchy_etp),
                            colmeans(cauchy_etp_d1)))
names(plotdata2)<-"etp"
plotdata2$type<-as.factor(rep(c('Single Study','Cauchy + Imputed','Cauchy + IRT','Cauchy','Cauchy d1'),each=length(rho2_vec)))
plotdata2$d<- rep(rho2_vec,5)
g2<-ggplot()+geom_line(data=plotdata2,aes(x=d,y=etp,color=type),size=1)+
  geom_point(data=plotdata2,aes(x=d,y=etp,color=type,shape=type),size=4,fill=NA)+
  ylab('ETP')+theme_bw()+xlab(TeX('$\\rho_1$'))+scale_x_continuous(limits = c(0, 0.9),
                                                                 breaks = rho2_vec) +
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

save.image('Example_2.Rdata')
