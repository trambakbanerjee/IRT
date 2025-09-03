
library(Rfast)
library(poolr)
library(harmonicmeanp)
library(heavytailcombtest)
library(Matrix)

setwd('C:/Users/Trambak Banerjee/Dropbox/Research/Federated multiple testing/simulations/biometrika major rev')
source('funcs.R')


nd = 15
d1 = 2
d2 = nd-d1
m = 1000
reps = 2000
alph = 0.005
null_prop = 0.95
rho1 = 0
rho2_vec = c(0,0.1,0.3,0.5,0.7,0.9)#c(-0.9,-0.7,-0.5,-0.3,-0.1,0,0.1,0.3,0.5,0.7,0.9)


onestudy_fdp = array(0,c(reps,length(rho2_vec),nd))
IRT_fdp =  cauchy_fdp = fisher_fdp_d1 = meanimputed_fdp = matrix(0,reps,length(rho2_vec))
onestudy_etp = array(0,c(reps,length(rho2_vec),nd))
IRT_etp = cauchy_etp = fisher_etp_d1 = meanimputed_etp = matrix(0,reps,length(rho2_vec))

alphd = rep(0.01,nd)


for(rr in 1:length(rho2_vec)){
  
  rho2=rho2_vec[rr]
  U = chol(rho1*matrix(1,m,m)+(1-rho1)*diag(m))
  S1 = diag(d1)
  S2 = rho2*matrix(1,d2,d2)+(1-rho2)*diag(d2)
  S = as.matrix(bdiag(S1,S2))
  V = chol(S)#AR1_chol(nd,rho2)#
  
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
    
    
    ### fisher (d1) + IRT (d2)
    #1. fisher combination p-values on the first d1 studies
    fisher_pvals_d1 = sapply(1:m,function(i) fisher(pvals[i,1:d1])$p)
    
    #2. IRT on the remaining studies
    estat = m*sapply((d1+1):nd,function(i) bh_decisions[,i]/(alphd[i]*max(sum(bh_decisions[,i],1))))
    fedeval_mean = rowmeans(estat)
    
    #3. Use ep-BH on fisher p-vals and IRT evals
  
    ebh = ep_BH(fisher_pvals_d1,fedeval_mean,alph)
    IRT_fdp[r,rr] = sum((1-theta)*ebh)/max(sum(ebh),1)
    IRT_etp[r,rr] = sum(theta*ebh)/max(sum(theta),1)
    
    ### All Cauchy pvals
    cauchy_pvals = sapply(1:m,function(i) cauchy.comb(pvals[i,],
                                                      rep(1/nd,nd)))
    bh_cauchy = bh.func(cauchy_pvals,alph)$de
    cauchy_fdp[r,rr] = sum((1-theta)*bh_cauchy)/max(sum(bh_cauchy),1)
    cauchy_etp[r,rr] = sum(theta*bh_cauchy)/max(sum(theta),1)
    
    ### Fisher pvals using only d1 studies
    bh_fisher_d1 = bh.func(fisher_pvals_d1,alph)$de
    fisher_fdp_d1[r,rr] = sum((1-theta)*bh_fisher_d1)/max(sum(bh_fisher_d1),1)
    fisher_etp_d1[r,rr] = sum(theta*bh_fisher_d1)/max(sum(theta),1)
    
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
    
    # ##### P2E ###############################
    # p2e = matrix(0,m,nd)
    # for(i in 1:m){
    #   id1 = (pvals[i,]==0)
    #   id2 = (pvals[i,]>0 & pvals[i,]<= exp(-2))
    #   p2e[i,id1] = Inf 
    #   p2e[i,id2] = 2/(pvals[i,id2]*((-log(pvals[i,id2]))^2))
    #   #p2e[i,] = 0.1*(pvals[i,]^{0.1-1})
    # }
    # p2e_mean = rowmeans(p2e)
    # ebh_p2e = ebh.func(p2e_mean,alph)
    # P2E_fdp[r,rr] = sum((1-theta)*ebh_p2e$de)/max(sum(ebh_p2e$de),1)
    # P2E_etp[r,rr] = sum(theta*ebh_p2e$de)/max(sum(theta),1)
    
    ##### Fisher ################################
    # pooled_pvals = sapply(1:m,function(i) fisher(pvals[i,])$p)
    # bh_pooled = bh.func(pooled_pvals,alph)$de
    # fisher_fdp[r,rr] = sum((1-theta)*bh_pooled)/max(sum(bh_pooled),1)
    # fisher_etp[r,rr] = sum(theta*bh_pooled)/max(sum(theta),1)
    # 
    # ### Cauchy combination p-values (Liu and Xie, 2020)
    # cauchy_pvals = sapply(1:m,function(i) cauchy.comb(pvals[i,],
    #                                                   rep(1/nd,nd)))
    # 
    # 
    # ### Truncated p-values (Gui et. al. 2023)
    # trunct_pvals = sapply(1:m,function(i) combination.test(pvals[i,], 
    #                                                        method = 't', tail.idx = 1, truncate = T, truncate.threshold = 0.9))
    # bh_trunct = bh.func(trunct_pvals,alph)$de
    # trunct_fdp[r,rr] = sum((1-theta)*bh_trunct)/max(sum(bh_trunct),1)
    # trunct_etp[r,rr] = sum(theta*bh_trunct)/max(sum(theta),1)
    # 
    # ### Harmonic mean p-values (Wilson 2019)
    # hm_pvals = sapply(1:m,function(i) p.hmp(pvals[i,],L=nd,
    #                                         multilevel = FALSE))
    # bh_hm = bh.func(hm_pvals,alph)$de
    # hm_fdp[r,rr] = sum((1-theta)*bh_hm)/max(sum(bh_hm),1)
    # hm_etp[r,rr] = sum(theta*bh_hm)/max(sum(theta),1)
    
    
    ##########################################################
    print(r)
  }
}

library(ggplot2)
library(ggpubr)
library(latex2exp)



plotdata1<- as.data.frame(c(rep(mean(sapply(1:nd,function(i) colmeans(onestudy_fdp[,,i]))),
                                length(rho2_vec)),colmeans(meanimputed_fdp),colmeans(IRT_fdp),colmeans(cauchy_fdp),colmeans(fisher_fdp_d1)))
names(plotdata1)<-"FDR"
plotdata1$type<-as.factor(rep(c('Single Study','Cauchy + Imputed','Fisher + IRT','Cauchy','Fisher d1'),each=length(rho2_vec)))
plotdata1$d<- rep(rho2_vec,5)
g1<-ggplot()+geom_line(data=plotdata1,aes(x=d,y=FDR,color=type),size=1)+
  geom_point(data=plotdata1,aes(x=d,y=FDR,color=type,shape=type),size=4,fill=NA)+
  geom_hline(yintercept = alph,linetype='dotted', size=1,col = 'black')+
  scale_y_continuous(limits = c(0,2*alph))+scale_x_continuous(limits = c(0, 0.9),
                          breaks = rho2_vec) +
  ylab('FDP')+theme_bw()+xlab(TeX('$\\rho_2$'))+
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
                            colmeans(fisher_etp_d1)))
names(plotdata2)<-"etp"
plotdata2$type<-as.factor(rep(c('Single Study','Cauchy + Imputed','Fisher + IRT','Cauchy','Fisher d1'),each=length(rho2_vec)))
plotdata2$d<- rep(rho2_vec,5)
g2<-ggplot()+geom_line(data=plotdata2,aes(x=d,y=etp,color=type),size=1)+
  geom_point(data=plotdata2,aes(x=d,y=etp,color=type,shape=type),size=4,fill=NA)+
  ylab('ETP')+theme_bw()+xlab(TeX('$\\rho_2$'))+scale_x_continuous(limits = c(0, 0.9),
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

save.image('Example_3.Rdata')
