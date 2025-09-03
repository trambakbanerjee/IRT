
library(Rfast)
library(poolr)
library(harmonicmeanp)
library(heavytailcombtest)

source('funcs.R')


#2 dependent agents (equi-correlated with negative to positive correlation)
nd = 2
m = 1000
reps = 500
alph = 0.02
null_prop = 0.8
null_lb = 0.5
rho1 = 0
rho2_vec = c(-0.9,-0.7,-0.5,-0.3,-0.1,0,0.1,0.3,0.5,0.7,0.9)

fisher_fdp = e_fdp = eprod_fdp = imputed_fdp = cauchy_fdp = trunct_fdp = hm_fdp = matrix(0,reps,length(rho2_vec))
fisher_etp = e_etp = eprod_etp = imputed_etp = cauchy_etp = trunct_etp = hm_etp = matrix(0,reps,length(rho2_vec))

for(rr in 1:length(rho2_vec)){
  
  rho2 = rho2_vec[rr]
  alphd = rep(0.01,nd)
  U = chol(rho1*matrix(1,m,m)+(1-rho1)*diag(m))
  V = chol(rho2*matrix(1,nd,nd)+(1-rho2)*diag(nd))#AR1_chol(nd,rho2)#
  
  for(r in 1:reps){
    set.seed(r)
    y = matrix(rnorm(m*nd),m,nd)
    p = runif(m)
    mu = sapply(1:m,function(i) 0*(p[i]<=null_prop)+rnorm(1,3,1)*(p[i]>null_prop))
    MU = t(sapply(1:m,function(i) rep(mu[i],nd)))
    theta = 1*(p>null_prop)
    x = MU+t(U)%*%y%*%V
    pvals = sapply(1:nd,function(i) pnorm(x[,i],0,1,lower.tail = FALSE))
    
    
    ### IRT
    bh_decisions =  sapply(1:nd,function(i) bh.func(pvals[,i],alphd[i])$de)
    estat = m*sapply(1:nd,function(i) bh_decisions[,i]/(alphd[i]*max(sum(bh_decisions[,i],1))))
    fedeval_mean = rowmeans(estat)
    ebh = ebh.func(fedeval_mean,alph)
    e_fdp[r,rr] = sum((1-theta)*ebh$de)/max(sum(ebh$de),1)
    e_etp[r,rr] = sum(theta*ebh$de)/max(sum(theta),1)
    
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
    
    ##########################################################
    ###### ---- calculating the product e-values ---- #######
    estat_prod = matrix(0,m,1)
    for(ii in 1:m){
      ee = null_lb*estat[ii,]
      if(max(ee)>0){
        tt = rep(0,nd)
        tt[1] = mean(ee)
        tt[nd] = prod(ee)
        for(j in 2:(nd-1)){
          
          combs = comb_n(nd,j)
          ncols = ncol(combs)
          tt[j] = mean(unlist(lapply(1:ncols,function(i) prod(ee[combs[,i]]))))
        }
        estat_prod[ii] = mean(tt[-1])
      }
    }
    
    ebh_prod = ebh.func(estat_prod,alph)#ebh.func(fedeval_mean,alph)
    eprod_fdp[r,rr] = sum((1-theta)*ebh_prod$de)/max(sum(ebh_prod$de),1)
    eprod_etp[r,rr] = sum(theta*ebh_prod$de)/max(sum(theta),1)
    ##########################################################
    
    ##########################################################
    #print(r)
  }
}

library(ggplot2)
library(ggpubr)
library(latex2exp)

plotdata1<- as.data.frame(c(colmeans(e_fdp),colmeans(eprod_fdp),colmeans(imputed_fdp),
                            colmeans(fisher_fdp),
                            colmeans(cauchy_fdp)#,colmeans(hm_fdp),
                            #colmeans(trunct_fdp)
                            ))
names(plotdata1)<-"FDR"
plotdata1$type<-as.factor(rep(c('IRT','IRT*','Cauchy + Imputed',
                                'Fisher','Cauchy'#,'HM','Trunc t'
                                ),each=length(rho2_vec)))
plotdata1$d<- rep(rho2_vec,5)
g1<-ggplot()+geom_line(data=plotdata1,aes(x=d,y=FDR,color=type),size=1)+
  geom_point(data=plotdata1,aes(x=d,y=FDR,color=type,shape=type),size=4,fill=NA)+
  geom_hline(yintercept = alph,linetype='dotted', size=1,col = 'black')+
  scale_y_continuous(limits = c(0,2*alph))+
  scale_x_continuous(breaks = rho2_vec)+
  ylab('FDP')+theme_bw()+xlab(TeX('$\\rho$'))+
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

plotdata2<- as.data.frame(c(colmeans(e_etp),colmeans(eprod_etp),colmeans(imputed_etp),
                            colmeans(fisher_etp),
                            colmeans(cauchy_etp)#,colmeans(hm_etp),
                            #colmeans(trunct_etp)
                            ))
names(plotdata2)<-"etp"
plotdata2$type<-as.factor(rep(c('IRT','IRT*','Cauchy + Imputed',
                                'Fisher','Cauchy'
                                #,'HM','Trunc t'
                                ),each=length(rho2_vec)))
plotdata2$d<- rep(rho2_vec,5)
g2<-ggplot()+geom_line(data=plotdata2,aes(x=d,y=etp,color=type),size=1)+
  geom_point(data=plotdata2,aes(x=d,y=etp,color=type,shape=type),size=4,fill=NA)+
  ylab('ETP')+theme_bw()+xlab(TeX('$\\rho$'))+
  scale_x_continuous(breaks = rho2_vec)+
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

save.image('scenario_3.Rdata')
