
library(Rfast)
library(poolr)
library(harmonicmeanp)
library(heavytailcombtest)
library(Matrix)

source('funcs.R')

# d1 studies are independent give p-values
# d2 studies are independent but only give 0/1 decisions
# d1 and d2 studies are independent


nd = 15
d2 = 13
m = 1000
reps = 1000
alph = 0.005
null_prop_vec = c(0.5,0.8,0.95,0.98)
rho1 = 0
null_lb = 0.5

onestudy_fdp = array(0,c(reps,length(null_prop_vec),nd))
IRTst_fdp = IRT_fdp =  fisher_fdp = fisher_fdp_d1 = meanimpute_fdp = matrix(0,reps,length(null_prop_vec))
onestudy_etp = array(0,c(reps,length(null_prop_vec),nd))
IRTst_etp =IRT_etp = fisher_etp = fisher_etp_d1 = meanimpute_etp = matrix(0,reps,length(null_prop_vec))

alphd = rep(0.01,nd)


for(rr in 1:length(null_prop_vec)){
  
  null_prop = null_prop_vec[rr]
  ei = null_prop+(1-null_prop)/2
  d1 = nd-d2
  U = chol(rho1*matrix(1,m,m)+(1-rho1)*diag(m))
  V = chol(diag(d1+d2))#AR1_chol(nd,rho2)#
  
  for(r in 1:reps){
    set.seed(r)
    y = matrix(rnorm(m*nd),m,nd)
    p = runif(m)
    mu = sapply(1:m,function(i) 0*(p[i]<=null_prop)+rnorm(1,3,1)*(p[i]>null_prop & p[i]<=ei)+
                  rnorm(1,-3,1)*(p[i]>ei))
    MU = t(sapply(1:m,function(i) rep(mu[i],nd)))
    theta = 1*(abs(mu)>0)
    x = MU+t(U)%*%y%*%V
    pvals = sapply(1:nd,function(i) 2*pnorm(-abs(x[,i]),0,1))
    
    bh_decisions =  sapply(1:nd,function(i) bh.func(pvals[,i],alphd[i])$de)
    
    ### Best single study
    bh_onestudy = sapply(1:nd,function(i) bh.func(pvals[,i],alphd[i])$de)
    onestudy_fdp[r,rr,] = sapply(1:nd,function(i) sum((1-theta)*bh_onestudy[,i])/max(sum(bh_onestudy[,i]),1))
    onestudy_etp[r,rr,] = sapply(1:nd,function(i) sum(theta*bh_onestudy[,i])/max(sum(theta),1))
    
    
    ### fisher (d1) + IRT* (d2)
    #1. fisher combination p-values on the first d1 studies
    fisher_pvals_d1 = sapply(1:m,function(i) fisher(pvals[i,1:d1])$p)
    
    #2. IRT* on the remaining studies
    ###### ---- calculating the product e-values ---- #######
    estat = m*sapply((d1+1):nd,function(i) bh_decisions[,i]/(alphd[i]*max(sum(bh_decisions[,i],1))))
    estat_prod = matrix(0,m,1)
    for(ii in 1:m){
      ee = estat[ii,]
      if(max(ee)>0){
        tt = rep(0,d2)
        tt[1] = null_lb*mean(ee)
        tt[d2] = (null_lb^d2)*prod(ee)
        for(j in 2:(d2-1)){
          
          combs = comb_n(d2,j)
          ncols = ncol(combs)
          tt[j] = (null_lb^j)*mean(unlist(lapply(1:ncols,function(i) prod(ee[combs[,i]]))))
          
        }
        estat_prod[ii] = mean(tt[-1])
      }
    }
    ##########################################################
    #3. Use ep-BH on fisher p-vals and IRT* evals
  
    ebh = ep_BH(fisher_pvals_d1,estat_prod,alph)
    IRTst_fdp[r,rr] = sum((1-theta)*ebh)/max(sum(ebh),1)
    IRTst_etp[r,rr] = sum(theta*ebh)/max(sum(theta),1)
    
    ### All Fisher pvals
    fisher_pvals = sapply(1:m,function(i) fisher(pvals[i,])$p)
    bh_fisher = bh.func(fisher_pvals,alph)$de
    fisher_fdp[r,rr] = sum((1-theta)*bh_fisher)/max(sum(bh_fisher),1)
    fisher_etp[r,rr] = sum(theta*bh_fisher)/max(sum(theta),1)
    
    ### Fisher pvals using only d1 studies
    bh_fisher_d1 = bh.func(fisher_pvals_d1,alph)$de
    fisher_fdp_d1[r,rr] = sum((1-theta)*bh_fisher_d1)/max(sum(bh_fisher_d1),1)
    fisher_etp_d1[r,rr] = sum(theta*bh_fisher_d1)/max(sum(theta),1)
    
    ### Mean imputation pvals for d2 studies
    imputed = matrix(0,m,d2)
    for(j in 1:d2){
      
      tj = (alphd[d1+j]/m)*sum(bh_decisions[,d1+j])
      imputed[,j] = sapply(1:m,function(i) (tj)*(bh_decisions[i,d1+j]==1)+
                             1*(bh_decisions[i,d1+j]==0))
    }
    fisher_pvals = sapply(1:m,function(i) fisher(c(pvals[i,1:d1],imputed[i,]))$p)
    bh_meanimpute = bh.func(fisher_pvals,alph)$de
    meanimpute_fdp[r,rr] = sum((1-theta)*bh_meanimpute)/max(sum(bh_meanimpute),1)
    meanimpute_etp[r,rr] = sum(theta*bh_meanimpute)/max(sum(theta),1)
    
    
    ### fisher (d1) + IRT (d2)
    #2. IRT on the d2 studies
    fedeval_mean = rowmeans(estat)
    ebh = ep_BH(fisher_pvals_d1,fedeval_mean,alph)
    IRT_fdp[r,rr] = sum((1-theta)*ebh)/max(sum(ebh),1)
    IRT_etp[r,rr] = sum(theta*ebh)/max(sum(theta),1)
   
    
    ##########################################################
    print(r)
  }
}

library(ggplot2)
library(ggpubr)
library(latex2exp)

plotdata1<- as.data.frame(c(rep(mean(sapply(1:nd,function(i) colmeans(onestudy_fdp[,,i]))),
                                length(null_prop_vec)),colmeans(meanimpute_fdp),colmeans(IRTst_fdp),
                            colmeans(IRT_fdp),colmeans(fisher_fdp),colmeans(fisher_fdp_d1)))
names(plotdata1)<-"FDR"
plotdata1$type<-as.factor(rep(c('Single Study','Fisher + Imputed','Fisher + IRT*','Fisher + IRT','Fisher','Fisher d1'),
                              each=length(null_prop_vec)))
plotdata1$d<- rep(null_prop_vec,6)
g1<-ggplot()+geom_line(data=plotdata1,aes(x=d,y=FDR,color=type),size=1)+
  geom_point(data=plotdata1,aes(x=d,y=FDR,color=type,shape=type),size=4,fill=NA)+
  geom_hline(yintercept = alph,linetype='dotted', size=1,col = 'black')+
  scale_y_continuous(limits = c(0,2*alph))+
  ylab('FDP')+theme_bw()+xlab(TeX('$d_2$'))+
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
                                length(null_prop_vec)),colmeans(meanimpute_etp),colmeans(IRTst_etp),
                            colmeans(IRT_etp),colmeans(fisher_etp),
                            colmeans(fisher_etp_d1)))
names(plotdata2)<-"etp"
plotdata2$type<-as.factor(rep(c('Single Study','Fisher + Imputed','Fisher + IRT*','Fisher + IRT','Fisher','Fisher d1'),
                              each=length(null_prop_vec)))
plotdata2$d<- rep(null_prop_vec,6)
g2<-ggplot()+geom_line(data=plotdata2,aes(x=d,y=etp,color=type),size=1)+
  geom_point(data=plotdata2,aes(x=d,y=etp,color=type,shape=type),size=4,fill=NA)+
  ylab('ETP')+theme_bw()+xlab(TeX('$d_2$'))+
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

save.image('Example_4.Rdata')
