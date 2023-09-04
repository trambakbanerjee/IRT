
source('funcs.R')

library(Rfast)
library(poolr)


d = c(5,10,15,20,25,30,35,40,45,50)
m = 1000
reps = 500
alph = 0.1
null_prop = 0.8
rho1 = 0.5
rho2 = 0.7

am_mean_fdp = fedeval_mean_fdp = p2e_mean_fdp = fed_naive_fdp = pooled_fdp = matrix(0,reps,length(d))
am_mean_etp = fedeval_mean_etp = p2e_mean_etp = fed_naive_etp = pooled_etp = matrix(0,reps,length(d))

for(dd in 1:length(d)){
  nd = d[dd]
  alphd = rep(0.01,nd)
  U = chol(rho1*matrix(1,m,m)+(1-rho1)*diag(m))
  V = chol(rho2*matrix(1,nd,nd)+(1-rho2)*diag(nd))
  
  for(r in 1:reps){
    set.seed(r)
    y = matrix(rnorm(m*nd),m,nd)
    p = runif(m)
    mu = sapply(1:m,function(i) 0*(p[i]<=null_prop)+rnorm(1,-3,1)*(p[i]>null_prop & p[i]<=0.9)+
                  rnorm(1,3,1)*(p[i]>0.9))
    MU = t(sapply(1:m,function(i) rep(mu[i],nd)))
    theta = 1*(p>null_prop)
    x = MU+t(U)%*%y%*%V
    pvals = sapply(1:nd,function(i) 2*pnorm(-abs(x[,i]),0,1))
    bh_decisions =  sapply(1:nd,function(i) bh.func(pvals[,i],alphd[i])$de)
    
    estat = m*sapply(1:nd,function(i) bh_decisions[,i]/(alphd[i]*max(sum(bh_decisions[,i],1))))
    
    fedeval_mean = rowmeans(estat)
    ebh = ebh.func(fedeval_mean,alph)
    fedeval_mean_fdp[r,dd] = sum((1-theta)*ebh$de)/max(sum(ebh$de),1)
    fedeval_mean_etp[r,dd] = sum(theta*ebh$de)/max(sum(theta),1)
    
    count_decisions = rowsums(bh_decisions)
    naive_decisions = sapply(1:m,function(i) 1*(count_decisions[i]>=(nd/2)))
    fed_naive_fdp[r,dd] = sum((1-theta)*naive_decisions)/max(sum(naive_decisions),1)
    fed_naive_etp[r,dd] = sum(theta*naive_decisions)/max(sum(theta),1)
    
    pooled_pvals = sapply(1:m,function(i) fisher(pvals[i,])$p)
    bh_pooled = bh.func(pooled_pvals,alph)$de
    pooled_fdp[r,dd] = sum((1-theta)*bh_pooled)/max(sum(bh_pooled),1)
    pooled_etp[r,dd] = sum(theta*bh_pooled)/max(sum(theta),1)
    
    ##### p2e ###############################
    p2e = matrix(0,m,nd)
    for(i in 1:m){
      id1 = (pvals[i,]==0)
      id2 = (pvals[i,]>0 & pvals[i,]<= exp(-2))
      p2e[i,id1] = Inf 
      p2e[i,id2] = 2/(pvals[i,id2]*((-log(pvals[i,id2]))^2))
    }
    p2e_mean = rowmeans(p2e)
    ebh_p2e = ebh.func(p2e_mean,alph)
    p2e_mean_fdp[r,dd] = sum((1-theta)*ebh_p2e$de)/max(sum(ebh_p2e$de),1)
    p2e_mean_etp[r,dd] = sum(theta*ebh_p2e$de)/max(sum(theta),1)
    ###########################################
    ######## Harmonic Mean p-value
    am_pvals = sapply(1:m,function(i) nd/sum(1/pvals[i,]))#log(nd)*sapply(1:m,function(i) nd/sum(1/pvals[i,]))
    ebh_am = ebh.func(am_pvals,alph)
    am_mean_fdp[r,dd] = sum((1-theta)*ebh_am$de)/max(sum(ebh_am$de),1)
    am_mean_etp[r,dd] = sum(theta*ebh_am$de)/max(sum(theta),1)
    
    ###################################
    print(r)
  }
}

library(ggplot2)
library(ggpubr)

plotdata1<- as.data.frame(c(colmeans(fedeval_mean_fdp),colmeans(p2e_mean_fdp),
                            colmeans(pooled_fdp),
                            colmeans(fed_naive_fdp),
                            colmeans(am_mean_fdp)))
names(plotdata1)<-"FDR"
plotdata1$type<-as.factor(rep(c('IRT','P2E','Fisher','Naive','HM'),each=length(d)))
plotdata1$d<- rep(d,5)
g1<-ggplot()+geom_line(data=plotdata1,aes(x=d,y=FDR,color=type),size=1)+
  geom_point(data=plotdata1,aes(x=d,y=FDR,color=type,shape=type),size=4,fill=NA)+
  geom_hline(yintercept = 0.1,linetype='dotted', size=1,col = 'black')+
  scale_y_continuous(limits = c(0,0.2))+
  ylab('FDP')+theme_bw()+xlab('d')+
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

plotdata2<- as.data.frame(c(colmeans(fedeval_mean_etp),colmeans(p2e_mean_etp),
                            colmeans(pooled_etp),
                            colmeans(fed_naive_etp),
                            colmeans(am_mean_etp)))
names(plotdata2)<-"etp"
plotdata2$type<-as.factor(rep(c('IRT','P2E','Fisher','Naive','HM'),each=length(d)))
plotdata2$d<- rep(d,5)
g2<-ggplot()+geom_line(data=plotdata2,aes(x=d,y=etp,color=type),size=1)+
  geom_point(data=plotdata2,aes(x=d,y=etp,color=type,shape=type),size=4,fill=NA)+
  ylab('ETP')+theme_bw()+xlab('d')+
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

save.image('scenario 3.Rdata')
