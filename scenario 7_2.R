

library(Rfast)
library(poolr)
library(locfdr)

source('funcs.R')

d = 10
df = c(5,10,15,30,50,100)
m = 1000
reps = 200
alph = 0.02
null_prop = 0.8
rho = 0
null_lb = 0.5
rho1 = 0
rho2 = 0.7

e_fdp = p2e_fdp = matrix(0,reps,length(df))
e_etp = p2e_etp = matrix(0,reps,length(df))

nd = d
alphd = rep(0.01,nd)
U = chol(rho1*matrix(1,m,m)+(1-rho1)*diag(m))
V = chol(rho2*matrix(1,nd,nd)+(1-rho2)*diag(nd))
for(dd in 1:length(df)){
  
  for(r in 1:reps){
    set.seed(r)
    y = matrix(rt(m*nd,df[dd]),m,nd)
    p = runif(m)
    mu = sapply(1:m,function(i) 0*(p[i]<=null_prop)+rnorm(1,-3,1)*(p[i]>null_prop & p[i]<=0.9)+
                  rnorm(1,3,1)*(p[i]>0.9))
    MU = t(sapply(1:m,function(i) rep(mu[i],nd)))
    theta = 1*(p>null_prop)
    x = MU+t(U)%*%y%*%V
    sig = sqrt(rep(df[dd]/(df[dd]-2),nd))
    pvals = sapply(1:nd,function(i) 2*pnorm(-abs(x[,i])/sig[i],0,1))
    bh_decisions =  sapply(1:nd,function(i) bh.func(pvals[,i],alphd[i])$de)

    estat = m*sapply(1:nd,function(i) bh_decisions[,i]/(alphd[i]*max(sum(bh_decisions[,i],1))))
    fedeval_mean = rowmeans(estat)
    ebh = ebh.func(fedeval_mean,alph)#uebh.func(fedeval_mean,alph,r)#ebh.func(fedeval_mean,alph)
    e_fdp[r,dd] = sum((1-theta)*ebh$de)/max(sum(ebh$de),1)
    e_etp[r,dd] = sum(theta*ebh$de)/max(sum(theta),1)
    
    
    ##### p2e ###############################
    p2e = matrix(0,m,nd)
    for(i in 1:m){
      id1 = (pvals[i,]==0)
      id2 = (pvals[i,]>0 & pvals[i,]<= exp(-2))
      p2e[i,id1] = Inf 
      p2e[i,id2] = 2/(pvals[i,id2]*((-log(pvals[i,id2]))^2))
    }
    p2e_mean = rowmeans(p2e)
    ebh_p2e = ebh.func(p2e_mean,alph)#uebh.func(p2e_mean,alph,r)
    p2e_fdp[r,dd] = sum((1-theta)*ebh_p2e$de)/max(sum(ebh_p2e$de),1)
    p2e_etp[r,dd] = sum(theta*ebh_p2e$de)/max(sum(theta),1)
    ###########################################
     print(r)
  }
}

library(ggplot2)
library(ggpubr)

plotdata1<- as.data.frame(c(colmeans(e_fdp),
                            colmeans(p2e_fdp)))
names(plotdata1)<-"FDR"
plotdata1$type<-as.factor(rep(c('IRT','P2E'),each=length(df)))
plotdata1$d<- rep(df,2)
g1<-ggplot()+geom_line(data=plotdata1,aes(x=d,y=FDR,color=type),size=1)+
  geom_point(data=plotdata1,aes(x=d,y=FDR,color=type,shape=type),size=4,fill=NA)+
  geom_hline(yintercept = alph,linetype='dotted', size=1,col = 'black')+
  scale_x_continuous(breaks = df)+
  ylab('FDP')+theme_bw()+xlab('Degrees of freedom')+
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

plotdata2<- as.data.frame(c(colmeans(e_etp),
                            colmeans(p2e_etp)))
names(plotdata2)<-"etp"
plotdata2$type<-as.factor(rep(c('IRT','P2E'),each=length(df)))
plotdata2$d<- rep(df,2)
g2<-ggplot()+geom_line(data=plotdata2,aes(x=d,y=etp,color=type),size=1)+
  geom_point(data=plotdata2,aes(x=d,y=etp,color=type,shape=type),size=4,fill=NA)+
  scale_x_continuous(breaks = df)+
  ylab('ETP')+theme_bw()+xlab('Degrees of freedom')+
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

save.image('scenario 7_2.Rdata')
