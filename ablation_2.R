
library(Rfast)
library(poolr)
library(Matrix)

source('funcs.R')


#Within each study, the p-values are not exchangeable
d = 10
m = 1000
reps = 2000
alph = 0.005
null_prop = 0.8
null_lb = 0.5
rho2_vec = c(0.1,0.3,0.5,0.7,0.9)
rho1_vec = c(0.3,0.5,0.7)

rhogrid = expand.grid(rho2_vec,rho1_vec)

ehybrid_fdp = eprod_fdp = matrix(0,reps,dim(rhogrid)[1])
ehybrid_etp = eprod_etp = matrix(0,reps,dim(rhogrid)[1])

for(dd in 1:dim(rhogrid)[1]){
  
  m1 = ceiling(m/3)
  m2 = m-(2*m1)
  nd = d
  d1 = floor(nd/2)
  d2 = nd-d1
  alphd = rep(0.01,nd)
  rho2 = rhogrid[dd,1]
  rho1 = rhogrid[dd,2]
  S1 = diag(m1)
  S2 = rho1*matrix(1,m1,m1)+(1-rho1)*diag(m1)
  S3 = 0.9*matrix(1,m2,m2)+(1-0.9)*diag(m2)
  S = as.matrix(bdiag(S1,S2,S3))
  U = chol(S)
#   U = AR1_chol(m,rho1)
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
    eprod_fdp[r,dd] = sum((1-theta)*ebh_prod$de)/max(sum(ebh_prod$de),1)
    eprod_etp[r,dd] = sum(theta*ebh_prod$de)/max(sum(theta),1)
    ##########################################################

    ###### ---- calculating the IRT H ---- #######
    estat_1 = m*sapply(1:d1,function(i) bh_decisions[,i]/(alphd[i]*max(sum(bh_decisions[,i],1))))
    fedeval_mean_1 = rowmeans(estat_1)
    estat_prod_2 = matrix(0,m,1)
    for(ii in 1:m){
      ee = null_lb*estat[ii,(d1+1):nd]
      if(max(ee)>0){
        if(d2 == 1){
          estat_prod_2[ii] = ee
        } else if(d2 == 2){
          estat_prod_2[ii] = mean(mean(ee),prod(ee))
        } else{
          tt = rep(0,d2)
          tt[1] = mean(ee)
          tt[d2] = prod(ee)
          for(j in 2:(d2-1)){
            
            combs = comb_n(d2,j)
            ncols = ncol(combs)
            tt[j] = mean(unlist(lapply(1:ncols,function(i) prod(ee[combs[,i]]))))
            
          }
          estat_prod_2[ii] = mean(tt[-1])
        }
      }
    }
    ehybrid = sapply(1:m,function(i) (d1/nd)*fedeval_mean_1[i]+(d2/nd)*estat_prod_2[i])
    ebh = ebh.func(ehybrid,alph)
    ehybrid_fdp[r,dd] = sum((1-theta)*ebh$de)/max(sum(ebh$de),1)
    ehybrid_etp[r,dd] = sum(theta*ebh$de)/max(sum(theta),1)
    ##########################################################
    print(r)
  }
}

library(ggplot2)
library(ggpubr)
library(latex2exp)

plotdata1<- as.data.frame(c(colmeans(eprod_fdp[,1:5]),
                            colmeans(ehybrid_fdp[,1:5])))
names(plotdata1)<-"FDR"
plotdata1$type<-as.factor(rep(c('IRT*','IRT H'),each=length(rho2_vec)))
plotdata1$d<- rep(rho2_vec,2)
g1<-ggplot()+geom_line(data=plotdata1,aes(x=d,y=FDR,color=type),size=1)+
  geom_point(data=plotdata1,aes(x=d,y=FDR,color=type,shape=type),size=4,fill=NA)+
  geom_hline(yintercept = alph,linetype='dotted', size=1,col = 'black')+
  scale_y_continuous(limits = c(0,2*alph))+
  scale_x_continuous(breaks=rho2_vec)+
  ylab('FDP')+theme_bw()+xlab(TeX('$\\rho_2$ when $\\rho_1=0.3'))+
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

plotdata1<- as.data.frame(c(colmeans(eprod_fdp[,6:10]),
                            colmeans(ehybrid_fdp[,6:10])))
names(plotdata1)<-"FDR"
plotdata1$type<-as.factor(rep(c('IRT*','IRT H'),each=length(rho2_vec)))
plotdata1$d<- rep(rho2_vec,2)
g2<-ggplot()+geom_line(data=plotdata1,aes(x=d,y=FDR,color=type),size=1)+
  geom_point(data=plotdata1,aes(x=d,y=FDR,color=type,shape=type),size=4,fill=NA)+
  geom_hline(yintercept = alph,linetype='dotted', size=1,col = 'black')+
  scale_y_continuous(limits = c(0,2*alph))+
  scale_x_continuous(breaks=rho2_vec)+
  ylab('FDP')+theme_bw()+xlab(TeX('$\\rho_2$ when $\\rho_1=0.5'))+
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

plotdata1<- as.data.frame(c(colmeans(eprod_fdp[,11:15]),
                            colmeans(ehybrid_fdp[,11:15])))
names(plotdata1)<-"FDR"
plotdata1$type<-as.factor(rep(c('IRT*','IRT H'),each=length(rho2_vec)))
plotdata1$d<- rep(rho2_vec,2)
g3<-ggplot()+geom_line(data=plotdata1,aes(x=d,y=FDR,color=type),size=1)+
  geom_point(data=plotdata1,aes(x=d,y=FDR,color=type,shape=type),size=4,fill=NA)+
  geom_hline(yintercept = alph,linetype='dotted', size=1,col = 'black')+
  scale_y_continuous(limits = c(0,2*alph))+
  scale_x_continuous(breaks=rho2_vec)+
  ylab('FDP')+theme_bw()+xlab(TeX('$\\rho_2$ when $\\rho_1=0.7'))+
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

# plotdata2<- as.data.frame(c(colmeans(eprod_etp),
#                             colmeans(ehybrid_etp)))
# names(plotdata2)<-"etp"
# plotdata2$type<-as.factor(rep(c('IRT*','IRT H'),each=length(rho2_vec)))
# plotdata2$d<- rep(rho2_vec,2)
# g2<-ggplot()+geom_line(data=plotdata2,aes(x=d,y=etp,color=type),size=1)+
#   geom_point(data=plotdata2,aes(x=d,y=etp,color=type,shape=type),size=4,fill=NA)+
#   ylab('ETP')+theme_bw()+xlab(TeX('$\\rho$'))+
#   scale_x_continuous(breaks=rho2_vec)+
#   theme(legend.position='top',legend.title=element_blank(),
#         legend.background = element_rect(fill="white",
#                                          size=1, linetype="solid", colour ="black"),
#         legend.text=element_text(size=12),
#         axis.text.x = element_text(size=15),
#         axis.text.y = element_text(size=15),
#         axis.title.x = element_text(size=15),
#         axis.title.y = element_text(size=15),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())

ggarrange(g1,g2,g3,ncol=3,nrow=1,common.legend = TRUE)

save.image('ablation_2.Rdata')
