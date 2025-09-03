

library(Rfast)
library(poolr)

source('funcs.R')

#case 1: studies have no FDR control
d = 5
m = 1000
reps = 2000
alph = 0.005
null_prop = 0.8
rho = 0
null_lb = 0.5
alphd_vec = c(0.02,0.05,0.08,0.1,0.12,0.15,0.18,0.2)

ehybrid_fdp = eprod_fdp = matrix(0,reps,length(alphd_vec))
ehybrid_etp = eprod_etp = matrix(0,reps,length(alphd_vec))

ei = null_prop+0.5*(1-null_prop)
for(dd in 1:length(alphd_vec)){
  nd = d
  d1 = floor(nd/2)
  d2 = nd-d1
  alphd = rep(alphd_vec[dd],nd)
  for(r in 1:reps){
    set.seed(r)#set.seed(r*dd)
    y = matrix(rnorm(m*nd,0,1),nd,m)
    sig = runif(nd,0.75,2)
    p = runif(m)
    mu = sapply(1:m,function(i) 0*(p[i]<=null_prop)+rnorm(1,-3,1)*(p[i]>null_prop & p[i]<=ei)+
                  rnorm(1,3,1)*(p[i]>ei))
    theta = 1*(p>null_prop)
    x = sapply(1:nd,function(i) sig[i]*y[i,]+mu)#sapply(1:nd,function(i) y[i,]+mu)
    pvals = sapply(1:nd,function(i) 2*pnorm(-abs(x[,i]/sig[i]),0,1))
    bh_decisions =  sapply(1:nd,function(i) bh.func(pvals[,i],alphd[i])$de)
    
    estat = m*sapply(1:nd,function(i) bh_decisions[,i]/(0.01*max(sum(bh_decisions[,i],1))))
    fedeval_mean = rowmeans(estat)
    ebh = ebh.func(fedeval_mean,alph)#ebh.func(fedeval_mean,alph)
    
    ###### ---- calculating the product e-values ---- #######
    estat_prod = matrix(0,m,1)
    for(ii in 1:m){
      ee = estat[ii,]
      if(max(ee)>0){
        tt = rep(0,nd)
        tt[1] = null_lb*mean(ee)
        tt[nd] = (null_lb^nd)*prod(ee)
        for(j in 2:(nd-1)){
          
          combs = comb_n(nd,j)
          ncols = ncol(combs)
          tt[j] = (null_lb^j)*mean(unlist(lapply(1:ncols,function(i) prod(ee[combs[,i]]))))
          
        }
        estat_prod[ii] = mean(tt[-1])
      }
    }
    
    ebh = ebh.func(estat_prod,alph)
    eprod_fdp[r,dd] = sum((1-theta)*ebh$de)/max(sum(ebh$de),1)
    eprod_etp[r,dd] = sum(theta*ebh$de)/max(sum(theta),1)
    ##########################################################
  
    ##########################################################
    ###### ---- calculating the hybrid e-values ---- #######
    estat_1 = m*sapply(1:d1,function(i) bh_decisions[,i]/(0.01*max(sum(bh_decisions[,i],1))))
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

plotdata1<- as.data.frame(c(colmeans(eprod_fdp),
                            colmeans(ehybrid_fdp)))
names(plotdata1)<-"FDR"
plotdata1$type<-as.factor(rep(c('IRT*','IRT H'),each=length(alphd_vec)))
plotdata1$d<- rep(alphd_vec,2)
g1<-ggplot()+geom_line(data=plotdata1,aes(x=d,y=FDR,color=type),size=1)+
  geom_point(data=plotdata1,aes(x=d,y=FDR,color=type,shape=type),size=4,fill=NA)+
  geom_hline(yintercept = alph,linetype='dotted', size=1,col = 'black')+
  scale_y_continuous(limits = c(0,4*alph))+
  scale_x_continuous(breaks=alphd_vec)+
  ylab('FDP')+theme_bw()+xlab(TeX('$\\alpha^*$'))+
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

plotdata2<- as.data.frame(c(colmeans(eprod_etp),
                            colmeans(ehybrid_etp)))
names(plotdata2)<-"etp"
plotdata2$type<-as.factor(rep(c('IRT*','IRT H'),each=length(alphd_vec)))
plotdata2$d<- rep(alphd_vec,2)
g2<-ggplot()+geom_line(data=plotdata2,aes(x=d,y=etp,color=type),size=1)+
  geom_point(data=plotdata2,aes(x=d,y=etp,color=type,shape=type),size=4,fill=NA)+
  ylab('ETP')+theme_bw()+xlab(TeX('$\\alpha^*$'))+
  scale_x_continuous(breaks=alphd_vec)+
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

save.image('ablation_1.Rdata')
