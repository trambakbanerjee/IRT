

library(Rfast)
library(poolr)

source('funcs.R')

d = 10
mixprop = c(0.1,0.3,0.5,0.7,0.9)
m = 1000
reps = 200
alph = 0.01
null_prop = 0.8
rho = 0
null_lb = 0.5

ehybrid_fdp = e_fdp = eprod_fdp = p2e_fdp = p2eprod_fdp = pooled_fdp = matrix(0,reps,length(mixprop))
ehybrid_etp = e_etp = eprod_etp = p2e_etp = p2eprod_etp = pooled_etp = matrix(0,reps,length(mixprop))

nd = d
d1 = floor(nd/2)
d2 = nd-d1
alphd = rep(0.05,nd)
for(dd in 1:length(mixprop)){
  for(r in 1:reps){
    set.seed(r)
    p = runif(m)
    mu = sapply(1:m,function(i) 0*(p[i]<=null_prop)+rnorm(1,-3,1)*(p[i]>null_prop & p[i]<=0.9)+
                  rnorm(1,3,1)*(p[i]>0.9))
    theta = 1*(p>null_prop)
    set.seed(r+1)
    p = runif(m*nd)
    y = matrix(sapply(1:(m*nd),function(i) 
      (p[i]<=mixprop[dd])*rnorm(1)+
        (p[i]>mixprop[dd])*rt(1,5)),nd,m)
    
    x = sapply(1:nd,function(i) y[i,]+mu)#sapply(1:nd,function(i) y[i,]+mu)
    sig = sqrt(rep(mixprop[dd]+(1-mixprop[dd])*(5/3),nd))
    pvals = sapply(1:nd,function(i) 2*pnorm(-abs(x[,i]/sig[i]),0,1))
    bh_decisions =  sapply(1:nd,function(i) bh.func(pvals[,i],alphd[i])$de)
    
    estat = m*sapply(1:nd,function(i) bh_decisions[,i]/(alphd[i]*max(sum(bh_decisions[,i],1))))
    fedeval_mean = rowmeans(estat)
    ebh = ebh.func(fedeval_mean,alph)
    e_fdp[r,dd] = sum((1-theta)*ebh$de)/max(sum(ebh$de),1)
    e_etp[r,dd] = sum(theta*ebh$de)/max(sum(theta),1)
    
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
    p2e_fdp[r,dd] = sum((1-theta)*ebh_p2e$de)/max(sum(ebh_p2e$de),1)
    p2e_etp[r,dd] = sum(theta*ebh_p2e$de)/max(sum(theta),1)
    ###########################################
    ###### ---- calculating the product P2E ---- #######
    p2e_prod = matrix(0,m,1)
    for(ii in 1:m){
      ee = p2e[ii,]
      if(max(ee)>0){
        tt = rep(0,nd)
        tt[1] = mean(ee)
        tt[nd] = prod(ee)
        for(j in 2:(nd-1)){
          
          combs = comb_n(nd,j)
          ncols = ncol(combs)
          tt[j] = mean(unlist(lapply(1:ncols,function(i) prod(ee[combs[,i]]))))
          
        }
        p2e_prod[ii] = mean(tt[-1])
      }
    }
    
    ebh_p2eprod = ebh.func(p2e_prod,alph)
    p2eprod_fdp[r,dd] = sum((1-theta)*ebh_p2eprod$de)/max(sum(ebh_p2eprod$de),1)
    p2eprod_etp[r,dd] = sum(theta*ebh_p2eprod$de)/max(sum(theta),1)
    ###########################################
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

plotdata1<- as.data.frame(c(colmeans(e_fdp),colmeans(eprod_fdp),
                            colmeans(ehybrid_fdp),
                            colmeans(p2e_fdp),colmeans(p2eprod_fdp),
                            colmeans(pooled_fdp)))
names(plotdata1)<-"FDR"
plotdata1$type<-as.factor(rep(c('IRT','IRT*','IRT H','P2E','P2E*','Fisher'),each=length(mixprop)))
plotdata1$d<- rep(mixprop,6)
g1<-ggplot()+geom_line(data=plotdata1,aes(x=d,y=FDR,color=type),size=1)+
  geom_point(data=plotdata1,aes(x=d,y=FDR,color=type,shape=type),size=4,fill=NA)+
  geom_hline(yintercept = alph,linetype='dotted', size=1,col = 'black')+
  scale_x_continuous(breaks = mixprop)+
  ylab('FDP')+theme_bw()+xlab('q')+
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

plotdata2<- as.data.frame(c(colmeans(e_etp),colmeans(eprod_etp),
                            colmeans(ehybrid_etp),
                            colmeans(p2e_etp),colmeans(p2eprod_etp),
                            colmeans(pooled_etp)))
names(plotdata2)<-"etp"
plotdata2$type<-as.factor(rep(c('IRT','IRT*','IRT H','P2E','P2E*','Fisher'),each=length(mixprop)))
plotdata2$d<- rep(mixprop,6)
g2<-ggplot()+geom_line(data=plotdata2,aes(x=d,y=etp,color=type),size=1)+
  geom_point(data=plotdata2,aes(x=d,y=etp,color=type,shape=type),size=4,fill=NA)+
  scale_x_continuous(breaks = mixprop)+
  ylab('ETP')+theme_bw()+xlab('q')+
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

save.image('scenario 8_1.Rdata')
