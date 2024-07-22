
source('funcs.R')

library(Rfast)
library(poolr)

d = 10
d1 = 1:floor(d/2)
d2 = (floor(d/2)+1):d
m = 1000
set.seed(1)
md = sample(c(2000,2500,3000),d,replace = TRUE)
alphd = sapply(1:d,function(i) 0.05*(md[i]<=1500)+0.03*(md[i]>1500 & md[i]<=2500)+
                 0.01*(md[i]>2500))
set.seed(2)
h_id = sample(1000:max(md),m,replace = FALSE)

alph = c(0.001,0.003,0.005,0.01,0.03,0.05,0.1)
null_prop = 0.8
null_lb = 0.5
rho = 0
reps = 200

ehybrid_fdp = e_fdp = estar_fdp = p2e_fdp = eprod_fdp = pooled_fdp = matrix(0,reps,length(alph))
ehybrid_etp = e_etp = estar_etp = p2e_etp = eprod_etp = pooled_etp = matrix(0,reps,length(alph))

############## sanity ######################
n_i = matrix(0,m,1)
md_id = sapply(1:d,function(i) 1:md[i])
for(ii in 1:m){
  n_i[ii] = sum(sapply(1:d,function(i) 1*(h_id[ii]%in%md_id[[i]])))
}
if(min(n_i)==0){
  break
}
###########################################


for(k in 1:length(alph)){
  for(r in 1:reps){
    set.seed(r)
    p = runif(max(md))
    mu = sapply(1:max(md),function(i) 0*(p[i]<=null_prop)+rnorm(1,-3,1)*(p[i]>null_prop & p[i]<=0.9)+
                  rnorm(1,3,1)*(p[i]>0.9))
    theta = 1*(p>null_prop)
    theta = theta[h_id]
    bh_decisions = pvals = list()
    for(i in 1:d){
      set.seed(r*i)
      y = rnorm(md[i],0,1)
      x = y+mu[1:md[i]]
      pvals[[i]] = 2*pnorm(-abs(x),0,1)
      bh_decisions[[i]] = bh.func(pvals[[i]],alphd[i])$de
    }
    
    estat = estat_prod = pooled_pvals = p2e = matrix(0,m,1)
    for(ii in 1:m){
      
      cc = max(m/sum(md),1/max(n_i))
      estat[ii]=sum(unlist(sapply(1:d,
                                  function(i) (md[i]*cc)*bh_decisions[[i]][which(md_id[[i]]==h_id[ii])]/(alphd[i]*max(sum(bh_decisions[[i]],1))))))
      
      pvals_ii = unlist(sapply(1:d,function(i) pvals[[i]][which(md_id[[i]]==h_id[ii])]))
      pooled_pvals[ii] = fisher(pvals_ii)$p
      
      temp1 = matrix(0,n_i[ii],1)
      temp2 = pvals_ii
      id1 = (temp2==0)
      id2 = (temp2>0 & temp2<= exp(-2))
      temp1[id1] = Inf 
      temp1[id2] = 2/(temp2[id2]*((-log(temp2[id2]))^2))
      p2e[ii] = sum(temp1)/n_i[ii]
      
      #---- calculating the product e-values ----------------------------------
      ee = null_lb*unlist(sapply(1:d,
                                 function(i) md[i]*bh_decisions[[i]][which(md_id[[i]]==h_id[ii])]/(alphd[i]*max(sum(bh_decisions[[i]],1)))))
      
      d_ii = n_i[ii]
      estat_prod[ii] = 0
      if(max(ee)>0){
        if(d_ii == 1){
          estat_prod[ii] = ee
        } else if(d_ii == 2){
          estat_prod[ii] = mean(mean(ee),prod(ee))
        } else{
          tt = rep(0,d_ii)
          tt[1] = mean(ee)
          tt[d_ii] = prod(ee)
          for(j in 2:(d_ii-1)){
            combs = comb_n(d_ii,j)
            ncols = ncol(combs)
            tt[j] = mean(sapply(1:ncols,function(i) prod(ee[combs[,i]])))
          }
          estat_prod[ii] = mean(tt[-1])
        }
      }
      
    }
    
    fed_ebh = ebh.func(estat,alph[k])$de
    e_fdp[r,k] = sum((1-theta)*fed_ebh)/max(sum(fed_ebh),1)
    e_etp[r,k] = sum(theta*fed_ebh)/max(sum(theta),1)
    
    fed_ebh_prod = ebh.func(estat_prod,alph[k])$de
    eprod_fdp[r,k] = sum((1-theta)*fed_ebh_prod)/max(sum(fed_ebh_prod),1)
    eprod_etp[r,k] = sum(theta*fed_ebh_prod)/max(sum(theta),1)
    
    bh_pooled = bh.func(pooled_pvals,alph[k])$de
    pooled_fdp[r,k] = sum((1-theta)*bh_pooled)/max(sum(bh_pooled),1)
    pooled_etp[r,k] = sum(theta*bh_pooled)/max(sum(theta),1)
    
    ebh_p2e_ni = ebh.func(p2e,alph[k])
    p2e_fdp[r,k] = sum((1-theta)*ebh_p2e_ni$de)/max(sum(ebh_p2e_ni$de),1)
    p2e_etp[r,k] = sum(theta*ebh_p2e_ni$de)/max(sum(theta),1)
    
    ## Calculating the hybrid e-values
    ehybrid = matrix(0,m,1)
    cc = m/sum(sapply(1:length(d1),function(i) md[i]))
    for(ii in 1:m){
      agents = which(sapply(1:d,function(i) i*(h_id[ii]%in%md_id[[i]]))>0)
      agents_1 = agents[which(agents%in%d1)]
      agents_2 = agents[which(agents%in%d2)]
      estat_1 = NA
      estat_prod_2 = NA
      if(length(agents_1)>0){
        estat_1 = sum(unlist(sapply(1:length(agents_1),
                                    function(i) (md[agents_1[i]]*cc)*bh_decisions[[agents_1[i]]][which(md_id[[agents_1[i]]]==h_id[ii])]/(alphd[agents_1[i]]*max(sum(bh_decisions[[agents_1[i]]],1))))))
        
      }
      if(length(agents_2)>0){
        estat_prod_2 = 0
        ee = null_lb*unlist(sapply(1:length(agents_2),
                                   function(i) md[agents_2[i]]*bh_decisions[[agents_2[i]]][which(md_id[[agents_2[i]]]==h_id[ii])]/(alphd[agents_2[i]]*max(sum(bh_decisions[[agents_2[i]]],1)))))
        if(max(ee)>0){
          if(length(agents_2) == 1){
            estat_prod_2 = ee
          } else if(length(agents_2) == 2){
            estat_prod_2 = mean(mean(ee),prod(ee))
          } else{
            tt = rep(0,length(agents_2))
            tt[1] = mean(ee)
            tt[length(agents_2)] = prod(ee)
            for(j in 2:(length(agents_2)-1)){
              
              combs = comb_n(length(agents_2),j)
              ncols = ncol(combs)
              tt[j] = mean(unlist(lapply(1:ncols,function(i) prod(ee[combs[,i]]))))
              
            }
            estat_prod_2 = mean(tt[-1])
          }
        }
      }
      ehybrid[ii] = NA
      if(length(agents_1)>0 & length(agents_2)>0){
        ehybrid[ii] = (1/length(agents))*(length(agents_1)*estat_1+length(agents_2)*estat_prod_2)
        
      }
      if(length(agents_1)==0){
        ehybrid[ii] = estat_prod_2
      }
      if(length(agents_2)==0){
        ehybrid[ii] = estat_1
      }
      
    }
    
    ebh = ebh.func(ehybrid,alph[k])
    ehybrid_fdp[r,k] = sum((1-theta)*ebh$de)/max(sum(ebh$de),1)
    ehybrid_etp[r,k] = sum(theta*ebh$de)/max(sum(theta),1)
    print(r)
  }
}
save.image('scenario 6.Rdata')
load("scenario 6_p2estar.Rdata")

library(ggplot2)
library(ggpubr)
library(latex2exp)

plotdata1<- as.data.frame(c(colmeans(e_fdp),
                            colmeans(eprod_fdp),
                            colmeans(ehybrid_fdp),
                            colmeans(p2e_fdp),colmeans(p2eprod_fdp),
                            colmeans(pooled_fdp)))
names(plotdata1)<-"FDR"
plotdata1$type<-as.factor(rep(c('IRT','IRT*','IRT H','P2E','P2E*','Fisher'),
                              each=length(alph)))
plotdata1$d<- rep(log10(alph),6)
g1<-ggplot()+geom_line(data=plotdata1,aes(x=d,y=log10(FDR),color=type),size=1)+
  geom_point(data=plotdata1,aes(x=d,y=log10(FDR),color=type,shape=type),size=4,fill=NA)+
  #  geom_abline(intercept = alph[1],slope =1)+
  geom_hline(yintercept = log10(alph[1]),linetype='dotted', size=1,col = 'black')+
  geom_hline(yintercept = log10(alph[2]),linetype='dotted', size=1,col = 'black')+
  geom_hline(yintercept = log10(alph[3]),linetype='dotted', size=1,col = 'black')+
  geom_hline(yintercept = log10(alph[4]),linetype='dotted', size=1,col = 'black')+
  geom_hline(yintercept = log10(alph[5]),linetype='dotted', size=1,col = 'black')+
  geom_hline(yintercept = log10(alph[6]),linetype='dotted', size=1,col = 'black')+
  geom_hline(yintercept = log10(alph[7]),linetype='dotted', size=1,col = 'black')+
  #scale_y_continuous(limits = c(0,alph[5]))+
  ylab(TeX('$\\log_{10}(FDP)$'))+theme_bw()+xlab(TeX('$\\log_{10}(\\alpha)$'))+
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
                            colmeans(eprod_etp),
                            colmeans(ehybrid_etp),
                            colmeans(p2e_etp),colmeans(p2eprod_etp),
                            colmeans(pooled_etp)))
names(plotdata2)<-"etp"
plotdata2$type<-as.factor(rep(c('IRT','IRT*','IRT H','P2E','P2E*','Fisher'),each=length(alph)))
plotdata2$d<- rep(log10(alph),6)
g2<-ggplot()+geom_line(data=plotdata2,aes(x=d,y=etp,color=type),size=1)+
  geom_point(data=plotdata2,aes(x=d,y=etp,color=type,shape=type),size=4,fill=NA)+
  ylab('ETP')+theme_bw()+xlab(TeX('$\\log_{10}(\\alpha)$'))+
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



