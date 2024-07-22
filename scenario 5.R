
source('funcs.R')

library(Rfast)
library(poolr)


d = 10
d1 = 1:floor(d/2)
d2 = (floor(d/2)+1):d
m = 1000
nmax = 900
rratio = c(0.1,0.2,0.3,0.4,0.5)#nmin/nmax
nmin = ceiling(rratio*nmax)

alph = 0.03
null_prop = 0.8
null_lb = 0.5
rho = 0
reps = 200

ehybrid_fdp = e_fdp = p2e_fdp = eprod_fdp = pooled_fdp = matrix(0,reps,length(rratio))
ehybrid_etp = e_etp = p2e_etp = eprod_etp = pooled_etp = matrix(0,reps,length(rratio))

for(k in 1:length(rratio)){
  kk = nmin[k]
  set.seed(1)
  md = sample(kk:nmax,d,replace = TRUE)
  md_id = vector("list", length = d)
  for(i in 1:d){
    set.seed(i)
    md_id[[i]] = sample(1:m,md[i],replace = FALSE)
  }
  alphd = sapply(1:d,function(i) 0.05*(md[i]<=600)+0.03*(md[i]>600 & md[i]<=800)+
                   0.01*(md[i]>800))
  ############## sanity ######################
  n_i = matrix(0,m,1)
  h_id = list()
  for(ii in 1:m){
    n_i[ii] = sum(sapply(1:d,function(i) sum(md_id[[i]]==ii)))
    h_id[[ii]] = sapply(1:d,function(i) ii%in%md_id[[i]])
  }
  if(min(n_i)==0){
    break
  }
  ###########################################
  
  for(r in 1:reps){
    set.seed(r*kk)
    y = matrix(rnorm(m*d,0,1),d,m)
    p = runif(m)
    mu = sapply(1:m,function(i) 0*(p[i]<=null_prop)+rnorm(1,-3,1)*(p[i]>null_prop & p[i]<=0.9)+
                  rnorm(1,3,1)*(p[i]>0.9))
    theta = 1*(p>null_prop)
    x = sapply(1:d,function(i) y[i,]+mu)
    pvals = sapply(1:d,function(i) 2*pnorm(-abs(x[,i]),0,1))
    bh_decisions =  sapply(1:d,function(i) bh.func(pvals[md_id[[i]],i],alphd[i])$de)
    
    estat = estat_hybrid = estat_prod = matrix(0,m,1)
    cc = max(1/max(n_i),m/sum(md))
    for(ii in 1:m){
      
      
      estat[ii]=sum(unlist(sapply(1:d,
                                  function(i) (md[i]*cc)*bh_decisions[[i]][which(md_id[[i]]==ii)]/(alphd[i]*max(sum(bh_decisions[[i]],1))))))
      
      
      #---- calculating the product e-values ----------------------------------
      ee = null_lb*unlist(sapply(1:d,
                                 function(i) md[i]*bh_decisions[[i]][which(md_id[[i]]==ii)]/(alphd[i]*max(sum(bh_decisions[[i]],1)))))
      
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
    
    fed_ebh = ebh.func(estat,alph)$de
    e_fdp[r,k] = sum((1-theta)*fed_ebh)/max(sum(fed_ebh),1)
    e_etp[r,k] = sum(theta*fed_ebh)/max(sum(theta),1)
    
    fed_ebh_prod = ebh.func(estat_prod,alph)$de
    eprod_fdp[r,k] = sum((1-theta)*fed_ebh_prod)/max(sum(fed_ebh_prod),1)
    eprod_etp[r,k] = sum(theta*fed_ebh_prod)/max(sum(theta),1)
    
    pooled_pvals = sapply(1:m,function(i) fisher(pvals[i,h_id[[i]]])$p)
    bh_pooled = bh.func(pooled_pvals,alph)$de
    pooled_fdp[r,k] = sum((1-theta)*bh_pooled)/max(sum(bh_pooled),1)
    pooled_etp[r,k] = sum(theta*bh_pooled)/max(sum(theta),1)
    
    ##### p2e ###############################
    p2e_ni = matrix(0,m,1)
    for(i in 1:m){
      agents = which(h_id[[i]]==TRUE)
      temp1 = matrix(0,length(agents),1)
      temp2 = pvals[i,agents]
      id1 = (temp2==0)
      id2 = (temp2>0 & temp2<= exp(-2))
      temp1[id1] = Inf 
      temp1[id2] = 2/(temp2[id2]*((-log(temp2[id2]))^2))
      p2e_ni[i] = sum(temp1)/n_i[i]
    }
    
    p2e_ni_mean = c(p2e_ni)
    ebh_p2e_ni = ebh.func(p2e_ni_mean,alph)
    p2e_fdp[r,k] = sum((1-theta)*ebh_p2e_ni$de)/max(sum(ebh_p2e_ni$de),1)
    p2e_etp[r,k] = sum(theta*ebh_p2e_ni$de)/max(sum(theta),1)
    ##########################################################
    ###### ---- calculating the hybrid e-values ---- #######
    ehybrid = matrix(0,m,1)
    cc = m/sum(sapply(1:length(d1),function(i) md[i]))
    for(ii in 1:m){
      agents = which(h_id[[ii]]==TRUE)
      agents_1 = agents[which(agents%in%d1)]
      agents_2 = agents[which(agents%in%d2)]
      
      estat_1 = NA
      estat_prod_2 = NA
      if(length(agents_1)>0){
        estat_1 = sum(unlist(sapply(1:length(agents_1),
                                    function(i) (md[agents_1[i]]*cc)*bh_decisions[[agents_1[i]]][which(md_id[[agents_1[i]]]==ii)]/(alphd[agents_1[i]]*max(sum(bh_decisions[[agents_1[i]]],1))))))
      }
      if(length(agents_2)>0){
        estat_prod_2 = 0
        ee = null_lb*unlist(sapply(1:length(agents_2),
                                   function(i) md[agents_2[i]]*bh_decisions[[agents_2[i]]][which(md_id[[agents_2[i]]]==ii)]/(alphd[agents_2[i]]*max(sum(bh_decisions[[agents_2[i]]],1)))))
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
    ebh = ebh.func(ehybrid,alph)
    ehybrid_fdp[r,k] = sum((1-theta)*ebh$de)/max(sum(ebh$de),1)
    ehybrid_etp[r,k] = sum(theta*ebh$de)/max(sum(theta),1)
    ##########################################################
    print(r)
  }
}
save.image('scenario 5.Rdata')

library(ggplot2)
library(ggpubr)
library(latex2exp)
load("scenario 5_p2estar.Rdata")
plotdata1<- as.data.frame(c(colmeans(e_fdp),
                            colmeans(eprod_fdp),
                            colmeans(ehybrid_fdp),
                            colmeans(p2e_fdp),colmeans(p2eprod_fdp),
                            colmeans(pooled_fdp)))
names(plotdata1)<-"FDR"
plotdata1$type<-as.factor(rep(c('IRT','IRT*','IRT H','P2E','P2E*','Fisher'),
                              each=length(rratio)))
plotdata1$d<- rep(rratio,6)
g1<-ggplot()+geom_line(data=plotdata1,aes(x=d,y=FDR,color=type),size=1)+
  geom_point(data=plotdata1,aes(x=d,y=FDR,color=type,shape=type),size=4,fill=NA)+
  geom_hline(yintercept = alph,linetype='dotted', size=1,col = 'black')+
  scale_y_continuous(limits = c(0,2*alph))+
  ylab('FDP')+theme_bw()+xlab(TeX('$\\eta$'))+
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
plotdata2$type<-as.factor(rep(c('IRT','IRT*','IRT H','P2E','P2E*','Fisher'),each=length(rratio)))
plotdata2$d<- rep(rratio,6)
g2<-ggplot()+geom_line(data=plotdata2,aes(x=d,y=etp,color=type),size=1)+
  geom_point(data=plotdata2,aes(x=d,y=etp,color=type,shape=type),size=4,fill=NA)+
  ylab('ETP')+theme_bw()+xlab(TeX('$\\eta$'))+
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




