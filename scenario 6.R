
source('funcs.R')

library(Rfast)
library(poolr)

d = 30
m = 1000
nmax = 900
rratio = c(0.1,0.2,0.3,0.5,0.6)
nmin = ceiling(rratio*nmax)

alph = 0.15
null_prop = 0.8
rho = 0
reps = 500

fedeval_star_mean_fdp = p2e_mean_fdp = p2e_ni_mean_fdp = fedeval_mean_fdp = fed_naive_fdp = pooled_fdp = matrix(0,reps,length(rratio))
fedeval_star_mean_etp = p2e_mean_etp = p2e_ni_mean_etp = fedeval_mean_etp = fed_naive_etp = pooled_etp = matrix(0,reps,length(rratio))

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
    
    estat_n = estat_ni = count_decisions = matrix(0,m,1)
    for(ii in 1:m){
      
      estat_n[ii]=sum(unlist(sapply(1:d,
                                    function(i) (md[i]/max(n_i))*bh_decisions[[i]][which(md_id[[i]]==ii)]/(alphd[i]*max(sum(bh_decisions[[i]],1))))))
      
      estat_ni[ii]=sum(unlist(sapply(1:d,
                                    function(i) (md[i]/n_i[ii])*bh_decisions[[i]][which(md_id[[i]]==ii)]/(alphd[i]*max(sum(bh_decisions[[i]],1))))))
      
      count_decisions[ii] = sum(unlist(sapply(1:d,
                                              function(i) bh_decisions[[i]][which(md_id[[i]]==ii)])))
      
    
    }
    
    fed_ebh = ebh.func(estat_n,alph)$de
    fedeval_mean_fdp[r,k] = sum((1-theta)*fed_ebh)/max(sum(fed_ebh),1)
    fedeval_mean_etp[r,k] = sum(theta*fed_ebh)/max(sum(theta),1)
    
    fed_star_ebh = ebh.func(estat_ni,alph)$de
    fedeval_star_mean_fdp[r,k] = sum((1-theta)*fed_star_ebh)/max(sum(fed_star_ebh),1)
    fedeval_star_mean_etp[r,k] = sum(theta*fed_star_ebh)/max(sum(theta),1)
    
    naive_decisions = sapply(1:m,function(i) 1*(count_decisions[i]>=(n_i[i]/2)))
    fed_naive_fdp[r,k] = sum((1-theta)*naive_decisions)/max(sum(naive_decisions),1)
    fed_naive_etp[r,k] = sum(theta*naive_decisions)/max(sum(theta),1)
    
    pooled_pvals = sapply(1:m,function(i) fisher(pvals[i,h_id[[i]]])$p)
    bh_pooled = bh.func(pooled_pvals,alph)$de
    pooled_fdp[r,k] = sum((1-theta)*bh_pooled)/max(sum(bh_pooled),1)
    pooled_etp[r,k] = sum(theta*bh_pooled)/max(sum(theta),1)
    
    ##### p2e ###############################
    p2e = p2e_ni = matrix(0,m,1)
    for(i in 1:m){
      agents = which(h_id[[i]]==TRUE)
      temp1 = matrix(0,length(agents),1)
      temp2 = pvals[i,agents]
      id1 = (temp2==0)
      id2 = (temp2>0 & temp2<= exp(-2))
      temp1[id1] = Inf 
      temp1[id2] = 2/(temp2[id2]*((-log(temp2[id2]))^2))
      p2e[i] = sum(temp1)/max(n_i)
      p2e_ni[i] = sum(temp1)/n_i[i]
    }
    p2e_mean = rowmeans(p2e)
    ebh_p2e = ebh.func(p2e_mean,alph)
    p2e_mean_fdp[r,k] = sum((1-theta)*ebh_p2e$de)/max(sum(ebh_p2e$de),1)
    p2e_mean_etp[r,k] = sum(theta*ebh_p2e$de)/max(sum(theta),1)
    
    p2e_ni_mean = rowmeans(p2e_ni)
    ebh_p2e_ni = ebh.func(p2e_ni_mean,alph)
    p2e_ni_mean_fdp[r,k] = sum((1-theta)*ebh_p2e_ni$de)/max(sum(ebh_p2e_ni$de),1)
    p2e_ni_mean_etp[r,k] = sum(theta*ebh_p2e_ni$de)/max(sum(theta),1)
    print(r)
  }
}
rbind(c(colmeans(fedeval_mean_fdp)),c(colmeans(pooled_fdp)),c(colmeans(fed_naive_fdp)))
rbind(c(colmeans(fedeval_mean_etp)),c(colmeans(pooled_etp)),c(colmeans(fed_naive_etp)))

library(ggplot2)
library(ggpubr)
library(latex2exp)

plotdata1<- as.data.frame(c(colmeans(fedeval_mean_fdp),colmeans(fedeval_star_mean_fdp),
                            colmeans(p2e_mean_fdp),colmeans(p2e_ni_mean_fdp),
                            colmeans(pooled_fdp),
                            colmeans(fed_naive_fdp)))
names(plotdata1)<-"FDR"
plotdata1$type<-as.factor(rep(c('IRT','IRT*','P2E','P2E*','Fisher','Naive'),each=length(rratio)))
plotdata1$d<- rep(rratio,6)
g1<-ggplot()+geom_line(data=plotdata1,aes(x=d,y=FDR,color=type),size=1)+
  geom_point(data=plotdata1,aes(x=d,y=FDR,color=type,shape=type),size=4,fill=NA)+
  geom_hline(yintercept = 0.15,linetype='dotted', size=1,col = 'black')+
  scale_y_continuous(limits = c(0,0.3))+
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

plotdata2<- as.data.frame(c(colmeans(fedeval_mean_etp),colmeans(fedeval_star_mean_etp),
                            colmeans(p2e_mean_etp),colmeans(p2e_ni_mean_etp),colmeans(pooled_etp),
                            colmeans(fed_naive_etp)))
names(plotdata2)<-"etp"
plotdata2$type<-as.factor(rep(c('IRT','IRT*','P2E','P2E*','Fisher','Naive'),each=length(rratio)))
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

save.image('scenario 6.Rdata')

