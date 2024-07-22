

source('funcs.R')

library(Rfast)
library(poolr)

d = 10
m = 1000
nmax = 900
rratio = c(0.1,0.2,0.3,0.4,0.5)#nmin/nmax
nmin = ceiling(rratio*nmax)

alph = 0.03
null_prop = 0.8
rho = 0
reps = 200

p2eprod_fdp = matrix(0,reps,length(rratio))
p2eprod_etp = matrix(0,reps,length(rratio))

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
    
    ##### p2e prod ###############################
    p2e_prod = matrix(0,m,1)
    for(ii in 1:m){
      agents = which(h_id[[ii]]==TRUE)
      temp1 = matrix(0,length(agents),1)
      temp2 = pvals[ii,agents]
      id1 = (temp2==0)
      id2 = (temp2>0 & temp2<= exp(-2))
      temp1[id1] = Inf 
      temp1[id2] = 2/(temp2[id2]*((-log(temp2[id2]))^2))
      
      ee = temp1
      d_ii = n_i[ii]
      p2e_prod[ii] = 0
      if(max(ee)>0){
        if(d_ii == 1){
          p2e_prod[ii] = ee
        } else if(d_ii == 2){
          p2e_prod[ii] = mean(mean(ee),prod(ee))
        } else{
          tt = rep(0,d_ii)
          tt[1] = mean(ee)
          tt[d_ii] = prod(ee)
          for(j in 2:(d_ii-1)){
            combs = comb_n(d_ii,j)
            ncols = ncol(combs)
            tt[j] = mean(sapply(1:ncols,function(i) prod(ee[combs[,i]])))
          }
          p2e_prod[ii] = mean(tt[-1])
        }
      }
    }
    
    ebh_p2eprod = ebh.func(p2e_prod,alph)
    p2eprod_fdp[r,k] = sum((1-theta)*ebh_p2eprod$de)/max(sum(ebh_p2eprod$de),1)
    p2eprod_etp[r,k] = sum(theta*ebh_p2eprod$de)/max(sum(theta),1)
    ##########################################################
    print(r)
  }
}

save.image('scenario 5_p2estar.Rdata')


