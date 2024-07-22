
source('funcs.R')

library(Rfast)
d = 10
m = 1000
set.seed(1)
md = sample(c(2000,2500,3000),d,replace = TRUE)
alphd = sapply(1:d,function(i) 0.05*(md[i]<=1500)+0.03*(md[i]>1500 & md[i]<=2500)+
                 0.01*(md[i]>2500))
set.seed(2)
h_id = sample(1000:max(md),m,replace = FALSE)#1000

alph = c(0.001,0.003,0.005,0.01,0.03,0.05,0.1)
null_prop = 0.8
rho = 0
reps = 200

p2eprod_fdp = matrix(0,reps,length(alph))
p2eprod_etp = matrix(0,reps,length(alph))

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
    pvals = list()
    for(i in 1:d){
      set.seed(r*i)
      y = rnorm(md[i],0,1)
      x = y+mu[1:md[i]]
      pvals[[i]] = 2*pnorm(-abs(x),0,1)
    }
    
    p2e_prod = matrix(0,m,1)
    for(ii in 1:m){
      
      pvals_ii = unlist(sapply(1:d,function(i) pvals[[i]][which(md_id[[i]]==h_id[ii])]))
      temp1 = matrix(0,n_i[ii],1)
      temp2 = pvals_ii
      id1 = (temp2==0)
      id2 = (temp2>0 & temp2<= exp(-2))
      temp1[id1] = Inf 
      temp1[id2] = 2/(temp2[id2]*((-log(temp2[id2]))^2))
      
      #---- calculating the product p2e-values ----------------------------------
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
    ebh_p2eprod = ebh.func(p2e_prod,alph[k])
    p2eprod_fdp[r,k] = sum((1-theta)*ebh_p2eprod$de)/max(sum(ebh_p2eprod$de),1)
    p2eprod_etp[r,k] = sum(theta*ebh_p2eprod$de)/max(sum(theta),1)
    print(r)
  }
}
save.image('scenario 6_p2estar.Rdata')


