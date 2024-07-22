
library(Rfast)
library(limma)
library(Biobase)
source('funcs.R')
load(file = "prostate8.rda")

pvalmat = list()
glist = list()
#Study 1
welsh = prostate8$data$Welsh
glist[[1]] = rownames(welsh)
y = ExpressionSet(welsh)
ourData <- y
ourData$Type <- factor(prostate8$dataLabel$Welsh)
design <- model.matrix(~ ourData$Type)
fit <- lmFit(ourData, design)
fit <- eBayes(fit)
tT <- topTable(fit, number=nrow(fit))
pvalmat[[1]] = tT$P.Value
#Study 2
yu = prostate8$data$Yu
glist[[2]] = rownames(yu)
y = ExpressionSet(yu)
ourData <- y
ourData$Type <- factor(prostate8$dataLabel$Yu)
design <- model.matrix(~ ourData$Type)
fit <- lmFit(ourData, design)
fit <- eBayes(fit)
tT <- topTable(fit, number=nrow(fit))
pvalmat[[2]] = tT$P.Value
#Study 3
lp = prostate8$data$Lapointe
glist[[3]] = rownames(lp)
y = ExpressionSet(lp)
ourData <- y
ourData$Type <- factor(prostate8$dataLabel$Lapointe)
design <- model.matrix(~ ourData$Type)
fit <- lmFit(ourData, design)
fit <- eBayes(fit)
tT <- topTable(fit, number=nrow(fit))
pvalmat[[3]] = tT$P.Value
#Study 4
varam = prostate8$data$Varambally
glist[[4]] = rownames(varam)
y = ExpressionSet(varam)
ourData <- y
ourData$Type <- factor(prostate8$dataLabel$Varambally)
design <- model.matrix(~ ourData$Type)
fit <- lmFit(ourData, design)
fit <- eBayes(fit)
tT <- topTable(fit, number=nrow(fit))
pvalmat[[4]] = tT$P.Value
#Study 5
singh = prostate8$data$Singh
glist[[5]] = rownames(singh)
y = ExpressionSet(singh)
ourData <- y
ourData$Type <- factor(prostate8$dataLabel$Singh)
design <- model.matrix(~ ourData$Type)
fit <- lmFit(ourData, design)
fit <- eBayes(fit)
tT <- topTable(fit, number=nrow(fit))
pvalmat[[5]] = tT$P.Value
#Study 6
wal = prostate8$data$Wallace
glist[[6]] = rownames(wal)
y = ExpressionSet(wal)
ourData <- y
ourData$Type <- factor(prostate8$dataLabel$Wallace)
design <- model.matrix(~ ourData$Type)
fit <- lmFit(ourData, design)
fit <- eBayes(fit)
tT <- topTable(fit, number=nrow(fit))
pvalmat[[6]] = tT$P.Value
#Study 7
nani = prostate8$data$Nanni
glist[[7]] = rownames(nani)
y = ExpressionSet(nani)
ourData <- y
ourData$Type <- factor(prostate8$dataLabel$Nanni)
design <- model.matrix(~ ourData$Type)
fit <- lmFit(ourData, design)
fit <- eBayes(fit)
tT <- topTable(fit, number=nrow(fit))
pvalmat[[7]] = tT$P.Value
#Study 8
tom = prostate8$data$Tomlins
glist[[8]] = rownames(tom)
y = ExpressionSet(tom)
ourData <- y
ourData$Type <- factor(prostate8$dataLabel$Tomlins)
design <- model.matrix(~ ourData$Type)
fit <- lmFit(ourData, design)
fit <- eBayes(fit)
tT <- topTable(fit, number=nrow(fit))
pvalmat[[8]] = tT$P.Value

#---------------------------Basic IRT and P2E ----------------------------------
estat = list()
n = data.frame(table(c(unlist(glist))))
n = n[n$Freq>=1,]
m = length(n$Var1)
g = as.character(unique(n$Var1))
nmax = max(n$Freq)
d = 8
md_id = sapply(1:d,function(i) glist[[i]][which(glist[[i]]%in%n$Var1)])
md = sapply(1:d,function(i) length(md_id[[i]]))
alphd = c(0.01,0.05,0.05,0.01,0.05,0.05,0.01,0.01)
alph = 0.1

bh_pvalscutoff = sapply(1:d,function(i) bh.func(pvalmat[[i]][which(glist[[i]]%in%n$Var1)],alphd[i])$th)
bh_decisions =  sapply(1:d,function(i) bh.func(pvalmat[[i]][which(glist[[i]]%in%n$Var1)],alphd[i])$de)

estat = sapply(1:d,function(i) (length(bh_decisions[[i]])/alphd[i])*bh_decisions[[i]]/max(sum(bh_decisions[[i]]),1))
e_agg = n_rejections = matrix(0,m,1)
for(i in 1:m){
  gi = g[i]
  e_agg[i]=(m/sum(md))*sum(unlist(sapply(1:d,
                             function(i) md[i]*bh_decisions[[i]][which(md_id[[i]]==gi)]/(alphd[i]*max(sum(bh_decisions[[i]],1))))))
  n_rejections[i] = sum(unlist(sapply(1:d,function(i) (bh_decisions[[i]][which(md_id[[i]]==gi)]))))
}
ebh = ebh.func(e_agg,alph)$de
imt_rejections =sum(ebh)
sum(1*(n_rejections>0))
max(imt_rejections)


p2e = matrix(0,m,1)
for(i in 1:m){
  gi = g[i]
  temp1 = unlist(sapply(1:d,function(i) pvalmat[[i]][which(md_id[[i]]==gi)]))
  temp2 = matrix(0,length(temp1),1)
  id1 = (temp1==0)
  id2 = (temp1>0 & temp1<= exp(-2))
  temp2[id1] = Inf 
  temp2[id2] = 2/(temp1[id2]*((-log(temp1[id2]))^2))
  p2e[i] = sum(temp2)/n$Freq[i]
}
ebh_p2e = ebh.func(p2e,alph)$de

hyp_p2e_alpha0.1 = g[which(ebh_p2e==1)]
hyp_alpha0.1 = g[which(ebh==1)]
cc_p2e= sapply(1:8,function(i) which(hyp_p2e_alpha0.1%in%glist[[i]][bh_decisions[[i]]==1]))
cc_estat= sapply(1:8,function(i) which(hyp_alpha0.1%in%glist[[i]][bh_decisions[[i]]==1]))

dd1 = sapply(1:d,function(i) sum(bh_decisions[[i]]))
dd3 = sapply(1:d,function(i) max(estat[[i]]))

###Plots
library(ggplot2)
library(ggpubr)
library(latex2exp)
df1 = data.frame('x'=n$Freq)
g1 = ggplot(df1)+geom_histogram(aes(x=x),color="black", fill="gray")+
  xlab(TeX(r'($n_i$)'))+
  theme_light()+
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))

study = matrix(FALSE,8,8)
for(j in 1:8){
  idd = n$Var1[n$Freq==j]
  study[j,] = (md/alphd)*sapply(1:8,function(i) length(which(glist[[i]]%in%idd))>0)
  
}
nn = 100
set.seed(1)
gene_sample = sample(g,nn,replace = FALSE)
estat_heatmap = matrix(NA,nn,8)
for(i in 1:length(gene_sample)){
  idd = gene_sample[i]
  for(j in 1:8){
    if(length(which(glist[[j]]%in%idd))>0){
      iy = which(glist[[j]]%in%idd)
      estat_heatmap[i,j] = max(0,log(estat[[j]][iy]))
    }
  }
}
library(lattice)
library(viridisLite)
x.scale <- list(at=seq(1,8,1),labels=c(2,3,4,5,1,8,7,6))#list(at=seq(1,8,1))

coul <- turbo(20,direction = -1)
p1<-levelplot(t(estat_heatmap),
              at=seq(0, 9, length.out=20),
              col.regions=coul,aspect="fill",region=TRUE,
              scales=list(x=x.scale),xlab="Studies",
              ylab="100 Genes")

plot(g1)
plot(p1)

top20 = order(e_agg,decreasing = TRUE)[1:25]
z = matrix(0,25,1)
for(i in 1:25){
  
  z[i] = n$Freq[which(as.character(n$Var1)%in%g[top20][i])]
  
}
#g[top100]
df3 = data.frame('y'=log(e_agg[top20]),'x'=1:25,'z'=as.factor(z))
g3 = ggplot(df3)+geom_point(aes(x=x,y=y,group=z,color=z,shape=z),size=3)+
  scale_shape_manual(values=1:nlevels(df3$z))+
  ylab(TeX(r'($\log(e^{agg})$)'))+xlab('Top 25 Genes')+
  scale_x_continuous(breaks=1:25,
                   labels=g[top20])+
  theme_light()+ 
  theme(axis.text.x = element_text(size=8,angle = 90,hjust=1,
                                   vjust = 1),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        legend.title = element_blank())

df4 = data.frame('x'=log(e_agg[e_agg>0]))
g4 = ggplot(df4)+geom_histogram(aes(x=x),color="black", fill="gray")+
  xlab(TeX(r'($\log(e^{agg})$ when $e^{agg}>0$)'))+
  theme_light()+
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))

ggarrange(g4,g3,ncol=2)

top20 = order(p2e,decreasing = TRUE)[1:25]
z = matrix(0,25,1)
for(i in 1:25){
  
  z[i] = n$Freq[which(as.character(n$Var1)%in%g[top20][i])]
  
}
#g[top100]
df3 = data.frame('y'=log(p2e[top20]),'x'=1:25,'z'=as.factor(z))
g3 = ggplot(df3)+geom_point(aes(x=x,y=y,group=z,color=z,shape=z),size=3)+
  scale_shape_manual(values=1:nlevels(df3$z))+
  ylab(TeX(r'($\log(e^{p2e})$)'))+xlab('Top 25 Genes')+
  scale_x_continuous(breaks=1:25,
                     labels=g[top20])+
  theme_light()+ 
  theme(axis.text.x = element_text(size=8,angle = 90,hjust=1,
                                   vjust = 1),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        legend.title = element_blank())

df4 = data.frame('x'=log(p2e[p2e>0]))
g4 = ggplot(df4)+geom_histogram(aes(x=x),color="black", fill="gray")+
  xlab(TeX(r'($\log(e^{p2e})$ when $e^{p2e}>0$)'))+
  theme_light()+
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))

ggarrange(g4,g3,ncol=2)



#### Table 2 calculations
cc1= sapply(1:8,function(i) glist[[i]][bh_decisions[[i]]==1])
sapply(1:8,function(i) length(which(cc1[[1]]%in%cc1[[i]])))
sapply(1:8,function(i) length(which(cc1[[2]]%in%cc1[[i]])))
sapply(1:8,function(i) length(which(cc1[[3]]%in%cc1[[i]])))
sapply(1:8,function(i) length(which(cc1[[4]]%in%cc1[[i]])))
sapply(1:8,function(i) length(which(cc1[[5]]%in%cc1[[i]])))
sapply(1:8,function(i) length(which(cc1[[6]]%in%cc1[[i]])))
sapply(1:8,function(i) length(which(cc1[[8]]%in%cc1[[i]])))

sapply(1:8, function(i) length(cc_estat[[i]])/sum(bh_decisions[[i]]))
sapply(1:8, function(i) length(cc_p2e[[i]])/sum(bh_decisions[[i]]))

#### Table 3 calculations

ff1 = n$Freq[which(as.character(n$Var1)%in%hyp_alpha0.1)]
ff2 = n$Freq[which(as.character(n$Var1)%in%hyp_p2e_alpha0.1)]
100*tabulate(ff1)/length(ff1)
100*tabulate(ff2)/length(ff2)
