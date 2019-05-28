# algae trial

library(cluster)
library(clue)
library(fpc)
library(ggplot2)
library(dendextend)
library(circlize)
library(vegan)
library(caret)
library(GGally)
library(dplyr)
library(reshape2)
library(ade4)
library(ggrepel)
library(adehabitatHR)
library(sf)
library(gridExtra)

# source clVal function
source('C:/coral_fish/scripts/coral_fish/clVal.R')
source('C:/coral_fish/scripts/coral_fish/functions.R')

dat<-read.csv('C:/coral_fish/data/Traits/algae_clean08May.csv')

dat1<-dat[,c(1:7, 9, 14, 11, 16, 17, 22, 23, 24, 25, 29)]

summary(dat1) # 23 missing vals in max depth and reporduction
str(dat1)

alg_out<-clVal(data=dat1[,8:17], runs=1000, min_cl=3, max_cl=20, subs_perc=0.95)


a_melt<-melt(alg_out$stats, id.vars=c( 'k', 'runs'))
a_sum<-a_melt%>%group_by(k, variable)%>%
  summarise(mean=mean(value), median=median(value))

ggplot()+
  geom_violin(data=a_melt, aes(x=k, y=value, group=k))+
  geom_point(data=a_sum, aes(x=k, y=mean), color='red', shape=1)+
  geom_line(data=a_sum, aes(x=k, y=mean), color='red')+
  geom_point(data=a_sum, aes(x=k, y=median), color='green', shape=1)+
  geom_line(data=a_sum, aes(x=k, y=median), color='green')+
  scale_x_continuous(breaks=3:20)+
  facet_wrap(~variable, scales='free_y')+
  geom_vline(xintercept = 7, color='cyan')

alg_out$clust_centres<-alg_out$clust_centres %>% group_by(kval, jc_match) %>%
  mutate_if(is.numeric, funs(replace(., is.na(.), mean(., na.rm=T)))) %>%
  as.data.frame()

alg_pca<-pca_vis(rundat=dat1[,8:17], clValresult=alg_out$clust_centres, kval=7)

grid.arrange(alg_pca[[4]], alg_pca[[5]], alg_pca[[6]])

# write clusters out

alg_7_wig<-rowMeans(alg_out$wiggle[[7]], na.rm=T)
row.names(dat1)<-paste(dat1$Genus, dat1$Species)
full_alg_clust<-cutree(hclust(daisy(dat1[,8:17], metric='gower', stand = FALSE), method='average'), k=7)

cl_sp<-data.frame(Species=attr(full_alg_clust, 'names'), cluster=full_alg_clust, wiggliness=alg_7_wig)
cl_sp<-cl_sp[order(cl_sp$cluster, cl_sp$wiggliness),]
#write.csv(cl_sp, 'C:/coral_fish/outputs/alg_clust_wigg.csv', quote=F, row.names=F)

