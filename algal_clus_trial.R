# algae trial

library(cluster)
library(fpc)
library(ggplot2)
library(dendextend)
library(dplyr)
library(reshape2)
library(gridExtra)

# source clVal function
source('C:/coral_fish/scripts/coral_fish/clVal.R')
source('C:/coral_fish/scripts/coral_fish/functions.R')

dat<-read.csv('C:/coral_fish/data/Traits/algae_db_0407_clean.csv')

# selects main columns and clips out NAs at end of sheet
dat1<-dat[1:88,c(1:7,  11, 18, 13,21, 26, 27, 28, 29, 33)]

summary(dat1) # 23 missing vals in max depth and reporduction
str(dat1)

table(dat1$Structure)
table(dat1$Holdfast.morphology)
table(dat1$Thallus.Height..cm....max)
table(dat1$Support.mechanism.KMC)
table(dat1$Min.depth..m.)
table(dat1$Max.depth..m.)
table(dat1$Substrate)
table(dat1$Tidal.Zone)
table(dat1$Reproduction)

# make DepthRange variable
dat1$Depthrange<-dat1$Max.depth..m.-dat1$Min.depth..m.

table(dat1$Depthrange)

# edit large Thallus Height from 1500 to 40
dat1[dat1$Thallus.Height..cm....max==1500 &
       !is.na(dat1$Thallus.Height..cm....max),]$Thallus.Height..cm....max<- 40

# make ordered variable
dat1$Tidal.Zone<-factor(dat1$Tidal.Zone,
                        levels=c("Intertidal", "Intertidal/ subtidal","Subtidal"),
                        ordered = T)

#remove NA row
dat1<-filter(dat1, Genus!='Sargassum sp.')

alg_out<-clVal(data=dat1[,c(8:11, 14:17)],
               runs=1000, min_cl=2, max_cl=20, subs_perc=0.95,
               calc_wigl = F)


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

