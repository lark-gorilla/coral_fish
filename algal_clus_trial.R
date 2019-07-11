# algae trial
rm(list=ls())

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

dat<-read.csv('C:/coral_fish/data/Traits/algae_db_1007_clean.csv')
row.names(dat)<-dat$species

dat1<-dat
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

dat1<-dat1[c('Structure', 'Holdfast.morphology', 'Thallus.Height..cm....max',
              'Support.mechanism.KMC', 'Depthrange', 'Substrate', 'Tidal.Zone',
              'Reproduction')]
# check NAs per row

table(apply(dat1, 1, function(x){length(which(is.na(x)==T))}))

which(apply(dat1, 1, function(x){length(which(is.na(x)==T))})>4) # 8,71, 90, 98, 111
#remove NA row
dat1<-dat1[-c(8,71, 90, 98, 111),]

# remove 'Dictyota adnata Zanardini'
dat1<-dat1[-3,]

alg_out<-clVal(data=dat1,
               runs=1000, min_cl=2, max_cl=20, subs_perc=0.95,
               calc_wigl = F)


a_melt<-melt(alg_out$stats, id.vars=c( 'k', 'runs'))
a_sum<-a_melt%>%group_by(k, variable)%>%
  summarise(mean=mean(value), median=median(value))

a_melt<-filter(a_melt, variable!='wig')
a_sum<-filter(a_sum, variable!='wig')

ggplot()+
  geom_violin(data=a_melt, aes(x=k, y=value, group=k))+
  geom_point(data=a_sum, aes(x=k, y=mean), color='red', shape=1)+
  geom_line(data=a_sum, aes(x=k, y=mean), color='red')+
  geom_point(data=a_sum, aes(x=k, y=median), color='green', shape=1)+
  geom_line(data=a_sum, aes(x=k, y=median), color='green')+
  scale_x_continuous(breaks=2:20)+
  facet_grid(~variable~., scales='free_y')+
  geom_vline(xintercept = 7, color='cyan')

alg_out$clust_centres<-alg_out$clust_centres %>% group_by(kval, jc_match) %>%
  mutate_if(is.numeric, funs(replace(., is.na(.), mean(., na.rm=T)))) %>%
  as.data.frame()

alg_pca<-pca_vis(rundat=dat1[,8:17], clValresult=alg_out$clust_centres, kval=7)

grid.arrange(alg_pca[[4]], alg_pca[[5]], alg_pca[[6]])

# write clusters out

# do 2 and 4 k solution

full_alg_clust<-cutree(hclust(daisy(dat1, metric='gower', stand = FALSE),
                              method='average'), k=5)

plot(hclust(daisy(dat1, metric='gower', stand = FALSE),
            method='average'))

dat1$group<-full_alg_clust
dat1$Species<-row.names(dat1)
write.csv(dat1, 'C:/coral_fish/outputs/alg_clust_k5.csv', quote=F, row.names=F)

# informal validation

library(gbm)
library(mice)

dat_mice<-mice(dat1, m=5, method=c('polyreg', 'polyreg','norm.predict',
                                                    rep('polyreg', 4),'norm.predict'))

dat_imp<-complete(dat_mice)
dat_imp<-cbind(dat_imp, group=full_alg_clust)

dat_imp$group<-factor(dat_imp$group)

brt_alg<-gbm(group~., distribution='multinomial', n.trees=1000,
             data=dat_imp) # other parameters default
summary(brt_alg)

