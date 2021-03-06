# Coral trial

library(cluster)
library(fpc)
library(ggplot2)
library(dendextend)
library(dplyr)
library(reshape2)
library(gridExtra)
library(readxl)

# source clVal function
source('C:/coral_fish/scripts/coral_fish/clVal.R')

dat<-read.csv('C:/coral_fish/data/Traits/echino_traitDB_jpn.csv')

summary(dat) # 23 missing vals in max depth and reporduction
str(dat)

dat[dat$Spines=='No ',]$Spines<-'No'
dat$Spines<-factor(dat$Spines)

table(dat$Spines)
table(dat$Max_Length)
table(dat$Depth_Range)
table(dat$Aggregation)
table(dat$Tidal_zone)
table(dat$Exposure)
table(dat$habitat_KMC)
table(dat$Diet)
table(dat$Mating_System_KMC)

#make ordered factor
dat$Exposure<-factor(dat$Exposure,
                          levels=c("Protected", "both","Exposed"), ordered = T,
                          exclude='NA')

str(dat)

table(dat$Spines)
table(dat$Max_Length)
table(dat$Depth_Range)
table(dat$Aggregation_KMC)
table(dat$Tidal_zone)
table(dat$Exposure)
table(dat$habitat_KMC)
table(dat$Diet)
table(dat$Mating_System_KMC)

crl_out<-clVal(data=dat[,c('Spines', 'Max_Length','Depth_Range',
                           'Aggregation_KMC','Tidal_zone',
                           'Exposure','habitat_KMC',
                           'Diet', 'Mating_System_KMC')],
               daisytypelist = list(symm=c(1,4,9)),
               runs=500, min_cl=2, max_cl=20, subs_perc=0.95,
               fast.k.h = 0.1, calc_wigl = F, daisyweights=rep(1, 9))
  

a_melt<-melt(crl_out$stats, id.vars=c( 'k', 'runs'))
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
  facet_wrap(~variable, scales='free_y')+
  geom_vline(xintercept = 8, color='cyan')+
  geom_vline(xintercept = 12, color='cyan')


# write clusters out

hc<-hclust(daisy(dat[,c('Spines', 'Max_Length','Depth_Range',
                        'Aggregation_KMC','Tidal_zone',
                        'Exposure','habitat_KMC',
                        'Diet', 'Mating_System_KMC')],
                  metric='gower', 
                 type = list(symm=c(1,4,9)),stand = FALSE), method='average')

plot(hc); rect.hclust(hc, k=8);rect.hclust(hc, k=12, border=4)


dat$groupk8<-cutree(hc, k=8)
dat$groupk12<-cutree(hc, k=12)


write.csv(dat, 'C:/coral_fish/outputs/echinoderm_clust_jpn.csv', quote=F, row.names=F)

