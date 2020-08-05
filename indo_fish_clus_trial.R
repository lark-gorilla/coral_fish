# Maarten's Indo fish clustering trial

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

dat<-read_xlsx('C:/coral_fish/data/Traits/Maarten_Wallacea_Fish_Traits.xlsx')

summary(dat) # 23 missing vals in max depth and reporduction
str(dat)
names(dat)[5]<-'PLD'
dat$MaxLength<-as.numeric(dat$MaxLength)
dat$DRange<-as.numeric(dat$DRange)
dat$PLD<-as.numeric(dat$PLD)

table(dat$Trophic)
table(dat$Aggregation)
dat[dat$Aggregation=='Groups',]$Aggregation<-'groups'
dat[dat$Aggregation=='Harems',]$Aggregation<-'harems'
dat[dat$Aggregation=='Solitary',]$Aggregation<-'solitary'
table(dat$Position)
table(dat$ParentalMode)
dat[dat$ParentalMode=='Brooders',]$ParentalMode<-'brooders'
table(dat$SpawnMode)
dat[dat$SpawnMode=='pairedSpawn'& !is.na(dat$SpawnMode),]$SpawnMode<-'PairedSpawn'
dat[dat$SpawnMode=='Pairedspawn'& !is.na(dat$SpawnMode),]$SpawnMode<-'PairedSpawn'

#make factor
dat$Trophic<-factor(dat$Trophic, exclude='NA')
dat$Aggregation<-factor(dat$Aggregation, exclude='NA')
dat$Position<-factor(dat$Position, exclude='NA')
dat$ParentalMode<-factor(dat$ParentalMode, exclude='NA')
dat$SpawnMode<-factor(dat$SpawnMode, exclude='NA')

str(dat)

crl_out<-clVal(data=dat[,c('MaxLength', 'DRange','PLD',
                           'Trophic','Aggregation',
                           'Position','ParentalMode',
                           'SpawnMode')],
               daisytypelist = list(logratio = c(1,2,3)),
               runs=100, min_cl=2, max_cl=70, subs_perc=0.95,
               fast.k.h = 0.1, calc_wigl = F, daisyweights=rep(1, 8))
  

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
  scale_x_continuous(breaks=2:70)+
  facet_wrap(~variable, scales='free_y')+
  geom_vline(xintercept = 20, color='cyan')+
  geom_vline(xintercept = 42, color='cyan')+
  geom_vline(xintercept = 56, color='cyan')+
  geom_vline(xintercept = 70, color='cyan')


# write clusters out

hc<-hclust(daisy(dat[,c('MaxLength', 'DRange','PLD',
                        'Trophic','Aggregation',
                        'Position','ParentalMode',
                        'SpawnMode')],
                  metric='gower', 
                 type = list(logratio = c(1,2,3)),stand = FALSE), method='average')

plot(hc); rect.hclust(hc, k=20);rect.hclust(hc, k=42, border=4)
rect.hclust(hc, k=56, border=5);rect.hclust(hc, k=70, border=6)

dat$FG_k20<-cutree(hc, k=20)
dat$FG_k42<-cutree(hc, k=42)
dat$FG_k56<-cutree(hc, k=56)
dat$FG_k70<-cutree(hc, k=70)

write.csv(dat, 'C:/coral_fish/outputs/Wallacea_Fish_Traits_clustered_FINE.csv', quote=F, row.names=F)

