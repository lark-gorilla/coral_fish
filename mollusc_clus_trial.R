# Coral trial
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

dat<-read.csv('C:/coral_fish/data/Traits/mollusc_traits_condensed.csv',h=T,
              na.strings = c('NA', ''))
                

summary(dat) # 23 missing vals in max depth and reporduction
str(dat)

table(dat$Tidal.Zone)
table(dat$Trophic)
table(dat$Habitat_KMC)
table(dat$position_KMC)
table(dat$adult.mobility..Paganelli.et.al..)
table(dat$max_size_mm)
table(dat$Shell_category_KMC)
table(dat$depth_range)

dat$max_size_mm<-as.numeric(dat$max_size_mm)

# fix naming issues

dat[dat$Tidal.Zone==8 & !is.na(dat$Tidal.Zone),]$Tidal.Zone<-NA
dat[dat$Tidal.Zone=='intertidal ' & !is.na(dat$Tidal.Zone),]$Tidal.Zone<-'intertidal'
dat[dat$Tidal.Zone=='Intertidal/subtidal' & !is.na(dat$Tidal.Zone),]$Tidal.Zone<-'intertidal/subtidal'
dat$Tidal.Zone<-factor(dat$Tidal.Zone)
table(dat$Tidal.Zone)

dat[dat$Trophic=='predator ' & !is.na(dat$Trophic),]$Trophic<-'predator'
dat[dat$Trophic=='Predator' & !is.na(dat$Trophic),]$Trophic<-'predator'
dat[dat$Trophic==' grazer' & !is.na(dat$Trophic),]$Trophic<-'grazer'
dat$Trophic<-factor(dat$Trophic)
table(dat$Trophic)

dat[dat$Habitat_KMC=='rocky ' & !is.na(dat$Habitat_KMC),]$Habitat_KMC<-'rocky'
dat[dat$Habitat_KMC=='rocky/coral ' & !is.na(dat$Habitat_KMC),]$Habitat_KMC<-'rocky/coral'
dat$Habitat_KMC<-factor(dat$Habitat_KMC)
table(dat$Habitat_KMC)

dat[dat$position_KMC=='benthic ' & !is.na(dat$position_KMC),]$position_KMC<-'benthic'
dat$position_KMC<-factor(dat$position_KMC)
table(dat$position_KMC)

dat[dat$adult.mobility..Paganelli.et.al..=='crawling ' & !is.na(dat$adult.mobility..Paganelli.et.al..),]$adult.mobility..Paganelli.et.al..<-'crawling'
dat$adult.mobility..Paganelli.et.al..<-factor(dat$adult.mobility..Paganelli.et.al..)
table(dat$adult.mobility..Paganelli.et.al..)

dat[dat$Shell_category_KMC=='no shell ' & !is.na(dat$Shell_category_KMC),]$Shell_category_KMC<-'no shell'
dat[dat$Shell_category_KMC=='No shell' & !is.na(dat$Shell_category_KMC),]$Shell_category_KMC<-'no shell'
dat$Shell_category_KMC<-factor(dat$Shell_category_KMC)
table(dat$Shell_category_KMC)

#make 2 ordered factors
dat$Tidal.Zone<-factor(dat$Tidal.Zone,
                          levels=c("intertidal", "intertidal/subtidal","subtidal"), ordered = T,
                          exclude='NA')

# strip na rows

st_rws<-apply(dat[,5:12], 1, function(x){length(which(is.na(x)))})

dat1<-dat[-which(st_rws==8),]

mol_out<-clVal(data=dat1[,5:12], runs=1000,
               min_cl=2, max_cl=20, subs_perc=0.95,
               fast.k.h = 0.1, calc_wigl = F)


a_melt<-melt(mol_out$stats, id.vars=c( 'k', 'runs'))
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
  geom_vline(xintercept = 4, color='cyan')

alg_out$clust_centres<-alg_out$clust_centres %>% group_by(kval, jc_match) %>%
  mutate_if(is.numeric, funs(replace(., is.na(.), mean(., na.rm=T)))) %>%
  as.data.frame()

alg_pca<-pca_vis(rundat=dat[,8:17], clValresult=alg_out$clust_centres, kval=7)

grid.arrange(alg_pca[[4]], alg_pca[[5]], alg_pca[[6]])

# write clusters out

plot(hclust(daisy(dat1[,5:12],
                  metric='gower', stand = FALSE), method='average'))


full_crl_clust<-cutree(hclust(daisy(dat1[,5:12],
                      metric='gower', stand = FALSE), method='average'), k=2)

dat1$group_k2<-full_crl_clust

# fix to remove commas

dat1$Trophic<-gsub(',', '/', dat1$Trophic)
write.csv(dat1, 'C:/coral_fish/outputs/mollusc_clust.csv', quote=F, row.names=F)

