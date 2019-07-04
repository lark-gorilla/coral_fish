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
source('C:/coral_fish/scripts/coral_fish/functions.R')

dat<-read_xlsx('C:/coral_fish/data/Traits/Coral_DB_condensed_KMC0107.xlsx',
                sheet = 2, na='NA')

summary(dat) # 23 missing vals in max depth and reporduction
str(dat)

dat$Corallite.width.maximum<-as.numeric(dat$Corallite.width.maximum)
dat$Depth.lower<-as.numeric(dat$Depth.lower)
dat$growth_rate<-as.numeric(dat$growth_rate)

table(dat$Coloniality)
table(dat$Corallite.width.maximum)
table(dat$Depth.lower)
table(dat$'Growth.form.typical (Brig revised)')
table(dat$growth_rate)
table(dat$larval_development)
table(dat$Sexual_system)
table(dat$Wave.exposure.preference)
table(dat$Water.clarity.preference)

#make 2 ordered factors
dat$Wave.exposure.preference<-factor(dat$Wave.exposure.preference,
                          levels=c("protected", "broad","exposed"), ordered = T,
                          exclude='NA')

dat$Water.clarity.preference<-factor(dat$Water.clarity.preference,
                                     levels=c("clear", "both","turbid"), ordered = T,
                                     exclude='NA')

# apparently daisy won't accept characters? change to factor

dat$Coloniality<-factor(dat$Coloniality, exclude='NA')
dat$'Growth.form.typical (Brig revised)'<-factor(dat$'Growth.form.typical (Brig revised)', exclude='NA')
dat$larval_development<-factor(dat$larval_development, exclude='NA')
dat$Sexual_system<-factor(dat$larval_development, exclude='NA')

str(dat)

crl_out<-clVal(data=dat[,2:10], runs=1000,
               min_cl=2, max_cl=20, subs_perc=0.95,
               fast.k.h = 0.1, calc_wigl = F)


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
  geom_vline(xintercept = 4, color='cyan')

alg_out$clust_centres<-alg_out$clust_centres %>% group_by(kval, jc_match) %>%
  mutate_if(is.numeric, funs(replace(., is.na(.), mean(., na.rm=T)))) %>%
  as.data.frame()

alg_pca<-pca_vis(rundat=dat[,8:17], clValresult=alg_out$clust_centres, kval=7)

grid.arrange(alg_pca[[4]], alg_pca[[5]], alg_pca[[6]])

# write clusters out

full_crl_clust<-cutree(hclust(daisy(dat[,2:10],
                      metric='gower', stand = FALSE), method='average'), k=4)

dat$group<-full_crl_clust

write.csv(dat, 'C:/coral_fish/outputs/coral_clust.csv', quote=F, row.names=F)

