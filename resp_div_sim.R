
# simulate differences in response diversity dependent on traits
# in this example one trait will be partially correlated with
# species response. It it hypothesised that clusters resulting from
# distances built featuring this trait will show lower response diversity
rm(list=ls())
library(ggplot2)
library(dplyr)
library(vegan)
library(cluster)


dat<-read.csv('C:/coral_fish/data/Traits/JPN_AUS_RMI_CHK_MLD_TMR_trait_master_opt2.csv', h=T)
# Running on simplist classification of Position trait

dat<-dat[which(dat$AUS_sp>0 | dat$JPN_sp>0),] # we will focus on Australia and Japan for this prelim

row.names(dat)<-dat$Species

# simple winners/losers classification linked to one trait
dat$fishing_imp<-'-'
dat[dat$BodySize<50,]$fishing_imp<-'+'

# with 2 class: winners/losers, small no divided by larger
diversity(table(dat$fishing_imp), 'simpson')# use if we have > 2 classes

w_size<-dat[,c(3, 6)]
wo_size<-dat[,c(6,7)]

dist1<-daisy(w_size, metric='gower', stand = FALSE)

out<-NULL
for(i in 1:100)
{
  cutz<-cutree(hclust(dist1, method='average'), k=i)
  
  out=rbind(out, 
            dat%>%mutate(cutz=cutz)%>%group_by(cutz)%>%
              summarize(n_sp=n(),fun_div=0.5*(1-(table(fishing_imp)[which.min(table(fishing_imp))]/
                                                   table(fishing_imp)[which.max(table(fishing_imp))])), 
                        fun_div2=diversity(table(fishing_imp), 'simpson')) %>%
              mutate(k=i)
  )
}

out2<-out %>% group_by(k) %>% summarise(wm=weighted.mean(fun_div, n_sp), mn=mean(fun_div),
                                        wm2=weighted.mean(fun_div2, n_sp), mn2=mean(fun_div2))

ggplot()+
  geom_jitter(data=out, aes(x=k, y=fun_div), height=0.001, shape=1, alpha=0.5, colour='orange')+
  geom_line(data=out2, aes(x=k, y=wm), colour='red')+
  geom_line(data=out2, aes(x=k, y=mn), colour='green')+
  geom_jitter(data=out, aes(x=k, y=fun_div2), height=0.001, shape=1, alpha=0.5, colour='purple')+
  geom_line(data=out2, aes(x=k, y=wm2), colour='dark red')+
  geom_line(data=out2, aes(x=k, y=mn2), colour='dark green')
# ok so metrics roughly track


