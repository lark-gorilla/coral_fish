
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
allv<-dat[,c(3:9)]

dist1<-daisy(w_size, metric='gower', stand = FALSE)
dist2<-daisy(wo_size, metric='gower', stand = FALSE)
dist3<-daisy(allv, metric='gower', stand = FALSE)

out<-NULL
for(i in 1:100)
{
  cutz_w<-cutree(hclust(dist1, method='average'), k=i)
  cutz_wo<-cutree(hclust(dist2, method='average'), k=i)
  cutz_all<-cutree(hclust(dist3, method='average'), k=i)
    
  out=rbind(out, cbind( 
            dat%>%mutate(cutz_w=cutz_w)%>%group_by(cutz_w)%>%
              summarize(n_sp_w=n(), 
                        fun_div2_w=diversity(table(fishing_imp), 'simpson')),
            dat%>%mutate(cutz_all=cutz_all)%>%group_by(cutz_all)%>%
              summarize(n_sp_a=n(), 
                        fun_div2_a=diversity(table(fishing_imp), 'simpson')),
              dat%>%mutate(cutz_wo=cutz_wo)%>%group_by(cutz_wo)%>%
              summarize(n_sp_wo=n(), 
                        fun_div2_wo=diversity(table(fishing_imp), 'simpson'))%>%
              mutate(k=i))
  )
  print(i)
}

out2<-out %>% group_by(k) %>% summarise(wm_w=weighted.mean(fun_div2_w, n_sp_w), mn_w=mean(fun_div2_w),
                                        wm_wo=weighted.mean(fun_div2_wo, n_sp_wo), mn_wo=mean(fun_div2_wo),
                                        wm_a=weighted.mean(fun_div2_a, n_sp_a), mn_a=mean(fun_div2_a))

ggplot()+
  geom_jitter(data=out, aes(x=k, y=fun_div_w), height=0.001, shape=1, alpha=0.5, colour='red')+
  geom_line(data=out2, aes(x=k, y=wm_w), colour='dark red')+

  geom_jitter(data=out, aes(x=k, y=fun_div_wo), height=0.001, shape=1, alpha=0.5, colour='green')+
  geom_line(data=out2, aes(x=k, y=wm_wo), colour='dark green')+
  
  geom_jitter(data=out, aes(x=k, y=fun_div_a), height=0.001, shape=1, alpha=0.5, colour='purple')+
  geom_line(data=out2, aes(x=k, y=wm_a), colour='purple')

ggplot()+
  geom_line(data=out2, aes(x=k, y=wm_w), colour='dark red')+
  geom_line(data=out2, aes(x=k, y=wm_wo), colour='dark green')+
  geom_line(data=out2, aes(x=k, y=wm_a), colour='purple')

# simulate null model of species randomly assigned to
# k groups and resp div calculated


sims<-NULL
for(i in 1:5)
{
for(j in 1:694)
  {
  # takes dendrogram allocation of clusters at cut k and randomises
  # the species into these groups
  cutz_w<-cutree(hclust(dist1, method='average'), k=j)
  
    sims<-rbind(sims,
    dat%>%mutate(cutz=sample(cutz_w, replace=F))%>%
      group_by(cutz)%>%summarize(n_sp=n(),
                fun_div2=diversity(table(fishing_imp), 'simpson'))%>%
      mutate(k=j, run=i)
    )
  }
print(i)
}

sims2<-sims %>% group_by(k, run) %>% summarise(wm=weighted.mean(fun_div2, n_sp), mn=mean(fun_div2))

# geom_jitter makes it slow                                        
ggplot()+
  geom_jitter(data=sims, aes(x=k, y=fun_div2), height=0.001, shape=1, alpha=0.5)+
  geom_line(data=sims2, aes(x=k, y=wm), colour='red')+
  geom_line(data=sims2, aes(x=k, y=mn), colour='green')

# Trial to plot null simulation and data together for
# with-fishing-pressure data (dist1), using Simpson index

out<-NULL
for(i in 1:694)
{
  cutz_w<-cutree(hclust(dist1, method='average'), k=i)

  out=rbind(out,
    dat%>%mutate(cutz_w=cutz_w)%>%group_by(cutz_w)%>%
      summarize(n_sp_w=n(), fun_div2_w=diversity(table(fishing_imp), 'simpson'))%>%
      mutate(k=i))
  print(i)
}

out2<-out %>% group_by(k) %>% summarise(wm_w=weighted.mean(fun_div2_w, n_sp_w))

ggplot()+
  geom_line(data=sims2, aes(x=k, y=wm, colour=factor(run)))+
  geom_line(data=out2, aes(x=k, y=wm_w), colour='red')
  
# Ok so conceptually it works..

## Trial with Japan data

# sort trait data
row.names(dat)<-dat$Species
## Edit to some trait values from MB 25/10/18
#dat[dat$Species=='Brotula multibarbata',]$DepthRange<-219
dat[dat$Species=='Mobula birostris',]$BodySize<-450
dat[dat$Species=='Amphiprion sandaracinos',]$BodySize<-14
## set ceiling for numeric variables for scaling purposes 04/03/19
dat[dat$PLD>=100 & !is.na(dat$PLD),]$PLD<-100
dat[dat$DepthRange>=200 & !is.na(dat$DepthRange),]$DepthRange<-200
# ORDER necessary categorical variables
dat$Aggregation<-factor(dat$Aggregation, levels=c("solitary", "pairs","groups","schools"), ordered = T)
# 

trop_class<-read.csv('C:/coral_fish/data/Japan/JPN_species_tropical_class.csv', h=T)
trop_class$variable<-gsub('\\.', ' ', trop_class$variable)

dat_jpn<-dat[which(dat$JPN_sp>0),]

dat_jpn<-left_join(dat_jpn, trop_class, by=c('Species'='variable'))

# Informal check to see how latitudinal classification lines up with lit.
table(dat_jpn$ThermalAffinity, dat_jpn$class) # not too bad, some lit classed tropical sp to to sub-tropical

# quick trial 
library(gbm)
g1<-gbm(class~., data=na.omit(dat_jpn[,c(3:9,16)]) )
summary(g1)


# with 2 class: winners/losers, small no divided by larger
diversity(table(dat_jpn$class), 'simpson')

t_res<-dat_jpn[,c("BodySize","DepthRange","PLD","ParentalMode")]
t_eff<-dat_jpn[,c("ThermalAffinity", "BodySize","Diet",  "Position", "Aggregation")]
t_all<-dat_jpn[,c(3:9)]

dist1<-daisy(t_res, metric='gower', stand = FALSE)
dist2<-daisy(t_eff, metric='gower', stand = FALSE)
dist3<-daisy(t_all, metric='gower', stand = FALSE)

out<-NULL
for(i in 1:nrow(dat_jpn))
{
  cutz_res<-cutree(hclust(dist1, method='average'), k=i)
  cutz_eff<-cutree(hclust(dist2, method='average'), k=i)
  cutz_all<-cutree(hclust(dist3, method='average'), k=i)
  
  out=rbind(out, cbind( 
    dat_jpn%>%mutate(cutz_res=cutz_res)%>%group_by(cutz_res)%>%
      summarize(n_sp_res=n(), 
                div_res=diversity(table(class), 'simpson')),
    dat_jpn%>%mutate(cutz_eff=cutz_eff)%>%group_by(cutz_eff)%>%
      summarize(n_sp_eff=n(), 
                div_eff=diversity(table(class), 'simpson')),
    dat_jpn%>%mutate(cutz_all=cutz_all)%>%group_by(cutz_all)%>%
      summarize(n_sp_all=n(), 
                div_all=diversity(table(class), 'simpson'))%>%
      mutate(k=i))
  )
  print(i)
}

out2<-out %>% group_by(k) %>% summarise(wm_res=weighted.mean(div_res, n_sp_res), mn_res=mean(div_res),
                                        wm_eff=weighted.mean(div_eff, n_sp_eff), mn_eff=mean(div_eff),
                                        wm_all=weighted.mean(div_all, n_sp_all), mn_all=mean(div_all))

ggplot()+
  geom_line(data=out2, aes(x=k, y=wm_res), colour='dark blue')+
  geom_line(data=out2, aes(x=k, y=wm_eff), colour='dark green')+
  geom_line(data=out2, aes(x=k, y=wm_all), colour='dark red')+
  geom_line(data=sims2, aes(x=k, y=wm, colour=factor(run)))+
  xlim(c(0, 20))+ylim(c(0.4, 0.6))




sims<-NULL
for(i in 1:5)
{
  for(j in 1:nrow(dat_jpn))
  {
    cutz_w<-cutree(hclust(dist1, method='average'), k=j)
    
    sims<-rbind(sims,
                dat_jpn%>%mutate(cutz=sample(cutz_w, replace=F))%>%
                  group_by(cutz)%>%summarize(n_sp=n(),
                                             fun_div2=diversity(table(class), 'simpson'))%>%
                  mutate(k=j, run=i)
    )
  }
  print(i)
}

sims2<-sims %>% group_by(k, run) %>% summarise(wm=weighted.mean(fun_div2, n_sp), mn=mean(fun_div2))



  