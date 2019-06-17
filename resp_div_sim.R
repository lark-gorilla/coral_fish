
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
  