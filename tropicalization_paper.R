# paper code 27/06/19

# picks up from tage where optimum number of clusters are known
#Steps
#1 calculate redundancy for each community and group
#2 split functional groups by thermal affinity to calculate 
# redundancy/vulnerability for each FG/TG group due to thermal disturbance
#3 calculate recovery diversity for each FG/TH

library(cluster)
library(dplyr)
library(ggplot2)
library(gridExtra)

#################### data clean ##########################
##########################################################

dat<-read.csv('C:/coral_fish/data/Traits/JPN_AUS_RMI_CHK_MLD_TMR_trait_master_opt2.csv', h=T)
# Running on simplist classification of Position trait

# Do we want to include sp. in the analyses? - yes.
# dat[grep('\\.', dat$Species),]

dat<-dat[which(dat$AUS_sp>0 | dat$JPN_sp>0),] # we will focus on Australia and Japan for this prelim

# remove functional duplicates
dup_trait<-paste(dat$BodySize, dat$DepthRange, dat$PLD, dat$Diet, dat$Aggregation, dat$Position, 
                 dat$ParentalMode)
dat<-dat[-which(duplicated(dup_trait)),]
dat$Species<-as.character(dat$Species)
dat[which(dat$Species=='Scarus psittacus'),]$Species<-'Scarus psittacus/spinus' # edit for one

# Including those that default to duplicates via NA after gower dist
dat<-dat[-which(dat$Species=='Caesio sp.'),]
dat<-dat[-which(dat$Species=='Ostracion immaculatus'),]

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
# Position doensn't follow a logical SINGLE order
#dat$Position<-factor(dat$Position, levels=c("SubBenthic", "Benthic","UpperBenthic",
#                                            "Demersal", "ReefPelagic","Pelagic"), ordered = T)

eff_both<-daisy(dat[,c("BodySize","Diet",  "Position", "Aggregation", 'DepthRange')], metric='gower', stand = FALSE)

both_FG<-cutree(hclust(eff_both, method='average'), k=12)

dat$FG<-both_FG

dat_aus<-dat[which(dat$AUS_sp>0),]

dat_jpn<-dat[which(dat$JPN_sp>0),]

dat_aus %>% group_by(FG, ThermalAffinity) %>% summarise(n()) %>% as.data.frame()
dat_jpn %>% group_by(FG, ThermalAffinity) %>% summarise(n()) %>% as.data.frame()

# probability that each FG has more non-tropical species
# than due to random chance

#Australia

out_aus<-NULL
for(i in 1:1000){
out_aus<-rbind(out_aus, dat_aus%>%mutate(FG2=sample(FG, replace=F))%>% 
          group_by(FG2)%>%summarize(prop_trop=length(which(ThermalAffinity=='tropical'))/n())%>%
          mutate(run=i))
}

#full
aus_full<-dat_aus%>%group_by(FG)%>%summarize(prop_trop=length(which(ThermalAffinity=='tropical'))/n())

pa<-ggplot()+geom_violin(data=out_aus, aes(x=factor(FG2), y=prop_trop))+
  geom_hline(yintercept = length(which(dat_aus$ThermalAffinity=='tropical'))/nrow(dat_aus),
             linetype='dashed', colour='dark green')+
  geom_point(data=aus_full, aes(x=factor(FG), y=prop_trop), colour='red')+
  geom_label(data=dat_aus%>%group_by(FG)%>%summarize(n=n()),
  aes(x=factor(FG), y=0.1, label=n), colour='blue')+
  xlab('Australia FGs')+ylab('Proportion of tropical sp.')

## Japan

out_jpn<-NULL
for(i in 1:1000){
  out_jpn<-rbind(out_jpn, dat_jpn%>%mutate(FG2=sample(FG, replace=F))%>% 
               group_by(FG2)%>%summarize(prop_trop=length(which(ThermalAffinity=='tropical'))/n())%>%
               mutate(run=i))
}

#full
jpn_full<-dat_jpn%>%group_by(FG)%>%summarize(prop_trop=length(which(ThermalAffinity=='tropical'))/n())

# hack to add extra group
out_jpn<-rbind(out_jpn, data.frame(FG2=c(11,12), prop_trop=c(NA, NA), run=c(1,1)))

pj<-ggplot()+geom_violin(data=out_jpn, aes(x=factor(FG2), y=prop_trop))+
  geom_hline(yintercept = length(which(dat_jpn$ThermalAffinity=='tropical'))/nrow(dat_jpn),
             linetype='dashed', colour='dark green')+
  geom_point(data=jpn_full, aes(x=factor(FG), y=prop_trop), colour='red')+
  geom_label(data=dat_jpn%>%group_by(FG)%>%summarize(n=n()),
  aes(x=as.factor(FG), y=0.1, label=n), colour='blue')+
  xlab('Japan FGs')+ylab('Proportion of tropical sp.')

grid.arrange(pa, pj, nrow=2)


