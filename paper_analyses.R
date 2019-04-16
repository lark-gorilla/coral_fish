## 11/04/19 
## fish trait cluster paper analyses

rm(list=ls())
library(cluster)
library(clue)
library(fpc)
library(ggplot2)
library(dendextend)
library(circlize)
library(vegan)
library(caret)

# source clVal function
source('C:/coral_fish/scripts/coral_fish/clVal.R')


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
dat<-dat[-which(dat$Species=='Choerodon cyanodus'),]
dat<-dat[-which(dat$Species=='Ostracion immaculatus'),]


row.names(dat)<-dat$Species


## Edit to some trait values from MB 25/10/18
#dat[dat$Species=='Brotula multibarbata',]$DepthRange<-219
dat[dat$Species=='Mobula birostris',]$BodySize<-450
dat[dat$Species=='Amphiprion sandaracinos',]$BodySize<-14

## set ceiling for numeric variables for scaling purposes 04/03/19
dat[dat$PLD>=100 & !is.na(dat$PLD),]$PLD<-100
dat[dat$DepthRange>=200 & !is.na(dat$DepthRange),]$DepthRange<-200

#################### run cluster validation function ##########################
###############################################################################

# analyses on Aus + Jpn combined

#trial to see difference between jaccad and distance measure
jac_trial<-clVal(data=dat[,3:9], runs=50, max_cl=20, cluster_match = 'jaccard')
dist_trial<-clVal(data=dat[,3:9], runs=50, max_cl=20, cluster_match = 'distance')

jac_trial$stats$method='jac'
dist_trial$stats$method='dist'

outy<-rbind(jac_trial$stats, dist_trial$stats)

library(GGally)
ggpairs(outy, columns = 3:7, ggplot2::aes(colour=method))


# see how many clusters in each regional community
dat_aus<-dat[which(dat$AUS_sp>0),]

dat_jpn<-dat[which(dat$JPN_sp>0),]

aus_out<-clVal(data=dat_aus[,3:9], runs=50, max_cl=20, cluster_match = 'jaccard')
jpn_out<-clVal(data=dat_jpn[,3:9], runs=50, max_cl=20, cluster_match = 'jaccard')

aus_out$stats$region='AUS'
jpn_out$stats$region='JPN'

library(reshape2)

aj_melt<-melt(rbind(aus_out$stats, jpn_out$stats), id.vars=c('region', 'k', 'runs'))

ggplot(data=aj_melt[aj_melt$region=='AUS',], aes(x=k, y=value, group=k))+
      geom_boxplot()+facet_grid(variable~., scales='free_y')

ggplot(data=aj_melt[aj_melt$region=='AUS' & aj_melt$variable=='kap',], aes(x=k, y=value, group=k))+
  geom_boxplot()# 7
ggplot(data=aj_melt[aj_melt$region=='AUS' & aj_melt$variable=='acc',], aes(x=k, y=value, group=k))+
  geom_boxplot()# 7
ggplot(data=aj_melt[aj_melt$region=='AUS' & aj_melt$variable=='sil',], aes(x=k, y=value, group=k))+
  geom_boxplot()# 9, 14
ggplot(data=aj_melt[aj_melt$region=='AUS' & aj_melt$variable=='jac',], aes(x=k, y=value, group=k))+
  geom_boxplot()# 7
ggplot(data=aj_melt[aj_melt$region=='AUS' & aj_melt$variable=='wig',], aes(x=k, y=value, group=k))+
  geom_boxplot()# 5, 7
## rand index ?

ggplot(data=aj_melt[aj_melt$region=='JPN' & aj_melt$variable=='kap',], aes(x=k, y=value, group=k))+
  geom_boxplot()# 10
ggplot(data=aj_melt[aj_melt$region=='JPN' & aj_melt$variable=='acc',], aes(x=k, y=value, group=k))+
  geom_boxplot()# 10
ggplot(data=aj_melt[aj_melt$region=='JPN' & aj_melt$variable=='sil',], aes(x=k, y=value, group=k))+
  geom_boxplot()# 10, 13
ggplot(data=aj_melt[aj_melt$region=='JPN' & aj_melt$variable=='jac',], aes(x=k, y=value, group=k))+
  geom_boxplot()# 10
ggplot(data=aj_melt[aj_melt$region=='JPN' & aj_melt$variable=='wig',], aes(x=k, y=value, group=k))+
  geom_boxplot()# 10, 13

# can now check individual cluster stability to see which is best
rowMeans(aus_out$jaccard[[7]])
rowMeans(aus_out$jaccard[[9]])
rowMeans(aus_out$jaccard[[14]])

summary(rowMeans(aus_out$jaccard[[7]])) #best
summary(rowMeans(aus_out$jaccard[[9]]))
summary(rowMeans(aus_out$jaccard[[14]]))

rowMeans(jpn_out$jaccard[[10]])
rowMeans(jpn_out$jaccard[[13]])


summary(rowMeans(jpn_out$jaccard[[10]])) # best (lowest), but pretty much even
summary(rowMeans(jpn_out$jaccard[[13]]))

