## 11/04/19 
## fish trait cluster paper analyses

#rm(list=ls())
library(cluster)
library(clue)
library(fpc)
library(ggplot2)
library(dendextend)
library(circlize)
library(vegan)
library(caret)
library(GGally)
library(dplyr)
library(reshape2)
library(ade4)
library(ggrepel)
library(adehabitatHR)
library(sf)
library(gridExtra)
library(gbm)
library(mice)

# source clVal function
source('C:/coral_fish/scripts/coral_fish/clVal.R')
source('C:/coral_fish/scripts/coral_fish/functions.R')

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

#################### run cluster validation function ##########################
###############################################################################

# analyses on Aus + Jpn combined

#trial to see difference between jaccad and distance measure
dat_aus<-dat[which(dat$AUS_sp>0),]

dat_jpn<-dat[which(dat$JPN_sp>0),]

jpn_out<-clVal(data=dat_jpn[,c("BodySize","Diet",  "Position", "Aggregation", 'DepthRange')],
               runs=100, min_cl=3,max_cl=20, subs_perc=0.95, fast.k.h = 0.2)

aus_out<-clVal(data=dat_aus[,c("BodySize","Diet",  "Position", "Aggregation", 'DepthRange')],
               runs=100, min_cl=3, max_cl=20, subs_perc=0.95, fast.k.h = 0.2)


jac_trial<-clVal(data=dat[,c("BodySize","Diet",  "Position", "Aggregation", 'DepthRange')],
                 runs=100, max_cl=20, subs_perc=0.95)


# outputs saved to memory

################## find optimal n clusters per dataset ########################
###############################################################################

# objects jac_trial (aus+jpn), aus_out and jpn_out loaded from memory

ggpairs(aus_out$stats, columns = 3:7)

#ggpairs(aus_out$stats, columns = 3:7, ggplot2::aes(colour=as.character(k)))

# remember silhouette won't necessarily correlate because it
# has a polynomial relationship with the others?

qplot(data=filter(aus_out$stats, k==6), x=sil, y=jac)

# Australia

a_melt<-melt(aus_out$stats, id.vars=c( 'k', 'runs'))
a_sum<-a_melt%>%group_by(k, variable)%>%
  summarise(mean=mean(value), median=median(value))

ggplot(data=a_melt, aes(x=k, y=value, group=k))+
  geom_boxplot()+facet_wrap(~variable, scales='free_y') # 7 apart from sil:9

ggplot()+
  geom_violin(data=a_melt, aes(x=k, y=value, group=k))+
  geom_point(data=a_sum, aes(x=k, y=mean), color='red', shape=1)+
  geom_line(data=a_sum, aes(x=k, y=mean), color='red')+
  geom_point(data=a_sum, aes(x=k, y=median), color='green', shape=1)+
  geom_line(data=a_sum, aes(x=k, y=median), color='green')+
  scale_x_continuous(breaks=3:20)+
  facet_wrap(~variable, scales='free_y')+
  geom_vline(xintercept = 7, color='cyan')
# 7 for most, 8 or 9 for silhouette
p

#ggsave('C:/coral_fish/outputs/aus_nclust2.png', width=8, height=8)

# Japan

j_melt<-melt(jpn_out$stats, id.vars=c( 'k', 'runs'))
j_sum<-j_melt%>%group_by(k, variable)%>%
  summarise(mean=mean(value), median=median(value))

p<-ggplot()+
  geom_violin(data=j_melt, aes(x=k, y=value, group=k))+
  geom_point(data=j_sum, aes(x=k, y=mean), color='red', shape=1)+
  geom_line(data=j_sum, aes(x=k, y=mean), color='red')+
  geom_point(data=j_sum, aes(x=k, y=median), color='green', shape=1)+
  geom_line(data=j_sum, aes(x=k, y=median), color='green')+
  scale_x_continuous(breaks=3:20)+
  facet_wrap(~variable, scales='free_y')+
  geom_vline(xintercept = 10, color='cyan')
# 10 for all

p
#ggsave('C:/coral_fish/outputs/jpn_nclust2.png', width=8, height=8)

# Combined Australia & Japan

c_melt<-melt(jac_trial$stats, id.vars=c( 'k', 'runs'))
c_sum<-c_melt%>%group_by(k, variable)%>%
  summarise(mean=mean(value), median=median(value))

p<-ggplot()+
  geom_violin(data=c_melt, aes(x=k, y=value, group=k))+
  geom_point(data=c_sum, aes(x=k, y=mean), color='red', shape=1)+
  geom_line(data=c_sum, aes(x=k, y=mean), color='red')+
  geom_point(data=c_sum, aes(x=k, y=median), color='green', shape=1)+
  geom_line(data=c_sum, aes(x=k, y=median), color='green')+
  scale_x_continuous(breaks=3:20)+
  facet_wrap(~variable, scales='free_y')+
  geom_vline(xintercept = 6, color='cyan')
# 6 for most, possibly 9

p
#ggsave('C:/coral_fish/outputs/combined_nclust.png', width=8, height=8)


# can now check individual cluster stability to see which is best

summary(aus_out)

rowMeans(aus_out$jaccard[[7]])
rowMeans(aus_out$jaccard[[8]])
rowMeans(aus_out$jaccard[[9]])

lapply(aus_out$jaccard, function(x){min(rowMeans(data.frame(x)))})

rowMeans(jpn_out$jaccard[[10]])
rowMeans(jpn_out$jaccard[[13]])
rowMeans(jpn_out$jaccard[[16]])

lapply(jpn_out$jaccard, function(x){min(rowMeans(data.frame(x)))})

rowMeans(jac_trial$jaccard[[6]])
rowMeans(jac_trial$jaccard[[9]])

lapply(jac_trial$jaccard, function(x){min(rowMeans(data.frame(x)))})


apply(jac_trial$jaccard[[9]], 1, var)
table(aus_out$clust_centres[aus_out$clust_centres$kval==7,]$jc_match)

################## Variable importance using BRT sensu Darling 2012 ###########
###############################################################################

eff_aus<-daisy(dat_aus[,3:9], metric='gower', stand = FALSE)

eff_jpn<-daisy(dat_jpn[,3:9], metric='gower', stand = FALSE)

eff_both<-daisy(dat[,c("BodySize","Diet",  "Position", "Aggregation", 'DepthRange')], metric='gower', stand = FALSE)

aus_FG<-cutree(hclust(eff_aus, method='average'), k=8)
jpn_FG<-cutree(hclust(eff_jpn, method='average'), k=10)
both_FG<-cutree(hclust(eff_both, method='average'), k=12)

## impute missing values using MICE (ladds et al. 2018)
dat_mice<-mice(dat[,c(3:9)], m=5, method=c(rep('norm.predict', 3), rep('polyreg', 4)))
dat_imp<-complete(dat_mice)
dat_imp<-cbind(dat[,1:2], dat_imp, dat[,10:14])

imp_aus<-dat_imp[which(dat_imp$AUS_sp>0),]
imp_jpn<-dat_imp[which(dat_imp$JPN_sp>0),]

imp_jpn$group<-factor(jpn_FG)
imp_aus$group<-factor(aus_FG)
dat_imp$group<-factor(both_FG)

brt_jpn<-gbm(group~., distribution='multinomial', n.trees=1000,
             data=imp_jpn[,c(3:9,15)]) # other parameters default
summary(brt_jpn)

brt_aus<-gbm(group~., distribution='multinomial', n.trees=1000,
             data=imp_aus[,c(3:9,15)]) # other parameters default
summary(brt_aus)

brt_both<-gbm(group~BodySize+Diet+Position+Aggregation+DepthRange, distribution='multinomial', n.trees=1000,
             data=dat_imp) # other parameters default
summary(brt_both)

## trial with party package
library(partykit)

party1<-ctree(group~BodySize+Diet+Position+Aggregation+DepthRange,
    data=imp_aus)
plot(party1)

party1<-ctree(group~Position,data=imp_aus)
plot(party1)


# Darling approach
for(i in 1:7)
{
  mydat<-imp_jpn[,c(3:9,15)]
  mydat<-mydat[,-i]
  brt_jpn<-gbm(group~., distribution='multinomial', n.trees=100,
               data=mydat)
  predBST = predict(brt_jpn,n.trees=100, newdata=imp_jpn,type='response')
  
  print(names(imp_jpn[,c(3:9,15)][i]))
 print(confusionMatrix(table(jpn_FG, apply(predBST, 1, which.max)))$overall[1])
}

adonis(eff_jpn~BodySize+Diet+Position+Aggregation+DepthRange+PLD+ParentalMode, data=imp_jpn)


for(i in 2:20)
{
both_FG<-cutree(hclust(eff_both, method='average'), k=i)
dat_imp$group<-factor(both_FG)
 brt_both<-gbm(group~BodySize+Diet+Position+Aggregation+DepthRange, distribution='multinomial', n.trees=1000,
 data=dat_imp)
print(i)
print(summary(brt_both))
}


