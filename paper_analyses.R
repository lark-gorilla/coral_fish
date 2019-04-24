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
jac_trial<-clVal(data=dat[,3:9], runs=1000, max_cl=20, subs_perc=0.95)

dat_aus<-dat[which(dat$AUS_sp>0),]

dat_jpn<-dat[which(dat$JPN_sp>0),]

aus_out<-clVal(data=dat_aus[,3:9], runs=1000, max_cl=20, subs_perc=0.95)
jpn_out<-clVal(data=dat_jpn[,3:9], runs=1000, max_cl=20, subs_perc=0.95)

# outputs saved to memory

################## find optimal n clusters per dataset ########################
###############################################################################

# objects jac_trial (aus+jpn), aus_out and jpn_out loaded from memory
jac_trial$stats$dis<-NULL
aus_out$stats$dis<-NULL
jpn_out$stats$dis<-NULL

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

#ggsave('C:/coral_fish/outputs/aus_nclust.png', width=12, height=8)

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
#ggsave('C:/coral_fish/outputs/jpn_nclust.png', width=12, height=8)

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
#ggsave('C:/coral_fish/outputs/combined_nclust.png', width=12, height=8)


# can now check individual cluster stability to see which is best

summary(aus_out)

rowMeans(aus_out$jaccard[[7]])
rowMeans(aus_out$jaccard[[8]])
rowMeans(aus_out$jaccard[[9]])

rowMeans(jpn_out$jaccard[[10]])
rowMeans(jpn_out$jaccard[[13]])

rowMeans(jac_trial$jaccard[[6]])
rowMeans(jac_trial$jaccard[[9]])

apply(jac_trial$jaccard[[9]], 1, var)

# Do PCA of optimal clustering solution

aus7<-aus_out$clust_centres[aus_out$clust_centres$kval==7,]

# setup weights for each column, factors penalised for n levels
my_wgt<-c(1,1,1,rep(1/7, 7), rep(1/4, 4), rep(1/6, 6), rep(1/5, 5))

xpand <- function(d) do.call("expand.grid", rep(list(1:nrow(d)), 2))
euc_norm <- function(x) sqrt(sum(x^2))
euc_dist <- function(mat, weights=1) {
  iter <- xpand(mat)
  vec <- mapply(function(i,j) euc_norm(weights*(mat[i,] - mat[j,])), 
                iter[,1], iter[,2])
  matrix(vec,nrow(mat), nrow(mat))
}

wgt_dist<-euc_dist(mat=aus7[,4:28], weights=my_wgt)

scl_wt<-corpcor::wt.scale(aus7[,4:28], w=my_wgt, center=T, scale=T)


#maybe need to stand and centre before distance calc to be sure..? (see vegan)

d1<-daisy(aus7[,4:28], metric='euclidean', stand=T, weights=rep(1, 25))


d1<-scalewt(aus7[,4:28])
