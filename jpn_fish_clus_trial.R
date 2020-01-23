# algae trial
rm(list=ls())

library(cluster)
library(fpc)
library(ggplot2)
library(dendextend)
library(dplyr)
library(reshape2)
library(gridExtra)
library(gbm)
library(mice)
library(clue)

# source clVal function
source('C:/coral_fish/scripts/coral_fish/clVal.R')

#################### data clean ##########################
##########################################################

dat<-read.csv('C:/coral_fish/data/Traits/JPN_AUS_RMI_CHK_MLD_TMR_trait_master_opt2.csv', h=T)
# Running on simplist classification of Position trait

# Do we want to include sp. in the analyses? - yes.
# dat[grep('\\.', dat$Species),]

dat<-dat[which(dat$JPN_sp>0),] # we will focus on JAPAN 

# remove functional duplicates
#dup_trait<-paste(dat$BodySize, dat$DepthRange, dat$PLD, dat$Diet, dat$Aggregation, dat$Position, 
#                 dat$ParentalMode)
#dat<-dat[-which(duplicated(dup_trait)),]
#dat$Species<-as.character(dat$Species)
#dat[which(dat$Species=='Scarus psittacus'),]$Species<-'Scarus psittacus/spinus' # edit for one
# Including those that default to duplicates via NA after gower dist
#dat<-dat[-which(dat$Species=='Caesio sp.'),]
#dat<-dat[-which(dat$Species=='Ostracion immaculatus'),]

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

### check cor between distance and hclust copophenetic and best algorithum ####
###############################################################################

# regular gower dist
distreg<-daisy(dat[,c("BodySize","Diet",  "Position", "Aggregation", 'DepthRange',
                      'PLD', 'ParentalMode')],
               metric='gower', stand = FALSE)

distlog<-daisy(dat[,c("BodySize","Diet",  "Position", "Aggregation", 'DepthRange',
                      'PLD', 'ParentalMode')],
               metric='gower', stand = FALSE, type = list(logratio = c(1,5, 6)))

plot(hclust(distreg, 'average'), labels=F);rect.hclust(hclust(distreg, 'average'),k=14)
plot(hclust(distlog, 'average'), labels=F);rect.hclust(hclust(distlog, 'average'),k=14)

hclust_methods <-c("ward.D", "ward.D2", "single",
                   "complete", "average", "mcquitty") # median and centroid dont work

hclustreg <- lapply(hclust_methods, function(m) hclust(distreg, m))
names(hclustreg) <- hclust_methods

## see how well the trees represent the original distance matrix.
## method=spectral is the 2-norm method proposed by Merigot et al. (2010)
alg_comp_reg<-cl_dissimilarity(hclustreg, distreg, method = "spectral")
alg_comp_reg

## get threshold 2-norm value for dendrgram to accurately represent dist. 
## code from within Mouchet et al. 2008 function: ("http://villeger.sebastien.free.fr/R%20scripts/GFD_matcomm.R")

thresh<-lapply(hclustreg, function(x){2*sqrt((var(distreg)+
                                                var(cl_ultrametric(x))))*
    sqrt(dim(dat)[1])})

# put together
data.frame(method=hclust_methods, two.norm=alg_comp_reg[,1],
           thresh=unlist(thresh), diff=alg_comp_reg[,1]-unlist(thresh))

# Nothing really in it but log marginally better

#################### run cluster validation function ##########################
###############################################################################

# analyses on Aus + Jpn combined


clus_out<-clVal(data=dat[,c("BodySize","Diet",  "Position", "Aggregation", 'DepthRange',
                            'PLD', 'ParentalMode')],
                runs=100, min_cl=3, max_cl=20, subs_perc=0.95,
                fast.k.h = 0.2, calc_wigl = F, logvars = F)


################## find optimal n clusters ########################
###############################################################################

# regular daisy
a_melt<-melt(clus_out$stats[c(1:4,6)], id.vars=c( 'k', 'runs'))
a_sum<-a_melt%>%group_by(k, variable)%>%
  summarise(mean=mean(value), median=median(value))

ggplot(data=a_melt, aes(x=k, y=value, group=k))+
  geom_boxplot()+facet_wrap(~variable, scales='free_y') # 9

ggplot()+
  geom_violin(data=a_melt, aes(x=k, y=value, group=k))+
  geom_point(data=a_sum, aes(x=k, y=mean), color='red', shape=1)+
  geom_line(data=a_sum, aes(x=k, y=mean), color='red')+
  geom_point(data=a_sum, aes(x=k, y=median), color='green', shape=1)+
  geom_line(data=a_sum, aes(x=k, y=median), color='green')+
  scale_x_continuous(breaks=3:20)+
  facet_wrap(~variable, scales='free_y')+
  geom_vline(xintercept = 14, color='cyan')


# can now check individual cluster stability to see which is best
# see Details in ?clusterboot and Henning (2008) Journal of Multivariate Analysis
# Jaccard guidelines
# <0.5 = 'dissolved cluster'
# 0.6-0.75 = 'patterns, but clusters highly doubtful'
# 0.75-0.85 = 'valid, stable cluster'
# >0.85 = 'highly stable'

# regular daisy

lapply(clus_out$jaccard, function(x){min(rowMeans(data.frame(x)))})

#mean jaccard similarity per cluster
rowMeans(clus_out$jaccard[[9]])
rowMeans(clus_out$jaccard[[14]])


# variance in jaccard similarity per cluster
apply(clus_out$jaccard[[9]], 1, var)
table(clus_out$clust_centres[clus_out$clust_centres$kval==9,]$jc_match)

################## Variable importance using BRT sensu Darling 2012 ###########
###############################################################################

clust14_reg<-cutree(hclust(distreg, method='average'), k=14)
clust14_log<-cutree(hclust(distlog, method='average'), k=14)


## impute missing values using MICE (ladds et al. 2018)
dat_mice<-mice(dat[,c(3:9)], m=5, method=c(rep('norm.predict', 3), rep('polyreg', 4)))
dat_imp<-complete(dat_mice)
dat_imp<-cbind(dat[,1:2], dat_imp, dat[,10:14])

dat_imp$groupreg<-factor(clust14_reg)
dat_imp$grouplog<-factor(clust14_log)

table(dat_imp$groupreg, dat_imp$grouplog)

brt_reg<-gbm(groupreg~BodySize+Diet+Position+Aggregation+DepthRange+PLD+ParentalMode,
             distribution='multinomial', n.trees=1000, data=dat_imp) # other parameters default
summary(brt_reg)
brt_log<-gbm(grouplog~BodySize+Diet+Position+Aggregation+DepthRange+PLD+ParentalMode,
             distribution='multinomial', n.trees=1000, data=dat_imp) # other parameters default
summary(brt_log)


# Write out data with n clusters attached
dat$group<-cutree(hclust(distreg, method='average'), k=14)
dat$group_log<-cutree(hclust(distlog, method='average'), k=14)
write.csv(dat, 'C:/coral_fish/outputs/jnp_fish_clust_k14.csv', quote=F, row.names=F)
