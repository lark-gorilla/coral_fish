## 28/05/19 
## functional entities and trophic groupings of AUS and JPN datasets
library(ggplot2)
library(dplyr)
library(mice)

dat<-read.csv('C:/coral_fish/data/Traits/JPN_AUS_RMI_CHK_MLD_TMR_trait_master_opt2.csv', h=T)
# Running on simplist classification of Position trait

dat<-dat[which(dat$AUS_sp>0 | dat$JPN_sp>0),] # we will focus on Australia and Japan for this prelim


dat_fine<-read.csv('C:/coral_fish/data/Traits/JPN_AUS_RMI_CHK_MLD_TMR_trait_master.csv', h=T)
# more levels in aggregation and position 

dat_fine<-dat_fine[which(dat_fine$AUS_sp>0 | dat_fine$JPN_sp>0),] # we will focus on Australia and Japan for this prelim


## impute missing values using MICE (ladds et al. 2018)

dat_mice<-mice(dat[,c(3:9)], m=5, method=c(rep('norm.predict', 3), rep('polyreg', 4)))
dat2<-complete(dat_mice)
dat2<-cbind(Species=dat[,1], dat2)

# then bin continuous data
summary(dat2)

# 6 classes on cuts from Mouillot et al (2014)
dat2$BodySize <- cut(dat2$BodySize ,breaks = c(-Inf, 7,15,30,50, 80, Inf),
                              labels = c("Tiny","VSmall", "Small", "Medium", "Large", "VLarge"))

ggplot(data=dat2, aes(x=Position, y=DepthRange))+geom_jitter(width=0.05, shape=1)+geom_violin()

dat2%>%group_by(Position)%>%summarise(depth=mean(DepthRange)) 

dat2$BodySize <- cut(dat2$BodySize ,breaks = c(-Inf, 7,15,30,50, 80, Inf),
                     labels = c("Reef", "Shallow", "Ocean", "Deep"))


fishes$maxLength_bin <- cut(fishes$maxLength,breaks = c(-Inf, 15,50,100,250,500, Inf),
                            labels = c("Tiny","VSmall", "Small", "Medium", "Large", "VLarge"))


# Do we want to include sp. in the analyses? - yes.
# dat[grep('\\.', dat$Species),]


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

