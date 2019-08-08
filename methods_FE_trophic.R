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
# DepthRange
# Is depthrange cor with Position after cut?
ggplot(data=dat2, aes(x=Position, y=DepthRange))+geom_jitter(width=0.05, shape=1)+geom_violin()

dat2%>%group_by(Position)%>%summarise(depth=mean(DepthRange)) 

ggplot(data=dat2, aes(x=DepthRange))+geom_histogram(binwidth=10)+
  scale_x_continuous(breaks=seq(20, 420, 20))

dat2$DepthRange <- cut(dat2$DepthRange ,breaks = c(-Inf, 10,20,40,100, 200, Inf),
                     labels = c("Tiny","VSmall", "Small", "Medium", "Large", "VLarge"))

table(dat2$DepthRange, dat2$Position) # not correlated

## PLD

ggplot(data=dat2, aes(x=PLD))+geom_histogram(binwidth=10)

dat2$PLD <- cut(dat2$PLD ,breaks = c(-Inf, 0, 20,40,60, 100, Inf),
                     labels = c("None", "VShort", "Short",'Medium', 'Long', "Huge"))

# ORDER necessary categorical variables
# in dat object
dat$Aggregation<-factor(dat$Aggregation, levels=c("solitary", "pairs","groups","schools"), ordered = T)
dat$Position<-factor(dat$Position, levels=c("SubBenthic", "Benthic","UpperBenthic",
                                            "Demersal", "ReefPelagic","Pelagic"), ordered = T)
# and dat2 object
dat2$Aggregation<-factor(dat2$Aggregation, levels=c("solitary", "pairs","groups","schools"), ordered = T)
dat2$Position<-factor(dat2$Position, levels=c("SubBenthic", "Benthic","UpperBenthic",
                                            "Demersal", "ReefPelagic","Pelagic"), ordered = T)
dat2$BodySize<-ordered(dat2$BodySize)
dat2$DepthRange<-ordered(dat2$DepthRange)
dat2$PLD<-ordered(dat2$PLD)

# Functional entities

#possible combinations
6*3*3*5*3*7 # 5670 Mouillot 2014
6*6*6*7*4*6*5 # 181440 This study
6*3*3*7*4*3*5 # 22680 possible refinement

dat2$FE<-paste(dat2$BodySize, dat2$DepthRange, dat2$PLD, dat2$Diet, dat2$Aggregation,
               dat2$Position, dat2$ParentalMode)

length(unique(dat2$FE)) #489
489/181440*100 #0.3% of possible combinations actually filled

FE_sum<-dat2 %>% group_by(FE) %>% summarise(n_sp=n()) %>% arrange(desc(n_sp))

ggplot(data=dat2 %>% group_by(FE) %>% summarise(n_sp=n()) %>% arrange(desc(n_sp)),
       aes(x=reorder(FE, -n_sp), y=n_sp))+
  geom_bar(stat='identity')+theme(axis.text.x = element_text(angle = 90, hjust = 1))

## Do per region not all together!!

ggplot(data=dat2, aes(x=reorder(dat2, - value)))+geom_bar(stat='count')

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





## calcing species occurence along latitudinal gradient - Japan trial

specs<-read.csv('C:/coral_fish/data/Traits/JPN_AUS_RMI_CHK_MLD_TMR_trait_master_opt2.csv', h=T)

specs<-specs[which(specs$JPN_sp>0),] 

dat<-read.csv('C:/coral_fish/data/Japan/FishData_JP_2016_final.csv', h=T)

dat$SpeciesFish<-as.character(dat$SpeciesFish)

nrow(specs);length(unique(dat$SpeciesFish))
specs$Species[which(!specs$Species %in% dat$SpeciesFish)]

dat$SpeciesFish[dat$SpeciesFish=="Apogon aureus"]<- "Ostorhinchus aureus"
dat$SpeciesFish[dat$SpeciesFish=='PLectroglyphidodon dickii']<-'Plectroglyphidodon dickii'
dat$SpeciesFish[dat$SpeciesFish=='Goniistius zebra']<-'Cheilodactylus zebra'
dat$SpeciesFish[dat$SpeciesFish=='Diagramma picta']<-'Diagramma pictum'
dat$SpeciesFish[dat$SpeciesFish=='Goniistius zonatus']<-'Cheilodactylus zonatus'
dat$SpeciesFish[dat$SpeciesFish=='Halichoeres poecilopterus']<-'Parajulis poecilepterus'
dat$SpeciesFish[dat$SpeciesFish=='Sebasticus marmoratus']<-'Sebastiscus marmoratus'
dat$SpeciesFish[dat$SpeciesFish=="Apogon doederleini"]<-'Ostorhinchus doederleini'
dat$SpeciesFish[dat$SpeciesFish=='Siganus stellatus']<-'Siganus punctatus'
dat$SpeciesFish[dat$SpeciesFish=="Apogon limenus"]<-'Ostorhinchus limenus'
dat$SpeciesFish[dat$SpeciesFish=="Chaetodon modestus"]<-'Roa modesta'


locs<-read.csv('C:/coral_fish/data/Japan/jp2015_16_waypoints.csv', h=T)

dat<-left_join(dat, locs, by=c('SiteID'= 'Site'))

# rm blanks
dat<-dat[-which(dat$SiteID==''),]
# abun to PA
dat$PA<-1

# sp in trait dataset filter?

library(labdsv)
library(ade4)

mat1<-matrify(data.frame(round(dat$lat, 2), dat$SpeciesFish, dat$PA))

table(dat[dat$Year==2015,]$SiteID)
table(dat[dat$Year==2016,]$SiteID)

bdist<-dist.binary(mat1, method=1) # jaccard dist

#shouldn't use ward centroid or median methods for jaccard dist
plot(hclust(bdist, 'single'))  #splits sites into 2 clusts from lat <31deg N

cor(cophenetic(hclust(bdist, 'average')), bdist) # best
cor(cophenetic(hclust(bdist, 'single')), bdist)
cor(cophenetic(hclust(bdist, 'complete')), bdist)

#make split and species list for each community
specs$JPN_trop<-ifelse(specs$Species%in% dat[dat$lat<31,]$SpeciesFish, 1, 0)
specs$JPN_temp<-ifelse(specs$Species%in% dat[dat$lat>31,]$SpeciesFish, 1, 0)

df4<-specs[,c(1,15, 16)]

write.csv(df4, 'C:/coral_fish/data/Japan/JPN_species_tropical_class.csv', quote=F, row.names=F)

# alternate approach n tropical and n temperate sp per site

dat<-left_join(dat, specs[,1:2], by=c('SpeciesFish'= 'Species'))
dat[dat$ThermalAffinity!='tropical',]$ThermalAffinity<-'subtropical'

d2<-dat%>%group_by(lat)%>%distinct(SpeciesFish, .keep_all=T)%>%
  summarize(n_trop=length(which(ThermalAffinity=='tropical')),
            n_subt=length(which(ThermalAffinity=='subtropical')))%>%
  as.data.frame()

row.names(d2)<-round(d2$lat, 4)
d3<-dist(d2[,2:3])
plot(hclust(d3, 'average'))

cor(cophenetic(hclust(d3, 'average')), d3) # best
cor(cophenetic(hclust(d3, 'single')), d3)
cor(cophenetic(hclust(d3, 'complete')), d3)

## calcing species occurence along latitudinal gradient - Australia data

specs<-read.csv('C:/coral_fish/data/Traits/JPN_AUS_RMI_CHK_MLD_TMR_trait_master_opt2.csv', h=T)

specs_aus<-specs[which(specs$AUS_sp>0),] 

#read in australia survey species
dat_aus<-read.csv('C:/coral_fish/data/Australia/LongTransect_Subtropical_fish_Aug2017_18Aug2018.csv')

#get list of aus_species names
aus_species<-aus_species_list$species

dat_aus$Fish<-as.character(dat_aus$Fish)
specs_aus$Species<-as.character(specs_aus$Species)

nrow(specs_aus);length(unique(dat_aus$Fish)) # hmm extra 50 sp in Maria's list?
specs_aus$Species[which(!specs_aus$Species %in% dat_aus$Fish)]
## Hmm 72 species not shared, not ~ 50.. naming issues? Hopefully there is updated sheet somewhere

# remove
dat_aus<-dat_aus[dat_aus$Site!='Siganus fuscescens',]

locs<-read.csv('C:/coral_fish/data/Australia/Australia_SitesMar2010toAug2017.csv', h=T)

dat_aus<-left_join(dat_aus, locs[,c(1, 4, 5)], by=c('Site'= 'Site.name'))

aggregate(Lat~Site, dat_aus, unique)
# sort different lats for Mudjima
dat_aus[dat_aus$Site=='Mudjimba Island',]$Lat<--26.61623

# site names that did not line up or had NA for lat, not fixing as subject to change
unique(dat_aus[is.na(dat_aus$Lat),]$Site)
paste(unique(locs$Site.name)[!unique(locs$Site.name)  %in%  unique(dat_aus$Site)])
# abun to PA
dat_aus$PA<-1

table(dat_aus[dat_aus$Year==2015,]$Site)

# sp in trait dataset filter?

library(labdsv)
library(ade4)

dat_aus%>%group_by(Site)%>%summarise(length(unique(Lat)))%>%View()
# ok a bit screwed up but proceed anyway

dat_aus_yr<-dat_aus[dat_aus$Year==2010,]
mat1<-matrify(data.frame(dat_aus_yr$Site, dat_aus_yr$Fish, dat_aus_yr$PA))

bdist<-dist.binary(mat1, method=1) # jaccard dist

#shouldn't use ward centroid or median methods for jaccard dist
plot(hclust(bdist, 'single'))  #splits sites into 3 clusts, make cut at Flat Rock 28 deg S


specs_aus$AUS_trop<-ifelse(specs_aus$Species%in% dat[dat$lat<31,]$SpeciesFish, 1, 0)
specs_aus$AUS_temp<-ifelse(specs_aus$Species%in% dat[dat$lat>31,]$SpeciesFish, 1, 0)

df4<-specs[,c(1,15, 16)]

write.csv(df4, 'C:/coral_fish/data/Japan/JPN_species_tropical_class.csv', quote=F, row.names=F)


write.csv(df4, 'C:/coral_fish/data/Japan/JPN_species_tropical_class.csv', quote=F, row.names=F)
