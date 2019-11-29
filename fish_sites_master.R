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

mat_jpn<-matrify(data.frame(dat$Name.x, dat$SpeciesFish, dat$PA))

table(dat[dat$Year==2015,]$SiteID)
table(dat[dat$Year==2016,]$SiteID)

bdist<-dist.binary(mat_jpn, method=1) # jaccard dist

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
dat_aus<-read.csv('C:/coral_fish/data/Australia/LongTransect_Subtropical_fish_Sep2018.csv')

#get list of aus_species names

dat_aus$Fish<-as.character(dat_aus$Fish)
specs_aus$Species<-as.character(specs_aus$Species)

nrow(specs_aus);length(unique(dat_aus$Fish)) # OK only 2 different
specs_aus$Species[which(!specs_aus$Species %in% dat_aus$Fish)]

#remove species that occur in sampling data but not species/traits list
badfish<-unique(dat_aus$Fish[which(!dat_aus$Fish %in% specs_aus$Species)])
dat_aus<-dat_aus[-which(dat_aus$Fish%in%badfish),]
dat_aus$Fish<-factor(dat_aus$Fish)


# remove and fix
dat_aus$Site<-as.character(dat_aus$Site)
dat_aus<-dat_aus[dat_aus$Site!='',]

locs<-read.csv('C:/coral_fish/data/Australia/Australia_SitesMar2010toSep2019.csv', h=T)

# edit locs Site.name so it lines up
locs$Site.name<-as.character(locs$Site.name)
locs[locs$Site.name=='Pialba',]$Site.name<-"Pialba Shallow"
locs[locs$Site.name=='Big Woody',]$Site.name<-"Big Woody Shallow"
locs[locs$Site.name=='Gatakers hig div site',]$Site.name<-"Gataker High Diversity Site"
locs[locs$Site.name=='Woolgoolga Headland Reef',]$Site.name<-"Woolgoolga Headland"
locs[locs$Site.name=='Inner Gneering Shoals',]$Site.name<-"Inner Gneerings"
locs[locs$Site.name=='Hendersons Shoals',]$Site.name<-"Hendersons Rock"
locs[locs$Site.name=='Latitude Rock, Forster',]$Site.name<-"Latitude Rock"
locs[locs$Site.name=='Julian Rocks False Trench',]$Site.name<-"Julian Rock False Trench"
locs[locs$Site.name=='Julian Rocks Nursery',]$Site.name<-"Julian Rock Nursery"
locs[locs$Site.name=='Lady Elliot',]$Site.name<-"Lady Elliot Island"
locs[locs$Site.name=='Lady Musgrave',]$Site.name<-"Lady Musgrave Island"
locs[locs$Site.name=='Heron - Tenemants',]$Site.name<-"Tenemants Buoy"
locs[locs$Site.name=='Heron - Turbistari',]$Site.name<-"Turbistari"
locs[locs$Site.name=="Heron - Libby's Lair",]$Site.name<-"Libbys Lair"

dat_aus<-left_join(dat_aus, locs[,c(1, 6,7)], by=c('Site'= 'Site.name'))

aggregate(Lat~Site, dat_aus, unique)
# left_join puts both values if sites are not unique in locs - fix manually
dat_aus[dat_aus$Site=='Mudjimba Island',]$Lat<--26.61623
dat_aus[dat_aus$Site=='Mudjimba Island',]$Long<-153.1130
dat_aus[dat_aus$Site=='Inner Gneerings',]$Lat<--26.64829
dat_aus[dat_aus$Site=='Inner Gneerings',]$Long<-153.1835

# site names that did not line up or had NA for lat - fixed
unique(dat_aus[is.na(dat_aus$Lat),]$Site)
paste(unique(locs$Site.name)[!unique(locs$Site.name)  %in%  unique(dat_aus$Site)])
# abun to PA
dat_aus$PA<-1
#rm unused cols
dat_aus<-dat_aus[,- c(10:12)]

# sp in trait dataset filter?

library(labdsv)
library(ade4)

table(dat_aus$Year)

#dat_aus_yr<-dat_aus[dat_aus$Year>2015,]

#Remove dodgy sites from upon Maria's advice
# actually keep in as clustering doesn't change
#dat_aus_sub<-dat_aus[-which(dat_aus$Site %in% c("Pialba Shallow",
#    "Gataker High Diversity Site",'Big Woody Shallow', "Mudjimba Island Shallow" )),]
#dat_aus_sub$Site<-factor(dat_aus_sub$Site)

mat_aus<-matrify(data.frame(dat_aus$Site, dat_aus$Fish, dat_aus$PA))

bdist<-dist.binary(mat_aus, method=1) # jaccard dist

#shouldn't use ward centroid or median methods for jaccard dist
plot(hclust(bdist, 'single'))  #splits sites into 3 clusts, make cut at Flat Rock 28 deg S

#plot with lats
bdist_lat<-bdist
ordlats<-dat_aus%>%group_by(Site)%>% summarize_all(first)
 

specs_aus$AUS_trop<-ifelse(specs_aus$Species%in% dat_aus[dat_aus$Lat<'-25.5',]$Fish, 1, 0)

specs_aus$AUS_temp<-ifelse(specs_aus$Species%in% dat_aus[dat_aus$Lat>'-25.5',]$Fish, 1, 0)

length(which(rowSums(specs_aus[,15:16])==0)) # here are the 72 missing sp

df4<-specs_aus[,c(1,15, 16)]

write.csv(df4, 'C:/coral_fish/data/Australia/AUS_species_tropical_class.csv', quote=F, row.names=F)

### Section to identify max latitude of each tropical species 

# mat_jpn and mat_aus being used from code above.
ordlats_aus<-dat_aus%>%group_by(Site)%>% summarize_all(first)
specs_aus$max_lat<-apply(mat_aus, 2, function(x){min(ordlats_aus[which(x==1),]$Lat)})# min for Aus

# could sort based on 95th percentile but is more conservative,
# with only few sites just take max
#sort(v1, decreasing=T)[0.95*length(v1)]

ordlats_jpn<-dat%>%group_by(Name.x)%>% summarize_all(first)
specs_jpn$max_lat<-apply(mat_jpn, 2, function(x){max(ordlats_jpn[which(x==1),]$lat)})# max for Jpn


