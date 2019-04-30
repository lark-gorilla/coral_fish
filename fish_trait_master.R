# 22/10/2018

# preparation of fish trait and abundance data from Japan, Australia and RMI
# for several pieces of comparative lab analyses 

library(dplyr)
library(readxl)


setwd("C:/coral_fish")

#****** RMI data ******#

# read in coral and fish abundance and trait data

fishdat<-read_excel('data/RMI/RMIfish_siteDepthMeans_ADepth.xlsx', sheet=1)

fishtrait<-read.csv('data/Traits/RMI_fish_traits4.csv', h=T)

#test names line up

n1<-names(fishdat)[-1] # without site col name
n2<-fishtrait$Species

# FYI more species in trait dataset than abundance dataset

ingrid<-expand.grid(n1, n2)
ingrid$Var1<-as.character(ingrid$Var1)
ingrid$Var2<-as.character(ingrid$Var2)

ingrid$test=ingrid$Var1==ingrid$Var2

ingrid_out<- ingrid %>% group_by(Var1) %>%
  summarise(has_trait=TRUE%in%test)

filter(ingrid_out, has_trait==F)
# only Bryanops natans
# called Bryaninops natans in trait data - correct
names(fishdat)[names(fishdat)=='Bryanops natans']<-'Bryaninops natans'

# subset rmi triat data to only include species seen in surveys
rmi_species<-names(fishdat)

fishtrait<-fishtrait[fishtrait$Species %in% rmi_species,]



#****** Japan and Australia data ******#
# code from Katie C

#read in fish survey data japan 
fish_survey<-read.csv('data/Japan/FishData_JP_2016_final.csv')

#read in australia survey species
aus_species_list<-read.csv('data/Australia/maria_aus_list_20feb.csv')

#get list of aus_species names
aus_species<-aus_species_list$species

#check Japan fish survey data and remove blank rows (7534 onwards) and columns (12:14)
glimpse(fish_survey)

fish_survey<-fish_survey[1:7533,1:11]

#get list of japan species names
fish_survey$SpeciesFish<-as.character(fish_survey$SpeciesFish)

#remove duplicate species
#species name is wrong so change (JPN data )
fish_survey$SpeciesFish[fish_survey$SpeciesFish=="Apogon aureus"]<- "Ostorhinchus aureus"
fish_survey$SpeciesFish[fish_survey$SpeciesFish=='PLectroglyphidodon dickii']<-'Plectroglyphidodon dickii'
fish_survey$SpeciesFish[fish_survey$SpeciesFish=='Goniistius zebra']<-'Cheilodactylus zebra'
fish_survey$SpeciesFish[fish_survey$SpeciesFish=='Diagramma picta']<-'Diagramma pictum'
fish_survey$SpeciesFish[fish_survey$SpeciesFish=='Goniistius zonatus']<-'Cheilodactylus zonatus'
fish_survey$SpeciesFish[fish_survey$SpeciesFish=='Halichoeres poecilopterus']<-'Parajulis poecilepterus'
fish_survey$SpeciesFish[fish_survey$SpeciesFish=='Sebasticus marmoratus']<-'Sebastiscus marmoratus'
fish_survey$SpeciesFish[fish_survey$SpeciesFish=="Apogon doederleini"]<-'Ostorhinchus doederleini'
fish_survey$SpeciesFish[fish_survey$SpeciesFish=='Siganus stellatus']<-'Siganus punctatus'
fish_survey$SpeciesFish[fish_survey$SpeciesFish=="Apogon limenus"]<-'Ostorhinchus limenus'
fish_survey$SpeciesFish[fish_survey$SpeciesFish=="Chaetodon modestus"]<-'Roa modesta'

#now clean aus species
aus_species<-as.character(aus_species)


#****** Chuuk data ******#

# read in Louise's 2016 Chuuk survey data, note these are mainly food-fish species so 
# do not fully represent community
chuuk_species_list<-read.csv('data/Chuuk/ChuukFishSpp.csv', h=T)
# and non-database species
chuuk_nonDB<-read.csv('data/Chuuk/UniqueFishChuuk.csv', h=T)

chuuk_species<-chuuk_species_list$x

chuuk_nonDB

chuuk_species<-c(as.character(chuuk_species), 'Parupeneus trifasciatus') # added fom non-bd species list
# note they may be in RMI trait database

#****** East Timor data ******#

timor_species_list<-read_excel('data/Etimor/FishData_Maldives_Mar2015_list.xlsx', sheet=1)

timor_species<-timor_species_list$Species

#****** Maldives data ******#

maldives_species_list<-read_excel('data/Maldives/FishData_Maldives_Mar2015_list.xlsx', sheet=1)

maldives_species<-maldives_species_list$Species
#preliminary test shows all maldives species are in the trait database

# read in Australia/Japan trait database

#library(readr)
# base read.csv fails so use readr's read_csv
#bigtrait<-read_csv('~/leeds_postdoc/data/Traits/fish_traits.csv')

# Updated database from Maria 20 Feb

bigtrait<-read_excel('data/Traits/_database_index11_19Feb2019.xlsx', sheet=1, skip=1)

#remove extra rows from fish_trait
#bigtrait<-bigtrait[1:1122,]

#remove unnecessary traits
bigtrait<- bigtrait[,c(2, 11,17,21,30,36,38,39,41)]

# combine with RMI triat data

names(bigtrait)
#rename to better names 
names(bigtrait)<-c('Species', 'ThermalAffinity', 'BodySize','DepthRange',
                   'PLD', 'Diet', 'Aggregation', 'Position', 'ParentalMode')


# check which RMI species are already in the japan/aus trait dataset

rmi_sp<- fishtrait$Species[-which(fishtrait$Species %in% bigtrait$Species)]

fishtrait[fishtrait$Species %in% rmi_sp,]
# all bar two

# checking eac against xlsx file and updating species list name from 
# survey data to be correct with bigtrait species name

rmi_species[rmi_species=='Chaetodon rafflesi']<-'Chaetodon rafflesii'
rmi_species[rmi_species=='Zebrasoma veliferum']<-'Zebrasoma velifer'

# note there are some discrepencies between traits (trophic and aggregation)
# between the aus/japan train master database (bigtrait) and the RMI trait database (fishtrait)
# however RMI traits do not include depth. For simplicity we take all species traits from 
# the master dataset (bigtrait) - furthermore this hs been manually corrected in places where the former
# has not.

# Add columns for each survey region

bigtrait$JPN_sp<-ifelse(bigtrait$Species %in% fish_survey$SpeciesFish, 1, 0)
bigtrait$AUS_sp<-ifelse(bigtrait$Species %in% aus_species, 1, 0)
bigtrait$RMI_sp<-ifelse(bigtrait$Species %in% rmi_species, 1, 0)
bigtrait$CHK_sp<-ifelse(bigtrait$Species %in% chuuk_species, 1, 0)
bigtrait$MLD_sp<-ifelse(bigtrait$Species %in% maldives_species, 1, 0)
bigtrait$TMR_sp<-ifelse(bigtrait$Species %in% timor_species, 1, 0)

# subset to only spcies in one of the three regions - change if want to include Chuuk, Maldives and Timor
bigtrait2<-bigtrait[which(rowSums(bigtrait[,10:12])>0),]

# Now to tidy the data

table(bigtrait2$ThermalAffinity)

table(bigtrait2$BodySize)
summary(bigtrait2$BodySize)

table(bigtrait2$DepthRange)
summary(bigtrait2$DepthRange)

table(bigtrait2$PLD)
summary(bigtrait2$PLD)

table(bigtrait2$Diet)

table(bigtrait2$Aggregation)
bigtrait2[which(bigtrait2$Aggregation=='harems'),]$Aggregation<-'groups'

table(bigtrait2$Position)
bigtrait2[which(bigtrait2$Position=='AlgaeAssociated'),]$Position<-'Benthic'
bigtrait2[which(bigtrait2$Position=='CnidarianAssociated'),]$Position<-'BenthicSpecialist'
bigtrait2[which(bigtrait2$Position=='EchinodermAssociated'),]$Position<-'BenthicSpecialist'
bigtrait2[which(bigtrait2$Position=='SandAssociated'),]$Position<-'BenthicSpecialist'

# option 2
bigtrait2[which(bigtrait2$Position=='BenthicSpecialist'),]$Position<-'Demersal'

table(bigtrait2$ParentalMode)
bigtrait2[which(bigtrait2$ParentalMode=='Brooders'),]$ParentalMode<-'brooders'
bigtrait2[which(bigtrait2$ParentalMode=='nesters'),]$ParentalMode<-'Nesters'
bigtrait2[which(bigtrait2$ParentalMode=='Viviparous'),]$ParentalMode<-'Live bearers'


# save file

write.csv(bigtrait2, 'data/Traits/JPN_AUS_RMI_CHK_MLD_TMR_trait_master_opt2.csv', quote=F, row.names=F) 




##########################################
# Check to see how many species are shared between Aus/JPN

dat<-read.csv('C:/coral_fish/data/Traits/JPN_AUS_RMI_CHK_MLD_TMR_trait_master_opt2.csv', h=T)

dat<-dat[which(dat$AUS_sp>0 | dat$JPN_sp>0),] # we will focus on Australia for this prelim

# remove functional duplicates
dup_trait<-paste(dat$BodySize, dat$DepthRange, dat$PLD, dat$Diet, dat$Aggregation, dat$Position, 
                 dat$ParentalMode)
dat<-dat[-which(duplicated(dup_trait)),]
dat$Species<-as.character(dat$Species)
dat[which(dat$Species=='Scarus psittacus'),]$Species<-'Scarus psittacus/spinus' # edit for one

# Including those that default to duplicates via NA after gower dist
dat<-dat[-which(dat$Species=='Caesio sp.'),]

row.names(dat)<-dat$Species

## Edit to some trait values from MB 25/10/18
dat[dat$Species=='Brotula multibarbata',]$DepthRange<-219
dat[dat$Species=='Mobula birostris',]$BodySize<-450
dat[dat$Species=='Amphiprion sandaracinos',]$BodySize<-14

nrow(dat) #689
nrow(dat[which(dat$AUS_sp>0 & dat$JPN_sp>0 ),]) #229
nrow(dat[which(dat$AUS_sp>0 & dat$JPN_sp==0 ),])#299
nrow(dat[which(dat$AUS_sp==0 & dat$JPN_sp>0 ),])#161

dist1<-daisy(dat[,3:9], metric='gower', stand = FALSE)

