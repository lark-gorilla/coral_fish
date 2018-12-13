# 22/10/2018

# preparation of fish trait and abundance data from Japan, Australia and RMI
# for several pieces of comparative lab analyses 

library(dplyr)
library(readxl)


setwd("M:/coral_fish")

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
aus_species_list<-read.csv('data/Australia/LongTransect_Subtropical_fish_list.csv')

#get list of aus_species names
aus_species<-aus_species_list$Fish

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

#now clean aus species
aus_species<-as.character(aus_species)

#add extra aus species 
aus_species<-c(aus_species, 'Scarus chameleon','Triaenodon obesus','Cheilinus aerenatus',
               'Scyliorhinidae','Priacanthus macracanthus','Priacanthus blochii','Cheilodipterus artus',
               'Pomacanthus imperator','Scarus flavipectoralis','Scarus rubroviolaceaus','Scarus fraenatus')

#species names wrong in aus data change/ edit
aus_species[aus_species== "Archamia zosterophora"]<-'Taeniamia zosterophora'
aus_species[aus_species=='Scarus rubroviolaceaus']<-'Scarus rubroviolaceus'
aus_species[aus_species=="Scarus fraenatus"]<-'Scarus frenatus'
aus_species[aus_species=='Cheilinus aerenatus']<-'Oxycheilinus arenatus'
aus_species[aus_species=='Archamia fucata']<-'Taeniamia fucata'
aus_species[aus_species=='Apogon limenus']<-'Ostorhinchus limenus'

#****** Chuuk data ******#

# read in Louise's 2016 Chuuk survey data, note these are mainly food-fish species so 
# do not fully represent community
chuuk_species_list<-read.csv('data/Chuuk/ChuukFishSpp.csv', h=T)
# and non-database species
chuuk_nonDB<-read.csv('data/Chuuk/UniqueFishChuuk.csv', h=T)
# note they may be in RMI trait database

#****** East Timor data ******#

#****** Maldives data ******#

# read in Australia/Japan trait database

#library(readr)
# base read.csv fails so use readr's read_csv
#bigtrait<-read_csv('~/leeds_postdoc/data/Traits/fish_traits.csv')

# Updated database from Maria 25 Oct

bigtrait<-read_excel('data/Traits/_database_index10_25Oct2018.xlsx', sheet=1, skip=1)

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

# subset to only spcies in one of the three regions

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

table(bigtrait2$ParentalMode)
bigtrait2[which(bigtrait2$ParentalMode=='Brooders'),]$ParentalMode<-'brooders'
bigtrait2[which(bigtrait2$ParentalMode=='nesters'),]$ParentalMode<-'Nesters'

# save file

write.csv(bigtrait2, 'data/Traits/JPN_AUS_RMI_trait_master.csv', quote=F, row.names=F) 



