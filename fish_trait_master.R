# 22/10/2018

# preparation of fish trait and abundance data from Japan, Australia and RMI
# for several pieces of comparative lab analyses

library(dplyr)
library(readxl)

if(Sys.info()['nodename']=="FBS5PCW223"){
  setwd("M:/coral_fish")}else{
    setwd("~/leeds_postdoc")}

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

# read in Australia/Japan trait database

#library(readr)
# base read.csv fails so use readr's read_csv
#bigtrait<-read_csv('~/leeds_postdoc/data/Traits/fish_traits.csv')

# Updated database from Maria 22 Oct


bigtrait<-read_excel('data/Traits/_database_index10_22Oct2018.xlsx', sheet=1, skip=1)

#remove extra rows from fish_trait
#bigtrait<-bigtrait[1:1122,]

#remove unnecessary traits
bigtrait<- bigtrait[,c(2, 11,17,21,30,36,38,39,41)]

# combine with RMI triat data

names(bigtrait)
names(fishtrait)
fishtrait$Family<-NULL
fishtrait$HomeRange<-NULL
fishtrait$PD50<-NULL
fishtrait$TrophLevel<-NULL
fishtrait$Food<-NULL
fishtrait$SpawnMode<-NULL
fishtrait$Active<-NULL
fishtrait$Envtemp<-'tropical'

# check which RMI species are already in the japan/aus trait dataset

rmi_sp<- fishtrait$Species[-which(fishtrait$Species %in% bigtrait$Species)]

fishtrait[fishtrait$Species %in% rmi_sp,]
# all bar two

# checking eac against xlsx file and updating species list name from 
# survey data to be correct with bigtrait species name

'Chaetodon rafflesii' %in% c(fish_survey$SpeciesFish, aus_species) # not in so not altered

'Zebrasoma velifer' %in% c(fish_survey$SpeciesFish, aus_species) # is in so changed

rmi_species[rmi_species=='Zebrasoma veliferum']<-'Zebrasoma velifer'

# note there are some discrepencies between traits (trophic and aggregation)
# between the aus/japan train master database (bigtrait) and the RMI trait database (fishtrait)
# for RMI species that DO NOT occur in the aus/japan data we use
# the RMI triats, for shared species we use the master trait dataset

# subset rmi_traits to only rmi species outstanding from jap/aus surveys
rmi_only<-rmi_species[-which(rmi_species %in% c(fish_survey$SpeciesFish, aus_species))]

fishtrait_rmi<-filter(fishtrait, Species %in% rmi_only)


#rename to better names 
names(bigtrait)<-c('Species', 'ThermalAffinity', 'BodySize','DepthRange',
                   'PLD', 'Diet', 'Aggregation', 'Position', 'ParentalMode')

#filling in the blanks#
#first make matrix 
rownames(fish_trait_c)<-fish_trait$Species

#fill in blanks with NA
fish_trait_c[fish_trait_c == '']<-NA
fish_trait_c[fish_trait_c == ' ']<-NA


#make into numeric/ factor/ character
sapply(fish_trait_c, class) #ok
########

###do traits## 

#filter for just japan species 
jp_species<-c(colnames(fish_matrix))
jp_species<-factor(jp_species)

fish_trait<- filter(fish_trait, Species %in% jp_species )

tropicalising_japan<- filter(tropicalising_species, species %in% jp_species)

#Mouillot traits + extra relevant 
colnames(fish_trait)





#check for errors in the factors
fish_trait_c$Diet<-factor(fish_trait_c$Diet)
levels(fish_trait_c$Diet) 

fish_trait_c$Position<-factor(fish_trait_c$Position)
levels(fish_trait_c$Position)

fish_trait_c$Aggregation<-factor(fish_trait_c$Aggregation)
levels(fish_trait_c$Aggregation)

fish_trait_c$Aggregation[fish_trait_c$Aggregation =='harems']<-'groups'
fish_trait_c$Aggregation<-factor(fish_trait_c$Aggregation)

fish_trait_c$ParentalMode<-factor(fish_trait_c$ParentalMode)
levels(fish_trait_c$ParentalMode)
which(fish_trait_c$ParentalMode=='Nesters') #right one
which(fish_trait_c$ParentalMode=='nesters')
fish_trait_c$ParentalMode[fish_trait_c$ParentalMode=='nesters']<-'Nesters'
fish_trait_c$ParentalMode<-factor(fish_trait_c$ParentalMode)

fish_trait_c$ThermalAffinity<-factor(fish_trait_c$ThermalAffinity)
levels(fish_trait_c$ThermalAffinity)


#check species are identical to fish matrix
which(rownames(fish_trait_c) %in% colnames(fish_matrix))

identical(rownames(fish_trait_c), colnames(fish_matrix)) #not identical

fish_trait_c<-fish_trait_c[ order(rownames(fish_trait_c)),]


identical(rownames(fish_trait_c), colnames(fish_matrix)) #true







# Subsetting to only include traits from species from the three regions

# filter to only include traits that are in abundance data

fishtrait2<-filter(fishtrait, Species%in% names(fishdat))

row.names(fishtrait2)<-fishtrait2$Species

# drop trait 'PD50' as seems to cock up calculaion
fishtrait2<-fishtrait2[,c(3:5,7:length(fishtrait2))]
# could also drop 'Function' as coarse representative of 'Food'?


##~~~