## 29/11/19

## Renamed script focused on fish and site data. Code clusters sites based on
## fish presence-absence and calculates maximum latitudes of each species
## results are added to the trait master file

library(dplyr)
library(labdsv)
library(ade4)

specs<-read.csv('C:/coral_fish/data/Traits/JPN_AUS_RMI_CHK_MLD_TMR_trait_master_opt2.csv', h=T)

# JAPAN

# survey data
dat<-read.csv('C:/coral_fish/data/Japan/FishData_JP_2016_final.csv', h=T)

dat$SpeciesFish<-as.character(dat$SpeciesFish)

nrow(specs[which(specs$JPN_sp>0),] );length(unique(dat$SpeciesFish))
specs[which(specs$JPN_sp>0),]$Species[which(!specs[which(specs$JPN_sp>0),]$Species %in% dat$SpeciesFish)]

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

# site data

locs<-read.csv('C:/coral_fish/data/Japan/jp2015_16_waypoints.csv', h=T)

dat<-left_join(dat, locs, by=c('SiteID'= 'Site'))

# rm blanks
dat<-dat[-which(dat$SiteID==''),]
# abun to PA
dat$PA<-1

# create p/a matrix

mat_jpn<-matrify(data.frame(dat$Name.x, dat$SpeciesFish, dat$PA))

table(dat[dat$Year==2015,]$SiteID)
table(dat[dat$Year==2016,]$SiteID)

bdist<-dist.binary(mat_jpn, method=1) # jaccard dist

#shouldn't use ward centroid or median methods for jaccard dist
plot(hclust(bdist, 'single'))  #splits sites into 2 clusts from lat <31deg N

#make split and fill in specs master dataframe

specs$JPN_trop<-ifelse(specs$JPN_sp==0, NA, ifelse(specs$Species%in% dat[dat$lat<31,]$SpeciesFish, 1, 0))
specs$JPN_tran<-ifelse(specs$JPN_sp==0, NA, ifelse(specs$Species%in% dat[dat$lat>31,]$SpeciesFish, 1, 0))

specs$JPN_maxlat<-NA
ordlats_jpn<-dat%>%group_by(Name.x)%>% summarize_all(first)
specs$JPN_maxlat[which(specs$JPN_sp==1)]<-apply(mat_jpn, 2, 
                  function(x){max(ordlats_jpn[which(x==1),]$lat)})# max for Jpn


# AUSTRLIA

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

specs$AUS_trop<-ifelse(specs$AUS_sp==0, NA, ifelse(specs$Species%in% dat_aus[dat_aus$Lat<'-25.5',]$Fish, 1, 0))
specs$AUS_tran<-ifelse(specs$AUS_sp==0, NA, ifelse(specs$Species%in% dat_aus[dat_aus$Lat>'-25.5',]$Fish, 1, 0))

specs$AUS_maxlat<-NA
ordlats_aus<-dat_aus%>%group_by(Site)%>% summarize_all(first)
specs$AUS_maxlat[which(specs$AUS_sp==1)]<-apply(mat_aus, 2, 
                        function(x){min(ordlats_aus[which(x==1),]$Lat)})# # min for Aus

# write out
write.csv(specs, 'C:/coral_fish/data/Traits/JPN_AUS_RMI_CHK_MLD_TMR_trait_master_opt2_lats.csv', quote=F, row.names=F) 


# could sort max lats based on 95th percentile but is more conservative,
# with only few sites just take max
#sort(v1, decreasing=T)[0.95*length(v1)]


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

