## 29/11/19

## Renamed script focused on fish and site data. Code clusters sites based on
## fish presence-absence and calculates maximum latitudes of each species
## results are added to the trait master file

library(dplyr)
library(labdsv)
library(ade4)
library(vegan)
library(reshape2)
library(readxl)
library(iNEXT) # installed from tar.gz from CRAN archive
#http://johnsonhsieh.github.io/iNEXT/
library(tidyverse)

# function to calc species acc curves for set threshold broken down by FG or thermalaffinity2
# requires site, species, abun matrix (dat), trait database with FGs assigned (fgs), species in dat and fg must match
spacBreakdown<-function(data=dat, fgs=fgs, FGZ=c(1,2), TM2=c('tropical', 'subtropical'), thresh=200)
{
  # FG or Tm2 can be 'all'
  runs=expand.grid(FGZ, TM2)  
  runs$Var2<-as.character(runs$Var2)
  return_df<-NULL
  
  for(i in 1:nrow(runs))
  {
    FGi=runs[i,]$Var1
    TMi=runs[i,]$Var2
    temp<-NULL 
    
    if(FGi=='all' & TMi!='all'){   
      temp<-t(data[,which(names(data) %in% 
      fgs[fgs$ThermalAffinity2==TMi,]$Species)])}
    if(TMi=='all'& FGi!='all'){   
      temp<-t(data[,which(names(data) %in% 
      fgs[fgs$FG==FGi,]$Species)])}
    if(TMi=='all'& FGi=='all'){   
      temp<-t(data)}
    if(TMi!='all'& FGi!='all'){
      temp<-t(data[,which(names(data) %in% 
      fgs[fgs$FG==FGi & fgs$ThermalAffinity2==TMi,]$Species)])}
    #remove sites with 0 individuals from rarefaction
    if(length(which(colSums(temp)==0))==0){
      temp2<-temp}else{
        temp2<-temp[,-which(colSums(temp)==0)]}
    
    out<-iNEXT(temp2, q=0, size=seq(10, 500, 10), datatype="abundance")
    plot(out, se=F)
    
    out<-do.call(rbind, lapply(out$iNextEst, function(x){data.frame(x[which(x$m==thresh), c(1,2,4:6)])}))
    out$Site=row.names(out)
    if(length(which(colSums(temp)==0))>0){
      out<-rbind(out, data.frame(m=thresh, method='zero.abun',
                                 qD=0, qD.LCL=0, qD.UCL=0, Site=names(which(colSums(temp)==0))))}
    out$FG=FGi
    out$ThermalAffinity2=TMi
    out$SPRICraw<-apply(t(temp), 1, function(x){length(x[which(x>0)])})[out$Site]
    return_df<-rbind(return_df, out)
    print(runs[i,])
  }
  return(return_df)
}

## Read in data

# Species trait data
specs<-read.csv('C:/coral_fish/data/Traits/JPN_AUS_RMI_CHK_MLD_TMR_trait_master_opt2.csv', h=T)

# species weight equation data
wtlen<-read_xlsx('C:/coral_fish/data/fish/fish_weight_length_calc_a_and_b_edited.xlsx', sheet=1)
#warning fine
# for some reason there are multiple/duplicate entries for some species
#wtlen[wtlen$SpeciesName%in%wtlen[which(duplicated(wtlen$SpeciesName)==TRUE),]$SpeciesName,]%>%View()
# remove
wtlen<-wtlen[-which(duplicated(wtlen$SpeciesName)),]
specs[-which(specs$Species %in% wtlen$SpeciesName),] %>% filter(JPN_sp==1 | AUS_sp==1)
# fixed naming issues
# read in FGs
fgs<-read.csv('C:/coral_fish/data/Traits/JPN_AUS_RMI_CHK_MLD_TMR_trait_master_opt2_lats_FG.csv')

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

# Species accumulation curve
#sumarise abundance per species per site, doing by transect sets the min abundance too low (6)
jpn_abun_mat<-dat %>% group_by(SiteID, SpeciesFish) %>% summarise(sum_abun=sum(Number))

jpn_abun_mat<-matrify(data.frame(jpn_abun_mat$SiteID, jpn_abun_mat$SpeciesFish, jpn_abun_mat$sum_abun))

jpn_spac<-spacBreakdown(data=jpn_abun_mat, fgs=fgs, FGZ=c(1,2,4,6), TM2=c('tropical', 'subtropical'), thresh=100)

jpn_spac2<-left_join(jpn_spac, locs[,2:4], by='Site')

ggplot(jpn_spac2, aes(x = lat, y = qD, colour=ThermalAffinity2)) + 
  geom_point()+geom_smooth(se=F)+geom_hline(yintercept=0)+geom_hline(yintercept=5, linetype='dotted')+
  geom_vline(xintercept=31, linetype='dotted')+facet_wrap(~FG, scales='free')+theme_bw()

# check correlation between spacc output and sprich from data in Japan

ggplot(jpn_spac2, aes(x = SPRICraw, y = qD)) + 
  geom_point()+geom_smooth(se=F)+facet_wrap(~FG+ThermalAffinity2, scales='free')+theme_bw()


# abun to PA
dat$PA<-1

# create p/a matrix

mat_jpn<-matrify(data.frame(dat$Name.x, dat$SpeciesFish, dat$PA))

table(dat[dat$Year==2015,]$SiteID)
table(dat[dat$Year==2016,]$SiteID)

bdist<-dist.binary(mat_jpn, method=1) # jaccard dist

#shouldn't use ward centroid or median methods for jaccard dist
plot(hclust(bdist, 'single'))  #splits sites into 2 clusts from lat <31deg N


# biomass calculations and conversion to biomass per transect effort
# add FG to data
dat<-left_join(dat, fgs[,c(1,21,22)], by=c('SpeciesFish'='Species'))
# add mass calc columns to data
dat<-left_join(dat, wtlen[,c(1:3)], by=c('SpeciesFish'='SpeciesName'))

biom1<-dat%>%group_by(SiteID, FG, ThermalAffinity2)%>%
  summarise(tot_biom=sum(Number*(a*SizeCm^b), na.rm=T))%>%ungroup()%>%
  complete(SiteID, FG, ThermalAffinity2, fill=list(tot_biom=0))

biom1<-left_join(biom1, dat%>%group_by(SiteID)%>%
                   summarise(totMsurv=length(unique(Transect))*25), 
                 by='SiteID')

# summarise total biomass (abundance * length) per species per transect
# NOT needed in Japan: then for transects that are surveyed multiple times take the mean biomass per species 
# then take the mean of averaged transects within each site per species
# might want to take min or max instead of mean?

biom1<-dat%>%group_by(SpeciesFish, Name.x, Transect)%>%summarise(tot_biom=sum(Number*(a*SizeCm^b), na.rm=T))

biom3<-biom1%>%group_by(SpeciesFish, Name.x)%>%summarise(mean_trans_biom=median(tot_biom))

mat_biom_jpn<-matrify(data.frame(biom3$Name.x, biom3$SpeciesFish, biom3$mean_trans_biom))

# curiosity site cluster plot based on biomass
plot(hclust(vegdist(mat_biom_jpn, 'bray', na.rm=T), 'average')) 


# setup data.frame to make geom_violin of p/a and biomass of each FG at each site
#pa data
mat_jpn_df<-as.data.frame(mat_jpn)
mat_jpn_df$Name.x<-row.names(mat_jpn_df)
mat_jpn_df<-left_join(mat_jpn_df, dat%>%group_by(Name.x)%>%summarise(Lat=first(lat)),
                      by='Name.x')
mat_jpn_df<-melt(mat_jpn_df, id.vars = c('Name.x', 'Lat'))

#biomass data
mat_biom_jpn_df<-as.data.frame(mat_biom_jpn)
mat_biom_jpn_df$Name.x<-row.names(mat_biom_jpn_df)
mat_biom_jpn_df<-melt(mat_biom_jpn_df, id.vars = c('Name.x'))

#combine

names(mat_jpn_df)[3]<-'Species'
names(mat_jpn_df)[4]<-'pa'
mat_jpn_df$biom<-mat_biom_jpn_df$value
mat_jpn_df$Name.x<-gsub(',', '@', mat_jpn_df$Name.x)

write.csv(mat_jpn_df, 'C:/coral_fish/data/Japan/Jpn_sites_pa_biomass_median.csv', quote=F, row.names=F) 


# AUSTRALIA

#read in australia survey species
dat_aus<-read.csv('C:/coral_fish/data/Australia/LongTransect_Subtropical_fish_Sep2018.csv')

#get list of aus_species names

dat_aus$Fish<-as.character(dat_aus$Fish)

nrow(specs[which(specs$AUS_sp>0),] );length(unique(dat_aus$Fish)) # OK only 2 different
specs[which(specs$AUS_sp>0),] [which(!specs[which(specs$AUS_sp>0),]  %in% dat_aus$Fish)]

#remove species that occur in sampling data but not species/traits list
badfish<-unique(dat_aus$Fish[which(!dat_aus$Fish %in% specs[specs$AUS_sp>0,]$Species)])
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

# Need to remove duplicate location names otherwise left_join inflates
# survey data with duplicate sites!

locs<-locs[-which(duplicated(locs$Site.name)),]

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

#Remove dodgy sites from upon Maria's advice 02/12/19

dat_aus_sub<-filter(dat_aus, !Site %in% "Mudjimba Island Shallow")

dat_aus_sub[which(is.na(dat_aus_sub$Number)),]$Number<-1 #fix for sp with na number

### Remove Summer sampling records, so all Aussie sites represent winter
dat_aus_sub<-filter(dat_aus_sub, grepl('Win', dat_aus_sub$Trip))

# Fix to fill 3 NA fish size measures 
dat_aus_sub[which(is.na(dat_aus_sub$Size)),]$Size<-c(7, 16, 7)



# Species accumulation curve
#sumarise abundance per species per site, doing by transect sets the min abundance too low (6)
aus_abun_mat<-dat_aus_sub %>% group_by(Site, Fish) %>% summarise(sum_abun=sum(Number))

aus_abun_mat<-matrify(data.frame(aus_abun_mat$Site, aus_abun_mat$Fish, aus_abun_mat$sum_abun))

aus_spac<-spacBreakdown(data=aus_abun_mat, FGZ=c(1,2,4,6), TM2=c('tropical', 'subtropical'), thresh=100)

aus_spac2<-left_join(aus_spac, locs[,1:6], by=c('Site'='Site.name'))

ggplot(aus_spac2, aes(x = Lat, y = qD, colour=ThermalAffinity2)) + 
  geom_point()+geom_smooth(se=F)+geom_hline(yintercept=0)+geom_hline(yintercept=5, linetype='dotted')+
  geom_vline(xintercept=-25.5, linetype='dotted')+scale_x_reverse()+
  facet_wrap(~FG, scales='free')+theme_bw()

# make hclust visualisation
mat_aus<-matrify(data.frame(paste(dat_aus_sub$Site,dat_aus_sub$Trip), dat_aus_sub$Fish, dat_aus_sub$PA))

bdist<-dist.binary(mat_aus, method=1) # jaccard dist

#shouldn't use ward centroid or median methods for jaccard dist
plot(hclust(bdist, 'single'))  #splits sites into 3 clusts, make cut at Flat Rock 28 deg S

specs$AUS_trop<-ifelse(specs$AUS_sp==0, NA, ifelse(specs$Species%in% dat_aus[dat_aus$Lat<'-25.5',]$Fish, 1, 0))
specs$AUS_tran<-ifelse(specs$AUS_sp==0, NA, ifelse(specs$Species%in% dat_aus[dat_aus$Lat>'-25.5',]$Fish, 1, 0))

specs$AUS_maxlat<-NA
ordlats_aus<-dat_aus_sub%>%group_by(Site)%>% summarize_all(first)
specs$AUS_maxlat[which(specs$AUS_sp==1)]<-apply(mat_aus, 2, 
                        function(x){min(ordlats_aus[which(x==1),]$Lat)})# # min for Aus
# write out
write.csv(specs, 'C:/coral_fish/data/Traits/JPN_AUS_RMI_CHK_MLD_TMR_trait_master_opt2_lats.csv', quote=F, row.names=F) 

# add mass calc columns to data

dat_aus_sub<-left_join(dat_aus_sub, wtlen[,c(1:3)], by=c('Fish'='SpeciesName'))

# biomass weighting

# summarise total biomass (abundance * length) per species per transect
# then for transects that are surveyed multiple times take the mean biomass per species
# then take the mean of averaged transects within each site per species
# might want to take min or max instead of mean?

biom1<-dat_aus_sub%>%group_by(Fish, Site, id)%>%summarise(tot_biom=sum(Number*(a*Size^b), na.rm=T))
biom1$id_nodate<-unlist(lapply(strsplit(as.character(biom1$id), '_'), function(x){paste(x[1], x[2], sep='_')}))
biom2<-biom1%>%group_by(Fish, Site, id_nodate)%>%summarise(mean_tot_biom=median(tot_biom),
              min_tot_biom=min(tot_biom),max_tot_biom=max(tot_biom))
biom3<-biom2%>%group_by(Fish, Site)%>%summarise(mean_trans_biom=median(mean_tot_biom))

mat_biom_aus<-matrify(data.frame(biom3$Site, biom3$Fish, biom3$mean_trans_biom))

# curiosity site cluster plot based on biomass
plot(hclust(vegdist(mat_biom_aus, 'bray', na.rm=T), 'average')) 


# setup data.frame to make geom_violin of p/a and biomass of each FG at each site
#pa data
mat_aus_df<-as.data.frame(mat_aus)
mat_aus_df$Site<-row.names(mat_aus_df)
mat_aus_df<-left_join(mat_aus_df, dat_aus_sub%>%group_by(Site)%>%summarise(Lat=first(Lat)),
                      by='Site')
mat_aus_df<-melt(mat_aus_df, id.vars = c('Site', 'Lat'))

#biomass data
mat_biom_aus_df<-as.data.frame(mat_biom_aus)
mat_biom_aus_df$Site<-row.names(mat_biom_aus_df)
mat_biom_aus_df<-melt(mat_biom_aus_df, id.vars = c('Site'))

#combine

names(mat_aus_df)[3]<-'Species'
names(mat_aus_df)[4]<-'pa'
mat_aus_df$biom<-mat_biom_aus_df$value

write.csv(mat_aus_df, 'C:/coral_fish/data/Australia/Aus_sites_pa_biomass_median.csv', quote=F, row.names=F) 


### Correcting for uneven sampling

# have a look at each transects min abun 
S <- specnumber(jpn_abun_mat) # observed number of species
(raremax <- min(rowSums(jpn_abun_mat)))
Srare <- rarefy(jpn_abun_mat, raremax)
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
rarecurve(jpn_abun_mat, step = 20, sample = raremax, col = "blue", cex = 0.6)


# could sort max lats based on 95th percentile but is more conservative,
# with only few sites just take max
#sort(v1, decreasing=T)[0.95*length(v1)]



