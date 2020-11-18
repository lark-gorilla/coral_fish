# paper code 23/01/2020

# picks up from stage where optimum number of clusters are known
#Steps
#1 calculate redundancy for each community and group and thermal preference
#2 plot and analyse changes in sprich and biomass of tropical/non-tropical FGs over latitude
#3 plot functional niche overlap of tropical/non-tropical components of each FG

#rm(list=ls())
library(dplyr)
library(ggplot2)
library(gridExtra)
library(mice)
library(scales)
library(ggwordcloud)
library(vegan)
library(ade4)
library(readxl)
library(mgcv)
library(emmeans)
library(ggResidpanel)
library(adehabitatHR)
library(sf)
library(cluster)

#################### read data ##########################
##########################################################

dat<-read.csv('C:/coral_fish/data/Traits/JPN_AUS_RMI_CHK_MLD_TMR_trait_master_opt2_clusters.csv', h=T)
# Running on simplist classification of Position trait

#Remove Aus summer only species
aus_summer<-read.csv('C:/coral_fish/data/Australia/sp_list_summer_only.csv')
dat[dat$Species %in% aus_summer$Fish,]$AUS_sp<-0

# Read in spacc-sprich and biomass data
spr_jpn<-read.csv('C:/coral_fish/data/Japan/Jpn_sites_sprich_combos.csv')
bio_jpn<-read.csv('C:/coral_fish/data/Japan/Jpn_transects_biomass.csv')

spr_aus<-read.csv('C:/coral_fish/data/Australia/Aus_sites_sprich_combos.csv')
bio_aus<-read.csv('C:/coral_fish/data/Australia/Aus_transects_biomass.csv')

# read in species per site biomass data for PCoA
aus_sp_site<-read.csv('C:/coral_fish/data/Australia/Aus_site_species_biomass.csv')
jpn_sp_site<-read.csv('C:/coral_fish/data/Japan/Jpn_site_species_biomass.csv')

# FG to factor
bio_jpn$FG<-factor(bio_jpn$FG)
bio_aus$FG<-factor(bio_aus$FG)

# Remove bad 2 sites from Jpn and Aus
bio_jpn<-filter(bio_jpn, !SiteID %in% c('JP27', 'JP32'))
bio_aus<-filter(bio_aus, !Site %in% c('Flat Rock', 'Wolf Rock'))

spr_jpn<-filter(spr_jpn, !Site %in% c('JP27', 'JP32'))
spr_aus<-filter(spr_aus, !Site %in% c('Flat Rock', 'Wolf Rock'))

jpn_sp_site<-filter(jpn_sp_site, !SiteID %in% c('JP27', 'JP32'))
aus_sp_site<-filter(aus_sp_site, !Site %in% c('Flat Rock', 'Wolf Rock'))

# correct biomass
bio_aus$cor_biom<-bio_aus$tot_biom/bio_aus$totMsurv
bio_jpn$cor_biom<-bio_jpn$tot_biom/bio_jpn$totMsurv

aus_sp_site$cor_biom<-aus_sp_site$tot_biom/aus_sp_site$totMsurv
jpn_sp_site$cor_biom<-jpn_sp_site$tot_biom/jpn_sp_site$totMsurv

# and summarise by site for sp site data
aus_sp_site<-aus_sp_site%>%group_by(Site, Lat, FG, ThermalAffinity2, Fish)%>%
  summarise(cor_biom=sum(cor_biom))

jpn_sp_site<-jpn_sp_site%>%group_by(SiteID, lat, FG, ThermalAffinity2, SpeciesFish)%>%
  summarise(cor_biom=sum(cor_biom))

#### Functional Entity creation and word clouds ####

dat_mice<-mice(dat[,c(3:9)], m=5, method=c(rep('norm.predict', 3), rep('polyreg', 4)))
dat_imp<-mice::complete(dat_mice)
dat_imp<-cbind(Species=dat[,1], ThermalAffinity=dat[,2], dat_imp, dat[,10:length(dat)])

#manual edit 
dat_imp[dat_imp$Species=='Cantheschenia grandisquamis',]$Diet<-'Omnivore'

# Bodysize

# 6 classes on cuts from Mouillot et al (2014)
dat_imp$BodySize <- cut(dat_imp$BodySize ,breaks = c(-Inf, 7,15,30,50, 80, Inf),
                     labels = c("Tiny","VSmall", "Small", "Medium", "Large", "VLarge"))
# DepthRange

ggplot(data=dat_imp, aes(x=DepthRange))+geom_histogram(binwidth=10, col='black')+
  scale_x_continuous(breaks=seq(0, 420, 10))
# split into those that can only dive to coral reef (~27m) depths, 
# the photic zone (~100m), and deeper > 100m. loosely similar to Mouillot et al (2014) 'Water column' trait
dat_imp$DepthRange <- cut(dat_imp$DepthRange ,breaks = c(-Inf, 30,100, Inf),
                       labels = c("shallow", "mid-depth", "deep"))
dat_imp$FE<-paste(dat_imp$BodySize, dat_imp$DepthRange,
                       dat_imp$Diet, dat_imp$Position, dat_imp$Aggregation)

# How to represent what each FG is using wordclouds of FEs

# reformat data
# split dataframe into list based on rows

FEdat<-dat_imp %>% group_by(FE) %>% summarise(FG=unique(groupk19),num=n())

FEdatlist<-split(FEdat, 1:nrow(FEdat))

FEwordlist<-lapply(FEdatlist, function(x){data.frame(FEcomp=unlist(strsplit(x$FE, ' ')),
                                                     n=x$num, FG=x$FG, FE=x$FE)})
FEword.df<-do.call('rbind', FEwordlist)

FEword.agg<-FEword.df %>% group_by(FG, FEcomp) %>% summarize(sum_word=sum(n)) 

# arranges FEcomp in descending order of n, by FG
FEword.agg<-FEword.agg %>% group_by(FG) %>% arrange(desc(sum_word), .by_group=TRUE)

# help with colouring https://stackoverflow.com/questions/18902485/colored-categories-in-r-wordclouds
FEword.agg$colorlist='blue'
FEword.agg$colorlist<-ifelse(FEword.agg$FEcomp%in%unique(dat_imp$Diet),
                             'orange', FEword.agg$colorlist)
FEword.agg$colorlist<-ifelse(FEword.agg$FEcomp%in%unique(dat_imp$Aggregation),
                             'green', FEword.agg$colorlist)
FEword.agg$colorlist<-ifelse(FEword.agg$FEcomp%in%unique(dat_imp$DepthRange),
                             'purple', FEword.agg$colorlist)
FEword.agg$colorlist<-ifelse(FEword.agg$FEcomp%in%unique(dat_imp$BodySize),
                             'red', FEword.agg$colorlist)

# preserves scaling i.e. sizes of FGs
ggplot(FEword.agg, aes(label=FEcomp, size=sum_word, colour=colorlist))+
  geom_text_wordcloud()+scale_size_area()+facet_wrap(~FG)

FEword.agg.cl<-split(FEword.agg, FEword.agg$FG)

# include seed
out<-lapply(FEword.agg.cl[as.numeric(names(sort(table(dat$groupk19), decreasing = T)))], function(x){
  ggplot(x, aes(label=FEcomp, size=sum_word, colour=colorlist))+
    geom_text_wordcloud(seed=300)+scale_size_area()+theme_minimal()+
    labs(title=paste('FG', unique(x$FG), 'n=', sum(x$sum_word)/5))})

do.call('grid.arrange', out) 

ggsave('C:/coral_fish/outputs/wordclouds.eps',
       plot=do.call('grid.arrange', out),width = 20, height = 20, units = "cm")

# get names for future plots
FEword.agg%>%group_by(FG)%>%
mutate(totword=sum(sum_word), cumword=cumsum(sum_word), wordperc=cumword/totword*100)%>%
filter(wordperc<59) ->FGnames

FGnames<-as.data.frame(aggregate(FEcomp~FG, FGnames, 
                                 function(x){paste(x, collapse='-')}))


#### Redundancy plots with thermal affinity split ####

red_dat<-rbind(dat[dat$AUS_sp>0,] %>% group_by(groupk19) %>% summarise(num=n()) %>%
  arrange(num) %>% mutate(val=18:1, dat='Australia'),
  dat[dat$JPN_sp>0,] %>% group_by(groupk19) %>% summarise(num=n()) %>%
    arrange(num) %>% mutate(val=18:1, dat='Japan'))

red_dat_therm<-rbind(dat[dat$AUS_sp>0,] %>% group_by(ThermalAffinity2, groupk19) %>%
                       summarise(num=n()) %>% mutate(dat='Australia'),
                     dat[dat$JPN_sp>0,] %>% group_by(ThermalAffinity2, groupk19) %>%
                       summarise(num=n())%>% mutate( dat='Japan'))

red_dat_therm$val<-left_join(red_dat_therm, red_dat, by=c('dat', 'groupk19'))$val

ggplot()+
  geom_bar(data=filter(red_dat_therm, dat=='Japan'), aes(x=val, y=num, fill=ThermalAffinity2),
           stat='identity',position=position_stack())+
  geom_hline(data=filter(red_dat_therm, dat=='Japan')%>%group_by(ThermalAffinity2)%>%summarise_all(mean),
             aes(yintercept=num, colour=ThermalAffinity2), linetype='dashed', colour=c('#2c7bb6', '#fdae61'))+
  scale_x_continuous(breaks=1:18, labels=FGnames[rev(filter(red_dat, dat=='Japan')$groupk19),]$FEcomp)+
  theme_bw()+theme(legend.justification = c(1, 1), legend.position = c(1, 1), 
        legend.background = element_rect(colour = "black"), 
        axis.text.x = element_text(angle = 60, hjust = 1))+
  guides(fill=guide_legend(title="Thermal Affinity"))+
  scale_fill_manual(values = c('#2c7bb6', '#fdae61'))+
  ylab('# of species per FG')+xlab('Functional niche')

ggplot()+
  geom_bar(data=filter(red_dat_therm, dat=='Australia'), aes(x=val, y=num, fill=ThermalAffinity2),
           stat='identity',position=position_stack())+
  geom_hline(data=filter(red_dat_therm, dat=='Australia')%>%group_by(ThermalAffinity2)%>%summarise_all(mean),
             aes(yintercept=num, colour=ThermalAffinity2), linetype='dashed', colour=c('#2c7bb6', '#fdae61'))+
  scale_x_continuous(breaks=1:18, labels=FGnames[rev(filter(red_dat, dat=='Australia')$groupk19),]$FEcomp)+
  theme_bw()+theme(legend.justification = c(1, 1), legend.position = c(1, 1), 
                   legend.background = element_rect(colour = "black"), 
                   axis.text.x = element_text(angle = 60, hjust = 1))+
  guides(fill=guide_legend(title="Thermal Affinity"))+
  scale_fill_manual(values = c('#2c7bb6', '#fdae61'))+
  ylab('# of species per FG')+xlab('Functional niche')

#### Tropical-substropical Biomass comparisons (log) ####

#sanity check to make sure per unit area biomass calc is correct

# If Logging: 0.01 chosen as min constant as min corr biomass val
# in Aus is 0.0016 and second min in Japan is 0.00098 but there are only 5

ggplot(filter(bio_jpn, FG %in% c(15, 10, 8, 2, 6, 12, 4, 1, 16)),
       aes(x = lat, y = log((tot_biom/totMsurv)+0.01), colour=ThermalAffinity2)) + 
  geom_point(shape=1)+geom_smooth(se=F)+facet_wrap(~FG, scales='free')+
  scale_colour_manual(values = c("#377eb8", "#ff7f00"))+
  theme_bw()+theme(legend.position = "none")+
  xlab('Latitude')+ylab('log standardised biomass')

ggplot(filter(bio_aus, FG %in% c(15, 10, 8, 2, 6, 12, 4, 1, 16)),
       aes(x = Lat, y = log((tot_biom/totMsurv)+0.01), colour=ThermalAffinity2)) + 
  geom_point(shape=1)+geom_smooth(se=F)+facet_wrap(~FG, scales='free')+
  scale_colour_manual(values = c("#377eb8", "#ff7f00"))+
  theme_bw()+theme(legend.position = "none")+scale_x_reverse()+
  xlab('Latitude')+ylab('log standardised biomass')
####
#### Tropical-substropical Biomass comparisons (4rt) ####

jpn_cor<-bio_jpn%>%filter(FG %in% c(15, 10, 8, 2, 6, 12, 4, 1, 16))%>%
  dplyr::select(FG, Site.trans.ID,ThermalAffinity2, cor_biom)%>%
  group_by(FG,Site.trans.ID)%>%tidyr::spread(ThermalAffinity2, cor_biom)%>%
  ungroup()%>%group_by(FG)%>%
  summarise(cor_p=cor.test((tropical^0.25), (subtropical^0.25), method = 'kendall')$p.value,
               cop_est=cor.test((tropical^0.25), (subtropical^0.25), method = 'kendall')$estimate)
                 
jpn_cor$txt=ifelse(jpn_cor$cor_p<0.001, paste0('tau=',round(jpn_cor$cop_est, 2), '***'),
                   ifelse(jpn_cor$cor_p<0.01, paste0('tau=',round(jpn_cor$cop_est, 2), '**'),
                          ifelse(jpn_cor$cor_p<0.05, paste0('tau=',round(jpn_cor$cop_est, 2), '*'),
                                 NA)))                     

aus_cor<-bio_aus%>%filter(FG %in% c(15, 10, 8, 2, 6, 12, 4, 1, 16))%>%
  dplyr::select(FG, Site.trans.ID,ThermalAffinity2, cor_biom)%>%
  group_by(FG,Site.trans.ID)%>%tidyr::spread(ThermalAffinity2, cor_biom)%>%
  ungroup()%>%group_by(FG)%>%
  summarise(cor_p=cor.test((tropical^0.25), (subtropical^0.25), method = 'kendall')$p.value,
            cop_est=cor.test((tropical^0.25), (subtropical^0.25), method = 'kendall')$estimate)

aus_cor$txt=ifelse(aus_cor$cor_p<0.001, paste0('tau=',round(aus_cor$cop_est, 2), '***'),
                   ifelse(aus_cor$cor_p<0.01, paste0('tau=',round(aus_cor$cop_est, 2), '**'),
                   ifelse(aus_cor$cor_p<0.05, paste0('tau=',round(aus_cor$cop_est, 2), '*'),
                   NA)))                     

sg_lat_spans<-data.frame(xmin=c(24.2, 26.2, 28.5, 31,   32.7,  33.38, 34.6), 
                         xmax=c(24.5, 28.4, 30.5, 31.6, 32.82, 33.5, 35 ))


bio_jpn$ThermalAffinity2<-factor(bio_jpn$ThermalAffinity2, levels=c('tropical', 'subtropical'))

p1<-ggplot(data=filter(bio_jpn, FG %in% c(15, 10, 8, 2, 6, 12, 4, 1, 16))) + 
  geom_rect(data=sg_lat_spans, aes(ymin=0, ymax=4, xmin=xmin, xmax=xmax), fill='grey', alpha=0.5)+
  geom_point(data=filter(bio_jpn, FG %in% c(15, 10, 8, 2, 6, 12, 4, 1, 16) & ((tot_biom/totMsurv)^0.25)<4),
             aes(x = lat, y = (tot_biom/totMsurv)^0.25, colour=ThermalAffinity2), shape=1)+
  geom_smooth(aes(x = lat, y = (tot_biom/totMsurv)^0.25, colour=ThermalAffinity2),se=F)+
  geom_label(data=jpn_cor, aes(x=28.5, y=3.5, label=txt))+
  facet_grid(FG~., scales='free_y')+
  theme_bw()+theme(legend.position = "none")+
  xlab('Latitude')+ylab('4rt-trans standardised biomass')+
  ggtitle('Japan')+
  theme(strip.background = element_blank(),strip.text.y = element_blank(),
        plot.title = element_text(hjust = 0.5))

sg_lat_spans<-data.frame(xmin=c(-23.4, -24.8, -26.61, -28.19, -29.9, -29.97), 
                           xmax=c(-24.1, -25.3, -26.98, -28.616, -30.96, -30.3))

bio_aus$ThermalAffinity2<-factor(bio_aus$ThermalAffinity2, levels=c('tropical', 'subtropical'))

p2<-ggplot(data=filter(bio_aus, FG %in% c(15, 10, 8, 2, 6, 12, 4, 1, 16))) + 
    geom_rect(data=sg_lat_spans, aes(ymin=0, ymax=4, xmin=xmin, xmax=xmax), fill='grey', alpha=0.5)+
    geom_point(data=filter(bio_aus, FG %in% c(15, 10, 8, 2, 6, 12, 4, 1, 16) & ((tot_biom/totMsurv)^0.25)<4),
               aes(x = Lat, y = (tot_biom/totMsurv)^0.25, colour=ThermalAffinity2), shape=1)+
    geom_smooth(aes(x = Lat, y = (tot_biom/totMsurv)^0.25, colour=ThermalAffinity2),se=F)+
  geom_label(data=aus_cor, aes(x=-26, y=3.5, label=txt))+  
  facet_grid(FG~., scales='free_y')+scale_x_reverse()+ylab(NULL)+
  ggtitle('Australia')+
    theme_bw()+theme(legend.position = "none")+
    xlab('Latitude')+theme(plot.title = element_text(hjust = 0.5))

#png('C:/coral_fish/outputs/biomass_thermal_plot.png',width = 8, height =12 , units ="in", res =600)
grid.arrange(p1, p2, ncol=2)
#dev.off()
#plots have outlier points removed (>4) but curves still fitted to full data
####
#### Tropical-substropical SPECIES RICHNESS comparisons (suppl) ####

# make same plot with spprich data while we're here

ggplot(filter(spr_jpn, FG %in% c(15, 10, 8, 2, 6, 12, 4, 1, 16)),
       aes(x = lat, y = qD, colour=ThermalAffinity2)) + 
  geom_point()+geom_smooth(se=F)+facet_wrap(~FG, scales='free')+
  scale_colour_manual(values = c("#377eb8", "#ff7f00"))+
  theme_bw()+theme(legend.position = "none")+
  xlab('Latitude')+ylab('Estimated species richness')

ggplot(filter(spr_aus, FG %in% c(15, 10, 8, 2, 6, 12, 4, 1, 16)),
       aes(x = Lat, y = qD, colour=ThermalAffinity2)) + 
  geom_point()+geom_smooth(se=F)+facet_wrap(~FG, scales='free')+
  scale_colour_manual(values = c("#377eb8", "#ff7f00"))+
  theme_bw()+theme(legend.position = "none")+scale_x_reverse()+
  xlab('Latitude')+ylab('Estimated species richness')

####

#### Setup tropical-only data for tropicalization analyses ####

# Calc standardisation
# filter out unwanted FGs for per FG objects but not for community total
jpn_trop_prop<-filter(bio_jpn, ThermalAffinity2=='tropical' & lat<25.5)%>%
  group_by(FG)%>%summarise(mean_biom=(mean(cor_biom^0.25, na.rm=T))^4)

aus_trop_prop<-filter(bio_aus, ThermalAffinity2=='tropical' & Lat> -24.5)%>%
  group_by(FG)%>%summarise(mean_biom=(mean(cor_biom^0.25, na.rm=T))^4)

jpn_trop_prop<-left_join(filter(bio_jpn, ThermalAffinity2=='tropical' & 
                                  FG %in% c(15, 10, 8, 2, 6, 12, 4, 1, 16)),
                         jpn_trop_prop, by='FG')

aus_trop_prop<-left_join(filter(bio_aus, ThermalAffinity2=='tropical' & 
                                  FG %in% c(15, 10, 8, 2, 6, 12, 4, 1, 16)),
                         aus_trop_prop, by='FG')

jpn_trop_comm<-filter(bio_jpn, ThermalAffinity2=='tropical')%>%
  group_by(lat, SiteID, Site.trans.ID)%>%summarise(cor_biom=sum(cor_biom)) 
jpn_trop_comm$mean_biom<-as.numeric(filter(jpn_trop_comm, lat<25.5)%>%ungroup()%>%
  summarise(mean_biom=(mean(cor_biom^0.25, na.rm=T))^4))

aus_trop_comm<-filter(bio_aus, ThermalAffinity2=='tropical')%>%
  group_by(Lat, Site, Site.trans.ID)%>%summarise(cor_biom=sum(cor_biom))
aus_trop_comm$mean_biom<-as.numeric(filter(aus_trop_comm, Lat> -24.5)%>%ungroup()%>%
  summarise(mean_biom=(mean(cor_biom^0.25, na.rm=T))^4))

jpn_trop_prop$FG<-factor(jpn_trop_prop$FG, levels=c(15, 10, 8, 2,6,12,4,1,16))
aus_trop_prop$FG<-factor(aus_trop_prop$FG, levels=c(15, 10, 8, 2,6,12,4,1,16))

# setup data for analyses
# add tropicalization_metric: >1 = more biomass at site compared to tropical site-group
# Using 4th root transformation to make data approx normality (still 0's though)
jpn_trop_comm$trop_met<-(jpn_trop_comm$cor_biom/jpn_trop_comm$mean_biom)^0.25
aus_trop_comm$trop_met<-(aus_trop_comm$cor_biom/aus_trop_comm$mean_biom)^0.25
jpn_trop_prop$trop_met<-(jpn_trop_prop$cor_biom/jpn_trop_prop$mean_biom)^0.25
aus_trop_prop$trop_met<-(aus_trop_prop$cor_biom/aus_trop_prop$mean_biom)^0.25

# create site-groups based on dendrogram clustering
jpn_trop_comm$FG<-'comm'
jpn_trop_prop$FG<-as.character(jpn_trop_prop$FG)
jpn_trop_tests<-rbind(data.frame(jpn_trop_comm), jpn_trop_prop[names(jpn_trop_comm)])
jpn_trop_tests$site.group<-'trop.base'
jpn_trop_tests[jpn_trop_tests$lat > 25 &  jpn_trop_tests$lat < 28.5,]$site.group<-'trop.island'
jpn_trop_tests[jpn_trop_tests$lat > 28.5 &  jpn_trop_tests$lat < 31,]$site.group<-'trans.island'
jpn_trop_tests[jpn_trop_tests$SiteID %in% c('JP28', 'JP29', 'JP30', 'JP31'),]$site.group<-'trans.inland'
jpn_trop_tests[jpn_trop_tests$SiteID %in% c('JP8', 'JP9', 'JP10', 'JP11', 'JP12'),]$site.group<-'trans.headld'
jpn_trop_tests[jpn_trop_tests$lat > 34,]$site.group<-'temp.headld'
table(jpn_trop_tests$SiteID, jpn_trop_tests$site.group)

aus_trop_comm$FG<-'comm'
aus_trop_prop$FG<-as.character(aus_trop_prop$FG)
aus_trop_tests<-rbind(data.frame(aus_trop_comm), aus_trop_prop[names(aus_trop_comm)])
aus_trop_tests$site.group<-'trop.base'
aus_trop_tests[aus_trop_tests$Lat > -25.6 &  aus_trop_tests$Lat < -24.5,]$site.group<-'trans.bay'
aus_trop_tests[aus_trop_tests$Lat > -28 &  aus_trop_tests$Lat < -25.6,]$site.group<-'trans.offshore'
aus_trop_tests[aus_trop_tests$Lat < -28,]$site.group<-'temp.offshore'
aus_trop_tests[aus_trop_tests$Site %in% c('Julian Rock False Trench', 'Julian Rock Nursery', 'Cook Island'),]$site.group<-'trans.temp'
aus_trop_tests[aus_trop_tests$Site %in% c('Muttonbird Island', 'Woolgoolga Reef', 'Woolgoolga Headland', 'North Rock'),]$site.group<-'temp.inshore'

table(aus_trop_tests$Site,aus_trop_tests$site.group)
####

#### Make all-FG tropicalization trends + comm ####

sg_lat_spans<-data.frame(xmin=c(24.2, 26.2, 28.5, 31,   32.7,  33.38, 34.6), 
                         xmax=c(24.5, 28.4, 30.5, 31.6, 32.82, 33.5, 35 ))

jpn_trend<-ggplot(jpn_trop_prop) + 
  geom_rect(data=sg_lat_spans, aes(ymin=0, ymax=2, xmin=xmin, xmax=xmax), fill='grey', alpha=0.5)+
  geom_smooth(aes(x = lat, y = (cor_biom/mean_biom)^0.25, colour=factor(FG)), se=F)+
  geom_smooth(data=jpn_trop_comm, aes(x = lat, y = (cor_biom/mean_biom)^0.25),
              se=F, colour='black', linetype='dashed')+
  theme_bw()+theme(legend.position = "none")+
  geom_hline(yintercept=0.05^0.25)+scale_x_continuous(breaks=24:35)+
  xlab('Latitude')+ylab('Proportion of tropical biomass')+
  scale_y_continuous(breaks=c(0, 0.05, 0.25, 0.5, 1, 2, 4, 20)^0.25, 
                     labels=c(0,0.05, 0.25, 0.5, 1, 2, 4, 20), minor_breaks = NULL)
  
sg_lat_spans<-data.frame(xmin=c(-23.4, -24.8, -26.61, -28.19, -29.9, -29.97), 
                         xmax=c(-24.1, -25.3, -26.98, -28.616, -30.96, -30.3))

aus_trend<-ggplot(aus_trop_prop) + 
  geom_rect(data=sg_lat_spans, aes(ymin=0, ymax=4^0.25, xmin=xmin, xmax=xmax), fill='grey', alpha=0.6)+
  geom_smooth(aes(x = Lat, y = (cor_biom/mean_biom)^0.25, colour=factor(FG)), se=F)+
  geom_smooth(data=aus_trop_comm, aes(x = Lat, y = (cor_biom/mean_biom)^0.25),
              se=F, colour='black', linetype='dashed')+
  theme_bw()+theme(legend.position = "none")+
  geom_hline(yintercept=0.05^0.25)+scale_x_reverse(breaks=-23:-31)+
  xlab('Latitude')+ylab('Proportion of tropical biomass')+
scale_y_continuous(breaks=c(0, 0.05, 0.25, 0.5, 1, 2, 4)^0.25, 
                   labels=c(0,0.05, 0.25, 0.5, 1, 2, 4), minor_breaks = NULL)
####


#### Visualise FG tropicalization using PCA (suppl) ####

#visualise using pca
jpn_pca<-jpn_trop_prop%>%group_by(FG, SiteID)%>%summarise(trop_met=mean(trop_met))%>%
         ungroup()%>% group_by(FG)%>%tidyr::spread(SiteID, trop_met)

aus_pca<-aus_trop_prop%>%group_by(FG, Site)%>%summarise(trop_met=mean(trop_met))%>%
  ungroup()%>% group_by(FG)%>%tidyr::spread(Site, trop_met)

# remove tropical base
jpn_pca<-jpn_pca[-c(2,13,22,25)]
aus_pca<-aus_pca[-c(12,13,14,23,24)]

jpnpc<-rda(jpn_pca[,c(2:length(jpn_pca)),], scale=T) # scaled so each transect equivelent
auspc<-rda(aus_pca[,c(2:length(aus_pca)),], scale=T)

jpnpc.dat<-as.data.frame(scores(jpnpc, choices=1:2, display='sites', scaling=1)) 
auspc.dat<-as.data.frame(scores(auspc, choices=1:2, display='sites', scaling=1)) 

jpnpc.dat$FG=factor(c('Upper-benthic Omnivores', 'Upper-benthic planktivores', 'Benthic planktivores',
               'Upper-benthic Herbivores', 'Benthic Herbivore/omnivores', 'Upper-benthic predators',
               'Demersal Predators', 'Benthic Predators', 'Corallivores'), 
               levels=c('Benthic Predators','Upper-benthic predators', 'Benthic Herbivore/omnivores',
              'Upper-benthic planktivores','Upper-benthic Herbivores','Demersal Predators',
              'Benthic planktivores','Upper-benthic Omnivores','Corallivores'))

auspc.dat$FG=factor(c('Upper-benthic Omnivores', 'Upper-benthic planktivores', 'Benthic planktivores',
               'Upper-benthic Herbivores', 'Benthic Herbivore/omnivores', 'Upper-benthic predators',
               'Demersal Predators', 'Benthic Predators', 'Corallivores'),
               levels=c('Benthic Predators','Upper-benthic predators', 'Benthic Herbivore/omnivores',
                        'Upper-benthic planktivores','Upper-benthic Herbivores','Demersal Predators',
                        'Benthic planktivores','Upper-benthic Omnivores','Corallivores'))
               

g<- ggplot()+
  geom_segment(data=NULL, aes(y=-Inf, x=0, yend=Inf, xend=0), linetype='dotted')+
  geom_segment(data=NULL, aes(y=0, x=-Inf, yend=0, xend=Inf), linetype='dotted')
  
g1<-g+geom_point(data=jpnpc.dat, aes(y=PC2, x=PC1, colour=FG),size=3)+
  geom_text_repel(data=jpnpc.dat, aes(y=PC2, x=PC1, colour=FG, label=FG),size=3)+
  theme_classic()+theme(legend.position = "none") 

eig<-eigenvals(jpnpc)
g1<- g1+scale_y_continuous(paste(names(eig[2]), sprintf('(%0.1f%% explained var.)', 100* eig[2]/sum(eig))))+
  scale_x_continuous(paste(names(eig[1]), sprintf('(%0.1f%% explained var.)', 100* eig[1]/sum(eig))))


g2<-g+geom_point(data=auspc.dat, aes(y=PC2, x=PC1, colour=FG),size=3)+
  geom_text_repel(data=auspc.dat, aes(y=PC2, x=PC1, colour=FG, label=FG),size=3)+
  theme_classic()+theme(legend.position = "none") 

eig<-eigenvals(auspc)
g2<- g2+scale_y_continuous(paste(names(eig[2]), sprintf('(%0.1f%% explained var.)', 100* eig[2]/sum(eig))))+
  scale_x_reverse(paste(names(eig[1]), sprintf('(%0.1f%% explained var.)', 100* eig[1]/sum(eig))))

grid.arrange(g1, g2)

# also using table 1 data

t1_dat<-read.csv('C:/coral_fish/outputs/paper_table_data.csv')
qplot(data=t1_dat, x=peak, y=edge, colour=FG, shape=regon)#nope
####



#### Japan FG tropicalization comparisons ####

# set comm as intercept
jpn_trop_tests$FG<-factor(jpn_trop_tests$FG, levels=c('comm', 15, 10, 8, 2,6,12,4,1,16))

# remove comm for diff between FG stat reporting
#jpn_trop_tests<-filter(jpn_trop_tests, FG!='comm')
#jpn_trop_tests$FG<-factor(jpn_trop_tests$FG, levels=c( 15, 10, 8, 2,6,12,4,1,16))

## trop.island

ggplot(data=filter(jpn_trop_tests, site.group=='trop.island'), aes(x=FG, y=trop_met))+
  geom_point(aes(colour=SiteID))+geom_boxplot(alpha=0.5) # remember boxplot = medians

m1<-lmer(trop_met~FG+(1|SiteID), data=filter(jpn_trop_tests, site.group=='trop.island'))
#resid_panel(m1)
boxplot(residuals(m1, type='pearson')~filter(jpn_trop_tests, site.group=='trop.island')$FG) #hmmm not good
summary(m1)

# try lme model
m1<-lme(trop_met~FG, random=~1|SiteID, data=filter(jpn_trop_tests, site.group=='trop.island'), 
        weights=varIdent(form=~1|FG))
boxplot(residuals(m1, type='pearson')~filter(jpn_trop_tests, site.group=='trop.island')$FG) # ok good
# pearson and normalized residuals incorporate the effect of the gls variance structure

em1<-emmeans(m1, specs='FG')
pairs(em1)
plot(em1, comparisons = TRUE)

trop_comps_out<-data.frame(site.group='trop.island',data.frame(em1), rbind(c(NA, NA), data.frame(pairs(em1))[1:9, c(1,6)]))

ggplot()+
  geom_jitter(data=filter(jpn_trop_tests, site.group=='trop.island'), aes(x=FG, y=trop_met), shape=1, alpha=0.5, width=0.2, colour='grey')+
  geom_hline(yintercept = trop_comps_out$emmean[1], linetype='dotted')+
  geom_errorbar(data=trop_comps_out, aes(x=FG, ymin=lower.CL, ymax=upper.CL))+
  geom_point(data=trop_comps_out, aes(x=FG, y=emmean, colour=ifelse(p.value>0.05| is.na(p.value), 'blue', 'red')), size=2)+
 theme_bw()+theme(legend.position = "none")

## trans.island

ggplot(data=filter(jpn_trop_tests, site.group=='trans.island'), aes(x=FG, y=trop_met))+
  geom_point(aes(colour=SiteID))+geom_boxplot(alpha=0.5) # remember boxplot = medians

m1<-lme(trop_met~FG, random=~1|SiteID, weights=varIdent(form=~1|FG), data=filter(jpn_trop_tests, site.group=='trans.island'))
#resid_panel(m1)
summary(m1)

em1<-emmeans(m1, specs='FG')
pairs(em1)
plot(em1, comparisons = TRUE)

trop_comps_out<-rbind(trop_comps_out,
                      data.frame(site.group='trans.island',data.frame(em1), rbind(c(NA, NA), data.frame(pairs(em1))[1:9, c(1,6)])))

ggplot()+
  geom_jitter(data=filter(jpn_trop_tests, site.group=='trans.island'), aes(x=FG, y=trop_met), shape=1, alpha=0.5, width=0.2, colour='grey')+
  geom_hline(yintercept =trop_comps_out[trop_comps_out$site.group=='trans.island',]$response[1], linetype='dotted')+
  geom_errorbar(data=filter(trop_comps_out, site.group=='trans.island'), aes(x=FG, ymin=lower.CL, ymax=upper.CL))+
  geom_point(data=filter(trop_comps_out, site.group=='trans.island'), aes(x=FG, y=emmean, colour=ifelse(p.value>0.05| is.na(p.value), 'blue', 'red')), size=2)+
  theme_bw()+theme(legend.position = "none")

## trans.inland

ggplot(data=filter(jpn_trop_tests, site.group=='trans.inland'), aes(x=FG, y=trop_met))+
  geom_point(aes(colour=SiteID))+geom_boxplot(alpha=0.5) # remember boxplot = medians

m1<-lme(trop_met~FG, random=~1|SiteID, weights=varIdent(form=~1|FG), data=filter(jpn_trop_tests, site.group=='trans.inland'))
#resid_panel(m1)
summary(m1)


em1<-emmeans(m1, specs='FG')
pairs(em1)
plot(em1, comparisons = TRUE)

trop_comps_out<-rbind(trop_comps_out,
                      data.frame(site.group='trans.inland',data.frame(em1), rbind(c(NA, NA), data.frame(pairs(em1))[1:9, c(1,6)])))

ggplot()+
  geom_jitter(data=filter(jpn_trop_tests, site.group=='trans.inland'), aes(x=FG, y=trop_met), shape=1, alpha=0.5, width=0.2, colour='grey')+
  geom_hline(yintercept =trop_comps_out[trop_comps_out$site.group=='trans.inland',]$response[1], linetype='dotted')+
  geom_errorbar(data=filter(trop_comps_out, site.group=='trans.inland'), aes(x=FG, ymin=lower.CL, ymax=upper.CL))+
  geom_point(data=filter(trop_comps_out, site.group=='trans.inland'), aes(x=FG, y=emmean, colour=ifelse(p.value>0.05| is.na(p.value), 'blue', 'red')), size=2)+
  theme_bw()+theme(legend.position = "none")


## trans.headld

ggplot(data=filter(jpn_trop_tests, site.group=='trans.headld'), aes(x=FG, y=trop_met))+
  geom_point(aes(colour=SiteID))+geom_boxplot(alpha=0.5) # remember boxplot = medians

m1<-lme(trop_met~FG, random=~1|SiteID, weights=varIdent(form=~1|FG),
        data=filter(jpn_trop_tests, site.group=='trans.headld'),
        control=lmeControl(opt = 'optim')) # different optimizer used to make run
#resid_panel(m1)
summary(m1)


em1<-emmeans(m1, specs='FG')
pairs(em1)
plot(em1, comparisons = TRUE)

trop_comps_out<-rbind(trop_comps_out,
                      data.frame(site.group='trans.headld',data.frame(em1), rbind(c(NA, NA),
                                data.frame(pairs(em1))[1:9, c(1,6)])))

ggplot()+
  geom_jitter(data=filter(jpn_trop_tests, site.group=='trans.headld'), aes(x=FG, y=trop_met), shape=1, alpha=0.5, width=0.2, colour='grey')+
  geom_hline(yintercept =trop_comps_out[trop_comps_out$site.group=='trans.headld',]$response[1], linetype='dotted')+
  geom_errorbar(data=filter(trop_comps_out, site.group=='trans.headld'), aes(x=FG, ymin=lower.CL, ymax=upper.CL))+
  geom_point(data=filter(trop_comps_out, site.group=='trans.headld'), aes(x=FG, y=emmean, colour=ifelse(p.value>0.05| is.na(p.value), 'blue', 'red')), size=2)+
  theme_bw()+theme(legend.position = "none")

## temp.headld

ggplot(data=filter(jpn_trop_tests, site.group=='temp.headld'), aes(x=FG, y=trop_met))+
  geom_point(aes(colour=SiteID))+geom_boxplot(alpha=0.5) # remember boxplot = medians
# remove FGs with only 0s - these will be sig diff 
zer_fgs<-filter(jpn_trop_tests, site.group=='temp.headld')%>%group_by(FG)%>%summarise(st=sum(trop_met))%>%filter(., st==0)%>%.$FG

m1<-lme(trop_met~FG, random=~1|SiteID, weights=varIdent(form=~1|FG), data=filter(jpn_trop_tests, site.group=='temp.headld' & !FG%in% zer_fgs))
#resid_panel(m1)
summary(m1)

em1<-emmeans(m1, specs='FG')
pairs(em1)
plot(em1, comparisons = TRUE)

trop_comps_out<-rbind(trop_comps_out,
                      data.frame(site.group='temp.headld',data.frame(em1), rbind(c(NA, NA), data.frame(pairs(em1))[1:(9-length(zer_fgs)), c(1,6)])))

ggplot()+
  geom_jitter(data=filter(jpn_trop_tests, site.group=='temp.headld'), aes(x=FG, y=trop_met), shape=1, alpha=0.5, width=0.2, colour='grey')+
  geom_hline(yintercept =trop_comps_out[trop_comps_out$site.group=='temp.headld',]$response[1], linetype='dotted')+
  geom_errorbar(data=filter(trop_comps_out, site.group=='temp.headld'), aes(x=FG, ymin=lower.CL, ymax=upper.CL))+
  geom_point(data=filter(trop_comps_out, site.group=='temp.headld'), aes(x=FG, y=emmean, colour=ifelse(p.value>0.05| is.na(p.value), 'blue', 'red')), size=2)+
  theme_bw()+theme(legend.position = "none")

#all plot
trop_comps_out<-rbind(trop_comps_out,
                      data.frame(site.group='temp.headld', FG=c(6,12,1,16), emmean=0, SE=0,
                      df=0, lower.CL=0, upper.CL=0, contrast=NA, p.value=0))

trop_comps_out$site.group<-factor(trop_comps_out$site.group,
              levels=c('trop.base', 'trop.island', 'trans.island', 'trans.inland',
                       'trans.headld', 'temp.headld'))
 

jpn_mods<-ggplot()+
  geom_hline(yintercept=0.05^0.25, linetype='dotted')+
  geom_errorbar(data=trop_comps_out, aes(x=FG, ymin=lower.CL, ymax=upper.CL))+
  geom_point(data=trop_comps_out, aes(x=FG, y=emmean, size=ifelse(is.na(p.value), 1, ifelse(p.value>0.05, 1, 2))), shape=1)+
  geom_point(data=trop_comps_out, aes(x=FG, y=emmean, colour=factor(FG)), size=2)+
  theme_bw()+theme(legend.position = "none")+facet_wrap(~site.group, nrow=1)+
  scale_colour_manual(values = c("black", "#F8766D", "#D39200" ,"#93AA00", "#00BA38",
 "#00C19F", "#00B9E3", "#619CFF", "#DB72FB", "#FF61C3"))+
  ylab('Proportion of tropical biomass')+
  scale_y_continuous(breaks=c(0, 0.05, 0.25, 0.5, 1, 2, 4, 20)^0.25, 
   labels=c(0,0.05, 0.25, 0.5, 1, 2, 4, 20), minor_breaks = NULL)


#### Australia FG tropicalization comparisons #### 

# drop some FGs and set comm as intercept
aus_trop_tests$FG<-factor(aus_trop_tests$FG, levels=c('comm', 15, 10, 8, 2,6,12,4,1,16))

# remove comm for diff between FG stat reporting
#aus_trop_tests<-filter(aus_trop_tests, FG!='comm')
#aus_trop_tests$FG<-factor(aus_trop_tests$FG, levels=c( 15, 10, 8, 2,6,12,4,1,16))


## trans.bay

ggplot(data=filter(aus_trop_tests, site.group=='trans.bay'), aes(x=FG, y=trop_met))+
  geom_point(aes(colour=Site))+geom_boxplot(alpha=0.5) # remember boxplot = medians

# drop 2 big outliers, actually don't
m1<-lme(trop_met~FG, random=~1|Site, weights=varIdent(form=~1|FG), data=filter(aus_trop_tests,
                         site.group=='trans.bay'))
#resid_panel(m1)
summary(m1)

em1<-emmeans(m1, specs='FG')
pairs(em1)

trop_comps_out_aus<-data.frame(site.group='trans.bay',data.frame(em1), rbind(c(NA, NA), data.frame(pairs(em1))[1:9, c(1,6)]))

ggplot()+
  geom_jitter(data=filter(aus_trop_tests, site.group=='trans.bay'), aes(x=FG, y=trop_met), shape=1, alpha=0.5, width=0.2, colour='grey')+
  geom_hline(yintercept =trop_comps_out_aus[trop_comps_out_aus$site.group=='trans.bay',]$response[1], linetype='dotted')+
  geom_errorbar(data=filter(trop_comps_out_aus, site.group=='trans.bay'), aes(x=FG, ymin=lower.CL, ymax=upper.CL))+
  geom_point(data=filter(trop_comps_out_aus, site.group=='trans.bay'), aes(x=FG, y=emmean, colour=ifelse(p.value>0.05| is.na(p.value), 'blue', 'red')), size=2)+
  theme_bw()+theme(legend.position = "none") # stupid massive outlier drags up 12, makes it look ns but it is sig diff


## trans.offshore

ggplot(data=filter(aus_trop_tests, site.group=='trans.offshore'), aes(x=FG, y=trop_met))+
  geom_point(aes(colour=Site))+geom_boxplot(alpha=0.5) # remember boxplot = medians

m1<-lme(trop_met~FG, random=~1|Site, weights=varIdent(form=~1|FG), data=filter(aus_trop_tests, site.group=='trans.offshore'))
#resid_panel(m1)
summary(m1)


em1<-emmeans(m1, specs='FG')
#pairs(em1)
plot(em1, comparisons = TRUE)

trop_comps_out_aus<-rbind(trop_comps_out_aus,
                      data.frame(site.group='trans.offshore',data.frame(em1), rbind(c(NA, NA),
                                                                                  data.frame(pairs(em1))[1:9, c(1,6)])))
ggplot()+
  geom_jitter(data=filter(aus_trop_tests, site.group=='trans.offshore'), aes(x=FG, y=trop_met), shape=1, alpha=0.5, width=0.2, colour='grey')+
  geom_hline(yintercept =trop_comps_out_aus[trop_comps_out_aus$site.group=='trans.offshore',]$response[1], linetype='dotted')+
  geom_errorbar(data=filter(trop_comps_out_aus, site.group=='trans.offshore'), aes(x=FG, ymin=lower.CL, ymax=upper.CL))+
  geom_point(data=filter(trop_comps_out_aus, site.group=='trans.offshore'), aes(x=FG, y=emmean, colour=ifelse(p.value>0.05| is.na(p.value), 'blue', 'red')), size=2)+
  theme_bw()+theme(legend.position = "none")


## trans.temp

ggplot(data=filter(aus_trop_tests, site.group=='trans.temp'), aes(x=FG, y=trop_met))+
  geom_point(aes(colour=Site))+geom_boxplot(alpha=0.5) # remember boxplot = medians

m1<-lme(trop_met~FG, random=~1|Site, weights=varIdent(form=~1|FG), data=filter(aus_trop_tests, site.group=='trans.temp'))
#resid_panel(m1)
summary(m1)

em1<-emmeans(m1, specs='FG')
#pairs(em1)
plot(em1, comparisons = TRUE)

trop_comps_out_aus<-rbind(trop_comps_out_aus,
                          data.frame(site.group='trans.temp',data.frame(em1), rbind(c(NA, NA),
                                                                                      data.frame(pairs(em1))[1:9, c(1,6)])))
ggplot()+
  geom_jitter(data=filter(aus_trop_tests, site.group=='trans.temp'), aes(x=FG, y=trop_met), shape=1, alpha=0.5, width=0.2, colour='grey')+
  geom_hline(yintercept =trop_comps_out_aus[trop_comps_out_aus$site.group=='trans.temp',]$response[1], linetype='dotted')+
  geom_errorbar(data=filter(trop_comps_out_aus, site.group=='trans.temp'), aes(x=FG, ymin=lower.CL, ymax=upper.CL))+
  geom_point(data=filter(trop_comps_out_aus, site.group=='trans.temp'), aes(x=FG, y=emmean, colour=ifelse(p.value>0.05| is.na(p.value), 'blue', 'red')), size=2)+
  theme_bw()+theme(legend.position = "none")

## temp.offshore

ggplot(data=filter(aus_trop_tests, site.group=='temp.offshore'), aes(x=FG, y=trop_met))+
  geom_point(aes(colour=Site))+geom_boxplot(alpha=0.5) # remember boxplot = medians

m1<-lme(trop_met~FG, random=~1|Site, weights=varIdent(form=~1|FG), data=filter(aus_trop_tests, 
                site.group=='temp.offshore' & trop_met<4)) # remove massive outlier
#resid_panel(m1)
summary(m1)


em1<-emmeans(m1, specs='FG')
#pairs(em1)
plot(em1, comparisons = TRUE)

trop_comps_out_aus<-rbind(trop_comps_out_aus,
                      data.frame(site.group='temp.offshore',data.frame(em1), rbind(c(NA, NA), data.frame(pairs(em1))[1:9, c(1,6)])))

ggplot()+
  geom_jitter(data=filter(aus_trop_tests, site.group=='temp.offshore'), aes(x=FG, y=trop_met), shape=1, alpha=0.5, width=0.2, colour='grey')+
  geom_hline(yintercept =trop_comps_out_aus[trop_comps_out_aus$site.group=='temp.offshore',]$response[1], linetype='dotted')+
  geom_errorbar(data=filter(trop_comps_out_aus, site.group=='temp.offshore'), aes(x=FG, ymin=lower.CL, ymax=upper.CL))+
  geom_point(data=filter(trop_comps_out_aus, site.group=='temp.offshore'), aes(x=FG, y=emmean, colour=ifelse(p.value>0.05| is.na(p.value), 'blue', 'red')), size=2)+
  theme_bw()+theme(legend.position = "none")


## temp.inshore

ggplot(data=filter(aus_trop_tests, site.group=='temp.inshore'), aes(x=FG, y=trop_met))+
  geom_point(aes(colour=Site))+geom_boxplot(alpha=0.5) # remember boxplot = medians
# remove FGs with only 0s - these will be sig diff 
zer_fgs<-filter(aus_trop_tests, site.group=='temp.inshore')%>%group_by(FG)%>%summarise(st=sum(trop_met))%>%filter(., st==0)%>%.$FG

m1<-lme(trop_met~FG, random=~1|Site, weights=varIdent(form=~1|FG), data=filter(aus_trop_tests, site.group=='temp.inshore' & !FG%in% zer_fgs))
#resid_panel(m1)
summary(m1)


em1<-emmeans(m1, specs='FG')
#pairs(em1)
#plot(em1, comparisons = TRUE)

trop_comps_out_aus<-rbind(trop_comps_out_aus,
                      data.frame(site.group='temp.inshore',data.frame(em1), rbind(c(NA, NA), data.frame(pairs(em1))[1:(9-length(zer_fgs)), c(1,6)])))

ggplot()+
  geom_jitter(data=filter(aus_trop_tests, site.group=='temp.inshore'), aes(x=FG, y=trop_met), shape=1, alpha=0.5, width=0.2, colour='grey')+
  geom_hline(yintercept =trop_comps_out_aus[trop_comps_out_aus$site.group=='temp.inshore',]$response[1], linetype='dotted')+
  geom_errorbar(data=filter(trop_comps_out_aus, site.group=='temp.inshore'), aes(x=FG, ymin=lower.CL, ymax=upper.CL))+
  geom_point(data=filter(trop_comps_out_aus, site.group=='temp.inshore'), aes(x=FG, y=emmean, colour=ifelse(p.value>0.05| is.na(p.value), 'blue', 'red')), size=2)+
  theme_bw()+theme(legend.position = "none")

#all plot
trop_comps_out_aus<-rbind(trop_comps_out_aus,
                      data.frame(site.group='temp.inshore', FG=c(6,12,16), emmean=0, SE=0,
                                 df=0, lower.CL=0, upper.CL=0, contrast=NA, p.value=0))


trop_comps_out_aus$site.group<-factor(trop_comps_out_aus$site.group, levels=c('trop.base', 'trans.bay', 'trans.offshore',
                                                                              'trans.temp',  'temp.offshore', 'temp.inshore'))

aus_mods<-ggplot()+
  geom_hline(yintercept=0.05^0.25, linetype='dotted')+
  geom_errorbar(data=trop_comps_out_aus, aes(x=FG, ymin=lower.CL, ymax=upper.CL))+
  geom_point(data=trop_comps_out_aus, aes(x=FG, y=emmean, size=ifelse(is.na(p.value), 1, ifelse(p.value>0.05, 1, 2))), shape=1)+
  geom_point(data=trop_comps_out_aus, aes(x=FG, y=emmean, colour=factor(FG)), size=2)+
  theme_bw()+theme(legend.position = "none")+facet_wrap(~site.group, nrow=1)+
  scale_colour_manual(values = c("black", "#F8766D", "#D39200" ,"#93AA00", "#00BA38",
  "#00C19F", "#00B9E3", "#619CFF", "#DB72FB", "#FF61C3"))+
  ylab('Proportion of tropical biomass')+
  scale_y_continuous(breaks=c(0, 0.05, 0.25, 0.5, 1, 2, 4, 20)^0.25, 
                     labels=c(0,0.05, 0.25, 0.5, 1, 2, 4, 20), minor_breaks = NULL)

# write out results
#write.csv(rbind(trop_comps_out, trop_comps_out_aus), 'C:/coral_fish/outputs/sitegroup_FG_tropicalization.csv', quote=F, row.names=F)


#### create tropicalization mega plot and write out ####

#png('C:/coral_fish/outputs/tropicalization_4plot.png',width = 12, height =12 , units ="in", res =600)

#grid.arrange(jpn_trend, jpn_mods, aus_trend, aus_mods, nrow=4)
#dev.off()
####


#### Look for sig diff between site.groups from comm level tropicalization ####

jpn_trop_tests$site.group<-factor(jpn_trop_tests$site.group,
                                  levels=c('trop.base', 'trop.island', 'trans.island', 'trans.inland',
                                           'trans.headld', 'temp.headld'))

ggplot(data=filter(jpn_trop_tests, FG=='comm'), aes(x=site.group, y=trop_met))+
  geom_point(aes(colour=SiteID))+geom_boxplot(alpha=0.5) # remember boxplot = medians

m1<-lme(trop_met~site.group, random=~1|SiteID, weights=varIdent(form=~1|site.group),
        data=filter(jpn_trop_tests, FG=='comm'))
resid_panel(m1)
summary(m1)
anova(m1)

em1<-emmeans(m1, specs='site.group')
pairs(em1)
plot(em1, comparisons = TRUE)

aus_trop_tests$site.group<-factor(aus_trop_tests$site.group, levels=c('trop.base', 'trans.bay', 'trans.offshore',
                                                                              'trans.temp',  'temp.offshore', 'temp.inshore'))
ggplot(data=filter(aus_trop_tests, FG=='comm'), aes(x=site.group, y=trop_met))+
  geom_point(aes(colour=Site))+geom_boxplot(alpha=0.5) # remember boxplot = medians

m1<-lme(trop_met~site.group, random=~1|Site, weights=varIdent(form=~1|site.group),
        data=filter(aus_trop_tests, FG=='comm'))
resid_panel(m1)
summary(m1)
anova(m1)

em1<-emmeans(m1, specs='site.group')
pairs(em1)
plot(em1, comparisons = TRUE)




#### stuart-smith FG tropicalization test ####

# compare tropicalization mean FG thermal midpoint data (stuart-smith)
therm_mid<-readxl::read_xlsx('C:/coral_fish/sourced_data/stuart_smith_thermal_midpoints/Thermal niche midpoints.xlsx')
tanika_match<-read.csv('C:/coral_fish/sourced_data/stuart_smith_thermal_midpoints/species_match_tanika.csv')
# bodge loop to fix 13 species names
for(i in tanika_match[tanika_match$Name.in.st.data!='',]$Name.in.st.data){
therm_mid[therm_mid$SPECIES_NAME==i,]$SPECIES_NAME<-as.character(tanika_match[tanika_match$Name.in.st.data==i,]$Not.matched)}
  
dat_ss<-left_join(dat, therm_mid, by=c('Species'='SPECIES_NAME'))
which(is.na(dat_ss$`95th SSTmax`)) # some naming mis-matches/ missing sp ~80
table(dat_ss$ThermalAffinity2, dat_ss$`Temp-Trop (23cutoff)`) # some differences

# check for coverage of FGs and compare within each region
# filter to FGs we're interested in
dat_ss<-filter(dat_ss, groupk19 %in% c(15, 10, 8, 2, 6, 12, 4, 1, 16))
dat_ss$groupk19<-factor(dat_ss$groupk19)
dat_ss$sst95<-dat_ss$`95th SSTmax`

#JPN
dat_ss%>%filter(JPN_sp>0 & ThermalAffinity2=='tropical')%>%group_by(groupk19)%>%
  summarise(n_sp=n(), n_sp_conf=length(confidence[which(confidence>1)])) # pretty good

conf_tm_jpn<-dat_ss%>%filter(JPN_sp>0 & ThermalAffinity2=='tropical' & confidence>1)

qplot(data=conf_tm_jpn, x=groupk19, y=sst95, geom='boxplot') # looks good but need closer 

m1<-lm(sst95~groupk19, data=conf_tm_jpn[-80,]) #remove outlier
par(mfrow=c(2,2));plot(m1)
summary(m1)
anova(m1) # ns
pairs(emmeans(m1, 'groupk19'))
boxplot(resid(m1, type='pearson')~factor(conf_tm_jpn$groupk19))

ggplot(data=data.frame(emmeans(m1, 'groupk19')), aes(x=factor(groupk19), y=emmean))+geom_pointrange(aes(ymin=lower.CL, ymax=upper.CL))+xlab('Functional Group')+ylab('Thermal Midpoint')

#try with GLS to be sure
m2<-gls(sst95~factor(groupk19), data=conf_tm_jpn,
        weights=varIdent(form=~1|groupk19))
boxplot(resid(m2, type='pearson')~factor(conf_tm_jpn$groupk19))
# no better, remember these groups have different n() so standard 
# error different anyway
m1<-ggplot(data=data.frame(emmeans(m1, 'groupk19')), aes(x=factor(groupk19), y=emmean))+
  geom_pointrange(aes(ymin=lower.CL, ymax=upper.CL))+xlab('Functional Group')+
  ylab('Thermal Midpoint (°C)')+scale_y_continuous(limits=c(30.5, 31.6),
  breaks=c(30.5,30.75, 31, 31.25, 31.5, 31.5))+labs(title ='Japan')

#AUS

dat_ss%>%filter(AUS_sp>0 & ThermalAffinity2=='tropical')%>%group_by(groupk19)%>%
  summarise(n_sp=n(), n_sp_conf=length(confidence[which(confidence>1)])) # pretty good

conf_tm_aus<-dat_ss%>%filter(AUS_sp>0 & ThermalAffinity2=='tropical' & confidence>1)

qplot(data=conf_tm_aus, x=groupk19, y=sst95, geom='boxplot') # looks good but need closer 

m1<-lm(sst95~groupk19, data=conf_tm_aus) #remove outlier
par(mfrow=c(2,2));plot(m1)
summary(m1)
anova(m1) # ns
pairs(emmeans(m1, 'groupk19'))
boxplot(resid(m1, type='pearson')~factor(conf_tm_aus$groupk19))

m2<-ggplot(data=data.frame(emmeans(m1, 'groupk19')), aes(x=factor(groupk19), y=emmean))+
  geom_pointrange(aes(ymin=lower.CL, ymax=upper.CL))+xlab('Functional Group')+
  ylab('Thermal Midpoint (°C)')+scale_y_continuous(limits=c(30.5, 31.6),
  breaks=c(30.5,30.75, 31, 31.25, 31.5, 31.5))+labs(title ='Australia')

#png('C:/coral_fish/outputs/themal_midpoint_suppl.png',width =8, height =4 , units ="in", res =600)
grid.arrange(m1, m2, ncol=2)
dev.off()


#### calculate functional distance/overlap between tropical invaders and

#### functional niche overlap, tropical vs higher latitude residents ####

# recreate distance matrix from clustering

distlog<-daisy(dat[,c("BodySize","Diet",  "Position", "Aggregation", 'DepthRange')],
               metric = "gower",stand = FALSE, type = list(logratio = c(1,5)))


func_dudi<-dudi.pco(d = cailliez(distlog, print=TRUE, cor.zero = F), scannf = FALSE, nf = 4)
#func_dudi<-dudi.pco(d = (distlog+3.04421), scannf = FALSE, nf =4)

screeplot(func_dudi)
hist(distlog$eig)

func_pco<-data.frame(func_dudi$li,dat)

# see how it looks

ppp<- ggplot()+
  geom_hline(yintercept=0, linetype="dotted") + 
  geom_vline(xintercept=0,  linetype="dotted")

ppp+geom_point(data=func_pco, aes(x=A1, y=A2))

# plot with variable contribution
#https://www.researchgate.net/post/how_can_i_produce_a_PCoA_biplot_using_R

efit <- envfit(func_dudi, dat[c("BodySize","Diet",  "Position", "Aggregation", 'DepthRange')], na.rm=T)
varibs<-data.frame(rbind(efit$vectors$arrows, efit$factors$centroids))
varibs$predictors=row.names(varibs)

ppp+geom_segment(data=varibs, aes(y=0, x=0, xend=A1, yend=A2),
                 arrow=arrow(length=unit(0.3,'lines')))+
  geom_text(data=varibs, aes(x=A1, y=A2, label=predictors))+
  
  scale_y_continuous(paste('PC2', sprintf('(%0.1f%% explained var.)',
                                          100* func_dudi$eig[2]/sum(func_dudi$eig[func_dudi$eig>0.007]))))+
  scale_x_continuous(paste('PC1', sprintf('(%0.1f%% explained var.)',
                                          100* func_dudi$eig[1]/sum(func_dudi$eig[func_dudi$eig>0.007]))))

# make plots per site-groups

jpn_sp_site$site.group<-'trop.base'
jpn_sp_site[jpn_sp_site$lat > 25 &  jpn_sp_site$lat < 28.5,]$site.group<-'trop.island'
jpn_sp_site[jpn_sp_site$lat > 28.5 &  jpn_sp_site$lat < 31,]$site.group<-'trans.island'
jpn_sp_site[jpn_sp_site$SiteID %in% c('JP28', 'JP29', 'JP30', 'JP31'),]$site.group<-'trans.inland'
jpn_sp_site[jpn_sp_site$SiteID %in% c('JP8', 'JP9', 'JP10', 'JP11', 'JP12'),]$site.group<-'trans.headld'
jpn_sp_site[jpn_sp_site$lat > 34,]$site.group<-'temp.headld'
table(jpn_sp_site$SiteID, jpn_sp_site$site.group)

aus_sp_site$site.group<-'trop.base'
aus_sp_site[aus_sp_site$Lat > -25.6 &  aus_sp_site$Lat < -24.5,]$site.group<-'trans.bay'
aus_sp_site[aus_sp_site$Lat > -28 &  aus_sp_site$Lat < -25.6,]$site.group<-'trans.offshore'
aus_sp_site[aus_sp_site$Lat < -28,]$site.group<-'temp.offshore'
aus_sp_site[aus_sp_site$Site %in% c('Julian Rock False Trench', 'Julian Rock Nursery', 'Cook Island'),]$site.group<-'trans.temp'
aus_sp_site[aus_sp_site$Site %in% c('Muttonbird Island', 'Woolgoolga Reef', 'Woolgoolga Headland', 'North Rock'),]$site.group<-'temp.inshore'

# join with PCoA axis data

# leave biomass summed per site

jpn_sp_site_pco<-left_join(jpn_sp_site, func_pco[,1:5], by=c('SpeciesFish'='Species'))
aus_sp_site_pco<-left_join(aus_sp_site, func_pco[,1:5], by=c('Fish'='Species'))

names(jpn_sp_site_pco)[1]<-'Site'
names(jpn_sp_site_pco)[2]<-'Lat'

funcOvl<-function(pcodat=mydat, FGz=c(10, 15), bywhat='site', mkern=F)
{  pcodat<-as.data.frame(pcodat)
  
  if(bywhat=='site'){
  pcodat$ID<-paste(pcodat$Site, pcodat$FG, pcodat$ThermalAffinity2, sep='@')}
  if(bywhat=='site.group'){
    pcodat$ID<-paste(pcodat$site.group, pcodat$FG, pcodat$ThermalAffinity2, sep='@')}
  
  # loop to find IDs with < 5 data points and replicate until n=5 so kernels can be made
  for(i in data.frame(table(pcodat$ID))[which(data.frame(table(pcodat$ID))$Freq<5),]$Var1)
  {
    pcodat<-rbind(pcodat, 
                  do.call(rbind,replicate(4, pcodat[pcodat$ID==i,], simplify=F)))
  }
  
  site_ovl_comp<-NULL
  site_ovl_filt<-NULL
  site_ovl_kernz<-NULL
  for(j in FGz)
  {
    pcFG<-pcodat[pcodat$FG==j,]
    
    if(bywhat=='site'){
    temp_base<-pcFG[pcFG$site.group=='trop.base',]
    temp_base$ID<-paste(temp_base$site.group, temp_base$FG, temp_base$ThermalAffinity2, sep='@')
    pcFG<-rbind(pcFG, temp_base)}  
    
    spdf<-SpatialPointsDataFrame(coords=pcFG[,c(8,9)],  data=data.frame(ID=pcFG$ID))
    
    spdf$ID<-factor(spdf$ID)
    
    KDE.Surface <- kernelUD(spdf,same4all = T, h=0.075, grid=500)
    
    if(mkern==T){
      KDE.99_site <- getverticeshr(KDE.Surface, percent = 99)
      KDE.99_site<-st_as_sf(KDE.99_site)
    temp<- do.call(rbind, lapply(strsplit(as.character(KDE.99_site$id), '@'),
    function(x) data.frame(Site=x[1], FG=x[2], ThermalAffinity2=x[3])))
    KDE.99_site$Site<-temp$Site
    KDE.99_site$FG<-temp$FG
    KDE.99_site$ThermalAffinity2<-temp$ThermalAffinity2
    }else{KDE.99_site<-NULL}

  ovl<-kerneloverlaphr(KDE.Surface, percent=99, meth="HR")
  ovl<-as.data.frame.table(ovl)
  ovl<-ovl[-which(ovl$Var1==ovl$Var2),]
  ovl<-cbind(ovl,do.call(rbind, lapply(strsplit(as.character(ovl$Var1), '@'),
                                       function(x) data.frame(V1Site=x[1], V1FG=x[2], V1ThermalAffinity2=x[3]))))
  ovl<-cbind(ovl,do.call(rbind, lapply(strsplit(as.character(ovl$Var2), '@'),
                                       function(x) data.frame(V2Site=x[1], V2FG=x[2], V2ThermalAffinity2=x[3]))))
  ovl<-ovl[ovl$V1FG==ovl$V2FG,]
  
  ovl_comp<-ovl[ovl$V1Site==ovl$V2Site,]
  ovl_comp<-ovl_comp%>%group_by(V1Site, V1FG)%>%summarise_all(first)%>%as.data.frame() # prop of tropical FG space covered by subtropical FG
  names(ovl_comp)[names(ovl_comp)=="V1Site"]<-'Site'
  names(ovl_comp)[names(ovl_comp)=="V1FG"]<-'FG'
  if(bywhat=='site'){ovl_comp<-ovl_comp[ovl_comp$Site!='trop.base',]}
  
  ovl_filt<-ovl[ovl$V1ThermalAffinity2==ovl$V2ThermalAffinity2,]
  ovl_filt<-ovl_filt[ovl_filt$V1ThermalAffinity2=='tropical',]
  ovl_filt<-ovl_filt[ovl_filt$V1Site=='trop.base',]
  names(ovl_filt)[names(ovl_filt)=="V2Site"]<-'Site'
  names(ovl_filt)[names(ovl_filt)=="V1FG"]<-'FG'
  
  site_ovl_comp<-rbind(site_ovl_comp, ovl_comp)
  site_ovl_filt<-rbind(site_ovl_filt, ovl_filt)
  site_ovl_kernz<-rbind(site_ovl_kernz, KDE.99_site)
  print(j)
    }
    return(list(site_ovl_comp, site_ovl_filt, site_ovl_kernz))
} # function end #

# run function per site for lat plots with tropicalization metric
ovl_site_aus<-funcOvl(pcodat=aus_sp_site_pco, FGz=c(15, 10, 8, 2,6, 12, 4,1, 16), bywhat='site', mkern=F)
ovl_site_jpn<-funcOvl(pcodat=jpn_sp_site_pco, FGz=c(15, 10, 8, 2,6, 12, 4,1, 16), bywhat='site', mkern=F)

# run function per site.group for kernel visualisation in FG space
ovl_site.group_aus<-funcOvl(pcodat=aus_sp_site_pco, FGz=c(15, 10, 8, 2,6, 12, 4,1, 16), bywhat='site.group', mkern=T)
ovl_site.group_jpn<-funcOvl(pcodat=jpn_sp_site_pco, FGz=c(15, 10, 8, 2,6, 12, 4,1, 16), bywhat='site.group', mkern=T)


#### Australia overlap by site ####

# trial plot for site level results
# fill overlap with 0's where no overlap
ovl_site_aus[[1]]<-ovl_site_aus[[1]]%>%tidyr::complete(Site, FG, fill=list(Freq=0))
ovl_site_aus[[2]]<-ovl_site_aus[[2]]%>%tidyr::complete(Site, FG, fill=list(Freq=0))

aus_ovl_comp<-left_join(ovl_site_aus[[1]],aus_sp_site%>%group_by(Site)%>%
                          summarise_all(first)%>%dplyr::select(Site, Lat), by='Site')
aus_ovl_filt<-left_join(ovl_site_aus[[2]],aus_sp_site%>%group_by(Site)%>%
                          summarise_all(first)%>%dplyr::select(Site, Lat), by='Site')

aus_ovl_comp<-aus_ovl_comp[-which(is.na(aus_ovl_comp$Site)),] # remove na row from FG 16 no overlap

sg_lat_spans<-data.frame(xmin=c(-23.4, -24.8, -26.61, -28.19, -29.9, -29.97), 
                         xmax=c(-24.1, -25.3, -26.98, -28.616, -30.96, -30.3))
p1<-ggplot()+
  geom_rect(data=sg_lat_spans, aes(ymin=0, ymax=2, xmin=xmin, xmax=xmax), fill='grey', alpha=0.5)+
  geom_smooth(data=aus_trop_prop, aes(x=Lat, y=trop_met), se=F, colour='black')+
  geom_smooth(data=aus_ovl_comp, aes(x=Lat, y=Freq), colour='red', se=F)+
  geom_smooth(data=aus_ovl_filt, aes(x=Lat, y=Freq), colour='blue', se=F)+
  geom_point(data=filter(aus_trop_prop, trop_met<2), aes(x=Lat, y=trop_met))+
  geom_point(data=aus_ovl_comp, aes(x=Lat, y=Freq), colour='red')+
  geom_point(data=aus_ovl_filt, aes(x=Lat, y=Freq), colour='blue')+
  theme_bw()+facet_grid(FG~.)+scale_x_reverse()
# note removal of outliers on geom_point: allows smooth tofit to full data
# but not show these outliers on plot

## calc correlation coef and significance

ovl_test_aus<-aus_trop_prop%>%group_by(FG, Site)%>%summarise(trop_met=mean(trop_met))%>%
  left_join(., ovl_site_aus[[1]][c(1,2,5)], by=c('FG', 'Site'))%>%
  left_join(., ovl_site_aus[[2]][c(1,2,5)], by=c('FG', 'Site'))

m1<-lm(trop_met~Freq.x, data=filter(ovl_test_aus, FG==1))

levz<-c('trop.base', 'trans.bay', 'trans.offshore', 'trans.temp',
                           'temp.inshore', 'temp.offshore')

aus_sp_site_pco$site.group<-factor(aus_sp_site_pco$site.group,levels=levz)
names(ovl_site.group_aus[[3]])[names(ovl_site.group_aus[[3]])=='Site']<-'site.group'
names(ovl_site.group_aus[[2]])[names(ovl_site.group_aus[[2]])=='Site']<-'site.group'
names(ovl_site.group_aus[[1]])[names(ovl_site.group_aus[[1]])=='Site']<-'site.group'

ovl_site.group_aus[[3]]$site.group<-factor(ovl_site.group_aus[[3]]$site.group,levels=levz)
ovl_site.group_aus[[1]]$site.group<-factor(ovl_site.group_aus[[1]]$site.group,levels=levz)
ovl_site.group_aus[[2]]$site.group<-factor(ovl_site.group_aus[[2]]$site.group,levels=levz)

ovl_site.group_aus[[1]]<-ovl_site.group_aus[[1]][-which(is.na(ovl_site.group_aus[[1]]$site.group)),]

ap1<-aus_sp_site_pco[aus_sp_site_pco$FG %in% c(15, 10, 8, 2,6, 12, 4,1, 16),]
ap1$FG<-factor(ap1$FG)

ap1$ThermalAffinity2<-factor(ap1$ThermalAffinity2, levels=c('tropical', 'subtropical'))
ovl_site.group_aus[[3]]$ThermalAffinity2<-factor(ovl_site.group_aus[[3]]$ThermalAffinity2, levels=c('tropical', 'subtropical'))

p2<-ppp+geom_point(data=ap1, aes(x=A1, y=A2, colour=ThermalAffinity2))+
  geom_sf(data=ovl_site.group_aus[[3]], aes(fill=ThermalAffinity2), alpha=0.5)+
  geom_text(data=ovl_site.group_aus[[1]], aes(x=0.9, y=0.9, label=paste(round(Freq, 2)*100,'%', sep='')))+
  geom_text(data=ovl_site.group_aus[[2]], aes(x=0.7, y=0.7, label=paste(round(Freq, 2)*100,'%', sep='')), colour='red')+
  scale_y_continuous(paste('PC2', sprintf('(%0.1f%% explained var.)',
                                          100* func_dudi$eig[2]/sum(func_dudi$eig[func_dudi$eig>0.007]))))+
  scale_x_continuous(paste('PC1', sprintf('(%0.1f%% explained var.)',
                                          100* func_dudi$eig[1]/sum(func_dudi$eig[func_dudi$eig>0.007]))))+
  facet_grid(FG~site.group)+theme_minimal()+theme(legend.position = "none")

p2a<-ppp+
  geom_sf(data=ovl_site.group_aus[[3]], aes(colour=ThermalAffinity2), fill=NA)+
  geom_text(data=ovl_site.group_aus[[1]], aes(x=0.9, y=0.9, label=paste(round(Freq, 2)*100,'%', sep='')))+
  geom_text(data=ovl_site.group_aus[[2]], aes(x=0.7, y=0.7, label=paste(round(Freq, 2)*100,'%', sep='')), colour='red')+
  scale_y_continuous(paste('PC2', sprintf('(%0.1f%% explained var.)',
                                          100* func_dudi$eig[2]/sum(func_dudi$eig[func_dudi$eig>0.007]))))+
  scale_x_continuous(paste('PC1', sprintf('(%0.1f%% explained var.)',
                                          100* func_dudi$eig[1]/sum(func_dudi$eig[func_dudi$eig>0.007]))))+
  facet_grid(FG~site.group)+theme_minimal()+theme(legend.position = "none")


#### Japan overlap by site ####

# trial plot for site level results
# fill overlap with 0's where no overlap
ovl_site_jpn[[1]]<-ovl_site_jpn[[1]]%>%tidyr::complete(Site, FG, fill=list(Freq=0))
ovl_site_jpn[[2]]<-ovl_site_jpn[[2]]%>%tidyr::complete(Site, FG, fill=list(Freq=0))

jpn_ovl_comp<-left_join(ovl_site_jpn[[1]],jpn_sp_site%>%group_by(SiteID)%>%summarise_all(first)%>%dplyr::select(SiteID, lat), by=c('Site'='SiteID'))
jpn_ovl_filt<-left_join(ovl_site_jpn[[2]],jpn_sp_site%>%group_by(SiteID)%>%summarise_all(first)%>%dplyr::select(SiteID, lat), by=c('Site'='SiteID'))

jpn_trop_prop$FG<-factor(jpn_trop_prop$FG,
                         levels=c(15, 10, 8, 2,6, 12, 4,1, 16))
jpn_ovl_comp$FG<-factor(jpn_ovl_comp$FG,
                        levels=c(15, 10, 8, 2,6, 12, 4,1, 16))
jpn_ovl_filt$FG<-factor(jpn_ovl_filt$FG,
                        levels=c(15, 10, 8, 2,6, 12, 4,1, 16))

jpn_ovl_comp<-jpn_ovl_comp[-which(is.na(jpn_ovl_comp$Site)),]

sg_lat_spans<-data.frame(xmin=c(24.2, 26.2, 28.5, 31,   32.7,  33.38, 34.6), 
                         xmax=c(24.5, 28.4, 30.5, 31.6, 32.82, 33.5, 35 ))

p3<-ggplot()+
  geom_rect(data=sg_lat_spans, aes(ymin=0, ymax=2, xmin=xmin, xmax=xmax), fill='grey', alpha=0.5)+
  geom_smooth(data=jpn_trop_prop, aes(x=lat, y=trop_met), se=F, colour='black')+
  geom_smooth(data=jpn_ovl_comp, aes(x=lat, y=Freq), colour='red', se=F)+
  geom_smooth(data=jpn_ovl_filt, aes(x=lat, y=Freq), colour='blue', se=F)+
  geom_point(data=jpn_trop_prop[jpn_trop_prop$trop_met<2,], aes(x=lat, y=trop_met))+
  geom_point(data=jpn_ovl_comp, aes(x=lat, y=Freq), colour='red')+
  geom_point(data=jpn_ovl_filt, aes(x=lat, y=Freq), colour='blue')+
  theme_bw()+facet_grid(FG~.)
# FYI tropicalization curve is same when using site-aggregated mean data
# jpn_trop_prop%>%group_by(FG, SiteID)%>%summarise_all(mean)

#### calc correlation coef and significance ####

ovl_test_jpn<-jpn_trop_prop%>%group_by(FG, SiteID)%>%summarise(trop_met=mean(trop_met))%>%
  left_join(., ovl_site_jpn[[1]][c(1,2,5)], by=c('FG', 'SiteID'='Site'))%>%
  left_join(., ovl_site_jpn[[2]][c(1,2,5)], by=c('FG', 'SiteID'='Site'))

qplot(data=ovl_test_jpn, x=Freq.y, y=trop_met)+facet_wrap(~FG, scales='free')

# tests for aus and japan then write out
spear_tests<-rbind(
ovl_test_jpn%>%group_by(FG)%>%summarise(comp_p=cor.test(trop_met, Freq.x, method = 'kendall')$p.value,
                                        comp_est=cor.test(trop_met, Freq.x, method = 'kendall')$estimate,
                                        filt_p=cor.test(trop_met, Freq.y, method = 'kendall')$p.value,
                                        filt_est=cor.test(trop_met, Freq.y, method = 'kendall')$estimate),
ovl_test_aus%>%group_by(FG)%>%summarise(comp_p=cor.test(trop_met, Freq.x, method = 'kendall')$p.value,
                                        comp_est=cor.test(trop_met, Freq.x, method = 'kendall')$estimate,
                                        filt_p=cor.test(trop_met, Freq.y, method = 'kendall')$p.value,
                                        filt_est=cor.test(trop_met, Freq.y, method = 'kendall')$estimate))
# actually using kendall's tau!
#write.csv(spear_tests, 'C:/coral_fish/outputs/func_overlap_correlation_tropicalization.csv', quote=F, row.names=F)

# Do correlation test between competition and subtropical biomass

ovl_test_jpn<-left_join(ovl_test_jpn, bio_jpn%>%filter(.,ThermalAffinity2=='subtropical')%>%
  group_by(FG, SiteID)%>%summarise(cor_biom=mean(cor_biom)), by=c('FG', 'SiteID'))
  
ovl_test_aus<-left_join(ovl_test_aus, bio_aus%>%filter(.,ThermalAffinity2=='subtropical')%>%
  group_by(FG, Site)%>%summarise(cor_biom=mean(cor_biom)), by=c('FG', 'Site'))

qplot(data=ovl_test_jpn, x=Freq.x, y=cor_biom^0.25)+facet_wrap(~FG, scales='free')
qplot(data=ovl_test_aus, x=Freq.x, y=cor_biom^0.25)+facet_wrap(~FG, scales='free')

ovl_test_aus%>%group_by(FG)%>%summarise(comp_p=cor.test(Freq.x, cor_biom^0.25, method = 'kendall')$p.value,
                                   comp_est=cor.test(Freq.x, cor_biom^0.25, method = 'kendall')$estimate)

ovl_test_jpn%>%group_by(FG)%>%summarise(comp_p=cor.test(Freq.x, cor_biom^0.25, method = 'kendall')$p.value,
                                   comp_est=cor.test(Freq.x, cor_biom^0.25, method = 'kendall')$estimate)
# not used in the end

# Do correlation test with subtropical biomass

aus_cor<-left_join(aus_trop_prop, bio_aus[bio_aus$ThermalAffinity2=='subtropical',c(1, 2, 11)],
          by=c('Site.trans.ID', 'FG'))
jpn_cor<-left_join(jpn_trop_prop, bio_jpn[bio_jpn$ThermalAffinity2=='subtropical',c(1, 2, 9)],
                   by=c('Site.trans.ID', 'FG'))

qplot(data=aus_cor, x=cor_biom.y^0.25, y=cor_biom.x^0.25)+facet_wrap(~FG)
qplot(data=aus_cor, x=cor_biom.y^0.25, y=trop_met)+facet_wrap(~FG)

aus_cor%>%group_by(FG)%>%summarise(comp_p=cor.test(trop_met, cor_biom.y^0.25, method = 'kendall')$p.value,
                                        comp_est=cor.test(trop_met, cor_biom.y^0.25, method = 'kendall')$estimate)

jpn_cor%>%group_by(FG)%>%summarise(comp_p=cor.test(trop_met, cor_biom.y^0.25, method = 'kendall')$p.value,
                                   comp_est=cor.test(trop_met, cor_biom.y^0.25, method = 'kendall')$estimate)
# not used in the end




levz<-c('trop.base', 'trop.island', 'trans.island', 'trans.inland',
        'trans.headld', 'temp.headld')

jpn_sp_site_pco$site.group<-factor(jpn_sp_site_pco$site.group,levels=levz)
names(ovl_site.group_jpn[[3]])[names(ovl_site.group_jpn[[3]])=='Site']<-'site.group'
names(ovl_site.group_jpn[[2]])[names(ovl_site.group_jpn[[2]])=='Site']<-'site.group'
names(ovl_site.group_jpn[[1]])[names(ovl_site.group_jpn[[1]])=='Site']<-'site.group'

ovl_site.group_jpn[[3]]$site.group<-factor(ovl_site.group_jpn[[3]]$site.group,levels=levz)
ovl_site.group_jpn[[1]]$site.group<-factor(ovl_site.group_jpn[[1]]$site.group,levels=levz)
ovl_site.group_jpn[[2]]$site.group<-factor(ovl_site.group_jpn[[2]]$site.group,levels=levz)

ovl_site.group_jpn[[1]]<-ovl_site.group_jpn[[1]][-which(is.na(ovl_site.group_jpn[[1]]$site.group)),]

ap1<-jpn_sp_site_pco[jpn_sp_site_pco$FG %in% c(15, 10, 8, 2,6, 12, 4,1, 16),]
ap1$FG<-factor(ap1$FG)

ap1$ThermalAffinity2<-factor(ap1$ThermalAffinity2, levels=c('tropical', 'subtropical'))
ovl_site.group_jpn[[3]]$ThermalAffinity2<-factor(ovl_site.group_jpn[[3]]$ThermalAffinity2, levels=c('tropical', 'subtropical'))

p4<-ppp+geom_point(data=ap1, aes(x=A1, y=A2, colour=ThermalAffinity2))+
  geom_sf(data=ovl_site.group_jpn[[3]], aes(fill=ThermalAffinity2), alpha=0.5)+
  geom_text(data=ovl_site.group_jpn[[1]], aes(x=0.9, y=0.9, label=paste(round(Freq, 2)*100,'%', sep='')))+
  geom_text(data=ovl_site.group_jpn[[2]], aes(x=0.7, y=0.7, label=paste(round(Freq, 2)*100,'%', sep='')), colour='red')+
  scale_y_continuous(paste('PC2', sprintf('(%0.1f%% explained var.)',
                                          100* func_dudi$eig[2]/sum(func_dudi$eig[func_dudi$eig>0.007]))))+
  scale_x_continuous(paste('PC1', sprintf('(%0.1f%% explained var.)',
                                          100* func_dudi$eig[1]/sum(func_dudi$eig[func_dudi$eig>0.007]))))+
  facet_grid(FG~site.group)+theme_minimal()+theme(legend.position = "none")

p4a<-ppp+
  geom_sf(data=ovl_site.group_jpn[[3]], aes(colour=ThermalAffinity2), fill=NA)+
  geom_text(data=ovl_site.group_jpn[[1]], aes(x=0.9, y=0.9, label=paste(round(Freq, 2)*100,'%', sep='')))+
  geom_text(data=ovl_site.group_jpn[[2]], aes(x=0.7, y=0.7, label=paste(round(Freq, 2)*100,'%', sep='')), colour='red')+
  scale_y_continuous(paste('PC2', sprintf('(%0.1f%% explained var.)',
                                          100* func_dudi$eig[2]/sum(func_dudi$eig[func_dudi$eig>0.007]))))+
  scale_x_continuous(paste('PC1', sprintf('(%0.1f%% explained var.)',
                                          100* func_dudi$eig[1]/sum(func_dudi$eig[func_dudi$eig>0.007]))))+
  facet_grid(FG~site.group)+theme_minimal()+theme(legend.position = "none")


library(rvg)
library(officer)

read_pptx('C:/coral_fish/outputs/portrait_template.pptx') %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  ph_with(dml(ggobj=p1), location = ph_location_fullsize()) %>% 
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  ph_with(dml(ggobj=p2), location = ph_location_fullsize()) %>% 
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  ph_with(dml(ggobj=p3), location = ph_location_fullsize()) %>% 
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  ph_with(dml(ggobj=p4), location = ph_location_fullsize()) %>% 
  print(target = 'C:/coral_fish/outputs/aus_lat_ovl.pptx')



ggsave('C:/coral_fish/outputs/fig5_Australia_trends.eps',
       plot=p1,width = 21, height = 30, units = "cm")
ggsave('C:/coral_fish/outputs/fig5_Australia_niche.eps',
       plot=p2,width = 21, height = 30, units = "cm")
ggsave('C:/coral_fish/outputs/fig5_Japan_trends.eps',
       plot=p3,width = 21, height = 30, units = "cm")
ggsave('C:/coral_fish/outputs/fig5_Japan_niche.eps',
       plot=p4,width = 21, height = 30, units = "cm")

ggsave('C:/coral_fish/outputs/fig5_Australia_nicheALT.eps',
       plot=p2a,width = 21, height = 30, units = "cm")
ggsave('C:/coral_fish/outputs/fig5_Japan_nicheALT.eps',
       plot=p4a,width = 21, height = 30, units = "cm")


##### single FG function

foverlapPlotz<-function(trop_metD=aus_trop_prop, site_comp=aus_ovl_comp, site_filt=aus_ovl_filt,
                        pco_pt=aus_sp_site_pco, pco_comp_text=ovl_site.group_aus[[1]],
                        pco_filt_text=ovl_site.group_aus[[2]], pco_kern=ovl_site.group_aus[[3]],
                        myFGz=15, aus.jpn='aus')
{
  #for(FGz in myFGz)
  #{  
  if(aus.jpn=='aus'){sg_lat_spans<-data.frame(xmin=c(-23.4, -24.8, -26.61, -28.19, -29.9, -29.97), 
                                              xmax=c(-24.1, -25.3, -26.98, -28.616, -30.96, -30.3))}
  if(aus.jpn=='jpn'){sg_lat_spans<-data.frame(xmin=c(24.2, 26.2, 28.5, 31, 32.7,  33.38, 34.6), 
                                              xmax=c(24.5, 28.4, 30.5, 31.6, 32.82, 33.5, 35 ))}
  
  
  p1<-ggplot()+
    geom_rect(data=sg_lat_spans, aes(ymin=0, ymax=2, xmin=xmin, xmax=xmax), fill='grey', alpha=0.5)+
    geom_smooth(data=filter(trop_metD, FG==FGz), aes(x=Lat, y=trop_met), se=F, colour='black')+
    geom_smooth(data=filter(site_comp, FG==FGz), aes(x=Lat, y=Freq), colour='red', se=F)+
    geom_smooth(data=filter(site_filt, FG==FGz), aes(x=Lat, y=Freq), colour='blue', se=F)+
    geom_point(data=filter(trop_metD, FG==FGz & trop_met<2), aes(x=Lat, y=trop_met))+
    geom_point(data=filter(site_comp, FG==FGz), aes(x=Lat, y=Freq), colour='red')+
    geom_point(data=filter(site_filt, FG==FGz), aes(x=Lat, y=Freq), colour='blue')+
    theme_bw()+theme(plot.margin=unit(c(1,1,-0.7,1), "cm"))
  
  if(aus.jpn=='aus'){p1<-p1+scale_x_reverse()}
  
  names(pco_kern)[names(pco_kern)=='Site']<-'site.group'
  names(pco_comp_text)[names(pco_comp_text)=='Site']<-'site.group'
  names(pco_filt_text)[names(pco_filt_text)=='Site']<-'site.group'
  
  pco_kern$ThermalAffinity2<-factor(pco_kern$ThermalAffinity2, levels=c('tropical', 'subtropical'))
  
  if(aus.jpn=='aus'){levz<-c('trop.base', 'trans.bay', 'trans.offshore', 'trans.temp',
                             'temp.inshore', 'temp.offshore')}
  if(aus.jpn=='jpn'){levz<-c('trop.base', 'trop.island', 'trans.island', 'trans.inland',
                             'trans.headld', 'temp.headld')}
  
  pco_pt$site.group<-factor(pco_pt$site.group,levels=levz)
  pco_kern$site.group<-factor(pco_kern$site.group,levels=levz)
  pco_comp_text$site.group<-factor(pco_comp_text$site.group,levels=levz)
  pco_filt_text$site.group<-factor(pco_filt_text$site.group,levels=levz)
  
  p2<-ggplot()+geom_point(data=filter(pco_pt, FG==FGz), aes(x=A1, y=A2, colour=ThermalAffinity2))+
    geom_sf(data=filter(pco_kern, FG==FGz), aes(fill=ThermalAffinity2), alpha=0.5)+
    geom_text(data=filter(pco_comp_text, FG==FGz), aes(x=max(filter(pco_pt, FG==FGz)$A1)+0.02, y=min(filter(pco_pt, FG==FGz)$A2)+0.2, label=paste(round(Freq, 2)*100,'%', sep='')), colour='red')+
    geom_text(data=filter(pco_filt_text, FG==FGz), aes(x=max(filter(pco_pt, FG==FGz)$A1)+0.02, y=min(filter(pco_pt, FG==FGz)$A2), label=paste(round(Freq, 2)*100,'%', sep='')), colour='blue')+
    facet_grid(~site.group)+
    scale_x_continuous(breaks=round(seq(min(filter(pco_pt, FG==FGz)$A1),
                                        max(filter(pco_pt, FG==FGz)$A1), by = 0.4), 1))+
    scale_y_continuous(breaks=round(seq(min(filter(pco_pt, FG==FGz)$A2),
                                        max(filter(pco_pt, FG==FGz)$A2), by = 0.4), 1))+
    theme_minimal()+xlab('PCoA1')+ylab('PCoA2')+
    theme(legend.position = "none",strip.text.x = element_blank(),plot.margin=unit(c(-3,1,1,1), "cm"))
  
  #return(list(p1, p2)) 
  
  
}

#### combine site-group tropicalization lme models with overlap results ####

# Add site.group classes to biomass data

bio_jpn$site.group<-'trop.base'
bio_jpn[bio_jpn$lat > 25 &  bio_jpn$lat < 28.5,]$site.group<-'trop.island'
bio_jpn[bio_jpn$lat > 28.5 &  bio_jpn$lat < 31,]$site.group<-'trans.island'
bio_jpn[bio_jpn$SiteID %in% c('JP28', 'JP29', 'JP30', 'JP31'),]$site.group<-'trans.inland'
bio_jpn[bio_jpn$SiteID %in% c('JP8', 'JP9', 'JP10', 'JP11', 'JP12'),]$site.group<-'trans.headld'
bio_jpn[bio_jpn$lat > 34,]$site.group<-'temp.headld'
table(bio_jpn$SiteID, bio_jpn$site.group)

bio_aus$site.group<-'trop.base'
bio_aus[bio_aus$Lat > -25.6 &  bio_aus$Lat < -24.5,]$site.group<-'trans.bay'
bio_aus[bio_aus$Lat > -28 &  bio_aus$Lat < -25.6,]$site.group<-'trans.offshore'
bio_aus[bio_aus$Lat < -28,]$site.group<-'temp.offshore'
bio_aus[bio_aus$Site %in% c('Julian Rock False Trench', 'Julian Rock Nursery', 'Cook Island'),]$site.group<-'trans.temp'
bio_aus[bio_aus$Site %in% c('Muttonbird Island', 'Woolgoolga Reef', 'Woolgoolga Headland', 'North Rock'),]$site.group<-'temp.inshore'

# summarise biomass data to site.group and FG, then filter to jsut subtropical
bio_jpn_sg<-bio_jpn%>%group_by(site.group, FG, ThermalAffinity2)%>%
  summarise(cor_biom=sum(cor_biom), lat=max(lat))%>%filter(., ThermalAffinity2=='subtropical')
bio_aus_sg<-bio_aus%>%group_by(site.group, FG, ThermalAffinity2)%>%
  summarise(cor_biom=sum(cor_biom), Lat=min(Lat))%>%filter(., ThermalAffinity2=='subtropical')

##### JAPAN ####

trop_comps_out<-trop_comps_out[trop_comps_out$FG!='comm',1:7]
trop_comps_out<-rbind(trop_comps_out, data.frame(site.group='trop.base', FG=c(15, 10, 8, 2,6, 12,4, 1, 16),
                   emmean=1, SE=0, df=0, lower.CL=1, upper.CL=1))

trop_expl<-left_join(trop_comps_out, jpn_ovl_comp[,c(1,2,5)], by=c('site.group', 'FG'))
names(trop_expl)[8]<-'comp_ovl'
trop_expl[is.na(trop_expl$comp_ovl),]$comp_ovl<-0
trop_expl<-left_join(trop_expl, jpn_ovl_filt[,c(5,7,3)], by=c('site.group', 'FG'))
names(trop_expl)[9]<-'filt_ovl'
trop_expl[is.na(trop_expl$filt_ovl),]$filt_ovl<-0
trop_expl[trop_expl$site.group=='trop.base',]$filt_ovl<-1
trop_expl<-left_join(trop_expl, bio_jpn_sg[,c(1,2,4,5)], by=c('site.group', 'FG'))
trop_expl$cor_biom<-trop_expl$cor_biom^0.25 #apply 4rt transformation
trop_expl%>%group_by(FG)%>%mutate(max_biom=max(cor_biom),
                                  res=resid(lm(emmean~poly(lat, 2))))->trop_expl

trop_expl$site.group<-factor(trop_expl$site.group,
                    levels=c('trop.base', 'trop.island', 'trans.island', 'trans.inland',
                     'trans.headld', 'temp.headld'))
trop_expl$FG<-factor(trop_expl$FG, levels=c(15, 10, 8, 2,6,12, 4,1, 16))

ggplot(data=trop_expl, aes(x=site.group))+
  geom_point(aes(y=comp_ovl, group=FG), stat='summary', fun.y=sum, colour='red') +
  stat_summary(aes(y=comp_ovl, group=FG),fun.y=sum, geom="line", colour='red', linetype='dotted')+
  geom_point(aes(y=filt_ovl, group=FG), stat='summary', fun.y=sum, colour='green', shape=17) +
  stat_summary(aes(y=filt_ovl, group=FG),fun.y=sum, geom="line", colour='green', linetype='dotted')+
  geom_pointrange(aes(y=emmean, ymin=lower.CL, ymax=upper.CL))+
  facet_wrap(~FG)+theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0))


ggplot(data=trop_expl, aes(x=emmean))+
  geom_point(aes(y=filt_ovl))+geom_smooth(aes(y=filt_ovl,colour=FG), method='lm', se=F)+
   geom_point(aes(y=comp_ovl), shape=1)+geom_smooth(aes(y=comp_ovl,colour=FG), method='lm', formula=y~poly(x,2),se=F, linetype='dashed')+theme_bw()+facet_wrap(~FG)

# models

m1<-lm(emmean~filt_ovl+ FG + FG:filt_ovl, data=trop_expl)
plot(m1)
summary(m1)
anova(m1)
m1<-lm(emmean~filt_ovl+ FG, data=trop_expl)
summary(m1)
anova(m1)

nd<-expand.grid(filt_ovl=seq(0, 1, 0.05), FG=factor(c(15, 10, 8, 2,6,12, 4,1, 16)))
nd$p1<-predict(m1, nd)

ggplot(data=trop_expl, aes(x=filt_ovl))+
  geom_point(aes(y=emmean,colour=FG))+
  geom_line(data=nd, aes(y=p1,colour=FG))+
  theme_bw()+facet_wrap(~FG)

m2<-lm(emmean~comp_ovl, data=trop_expl[trop_expl$FG!=16,])
plot(m2)
summary(m2)
anova(m2)
nd<-expand.grid(comp_ovl=seq(0, 1, 0.05), FG=factor(c(15, 10, 8, 2,6,12, 4,1)))
nd$p1<-predict(m2, nd)

ggplot(data=trop_expl, aes(x=comp_ovl))+
  geom_point(aes(y=emmean,colour=FG))+
  geom_line(data=nd, aes(y=p1,colour=FG))+
  theme_bw()+facet_wrap(~FG)


#### Australia ####

trop_comps_out_aus<-trop_comps_out_aus[trop_comps_out_aus$FG!='comm',1:7]
trop_comps_out_aus<-rbind(trop_comps_out_aus, data.frame(site.group='trop.base', FG=c(15, 10, 8, 2,6, 12,4, 1, 16),
                                                 emmean=1, SE=0, df=0, lower.CL=1, upper.CL=1))

trop_expl_aus<-left_join(trop_comps_out_aus, aus_ovl_comp[,c(1,2,5)], by=c('site.group', 'FG'))
names(trop_expl_aus)[8]<-'comp_ovl'
trop_expl_aus[is.na(trop_expl_aus$comp_ovl),]$comp_ovl<-0
trop_expl_aus<-left_join(trop_expl_aus, aus_ovl_filt[,c(5,7,3)], by=c('site.group', 'FG'))
names(trop_expl_aus)[9]<-'filt_ovl'
trop_expl_aus[is.na(trop_expl_aus$filt_ovl),]$filt_ovl<-0
trop_expl_aus[trop_expl_aus$site.group=='trop.base',]$filt_ovl<-1
trop_expl_aus<-left_join(trop_expl_aus, bio_aus_sg[,c(1,2,4,5)], by=c('site.group', 'FG'))
trop_expl_aus$cor_biom<-trop_expl_aus$cor_biom^0.25 #apply 4rt transformation
trop_expl_aus%>%group_by(FG)%>%mutate(max_biom=max(cor_biom),
                                  res=resid(lm(emmean~poly(Lat, 2))))->trop_expl_aus

trop_expl_aus$site.group<-factor(trop_expl_aus$site.group,
                          levels=c('trop.base', 'trans.bay', 'trans.offshore', 'trans.temp',
                           'temp.inshore', 'temp.offshore'))
trop_expl_aus$FG<-factor(trop_expl_aus$FG, levels=c(15, 10, 8, 2,6,12, 4,1, 16))


ggplot(data=trop_expl_aus, aes(x=site.group))+
  geom_point(aes(y=comp_ovl, group=FG), stat='summary', fun.y=sum, colour='red') +
  stat_summary(aes(y=comp_ovl, group=FG),fun.y=sum, geom="line", colour='red', linetype='dotted')+
  geom_point(aes(y=filt_ovl, group=FG), stat='summary', fun.y=sum, colour='green', shape=17) +
  stat_summary(aes(y=filt_ovl, group=FG),fun.y=sum, geom="line", colour='green', linetype='dotted')+
  geom_pointrange(aes(y=emmean, ymin=lower.CL, ymax=upper.CL))+
  facet_wrap(~FG)+theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust=0))


ggplot(data=trop_expl_aus, aes(x=emmean))+
  geom_point(aes(y=filt_ovl))+geom_smooth(aes(y=filt_ovl,colour=FG), method='lm', se=F)+
  geom_point(aes(y=comp_ovl), shape=1)+geom_smooth(aes(y=comp_ovl,colour=FG), method='lm', formula=y~poly(x,2),se=F, linetype='dashed')+theme_bw()+facet_wrap(~FG)






dat_jpn$ThermalAffinity2<-factor(dat_jpn$ThermalAffinity2, levels=c('tropical', 'subtropical'))



ppp+geom_point(data=dat_jpn%>% filter(JPN_temp==1 & FG %in% c(4,6,1,2)),
               aes(x=A1, y=A2, colour=ThermalAffinity2))+
  geom_polygon(data=dat_jpn %>% filter(JPN_temp==1 & FG %in% c(4,6,1,2))%>%
                 group_by(FG, ThermalAffinity2) %>% slice(chull(A1, A2)),
               aes(x=A1, y=A2, fill=ThermalAffinity2), alpha=0.2)+
  geom_point(data=dat_jpn%>% filter(JPN_temp==1 & FG %in% c(4,6,1,2))%>%
               group_by(FG, ThermalAffinity2)%>%summarize(A1=mean(A1), A2=mean(A2)),
             aes(x=A1, y=A2, colour=ThermalAffinity2), size=2)+
  geom_point(data=dat_jpn%>% filter(JPN_temp==1 & FG %in% c(4,6,1,2))%>%
               group_by(FG, ThermalAffinity2)%>%summarize(A1=mean(A1), A2=mean(A2)),
             aes(x=A1, y=A2, group=ThermalAffinity2), colour='black', size=2, shape=1)+
  scale_y_continuous(paste('PC2', sprintf('(%0.1f%% explained var.)',
                                          100* func_dudi$eig[2]/sum(func_dudi$eig))))+
  scale_x_continuous(paste('PC1', sprintf('(%0.1f%% explained var.)',
                                          100* func_dudi$eig[1]/sum(func_dudi$eig))))+
  facet_wrap(~FG)+theme_minimal()


# euc_distance vals
fd_aus<-dat_aus%>% filter(AUS_temp==1 & FG %in% c(4,6,1,2))%>%
  group_by(FG,  ThermalAffinity2)%>%summarize(A1=mean(A1), A2=mean(A2))%>%ungroup()%>%
  group_by(FG)%>%summarise(fn_d=ifelse(n()>1, dist(cbind(A1, A2)), 0.5))%>%
  mutate(comm='temp', region='Australia')

fd_jpn<-dat_jpn%>% filter(JPN_temp==1 & FG %in% c(4,6,1,2))%>%
  group_by(FG,  ThermalAffinity2)%>%summarize(A1=mean(A1), A2=mean(A2))%>%ungroup()%>%
  group_by(FG)%>%summarise(fn_d=ifelse(n()>1, dist(cbind(A1, A2)), 0.5))%>%
  mutate(comm='temp', region='Japan')

fdz<-rbind(fd_aus, fd_jpn)

trans_testdat<-left_join(filter(delta_trop, comm=='temp'& FG %in% c(1,2,4,6)),
                         fdz, by=c('region', 'FG'))

# calc anbsolute difference (n species different from expected rather than proportion)
# this adds in differences between groups in terms of size

trans_testdat$comm_trop_absol<-round(trans_testdat$n*trans_testdat$comm_trop)
trans_testdat$n_trop<-trans_testdat$n*trans_testdat$prop_trop
trans_testdat$delta_trop_absol<-trans_testdat$n_trop-trans_testdat$comm_trop_absol

qplot(data=trans_testdat, x=delta_trop, y=delta_trop_absol, colour=factor(FG))+
  geom_point(data=filter(trans_testdat, n==1), aes(x=delta_trop, y=delta_trop_absol),
             shape=1, colour='black')+facet_wrap(~region)

qplot(data=trans_testdat, x=delta_trop, y=fn_d, colour=factor(FG))
qplot(data=trans_testdat, x=delta_trop, y=fn_d, colour=factor(FG), shape=region)
qplot(data=trans_testdat, x=delta_trop_absol, y=fn_d, colour=factor(FG), shape=region)


# projection of FGs into 'colonisation space'

d_colon<-daisy(dat[c('BodySize', 'PLD', 'ParentalMode')], metric='gower', stand = FALSE)

colon_dudi<-dudi.pco(d = sqrt(d_colon), scannf = FALSE, nf = 4)


efit <- envfit(colon_dudi, dat[c('BodySize', 'PLD', 'ParentalMode')], na.rm=T)
varibs<-data.frame(rbind(efit$vectors$arrows, efit$factors$centroids))
varibs$predictors=row.names(varibs)
varibs$predictors<-gsub('ParentalMode', 'PM.',varibs$predictors)

ppp+geom_segment(data=varibs, aes(y=0, x=0, xend=A1, yend=A2),
                 arrow=arrow(length=unit(0.3,'lines')), colour='red')+
  geom_text(data=varibs, aes(x=A1, y=A2, label=predictors))+
  theme_minimal()+
  scale_y_continuous(paste('PC2', sprintf('(%0.1f%% explained var.)',
                                          100* colon_dudi$eig[2]/sum(colon_dudi$eig))))+
  scale_x_continuous(paste('PC1', sprintf('(%0.1f%% explained var.)',
                                          100* colon_dudi$eig[1]/sum(colon_dudi$eig))))

dat_aus$ThermalAffinity2<-factor(dat_aus$ThermalAffinity2, levels=c('tropical', 'subtropical'))

ppp+geom_point(data=dat_aus%>% filter(AUS_temp==1 & FG %in% c(4,6,1,2)),
               aes(x=A1.colon, y=A2.colon, colour=ThermalAffinity2))+
  geom_polygon(data=dat_aus %>% filter(AUS_temp==1 & FG %in% c(4,6,1,2))%>%
                 group_by(FG, ThermalAffinity2) %>% slice(chull(A1.colon, A2.colon)),
               aes(x=A1.colon, y=A2.colon, fill=ThermalAffinity2), alpha=0.2)+
  geom_point(data=dat_aus%>% filter(AUS_temp==1 & FG %in% c(4,6,1,2))%>%
               group_by(FG, ThermalAffinity2)%>%summarize(A1.colon=mean(A1.colon), A2.colon=mean(A2.colon)),
             aes(x=A1.colon, y=A2.colon, colour=ThermalAffinity2), size=2)+
  geom_point(data=dat_aus%>% filter(AUS_temp==1 & FG %in% c(4,6,1,2))%>%
               group_by(FG, ThermalAffinity2)%>%summarize(A1.colon=mean(A1.colon), A2.colon=mean(A2.colon)),
             aes(x=A1.colon, y=A2.colon, group=ThermalAffinity2), colour='black', size=2, shape=1)+
  scale_y_continuous(paste('PC2', sprintf('(%0.1f%% explained var.)',
                                          100* colon_dudi$eig[2]/sum(colon_dudi$eig))))+
  scale_x_continuous(paste('PC1', sprintf('(%0.1f%% explained var.)',
                                          100* colon_dudi$eig[1]/sum(colon_dudi$eig))))+
  facet_wrap(~FG)+theme_minimal()

ppp+geom_point(data=dat_jpn%>% filter(JPN_temp==1 & FG %in% c(4,6,1,2)),
               aes(x=A1.colon, y=A2.colon, colour=ThermalAffinity2))+
  geom_polygon(data=dat_jpn %>% filter(JPN_temp==1 & FG %in% c(4,6,1,2))%>%
                 group_by(FG, ThermalAffinity2) %>% slice(chull(A1.colon, A2.colon)),
               aes(x=A1.colon, y=A2.colon, fill=ThermalAffinity2), alpha=0.2)+
  geom_point(data=dat_jpn%>% filter(JPN_temp==1 & FG %in% c(4,6,1,2))%>%
               group_by(FG, ThermalAffinity2)%>%summarize(A1.colon=mean(A1.colon), A2.colon=mean(A2.colon)),
             aes(x=A1.colon, y=A2.colon, colour=ThermalAffinity2), size=2)+
  geom_point(data=dat_jpn%>% filter(JPN_temp==1 & FG %in% c(4,6,1,2))%>%
               group_by(FG, ThermalAffinity2)%>%summarize(A1.colon=mean(A1.colon), A2.colon=mean(A2.colon)),
             aes(x=A1.colon, y=A2.colon, group=ThermalAffinity2), colour='black', size=2, shape=1)+
  scale_y_continuous(paste('PC2', sprintf('(%0.1f%% explained var.)',
                                          100* colon_dudi$eig[2]/sum(colon_dudi$eig))))+
  scale_x_continuous(paste('PC1', sprintf('(%0.1f%% explained var.)',
                                          100* colon_dudi$eig[1]/sum(colon_dudi$eig))))+
  facet_wrap(~FG)+theme_minimal()



colon_pco<-data.frame(colon_dudi$li,dat)
names(colon_pco)[1]<-'A1.colon'
names(colon_pco)[2]<-'A2.colon'
names(colon_pco)[3]<-'A3.colon'
names(colon_pco)[4]<-'A4.colon'

# add to regional datasets
dat_aus<-left_join(dat_aus, colon_pco[,c(1:5)], by='Species')
dat_jpn<-left_join(dat_jpn, colon_pco[,c(1:5)], by='Species')

colon_lh<-rbind(
dat_aus%>%filter(AUS_temp==1)%>% group_by(FG, ThermalAffinity2)%>%
  summarize(Investment=mean(A1.colon), Dispersal=mean(A2.colon))%>%mutate(region='Australia'),
dat_jpn%>%filter(JPN_temp==1)%>% group_by(FG,ThermalAffinity2)%>%
  summarize(Investment=mean(A1.colon), Dispersal=mean(A2.colon))%>%mutate(region='Japan'))

trans_testdat<-left_join(trans_testdat, filter(colon_lh,
               ThermalAffinity2=='tropical'), by=c('region', 'FG'))
trans_testdat$comm.y<-NULL
trans_testdat$ThermalAffinity2<-NULL
names(trans_testdat)[13]<-'Investment.trop'
names(trans_testdat)[14]<-'Dispersal.trop'
trans_testdat<-left_join(trans_testdat, filter(colon_lh,
 ThermalAffinity2=='subtropical'), by=c('region', 'FG'))
trans_testdat$ThermalAffinity2<-NULL
names(trans_testdat)[15]<-'Investment.resi'
names(trans_testdat)[16]<-'Dispersal.resi'

qplot(data=trans_testdat, x=Investment.trop, y=delta_trop, colour=factor(FG), shape=region)
qplot(data=trans_testdat, x=Dispersal.trop, y=delta_trop, colour=factor(FG), shape=region)

qplot(data=trans_testdat, x=Investment.resi, y=delta_trop, colour=factor(FG), shape=region)
qplot(data=trans_testdat, x=Dispersal.resi, y=delta_trop, colour=factor(FG), shape=region)

m1<-lm(delta_trop~fn_d+Investment.trop+Dispersal.trop+Investment.resi+Dispersal.resi+
         region+factor(FG),data=trans_testdat)

# Modelling resident vs invader dfferences within group

# reclass Parental mode to binary variable
dat_aus$invest_bin<-1
dat_aus[which(dat_aus$ParentalMode=='scatterers'),]$invest_bin<-0
dat_jpn$invest_bin<-1
dat_jpn[which(dat_jpn$ParentalMode=='scatterers'),]$invest_bin<-0

# Functional space

# Australia
aus_mod<-dat_aus[dat_aus$FG %in% c(4,6,2,1) & 
                   dat_aus$AUS_temp==1,]
aus_mod$FG<-factor(aus_mod$FG, levels=c(4,6,1,2))

aus_fgD<-lm(A1~ThermalAffinity2:FG, data=aus_mod)
anova(aus_fgD)

library(emmeans)
#https://cran.r-project.org/web/packages/emmeans/vignettes/interactions.html
emmip(aus_fgD, FG~ThermalAffinity2)
emmeans(aus_fgD, pairwise ~ ThermalAffinity2 | FG)

man_aus<-manova(cbind(A1, A2)~ThermalAffinity2:FG, data=aus_mod)
anova(man_aus)
emmip(man_aus, FG~ThermalAffinity2)
emmeans(man_aus, pairwise ~ ThermalAffinity2 | FG)

# Japan
jpn_mod<-dat_jpn[dat_jpn$FG %in% c(4,6,2,1) & 
                   dat_jpn$JPN_temp==1,]
jpn_mod$FG<-factor(jpn_mod$FG, levels=c(4,6,1,2))

jpn_fgD<-lm(A1~ThermalAffinity2:FG, data=jpn_mod)
anova(jpn_fgD)

emmip(jpn_fgD, FG~ThermalAffinity2)
emmeans(jpn_fgD, pairwise ~ ThermalAffinity2 | FG)

man_jpn<-manova(cbind(A1, A2)~ThermalAffinity2:FG, data=jpn_mod)
anova(man_jpn)
emmip(man_jpn, FG~ThermalAffinity2)
emmeans(man_jpn, pairwise ~ ThermalAffinity2 | FG)


# Life-history space

aus_lhD<-glm(invest_bin~ThermalAffinity2:FG, data=aus_mod, family='binomial')
library(car)
Anova(aus_lhD)
sum(resid(aus_lhD, type='pearson')^2)/df.residual(aus_lhD)

emmip(aus_lhD, FG~ThermalAffinity2)
emmeans(aus_lhD, pairwise ~ ThermalAffinity2 | FG)
table(aus_mod$FG, aus_mod$ThermalAffinity2, aus_mod$invest_bin )

jpn_lhD<-glm(invest_bin~ThermalAffinity2:FG, data=jpn_mod, family='binomial')
library(car)
Anova(jpn_lhD)
sum(resid(jpn_lhD, type='pearson')^2)/df.residual(jpn_lhD)

emmip(jpn_lhD, FG~ThermalAffinity2)
emmeans(jpn_lhD, pairwise ~ ThermalAffinity2 | FG)
table(jpn_mod$FG, jpn_mod$ThermalAffinity2, jpn_mod$invest_bin )

# Modelling inter-group differences in PLD and Life-history

jpn_diff_dat<-dat_jpn%>% filter(JPN_temp==1 & FG %in% c(4,6,1,2) & ThermalAffinity2=='tropical')
aus_diff_dat<-dat_aus%>% filter(AUS_temp==1 & FG %in% c(4,6,1,2) & ThermalAffinity2=='tropical')

jpn_diff_dat$FG<-factor(jpn_diff_dat$FG, levels=c(4,6,1,2))
aus_diff_dat$FG<-factor(aus_diff_dat$FG, levels=c(4,6,1,2))

library(agricolae)
pld_aus<-lm(log(PLD)~FG, data=aus_diff_dat[aus_diff_dat$PLD>0,])
anova(pld_aus)
HSD.test(pld_aus, 'FG', console=T)

aus_lh<-glm(invest_bin~FG, data=aus_diff_dat, family='binomial')
Anova(aus_lh)
sum(resid(aus_lh, type='pearson')^2)/df.residual(aus_lh)

emmeans(aus_lh, pairwise ~ FG)
plot(emmeans(aus_lh, pairwise ~ FG), comparisons=T)


pld_jpn<-lm(log(PLD)~FG, data=jpn_diff_dat[jpn_diff_dat$PLD>0,])
anova(pld_jpn)
HSD.test(pld_jpn, 'FG', console=T)

jpn_lh<-glm(invest_bin~FG, data=jpn_diff_dat, family='binomial')
Anova(jpn_lh)
sum(resid(jpn_lh, type='pearson')^2)/df.residual(jpn_lh)

emmeans(jpn_lh, pairwise ~ FG)
plot(emmeans(jpn_lh, pairwise ~ FG), comparisons=T)

# check ratio of tropical to sub-tropical
subt_rat<-rbind(
dat_jpn%>% filter(JPN_temp==1 & FG %in% c(4,6,1,2))%>%
  group_by(FG) %>% summarise(n_nontrop=length(ThermalAffinity),
 n_subt=length(which(ThermalAffinity=='subtropical')),
 n_temp=length(which(ThermalAffinity=='temperate')))%>%
  mutate(subt_2_temp=n_subt/(n_subt+n_temp), region='Japan'),

dat_aus%>% filter(AUS_temp==1 & FG %in% c(4,6,1,2))%>%
  group_by(FG) %>% summarise(n_nontrop=length(ThermalAffinity),
  n_subt=length(which(ThermalAffinity=='subtropical')),
  n_temp=length(which(ThermalAffinity=='temperate')))%>%
  mutate(subt_2_temp=n_subt/(n_subt+n_temp), region='Australia'))

trans_testdat<-left_join(trans_testdat, subt_rat, by=c('region', 'FG'))

qplot(data=trans_testdat, x=subt_2_temp, y=delta_trop, colour=factor(FG), shape=region)
qplot(data=trans_testdat, x=subt_2_temp, y=delta_trop, colour=factor(FG))+facet_wrap(~region)

m_aus<-lm(delta_trop~fn_d+Investment.trop+Dispersal.trop+Investment.resi+Dispersal.resi+
         subt_2_temp+factor(FG),data=trans_testdat[trans_testdat$region=='Australia',])

summary(m_aus)

m_jpn<-lm(delta_trop~fn_d+Investment.trop+Dispersal.trop+Investment.resi+Dispersal.resi+
            subt_2_temp,data=trans_testdat[trans_testdat$region=='Japan',])

summary(m_jpn)


# plot with variable contribution
#https://www.researchgate.net/post/how_can_i_produce_a_PCoA_biplot_using_R

efit <- envfit(colon_dudi, dat[c('BodySize', 'PLD', 'ParentalMode')], na.rm=T)
varibs<-data.frame(rbind(efit$vectors$arrows, efit$factors$centroids))
varibs$predictors=row.names(varibs)
varibs$predictors<-gsub('ParentalMode', 'PM.',varibs$predictors)


p1<-ppp+geom_point(data=colon_pco%>% filter(JPN_trop==1),
                   aes(x=A1, y=A2), colour='red')+
  geom_point(data=colon_pco%>% filter(JPN_temp==1),
             aes(x=A1, y=A2), colour='blue')+
  geom_polygon(data=colon_pco %>% filter(JPN_trop==1)%>%
                 group_by(FG, ThermalAffinity2) %>% slice(chull(A1, A2)),
               aes(x=A1, y=A2), fill='red', alpha=0.2)+
  geom_polygon(data=colon_pco %>% filter(JPN_temp==1)%>%
                 group_by(FG, ThermalAffinity2) %>% slice(chull(A1, A2)),
               aes(x=A1, y=A2), fill='blue', alpha=0.2)+
  geom_point(data=colon_pco%>% filter(JPN_trop==1)%>%
               group_by(FG, ThermalAffinity2)%>%summarize(A1=mean(A1), A2=mean(A2)), 
             aes(x=A1, y=A2), colour='orange', size=2)+
  geom_point(data=colon_pco%>% filter(JPN_temp==1)%>%
               group_by(FG, ThermalAffinity2)%>%summarize(A1=mean(A1), A2=mean(A2)), 
             aes(x=A1, y=A2), colour='cyan', size=2)+

  scale_y_continuous(paste('PC2', sprintf('(%0.1f%% explained var.)',
                                          100* colon_dudi$eig[2]/sum(colon_dudi$eig))))+
  scale_x_continuous(paste('PC1', sprintf('(%0.1f%% explained var.)',
                                          100* colon_dudi$eig[1]/sum(colon_dudi$eig))))

p1+facet_wrap(~FG+ThermalAffinity2, ncol=2)+theme_minimal()

ppp+geom_segment(data=varibs, aes(y=0, x=0, xend=A1, yend=A2),
                arrow=arrow(length=unit(0.3,'lines')))+
  geom_text(data=varibs, aes(x=A1, y=A2, label=predictors))+
  
  scale_y_continuous(paste('PC2', sprintf('(%0.1f%% explained var.)',
                                          100* colon_dudi$eig[2]/sum(colon_dudi$eig))))+
  scale_x_continuous(paste('PC1', sprintf('(%0.1f%% explained var.)',
                                          100* colon_dudi$eig[1]/sum(colon_dudi$eig))))

# varible contrib plot
p1<-ppp+
  geom_point(data=colon_pco%>% filter(JPN_trop==1)%>%group_by(FG, ThermalAffinity2)%>%
               summarize(A1=mean(A1), A2=mean(A2)), 
             aes(x=A1, y=A2, colour=factor(FG), shape=ThermalAffinity2), size=2)+
  geom_segment(data=varibs, aes(y=0, x=0, yend=A1, xend=A2),
               arrow=arrow(length=unit(0.3,'lines')))+
  geom_text(data=varibs, aes(y=A1, x=A2, label=predictors))+

  scale_y_continuous(paste('PC2', sprintf('(%0.1f%% explained var.)',
                                          100* colon_dudi$eig[2]/sum(colon_dudi$eig))))+
  scale_x_continuous(paste('PC1', sprintf('(%0.1f%% explained var.)',
                                          100* colon_dudi$eig[1]/sum(colon_dudi$eig))))
  p2<-ppp+
  geom_point(data=colon_pco%>% filter(JPN_temp==1)%>%group_by(FG, ThermalAffinity2)%>%
               summarize(A1=mean(A1), A2=mean(A2)), 
             aes(x=A1, y=A2, colour=factor(FG), shape=ThermalAffinity2), size=2)+
  geom_segment(data=varibs, aes(y=0, x=0, yend=A1, xend=A2),
               arrow=arrow(length=unit(0.3,'lines')))+
  geom_text(data=varibs, aes(y=A1, x=A2, label=predictors))+
  
  scale_y_continuous(paste('PC2', sprintf('(%0.1f%% explained var.)',
                                          100* colon_dudi$eig[2]/sum(colon_dudi$eig))))+
  scale_x_continuous(paste('PC1', sprintf('(%0.1f%% explained var.)',
                                          100* colon_dudi$eig[1]/sum(colon_dudi$eig))))
  
  grid.arrange(p1, p2)


# non-multidimensional approach

ggplot(data=dat, aes(x=ParentalMode, y=PLD))+geom_boxplot()+facet_wrap(~FG)

ggplot(data=dat, aes(x=factor(FG), y=PLD))+geom_boxplot()+facet_wrap(~ParentalMode)

dat$ord_ParentalMode<-factor(dat$ParentalMode, 
    levels=c("Live bearers", "brooders","Nesters", "demersal","scatterers"),
    ordered = T)
ggplot(data=dat, aes(x=ord_ParentalMode, y=PLD, colour=factor(FG)))+
  geom_point()+
  geom_point(data=dat%>%group_by(FG)%>%summarize(PM=names(which.max(table(ParentalMode))), 
            PLD=mean(PLD, na.rm=T)),aes(x=PM, y=PLD), size=2, colour='black')+
              facet_wrap(~FG)+theme(axis.text.x = element_text(angle = 90, hjust = 1))

