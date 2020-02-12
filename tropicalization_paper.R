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

# FG to factor
bio_jpn$FG<-factor(bio_jpn$FG)
bio_aus$FG<-factor(bio_aus$FG)

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

#### BIOMASS tropicalization ####

#sanity check to make sure per unit area biomass calc is correct

ggplot(bio_jpn, aes(x = lat, y = tot_biom/totMsurv, colour=ThermalAffinity2)) + 
  geom_point()+geom_smooth(se=F)+facet_wrap(~FG, scales='free')+theme_bw()

ggplot(bio_aus, aes(x = Lat, y = tot_biom/totMsurv, colour=ThermalAffinity2)) + 
  geom_point()+geom_smooth(se=F)+facet_wrap(~FG, scales='free')+
  scale_x_reverse()+theme_bw()

ggplot(spr_jpn, aes(x = lat, y = qD, colour=ThermalAffinity2)) + 
  geom_point()+geom_smooth(se=F)+facet_wrap(~FG, scales='free')+theme_bw()

ggplot(spr_aus, aes(x = Lat, y = qD, colour=ThermalAffinity2)) + 
  geom_point()+geom_smooth(se=F)+facet_wrap(~FG, scales='free')+
  scale_x_reverse()+theme_bw()

# make same plot with spprich data while we're here

ggplot(bio_jpn, aes(x = lat, y = tot_biom/totMsurv, colour=ThermalAffinity2)) + 
  geom_point()+geom_smooth(se=F)+facet_wrap(~FG, scales='free')+theme_bw()

ggplot(bio_aus, aes(x = Lat, y = tot_biom/totMsurv, colour=ThermalAffinity2)) + 
  geom_point()+geom_smooth(se=F)+facet_wrap(~FG, scales='free')+
  scale_x_reverse()+theme_bw()

# If Logging: 0.001 chosen as min constant as min corr biomass val
# in Aus is 0.0016 and second min in Japan is 0.00098

# ok biomass looks ok
bio_aus$cor_biom<-bio_aus$tot_biom/bio_aus$totMsurv
bio_jpn$cor_biom<-bio_jpn$tot_biom/bio_jpn$totMsurv

# setup tropical only standardisation

# Going to standardise using biomass of each FG in tropical site group
# check for biomass outliers in these site groups

ggplot(data=filter(bio_jpn, ThermalAffinity2=='tropical' ),
       aes(x=lat, y=cor_biom))+geom_point(aes(colour=SiteID))+
  geom_hline(data=filter(bio_jpn, ThermalAffinity2=='tropical' & lat<29)%>%
               group_by(FG)%>%summarise(mean_biom=mean(cor_biom, na.rm=T)),
             aes(yintercept = mean_biom))+
  geom_hline(data=filter(bio_jpn, ThermalAffinity2=='tropical' & lat<25.5)%>%
               group_by(FG)%>%summarise(mean_biom=mean(cor_biom, na.rm=T)),
             aes(yintercept = mean_biom), col='red')+
                facet_wrap(~FG, scales='free')
# go for Irimote, lat<25.5 option in Japan

ggplot(data=filter(bio_aus, ThermalAffinity2=='tropical'),
       aes(x=Lat, y=cor_biom))+geom_point(aes(colour=Site))+
  geom_hline(data=filter(bio_aus, ThermalAffinity2=='tropical' & Lat> -24.5)%>%
               group_by(FG)%>%summarise(mean_biom=mean(cor_biom, na.rm=T)),aes(yintercept = mean_biom))+
  facet_wrap(~FG, scales='free')+scale_x_reverse()

# Calc standardisation
jpn_trop_prop<-filter(bio_jpn, ThermalAffinity2=='tropical' & lat<25.5)%>%
  group_by(FG)%>%summarise(mean_biom=mean(cor_biom, na.rm=T))

aus_trop_prop<-filter(bio_aus, ThermalAffinity2=='tropical' & Lat> -24.5)%>%
  group_by(FG)%>%summarise(mean_biom=mean(cor_biom, na.rm=T))

jpn_trop_prop<-left_join(filter(bio_jpn, ThermalAffinity2=='tropical'),
                         jpn_trop_prop, by='FG')

aus_trop_prop<-left_join(filter(bio_aus, ThermalAffinity2=='tropical'),
                         aus_trop_prop, by='FG')

jpn_trop_comm<-filter(bio_jpn, ThermalAffinity2=='tropical')%>%
  group_by(lat, SiteID, Site.trans.ID)%>%summarise(cor_biom=sum(cor_biom))
jpn_trop_comm$mean_biom<-as.numeric(filter(jpn_trop_comm, lat<25.5)%>%ungroup()%>%
  summarise(mean_biom=mean(cor_biom, na.rm=T)))

aus_trop_comm<-filter(bio_aus, ThermalAffinity2=='tropical')%>%
  group_by(Lat, Site, Site.trans.ID)%>%summarise(cor_biom=sum(cor_biom))
aus_trop_comm$mean_biom<-as.numeric(filter(aus_trop_comm, Lat> -24.5)%>%ungroup()%>%
                                     summarise(mean_biom=mean(cor_biom, na.rm=T)))

# visual check to make sure things look right
ggplot(data=rbind(data.frame(jpn_trop_comm, FG=99), jpn_trop_prop[c(7,1,9,10,2)]),
       aes(x=lat, y=cor_biom))+geom_point()+geom_smooth()+facet_wrap(~FG, scales='free')

ggplot(data=rbind(data.frame(jpn_trop_comm, FG=99), jpn_trop_prop[c(6,8,9,2)]),
       aes(x=lat, y=cor_biom/max_biom))+geom_point()+geom_smooth()+facet_wrap(~FG, scales='free')


# biomass
outlz<-which(jpn_trop_prop$cor_biom/jpn_trop_prop$mean_biom>10)

ggplot(jpn_trop_prop[-outlz,], aes(x = lat, y = cor_biom/mean_biom)) + 
  geom_point(aes(colour=factor(FG)))+
  geom_smooth(aes(colour=factor(FG)),se=F)+
  geom_smooth(data=jpn_trop_comm, se=F, colour='black', linetype='dashed')+
  scale_x_continuous(breaks=24:35)+theme_bw()+
  facet_wrap(~FG, scales='free')+theme(legend.position = "none")

outlz<-which(aus_trop_prop$cor_biom/aus_trop_prop$mean_biom>10)

ggplot(aus_trop_prop[-outlz,], aes(x = Lat, y = cor_biom/mean_biom)) + 
  geom_point(aes(colour=factor(FG)))+
  geom_smooth(aes(colour=factor(FG)),se=F)+
  geom_smooth(data=aus_trop_comm, se=F, colour='black', linetype='dashed')+
 scale_x_reverse(breaks=-23:-33)+theme_bw()+
  facet_wrap(~FG, scales='free')+theme(legend.position = "none")
# might need more wiggliness for Aussie GAM

# calc statistics 
#1) make gams per FG and for comm to visualise tropicalization
#2) make mixed anovas to compare FG tropicalization 'levels' within each site-group 

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
jpn_trop_tests[jpn_trop_tests$SiteID %in% c('JP27', 'JP32'),]$site.group<-'trans.bay'
jpn_trop_tests[jpn_trop_tests$lat > 34,]$site.group<-'temp.headld'
table(jpn_trop_tests$site.group, jpn_trop_tests$SiteID)

aus_trop_comm$FG<-'comm'
aus_trop_prop$FG<-as.character(aus_trop_prop$FG)
aus_trop_tests<-rbind(data.frame(aus_trop_comm), aus_trop_prop[names(aus_trop_comm)])
aus_trop_tests$site.group<-'trop.base'
aus_trop_tests[aus_trop_tests$Lat > -25.6 &  aus_trop_tests$Lat < -24.5,]$site.group<-'trans.bay'
aus_trop_tests[aus_trop_tests$Lat > -28 &  aus_trop_tests$Lat < -25.6,]$site.group<-'trans.offshore'
aus_trop_tests[aus_trop_tests$Lat < -28,]$site.group<-'temp.offshore'
aus_trop_tests[aus_trop_tests$Site %in% c('Muttonbird Island', 'Woolgoolga Reef', 'Woolgoolga Headland', 'North Rock'),]$site.group<-'temp.inshore'

table(aus_trop_tests$site.group, aus_trop_tests$Site)

#### run FG tropicalization comparisons ####

#### Japan FG tropicalization comps ####

# drop some FGs and set comm as intercept
jpn_trop_tests<-jpn_trop_tests[-which(jpn_trop_tests$FG %in% c(17,18,19,5,7,11)),]
jpn_trop_tests$FG<-factor(jpn_trop_tests$FG, levels=c('comm', 15, 10, 8, 2,6,12,4,1,16,13,14,3,9))

## trop.island

ggplot(data=filter(jpn_trop_tests, site.group=='trop.island'), aes(x=FG, y=trop_met))+
  geom_point(aes(colour=SiteID))+geom_boxplot(alpha=0.5) # remember boxplot = medians

m1<-lme4::lmer(trop_met~FG+(1|SiteID), data=filter(jpn_trop_tests, site.group=='trop.island'))
resid_panel(m1)
summary(m1)

mod.rg <- update(ref_grid(m1), tran = make.tran("power", 0.25))
em1<-emmeans(mod.rg, specs='FG', type='response')
pairs(em1)
plot(em1, comparisons = TRUE)

trop_comps_out<-data.frame(site.group='trop.island',data.frame(em1), rbind(c(NA, NA), data.frame(pairs(em1))[1:13, c(1,6)]))

ggplot()+
  geom_jitter(data=filter(jpn_trop_tests, site.group=='trop.island'), aes(x=FG, y=trop_met^4), shape=1, alpha=0.5, width=0.2, colour='grey')+
  geom_hline(yintercept = trop_comps_out$response[1], linetype='dotted')+
  geom_errorbar(data=trop_comps_out, aes(x=FG, ymin=lower.CL, ymax=upper.CL))+
  geom_point(data=trop_comps_out, aes(x=FG, y=response, colour=ifelse(p.value>0.05| is.na(p.value), 'blue', 'red')), size=2)+
 theme_bw()+theme(legend.position = "none")

## trans.island

ggplot(data=filter(jpn_trop_tests, site.group=='trans.island'), aes(x=FG, y=trop_met))+
  geom_point(aes(colour=SiteID))+geom_boxplot(alpha=0.5) # remember boxplot = medians

m1<-lme4::lmer(trop_met~FG+(1|SiteID), data=filter(jpn_trop_tests, site.group=='trans.island'))
resid_panel(m1)
summary(m1)

mod.rg <- update(ref_grid(m1), tran = make.tran("power", 0.25))
em1<-emmeans(mod.rg, specs='FG', type='response')
pairs(em1)
plot(em1, comparisons = TRUE)

trop_comps_out<-rbind(trop_comps_out,
                      data.frame(site.group='trans.island',data.frame(em1), rbind(c(NA, NA), data.frame(pairs(em1))[1:13, c(1,6)])))

ggplot()+
  geom_jitter(data=filter(jpn_trop_tests, site.group=='trans.island'), aes(x=FG, y=trop_met^4), shape=1, alpha=0.5, width=0.2, colour='grey')+
  geom_hline(yintercept =trop_comps_out[trop_comps_out$site.group=='trans.island',]$response[1], linetype='dotted')+
  geom_errorbar(data=filter(trop_comps_out, site.group=='trans.island'), aes(x=FG, ymin=lower.CL, ymax=upper.CL))+
  geom_point(data=filter(trop_comps_out, site.group=='trans.island'), aes(x=FG, y=response, colour=ifelse(p.value>0.05| is.na(p.value), 'blue', 'red')), size=2)+
  ylim(c(0,5))+theme_bw()+theme(legend.position = "none")

## trans.inland

ggplot(data=filter(jpn_trop_tests, site.group=='trans.inland'), aes(x=FG, y=trop_met))+
  geom_point(aes(colour=SiteID))+geom_boxplot(alpha=0.5) # remember boxplot = medians

zer_fgs<-filter(jpn_trop_tests, site.group=='trans.inland')%>%group_by(FG)%>%summarise(st=sum(trop_met))%>%filter(., st==0)%>%.$FG


m1<-lme4::lmer(trop_met~FG+(1|SiteID), data=filter(jpn_trop_tests, site.group=='trans.inland' & !FG%in% zer_fgs))
resid_panel(m1)
summary(m1)

mod.rg <- update(ref_grid(m1), tran = make.tran("power", 0.25))
em1<-emmeans(mod.rg, specs='FG', type='response')
pairs(em1)
plot(em1, comparisons = TRUE)

trop_comps_out<-rbind(trop_comps_out,
                      data.frame(site.group='trans.inland',data.frame(em1), rbind(c(NA, NA), data.frame(pairs(em1))[1:(13-length(zer_fgs)), c(1,6)])))

ggplot()+
  geom_jitter(data=filter(jpn_trop_tests, site.group=='trans.inland'), aes(x=FG, y=trop_met^4), shape=1, alpha=0.5, width=0.2, colour='grey')+
  geom_hline(yintercept =trop_comps_out[trop_comps_out$site.group=='trans.inland',]$response[1], linetype='dotted')+
  geom_errorbar(data=filter(trop_comps_out, site.group=='trans.inland'), aes(x=FG, ymin=lower.CL, ymax=upper.CL))+
  geom_point(data=filter(trop_comps_out, site.group=='trans.inland'), aes(x=FG, y=response, colour=ifelse(p.value>0.05| is.na(p.value), 'blue', 'red')), size=2)+
  theme_bw()+theme(legend.position = "none")


## trans.headld

ggplot(data=filter(jpn_trop_tests, site.group=='trans.headld'), aes(x=FG, y=trop_met))+
  geom_point(aes(colour=SiteID))+geom_boxplot(alpha=0.5) # remember boxplot = medians

zer_fgs<-filter(jpn_trop_tests, site.group=='trans.headld')%>%group_by(FG)%>%summarise(st=sum(trop_met))%>%filter(., st==0)%>%.$FG


m1<-lme4::lmer(trop_met~FG+(1|SiteID), data=filter(jpn_trop_tests, site.group=='trans.headld'& !FG%in% zer_fgs))
resid_panel(m1)
summary(m1)

mod.rg <- update(ref_grid(m1), tran = make.tran("power", 0.25))
em1<-emmeans(mod.rg, specs='FG', type='response')
pairs(em1)
plot(em1, comparisons = TRUE)

trop_comps_out<-rbind(trop_comps_out,
                      data.frame(site.group='trans.headld',data.frame(em1), rbind(c(NA, NA),
                                                                                  data.frame(pairs(em1))[1:(13-length(zer_fgs)), c(1,6)])))

ggplot()+
  geom_jitter(data=filter(jpn_trop_tests, site.group=='trans.headld'), aes(x=FG, y=trop_met^4), shape=1, alpha=0.5, width=0.2, colour='grey')+
  geom_hline(yintercept =trop_comps_out[trop_comps_out$site.group=='trans.headld',]$response[1], linetype='dotted')+
  geom_errorbar(data=filter(trop_comps_out, site.group=='trans.headld'), aes(x=FG, ymin=lower.CL, ymax=upper.CL))+
  geom_point(data=filter(trop_comps_out, site.group=='trans.headld'), aes(x=FG, y=response, colour=ifelse(p.value>0.05| is.na(p.value), 'blue', 'red')), size=2)+
  theme_bw()+ylim(c(0,5))+theme(legend.position = "none")

## trans.bay

ggplot(data=filter(jpn_trop_tests, site.group=='trans.bay'), aes(x=FG, y=trop_met))+
  geom_point(aes(colour=SiteID))+geom_boxplot(alpha=0.5) # remember boxplot = medians
# remove FGs with only 0s - these will be sig diff 
zer_fgs<-filter(jpn_trop_tests, site.group=='trans.bay')%>%group_by(FG)%>%summarise(st=sum(trop_met))%>%filter(., st==0)%>%.$FG

m1<-lme4::lmer(trop_met~FG+(1|SiteID), data=filter(jpn_trop_tests, site.group=='trans.bay' & !FG%in% zer_fgs))
resid_panel(m1)
summary(m1)

mod.rg <- update(ref_grid(m1), tran = make.tran("power", 0.25))
em1<-emmeans(mod.rg, specs='FG', type='response')
pairs(em1)
plot(em1, comparisons = TRUE)

trop_comps_out<-rbind(trop_comps_out,
                      data.frame(site.group='trans.bay',data.frame(em1), rbind(c(NA, NA), data.frame(pairs(em1))[1:(13-length(zer_fgs)), c(1,6)])))

ggplot()+
  geom_jitter(data=filter(jpn_trop_tests, site.group=='trans.bay'), aes(x=FG, y=trop_met^4), shape=1, alpha=0.5, width=0.2, colour='grey')+
  geom_hline(yintercept =trop_comps_out[trop_comps_out$site.group=='trans.bay',]$response[1], linetype='dotted')+
  geom_errorbar(data=filter(trop_comps_out, site.group=='trans.bay'), aes(x=FG, ymin=lower.CL, ymax=upper.CL))+
  geom_point(data=filter(trop_comps_out, site.group=='trans.bay'), aes(x=FG, y=response, colour=ifelse(p.value>0.05| is.na(p.value), 'blue', 'red')), size=2)+
  theme_bw()+ylim(c(0,1))+theme(legend.position = "none")

## temp.headld

ggplot(data=filter(jpn_trop_tests, site.group=='temp.headld'), aes(x=FG, y=trop_met))+
  geom_point(aes(colour=SiteID))+geom_boxplot(alpha=0.5) # remember boxplot = medians
# remove FGs with only 0s - these will be sig diff 
zer_fgs<-filter(jpn_trop_tests, site.group=='temp.headld')%>%group_by(FG)%>%summarise(st=sum(trop_met))%>%filter(., st==0)%>%.$FG

m1<-lme4::lmer(trop_met~FG+(1|SiteID), data=filter(jpn_trop_tests, site.group=='temp.headld' & !FG%in% zer_fgs))
resid_panel(m1)
summary(m1)

mod.rg <- update(ref_grid(m1), tran = make.tran("power", 0.25))
em1<-emmeans(mod.rg, specs='FG', type='response')
pairs(em1)
plot(em1, comparisons = TRUE)

trop_comps_out<-rbind(trop_comps_out,
                      data.frame(site.group='temp.headld',data.frame(em1), rbind(c(NA, NA), data.frame(pairs(em1))[1:(13-length(zer_fgs)), c(1,6)])))

ggplot()+
  geom_jitter(data=filter(jpn_trop_tests, site.group=='temp.headld'), aes(x=FG, y=trop_met^4), shape=1, alpha=0.5, width=0.2, colour='grey')+
  geom_hline(yintercept =trop_comps_out[trop_comps_out$site.group=='temp.headld',]$response[1], linetype='dotted')+
  geom_errorbar(data=filter(trop_comps_out, site.group=='temp.headld'), aes(x=FG, ymin=lower.CL, ymax=upper.CL))+
  geom_point(data=filter(trop_comps_out, site.group=='temp.headld'), aes(x=FG, y=response, colour=ifelse(p.value>0.05| is.na(p.value), 'blue', 'red')), size=2)+
  theme_bw()+ylim(c(0,2))+theme(legend.position = "none")

#all plot

trop_comps_out$site.group<-factor(trop_comps_out$site.group, levels=c('trop.base', 'trop.island', 'trans.island', 'trans.inland',
                                                                      'trans.headld', 'trans.bay', 'temp.headld'))

ggplot()+
  geom_hline(yintercept =filter(trop_comps_out, FG=='comm')$response, linetype='dotted')+
  geom_errorbar(data=trop_comps_out, aes(x=FG, ymin=lower.CL, ymax=upper.CL))+
  geom_point(data=trop_comps_out, aes(x=FG, y=response, colour=ifelse(p.value>0.05| is.na(p.value), 'blue', 'red')), size=2)+
  theme_bw()+theme(legend.position = "none")+facet_wrap(~site.group, scales='free')

#### Australia FG tropicalization comparisons #### 

# drop some FGs and set comm as intercept
aus_trop_tests<-aus_trop_tests[-which(aus_trop_tests$FG %in% c(17,18,19,11)),]
aus_trop_tests$FG<-factor(aus_trop_tests$FG, levels=c('comm', 15, 10, 8, 2,6,12,4,5,1,16,13,14,3,7,9))

## trans.bay

ggplot(data=filter(aus_trop_tests, site.group=='trans.bay'), aes(x=FG, y=trop_met))+
  geom_point(aes(colour=Site))+geom_boxplot(alpha=0.5) # remember boxplot = medians

zer_fgs<-filter(aus_trop_tests, site.group=='trans.bay')%>%group_by(FG)%>%summarise(st=sum(trop_met))%>%filter(., st==0)%>%.$FG


m1<-lme4::lmer(trop_met~FG+(1|Site), data=filter(aus_trop_tests, site.group=='trans.bay' & !FG%in% zer_fgs))
resid_panel(m1)
summary(m1)

mod.rg <- update(ref_grid(m1), tran = make.tran("power", 0.25))
em1<-emmeans(mod.rg, specs='FG', type='response')
#pairs(em1)

trop_comps_out_aus<-data.frame(site.group='trans.bay',data.frame(em1), rbind(c(NA, NA), data.frame(pairs(em1))[1:(15-length(zer_fgs)), c(1,6)]))

ggplot()+
  geom_jitter(data=filter(aus_trop_tests, site.group=='trans.bay'), aes(x=FG, y=trop_met^4), shape=1, alpha=0.5, width=0.2, colour='grey')+
  geom_hline(yintercept =trop_comps_out_aus[trop_comps_out_aus$site.group=='trans.bay',]$response[1], linetype='dotted')+
  geom_errorbar(data=filter(trop_comps_out_aus, site.group=='trans.bay'), aes(x=FG, ymin=lower.CL, ymax=upper.CL))+
  geom_point(data=filter(trop_comps_out_aus, site.group=='trans.bay'), aes(x=FG, y=response, colour=ifelse(p.value>0.05| is.na(p.value), 'blue', 'red')), size=2)+
  theme_bw()+theme(legend.position = "none")+ylim(c(0,2))


## trans.offshore

ggplot(data=filter(aus_trop_tests, site.group=='trans.offshore'), aes(x=FG, y=trop_met))+
  geom_point(aes(colour=Site))+geom_boxplot(alpha=0.5) # remember boxplot = medians

zer_fgs<-filter(aus_trop_tests, site.group=='trans.offshore')%>%group_by(FG)%>%summarise(st=sum(trop_met))%>%filter(., st==0)%>%.$FG


m1<-lme4::lmer(trop_met~FG+(1|Site), data=filter(aus_trop_tests, site.group=='trans.offshore'& !FG%in% zer_fgs))
resid_panel(m1)
summary(m1)

mod.rg <- update(ref_grid(m1), tran = make.tran("power", 0.25))
em1<-emmeans(mod.rg, specs='FG', type='response')
#pairs(em1)
#plot(em1, comparisons = TRUE)

trop_comps_out_aus<-rbind(trop_comps_out_aus,
                      data.frame(site.group='trans.offshore',data.frame(em1), rbind(c(NA, NA),
                                                                                  data.frame(pairs(em1))[1:(15-length(zer_fgs)), c(1,6)])))
ggplot()+
  geom_jitter(data=filter(aus_trop_tests, site.group=='trans.offshore'), aes(x=FG, y=trop_met^4), shape=1, alpha=0.5, width=0.2, colour='grey')+
  geom_hline(yintercept =trop_comps_out_aus[trop_comps_out_aus$site.group=='trans.offshore',]$response[1], linetype='dotted')+
  geom_errorbar(data=filter(trop_comps_out_aus, site.group=='trans.offshore'), aes(x=FG, ymin=lower.CL, ymax=upper.CL))+
  geom_point(data=filter(trop_comps_out_aus, site.group=='trans.offshore'), aes(x=FG, y=response, colour=ifelse(p.value>0.05| is.na(p.value), 'blue', 'red')), size=2)+
  theme_bw()+ylim(c(0,5))+theme(legend.position = "none")

## temp.offshore

ggplot(data=filter(aus_trop_tests, site.group=='temp.offshore'), aes(x=FG, y=trop_met))+
  geom_point(aes(colour=Site))+geom_boxplot(alpha=0.5) # remember boxplot = medians
# remove FGs with only 0s - these will be sig diff 
zer_fgs<-filter(aus_trop_tests, site.group=='temp.offshore')%>%group_by(FG)%>%summarise(st=sum(trop_met))%>%filter(., st==0)%>%.$FG

m1<-lme4::lmer(trop_met~FG+(1|Site), data=filter(aus_trop_tests, site.group=='temp.offshore' & !FG%in% zer_fgs))
resid_panel(m1)
summary(m1)

mod.rg <- update(ref_grid(m1), tran = make.tran("power", 0.25))
em1<-emmeans(mod.rg, specs='FG', type='response')
#pairs(em1)
#plot(em1, comparisons = TRUE)

trop_comps_out_aus<-rbind(trop_comps_out_aus,
                      data.frame(site.group='temp.offshore',data.frame(em1), rbind(c(NA, NA), data.frame(pairs(em1))[1:(15-length(zer_fgs)), c(1,6)])))

ggplot()+
  geom_jitter(data=filter(aus_trop_tests, site.group=='temp.offshore'), aes(x=FG, y=trop_met^4), shape=1, alpha=0.5, width=0.2, colour='grey')+
  geom_hline(yintercept =trop_comps_out_aus[trop_comps_out_aus$site.group=='temp.offshore',]$response[1], linetype='dotted')+
  geom_errorbar(data=filter(trop_comps_out_aus, site.group=='temp.offshore'), aes(x=FG, ymin=lower.CL, ymax=upper.CL))+
  geom_point(data=filter(trop_comps_out_aus, site.group=='temp.offshore'), aes(x=FG, y=response, colour=ifelse(p.value>0.05| is.na(p.value), 'blue', 'red')), size=2)+
  theme_bw()+ylim(c(0,1))+theme(legend.position = "none")

## temp.inshore

ggplot(data=filter(aus_trop_tests, site.group=='temp.inshore'), aes(x=FG, y=trop_met))+
  geom_point(aes(colour=Site))+geom_boxplot(alpha=0.5) # remember boxplot = medians
# remove FGs with only 0s - these will be sig diff 
zer_fgs<-filter(aus_trop_tests, site.group=='temp.inshore')%>%group_by(FG)%>%summarise(st=sum(trop_met))%>%filter(., st==0)%>%.$FG

m1<-lme4::lmer(trop_met~FG+(1|Site), data=filter(aus_trop_tests, site.group=='temp.inshore' & !FG%in% zer_fgs))
resid_panel(m1)
summary(m1)

mod.rg <- update(ref_grid(m1), tran = make.tran("power", 0.25))
em1<-emmeans(mod.rg, specs='FG', type='response')
#pairs(em1)
#plot(em1, comparisons = TRUE)

trop_comps_out_aus<-rbind(trop_comps_out_aus,
                      data.frame(site.group='temp.inshore',data.frame(em1), rbind(c(NA, NA), data.frame(pairs(em1))[1:(15-length(zer_fgs)), c(1,6)])))

ggplot()+
  geom_jitter(data=filter(aus_trop_tests, site.group=='temp.inshore'), aes(x=FG, y=trop_met^4), shape=1, alpha=0.5, width=0.2, colour='grey')+
  geom_hline(yintercept =trop_comps_out_aus[trop_comps_out_aus$site.group=='temp.inshore',]$response[1], linetype='dotted')+
  geom_errorbar(data=filter(trop_comps_out_aus, site.group=='temp.inshore'), aes(x=FG, ymin=lower.CL, ymax=upper.CL))+
  geom_point(data=filter(trop_comps_out_aus, site.group=='temp.inshore'), aes(x=FG, y=response, colour=ifelse(p.value>0.05| is.na(p.value), 'blue', 'red')), size=2)+
  theme_bw()+ylim(c(0,2))+theme(legend.position = "none")

#all plot

trop_comps_out_aus$site.group<-factor(trop_comps_out_aus$site.group, levels=c('trop.base', 'trans.bay', 'trans.offshore',
                                                                              'temp.offshore', 'temp.inshore'))

ggplot()+
  geom_hline(yintercept =filter(trop_comps_out_aus, FG=='comm')$response, linetype='dotted')+
  geom_errorbar(data=trop_comps_out_aus, aes(x=FG, ymin=lower.CL, ymax=upper.CL))+
  geom_point(data=trop_comps_out_aus, aes(x=FG, y=response, colour=ifelse(p.value>0.05| is.na(p.value), 'blue', 'red')), size=2)+
  theme_bw()+theme(legend.position = "none")+facet_wrap(~site.group, scales='free')


# 1) GAMS

# Japan all FGs
jpn_comm_m2<-gam(trop_met~s(lat), data=jpn_trop_comm, method = 'REML')
par(mfrow=c(2,2))
gam.check(jpn_comm_m2)
plot(jpn_comm_m2, residuals = T, pch=1, cex=1, rug=T,
     shade=T, seWithMean = T, shift = coef(jpn_comm_m2)[1],
     xlab='Latitude', ylab='delta Biomass relative to tropical site')
jpn_comm_m2$sp
jpn_comm_m3<-gam(trop_met~s(lat, k=7)+s(SiteID, bs='re'), data=jpn_trop_comm, method = 'REML')
par(mfrow=c(2,2));gam.check(jpn_comm_m3)
plot(jpn_comm_m3, residuals = T, pch=1, cex=1, rug=T,
     shade=T, seWithMean = T, shift = coef(jpn_comm_m3)[1],
     xlab='Latitude', ylab='delta Biomass relative to tropical site')

jpn_gam_pred<-expand.grid(lat=seq(24.3, 35, 0.1), FG='all')
jpn_gam_pred<-cbind(jpn_gam_pred,predict.gam(jpn_comm_m3, newdata =jpn_gam_pred, 
                     exclude='s(SiteID)', newdata.guaranteed = T, type='link', se.fit=T ))

qplot(data=jpn_gam_pred, x=lat, y=fit^4, geom='line')+
  geom_point(data=jpn_trop_comm, aes(x=lat, y=trop_met^4))

# Australia all FGs
# exclude 1 outlier
aus_comm_m3<-gam(trop_met~s(Lat, k=5 )+s(Site, bs='re'),data=aus_trop_comm, method = 'REML')
summary(aus_comm_m3)
par(mfrow=c(2,2));gam.check(aus_comm_m3)
plot(aus_comm_m3, residuals = T, pch=1, cex=1, rug=T,
     shade=T, seWithMean = T, shift = coef(aus_comm_m3)[1],
     xlab='Latitude', ylab='delta Biomass relative to tropical site')

aus_gam_pred<-expand.grid(Lat=seq(-31, -23.4, 0.1), FG='all')
aus_gam_pred<-cbind(aus_gam_pred,predict.gam(aus_comm_m3, newdata =aus_gam_pred,
                    type='link', se.fit=T, exclude='s(Site)', newdata.guaranteed = T ))

qplot(data=aus_gam_pred, x=Lat, y=fit^4, geom='line')+
  geom_point(data=aus_trop_comm, aes(x=Lat, y=trop_met^4))

#Individual FGs modelled in one GAM

#Japan

jpn_trop_prop2<-filter(jpn_trop_prop, FG %in% c(1,2,4,6,9,10,12,15,16)) # drop some FGs
jpn_trop_prop2$FG<-factor(jpn_trop_prop2$FG)

jpn_trop_prop2[jpn_trop_prop2$lat<25.5,]$trop_met<-1 # set tropical group-site
# to 1 rather than divide by group-site mean = forces intercept thru 1

jpn_fg_m1<-gam(trop_met~s(lat, by=FG, k=7)+FG +s(SiteID, bs='re'), 
               data=jpn_trop_prop2[jpn_trop_prop2$trop_met<10^0.25,], method = 'REML')
summary(jpn_fg_m1)
coef(jpn_fg_m1)
gam.check(jpn_fg_m1)
plot(jpn_fg_m1, residuals = F, rug=T, pages=1, all.terms = T)

jpn_gam_pred2<-expand.grid(lat=seq(24.3, 35, 0.1), FG=c(1,2,4, 6, 9, 10, 12, 15, 16))
jpn_gam_pred2<-cbind(jpn_gam_pred2,predict.gam(jpn_fg_m1, newdata =jpn_gam_pred2,
                                               type='link', se.fit=T,exclude='s(SiteID)', newdata.guaranteed = T ))

jpn_gam_pred2$FG<-as.character(jpn_gam_pred2$FG)
jpn_gam_pred_all<-rbind(jpn_gam_pred, jpn_gam_pred2)

ggplot(jpn_gam_pred_all[jpn_gam_pred_all$FG!='all',], aes(x = lat, y = fit^4,
                                                          ymax=(fit+se.fit*1.96)^4, ymin=(fit-se.fit*1.96)^4))+  
  geom_ribbon(data=jpn_gam_pred_all[jpn_gam_pred_all$FG=='all',]%>%rename(FG2=FG), fill='yellow', alpha=0.5)+
  geom_line(data=jpn_gam_pred_all[jpn_gam_pred_all$FG=='all',]%>%rename(FG2=FG))+
  geom_ribbon(fill='red', alpha=0.5)+geom_line(colour='red')+
  scale_x_continuous(breaks=24:35)+theme_bw()+
  facet_wrap(~FG, scales='free')

# Australia

aus_trop_prop2<-filter(aus_trop_prop, FG %in% c(1,2,4,6,9,10,12,15,16)) # drop some FGs
aus_trop_prop2$FG<-factor(aus_trop_prop2$FG)

aus_trop_prop2[aus_trop_prop2$Lat> -24.5,]$trop_met<-1 # set tropical group-site
# to 1 rather than divide by group-site mean = forces intercept thru 1


aus_fg_m1<-gam(trop_met~s(Lat, by=FG, k=5)+FG +s(Site, bs='re'), 
               data=aus_trop_prop2[aus_trop_prop2$trop_met<10^0.25,], method = 'REML')
summary(aus_fg_m1)
coef(aus_fg_m1)
gam.check(aus_fg_m1)
plot(aus_fg_m1, residuals = F, rug=T, pages=1, all.terms = T)

aus_gam_pred2<-expand.grid(Lat=seq(-31, -23, 0.1), FG=c(1,2,4, 6, 9, 10, 12, 15, 16))
aus_gam_pred2<-cbind(aus_gam_pred2,predict.gam(aus_fg_m1, newdata =aus_gam_pred2,
                                               type='link', se.fit=T,exclude='s(Site)', newdata.guaranteed = T ))

aus_gam_pred2$FG<-as.character(aus_gam_pred2$FG)
aus_gam_pred_all<-rbind(aus_gam_pred, aus_gam_pred2)

ggplot(aus_gam_pred_all[aus_gam_pred_all$FG!='all',], aes(x = Lat, y = fit^4,
                                                          ymax=(fit+se.fit*1.96)^4, ymin=(fit-se.fit*1.96)^4))+  
  geom_ribbon(data=aus_gam_pred_all[aus_gam_pred_all$FG=='all',]%>%rename(FG2=FG), fill='yellow', alpha=0.5)+
  geom_line(data=aus_gam_pred_all[aus_gam_pred_all$FG=='all',]%>%rename(FG2=FG))+
  geom_ribbon(fill='red', alpha=0.5)+geom_line(colour='red')+
  theme_bw()+scale_x_reverse()+
  facet_wrap(~FG, scales='free')



# Japan indivdiual FGs

# FG1
jpn.fg1<-gam(trop_met~s(lat, k=7)+s(SiteID, bs='re'), data=filter(jpn_trop_prop, FG==1), method = 'REML')
modl<-jpn.fg1
summary(modl)
par(mfrow=c(2,2));gam.check(modl)
plot(modl, residuals = T, pch=1, cex=1, rug=T,
     shade=T, seWithMean = T, shift = coef(modl)[1],xlab='Latitude', ylab='delta Biomass relative to tropical site')
qplot(data=data.frame(lat=seq(24.3, 35, 0.1),predict.gam(modl, newdata =data.frame(lat=seq(24.3, 35, 0.1)),type='link', se.fit=T, exclude='s(SiteID)', newdata.guaranteed = T )),x=lat, y=fit^4, geom='line')+
  geom_jitter(data=filter(jpn_trop_prop, FG==1),aes(x=lat, y=trop_met^4), shape=1, height=0.05,width=0.1)

jpn.fg2<-gam(trop_met~s(lat, k=7)+s(SiteID, bs='re'), data=filter(jpn_trop_prop, FG==2), method = 'REML')
modl<-jpn.fg2
summary(modl)
par(mfrow=c(2,2));gam.check(modl)
plot(modl, residuals = T, pch=1, cex=1, rug=T,
                       shade=T, seWithMean = T, shift = coef(modl)[1],xlab='Latitude', ylab='delta Biomass relative to tropical site')
qplot(data=data.frame(lat=seq(24.3, 35, 0.1),predict.gam(modl, newdata =data.frame(lat=seq(24.3, 35, 0.1)),type='link', se.fit=T, exclude='s(SiteID)', newdata.guaranteed = T )),x=lat, y=fit^4, geom='line')+
  geom_jitter(data=filter(jpn_trop_prop, FG==2),aes(x=lat, y=trop_met^4), shape=1, height=0.05,width=0.1)

jpn.fg3<-gam(trop_met~s(lat, k=5, sp=0.1)+s(SiteID, bs='re'), data=filter(jpn_trop_prop, FG==3), method = 'REML')
modl<-jpn.fg2
summary(modl)
par(mfrow=c(2,2));gam.check(modl)
plot(modl, residuals = T, pch=1, cex=1, rug=T,
     shade=T, seWithMean = T, shift = coef(modl)[1],xlab='Latitude', ylab='delta Biomass relative to tropical site')
qplot(data=data.frame(lat=seq(24.3, 35, 0.1),predict.gam(modl, newdata =data.frame(lat=seq(24.3, 35, 0.1)),type='link', se.fit=T, exclude='s(SiteID)', newdata.guaranteed = T )),x=lat, y=fit^4, geom='line')+
  geom_jitter(data=filter(jpn_trop_prop, FG==3),aes(x=lat, y=trop_met^4), shape=1, height=0.05,width=0.1)

# Aus individual FGs

# FG1
aus.fg1<-gam(trop_met~s(Lat, k=5, sp=0.001)+s(Site, bs='re'),
             data=filter(aus_trop_prop2, FG==1 & trop_met< 10^0.25), method = 'REML')
modl<-aus.fg1
summary(modl)
par(mfrow=c(2,2));gam.check(modl)
plot(modl, residuals = T, pch=1, cex=1, rug=T,
     shade=T, seWithMean = T, shift = coef(modl)[1],xlab='Latitude', ylab='delta Biomass relative to tropical site')
qplot(data=data.frame(Lat=seq(-31, -23, 0.1),predict.gam(modl, newdata =data.frame(Lat=seq(-31, -23, 0.1)),type='link', se.fit=T, exclude='s(Site)', newdata.guaranteed = T )),x=Lat, y=fit^4, geom='line')+
  geom_point(data=filter(aus_trop_prop, FG==1),aes(x=Lat, y=trop_met^4), shape=1)+ylim(c(0,10))

# FG2

aus.fg2<-gam(trop_met~s(Lat, k=5, sp=0.001)+s(Site, bs='re'),
             data=filter(aus_trop_prop2, FG==2 & trop_met< 10^0.25), method = 'REML')
modl<-aus.fg2
summary(modl)
par(mfrow=c(2,2));gam.check(modl)
plot(modl, residuals = T, pch=1, cex=1, rug=T,
     shade=T, seWithMean = T, shift = coef(modl)[1],xlab='Latitude', ylab='delta Biomass relative to tropical site')
qplot(data=data.frame(Lat=seq(-31, -23, 0.1),predict.gam(modl, newdata =data.frame(Lat=seq(-31, -23, 0.1)),type='link', se.fit=T, exclude='s(Site)', newdata.guaranteed = T )),x=Lat, y=fit^4, geom='line')+
  geom_point(data=filter(aus_trop_prop, FG==2),aes(x=Lat, y=trop_met^4), shape=1)

# FG4

aus.fg4<-gam(trop_met~s(Lat, k=5)+s(Site, bs='re'),
             data=filter(aus_trop_prop2, FG==4 & trop_met< 10^0.25), method = 'REML') # sp=2730
modl<-aus.fg4
summary(modl)
par(mfrow=c(2,2));gam.check(modl)
plot(modl, residuals = T, pch=1, cex=1, rug=T,
     shade=T, seWithMean = T, shift = coef(modl)[1],xlab='Latitude', ylab='delta Biomass relative to tropical site')
qplot(data=data.frame(Lat=seq(-31, -23, 0.1),predict.gam(modl, newdata =data.frame(Lat=seq(-31, -23, 0.1)),type='link', se.fit=T, exclude='s(Site)', newdata.guaranteed = T )),x=Lat, y=fit^4, geom='line')+
  geom_point(data=filter(aus_trop_prop, FG==4),aes(x=Lat, y=trop_met^4), shape=1)+ylim(c(0,10))+geom_hline(yintercept=1, colour='red')


# FG6

aus.fg6<-gam(trop_met~s(Lat, k=5)+s(Site, bs='re'),
             data=filter(aus_trop_prop2, FG==6 & trop_met< 10^0.25), method = 'REML') # sp=0.0019
modl<-aus.fg6
summary(modl)
par(mfrow=c(2,2));gam.check(modl)
plot(modl, residuals = T, pch=1, cex=1, rug=T,
     shade=T, seWithMean = T, shift = coef(modl)[1],xlab='Latitude', ylab='delta Biomass relative to tropical site')
qplot(data=data.frame(Lat=seq(-31, -23, 0.1),predict.gam(modl, newdata =data.frame(Lat=seq(-31, -23, 0.1)),type='link', se.fit=T, exclude='s(Site)', newdata.guaranteed = T )),x=Lat, y=fit^4, geom='line')+
  geom_point(data=filter(aus_trop_prop, FG==6),aes(x=Lat, y=trop_met^4), shape=1)+geom_hline(yintercept=1, colour='red')

# FG9

aus.fg9<-gam(trop_met~s(Lat, k=5)+s(Site, bs='re'),
             data=filter(aus_trop_prop2, FG==9 & trop_met< 10^0.25), method = 'REML') # sp=0.0019
modl<-aus.fg9
summary(modl)
par(mfrow=c(2,2));gam.check(modl)
plot(modl, residuals = T, pch=1, cex=1, rug=T,
     shade=T, seWithMean = T, shift = coef(modl)[1],xlab='Latitude', ylab='delta Biomass relative to tropical site')
qplot(data=data.frame(Lat=seq(-31, -23, 0.1),predict.gam(modl, newdata =data.frame(Lat=seq(-31, -23, 0.1)),type='link', se.fit=T, exclude='s(Site)', newdata.guaranteed = T )),x=Lat, y=fit^4, geom='line')+
  geom_point(data=filter(aus_trop_prop, FG==9),aes(x=Lat, y=trop_met^4), shape=1)+geom_hline(yintercept=1, colour='red')



                
ggplot(jpn_gam_pred[jpn_gam_pred$FG!='all',], aes(x = lat))+  
  geom_ribbon(data=jpn_gam_pred[jpn_gam_pred$FG=='all',]%>%rename(FG2=FG), 
              aes(ymax=(fit+se.fit*1.96)^4, ymin=(fit-se.fit*1.96)^4),fill='yellow', alpha=0.5)+
  geom_line(data=jpn_gam_pred[jpn_gam_pred$FG=='all',]%>%rename(FG2=FG), aes(y = fit^4))+
  geom_ribbon(aes(ymax=(fit+se.fit*1.96)^4, ymin=(fit-se.fit*1.96)^4),fill='red', alpha=0.5)+geom_line(aes(y = fit^4),colour='red')+
  scale_x_continuous(breaks=24:35)+theme_bw()+
  geom_point(data=jpn_trop_prop2, aes(x=lat, y=cor_biom/mean_biom), shape=1, size=0.6)+
  facet_wrap(~FG, scales='free')






# compare delta values with mean FG thermal midpoint data (stuart-smith)
therm_mid<-read_xlsx('C:/coral_fish/sourced_data/stuart_smith_thermal_midpoints/Thermal niche midpoints.xlsx')

dat_ss<-left_join(dat, therm_mid, by=c('Species'='SPECIES_NAME'))
which(is.na(dat_ss$`95th SSTmax`)) # some naming mis-matches/ missing sp ~90
table(dat_ss$ThermalAffinity2, dat_ss$`Temp-Trop (23cutoff)`)

# confidence > 1 to filter out low confidence midpoints
qplot(data=dat_ss[dat_ss$confidence>1,], x=ThermalAffinity, y=`MP(5min-95max)`)+geom_violin()
dat_ss%>%group_by(ThermalAffinity)%>%summarize(n(), mn=mean(`MP(5min-95max)`, na.rm=T),
     md=median(`MP(5min-95max)`, na.rm=T), sd=sd(`MP(5min-95max)`, na.rm=T))

# check for sig dif between classes
library(agricolae)
m1<-lm(`MP(5min-95max)`~ThermalAffinity, data=dat_ss[dat_ss$confidence>1 &
                         dat_ss$ThermalAffinity!='nonarctic',])
t1<-HSD.test(m1, 'ThermalAffinity')

# could us stu smith class or just keep mine. or blend?

# do by region
dat_aus$thermal_mid<-left_join(dat_aus, therm_mid, by=c('Species'='SPECIES_NAME'))$'MP(5min-95max)'
dat_jpn$thermal_mid<-left_join(dat_jpn, therm_mid, by=c('Species'='SPECIES_NAME'))$'MP(5min-95max)'
dat_aus$confidence<-left_join(dat_aus, therm_mid, by=c('Species'='SPECIES_NAME'))$confidence
dat_jpn$confidence<-left_join(dat_jpn, therm_mid, by=c('Species'='SPECIES_NAME'))$confidence

therms<-rbind(
dat_aus%>%filter(AUS_trop==1 & confidence>1)%>%group_by(FG)%>%
  summarise(mean_tm=mean(thermal_mid, na.rm=T),median_tm=median(thermal_mid, na.rm=T), 
            sd_tm=sd(thermal_mid, na.rm=T))%>%mutate(comm='trop', region='Australia'),
dat_aus%>%filter(AUS_temp==1& confidence>1)%>%group_by(FG)%>%
  summarise(mean_tm=mean(thermal_mid, na.rm=T),median_tm=median(thermal_mid, na.rm=T), 
            sd_tm=sd(thermal_mid, na.rm=T))%>%mutate(comm='temp', region='Australia'),
dat_jpn%>%filter(JPN_trop==1& confidence>1)%>%group_by(FG)%>%
  summarise(mean_tm=mean(thermal_mid, na.rm=T),median_tm=median(thermal_mid, na.rm=T), 
            sd_tm=sd(thermal_mid, na.rm=T))%>%mutate(comm='trop', region='Japan'),
dat_jpn%>%filter(JPN_temp==1& confidence>1)%>%group_by(FG)%>%
  summarise(mean_tm=mean(thermal_mid, na.rm=T),median_tm=median(thermal_mid, na.rm=T), 
            sd_tm=sd(thermal_mid, na.rm=T))%>%mutate(comm='temp', region='Japan'))

delta_trop_therm<-left_join(delta_trop, therms, by=c('region', 'comm', 'FG'))

ggplot(data=delta_trop_therm, aes(x=delta_trop, y=mean_tm))+
  geom_pointrange(aes(ymin=mean_tm-sd_tm, ymax=mean_tm+sd_tm))+
  facet_wrap(~region+comm)


ggplot(data=filter(delta_trop_therm, FG%in%c(4,6,1,2)),
       aes(x=delta_trop, y=mean_tm))+geom_vline(xintercept=0, linetype='dotted')+
  geom_pointrange(aes(ymin=mean_tm-sd_tm, ymax=mean_tm+sd_tm))+
  facet_wrap(~region+comm)+theme_bw()+xlab('Delta community prop. tropical')+
  ylab('Thermal midpoint C')

# plots per FG

ggplot(data=dat_jpn%>%filter(JPN_temp==1 &  FG %in%c(4,6,1,2)),
       aes(x=factor(FG), y=thermal_mid))+geom_violin()

### calculate functional distance/overlap between tropical invaders and
## higher latitude residents in transititon zone

func_dudi<-dudi.pco(d = sqrt(eff_both), scannf = FALSE, nf = 4)
# could do func distances for only 4 FGs?

func_pco<-data.frame(func_dudi$li,dat)

# add to regional datasets
dat_aus<-left_join(dat_aus, func_pco[,c(1:5)], by='Species')
dat_jpn<-left_join(dat_jpn, func_pco[,c(1:5)], by='Species')

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
                                          100* func_dudi$eig[2]/sum(func_dudi$eig))))+
  scale_x_continuous(paste('PC1', sprintf('(%0.1f%% explained var.)',
                                          100* func_dudi$eig[1]/sum(func_dudi$eig))))
dat_aus$ThermalAffinity2<-factor(dat_aus$ThermalAffinity2, levels=c('tropical', 'subtropical'))

ppp+geom_point(data=dat_aus%>% filter(AUS_temp==1 & FG %in% c(4,6,1,2)),
                   aes(x=A1, y=A2, colour=ThermalAffinity2))+
  geom_polygon(data=dat_aus %>% filter(AUS_temp==1 & FG %in% c(4,6,1,2))%>%
                 group_by(FG, ThermalAffinity2) %>% slice(chull(A1, A2)),
               aes(x=A1, y=A2, fill=ThermalAffinity2), alpha=0.2)+
  geom_point(data=dat_aus%>% filter(AUS_temp==1 & FG %in% c(4,6,1,2))%>%
               group_by(FG, ThermalAffinity2)%>%summarize(A1=mean(A1), A2=mean(A2)),
             aes(x=A1, y=A2, colour=ThermalAffinity2), size=2)+
  geom_point(data=dat_aus%>% filter(AUS_temp==1 & FG %in% c(4,6,1,2))%>%
               group_by(FG, ThermalAffinity2)%>%summarize(A1=mean(A1), A2=mean(A2)),
             aes(x=A1, y=A2, group=ThermalAffinity2), colour='black', size=2, shape=1)+
  scale_y_continuous(paste('PC2', sprintf('(%0.1f%% explained var.)',
                                          100* func_dudi$eig[2]/sum(func_dudi$eig))))+
  scale_x_continuous(paste('PC1', sprintf('(%0.1f%% explained var.)',
                                          100* func_dudi$eig[1]/sum(func_dudi$eig))))+
  facet_wrap(~FG)+theme_minimal()

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

