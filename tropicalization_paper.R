# paper code 27/06/19

# picks up from stage where optimum number of clusters are known
#Steps
#1 calculate redundancy for each community and group
#2 split functional groups by thermal affinity to calculate 
# redundancy/vulnerability for each FG/TG group due to thermal disturbance
#3 calculate recovery diversity for each FG/TH

library(cluster)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(mice)
library(scales)
library(ggwordcloud)
library(vegan)
library(ade4)
library(readxl)

#################### data clean ##########################
##########################################################

dat<-read.csv('C:/coral_fish/data/Traits/JPN_AUS_RMI_CHK_MLD_TMR_trait_master_opt2_lats.csv', h=T)
# Running on simplist classification of Position trait

# Do we want to include sp. in the analyses? - yes.
# dat[grep('\\.', dat$Species),]

dat<-dat[which(dat$AUS_sp>0 | dat$JPN_sp>0),] # we will focus on Australia and Japan for this prelim

# remove functional duplicates
#dup_trait<-paste(dat$BodySize, dat$DepthRange, dat$PLD, dat$Diet, dat$Aggregation, dat$Position, 
#                 dat$ParentalMode)
#dat<-dat[-which(duplicated(dup_trait)),]
#dat$Species<-as.character(dat$Species)
#dat[which(dat$Species=='Scarus psittacus'),]$Species<-'Scarus psittacus/spinus' # edit for one

# Including those that default to duplicates via NA after gower dist
#dat<-dat[-which(dat$Species=='Caesio sp.'),]
#dat<-dat[-which(dat$Species=='Ostracion immaculatus'),]

row.names(dat)<-dat$Species

## Edit to some trait values from MB 25/10/18
#dat[dat$Species=='Brotula multibarbata',]$DepthRange<-219
dat[dat$Species=='Mobula birostris',]$BodySize<-450
dat[dat$Species=='Amphiprion sandaracinos',]$BodySize<-14

## set ceiling for numeric variables for scaling purposes 04/03/19
dat[dat$PLD>=100 & !is.na(dat$PLD),]$PLD<-100
dat[dat$DepthRange>=200 & !is.na(dat$DepthRange),]$DepthRange<-200

# ORDER necessary categorical variables

dat$Aggregation<-factor(dat$Aggregation, levels=c("solitary", "pairs","groups","schools"), ordered = T)
# Position doensn't follow a logical SINGLE order
#dat$Position<-factor(dat$Position, levels=c("SubBenthic", "Benthic","UpperBenthic",
#                                            "Demersal", "ReefPelagic","Pelagic"), ordered = T)

# Set non-arctic to tropical ThermalAffinuty
dat[dat$ThermalAffinity=='nonarctic',]$ThermalAffinity<-'tropical'
dat$ThermalAffinity<-factor(dat$ThermalAffinity)
# create thermal affinity variable with just tropical/non-tropical
dat$ThermalAffinity2<-as.character(dat$ThermalAffinity)
dat[dat$ThermalAffinity2!='tropical',]$ThermalAffinity2<-'subtropical'

eff_both<-daisy(dat[,c("BodySize","Diet",  "Position", "Aggregation", 'DepthRange')], metric='gower', stand = FALSE)

both_FG<-cutree(hclust(eff_both, method='average'), k=12)

dat$FG<-both_FG
# some finer cuts to increase statistical power
dat$FG_22<-cutree(hclust(eff_both, method='average'), k=22)
dat$FG_36<-cutree(hclust(eff_both, method='average'), k=36)

# Split to regions

dat_aus<-dat[which(dat$AUS_sp>0),]
dat_jpn<-dat[which(dat$JPN_sp>0),]

# plot max lat of tropical species in FGs

ggplot(data=filter(dat, AUS_sp>0 & ThermalAffinity2=='tropical'), 
      aes(x=factor(FG), y= AUS_maxlat))+geom_violin()+
  geom_point(
    data=filter(dat, AUS_sp>0 & ThermalAffinity2=='tropical')%>%
      group_by(FG)%>%summarize(medlat=median(AUS_maxlat)), 
    aes(x=factor(FG), y= medlat), colour='red')

ggplot(data=filter(dat, JPN_sp>0 & ThermalAffinity2=='tropical'), 
       aes(x=factor(FG), y= JPN_maxlat))+geom_violin()+
  geom_point(
    data=filter(dat, JPN_sp>0 & ThermalAffinity2=='tropical')%>%
      group_by(FG)%>%summarize(medlat=median(JPN_maxlat)), 
    aes(x=factor(FG), y= medlat), colour='red')

### Functional Entity creation

dat_mice<-mice(dat[,c(2:9)], m=5, method=c('polyreg',rep('norm.predict', 3), rep('polyreg', 4)))
dat_imp<-complete(dat_mice)
dat_imp<-cbind(Species=dat[,1], dat_imp, dat[,10:length(dat)])

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

FEdat<-dat_imp %>% group_by(FE) %>% summarise(FG=unique(FG),num=n())

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
out<-lapply(FEword.agg.cl, function(x){
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

## first Redundancy plot with thermal affinity split

red_dat<-rbind(dat_aus %>% group_by(FG) %>% summarise(num=n()) %>%
  arrange(num) %>% mutate(val=12:1, dat='Australia'),
  dat_jpn %>% group_by(FG) %>% summarise(num=n()) %>%
    arrange(num) %>% mutate(val=10:1, dat='Japan'))

red_dat_therm<-rbind(dat_aus %>% group_by(ThermalAffinity2, FG) %>%
                       summarise(num=n()) %>% mutate(dat='Australia'),
                     dat_jpn %>% group_by(ThermalAffinity2, FG) %>%
                       summarise(num=n())%>% mutate( dat='Japan'))

red_dat_therm$val<-left_join(red_dat_therm, red_dat, by=c('dat', 'FG'))$val

ggplot()+
  geom_bar(data=red_dat_therm, aes(x=val, y=num, fill=ThermalAffinity2),
           stat='identity',position=position_dodge(preserve = 'single'))+
  geom_hline(data=red_dat_therm%>%group_by(dat, ThermalAffinity2)%>%summarise_all(mean),
             aes(yintercept=num, colour=ThermalAffinity2), linetype='dashed')+
  facet_grid(dat~.)+
  scale_x_continuous(breaks=1:12)+
  xlab('Rank of functional group')+ylab('# of species per FG')+
  theme_bw()

## Detailed Redundancy plot with thermal affinity and community splits - JAPAN

jpn_dat_therm<-rbind(dat_jpn %>% filter(JPN_trop==1) %>% 
                       group_by(ThermalAffinity2, FG) %>%
                       summarise(num=n()) %>% mutate(dat='Tropical'),
                     dat_jpn %>% filter(JPN_temp==1) %>%
                       group_by(ThermalAffinity2, FG) %>%
                       summarise(num=n())%>% mutate( dat='Transition'))

jpn_dat_therm$val<-left_join(jpn_dat_therm, filter(red_dat,dat=='Japan'), by='FG')$val

ggplot()+
  geom_bar(data=jpn_dat_therm, aes(x=val, y=num,
  fill=interaction(ThermalAffinity2, dat)),stat='identity',
  position=position_dodge2(preserve = 'single'))+
  geom_hline(data=jpn_dat_therm%>%group_by(dat)%>%summarise_all(mean),
             aes(yintercept=num),colour=c('#abd9e9','#d7191c'))+
  geom_hline(data=jpn_dat_therm%>%group_by(dat)%>%summarise_all(mean),
             aes(yintercept=num),colour=c('#2c7bb6', '#fdae61'), linetype='dotted', size=1)+
  scale_x_continuous(breaks=1:10, labels=FGnames[c(4,6,1,2,8,3,5,10,9, 7),]$FEcomp)+
  ylab('# of species per FG')+xlab('Functional niche')+
  theme_bw()+ guides(fill=guide_legend(title="ThermAffin.Comm"))+
  scale_fill_manual(values = c('#2c7bb6','#abd9e9', '#fdae61', '#d7191c'))+
  theme(legend.justification = c(1, 1), legend.position = c(1, 1), 
        legend.background = element_rect(colour = "black"), 
        axis.text.x = element_text(angle = 60, hjust = 1))

## Redundancy plot with thermal affinity and community splits - AUSTRALIA

aus_dat_therm<-rbind(dat_aus %>% filter(AUS_trop==1) %>% 
                       group_by(ThermalAffinity2, FG) %>%
                       summarise(num=n()) %>% mutate(dat='Tropical'),
                     dat_aus %>% filter(AUS_temp==1) %>%
                       group_by(ThermalAffinity2, FG) %>%
                       summarise(num=n())%>% mutate( dat='Transition'))

aus_dat_therm$val<-left_join(aus_dat_therm, filter(red_dat,dat=='Australia'), by='FG')$val

ggplot()+
  geom_bar(data=aus_dat_therm, aes(x=val, y=num,
                                   fill=interaction(ThermalAffinity2, dat)),stat='identity',
           position=position_dodge2(preserve = 'single'))+
  geom_hline(data=aus_dat_therm%>%group_by(dat)%>%summarise_all(mean),
             aes(yintercept=num),colour=c('#abd9e9','#d7191c'))+
  geom_hline(data=aus_dat_therm%>%group_by(dat)%>%summarise_all(mean),
             aes(yintercept=num),colour=c('#2c7bb6', '#fdae61'), linetype='dotted', size=1)+
  scale_x_continuous(breaks=1:12, labels=
                       FGnames[c(4,6,1,2,5,8,7,3,9,10,12,11),]$FEcomp)+
  ylab('# of species per FG')+xlab('Rank of functional group')+
  theme_bw()+ guides(fill=guide_legend(title="ThermAffin.Comm"))+
  scale_fill_manual(values = c('#2c7bb6','#abd9e9', '#fdae61', '#d7191c'))+
  theme(legend.justification = c(1, 1), legend.position = c(1, 1), 
        legend.background = element_rect(colour = "black"), 
        axis.text.x = element_text(angle = 60, hjust = 1))+
  scale_y_continuous(breaks=c(0,25,50,75,100,125))


##### Proportion of tropical species plot


dat_aus %>% group_by(FG, ThermalAffinity) %>% summarise(n()) %>% as.data.frame()
dat_jpn %>% group_by(FG, ThermalAffinity) %>% summarise(n()) %>% as.data.frame()

# probability that each FG has more non-tropical species
# than due to random chance

#Australia

out_aus<-NULL
for(i in 1:1000){
out_aus<-rbind(out_aus, dat_aus%>%mutate(FG2=sample(FG, replace=F))%>% 
          group_by(FG2)%>%summarize(prop_trop=length(which(ThermalAffinity=='tropical'))/n())%>%
          mutate(run=i))
}

#full
aus_full<-dat_aus%>%group_by(FG)%>%summarize(prop_trop=length(which(ThermalAffinity=='tropical'))/n())

pa<-ggplot()+geom_violin(data=out_aus, aes(x=factor(FG2), y=prop_trop))+
  geom_hline(yintercept = length(which(dat_aus$ThermalAffinity=='tropical'))/nrow(dat_aus),
             linetype='dashed', colour='dark green')+
  geom_point(data=aus_full, aes(x=factor(FG), y=prop_trop), colour='red')+
  geom_label(data=dat_aus%>%group_by(FG)%>%summarize(n=n()),
  aes(x=factor(FG), y=0.1, label=n), colour='blue')+
  xlab('Australia FGs')+ylab('Proportion of tropical sp.')

## Japan

out_jpn<-NULL
for(i in 1:9999){
  out_jpn<-rbind(out_jpn, dat_jpn%>%mutate(FG2=sample(FG, replace=F))%>% 
               group_by(FG2)%>%summarize(prop_trop=length(which(ThermalAffinity=='tropical'))/n())%>%
               mutate(run=i))
}

#full
jpn_full<-dat_jpn%>%group_by(FG)%>%summarize(prop_trop=length(which(ThermalAffinity=='tropical'))/n())

# hack to add extra group
out_jpn<-rbind(out_jpn, data.frame(FG2=c(11,12), prop_trop=c(NA, NA), run=c(1,1)))

pj<-ggplot()+geom_violin(data=out_jpn, aes(x=factor(FG2), y=prop_trop))+
  geom_hline(yintercept = length(which(dat_jpn$ThermalAffinity=='tropical'))/nrow(dat_jpn),
             linetype='dashed', colour='dark green')+
  geom_point(data=jpn_full, aes(x=factor(FG), y=prop_trop), colour='red')+
  geom_label(data=dat_jpn%>%group_by(FG)%>%summarize(n=n()),
  aes(x=as.factor(FG), y=0.1, label=n), colour='blue')+
  xlab('Japan FGs')+ylab('Proportion of tropical sp.')

grid.arrange(pa, pj, nrow=2)

## do any FEs have species with different thermal affinity?

imp_aus %>% group_by(FE) %>% summarise(n(), n_therm=length(unique(ThermalAffinity))) %>% as.data.frame()
dat_jpn %>% group_by(FG, ThermalAffinity) %>% summarise(n()) %>% as.data.frame()

# Compare ThermalAffinity class vs latitudinal sampling class JAPAN

# run with latitude data
out_jpn_trop<-NULL
for(i in 1:1000){
  out_jpn_trop<-rbind(out_jpn_trop, dat_jpn%>% filter(JPN_trop==1)%>%
                        mutate(FG2=sample(FG, replace=F))%>% group_by(FG2)%>%
                        summarize(prop_trop=length(which(ThermalAffinity2=='tropical'))/n())%>%
                   mutate(run=i))}
out_jpn_temp<-NULL
for(i in 1:1000){
  out_jpn_temp<-rbind(out_jpn_temp, dat_jpn%>% filter(JPN_temp==1)%>%
                        mutate(FG2=sample(FG, replace=F))%>% group_by(FG2)%>%
                        summarize(prop_trop=length(which(ThermalAffinity2=='tropical'))/n())%>%
                        mutate(run=i))}

#full latitude data
jpn_full_trop<-dat_jpn%>%filter(JPN_trop==1)%>%group_by(FG)%>%
  summarize(prop_trop=length(which(ThermalAffinity2=='tropical'))/n(), n=n())
jpn_full_temp<-dat_jpn%>%filter(JPN_temp==1)%>%group_by(FG)%>%
  summarize(prop_trop=length(which(ThermalAffinity2=='tropical'))/n(), n=n())

# plot ThermalAffin vs Latitudeinal classification

ggplot()+geom_violin(data=out_jpn_trop, aes(x=factor(FG2), y=prop_trop), fill='red', alpha=0.3, colour=NA)+
  geom_violin(data=out_jpn_temp, aes(x=factor(FG2), y=prop_trop), fill='blue', alpha=0.3, colour=NA)+
  geom_hline(yintercept = length(which(dat_jpn[dat_jpn$JPN_trop==1,]$ThermalAffinity2=='tropical'))/nrow(dat_jpn[dat_jpn$JPN_trop==1,]),
             linetype='dashed', colour='dark red')+
  geom_point(data=jpn_full_trop, aes(x=factor(FG), y=prop_trop), colour='red', size=2.2)+
 
  geom_hline(yintercept = length(which(dat_jpn[dat_jpn$JPN_temp==1,]$ThermalAffinity2=='tropical'))/nrow(dat_jpn[dat_jpn$JPN_temp==1,]),
             linetype='dashed', colour='dark blue')+
  geom_point(data=jpn_full_temp, aes(x=factor(FG), y=prop_trop), colour='blue', size=2)+
  geom_label(data=dat_jpn%>%filter(JPN_trop==1)%>%group_by(FG)%>%summarize(n=n()),
             aes(x=as.factor(FG), y=0.15, label=n), colour='red')+
  geom_label(data=dat_jpn%>%filter(JPN_temp==1)%>%group_by(FG)%>%summarize(n=n()),
             aes(x=as.factor(FG), y=0.05, label=n), colour='blue')+
  xlab('Functional Niche')+ylab('Proportion of tropical sp.')+
  scale_x_discrete(limits=c(4,6, 1, 2, 8, 3, 5, 10, 9, 7),
                  labels=c('1'=FGnames[4,]$FEcomp,'2'=FGnames[6,]$FEcomp,
                           '3'=FGnames[1,]$FEcomp, '4'=FGnames[2,]$FEcomp,
                           '5'=FGnames[8,]$FEcomp, '6'=FGnames[3,]$FEcomp,
                           '7'=FGnames[5,]$FEcomp, '8'=FGnames[10,]$FEcomp,
                           '9'=FGnames[9,]$FEcomp, '10'=FGnames[7,]$FEcomp))+
  theme_bw()+theme(axis.text.x = element_text(angle = 55, hjust = 1))

# Compare ThermalAffinity class vs latitudinal sampling class AUSTRALIA

# run with latitude data
out_aus_trop<-NULL
for(i in 1:1000){
  out_aus_trop<-rbind(out_aus_trop, dat_aus%>% filter(AUS_trop==1)%>%
                        mutate(FG2=sample(FG, replace=F))%>% group_by(FG2)%>%
                        summarize(prop_trop=length(which(ThermalAffinity2=='tropical'))/n())%>%
                        mutate(run=i))}
out_aus_temp<-NULL
for(i in 1:1000){
  out_aus_temp<-rbind(out_aus_temp, dat_aus%>% filter(AUS_temp==1)%>%
                        mutate(FG2=sample(FG, replace=F))%>% group_by(FG2)%>%
                        summarize(prop_trop=length(which(ThermalAffinity2=='tropical'))/n())%>%
                        mutate(run=i))}

#full latitude data
aus_full_trop<-dat_aus%>%filter(AUS_trop==1)%>%group_by(FG)%>%
  summarize(prop_trop=length(which(ThermalAffinity2=='tropical'))/n(), n=n())
aus_full_temp<-dat_aus%>%filter(AUS_temp==1)%>%group_by(FG)%>%
  summarize(prop_trop=length(which(ThermalAffinity2=='tropical'))/n(), n=n())

# plot ThermalAffin vs Latitudeinal classification

ggplot()+geom_violin(data=out_aus_trop, aes(x=factor(FG2), y=prop_trop), fill='red', alpha=0.3, colour=NA)+
  geom_violin(data=out_aus_temp, aes(x=factor(FG2), y=prop_trop), fill='blue', alpha=0.3, colour=NA)+
  geom_hline(yintercept = length(which(dat_aus[dat_aus$AUS_trop==1,]$ThermalAffinity2=='tropical'))/nrow(dat_aus[dat_aus$AUS_trop==1,]),
             linetype='dashed', colour='dark red')+
  geom_point(data=aus_full_trop, aes(x=factor(FG), y=prop_trop), colour='red', size=2.2)+
  
  geom_hline(yintercept = length(which(dat_aus[dat_aus$AUS_temp==1,]$ThermalAffinity2=='tropical'))/nrow(dat_aus[dat_aus$AUS_temp==1,]),
             linetype='dashed', colour='dark blue')+
  geom_point(data=aus_full_temp, aes(x=factor(FG), y=prop_trop), colour='blue', size=2)+
  geom_label(data=dat_aus%>%filter(AUS_trop==1)%>%group_by(FG)%>%summarize(n=n()),
             aes(x=as.factor(FG), y=0.15, label=n), colour='red')+
  geom_label(data=dat_aus%>%filter(AUS_temp==1)%>%group_by(FG)%>%summarize(n=n()),
             aes(x=as.factor(FG), y=0.05, label=n), colour='blue')+
  xlab('Functional Niche')+ylab('Proportion of tropical sp.')+
  scale_x_discrete(limits=c(4,6,1,2,5,8,7,3,9,10,12,11),
                   labels=c('1'=FGnames[4,]$FEcomp,'2'=FGnames[6,]$FEcomp,
                            '3'=FGnames[1,]$FEcomp, '4'=FGnames[2,]$FEcomp,
                            '5'=FGnames[5,]$FEcomp, '6'=FGnames[8,]$FEcomp,
                            '7'=FGnames[7,]$FEcomp, '8'=FGnames[3,]$FEcomp,
                            '9'=FGnames[9,]$FEcomp, '10'=FGnames[10,]$FEcomp, 
                            '11'=FGnames[12,]$FEcomp, '12'=FGnames[11,]$FEcomp))+
  theme_bw()+theme(axis.text.x = element_text(angle = 55, hjust = 1))

# Compare ThermalAffinity within TRANSITION ZONE class JAPAN

# run with latitude data
out_jpn_temp1<-NULL
for(i in 1:1000){
  out_jpn_temp1<-rbind(out_jpn_temp1, dat_jpn%>% filter(JPN_temp==1 &
                         FG %in% c(4,6,1,2))%>%
                        mutate(FG2=sample(FG, replace=F))%>% group_by(FG2)%>%
                        summarize(prop=length(which(ThermalAffinity=='subtropical'))/n())%>%
                        mutate(run=i, therm.niche='subtropical'))}
out_jpn_temp2<-NULL
for(i in 1:1000){
  out_jpn_temp2<-rbind(out_jpn_temp2, dat_jpn%>% filter(JPN_temp==1 &
                        FG %in% c(4,6,1,2))%>%
                        mutate(FG2=sample(FG, replace=F))%>% group_by(FG2)%>%
                        summarize(prop=length(which(ThermalAffinity=='temperate'))/n())%>%
                        mutate(run=i, therm.niche='temperate'))}

out_jpn_temp<-rbind(out_jpn_temp1, out_jpn_temp2)

#full latitude data
jpn_full_temp_trans<-rbind(dat_jpn%>%filter(JPN_temp==1&
   FG %in% c(4,6,1,2))%>%group_by(FG)%>%
  summarize(prop=length(which(ThermalAffinity=='subtropical'))/n())%>%
    mutate(therm.niche='subtropical'),
  dat_jpn%>%filter(JPN_temp==1&
  FG %in% c(4,6,1,2))%>%group_by(FG)%>%summarize(
  prop=length(which(ThermalAffinity=='temperate'))/n())%>%
  mutate(therm.niche='temperate'))

# plot ThermalAffin vs Latitudeinal classification

out_jpn_temp$FG2<-factor(out_jpn_temp$FG2, levels=c(4,6,1,2))
jpn_full_temp_trans$FG<-factor(jpn_full_temp_trans$FG, levels=c(4,6,1,2))

ggplot()+
  geom_violin(data=out_jpn_temp, aes(x=FG2, y=prop, fill=therm.niche),
              alpha=0.3, colour=NA, position=position_dodge())+
    scale_fill_manual(values = c('cyan','blue'))+
  geom_hline(yintercept = length(which(dat_jpn[dat_jpn$JPN_temp==1,]$ThermalAffinity=='temperate'))/nrow(dat_jpn[dat_jpn$JPN_temp==1,]),
             linetype='dashed', colour='dark blue')+
  geom_hline(yintercept = length(which(dat_jpn[dat_jpn$JPN_temp==1,]$ThermalAffinity=='subtropical'))/nrow(dat_jpn[dat_jpn$JPN_temp==1,]),
             linetype='dashed', colour='#238b45')+
  geom_point(data=jpn_full_temp_trans, aes(x=FG, y=prop, colour=therm.niche), size=2,
             position=position_dodge(width=0.9))+
    scale_colour_manual(values = c('#238b45','blue'))+

  xlab('Functional Niche')+ylab('Proportion of subtropical/temperate sp.')+
  scale_x_discrete(labels=c(FGnames[4,]$FEcomp,FGnames[6,]$FEcomp,
                            FGnames[1,]$FEcomp,FGnames[2,]$FEcomp))+
  theme_bw()+theme(axis.text.x = element_text(angle = 55, hjust = 1),
                   legend.justification = c(0, 1), legend.position = c(0, 1), 
                   legend.background = element_rect(colour = "black"))+
  labs(fill='Thermal Niche', colour='Thermal Niche')+ylim(c(0, 0.4))

# Compare ThermalAffinity within TRANSITION ZONE class Australia

# run with latitude data
out_aus_temp1<-NULL
for(i in 1:1000){
  out_aus_temp1<-rbind(out_aus_temp1, dat_aus%>% filter(AUS_temp==1 &
                                                          FG %in% c(4,6,1,2))%>%
                         mutate(FG2=sample(FG, replace=F))%>% group_by(FG2)%>%
                         summarize(prop=length(which(ThermalAffinity=='subtropical'))/n())%>%
                         mutate(run=i, therm.niche='subtropical'))}
out_aus_temp2<-NULL
for(i in 1:1000){
  out_aus_temp2<-rbind(out_aus_temp2, dat_aus%>% filter(AUS_temp==1 &
                                                          FG %in% c(4,6,1,2))%>%
                         mutate(FG2=sample(FG, replace=F))%>% group_by(FG2)%>%
                         summarize(prop=length(which(ThermalAffinity=='temperate'))/n())%>%
                         mutate(run=i, therm.niche='temperate'))}

out_aus_temp<-rbind(out_aus_temp1, out_aus_temp2)

#full latitude data
aus_full_temp_trans<-rbind(dat_aus%>%filter(AUS_temp==1&
                     FG %in% c(4,6,1,2))%>%group_by(FG)%>%
                       summarize(prop=length(which(ThermalAffinity=='subtropical'))/n())%>%
                       mutate(therm.niche='subtropical'),
                     dat_aus%>%filter(AUS_temp==1&
                     FG %in% c(4,6,1,2))%>%group_by(FG)%>%summarize(
                     prop=length(which(ThermalAffinity=='temperate'))/n())%>%
                       mutate(therm.niche='temperate'))

# plot ThermalAffin vs Latitudeinal classification

out_aus_temp$FG2<-factor(out_aus_temp$FG2, levels=c(4,6,1,2))
aus_full_temp_trans$FG<-factor(aus_full_temp_trans$FG, levels=c(4,6,1,2))

ggplot()+
  geom_violin(data=out_aus_temp, aes(x=FG2, y=prop, fill=therm.niche),
              alpha=0.3, colour=NA, position=position_dodge())+
  scale_fill_manual(values = c('cyan','blue'))+
  geom_hline(yintercept = length(which(dat_aus[dat_aus$AUS_temp==1,]$ThermalAffinity=='temperate'))/nrow(dat_aus[dat_aus$AUS_temp==1,]),
             linetype='dashed', colour='dark blue')+
  geom_hline(yintercept = length(which(dat_aus[dat_aus$AUS_temp==1,]$ThermalAffinity=='subtropical'))/nrow(dat_aus[dat_aus$AUS_temp==1,]),
             linetype='dashed', colour='#238b45')+
  geom_point(data=aus_full_temp_trans, aes(x=FG, y=prop, colour=therm.niche), size=2,
             position=position_dodge(width=0.9))+
  scale_colour_manual(values = c('#238b45','blue'))+
  
  xlab('Functional Niche')+ylab('Proportion of subtropical/temperate sp.')+
  scale_x_discrete(labels=c(FGnames[4,]$FEcomp,FGnames[6,]$FEcomp,
                            FGnames[1,]$FEcomp,FGnames[2,]$FEcomp))+
  theme_bw()+theme(axis.text.x = element_text(angle = 55, hjust = 1),
                   legend.justification = c(0, 1), legend.position = c(0, 1), 
                   legend.background = element_rect(colour = "black"))+
  labs(fill='Thermal Niche', colour='Thermal Niche')+ylim(c(0, 0.4))



## Get delta values between community prop trop and FG prop trop

aus_full_temp<-dat_aus%>%filter(AUS_temp==1)%>%group_by(FG)%>%
  summarize(prop_trop=length(which(ThermalAffinity2=='tropical'))/n(), n=n())
aus_full_trop<-dat_aus%>%filter(AUS_trop==1)%>%group_by(FG)%>%
  summarize(prop_trop=length(which(ThermalAffinity2=='tropical'))/n(), n=n())
jpn_full_temp<-dat_jpn%>%filter(JPN_temp==1)%>%group_by(FG)%>%
  summarize(prop_trop=length(which(ThermalAffinity2=='tropical'))/n(), n=n())
jpn_full_trop<-dat_jpn%>%filter(JPN_trop==1)%>%group_by(FG)%>%
  summarize(prop_trop=length(which(ThermalAffinity2=='tropical'))/n(), n=n())


aus_full_temp$comm_trop<-length(which(dat_aus[dat_aus$AUS_temp==1,]$ThermalAffinity2=='tropical'))/nrow(dat_aus[dat_aus$AUS_temp==1,])
aus_full_temp$delta_trop<-aus_full_temp$prop_trop-aus_full_temp$comm_trop
aus_full_temp$comm='temp'
aus_full_temp$region='Australia'
aus_full_trop$comm_trop<-length(which(dat_aus[dat_aus$AUS_trop==1,]$ThermalAffinity2=='tropical'))/nrow(dat_aus[dat_aus$AUS_trop==1,])
aus_full_trop$delta_trop<-aus_full_trop$prop_trop-aus_full_trop$comm_trop
aus_full_trop$comm='trop'
aus_full_trop$region='Australia'
jpn_full_temp$comm_trop<-length(which(dat_jpn[dat_jpn$JPN_temp==1,]$ThermalAffinity2=='tropical'))/nrow(dat_jpn[dat_jpn$JPN_temp==1,])
jpn_full_temp$delta_trop<-jpn_full_temp$prop_trop-jpn_full_temp$comm_trop
jpn_full_temp$comm='temp'
jpn_full_temp$region='Japan'
jpn_full_trop$comm_trop<-length(which(dat_jpn[dat_jpn$JPN_trop==1,]$ThermalAffinity2=='tropical'))/nrow(dat_jpn[dat_jpn$JPN_trop==1,])
jpn_full_trop$delta_trop<-jpn_full_trop$prop_trop-jpn_full_trop$comm_trop
jpn_full_trop$comm='trop'
jpn_full_trop$region='Japan'

delta_trop<-rbind(aus_full_trop, aus_full_temp, jpn_full_trop, jpn_full_temp)

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

