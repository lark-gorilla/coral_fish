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

#################### data clean ##########################
##########################################################

dat<-read.csv('C:/coral_fish/data/Traits/JPN_AUS_RMI_CHK_MLD_TMR_trait_master_opt2.csv', h=T)
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

# create thermal affinity variable with just tropical/non-tropical
dat$ThermalAffinity2<-as.character(dat$ThermalAffinity)
dat[dat$ThermalAffinity2!='tropical',]$ThermalAffinity2<-'subtropical'

eff_both<-daisy(dat[,c("BodySize","Diet",  "Position", "Aggregation", 'DepthRange')], metric='gower', stand = FALSE)

both_FG<-cutree(hclust(eff_both, method='average'), k=12)

dat$FG<-both_FG

# Add regional tropical/sub-tropical comm data and split to regions
jpn_trop<-read.csv('C:/coral_fish/data/Japan/JPN_species_tropical_class.csv')
aus_trop<-read.csv('C:/coral_fish/data/Australia/AUS_species_tropical_class.csv')


dat_aus<-dat[which(dat$AUS_sp>0),]
dat_jpn<-dat[which(dat$JPN_sp>0),]

dat_aus<-left_join(dat_aus, aus_trop, by='Species')
dat_jpn<-left_join(dat_jpn, jpn_trop, by='Species')

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
## Redundancy plot with thermal affinity split

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

## Redundancy plot with thermal affinity and community splits - JAPAN

jpn_dat_therm<-rbind(dat_jpn %>% filter(JPN_trop==1) %>% 
                       group_by(ThermalAffinity2, FG) %>%
                       summarise(num=n()) %>% mutate(dat='Trop'),
                     dat_jpn %>% filter(JPN_temp==1) %>%
                       group_by(ThermalAffinity2, FG) %>%
                       summarise(num=n())%>% mutate( dat='Temp'))


jpn_dat_therm$val<-left_join(jpn_dat_therm, filter(red_dat,dat=='Japan'), by='FG')$val


ggplot()+
  geom_bar(data=jpn_dat_therm, aes(x=val, y=num,
  fill=interaction(ThermalAffinity2, dat)),stat='identity',
  position=position_dodge2(preserve = 'single'))+
  geom_hline(data=jpn_dat_therm%>%group_by(dat)%>%summarise_all(mean),
             aes(yintercept=num),colour=c('#abd9e9','#d7191c'))+
  geom_hline(data=jpn_dat_therm%>%group_by(dat)%>%summarise_all(mean),
             aes(yintercept=num),colour=c('#2c7bb6', '#fdae61'), linetype='dotted', size=1)+
  scale_x_continuous(breaks=1:10, labels=c(4,6,1,2,8,3,5,10,9, 7))+
  ylab('# of species per FG')+xlab('Rank of functional group')+
  theme_bw()+ guides(fill=guide_legend(title="ThermAffin.Comm"))+
  scale_fill_manual(values = c('#abd9e9','#2c7bb6', '#fdae61', '#d7191c'))+
  theme(legend.justification = c(1, 1), legend.position = c(1, 1), 
        legend.background = element_rect(colour = "black"))

## Redundancy plot with thermal affinity and community splits - AUSTRALIA

aus_dat_therm<-rbind(dat_aus %>% filter(AUS_trop==1) %>% 
                       group_by(ThermalAffinity2, FG) %>%
                       summarise(num=n()) %>% mutate(dat='Trop'),
                     dat_aus %>% filter(AUS_temp==1) %>%
                       group_by(ThermalAffinity2, FG) %>%
                       summarise(num=n())%>% mutate( dat='Temp'))


aus_dat_therm$val<-left_join(aus_dat_therm, filter(red_dat,dat=='Australia'), by='FG')$val


ggplot()+
  geom_bar(data=aus_dat_therm, aes(x=val, y=num,
                                   fill=interaction(ThermalAffinity2, dat)),stat='identity',
           position=position_dodge2(preserve = 'single'))+
  geom_hline(data=aus_dat_therm%>%group_by(dat)%>%summarise_all(mean),
             aes(yintercept=num),colour=c('#abd9e9','#d7191c'))+
  geom_hline(data=aus_dat_therm%>%group_by(dat)%>%summarise_all(mean),
             aes(yintercept=num),colour=c('#2c7bb6', '#fdae61'), linetype='dotted', size=1)+
  scale_x_continuous(breaks=1:12, labels=unique(aus_dat_therm[order(aus_dat_therm$val),]$FG))+
  ylab('# of species per FG')+xlab('Rank of functional group')+
  theme_bw()+ guides(fill=guide_legend(title="ThermAffin.Comm"))+
  scale_fill_manual(values = c('#abd9e9','#2c7bb6', '#fdae61', '#d7191c'))+
  theme(legend.justification = c(1, 1), legend.position = c(1, 1), 
        legend.background = element_rect(colour = "black"))


  
## Redundancy/complimentarity plot

imp_aus<-dat_imp[which(dat_imp$AUS_sp>0),]
imp_jpn<-dat_imp[which(dat_imp$JPN_sp>0),]

#checks
length(unique(imp_aus$FE))
length(unique(imp_jpn$FE))
# All FEs are uniquely nested within FGs
max(aggregate(FG~FE, imp_aus, function(x){table(unique(x))})$FG)
max(aggregate(FG~FE, imp_jpn, function(x){table(unique(x))})$FG)

com_dat<-rbind(imp_aus %>% group_by(FG) %>% summarise(n_FE=length(unique(FE)),num=n()) %>%
                 arrange(num) %>% mutate(val=12:1, dat='Australia'),
                 imp_jpn %>% group_by(FG) %>% summarise(n_FE=length(unique(FE)),num=n()) %>%
                 arrange(num) %>% mutate(val=10:1, dat='Japan'))

mn_prop<-com_dat%>%mutate(prop_fefg=n_FE/num)%>%
             group_by(dat)%>%summarise_all(mean)

com_dat$prop<-mn_prop$prop_fefg[1]
com_dat[com_dat$dat=='Japan',]$prop<-mn_prop$prop_fefg[2]

p1<-ggplot()+
  geom_bar(data=com_dat, aes(x=val, y=num, fill=factor(FG)), stat='identity')+
  geom_bar(data=com_dat, aes(x=val, y=n_FE),colour='black', fill='black', stat='identity', alpha=0.6)+
  geom_hline(data=com_dat%>%group_by(dat)%>%summarise_all(mean),
             aes(yintercept=num), linetype='dashed')+
  geom_point(data=com_dat, aes(x=val, y=num*prop),colour='red')+
  facet_wrap(~dat)+
  scale_x_continuous(breaks=1:12)+
  xlab('Rank of functional group')+ylab('Richness per functional group')+
  theme_bw()+theme(legend.position = 'none')

# sp per FE mean+sd plot
com_dat2<-rbind(imp_aus %>% group_by(FE) %>% summarise(FG=unique(FG),num=n()) %>%
                  arrange(num, FG) %>% mutate(val=260:1, dat='Australia'),
                imp_jpn %>% group_by(FE) %>% summarise(FG=unique(FG),num=n()) %>%
                  arrange(num, FG) %>% mutate(val=204:1, dat='Japan'))

fg_mn<-com_dat2 %>% group_by(dat, FG) %>% summarise(nm_mn=mean(num), nm_sd=sd(num),
                                                    val_mn=mean(val), val_sd=sd(val))

fg_mn$val<-com_dat[order(com_dat$dat, com_dat$FG),]$val

# for mean level different vals depending on mean - check
fg_mn%>%group_by(dat)%>%summarise(mean(nm_mn))
com_dat2%>%group_by(dat)%>%summarise(mean(num))

# running hline on mean og FG means
p2<-ggplot(data=fg_mn, aes(x=val, y=nm_mn))+
  geom_errorbar(aes(ymin=nm_mn-nm_sd, ymax=nm_mn+nm_sd))+
  geom_bar(aes(fill=factor(FG)),colour='black',stat='identity')+
  geom_hline(data=fg_mn%>%group_by(dat)%>%summarise(mean(nm_mn)),
             aes(yintercept=`mean(nm_mn)`), linetype='dashed')+
  facet_wrap(~dat)+
  scale_x_continuous(breaks=1:12)+
  scale_y_continuous(limits=c(0, 6.5), breaks=(0:6), oob=rescale_none)+
  xlab('Rank of functional group')+ylab('Species per functional entity')+
  theme_bw()+theme(legend.position = 'none')

com_dat3<-com_dat2
com_dat3<-left_join(com_dat3, com_dat[,c(1, 4, 5)], by=c('dat', 'FG'))


p2.5<-ggplot(data=com_dat3, aes(x=val.y, y=num))+
  
  geom_violin(aes(fill=factor(FG)), position='dodge')+
  facet_grid(dat~.)+scale_x_continuous(breaks=1:12)
# ig nore colouring, it has reordered cos has missed FGs 11 & 12
               
grid.arrange(p1, p2)

## Mouillot redundancy plot ??

com_dat2<-rbind(imp_aus %>% group_by(FE) %>% summarise(FG=unique(FG),num=n()) %>%
                 arrange(num, FG) %>% mutate(val=260:1, dat='Australia'),
               imp_jpn %>% group_by(FE) %>% summarise(FG=unique(FG),num=n()) %>%
                 arrange(num, FG) %>% mutate(val=204:1, dat='Japan'))

fg_mn<-com_dat2 %>% group_by(dat, FG) %>% summarise(nm_mn=mean(num), nm_sd=sd(num),
                                                   val_mn=mean(val), val_sd=sd(val))

ggplot()+
  geom_bar(data=com_dat, aes(x=val, y=num, fill=factor(FG)), stat='identity')+
  geom_boxplot(data=com_dat, aes(x=FG+248, y=num, colour=factor(FG)))+

  geom_hline(data=com_dat%>%group_by(dat)%>%summarise_all(mean),
             aes(yintercept=num), linetype='dashed')+
  facet_grid(dat~.)+
  scale_x_continuous(breaks=1:260)+
  xlab('Rank of functional group')+ylab('# of species')+
  theme_bw()

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


##### Proportion of tropical species plot


dat_aus %>% group_by(FG, ThermalAffinity) %>% summarise(n()) %>% as.data.frame()
dat_jpn %>% group_by(FG, ThermalAffinity) %>% summarise(n()) %>% as.data.frame()

# probability that each FG has more non-tropical species
# than due to random chance

#Australia

out_aus<-NULL
for(i in 1:9999){
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

# Compare ThermalAffinity class vs latitudinal sampling class

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

#full latuitude data
jpn_full_trop<-dat_jpn%>%filter(JPN_trop==1)%>%group_by(FG)%>%
  summarize(prop_trop=length(which(ThermalAffinity2=='tropical'))/n())
jpn_full_temp<-dat_jpn%>%filter(JPN_temp==1)%>%group_by(FG)%>%
  summarize(prop_trop=length(which(ThermalAffinity2=='tropical'))/n())


# plot ThermalAffin vs Latitudeinal classification

ggplot()+geom_violin(data=out_jpn_trop, aes(x=factor(FG2), y=prop_trop), fill='red', alpha=0.3)+
  geom_violin(data=out_jpn_temp, aes(x=factor(FG2), y=prop_trop), fill='blue', alpha=0.3)+
  geom_hline(yintercept = length(which(dat_jpn[dat_jpn$JPN_trop==1,]$ThermalAffinity2=='tropical'))/nrow(dat_jpn[dat_jpn$JPN_trop==1,]),
             linetype='dashed', colour='dark red')+
  geom_point(data=jpn_full_trop, aes(x=factor(FG), y=prop_trop), colour='red', size=2.2)+
 
  geom_hline(yintercept = length(which(dat_jpn[dat_jpn$JPN_temp==1,]$ThermalAffinity2=='tropical'))/nrow(dat_jpn[dat_jpn$JPN_temp==1,]),
             linetype='dashed', colour='dark blue')+
  geom_point(data=jpn_full_temp, aes(x=factor(FG), y=prop_trop), colour='blue', size=2)+
  geom_label(data=dat_jpn%>%group_by(FG)%>%summarize(n=n()),
             aes(x=as.factor(FG), y=0.1, label=n), colour='blue')+
  xlab('Japan FGs')+ylab('Proportion of tropical sp.')+
  scale_x_discrete(limits=c(4,6, 1, 2, 8, 3, 5, 10, 9, 7),
                  labels=c('1'='4','2'='6', '3'='1', '4'='2', '5'='8', '6'='3', '7'='5', '8'='10', '9'='9', '10'='7'))

#check for subtropical thermalaffin sp that have become generalists?
table(dat_jpn[dat_jpn$ThermalAffinity!='tropical',]$class)



# projection of FGs into 'colonisation space'

d_colon<-daisy(dat[c('BodySize', 'PLD', 'ParentalMode')], metric='gower', stand = FALSE)

ppp<- ggplot()+
  geom_hline(yintercept=0, linetype="dotted") + 
  geom_vline(xintercept=0,  linetype="dotted")

colon_dudi<-dudi.pco(d = sqrt(d_colon), scannf = FALSE, nf = 4)

# plot with variable contribution
#https://www.researchgate.net/post/how_can_i_produce_a_PCoA_biplot_using_R

efit <- envfit(colon_dudi, dat[c('BodySize', 'PLD', 'ParentalMode')], na.rm=T)
varibs<-data.frame(rbind(efit$vectors$arrows, efit$factors$centroids))
varibs$predictors=row.names(varibs)
varibs$predictors<-gsub('ParentalMode', 'PM.',varibs$predictors)

colon_pco<-data.frame(colon_dudi$li,dat)
# filter to see only main groups
colon_pco<-filter(colon_pco, FG %in% c(4,6,1,2))

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

