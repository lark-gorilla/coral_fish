## 11/04/19 
## fish trait cluster paper analyses

#rm(list=ls())
library(cluster)
library(clue)
library(fpc)
library(ggplot2)
library(dendextend)
library(circlize)
library(vegan)
library(caret)
library(GGally)
library(dplyr)
library(reshape2)
library(ade4)
library(ggrepel)
library(adehabitatHR)
library(sf)
library(gridExtra)

# source clVal function
source('C:/coral_fish/scripts/coral_fish/clVal.R')
source('C:/coral_fish/scripts/coral_fish/functions.R')

#################### data clean ##########################
##########################################################

dat<-read.csv('C:/coral_fish/data/Traits/JPN_AUS_RMI_CHK_MLD_TMR_trait_master_opt2.csv', h=T)
# Running on simplist classification of Position trait

# Do we want to include sp. in the analyses? - yes.
# dat[grep('\\.', dat$Species),]

dat<-dat[which(dat$AUS_sp>0 | dat$JPN_sp>0),] # we will focus on Australia and Japan for this prelim

# remove functional duplicates
dup_trait<-paste(dat$BodySize, dat$DepthRange, dat$PLD, dat$Diet, dat$Aggregation, dat$Position, 
                 dat$ParentalMode)
dat<-dat[-which(duplicated(dup_trait)),]
dat$Species<-as.character(dat$Species)
dat[which(dat$Species=='Scarus psittacus'),]$Species<-'Scarus psittacus/spinus' # edit for one

# Including those that default to duplicates via NA after gower dist
dat<-dat[-which(dat$Species=='Caesio sp.'),]
dat<-dat[-which(dat$Species=='Choerodon cyanodus'),]
dat<-dat[-which(dat$Species=='Ostracion immaculatus'),]


row.names(dat)<-dat$Species


## Edit to some trait values from MB 25/10/18
#dat[dat$Species=='Brotula multibarbata',]$DepthRange<-219
dat[dat$Species=='Mobula birostris',]$BodySize<-450
dat[dat$Species=='Amphiprion sandaracinos',]$BodySize<-14

## set ceiling for numeric variables for scaling purposes 04/03/19
dat[dat$PLD>=100 & !is.na(dat$PLD),]$PLD<-100
dat[dat$DepthRange>=200 & !is.na(dat$DepthRange),]$DepthRange<-200

#################### run cluster validation function ##########################
###############################################################################

# analyses on Aus + Jpn combined

#trial to see difference between jaccad and distance measure
dat_aus<-dat[which(dat$AUS_sp>0),]

dat_jpn<-dat[which(dat$JPN_sp>0),]

aus_out<-clVal(data=dat_aus[,3:9], runs=1000, max_cl=20, subs_perc=0.95)
jpn_out<-clVal(data=dat_jpn[,3:9], runs=1000, max_cl=20, subs_perc=0.95)

jac_trial<-clVal(data=dat[,3:9], runs=1000, max_cl=20, subs_perc=0.95)


# outputs saved to memory

################## find optimal n clusters per dataset ########################
###############################################################################

# objects jac_trial (aus+jpn), aus_out and jpn_out loaded from memory

ggpairs(aus_out$stats, columns = 3:7)

#ggpairs(aus_out$stats, columns = 3:7, ggplot2::aes(colour=as.character(k)))

# remember silhouette won't necessarily correlate because it
# has a polynomial relationship with the others?

qplot(data=filter(aus_out$stats, k==6), x=sil, y=jac)

# Australia

a_melt<-melt(aus_out$stats, id.vars=c( 'k', 'runs'))
a_sum<-a_melt%>%group_by(k, variable)%>%
  summarise(mean=mean(value), median=median(value))

ggplot(data=a_melt, aes(x=k, y=value, group=k))+
  geom_boxplot()+facet_wrap(~variable, scales='free_y') # 7 apart from sil:9

ggplot()+
  geom_violin(data=a_melt, aes(x=k, y=value, group=k))+
  geom_point(data=a_sum, aes(x=k, y=mean), color='red', shape=1)+
  geom_line(data=a_sum, aes(x=k, y=mean), color='red')+
  geom_point(data=a_sum, aes(x=k, y=median), color='green', shape=1)+
  geom_line(data=a_sum, aes(x=k, y=median), color='green')+
  scale_x_continuous(breaks=3:20)+
  facet_wrap(~variable, scales='free_y')+
  geom_vline(xintercept = 7, color='cyan')
# 7 for most, 8 or 9 for silhouette
p

#ggsave('C:/coral_fish/outputs/aus_nclust2.png', width=8, height=8)

# Japan

j_melt<-melt(jpn_out$stats, id.vars=c( 'k', 'runs'))
j_sum<-j_melt%>%group_by(k, variable)%>%
  summarise(mean=mean(value), median=median(value))

p<-ggplot()+
  geom_violin(data=j_melt, aes(x=k, y=value, group=k))+
  geom_point(data=j_sum, aes(x=k, y=mean), color='red', shape=1)+
  geom_line(data=j_sum, aes(x=k, y=mean), color='red')+
  geom_point(data=j_sum, aes(x=k, y=median), color='green', shape=1)+
  geom_line(data=j_sum, aes(x=k, y=median), color='green')+
  scale_x_continuous(breaks=3:20)+
  facet_wrap(~variable, scales='free_y')+
  geom_vline(xintercept = 10, color='cyan')
# 10 for all

p
#ggsave('C:/coral_fish/outputs/jpn_nclust2.png', width=8, height=8)

# Combined Australia & Japan

c_melt<-melt(jac_trial$stats, id.vars=c( 'k', 'runs'))
c_sum<-c_melt%>%group_by(k, variable)%>%
  summarise(mean=mean(value), median=median(value))

p<-ggplot()+
  geom_violin(data=c_melt, aes(x=k, y=value, group=k))+
  geom_point(data=c_sum, aes(x=k, y=mean), color='red', shape=1)+
  geom_line(data=c_sum, aes(x=k, y=mean), color='red')+
  geom_point(data=c_sum, aes(x=k, y=median), color='green', shape=1)+
  geom_line(data=c_sum, aes(x=k, y=median), color='green')+
  scale_x_continuous(breaks=3:20)+
  facet_wrap(~variable, scales='free_y')+
  geom_vline(xintercept = 6, color='cyan')
# 6 for most, possibly 9

p
#ggsave('C:/coral_fish/outputs/combined_nclust.png', width=8, height=8)


# can now check individual cluster stability to see which is best

summary(aus_out)

rowMeans(aus_out$jaccard[[7]])
rowMeans(aus_out$jaccard[[8]])
rowMeans(aus_out$jaccard[[9]])

lapply(aus_out$jaccard, function(x){min(rowMeans(data.frame(x)))})

rowMeans(jpn_out$jaccard[[10]])
rowMeans(jpn_out$jaccard[[13]])
rowMeans(jpn_out$jaccard[[16]])

lapply(jpn_out$jaccard, function(x){min(rowMeans(data.frame(x)))})

rowMeans(jac_trial$jaccard[[6]])
rowMeans(jac_trial$jaccard[[9]])

lapply(jac_trial$jaccard, function(x){min(rowMeans(data.frame(x)))})


apply(jac_trial$jaccard[[9]], 1, var)


################## PCA of optimal clustering solution  ########################
###############################################################################

# edit NA values in JPN DepthRange & PLD
#apply(jpn_out$clust_cent, 2, function(x){TRUE %in% is.na(x)})

jpn_out$clust_cent<-jpn_out$clust_cent %>% group_by(kval, jc_match) %>%
  mutate(PLD=replace(PLD, is.na(PLD), mean(PLD, na.rm=T)),
         DepthRange=replace(DepthRange, is.na(DepthRange), mean(DepthRange, na.rm=T))) %>%
  as.data.frame()

# Australia

aus_pca<-pca_vis(rundat=dat_aus[,3:9], clValresult=aus_out, kval=7)

jpn_pca<-pca_vis(rundat=dat_jpn[,3:9], clValresult=jpn_out, kval=10)

png('C:/coral_fish/outputs/aus_pca3.png', width=4, height=12, units ="in", res=600)
grid.arrange(aus_pca[[4]], aus_pca[[5]], aus_pca[[6]], nrow=3)
dev.off()

# add cluster centre for full dataset ~ cluster centroid

enviro.sites.scores<-as.data.frame(scores(aus7_pca, display='sites', scaling=1)) 
enviro.sites.scores$jc_match<-aus7$jc_match
# i've put scaling to 1 for the sites to fit better on the plot

# Now make some plots
enviro.species.scores<-as.data.frame(scores(aus7_pca, display='species'))
enviro.species.scores$Predictors<-colnames(aus7[,4:28])
#enviro.species.scores$Pred_codes<-codez
head(enviro.species.scores)

g<- ggplot()+
  geom_segment(data=NULL, aes(y=-Inf, x=0, yend=Inf, xend=0), linetype='dotted')+
  geom_segment(data=NULL, aes(y=0, x=-Inf, yend=0, xend=Inf), linetype='dotted')+
  geom_point(data=enviro.sites.scores, aes(y=Axis2, x=Axis1, colour=factor(jc_match)))+
  geom_segment(data=enviro.species.scores, aes(y=0, x=0, yend=Comp2, xend=Comp1), arrow=arrow(length=unit(0.3,'lines')), colour='red')+
  theme_classic()+theme(legend.position = "none") 

g<-g+geom_text_repel(data=enviro.species.scores, aes(y=Comp2*2, x=Comp1*2, label=Predictors), segment.size=NA, colour='red')

eig<-aus7_pca$eig
g<- g+scale_y_continuous(paste('PC2', sprintf('(%0.1f%% explained var.)', 100* eig[2]/sum(eig))))+
  scale_x_continuous(paste('PC1', sprintf('(%0.1f%% explained var.)', 100* eig[1]/sum(eig))))

g

## kernel areas

jpn10<-jpn_out$clust_centres[jpn_out$clust_centres$kval==3,]

my_wgt<-c(1,1,1,rep(1/7, 7), rep(1/4, 4), rep(1/6, 6), rep(1/5, 5))

jpn10_pca<-dudi.pca(jpn10[,4:28], col.w=my_wgt, center=T, scale=T, scannf = FALSE, nf = 3)

jpn10<-cbind(jpn10, jpn10_pca$li)

spdf<-SpatialPointsDataFrame(coords=cbind(aus7$Axis1, aus7$Axis2),
                       data=data.frame(jc_match=factor(aus7$jc_match)))

KDE.Surface <- kernelUD(spdf,same4all = T, h=0.05, grid=1000)
kernel.area(KDE.Surface, percent = c(25, 50, 75, 95), standardize = T)

KDE.UD <- getverticeshr(KDE.Surface, percent = 99)

KDE.Surface <- kernelUD(spdf, h='LSCV', same4all = T)

KDE.Surface <- kernelUD(data.frame(TripCoords[,1], TripCoords[,2]), id=Trip$ID, h=(Scale * 1000), grid=100, extent=BExt, same4all=FALSE)
KDE.UD <- getverticeshr(KDE.Surface, lev = UDLev)
KDE.Sp1 <- kver2spol(KDE.UD)


# plots

ppp+geom_point(data=aus7[1:(nrow(aus7)-kval),], aes(x=Axis1, y=Axis2, colour=factor(jc_match)), shape=3)+geom_sf(data=pols, colour=rainbow(7), alpha=0.5)+scale_color_manual(values=rainbow(7))