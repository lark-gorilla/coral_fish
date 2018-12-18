# 02/11/18
# Code to determine optimal clustering method and n clusters for Japan-Australia trait dataset

library(dplyr)
library(ggplot2)
library(factoextra)
library(NbClust)
library(car)
library(FD)
library(cluster)
library(gridExtra)
library(fpc)

dat<-read.csv('C:/coral_fish/data/Traits/JPN_AUS_RMI_trait_master.csv', h=T)

dat<-dat[which((dat$JPN_sp+dat$AUS_sp)>0),] # keeps AUS and JPN combined, can splt and test later

## Edit to some trait values from MB 25/10/18

dat[dat$Species=='Brotula multibarbata',]$DepthRange<-219
dat[dat$Species=='Manta birostris',]$BodySize<-450
dat[dat$Species=='Amphiprion sandaracinos',]$BodySize<-14

########

### The warnings of 0 dist in is.euclid are due to functionally identical species

dat_temp<-dat
dat_temp<-mutate(dat_temp, key=paste(ThermalAffinity, BodySize,DepthRange, PLD, Diet, Aggregation, Position, ParentalMode))
dat_temp[which(duplicated(dat_temp$key)),]
#have a look at each
filter(dat_temp, key=='tropical 290 100 0 Piscivore solitary Demersal Live bearers' |
         key=='tropical 30 23 29 Herbivore solitary UpperBenthic scatterers' |
         key=='tropical 83 297 25.2 Predator solitary UpperBenthic scatterers')

# remove duplicates
dat<-dat[-which(dat$Species == 'Orectolobus sp.' |dat$Species == 'Scarus spinus' |
                  dat$Species == 'Variola sp.' ),]


#check if scaling is needed for numeric variables

summary(dat[,3:5]) # looks ok, but some big values

pairs(dat[,2:9]) # some outliers need to look at

qplot(data=dat, x=BodySize, geom='histogram')

dat[dat$BodySize>350,] # a ray and Manta, and something else

# dat2 is outlier-removed data
dat2<-filter(dat, BodySize<350 | is.na(BodySize))

qplot(data=dat, x=DepthRange, geom='histogram')

dat[which(dat$DepthRange>600),] # 3 species
dat2<-filter(dat2, DepthRange<=600 | is.na(DepthRange)) # allows species with NA depth through

qplot(data=dat, x=PLD, geom='histogram')

dat[which(dat$PLD>200),] # 10 species
dat2<-filter(dat2, PLD<200| is.na(PLD))

vif(lm(1:nrow(dat)~ThermalAffinity+BodySize+DepthRange+PLD+Diet+
         Aggregation+Position+ParentalMode, data=dat, na.action=na.omit))
#Collinearity looks ok

# Ok so proceed to clustering 

# check if distances are the same calced via gowdis() and daisy()

dist1<-gowdis(dat[,2:9])
dist2<-daisy(dat[,2:9],metric = "gower", stand = FALSE)
# there is essentially no difference between theses functions
# gowdis is more flexible as you can include weights and it better handles nominal data
# but neither apply to this analyses.

dist1==dist2 # not the same

diff(c(as.matrix(dist1)[2],as.matrix(dist2)[2]))
# differences tiny, use gowdis

# now check differences in distance matrices when removing outliers

dist1<-gowdis(dat[,2:9])
dist2<-gowdis(dat2[,2:9])

# what about transformation
dat4<-dat
dat4$BodySize<-log10(dat4$BodySize)
dat4$DepthRange<-log10(dat4$DepthRange)
dat4$PLD<-log10(dat4$PLD)
dist4<-gowdis(dat4[,2:9])
# can be done with daisy
dist4d<-daisy(dat[,2:9],metric = "gower", type = list(logratio = c(2,3,4)))

#optimal solution would be transformation on outlier removed data
dat5<-dat2
dat5$BodySize<-log10(dat5$BodySize)
dat5$DepthRange<-log10(dat5$DepthRange)
dat5$PLD<-log10(dat5$PLD)
dist5<-gowdis(dat5[,2:9])


pd1<-fviz_dist(dist1, order = TRUE, show_labels = TRUE, lab_size = 4)
pd2<-fviz_dist(dist2, order = TRUE, show_labels = TRUE, lab_size = 4)
pd4<-fviz_dist(dist4, order = TRUE, show_labels = TRUE, lab_size = 4)
pd4d<-fviz_dist(dist4d, order = TRUE, show_labels = TRUE, lab_size = 4)

grid.arrange(pd1, pd2, pd4) # Definately some differences. - but remember this uses hclust internally

grid.arrange(pd4, pd4d) # Just to check differnce in logging data

# OK so outliers have an impact but so does normality, see: https://www.r-bloggers.com/clustering-mixed-data-types-in-r/
# scaling is done internally by Gower, also continuous variables are already on similar scale.
# If we want to keep outlier species in then we need to think about best
# clustering method and/or algorithm OR try variable transformation.

#Now look at categorial levels and NA species rows

table(dat$ThermalAffinity)
table(dat$Diet)
table(dat$Aggregation)
table(dat$Position)
table(dat$ParentalMode)

# some variables have levels with only a few species which could
# be problematic?

dat[which(apply(is.na(dat[,2:9]), 1, sum)>1),]
# only ~ 10 species with two NAs, only 1, (Ostracion) that has 3.

# Ok lets make a quick PCoA to visulise the trait space

dist_euc <- cailliez(dist1)
dist_euc2 <- lingoes(dist1)

pd2<-fviz_dist(dist_euc, order = TRUE, show_labels = TRUE, lab_size = 4)

grid.arrange(pd1, pd2) # check if euclid corrected distance 
# still is same as Gower product - kinda - not sure if this is necessary

ppp <- ggplot() + coord_fixed() + 
  labs(x="Comp1, Axis1", y="Comp2, Axis2") +
  geom_hline(yintercept=0, col="darkgrey") + 
  geom_vline(xintercept=0, col="darkgrey")

# compare the different distance matrixes in PCoA space

p1<-ppp+geom_point(data=dudi.pco(d = cailliez(dist1),
     scannf = FALSE, nf = 2)$li, aes(x=A1, y=A2))+ggtitle('raw')
p2<-ppp+geom_point(data=dudi.pco(d = cailliez(dist2),
    scannf = FALSE, nf = 2)$li, aes(x=A1, y=A2))+ggtitle('outlier removed')
p3<-ppp+geom_point(data=dudi.pco(d = cailliez(dist4),
    scannf = FALSE, nf = 2)$li, aes(x=A1, y=A2))+ggtitle('log gowdis')
p4<-ppp+geom_point(data=dudi.pco(d = cailliez(dist4d),
    scannf = FALSE, nf = 2)$li, aes(x=A1, y=A2))+ggtitle('log daisy')
p5<-ppp+geom_point(data=dudi.pco(d = cailliez(dist5),
    scannf = FALSE, nf = 2)$li, aes(x=A1, y=A2))+ggtitle('log outlier removed')

grid.arrange(p1,p2,p3,p4,p5) # ok so differences are shown

#doesn't seem to be much in the way of clustering - possibly an impact 
# of correcting the distance data?

# check clustering in each variable individually

p1<-ppp+geom_point(data=dudi.pco(d = cailliez(gowdis(dat[,3:9])),
    scannf = FALSE, nf = 2)$li, aes(x=A1, y=A2))+ggtitle(paste('dropped',names(dat)[2]))
p2<-ppp+geom_point(data=dudi.pco(d = cailliez(gowdis(dat[,c(2, 4:9)])),
    scannf = FALSE, nf = 2)$li, aes(x=A1, y=A2))+ggtitle(paste('dropped',names(dat)[3]))
p3<-ppp+geom_point(data=dudi.pco(d = cailliez(gowdis(dat[,c(2:3, 5:9)])),
    scannf = FALSE, nf = 2)$li, aes(x=A1, y=A2))+ggtitle(paste('dropped',names(dat)[4]))
p4<-ppp+geom_point(data=dudi.pco(d = cailliez(gowdis(dat[,c(2:4, 6:9)])),
    scannf = FALSE, nf = 2)$li, aes(x=A1, y=A2))+ggtitle(paste('dropped',names(dat)[5]))
p5<-ppp+geom_point(data=dudi.pco(d = cailliez(gowdis(dat[,c(2:5, 7:9)])),
    scannf = FALSE, nf = 2)$li, aes(x=A1, y=A2))+ggtitle(paste('dropped',names(dat)[6]))
p6<-ppp+geom_point(data=dudi.pco(d = cailliez(gowdis(dat[,c(2:6, 8:9)])),
    scannf = FALSE, nf = 2)$li, aes(x=A1, y=A2))+ggtitle(paste('dropped',names(dat)[7]))
p7<-ppp+geom_point(data=dudi.pco(d = cailliez(gowdis(dat[,c(2:7, 9)])),
    scannf = FALSE, nf = 2)$li, aes(x=A1, y=A2))+ggtitle(paste('dropped',names(dat)[8]))
p8<-ppp+geom_point(data=dudi.pco(d = cailliez(gowdis(dat[,2:8])),
    scannf = FALSE, nf = 2)$li, aes(x=A1, y=A2))+ggtitle(paste('dropped',names(dat)[9]))

grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8, ncol=4, nrow=2)
# dropping parental mode looks nice..

### Now look at outlier removed data

# check clustering in each variable individually

p1<-ppp+geom_point(data=dudi.pco(d = cailliez(gowdis(dat2[,3:9])),
scannf = FALSE, nf = 2)$li, aes(x=A1, y=A2))+ggtitle(paste('dropped',names(dat2)[2]))
p2<-ppp+geom_point(data=dudi.pco(d = cailliez(gowdis(dat2[,c(2, 4:9)])),
scannf = FALSE, nf = 2)$li, aes(x=A1, y=A2))+ggtitle(paste('dropped',names(dat2)[3]))
p3<-ppp+geom_point(data=dudi.pco(d = cailliez(gowdis(dat2[,c(2:3, 5:9)])),
scannf = FALSE, nf = 2)$li, aes(x=A1, y=A2))+ggtitle(paste('dropped',names(dat2)[4]))
p4<-ppp+geom_point(data=dudi.pco(d = cailliez(gowdis(dat2[,c(2:4, 6:9)])),
scannf = FALSE, nf = 2)$li, aes(x=A1, y=A2))+ggtitle(paste('dropped',names(dat2)[5]))
p5<-ppp+geom_point(data=dudi.pco(d = cailliez(gowdis(dat2[,c(2:5, 7:9)])),
scannf = FALSE, nf = 2)$li, aes(x=A1, y=A2))+ggtitle(paste('dropped',names(dat2)[6]))
p6<-ppp+geom_point(data=dudi.pco(d = cailliez(gowdis(dat2[,c(2:6, 8:9)])),
scannf = FALSE, nf = 2)$li, aes(x=A1, y=A2))+ggtitle(paste('dropped',names(dat2)[7]))
p7<-ppp+geom_point(data=dudi.pco(d = cailliez(gowdis(dat2[,c(2:7, 9)])),
scannf = FALSE, nf = 2)$li, aes(x=A1, y=A2))+ggtitle(paste('dropped',names(dat2)[8]))
p8<-ppp+geom_point(data=dudi.pco(d = cailliez(gowdis(dat2[,2:8])),
scannf = FALSE, nf = 2)$li, aes(x=A1, y=A2))+ggtitle(paste('dropped',names(dat2)[9]))

grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8, ncol=4, nrow=2)

# same story as raw data - ok so check what heppens if we reclass the categorial outliers

# quick comparison of dropped parental model between different datasets

p1<-ppp+geom_point(data=dudi.pco(d = cailliez(gowdis(dat[,2:8])),
    scannf = FALSE, nf = 2)$li, aes(x=A1, y=A2))+ggtitle('raw')
p2<-ppp+geom_point(data=dudi.pco(d = cailliez(gowdis(dat2[,2:8])),
    scannf = FALSE, nf = 2)$li, aes(x=A1, y=A2))+ggtitle('outlier removed')
p3<-ppp+geom_point(data=dudi.pco(d = cailliez(gowdis(dat4[,2:8])),
    scannf = FALSE, nf = 2)$li, aes(x=A1, y=A2))+ggtitle('log gowdis')
p5<-ppp+geom_point(data=dudi.pco(d = cailliez(gowdis(dat5[,2:8])),
    scannf = FALSE, nf = 2)$li, aes(x=A1, y=A2))+ggtitle('log outlier removed')

grid.arrange(p1,p2,p3,p5) #maybe don't use log transformations
# Maybe don't get too hung up on PCO visualisation

# Update to PCO visusalisation

# Investigating negative eigenvalues

pco1<-dudi.pco(d = gowdis(dat[, c(2:5, 7, 9)]), scannf = T)
pco2<-dudi.pco(d = sqrt(gowdis(dat[, c(2:5, 7, 9)])), scannf = T)
pco3<-dudi.pco(d = cailliez(gowdis(dat[, c(2:5, 7, 9)])), scannf = T)
pco4<-dudi.pco(d = lingoes(gowdis(dat[, c(2:5, 7, 9)])), scannf = T)

# dudi.pco only returns positive eigenvals
nrow(dat[, c(2:5, 7, 9)])
length(eigenvals(pco1))
length(eigenvals(pco2))
length(eigenvals(pco3))
length(eigenvals(pco4))

# plot with cmdscale to get all eigenvals (even the negative ones)

pco1_c<-cmdscale( gowdis(dat[, c(2:5, 7, 9)]), k=652, eig = TRUE)
pco2_c<-cmdscale( sqrt(gowdis(dat[, c(2:5, 7, 9)])), k=652, eig = TRUE)
pco3_c<-cmdscale( cailliez(gowdis(dat[, c(2:5, 7, 9)])), k=652, eig = TRUE)
pco4_c<-cmdscale( lingoes(gowdis(dat[, c(2:5, 7, 9)])), k=652, eig = TRUE)

df<-data.frame(eigenvals=rep(1:653,4), eigenval_unit=c(pco1_c$eig, 
                                                       pco2_c$eig,
                                                       pco3_c$eig,
                                                       pco4_c$eig),
               Method=c(rep('raw', 653), rep('sqrt', 653), 
                        rep('cailliez', 653), rep('lingoes', 653)))

qplot(data=df, x=eigenvals, y=eigenval_unit, colour=Method, geom='line')+
  facet_wrap(~Method, scales='free_y')+
  geom_point(data=df[df$eigenval_unit<0,],  aes(x=eigenvals, y=eigenval_unit), shape=1)
# same as above, just different scaling
df<-data.frame(eigenvals=c(1:176, 1:436, 1:647, 1:642), eigenval_unit=c(pco1$eig, 
                                                       pco2$eig,
                                                       pco3$eig,
                                                       pco4$eig),
               Method=c(rep('raw', 176), rep('sqrt', 436), 
                        rep('cailliez', 647), rep('lingoes', 642)))

qplot(data=df, x=eigenvals, y=eigenval_unit, colour=Method, geom='line')+
  facet_wrap(~Method, scales='free_y')

# No correction method fully removes negative eigenvals

# check proportion of variance explained by first 4 axes

sum(eigenvals(pco1_c)[1:4])/sum(eigenvals(pco1_c)[which(eigenvals(pco1_c)>0)])
sum(eigenvals(pco1)[1:4])/sum(eigenvals(pco1))

sum(eigenvals(pco2)[1:4])/sum(eigenvals(pco2))
sum(eigenvals(pco3)[1:4])/sum(eigenvals(pco3))
sum(eigenvals(pco4)[1:4])/sum(eigenvals(pco4))

# basically the cailliez and lingoez correction aren't really doing anything.
# best stick with raw or sqrt for PCO plotting.
# FYI plots aren't that much different for full dataset

# See here for info 
#http://r-sig-ecology.471788.n2.nabble.com/Variability-explanations-for-the-PCO-axes-as-in-Anderson-and-Willis-2003-td6429547.html

# validation of n PCoA axis to represent original data, as per Mouillot et al. (2014)

d1<-gowdis(dat[, c(2:5, 7, 9)])
d2<-dist(pco1$li) #(assuming pco1 has k=4)
qplot(x=d1, y=d2)
# Do mantel test to get correlation coefficient and p val
mantel(d1, d2)

# could compare results between other corrections to see best r2

table(dat2$ThermalAffinity)
table(dat2$Diet)
table(dat2$Aggregation)
table(dat2$Position)
table(dat2$ParentalMode)

dat2_trial<-dat2

dat2_trial$Position<-as.character(dat2_trial$Position)

dat2_trial[dat2_trial$Position=='AlgaeAssociated' |
             dat2_trial$Position== 'EchinodermAssociated'|
             dat2_trial$Position== 'SandAssociated',]$Position<-'other'

dat2_trial$ParentalMode<-as.character(dat2_trial$ParentalMode)

dat2_trial[which(dat2_trial$ParentalMode=='demersal' |
             dat2_trial$ParentalMode== 'Viviparous'),]$ParentalMode<-'brooders'

p9<-ppp+geom_point(data=dudi.pco(d = cailliez(gowdis(dat2_trial[,c(2:9)])),
                                 scannf = FALSE, nf = 2)$li, aes(x=A1, y=A2))

grid.arrange(p7,p8,p9) # Hmm not sure whats going on here, should be better than that..

## maybe Position and Parental Mode just have lots of different combinations

sort(table(paste(dat2$Position, dat2$ParentalMode)))


length(table(paste(dat2$Position, dat2$ParentalMode)))

nrow(expand.grid(paste(dat2$Position, dat2$ParentalMode)))

#####################################
#    Hierarchical clustering
#####################################

mydist<-gowdis(dat[,c(2:3)])

hc_av<-hclust(d=mydist, method='average')
hc_si<-hclust(d=mydist,  method='single')
hc_co<-hclust(d=mydist,  method='complete')
hc_w1<-hclust(d=mydist,  method='ward.D')
hc_w2<-hclust(d=mydist,  method='ward.D2')
hc_mc<-hclust(d=mydist,  method='mcquitty')
hc_me<-hclust(d=mydist,  method='median')
hc_ce<-hclust(d=mydist,  method='centroid')

# number of clusters

# nbclust not useful for gower distance data

outdat<-NULL

for(i in 2:20){
  
  mnz<-list(cluster.stats(mydist, cutree(hc_av, k=i)),
            cluster.stats(mydist, cutree(hc_si, k=i)),
            cluster.stats(mydist, cutree(hc_co, k=i)),
            cluster.stats(mydist, cutree(hc_w1, k=i)),
            cluster.stats(mydist, cutree(hc_w2, k=i)),
            cluster.stats(mydist, cutree(hc_mc, k=i)),
            cluster.stats(mydist, cutree(hc_me, k=i)),
            cluster.stats(mydist, cutree(hc_ce, k=i)))
    
  
  out<-data.frame(nclust=i, method=c('average','single','complete','ward.D',
                                     'ward.D2','mcquitty','median','centroid'),
                   av_sil_width=c(mnz[[1]]$avg.silwidth,mnz[[2]]$avg.silwidth,mnz[[3]]$avg.silwidth,
                                  mnz[[4]]$avg.silwidth,mnz[[5]]$avg.silwidth,mnz[[6]]$avg.silwidth,
                                  mnz[[7]]$avg.silwidth,mnz[[8]]$avg.silwidth),
                  wb_ratio=c(mnz[[1]]$wb.ratio,mnz[[2]]$wb.ratio,mnz[[3]]$wb.ratio,
                             mnz[[4]]$wb.ratio,mnz[[5]]$wb.ratio,mnz[[6]]$wb.ratio,
                             mnz[[7]]$wb.ratio,mnz[[8]]$wb.ratio))
  
  outdat<-rbind(outdat, out)
  print(i)
}

# Plot sihouette width (higher is better)

qplot(data=outdat, x=nclust, y=av_sil_width, colour=method)+labs(x="Number of clusters", y= "Silhouette Width")+
  geom_line() +scale_x_continuous(breaks=2:20)


ggplot(data=outdat, aes(x=nclust))+labs(x="Number of clusters", y= "Silhouette Width/WB ratio")+
  geom_line(aes(y=av_sil_width), colour='red')+geom_point(aes(y=av_sil_width), colour='red')+
    geom_line(aes(y=wb_ratio), colour='black')+geom_point(aes(y=wb_ratio), colour='black')+
    scale_x_continuous(breaks=2:20)+facet_wrap(~method, ncol=2)

fviz_dend(hc_av, cex = 0.6, rect = TRUE)

fviz_silhouette(silhouette(kmeans(mydist, 8)$cluster, dist=mydist))
fviz_silhouette(silhouette(dist=mydist,cutree(hclust(d=mydist, method='average'), k=8)))

#trial of different combinations of variables


bigout<-NULL
for(i in 2:8)
{
  dfz<-combn(dat[c(2:9)], i, simplify = F)
  
  for(j in 1:length(dfz))
  {  
  mydist<-gowdis(dfz[[j]])
  
  if(class(try(hclust(d=mydist, method='average'), silent = T))== "try-error"){
    print(paste('ERROR', names(dfz[[j]])));next}
    
  hc_av<-hclust(d=mydist, method='average')
  hc_si<-hclust(d=mydist,  method='single')
  hc_co<-hclust(d=mydist,  method='complete')
  hc_w1<-hclust(d=mydist,  method='ward.D')
  hc_w2<-hclust(d=mydist,  method='ward.D2')
  hc_mc<-hclust(d=mydist,  method='mcquitty')
  hc_me<-hclust(d=mydist,  method='median')
  hc_ce<-hclust(d=mydist,  method='centroid')
  
  out1<-NULL
  for(m in 4:14)
  {
    mnz<-list(cluster.stats(mydist, cutree(hc_av, k=m)),
              cluster.stats(mydist, cutree(hc_si, k=m)),
              cluster.stats(mydist, cutree(hc_co, k=m)),
              cluster.stats(mydist, cutree(hc_w1, k=m)),
              cluster.stats(mydist, cutree(hc_w2, k=m)),
              cluster.stats(mydist, cutree(hc_mc, k=m)),
              cluster.stats(mydist, cutree(hc_me, k=m)),
              cluster.stats(mydist, cutree(hc_ce, k=m)))
    
    corz<-c(cor(mydist, dist(cutree(hc_av, k=m))),
           cor(mydist, dist(cutree(hc_si, k=m))),
           cor(mydist, dist(cutree(hc_co, k=m))),
           cor(mydist, dist(cutree(hc_w1, k=m))),
           cor(mydist, dist(cutree(hc_w2, k=m))),
           cor(mydist, dist(cutree(hc_mc, k=m))),
           cor(mydist, dist(cutree(hc_me, k=m))),
           cor(mydist, dist(cutree(hc_ce, k=m))))
    
    
    out<-data.frame(nvar=i,nclust=m,
                    method=c('average','single','complete','ward.D',
                                       'ward.D2','mcquitty','median','centroid'),
                    av_sil_width=c(mnz[[1]]$avg.silwidth,mnz[[2]]$avg.silwidth,mnz[[3]]$avg.silwidth,
                                   mnz[[4]]$avg.silwidth,mnz[[5]]$avg.silwidth,mnz[[6]]$avg.silwidth,
                                   mnz[[7]]$avg.silwidth,mnz[[8]]$avg.silwidth),
                    wb_ratio=c(mnz[[1]]$wb.ratio,mnz[[2]]$wb.ratio,mnz[[3]]$wb.ratio,
                               mnz[[4]]$wb.ratio,mnz[[5]]$wb.ratio,mnz[[6]]$wb.ratio,
                               mnz[[7]]$wb.ratio,mnz[[8]]$wb.ratio),
                    clust_d_cor=corz,
                    vars=paste(names(dfz[[j]]), collapse=" "))
    
  out1<-rbind(out1, out)
  }
  
  #print(names(dfz[[j]]))
  #print(out1)
  bigout<-rbind(bigout, out1)
  }
  print(i)
}    

#write out
#write.csv(bigout, 'c:/coral_fish/data/Traits/cluster_combinations.csv', quote=F, row.names=F)

bigout<-read.csv('c:/coral_fish/data/Traits/cluster_combinations.csv', h=T)

# visualise
# 7 variables...
qplot(data=bigout[bigout$nvar==7,], y=av_sil_width, x=nclust, colour=vars, geom='line' )+
  guides(colour=guide_legend(ncol=1))+facet_wrap(~method)

mydist<-gowdis(dat[,c(2:7,9)])

fviz_silhouette(silhouette(dist=mydist,cutree(hclust(d=mydist, method='ward.D2'), k=10)))

## 6 variables
qplot(data=bigout[bigout$nvar==6,], y=av_sil_width, x=nclust, colour=vars, geom='line' )+
  guides(colour=guide_legend(ncol=1))+facet_wrap(~method)


mydist<-gowdis(dat[,c(2:5,7,9)]) # dropped diet and position

# best option mcquitty =0.588 @ 10 clusters
fviz_silhouette(silhouette(dist=mydist,cutree(hclust(d=mydist, method='mcquitty'), k=10)))

# Identical with Ward.D2 @ 10 clusters but no clear peak at k=10
fviz_silhouette(silhouette(dist=mydist,cutree(hclust(d=mydist, method='ward.D2'), k=10)))

# another contender, 'better' variables but weaker

mydist<-gowdis(dat[,c(3:7,9)]) # dropped diet and position

# best option ward.D2 =0.49 @ 6 clusters
fviz_silhouette(silhouette(dist=mydist,cutree(hclust(d=mydist, method='ward.D2'), k=6)))

# second option average =0.42 @ 7 clusters
fviz_silhouette(silhouette(dist=mydist,cutree(hclust(d=mydist, method='average'), k=7)))

# joint second option average =0.43 @ 8 clusters
fviz_silhouette(silhouette(dist=mydist,cutree(hclust(d=mydist, method='complete'), k=8)))

# have a look at option 1
fviz_dend(hclust(d=mydist, method='ward.D2'), cex = 0.6, rect = TRUE)

fviz_dend(hclust(d=mydist, method='ward.D2'), k=6, cex = 0.6, rect = TRUE)

fviz_dend(hclust(d=mydist, method='ward.D2'), k=6, color_labels_by_k = TRUE, type='circular')

# Validation using cophenetic distances of cluster solutions
# See here for explanation 
#https://people.revoledu.com/kardi/tutorial/Clustering/Cophenetic.htm

hc_av<-hclust(d=mydist, method='average')
hc_si<-hclust(d=mydist,  method='single')

fviz_dend(hc_av, cex = 0.6, rect = TRUE)

# cophenetic distances of species based on cluster nodes
hc_av_cop<-cophenetic(hc_av)

# distances calculated between species based on species assigned to 1 of 10 clusters
cutty<-cutree(hc_av, k=10)
hc_av_group<-dist(cutty)
# so I think this is the correct approach. By calculating dist on each cluster
# as a numeirc object it gives clusters with bigger differences in their 'name'/number
# a bigger distance e.g. further apart. I think this is correct, - no its not!!

mydist<-gowdis(dat[,c(3:7,9)]) 
attr(mydist, 'Labels')<-paste(cutty) # bodge
hc_co<-hclust(d=mydist,  method='complete')
cutty<-cutree(hc_co, k=6)
plot(hc_co)
rect.hclust(hc_co, k=6, border=c(6,3,5,2,4,1)) # make colours line up
# ok so labelling from cutree doesnt follow dendrogram plot clustering structure

mantel(mydist, hc_av_cop) # gives p value to r2 correlation metric

mantel(mydist, hc_av_group) # unsurprisingly grouped data does not correlate as well

####### Validation attempt / alternate approach  

library(dendextend) # watch out for functions sharing the same name e.g. cutree

mydist<-gowdis(dat[,c(3:7,9)]) # dropped diet and position

hc_av<-hclust(d=mydist, method='average') %>% as.dendrogram()
hc_si<-hclust(d=mydist,  method='single')%>% as.dendrogram()
hc_co<-hclust(d=mydist,  method='complete')%>% as.dendrogram()
hc_w1<-hclust(d=mydist,  method='ward.D')%>% as.dendrogram()
hc_w2<-hclust(d=mydist,  method='ward.D2')%>% as.dendrogram()
hc_mc<-hclust(d=mydist,  method='mcquitty')%>% as.dendrogram()
hc_me<-hclust(d=mydist,  method='median')%>% as.dendrogram()
hc_ce<-hclust(d=mydist,  method='centroid')%>% as.dendrogram()

dend1234 <- dendlist('average'=hc_av,'single'=hc_si,'complete'=hc_co,'ward.D'=hc_w1,
'ward.D2'=hc_w2,'mcquitty'=hc_mc,'median'=hc_me,'centroid'=hc_ce)

cor_d<-cor.dendlist(dend1234)
# cophenetic correlation between different approaches

library(corrplot)
corrplot(cor_d, "pie", "lower")
# Ward alorithms seem least correlated with other methods

dend_FM <- dendlist('average'=cutree(hc_av, k=m),'single'=cutree(hc_si, k=m),'complete'=cutree(hc_co, k=m),'ward.D'=cutree(hc_w1, k=m),
                     'ward.D2'=cutree(hc_w2, k=m),'mcquitty'=cutree(hc_mc, k=m),'median'=cutree(hc_me, k=m),'centroid'=cutree(hc_ce, k=m))


# wow..
tanglegram(hc_av, hc_si)

FM_index(cutree(hc_av, k=3), cutree(hc_si, k=3), assume_sorted_vectors = T) 



### To view clusters with env data gradients

some_col_func <- function(n) rev(colorspace::heat_hcl(n, c = c(80, 30), l = c(30, 90), power = c(1/5, 1.5)))

library(gplots)
heatmap.2(as.matrix(iris2), 
          main = "Heatmap for the Iris data set",
          srtCol = 20,
          dendrogram = "row",
          Rowv = dend,
          Colv = "NA", 
          trace="none",          
          margins =c(5,0.1),      
          key.xlab = "Cm",
          denscol = "grey",
          density.info = "density",
          RowSideColors = rev(labels_colors(dend)), 
          col = some_col_func)

# http://www.ams.med.uni-goettingen.de/download/Steffen-Unkel/cluster1.html



## Or try this
temp<-dat[,c(3:5)]
row.names(temp)<-dat$Species
x  <- as.matrix(temp)

# d3heatmap(x)
# now let's spice up the dendrograms a bit:
Rowv  <- x %>% dist %>% hclust %>% as.dendrogram %>%
  set("branches_k_color", k = 3) %>% set("branches_lwd", 4) %>%
  ladderize
#    rotate_DendSer(ser_weight = dist(x))
Colv  <- x %>% t %>% dist %>% hclust %>% as.dendrogram %>%
  set("branches_k_color", k = 2) %>% set("branches_lwd", 4) %>%
  ladderize
#    rotate_DendSer(ser_weight = dist(t(x)))

  library(d3heatmap)
d3heatmap(x, Rowv = Rowv, Colv = Colv)

