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

dat<-read.csv('C:/coral_fish/data/Traits/JPN_AUS_RMI_trait_master.csv', h=T)

dat<-dat[which((dat$JPN_sp+dat$AUS_sp)>0),] # keeps AUS and JPN combined, can splt and test later

#check if scaling is needed for numeric variables

summary(dat[,3:5]) # looks ok, but some big values

pairs(dat[,2:9]) # some outliers need to look at

qplot(data=dat, x=BodySize, geom='histogram')

dat[dat$BodySize>350,] # a ray and Manta, and something else

# dat2 is outlier-removed data
dat2<-filter(dat, BodySize<350)

qplot(data=dat, x=DepthRange, geom='histogram')

dat[which(dat$DepthRange>600),] # 3 species
dat2<-filter(dat2, DepthRange<=600)

qplot(data=dat, x=PLD, geom='histogram')

dat[which(dat$PLD>200),] # 10 species
dat2<-filter(dat2, PLD<200)

vif(lm(1:nrow(dat)~ThermalAffinity+BodySize+DepthRange+PLD+Diet+
         Aggregation+Position+ParentalMode, data=dat, na.action=na.omit))
#Collinearity looks ok

# Ok so proceed to clustering 

# check if distances are the same calced via gowdis() and daisy()

dist1<-gowdis(dat[,2:9])
dist2<-daisy(dat[,2:9],metric = "gower", stand = FALSE)

dist1==dist2 # not the same

diff(c(as.matrix(dist1)[2],as.matrix(dist2)[2]))
# differences tiny, use gowdis

# now check differences in distance matrices when removing outliers

dist1<-gowdis(dat[,2:9])
dist2<-gowdis(dat2[,2:9])
# and what happens if we scale continuous variables also
dat3<-dat2
dat3[,3:5]<-scale(dat3[,3:5])
dist3<-gowdis(dat3[,2:9])

pd1<-fviz_dist(dist1, order = TRUE, show_labels = TRUE, lab_size = 4)
pd2<-fviz_dist(dist2, order = TRUE, show_labels = TRUE, lab_size = 4)
pd3<-fviz_dist(dist3, order = TRUE, show_labels = TRUE, lab_size = 4)

grid.arrange(pd1, pd2, pd3) # look quite different.

# OK so outliers have an impact but scaling doesn't
# seem to that much (as continouos variables are already on similar scale)
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

grid.arrange(pd1, pd2) # check if euclid corrected bdistance 
# still is same as Gower product - kinda

pc2<-dudi.pco(d = dist_euc, scannf = FALSE, nf = 4)

ppp <- ggplot() + coord_fixed() + 
  labs(x="Comp1, Axis1", y="Comp2, Axis2") +
  geom_hline(yintercept=0, col="darkgrey") + 
  geom_vline(xintercept=0, col="darkgrey")

pc2.dfs <- data.frame(pc2$li, dat)

#plot pcoa
ppp + geom_point(data=pc2.dfs, aes(x=A1, y=A2))

ppp + geom_point(data=pc2.dfs, aes(x=A3, y=A4))
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

