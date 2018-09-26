# 20/09/2018

# Initial explore of fish and coral survey and trait data

library(dplyr)
library(readxl)
library(ggplot2)
library(reshape2)
library(gridExtra)
#https://daijiang.name/en/2014/05/11/functional-diversity-in-r/
library(FD)
library(ade4)

#setwd("~/leeds_postdoc")
setwd("M:/coral_fish")


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

# Now for coral data

# read in coral and fish abundance and trait data

coraldat<-read.csv('data/RMI/RMI_corals_A.csv', h=T)

coraltrait<-read.csv('data/Traits/RMI_coralTraits4May15.csv', h=T)

# replace '.' between family & species name in abundance data
names(coraldat)<-gsub('\\.', '\\ ', names(coraldat))

#test names line up

n1<-names(coraldat)[-1] # without site col name
n2<-coraltrait$Species

# FYI more species in trait dataset than abundance dataset

ingrid<-expand.grid(n1, n2)
ingrid$Var1<-as.character(ingrid$Var1)
ingrid$Var2<-as.character(ingrid$Var2)

ingrid$test=ingrid$Var1==ingrid$Var2

ingrid_out<- ingrid %>% group_by(Var1) %>%
  summarise(has_trait=TRUE%in%test)

filter(ingrid_out, has_trait==F)
# two 'sp.' species had the '.' removed by gsub above
names(coraldat)[names(coraldat)=='Acropora sp ']<-'Acropora sp.'
names(coraldat)[names(coraldat)=='Alveopora sp ']<-'Alveopora sp.'

# prepare trait data for Gower distance calculation

# filter to only include traits that are in abundance data

fishtrait2<-filter(fishtrait, Species%in% names(fishdat))

row.names(fishtrait2)<-fishtrait2$Species

# drop trait 'PD50' as seems to cock up calculaion
fishtrait2<-fishtrait2[,c(3:5,7:length(fishtrait2))]
# could also drop 'Function' as coarse representative of 'Food'?

fishtrait_dist<-gowdis(fishtrait2)

is.euclid(fishtrait_dist) # not euclidean so needs transforming for use in pcoa


fishtrait_dist2 <- cailliez(fishtrait_dist)

pc1<-pcoa(fishtrait_dist2)

pc2<-dudi.pco(d = fishtrait_dist2, scannf = FALSE, nf = 4)

screeplot(pc2)
biplot(pc2) # yikes

#https://www.r-bloggers.com/7-functions-to-do-metric-multidimensional-scaling-in-r/
ppp <- ggplot() + coord_fixed() + 
  labs(x="Comp1, Axis1", y="Comp2, Axis2") +
  geom_hline(yintercept=0, col="darkgrey") + 
  geom_vline(xintercept=0, col="darkgrey")


pc2.dfs <- data.frame(pc2$li, fishtrait2)
pc2.dfs$Species<-row.names(pc2.dfs)
pc2.dfs<-left_join(pc2.dfs, fishtrait[,1:2], by='Species')

#plot gloab trait space by family
ppp + geom_point(data=pc2.dfs, aes(x=A1, y=A2, colour=Family))

# and for each trait
p1<-ppp + geom_point(data=pc2.dfs, aes(x=A1, y=A2, colour=HomeRange))+scale_colour_gradientn(colors=topo.colors(10))
p2<-ppp + geom_point(data=pc2.dfs, aes(x=A1, y=A2, colour=MaxLength))+scale_colour_gradientn(colors=topo.colors(10))
p3<-ppp + geom_point(data=pc2.dfs, aes(x=A1, y=A2, colour=PLD_Mean))+scale_colour_gradientn(colors=topo.colors(10))
p4<-ppp + geom_point(data=pc2.dfs, aes(x=A1, y=A2, colour=TrophLevel))+scale_colour_gradientn(colors=topo.colors(10))

p5<-ppp + geom_point(data=pc2.dfs, aes(x=A1, y=A2, colour=Food))
p6<-ppp + geom_point(data=pc2.dfs, aes(x=A1, y=A2, colour=Function))
p7<-ppp + geom_point(data=pc2.dfs, aes(x=A1, y=A2, colour=Aggregation))
p8<-ppp + geom_point(data=pc2.dfs, aes(x=A1, y=A2, colour=Position))
p9<-ppp + geom_point(data=pc2.dfs, aes(x=A1, y=A2, colour=SpawnMode))
p10<-ppp + geom_point(data=pc2.dfs, aes(x=A1, y=A2, colour=ParentalMode))
p11<-ppp + geom_point(data=pc2.dfs, aes(x=A1, y=A2, colour=Active))

# plot easily visualised traits together
grid.arrange(p2,p4,p6,p7,p8,p9,p10,p11, ncol=4)

# attribute PCOA species with counts per RMI region
fishdat$region=substr(fishdat$X__1, 1,2)

reg_count<-fishdat[,-1]%>%group_by(region)%>%
  summarise_all(sum, na.rm=T)

reg_count2<-data.frame(t(data.frame(reg_count)[,-1]))

names(reg_count2)<-c('Ailuk', 'Majuro', 'Rongerik', 'Rongelap')
reg_count2$Species<-names(fishdat)[2:371]
row.names(reg_count2)<-reg_count2$Species

pc2.dfs<-left_join(pc2.dfs, reg_count2, by='Species')

#subset for regional points with > 0 fish species
ailuk_df<-pc2.dfs[pc2.dfs$Ailuk>0,]
majuro_df<-pc2.dfs[pc2.dfs$Majuro>0,]
rongerik_df<-pc2.dfs[pc2.dfs$Rongerik>0,]
rongelap_df<-pc2.dfs[pc2.dfs$Rongelap>0,]

#hulls
glob_hull<-pc2.dfs[chull(pc2.dfs$A1, pc2.dfs$A2),]

ailuk_hull<-ailuk_df[chull(ailuk_df$A1, ailuk_df$A2),]
majuro_hull<-majuro_df[chull(majuro_df$A1, majuro_df$A2),]
rongerik_hull<-rongerik_df[chull(rongerik_df$A1, rongerik_df$A2),]
rongelap_hull<-rongelap_df[chull(rongelap_df$A1, rongelap_df$A2),]

p1<-ppp+
  geom_polygon(data=glob_hull,aes(x=A1,y=A2),fill=NA,colour="grey70")+
  geom_point(data=pc2.dfs, aes(x=A1, y=A2), colour='grey70')+theme_bw()+
  geom_polygon(data=ailuk_hull,aes(x=A1,y=A2),alpha=0.09, fill='red',colour="black")+
  geom_point(data=ailuk_df, aes(x=A1, y=A2), colour='red')+theme_bw()+labs(title='Ailuk')

p2<-ppp+
  geom_polygon(data=glob_hull,aes(x=A1,y=A2),fill=NA,colour="grey70")+
  geom_point(data=pc2.dfs, aes(x=A1, y=A2), colour='grey70')+theme_bw()+
  geom_polygon(data=majuro_hull,aes(x=A1,y=A2),alpha=0.09, fill='green',colour="black")+
  geom_point(data=majuro_df, aes(x=A1, y=A2), colour='green')+theme_bw()+labs(title='Majuro')


p3<-ppp+
  geom_polygon(data=glob_hull,aes(x=A1,y=A2),fill=NA,colour="grey70")+
  geom_point(data=pc2.dfs, aes(x=A1, y=A2), colour='grey70')+theme_bw()+
  geom_polygon(data=rongerik_hull,aes(x=A1,y=A2),alpha=0.09, fill='blue',colour="black")+
  geom_point(data=rongerik_df, aes(x=A1, y=A2), colour='blue')+theme_bw()+labs(title='Rongerik')


p4<-ppp+
  geom_polygon(data=glob_hull,aes(x=A1,y=A2),fill=NA,colour="grey70")+
  geom_point(data=pc2.dfs, aes(x=A1, y=A2), colour='grey70')+theme_bw()+
  geom_polygon(data=rongelap_hull,aes(x=A1,y=A2),alpha=0.09, fill='purple',colour="black")+
  geom_point(data=rongelap_df, aes(x=A1, y=A2), colour='purple')+theme_bw()+labs(title='Rongelap')

grid.arrange(p1, p2, p3, p4, ncol=2)


# OLD nMDS stuff

ord<-metaMDS(fishtrait_dist) # 

ordiplot(ord)
orditorp(ord,display="species",col="red",air=0.01)
orditorp(ord,display="sites",cex=1.25,air=0.01)
ordiellipse(ord, fish_trait.EnvTemp, col=1:4, draw='polygon', label=TRUE)

