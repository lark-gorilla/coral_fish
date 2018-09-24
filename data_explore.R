# 20/09/2018

# Initial explore of fish and coral survey and trait data

library(dplyr)
library(readxl)

#https://daijiang.name/en/2014/05/11/functional-diversity-in-r/

library(FD)

setwd("~/leeds_postdoc")

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

fishtrait_dist<-gowdis(fishtrait2)
  
library(ade4)

fishtrait_dist2 <- cailliez(fishtrait_dist)

pcoa<-dudi.pco(fishtrait_dist2)

plot(pcoa$li[,1], pcoa$li[,2], type = "n", xlab = "", ylab = "",
     axes = FALSE, main = "dudi.pco (ade4)")
text(pcoa$li[,1], pcoa$li[,2], labels(fishtrait_dist2), cex = 0.9)

fish_clust<-hclust(fishtrait_dist)

#https://www.r-bloggers.com/7-functions-to-do-metric-multidimensional-scaling-in-r/
