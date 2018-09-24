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

