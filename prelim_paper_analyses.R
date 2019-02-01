# 25/01/18
# Preliminary analyses for paper
# Getting method into workable script

library(cluster)
library(clue)

# read in data
dat<-read.csv('C:/coral_fish/data/Traits/JPN_AUS_RMI_trait_master.csv', h=T)

dat<-dat[which(dat$AUS_sp>0),] # we will focus on Australia for this prelim

row.names(dat)<-unlist(strsplit(as.character(dat$Species), ' '))[seq(1,(nrow(dat)-1),2)]


## Edit to some trait values from MB 25/10/18
#dat[dat$Species=='Brotula multibarbata',]$DepthRange<-219
dat[dat$Species=='Manta birostris',]$BodySize<-450
#dat[dat$Species=='Amphiprion sandaracinos',]$BodySize<-14

dist1<-daisy(dat[,3:9])

# test different algorithums and make consensus tree

hclust_methods <-
  c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty") # median and centroid dont work
hclust_results <- lapply(hclust_methods, function(m) hclust(dist1, m))
names(hclust_results) <- hclust_methods 
## Now create an ensemble from the results.
hens <- cl_ensemble(list = hclust_results)
hens

## Plotting.
plot(hens, main = names(hens))
## And continue to analyze the ensemble, e.g.
tree_agg<-cl_agreement(hens, dist1, method = "cophenetic") # this just calcs cophenetic corr 
tree_agg
# between deondrograms - but there are lots of different methods

lapply(hens, function(x){cor(cophenetic(x), dist1)})
# also between each other
cl_agreement(hens, method = "cophenetic")
# or disagreement
cl_dissimilarity(hens, method = "cophenetic")
# get consensus tree
constree<-cl_consensus(hens, 'manhattan',control=list(verbose=T))
# add to list
hens[[7]]<-constree
# not better than all
hens
# try weighting by agreement
constree2<-cl_consensus(hens[1:6], weights = as.vector(tree_agg))
hens[[8]]<-constree2#
cl_agreement(hens, dist1, method = "cophenetic") # better but not by much
plot(hens,main = names(hens))