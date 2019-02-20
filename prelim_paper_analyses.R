# 25/01/18
# Preliminary analyses for paper
# Getting method into workable script

library(cluster)
library(clue)
library(fpc)
library(ggplot2)
library(dendextend)
library(circlize)
library(vegan)

#################### data clean and distance matrix calculation ##########################
##########################################################################################

dat<-read.csv('C:/coral_fish/data/Traits/JPN_AUS_RMI_CHK_MLD_TMR_trait_master_opt2.csv', h=T)

dat<-dat[which(dat$AUS_sp>0),] # we will focus on Australia for this prelim

row.names(dat)<-paste(unlist(strsplit(as.character(dat$Species), ' '))[seq(1,((nrow(dat)*2)-1),2)], 1:nrow(dat), sep='@')


## Edit to some trait values from MB 25/10/18
dat[dat$Species=='Brotula multibarbata',]$DepthRange<-219
dat[dat$Species=='Manta birostris',]$BodySize<-450
dat[dat$Species=='Amphiprion sandaracinos',]$BodySize<-14

# remove duplicates
dat<-dat[-which(dat$Species == 'Orectolobus sp.' |dat$Species == 'Scarus spinus' |
                  dat$Species == 'Variola sp.' ),]

dist1<-daisy(dat[,3:9], metric='gower', stand = FALSE)

############### test different algorithums and make consensus tree #######################
##########################################################################################

hclust_methods <-
  c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty") # median and centroid dont work
hclust_results <- lapply(hclust_methods, function(m) hclust(dist1, m))
names(hclust_results) <- hclust_methods 

## Now create an ensemble from the results.
hens <- cl_ensemble(list = hclust_results)

## Plotting.
#plot(hens, main = names(hens))

## see how well the trees represent the original distance matrix.
## method=spectral is the 2-norm method proposed by Merigot et al. (2010)
tree_agg<-cl_dissimilarity(hens, dist1, method = "spectral")
tree_agg

## get threshold 2-norm value for dendrgram to accurately represent dist. 
## code from within Mouchet et al. 2008 function: ("http://villeger.sebastien.free.fr/R%20scripts/GFD_matcomm.R")

# best tree as hclust object
btree<-hclust(dist1, method='average')
sigma2 <- var(dist1) + var(cl_ultrametric(btree))
threshold <- 2*sqrt(sigma2)*sqrt(dim(dat)[1])
threshold # still too high

# This is the consensus tree calc from Mouchet - I wouldnt get to work
# cl_consensus(hens, method="majority", weights=1, control=list(p=1))

# get consensus tree
constree<-cl_consensus(hens[5:6], 'manhattan',control=list(verbose=T))
# add to list
hens[[7]]<-constree
# try weighting by agreement
constree2<-cl_consensus(hens[5:6], 'manhattan', weights = as.vector(tree_agg)[5:6])
hens[[8]]<-constree2
cl_dissimilarity(hens, dist1, method = "spectral")
# both worse

## we take average algorithum tree forward as the best
btree<-hclust(dist1, method='average')

btree %>% as.dendrogram() %>% plot

############### Work out optimal n clusters and define h value #######################
##########################################################################################

clust_out<-NULL
for(i in 2:20){
  mnz<-cluster.stats(dist1, cutree(btree, k=i))
  out<-data.frame(nclust=i, av.sil=mnz$avg.silwidth, wb.ratio=mnz$wb.ratio) 
  clust_out<-rbind(clust_out, out)}

par(mfrow=c(1,2))
plot(av.sil~nclust, data=clust_out, type='b')
plot(wb.ratio~nclust, data=clust_out, type='b')

#potential other approach
lp1<-lapply(2:20, function(m) cutree(btree, m))
dat.t<-cbind(dat, as.data.frame(do.call('cbind', lp1)))
library(vegan)
v_runs<-lapply(names(dat.t[13:31]), function(m) adonis(as.formula(paste('dist1~', m, sep="")), data=dat.t))
# see also: https://github.com/pmartinezarbizu/pairwiseAdonis
# not sure is this is valid

cbind(seq(0.25, 0.45, 0.01), unlist(lapply(seq(0.25, 0.45, 0.01),
                                           function(m) length(unique(cutree(btree, h=m))))))
par(mfrow=c(1,1))
btree %>% as.dendrogram() %>% plot
btree %>% as.dendrogram() %>% rect.dendrogram(h=0.37, 
                                              border = 1:8)

# View better in circle


bt<-btree %>% as.dendrogram()
#labels(bt)<-unlist(strsplit(labels(btree %>% as.dendrogram()), '@'))[seq(1, (nrow(dat)*2), 2)]

circlize_dendrogram(bt %>%  set("branches_k_color", h=0.37) %>%
                      set("branches_k_color", h=0.37) %>% set("labels_cex", 0.5) )

## OK number of clusters decided based on optimal number of groups that informs h value

bcut<-cutree(btree, h=0.37)

############### Run sensitivity analyses #######################
##########################################################################################


mod<-betadisper(dist1, bcut)
mod
plot(mod)
mod$centroids

as.matrix(mod$centroids)[1:8, 1:5]
# doesnt always make a centroid if cluster is tiny - apparently not. Seems ok.

## OK trial run

dsub<-dat[sort(sample(1:nrow(dat), (nrow(dat)*0.95))),]

distsub<-daisy(dsub[,3:9], metric='gower', stand = FALSE)

subtree<-hclust(distsub, method='average')

dend_diff(btree %>% as.dendrogram(), subtree %>% as.dendrogram())

par(mfrow=c(2,1))
btree %>% as.dendrogram()  %>% set("branches_k_color", h=0.37, col=bcut)  %>% plot
subtree %>% as.dendrogram()  %>% set("branches_k_color", h=0.37)  %>% plot
# colouring just works left to right
# wow in this example it doesnt matter about k or h.
# need to figure a solution.. work back up the original tree to assign to cultser?

mod2<-betadisper(distsub, cutree(subtree, k=8)) # tis time we use k rather than h..


# plot betadisper outputs alongside
par(mfrow=c(1,2))
plot(mod); plot(mod2)

dist(rbind(mod$centroids[1,], mod2$centroids[1,]))

sum(eigenvals(mod)[1:25])/sum(eigenvals(mod)[which(eigenvals(mod)>0)])
sum(eigenvals(mod)[1:4])/sum(eigenvals(mod)[which(eigenvals(mod)>0)])
#first 25 PCoA axes explain 99% of variation
# check if need to weight by var expl, or just select n eigenvals to elbow
# OK so we can select first 25 eigenvals then weight by sqrt of eigen value
# to diminish impact of less important axes following: https://www.pnas.org/content/110/38/15307#ref-57
# But this would mean mod and mod2 have differnt weights - not sure if this is a good idea e.g.

sqrt(eigenvals(mod)[1:25])
sqrt(eigenvals(mod2)[1:25])

# best ideas 
# 1) weight both mod and mod2 by mod's sqrt(eigenvals) - this assumes that they are basically the same,
# this holds for 95% subsample but maybe not if big diffs to gower distb matrix e.g. -1 trait variable
# 2) Take fewer eigenvals 1:10 gives ~ 87% var expl, or even 1:4 is 60%. tradeoff between info and crap

# test with ability to assign correct clusters


as.matrix(mod$centroids)[1:8, 1:25]

mod_cent<-as.matrix(mod$centroids)[1:8, 1:25] # get the centroids for 8 clusters in the first 25 dimensions
dimnames(mod_cent)[[1]]<-paste('full', 1:8, sep='')

mod2_cent<-as.matrix(mod2$centroids)[1:8, 1:25] # get the centroids for 8 clusters in the first 25 dimensions
dimnames(mod2_cent)[[1]]<-paste('subs', 1:8, sep='')

dist(rbind(mod_cent, mod2_cent))




