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

# Do we want to include sp. in the analyses?
dat[grep('\\.', dat$Species),]

dat<-dat[which(dat$AUS_sp>0),] # we will focus on Australia for this prelim

# remove functional duplicates
dup_trait<-paste(dat$BodySize, dat$DepthRange, dat$PLD, dat$Diet, dat$Aggregation, dat$Position, 
                 dat$ParentalMode)
dat<-dat[-which(duplicated(dup_trait)),]
dat$Species<-as.character(dat$Species)
dat[which(dat$Species=='Scarus psittacus'),]$Species<-'Scarus psittacus/spinus' # edit for one

# Including those that default to duplicates via NA after gower dist
dat<-dat[-which(dat$Species=='Caesio sp.'),]

row.names(dat)<-dat$Species


## Edit to some trait values from MB 25/10/18
dat[dat$Species=='Brotula multibarbata',]$DepthRange<-219
dat[dat$Species=='Mobula birostris',]$BodySize<-450
dat[dat$Species=='Amphiprion sandaracinos',]$BodySize<-14


dist1<-daisy(dat[,3:9], metric='gower', stand = FALSE)

# checking for more functionally identical species by distance
dm<-as.matrix(dist1)
diag(dm)<-999
which(apply(dm, 1, function(x){min(x)==0}))
# only for Caesio sp. in AUS (removed above) !!! BUT more in other datasets


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
constree<-cl_consensus(hens, method="majority", weights=1, control=list(p=1, verbose=T))
# try with ape
#library(ape)
#hens[[7]]<-as.hclust.phylo(ape_c)
#ape_c<-consensus(lapply(hclust_methods, function(m) as.phylo(hclust(dist1, m))), p=0.5)

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
btree %>% as.dendrogram() %>% rect.dendrogram(h=0.38, 
                                              border = 1:8)

# View better in circle


bt<-btree %>% as.dendrogram()
#labels(bt)<-unlist(strsplit(labels(btree %>% as.dendrogram()), '@'))[seq(1, (nrow(dat)*2), 2)]

circlize_dendrogram(bt %>%  set("branches_k_color", h=0.38) %>%
                      set("branches_k_color", h=0.38) %>% set("labels_cex", 0.5) )

## OK number of clusters decided based on optimal number of groups that informs h value

bcut<-cutree(btree, h=0.38)

############### Run sensitivity analyses #######################
##########################################################################################


mod<-betadisper(dist1, bcut, sqrt.dist = T)
mod
plot(mod)
mod$centroids

as.matrix(mod$centroids)[1:8, 1:5]
# doesnt always make a centroid if cluster is tiny - apparently not. Seems ok.

## OK trial run

# remake distance matrix based on original scaling - only works for species resmapling
dist1_matx<-as.matrix(dist1)

sub_index<-sort(sample(1:nrow(dat), (nrow(dat)*0.95)))

distsub<-dist1_matx[sub_index, sub_index]

distsub<-as.dist(distsub)

subtree<-hclust(distsub, method='average')
subcut<-cutree(subtree, h=0.38) # tis time we use k rather than h.. NO

dend_diff(btree %>% as.dendrogram(), subtree %>% as.dendrogram())

par(mfrow=c(2,1))

btree %>% as.dendrogram() %>% set("labels", bcut) %>%
  set("branches_k_color", h=0.38, unique(bcut[btree$order])) %>% rotate(bcut)%>% plot
subtree %>% as.dendrogram() %>% set("labels", subcut) %>%
  set("branches_k_color", h=0.38, unique(subcut[subtree$order])) %>%  rotate(subcut) %>%plot

par(mfrow=c(1,1))
circlize_dendrogram(btree %>% as.dendrogram() %>% set("labels", bcut) %>%
                  set("branches_k_color", h=0.38) %>% set("labels_cex", 0.5) )

# colouring just works left to right
# wow in this example it doesnt matter about k or h.
# need to figure a solution.. work back up the original tree to assign to cultser?

mod2<-betadisper(distsub, subcut, sqrt.dist = T) 


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

maxcl_mod<-dim(as.matrix(mod$centroids))[1]
maxcl_mod2<-dim(as.matrix(mod2$centroids))[1]

mod_cent<-as.matrix(mod$centroids)[1:maxcl_mod, 1:25] # get the centroids for 8 clusters in the first 25 dimensions
dimnames(mod_cent)[[1]]<-paste('full', 1:maxcl_mod, sep='')

mod2_cent<-as.matrix(mod2$centroids)[1:maxcl_mod2, 1:25] # get the centroids for 8 clusters in the first 25 dimensions
dimnames(mod2_cent)[[1]]<-paste('subs', 1:maxcl_mod2, sep='')

# add correction for eig importance - currently apply full mod importance to both as similar
cor_mod<-matrix(data=rep(eigenvals(mod)[1:25], maxcl_mod),
                ncol=25, nrow=maxcl_mod, byrow=T)

cor_mod2<-matrix(data=rep(eigenvals(mod)[1:25], maxcl_mod2),
                 ncol=25, nrow=maxcl_mod2, byrow=T)

dist(rbind((mod_cent*cor_mod), (mod2_cent*cor_mod2)))

### Option 2, via dist object

dat_val<-dat
dsub_val<-dat[sub_index,]

distval<-daisy(rbind(dat_val[,3:9], dsub_val[,3:9]), metric='gower', stand = FALSE)
# doesnt matter that we re-create gower distance object here as scaling is set by global 'full model' parameters

vm<-as.matrix(distval)

dimnames(vm)[[1]]<-c(paste('dat', bcut, sep=''),paste('sub', subcut, sep=''))
dimnames(vm)[[2]]<-c(paste('dat', bcut, sep=''),paste('sub', subcut, sep=''))

diag(vm)<-NA

vm<-vm[(nrow(dat)+1):nrow(vm), 1:nrow(dat)] # remove dat to dat and sub to sub distances

#empty matrix
vm2<-matrix(data=NA, nrow=8, ncol=8, dimnames=list(paste('sub', 1:8, sep=''), paste('dat', 1:8, sep='')))

for(i in unique(dimnames(vm)[[2]]))
{
  for(j in unique(dimnames(vm)[[1]])) 
  {
  mz<-median(vm[which(dimnames(vm)[[1]]==j),which(dimnames(vm)[[2]]==i)], na.rm=T)
  vm2[which(dimnames(vm2)[[1]]==j),which(dimnames(vm2)[[2]]==i)]<-mz
  }
}
  
#combination of options

mod3<-betadisper(distval, c(paste('dat', bcut, sep=''),paste('sub', subcut, sep='')))

dist(as.matrix(mod3$centroids)[1:16, 1:8])

  
vm4<-as.matrix(dist(as.matrix(mod3$centroids)[1:16, 1:8]))[9:16, 1:8]

## LOOP to test analyses concept using random 5% subsample impact

# prep steps making full dataset clusters

dist1<-daisy(dat[,3:9], metric='gower', stand = FALSE)
btree<-hclust(dist1, method='average')
bcut<-cutree(btree, h=0.38)
mod<-betadisper(dist1, bcut, sqrt.dist = T)

clust_wobble_pcoa<-data.frame(full_clust=1:max(bcut))
clust_wobble_raw<-data.frame(full_clust=1:max(bcut))
my_out<-list()
for ( i in 1:100)
{
  #resample original distance matrix
  dist1_matx<-as.matrix(dist1)
  sub_index<-sort(sample(1:nrow(dat), (nrow(dat)*0.95)))
  distsub<-dist1_matx[sub_index, sub_index]
  distsub<-as.dist(distsub)
  #make subsample tree and cut 
  subtree<-hclust(distsub, method='average')
  subcut<-cutree(subtree, h=0.38) # cut using original h value
  
  mod2<-betadisper(distsub, subcut, sqrt.dist = T) 
  
  ## link clusters to originals and calc distance
  
  #option 1
  maxcl_mod<-dim(as.matrix(mod$centroids))[1]
  maxcl_mod2<-dim(as.matrix(mod2$centroids))[1]
  
  mod_cent<-as.matrix(mod$centroids)[1:maxcl_mod, 1:25] 
  dimnames(mod_cent)[[1]]<-paste('full', 1:maxcl_mod, sep='')
  
  mod2_cent<-as.matrix(mod2$centroids)[1:maxcl_mod2, 1:25] 
  dimnames(mod2_cent)[[1]]<-paste('subs', 1:maxcl_mod2, sep='')
  
  # add correction for eig importance - currently apply full mod importance to both as similar
  cor_mod<-matrix(data=rep(eigenvals(mod)[1:25], maxcl_mod),
                  ncol=25, nrow=maxcl_mod, byrow=T)
  
  cor_mod2<-matrix(data=rep(eigenvals(mod)[1:25], maxcl_mod2),
                   ncol=25, nrow=maxcl_mod2, byrow=T)
  
  pcoa_dist<-as.matrix(dist(rbind((mod_cent*cor_mod), (mod2_cent*cor_mod2))))
  pcoa_dist<-pcoa_dist[(maxcl_mod+1):nrow(pcoa_dist), 1:maxcl_mod] #
  
  pcoa_match<-apply(pcoa_dist, 1, function(x){which.min(x)})
  pcoa_match2<-apply(pcoa_dist, 2, function(x){which.min(x)})
  
  #get 1:25 pcoa dists without eig weight
  
  pcoa_dist_unweight<-as.matrix(dist(rbind((mod_cent), (mod2_cent))))
  pcoa_dist_unweight<-pcoa_dist_unweight[(maxcl_mod+1):nrow(pcoa_dist_unweight), 1:maxcl_mod]
  
  # Assumes eig weight versio is better at matching clusters than raw  1:25 pcoa axes
  #pcoa_vals<-apply(pcoa_dist_unweight, 2, function(x){min(x)})
  pcoa_vals<-pcoa_dist_unweight[pcoa_match2+c(0, maxcl_mod2*1:(maxcl_mod2-1))]
  
  clust_wobble_pcoa<-cbind(clust_wobble_pcoa, x=pcoa_match2, y=pcoa_vals)
  names(clust_wobble_pcoa)[names(clust_wobble_pcoa)=='x']<-paste('m', i, sep="")
  names(clust_wobble_pcoa)[names(clust_wobble_pcoa)=='y']<-paste('d', i, sep="")
  
  #option 2
  
  #setup combined distance matrix of full and subsampled raw data
  dat_val<-dat
  dsub_val<-dat[sub_index,]
  distval<-daisy(rbind(dat_val[,3:9], dsub_val[,3:9]), metric='gower', stand = FALSE)
  vm<-as.matrix(distval)
  
  dimnames(vm)[[1]]<-c(paste('full', bcut, sep=''),paste('subs', subcut, sep=''))
  dimnames(vm)[[2]]<-c(paste('full', bcut, sep=''),paste('subs', subcut, sep=''))
  
  diag(vm)<-NA
  
  vm<-vm[(nrow(dat)+1):nrow(vm), 1:nrow(dat)] # remove dat to dat and sub to sub distances
  
  #empty matrix
  vm2<-matrix(data=NA, nrow=maxcl_mod2, ncol=maxcl_mod,
              dimnames=list(paste('subs', 1:maxcl_mod2, sep=''), paste('full', 1:maxcl_mod, sep='')))
  
  for(k in unique(dimnames(vm)[[2]]))
  {
    for(j in unique(dimnames(vm)[[1]])) 
    {
      mz<-median(vm[which(dimnames(vm)[[1]]==j),which(dimnames(vm)[[2]]==k)], na.rm=T)
      vm2[which(dimnames(vm2)[[1]]==j),which(dimnames(vm2)[[2]]==k)]<-mz
    }
  }
  
  raw_match<-apply(vm2, 1, function(x){which.min(x)})
  raw_match2<-apply(vm2, 2, function(x){which.min(x)})
  
  raw_vals<-apply(vm2, 2, function(x){min(x)})
  
  clust_wobble_raw<-cbind(clust_wobble_raw, x=raw_match2, y=raw_vals)
  names(clust_wobble_raw)[names(clust_wobble_raw)=='x']<-paste('m', i, sep="")
  names(clust_wobble_raw)[names(clust_wobble_raw)=='y']<-paste('d', i, sep="")
  
 print(i) 
 my_out[[i]]<-list(pcoa_dist, vm2)
 
}



#cophenetic stuff
copo1<-as.matrix(cophenetic(btree))
dimnames(copo1)[[1]]<-bcut
dimnames(copo1)[[2]]<-bcut
copo1[copo1<0.38]<-0
copo1[-which(duplicated(dimnames(copo1)[[1]])), -which(duplicated(dimnames(copo1)[[2]]))]

# testing how well kmeans and mediods can recreate hieracrical clusters.. not that well
test1<-data.frame( hclust=as.character(cutree(btree, k=8)),
                   kmeans=as.character(kmeans(dist1, 8)$cluster),
                   pam=as.character(pam(dist1, 8)$cluster))

table(test1[,c(2,1)])

table(test1[,c(3,1)])