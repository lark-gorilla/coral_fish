# 25/01/18
# Preliminary analyses for paper
# Getting method into workable script
rm(list=ls())
library(cluster)
library(clue)
library(fpc)
library(ggplot2)
library(dendextend)
library(circlize)
library(vegan)
library(caret)



#################### data clean and distance matrix calculation ##########################
##########################################################################################

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


dist1<-daisy(dat[,3:9], metric='gower', stand = FALSE)

# checking for more functionally identical species by distance
dm<-as.matrix(dist1)
diag(dm)<-999
which(apply(dm, 1, function(x){min(x)==0}))
# only for 3 sp. in AUS+JPN (removed above) !!! BUT more in other datasets?


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

cbind(seq(0.25, 0.5, 0.01), unlist(lapply(seq(0.25, 0.5, 0.01),
                                           function(m) length(unique(cutree(btree, h=m))))))
par(mfrow=c(1,1))
btree %>% as.dendrogram() %>% plot
btree %>% as.dendrogram() %>% rect.dendrogram(h=0.45, 
                                              border = 1:6)

# View better in circle


bt<-btree %>% as.dendrogram()
#labels(bt)<-unlist(strsplit(labels(btree %>% as.dendrogram()), '@'))[seq(1, (nrow(dat)*2), 2)]

circlize_dendrogram(bt %>%  set("branches_k_color", h=0.38) %>%
                      set("branches_k_color", h=0.38) %>% set("labels_cex", 0.5) )

## OK number of clusters decided based on optimal number of groups that informs h value
cutval=0.46
bcut<-cutree(btree, h=cutval)

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
subcut<-cutree(subtree, h=cutval) # tis time we use k rather than h.. NO

dend_diff(btree %>% as.dendrogram(), subtree %>% as.dendrogram())

par(mfrow=c(2,1))

btree %>% as.dendrogram() %>% set("labels", bcut) %>%
  set("branches_k_color", h=cutval, unique(bcut[btree$order])) %>% rotate(bcut)%>% plot
subtree %>% as.dendrogram() %>% set("labels", subcut) %>%
  set("branches_k_color", h=cutval, unique(subcut[subtree$order])) %>%  rotate(subcut) %>%plot

par(mfrow=c(1,1))
circlize_dendrogram(btree %>% as.dendrogram() %>% set("labels", bcut) %>%
                  set("branches_k_color", h=cutval) %>% set("labels_cex", 0.5) )

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

# Linking full and subsample clusters using species identity

sp2<-matrix(data=NA, nrow=maxcl_mod2, ncol=maxcl_mod,
            dimnames=list(paste('subs', 1:maxcl_mod2, sep=''), paste('full', 1:maxcl_mod, sep='')))

for(k in 1:maxcl_mod2)
{
  sp_cl<-names(subcut)[subcut==k]
  
  full_sp<-list()
  for(j in 1:maxcl_mod){full_sp[[j]]<-names(bcut)[bcut==j]}
  
  sp2[k,]<-unlist(lapply(full_sp, function(x){ceiling(length(which(sp_cl %in% x)))}))
}



## LOOP to test analyses concept using random 5% subsample impact

# prep steps making full dataset clusters

dist1<-daisy(dat[,3:9], metric='gower', stand = FALSE)
btree<-hclust(dist1, method='average')
cutval=0.46
bcut<-cutree(btree, h=cutval)
mod<-betadisper(dist1, bcut, sqrt.dist = T)


out_pcoD<-list()
out_pcoC<-list()
out_rawD<-list()
out_rawC<-list()
for ( i in 1:100)
{
  #resample original distance matrix
  dist1_matx<-as.matrix(dist1)
  sub_index<-sort(sample(1:nrow(dat), (nrow(dat)*0.95)))
  distsub<-dist1_matx[sub_index, sub_index]
  distsub<-as.dist(distsub)
  #make subsample tree and cut 
  subtree<-hclust(distsub, method='average')
  subcut<-cutree(subtree, h=cutval) # cut using original h value
  
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
  
  # fix for clusters in full dataset not found in subset via pcoa_match but found via pcoa_match2
  if(maxcl_mod2>maxcl_mod) # if equal to then can replace correct assinments..
    {
    ms_clust<-which(!1:maxcl_mod %in% pcoa_match)
    print(paste('missing full clust', ms_clust))
    print(pcoa_match)
    print(pcoa_match2)
  
    pcoa_match[pcoa_match2[which(!1:maxcl_mod %in% pcoa_match)]]<-ms_clust
    print('edited to..');print(pcoa_match)
    }
    
  #get 1:25 pcoa dists without eig weight
  
  pcoa_dist_unweight<-as.matrix(dist(rbind((mod_cent), (mod2_cent))))
  pcoa_dist_unweight<-pcoa_dist_unweight[(maxcl_mod+1):nrow(pcoa_dist_unweight), 1:maxcl_mod]
  
  # Assumes eig weight versio is better at matching clusters than raw  1:25 pcoa axes
  #pcoa_vals<-apply(pcoa_dist_unweight, 2, function(x){min(x)}) get distance per subset cluster (eig weighted data)
  pcoa_vals<-pcoa_dist_unweight[cbind(1:maxcl_mod2,pcoa_match)] #get distance per full cluster (unweighted data)
  #pcoa_vals<-pcoa_dist_unweight[pcoa_match2+c(0, maxcl_mod2*1:(maxcl_mod2-1))] get distance per full cluster (unweighted data)
  
  lz<-as.list(1:maxcl_mod)
  for(m in lz)
  {
    lz[[m]]<-pcoa_vals[which(pcoa_match==m)]
  }
    
  
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
  
  # fix for clusters in full dataset not found in subset via raw_match but found via raw_match2
  if(maxcl_mod2>maxcl_mod)
    {
    ms_clust<-which(!1:maxcl_mod %in% raw_match)
    print(paste('missing full clust RAW', ms_clust))
    print(raw_match)
    print(raw_match2)
  
    raw_match[raw_match2[which(!1:maxcl_mod %in% raw_match)]]<-ms_clust
    print('RAW edited to..');print(raw_match)
    }
  
  raw_vals<-apply(vm2, 1, function(x){min(x)})
  
  lz2<-as.list(1:maxcl_mod)
  for(m in lz2)
  {
    lz2[[m]]<-raw_vals[which(raw_match==m)]
  }
  
  
 print(i) 
 out_pcoD[[i]]<-pcoa_dist
 out_pcoC[[i]]<-lz
 out_rawD[[i]]<-vm2
 out_rawC[[i]]<-lz2
 
 
 # debug mode
 print(pcoa_dist)
 print(pcoa_match)
 print(pcoa_match2)
 print(vm2)
 print(raw_match)
 print(raw_match2)
 
 sp2<-matrix(data=NA, nrow=maxcl_mod2, ncol=maxcl_mod,
             dimnames=list(paste('subs', 1:maxcl_mod2, sep=''), paste('full', 1:maxcl_mod, sep='')))
 
 for(k in 1:maxcl_mod2)
 {
   sp_cl<-names(subcut)[subcut==k]
   
   full_sp<-list()
   for(j in 1:maxcl_mod){full_sp[[j]]<-names(bcut)[bcut==j]}
   
   sp2[k,]<-unlist(lapply(full_sp, function(x){ceiling(length(which(sp_cl %in% x)))}))
 }
 
 print(sp2)
 
 par(mfrow=c(2,1))
 
 btree %>% as.dendrogram() %>% set("labels", bcut) %>%
   set("branches_k_color", h=cutval, unique(bcut[btree$order])) %>% rotate(bcut)%>% plot
 subtree %>% as.dendrogram() %>% set("labels", subcut) %>%
   set("branches_k_color", h=cutval, unique(subcut[subtree$order])) %>%  rotate(subcut) %>%plot
 
 readline('')
 
 
}

do.call(cbind, out_pcoC)

#cophenetic stuff
copo1<-as.matrix(cophenetic(btree))
dimnames(copo1)[[1]]<-bcut
dimnames(copo1)[[2]]<-bcut
copo1[copo1<= cutval]<-0 # below AND equal

#copo distance between clusters
copo1[-which(duplicated(dimnames(copo1)[[1]])), -which(duplicated(dimnames(copo1)[[2]]))]

copo1<-as.matrix(cophenetic(btree))
dimnames(copo1)[[1]]<-names(bcut)
dimnames(copo1)[[2]]<-bcut
copo1[copo1<= cutval]<-0 # below AND equal

#copo distance of indivdiual species between clusters, same as above but species level 
copo1[1:10, -which(duplicated(dimnames(copo1)[[2]]))]

## for resampled matrix

copo2<-as.matrix(cophenetic(subtree))
dimnames(copo2)[[1]]<-names(subcut)
dimnames(copo2)[[2]]<-subcut
copo2[copo2<= cutval]<-0 # below AND equal

sp2_copod<-copo2[, -which(duplicated(dimnames(copo2)[[2]]))]

# proof of concept for cluster 1
full_n_clust<-bcut[sub_index] # only take species present in subsample for comparison
full_n_clust<-names(full_n_clust)[full_n_clust==1] 
sp2_copod[which(dimnames(sp2_copod)[[1]] %in% full_n_clust),]

# loop for all clusters
sp3<-matrix(data=NA, nrow=maxcl_mod2, ncol=maxcl_mod,
            dimnames=list(paste('subs', 1:maxcl_mod2, sep=''), paste('full', 1:maxcl_mod, sep='')))

for(k in 1:maxcl_mod)
{
  full_n_clust<-names(bcut)[bcut==k]
  
  cl_cop_mv<-sp2_copod[which(dimnames(sp2_copod)[[1]] %in% full_n_clust), k]
  
 print(k)
 print(sp2_copod[which(dimnames(sp2_copod)[[1]] %in% full_n_clust),])
 print(mean(cl_cop_mv))
}
 
# loop to calc copo movement value for each species 'wiggliness'

wigl_out<-data.frame(Species=dat$Species)

for(k in 1:maxcl_mod)
{
  full_n_clust<-names(bcut)[bcut==k]
  
  cl_cop_mv<-sp2_copod[which(dimnames(sp2_copod)[[1]] %in% full_n_clust), k]
  
  wigl_out<-cbind(wigl_out, runni=NA)
  names(wigl_out)[names(wigl_out)=='runni']<-paste('run', i, sep='')
  
  wigl_out[which(wigl_out$Species %in% names(cl_cop_mv)),]$run1<-cl_cop_mv

} # doesnt' work but just concept




# accuracy for sp2 table
library(caret)

sp_tab<-sp2
sp_tab<-cbind(sp_tab, rep(0, 7)) # make dummy for extra cluster in subs
dimnames(sp_tab)[[2]][7]<-'full7'
row.names(sp_tab)<-colnames(sp_tab)

confusionMatrix(as.table(sp_tab))

# Loop to run species wiggliness and accuracy of clustering classification  using random 5% subsample impact

# prep steps making full dataset clusters

dist1<-daisy(dat[,3:9], metric='gower', stand = FALSE)
btree<-hclust(dist1, method='average')
runz=50

wigl_list<-list()
sens_list<-list()
spec_list<-list()
stats_out<-NULL
for(o in 3:20)
{
kval=o
bcut<-cutree(btree, k=kval)

wigl_out<-matrix(data=NA, nrow=nrow(dat), ncol=runz, dimnames=list(dat$Species, 1:runz))
clust_sens<-matrix(data=NA, nrow=kval, ncol=runz, dimnames=list(paste('clust', 1:kval, sep=''), 1:runz))
clust_spec<-matrix(data=NA, nrow=kval, ncol=runz, dimnames=list(paste('clust', 1:kval, sep=''), 1:runz))

out_kap<-NULL
out_acc<-NULL
for ( i in 1:runz)
{
  #resample original distance matrix
  dist1_matx<-as.matrix(dist1)
  sub_index<-sort(sample(1:nrow(dat), (nrow(dat)*0.95)))
  distsub<-dist1_matx[sub_index, sub_index]
  distsub<-as.dist(distsub)
  #make subsample tree and cut 
  subtree<-hclust(distsub, method='average')
  subcut<-cutree(subtree, k=kval) # cut using k value
  
  # get distances between clusters, option 2!
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
  vm2<-matrix(data=NA, nrow=kval, ncol=kval,
              dimnames=list(paste('subs', 1:kval, sep=''), paste('full', 1:kval, sep='')))
  
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
  
  
  # make confusion matrix of clustering. COLUMNS= FULL CLUSTER, ROWS=SUBSAMPLE CLUSTERS
  sp2<-matrix(data=NA, nrow=kval, ncol=kval)
  
  for(k in 1:kval)
  {
    sp_cl<-names(subcut)[subcut==k]
    
    full_sp<-list()
    for(j in 1:kval){full_sp[[j]]<-names(bcut)[bcut==j]}
    
    sp2[k,1:kval]<-unlist(lapply(full_sp, function(x){length(which(sp_cl %in% x))}))
  }
  
  print(sp2)
  
  # re-order confusion matrix using vm2 min distances
  
  sp3<-sp2[raw_match2,] # reorders rows based on most likely match
  
  # give lumped clusters a 0 for row
  sp3[ which(duplicated(raw_match2)),]<-0
  
  # make new col and row for split clusters
  
  sp3<-rbind(sp3, sp2[which(!1:kval %in% raw_match2),])
  while(ncol(sp3)<nrow(sp3)){sp3<-cbind(sp3, 0)}

  # confusion matrix
  my_conf<-confusionMatrix(as.table(sp3)) 
  
  out_acc<-c(out_acc, my_conf$overall[1])
  out_kap<-c(out_kap, my_conf$overall[2])
  # although we have inflated the confusion matrix with splitters we do not output
  # these as only interested in changes to original clusters. The edits to the confusion
  # matrix via lumpers and splitters change the sens/spec of each original cluster - which is what we want
  
  clust_sens[,i]<-my_conf$byClass[1:kval,1] # this tells how accurately species in original
  # cluster are still in subs cluster (identical/splitting) 
  
  clust_spec[,i]<-my_conf$byClass[1:kval,2] # this tells how accurately species in sub
  # cluster are from the same original cluster (identical/lumped) 
  
  # get h val from k val
  # hacked function from dendextend function heights_per_k.dendrogram()
  
  dend=as.dendrogram(subtree)
  our_dend_heights <- sort(unique(get_branches_heights(dend, 
                                                       sort = FALSE)), TRUE)
  heights_to_remove_for_A_cut <- min(-diff(our_dend_heights))/2
  heights_to_cut_by <- c((max(our_dend_heights) + heights_to_remove_for_A_cut), 
                         (our_dend_heights - heights_to_remove_for_A_cut))
  heights_to_cut_by<-heights_to_cut_by[heights_to_cut_by>0.3] # hack to only do larger clusters
  names(heights_to_cut_by) <- sapply(heights_to_cut_by, function(h) {
    length(cut(dend, h = h)$lower)})
  
  cutval<-heights_to_cut_by[names(heights_to_cut_by)==kval]
  
  # copo distance per species
  copo2<-as.matrix(cophenetic(subtree))
  dimnames(copo2)[[1]]<-names(subcut)
  dimnames(copo2)[[2]]<-subcut
  copo2[copo2<= cutval]<-0 # below AND equal
  
  sp2_copod<-copo2[, -which(duplicated(dimnames(copo2)[[2]]))]
  
  for(k in 1:kval)
  {
    full_n_clust<-names(bcut)[bcut==k]
    
    cl_cop_mv<-sp2_copod[which(dimnames(sp2_copod)[[1]] %in% full_n_clust), raw_match2[k]] 
    # using raw_match2 here lines up original cluster with most similar subs cluster 
    
    wigl_out[which(dimnames(wigl_out)[[1]] %in% names(cl_cop_mv)),i]<-cl_cop_mv
    # allows for 5% species dropped to remain NA
    } 
  
  print(o)
  
  wigl_list[[o]]<-wigl_out
  sens_list[[o]]<-clust_sens
  spec_list[[o]]<-clust_spec
  stats_out<-rbind(stats_out, data.frame(k=o, acc=out_acc, kap=out_kap))
  
}
}

# summary stats

sort(rowMeans(wigl_out, na.rm=T))

df_mean<-data.frame(Species=row.names(wigl_out), av=rowMeans(wigl_out, na.rm=T))
df_mean<-df_mean[order(df_mean$av),]

df_mean$Species<-as.character(df_mean$Species)
df_mean$Species <- factor(df_mean$Species, levels=order(df_mean$av))


library(reshape2)
dfr<-data.frame(wigl_out[order(df_mean$av),])
dfr$Species<-row.names(wigl_out)

dfr<-melt(dfr, id.vars = 'Species')

dfr$Species<-as.character(dfr$Species)
dfr$Species <- factor(dfr$Species, levels=df_mean[order(df_mean$av, order=-1),]$Species)

df_mean[which(df_mean$av==0),]$Species

mv_sp<-df_mean[which(df_mean$av>0),]$Species

ggplot(data=dfr[dfr$Species%in% mv_sp,], aes(x=Species, y=value))+
  geom_jitter(height=0.01, width=0.0001, shape=1, alpha=0.05)+
  geom_point(data=df_mean[df_mean$Species%in% mv_sp,], aes(x=Species, y=av), colour='red')+
  theme_bw()+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))

#acc
qplot(data=stats_out, x=k, y=acc)+geom_smooth()
ggplot(data=stats_out, aes(x=k, y=acc, group=k))+geom_boxplot()+scale_x_continuous(breaks=3:20)

ggplot(data=stats_out, aes(x=k, y=kap, group=k))+geom_boxplot()+scale_x_continuous(breaks=3:20)+
  geom_line(data=clust_out, aes(x=nclust, y=av.sil), colour=2)

clust_out<-NULL
for(i in 3:20){
  mnz<-cluster.stats(dist1, cutree(btree, k=i))
  out<-data.frame(nclust=i, av.sil=mnz$avg.silwidth, wb.ratio=mnz$wb.ratio) 
  clust_out<-rbind(clust_out, out)}

par(mfrow=c(1,2))
plot(av.sil~nclust, data=clust_out, type='b')
