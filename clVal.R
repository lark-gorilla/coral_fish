clVal<-function(data=data, runs=10, min_cl=3,  max_cl=20, subs_perc=0.95, fast.k.h=0.3, calc_wigl=T)
 
{
  require(cluster)
  require(fpc)
  require(dendextend)
  
  dist1<-daisy(data, metric='gower', stand = FALSE)
  btree<-hclust(dist1, method='average')
  if(is.null(btree$labels)){btree$labels<-1:nrow(data)} # edit to make sure tree has labels
  
  
  #create output objects
  wigl_list<-as.list(rep(list(matrix(data=NA, nrow=nrow(data), ncol=runs,
                                     dimnames=list(row.names(data), 1:runs))), max_cl))
  jacc_list<-list()
  for(i in 1:max_cl)
  {jacc_list[[i]]<-matrix(data=NA, nrow=i, ncol=runs, 
                   dimnames=list(paste('clust', 1:i, sep=''), 1:runs))}
  
  jacc_out<-expand.grid(k=min_cl:max_cl, runs=1:runs)
  out_centres<-NULL
  stats_out<-expand.grid(k=min_cl:max_cl, runs=1:runs, rnd=NA, jac=NA, wig=NA, sil=NA)
  


for ( i in 1:runs)
{
  #resample original distance matrix
  dist1_matx<-as.matrix(dist1)
  sub_index<-sort(sample(1:nrow(data), (nrow(data)*subs_perc)))
  distsub<-dist1_matx[sub_index, sub_index]
  distsub<-as.dist(distsub)
  
      for(kval in min_cl:max_cl)
      {
      bcut<-cutree(btree, k=kval)
        
      #make subsample tree and cut 
      subtree<-hclust(distsub, method='average')
      subcut<-cutree(subtree, k=kval) # cut using k value
      subdata<-data[sub_index,]
      
      # Jaccard index of similarity between clusters based on species presence/absence
      # Used in fpc::clusterboot
      
      fpc_btree<-disthclustCBI(dist1,  method="average", 
                               cut="number", k=kval)
      
      fpc_subtree<-disthclustCBI(distsub,  method="average", 
                                 cut="number", k=kval)
      
      #empty matrix
      jc2<-matrix(data=NA, nrow=kval, ncol=kval,
                  dimnames=list(paste('subs', 1:kval, sep=''), paste('full', 1:kval, sep='')))
     
      
      for(k in 1:dim(jc2)[2])
      {
        for(j in 1:dim(jc2)[1]) 
        {
          
          jc2[j,k]<-clujaccard(fpc_btree$clusterlist[[k]][sub_index], 
                         fpc_subtree$clusterlist[[j]], zerobyzero = 0)
        }
      }
  
      jc_match<-apply(jc2, 1, function(x){which.max(x)}) 
      jc_match2<-apply(jc2, 2, function(x){which.max(x)})
      
      jacc_list[[kval]][,i]<-apply(jc2, 2, max) # gives the jaccard similarity of each original cluster to
      # the MOST similar cluster in the resampled data
      
      stats_out[stats_out$k==kval & stats_out$runs==i,]$jac<-mean(apply(jc2, 2, max))
      
      # loop to create cluster centres
      
      clust_cent<-NULL
      for(h in 1:kval)
        
      {
        clust_sp<-subdata[which(subcut==h),]
        my_out<-do.call('c', lapply(as.list(clust_sp), function(x){if(is.factor(x)) {table(x)}else{median(x, na.rm=T)}}))
        out2<-data.frame(run=i, kval=kval, cluster= h, as.list(my_out))
        clust_cent<-rbind(clust_cent, out2)
      }  
      clust_cent$jc_match<-jc_match
      out_centres<-rbind(out_centres, clust_cent)
      

      ### adjusted rand index, code hacked from mclust::adjustedRandIndex
      
      tab <- table(subcut, bcut[which(names(bcut) %in% names(subcut))])
      if (all(dim(tab) == c(1, 1))) 
        return(1)
      a <- sum(choose(tab, 2))
      b <- sum(choose(rowSums(tab), 2)) - a
      c <- sum(choose(colSums(tab), 2)) - a
      d <- choose(sum(tab), 2) - a - b - c
      ARI <- (a - (a + b) * (a + c)/(a + b + c + d))/((a + b + 
                   a + c)/2 - (a + b) * (a + c)/(a + b + c + d))
      
      stats_out[stats_out$k==kval & stats_out$runs==i,]$rnd<- ARI
      
      if(calc_wigl==T)
          {
          # get h val from k val
          # hacked function from dendextend function heights_per_k.dendrogram()
          
          dend=as.dendrogram(subtree)
          our_dend_heights <- sort(unique(get_branches_heights(dend, 
                                                               sort = FALSE)), TRUE)
          heights_to_remove_for_A_cut <- min(-diff(our_dend_heights))/2
          heights_to_cut_by <- c((max(our_dend_heights) + heights_to_remove_for_A_cut), 
                                 (our_dend_heights - heights_to_remove_for_A_cut))
          heights_to_cut_by<-heights_to_cut_by[heights_to_cut_by>fast.k.h] # hack to only do larger clusters = reduces time
          names(heights_to_cut_by) <- sapply(heights_to_cut_by, function(h) {
            length(cut(dend, h = h)$lower)})
          
          cutval<-heights_to_cut_by[names(heights_to_cut_by)==kval]
          
          if(length(cutval)==0){
            print(paste('looking for', kval, 'clusters in dendrogram went below fast cutoff search limit, set lower value for fast.k.h', sep=' '))
            break}
          
          # copo distance per species
          copo2<-as.matrix(cophenetic(subtree))
          dimnames(copo2)[[1]]<-names(subcut)
          dimnames(copo2)[[2]]<-subcut
          copo2[copo2<= cutval]<-0 # below AND equal
          
          sp2_copod<-copo2[, -which(duplicated(dimnames(copo2)[[2]]))]
          
          for(k in 1:kval)
              {
                full_n_clust<-names(bcut)[bcut==k]
                
                cl_cop_mv<-sp2_copod[which(dimnames(sp2_copod)[[1]] %in% full_n_clust), jc_match2[k]] 
                # using jc2_match2 here lines up original cluster with most similar subs cluster 
                
                wigl_list[[kval]][which(dimnames(wigl_list[[kval]])[[1]] %in% names(cl_cop_mv)),i]<-cl_cop_mv
                # allows for 5% species dropped to remain NA
              } 
          
          stats_out[stats_out$k==kval & stats_out$runs==i,]$wig<- mean(wigl_list[[kval]][,i], na.rm=T)
          }
      
      stats_out[stats_out$k==kval & stats_out$runs==i,]$sil<-cluster.stats(distsub, subcut)$avg.silwidth
      
      }#close kval loop
print(i)
print(Sys.time())
}# close runs loop
  
all_out<-list(n_runs=runs, stats=stats_out, clust_centres=out_centres, wiggle=wigl_list, jaccard=jacc_list)
return(all_out)
}