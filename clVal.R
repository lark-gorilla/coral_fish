clVal<-function(data=data, runs=10,  max_cl=20, subs_perc=0.95, fast.k.h=0.3)
 
{
  require(cluster)
  require(fpc)
  require(dendextend)
  
  dist1<-daisy(data, metric='gower', stand = FALSE)
  btree<-hclust(dist1, method='average')
  
  wigl_list<-list()
  sens_list<-list()
  spec_list<-list()
  jacc_list<-list()
  
  out_centres<-NULL
  stats_out<-NULL
  for(o in 3:max_cl)
  {
    kval=o
    bcut<-cutree(btree, k=kval)
    
    wigl_out<-matrix(data=NA, nrow=nrow(data), ncol=runs, dimnames=list(row.names(data), 1:runs))
    clust_jacc<-matrix(data=NA, nrow=kval, ncol=runs, dimnames=list(paste('clust', 1:kval, sep=''), 1:runs))
    
    out_metrics<-data.frame(k=o, runs=1:runs, rnd=NA, jac=NA, wig=NA, sil=NA)

    for ( i in 1:runs)
    {
      #resample original distance matrix
      dist1_matx<-as.matrix(dist1)
      sub_index<-sort(sample(1:nrow(data), (nrow(data)*subs_perc)))
      distsub<-dist1_matx[sub_index, sub_index]
      distsub<-as.dist(distsub)
      #make subsample tree and cut 
      subtree<-hclust(distsub, method='average')
      subcut<-cutree(subtree, k=kval) # cut using k value
      subdata<-data[sub_index,]
      
      # Jaccard index of similarity between clusters based on species presence/absence
      # Used in fpc::clusterboot
      
      fpc_btree<-disthclustCBI(dist1,  method="average", 
                               cut="number", k=o)
      
      fpc_subtree<-disthclustCBI(distsub,  method="average", 
                                 cut="number", k=o)
      
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
      
      clust_jacc[,i]<-apply(jc2, 2, max) # gives the jaccard similarity of each original cluster to
      # the MOST similar cluster in the resampled data
      
      out_metrics[i,]$jac<-mean(apply(jc2, 2, max))
      
      # loop to create cluster centres
      
      clust_cent<-NULL
      for(h in 1:kval)
        
      {
        clust_sp<-subdata[which(subcut==h),]
        my_out<-do.call('c', lapply(as.list(clust_sp), function(x){if(is.factor(x)) {table(x)}else{median(x, na.rm=T)}}))
        out2<-data.frame(run=i, kval=o, cluster= h, as.list(my_out))
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
      
      out_metrics[i,]$rnd<- ARI
      
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
        print(paste('looking for', o, 'clusters in dendrogram went below fast cutoff search limit, set lower value for fast.k.h', sep=' '))
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
        
        wigl_out[which(dimnames(wigl_out)[[1]] %in% names(cl_cop_mv)),i]<-cl_cop_mv
        # allows for 5% species dropped to remain NA
      } 
      
      out_metrics[i,]$wig<- mean(wigl_out[,i], na.rm=T)
  
      out_metrics[i,]$sil<-cluster.stats(distsub, subcut)$avg.silwidth
      
      print(o)
      
      wigl_list[[o]]<-wigl_out
      jacc_list[[o]]<-clust_jacc
      
    }
    stats_out<-rbind(stats_out, out_metrics)
  }
  
all_out<-list(n_runs=runs, stats=stats_out, clust_centres=out_centres, wiggle=wigl_list, jaccard=jacc_list)
return(all_out)
}