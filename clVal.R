clVal<-function(data=data, runs=10,  max_cl=20, cluster_match='distance')
# cluster_match can be 1 if 'distance' or 'jaccard'  
{
  
  dist1<-daisy(data, metric='gower', stand = FALSE)
  btree<-hclust(dist1, method='average')
  
  wigl_list<-list()
  sens_list<-list()
  spec_list<-list()
  jacc_list<-list()
  
  stats_out<-NULL
  for(o in 3:max_cl)
  {
    kval=o
    bcut<-cutree(btree, k=kval)
    
    wigl_out<-matrix(data=NA, nrow=nrow(data), ncol=runs, dimnames=list(row.names(data), 1:runs))
    clust_sens<-matrix(data=NA, nrow=kval, ncol=runs, dimnames=list(paste('clust', 1:kval, sep=''), 1:runs))
    clust_spec<-matrix(data=NA, nrow=kval, ncol=runs, dimnames=list(paste('clust', 1:kval, sep=''), 1:runs))
    clust_jacc<-matrix(data=NA, nrow=kval, ncol=runs, dimnames=list(paste('clust', 1:kval, sep=''), 1:runs))
    
    
    out_kap<-NULL
    out_acc<-NULL
    out_jac<-NULL
    
    for ( i in 1:runs)
    {
      #resample original distance matrix
      dist1_matx<-as.matrix(dist1)
      sub_index<-sort(sample(1:nrow(data), (nrow(data)*0.95)))
      distsub<-dist1_matx[sub_index, sub_index]
      distsub<-as.dist(distsub)
      #make subsample tree and cut 
      subtree<-hclust(distsub, method='average')
      subcut<-cutree(subtree, k=kval) # cut using k value
      
      if( cluster_match=='distance')
          {  
          # get distances between clusters, option 2!
          #setup combined distance matrix of full and subsampled raw data
          dat_val<-data
          dsub_val<-data[sub_index,]
          distval<-daisy(rbind(dat_val, dsub_val), metric='gower', stand = FALSE)
          vm<-as.matrix(distval)
          
          dimnames(vm)[[1]]<-c(paste('full', bcut, sep=''),paste('subs', subcut, sep=''))
          dimnames(vm)[[2]]<-c(paste('full', bcut, sep=''),paste('subs', subcut, sep=''))
          
          diag(vm)<-NA
          
          vm<-vm[(nrow(data)+1):nrow(vm), 1:nrow(data)] # remove dat to dat and sub to sub distances
          
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
          
          #! running on median rather than mean !#
          
          raw_match<-apply(vm2, 1, function(x){which.min(x)})
          raw_match2<-apply(vm2, 2, function(x){which.min(x)})
      }
      
      # Jaccard index of similarity between clusters based on species presence/absence
      # Used in fpc::clusterboot
      
      fpc_btree<-disthclustCBI(dist1,  method="average", 
                               cut="number", k=o)
      
      fpc_subtree<-disthclustCBI(distsub,  method="average", 
                                 cut="number", k=o)
      
      jc2<-matrix(data=NA, nrow=kval, ncol=kval,
                  dimnames=list(paste('subs', 1:kval, sep=''), paste('full', 1:kval, sep='')))
      
      
      for(k in 1:dim(jc2)[2])
      {
        for(j in 1:dim(jc2)[1]) 
        {
          mz<-clujaccard(fpc_btree$clusterlist[[k]][sub_index], 
                         fpc_subtree$clusterlist[[j]], zerobyzero = 0)
          jc2[j,k]<-mz
        }
      }
      
      jc_match<-apply(jc2, 1, function(x){which.max(x)})
      jc_match2<-apply(jc2, 2, function(x){which.max(x)})
      
      clust_jacc[,i]<-apply(jc2, 2, max) # gives the jaccard similarity of each original cluster to
      # the MOST similar cluster in the resampled data
      
      out_jac<-c(out_jac, mean(apply(jc2, 2, max)))
      
      # make confusion matrix of clustering. COLUMNS= FULL CLUSTER, ROWS=SUBSAMPLE CLUSTERS
      sp2<-matrix(data=NA, nrow=kval, ncol=kval)
      
      for(k in 1:kval)
      {
        sp_cl<-names(subcut)[subcut==k]
        
        full_sp<-list()
        for(j in 1:kval){full_sp[[j]]<-names(bcut)[bcut==j]}
        
        sp2[k,1:kval]<-unlist(lapply(full_sp, function(x){length(which(sp_cl %in% x))}))
      }
      
      
      # re-order confusion matrix using vm2 min distances of jc2 max similarities 
      
      if(cluster_match=='distance'){mymatch=raw_match2}else{
                                    mymatch=jc_match2}
      
      sp3<-sp2[mymatch,] # reorders rows based on most likely match
      
      # give lumped clusters a 0 for row
      sp3[ which(duplicated(mymatch)),]<-0
      
      # make new col and row for split clusters
      
      sp3<-rbind(sp3, sp2[which(!1:kval %in% mymatch),])
      while(ncol(sp3)<nrow(sp3)){sp3<-cbind(sp3, 0)}
      
      # confusion matrix
      my_conf<-confusionMatrix(as.table(sp3)) 
      
      #sanity check
      #http://www.marcovanetti.com/pages/cfmatrix
      ## acc<-sum(diag(my_conf$table))/sum(my_conf$table)
      ## exp_acc<-sum(my_conf$byClass[,'Prevalence']*my_conf$byClass[,'Detection Prevalence'])
      ## kap<- (acc - exp_acc)/(1 - exp_acc)
      
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
      heights_to_cut_by<-heights_to_cut_by[heights_to_cut_by>0.3] # hack to only do larger clusters = reduces time
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
        
        cl_cop_mv<-sp2_copod[which(dimnames(sp2_copod)[[1]] %in% full_n_clust), mymatch[k]] 
        # using raw_match2 here lines up original cluster with most similar subs cluster 
        
        wigl_out[which(dimnames(wigl_out)[[1]] %in% names(cl_cop_mv)),i]<-cl_cop_mv
        # allows for 5% species dropped to remain NA
      } 
      
      print(o)
      
      wigl_list[[o]]<-wigl_out
      sens_list[[o]]<-clust_sens
      spec_list[[o]]<-clust_spec
      jacc_list[[o]]<-clust_jacc
      
    }
    stats_out<-rbind(stats_out, data.frame(k=o, acc=out_acc, kap=out_kap, jac=out_jac))
  }
  
all_out<-list(n_runs=runs, stats=stats_out, wiggle=wigl_list, sensitivity=sens_list, 
                specificity=spec_list, jaccard=jacc_list)
return(all_out)
}