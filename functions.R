## 26/04/19
## Functions to assist analyses

## function to do pca, get axis 1:3, calc kernel density for UD 50 & 99 and make plots

pca_vis<-function(rundat=dat, clValresult=clv, kval=7)
{  
  require(ade4)
  require(ggplot2)
  require(adehabitatHR)
  require(sf) 
  require(cluster)
  
  dist1<-daisy(rundat, metric='gower', stand = FALSE)
  btree<-hclust(dist1, method='average')
  bcut<-cutree(btree, k=kval)
  
  clust_cent<-NULL
  for(h in 1:kval)
  {
    clust_sp<-rundat[which(bcut==h),]
    my_out<-do.call(c, lapply(as.list(clust_sp), function(x){if(is.factor(x)) {table(x)}else{median(x, na.rm=T)}}))
    out2<-data.frame(run=0, kval=kval, cluster= h, as.list(my_out), jc_match=h)
    clust_cent<-rbind(clust_cent, out2)
  }  
  
  # get cluster centres for optimal clusters from subsampling then bind full cluster centres
  aus7<-clValresult[clValresult$kval==kval,]
  
  aus7<-rbind(aus7, clust_cent)
  
  ## Run PCA
  
  # setup weights for each column, factors penalised for n levels
  my_wgt<-c(1,1,1,rep(1/7, 7), rep(1/4, 4), rep(1/6, 6), rep(1/5, 5))
  
  aus7_pca<-dudi.pca(aus7[,4:28], col.w=my_wgt, center=T, scale=T, scannf = FALSE, nf = 3)
  
  aus7<-cbind(aus7, aus7_pca$li)
  
  ## plots and kernel areas for 1:3 PCAs
  
  ### PCA 1- PCA2
  
  #kernel calc
  
  spdf<-SpatialPointsDataFrame(coords=cbind(aus7$Axis1, aus7$Axis2),
                               data=data.frame(jc_match=factor(aus7$jc_match)))
  
  KDE.Surface <- kernelUD(spdf,same4all = T, h=0.05, grid=1000)
  
  k1_2<-kernel.area(KDE.Surface, percent = c(50, 99))
  
  KDE.UD1_2 <- getverticeshr(KDE.Surface, percent = 99)
  
  ## plotting
  
  
  g1_2<- ggplot()+
    geom_hline(yintercept=0, linetype="dotted") + 
    geom_vline(xintercept=0,  linetype="dotted")+
    geom_sf(data=st_as_sf(KDE.UD1_2), colour=rainbow(kval), alpha=0.5)+
    geom_point(data=aus7[1:(nrow(aus7)-kval),], aes(x=Axis1, y=Axis2, colour=factor(jc_match)), shape=3)+
    geom_point(data=aus7[(nrow(aus7)-(kval-1)):nrow(aus7),],
               aes(x=Axis1, y=Axis2, fill=factor(jc_match)), shape=21, colour='black', size=3)+
    theme_bw()+theme(legend.position = "none")+
    scale_color_manual(values=rainbow(kval))
  
  eig<-aus7_pca$eig
  g1_2<- g1_2+scale_y_continuous(paste('PC2', sprintf('(%0.1f%% explained var.)', 100* eig[2]/sum(eig))))+
    scale_x_continuous(paste('PC1', sprintf('(%0.1f%% explained var.)', 100* eig[1]/sum(eig))))
  
  
  ### PCA 1- PCA3
  
  #kernel calc
  
  spdf<-SpatialPointsDataFrame(coords=cbind(aus7$Axis1, aus7$Axis3),
                               data=data.frame(jc_match=factor(aus7$jc_match)))
  
  KDE.Surface <- kernelUD(spdf,same4all = T, h=0.05, grid=1000)
  
  k1_3<-kernel.area(KDE.Surface, percent = c(50, 99))
  
  KDE.UD1_3 <- getverticeshr(KDE.Surface, percent = 99)
  
  ## plotting
  
  
  g1_3<- ggplot()+
    geom_hline(yintercept=0, linetype="dotted") + 
    geom_vline(xintercept=0,  linetype="dotted")+
    geom_sf(data=st_as_sf(KDE.UD1_3), colour=rainbow(kval), alpha=0.5)+
    geom_point(data=aus7[1:(nrow(aus7)-kval),], aes(x=Axis1, y=Axis3, colour=factor(jc_match)), shape=3)+
    geom_point(data=aus7[(nrow(aus7)-(kval-1)):nrow(aus7),],
               aes(x=Axis1, y=Axis3, fill=factor(jc_match)), shape=21, colour='black', size=3)+
    theme_bw()+theme(legend.position = "none")+
    scale_color_manual(values=rainbow(kval))
  
  eig<-aus7_pca$eig
  g1_3<- g1_3+scale_y_continuous(paste('PC3', sprintf('(%0.1f%% explained var.)', 100* eig[3]/sum(eig))))+
    scale_x_continuous(paste('PC1', sprintf('(%0.1f%% explained var.)', 100* eig[1]/sum(eig))))
  
  ### PCA 2- PCA3
  
  #kernel calc
  
  spdf<-SpatialPointsDataFrame(coords=cbind(aus7$Axis2, aus7$Axis3),
                               data=data.frame(jc_match=factor(aus7$jc_match)))
  
  KDE.Surface <- kernelUD(spdf,same4all = T, h=0.05, grid=1000)
  
  k2_3<-kernel.area(KDE.Surface, percent = c(50, 99))
  
  KDE.UD2_3 <- getverticeshr(KDE.Surface, percent = 99)
  
  ## plotting
  
  
  g2_3<- ggplot()+
    geom_hline(yintercept=0, linetype="dotted") + 
    geom_vline(xintercept=0,  linetype="dotted")+
    geom_sf(data=st_as_sf(KDE.UD2_3), colour=rainbow(kval), alpha=0.5)+
    geom_point(data=aus7[1:(nrow(aus7)-kval),], aes(x=Axis2, y=Axis3, colour=factor(jc_match)), shape=3)+
    geom_point(data=aus7[(nrow(aus7)-(kval-1)):nrow(aus7),],
               aes(x=Axis2, y=Axis3, fill=factor(jc_match)), shape=21, colour='black', size=3)+
    theme_bw()+theme(legend.position = "none")+
    scale_color_manual(values=rainbow(kval))
  
  eig<-aus7_pca$eig
  g2_3<- g2_3+scale_y_continuous(paste('PC3', sprintf('(%0.1f%% explained var.)', 100* eig[3]/sum(eig))))+
    scale_x_continuous(paste('PC2', sprintf('(%0.1f%% explained var.)', 100* eig[2]/sum(eig))))
  
  return(list(area1.2=k1_2, area1.3=k1_3, area2.3=k2_3, g1.2=g1_2, g1.3=g1_3, g2.3=g2_3,
              UD1_2=KDE.UD1_2, UD1_3=KDE.UD1_3,UD2_3=KDE.UD2_3, full_cent=aus7[(nrow(aus7)-(kval-1)):nrow(aus7),]))
}
