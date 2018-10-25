#

library(shiny)
library(dplyr)
library(ggplot2)
#https://daijiang.name/en/2014/05/11/functional-diversity-in-r/
library(ade4)
library(cluster)
library(factoextra)
library(FD)


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]}
  
# read in coral and fish abundance and trait data

fishtrait<-read.csv('./data/JPN_AUS_RMI_trait_master.csv', h=T)
row.names(fishtrait)<-fishtrait$Species
fishtrait$Species<-NULL


shinyServer(function(input, output) {
   
  mydata=reactive({

    fishtrait_sub<-fishtrait[which(rowSums(fishtrait[,input$reg_dat, drop=FALSE])>0),]

    md=fishtrait_sub[,input$trait_dat, drop=FALSE]
    return(md)
  })
  
  fdist=reactive({
    gowdis(mydata())
  })
  
  
  clusters = reactive ({hclust (fdist(), method = "ward.D")}) 
  
  phyly= reactive({
    FDisClus_phy = as.phylo(clusters())
    labelClust = as.character (row.names(mydata()))
    FDisClus_phy$tip.label = labelClust
    return(FDisClus_phy)
  })
  
  cutty=reactive ({hcut(fdist(), k=input$clusters, isdiss=T, hc_method='ward.D')})
  
  cutty2<-reactive({cutree(clusters(), k=input$clusters)})
  
  
  output$plotgraph1 = renderPlot({
    colz<-gg_color_hue(input$clusters)
    cutz<-cutty2()
    plot(phyly(), type="fan", use.edge.length = TRUE, node.pos = NULL,
         show.tip.label = TRUE, show.node.label = FALSE, tip.color = colz[cutz] ,
         edge.width = 1, edge.lty = 1, font = 2, cex = 0.5, label.offset = 0.1)
    
    
  })
  
  output$plotgraph2 = renderPlot({
    fviz_silhouette(cutty())
  })
  
})
