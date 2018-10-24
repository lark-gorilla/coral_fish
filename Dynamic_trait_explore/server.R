#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(dplyr)
library(readxl)
library(ggplot2)
#https://daijiang.name/en/2014/05/11/functional-diversity-in-r/
library(ade4)
library(cluster)
library(factoextra)
library(FD)
library(shiny)


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]}
  
# read in coral and fish abundance and trait data

fishdat<-read_excel('~/leeds_postdoc/data/RMI/RMIfish_siteDepthMeans_ADepth.xlsx', sheet=1)

fishtrait<-read.csv('~/leeds_postdoc/data/Traits/RMI_fish_traits4.csv', h=T)

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

#filter(ingrid_out, has_trait==F)
# only Bryanops natans
# called Bryaninops natans in trait data - correct
names(fishdat)[names(fishdat)=='Bryanops natans']<-'Bryaninops natans'

fishtrait2<-filter(fishtrait, Species%in% names(fishdat))

row.names(fishtrait2)<-fishtrait2$Species


# Define server logic required to draw a histogram
shinyServer(function(input, output) {
   
  mydata=reactive({
    selnames=NULL
    if(input$HomeRange){selnames=c(selnames, 1)}
    if(input$MaxLength){selnames=c(selnames, 2)}
    if(input$PLD_Mean){selnames=c(selnames, 3)}
    if(input$TrophLevel){selnames=c(selnames, 4)}
    if(input$Food){selnames=c(selnames, 5)}
    if(input$Function){selnames=c(selnames, 6)}
    if(input$Aggregation){selnames=c(selnames, 7)}
    if(input$Position){selnames=c(selnames, 8)}
    if(input$SpawnMode){selnames=c(selnames, 9)}
    if(input$ParentalMode){selnames=c(selnames, 10)}
    if(input$Active){selnames=c(selnames, 11)}
    
    md=fishtrait2[,selnames]
    return(md)
  })
  
  fdist=reactive({
    gowdis(mydata())
  })
  
  
  clusters = reactive ({hclust (fdist(), method = "ward.D")}) 
  
  phyly= reactive({
    FDisClus_phy = as.phylo(clusters())
    labelClust = as.character (row.names(fishtrait2))
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
