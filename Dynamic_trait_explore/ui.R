#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  titlePanel('Dynamic clustering'),
  
  sidebarLayout(position='left',
                sidebarPanel('Select Traits',
                             checkboxInput("HomeRange", "HomeRange", FALSE), 
                             checkboxInput("MaxLength", "MaxLength", FALSE), 
                             checkboxInput("PLD_Mean", "PLD_Mean", FALSE), 
                             checkboxInput("TrophLevel", "TrophLevel", FALSE), 
                             checkboxInput("Food", "Food", FALSE), 
                             checkboxInput("Function", "Function", FALSE), 
                             checkboxInput("Aggregation", "Aggregation", FALSE), 
                             checkboxInput("Position", "Position", FALSE), 
                             checkboxInput("SpawnMode", "SpawnMode", FALSE), 
                             checkboxInput("ParentalMode", "ParentalMode", FALSE), 
                             checkboxInput("Active", "Active", FALSE), 
                             sliderInput("clusters","No. of clusters",min=1,max=10,value=1)),
                
                mainPanel("main panel",
                          fluidRow(
                            plotOutput("plotgraph1", height='550px'), plotOutput("plotgraph2", height = '250px')
                            ))
               )
  ))

