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
                sidebarPanel( checkboxGroupInput('reg_dat', 'Select region',
                                                 c('Japan' = 'JPN_sp',
                                                   'Australia'= 'AUS_sp',
                                                   'RMI' = 'RMI_sp')),
                             checkboxGroupInput('trait_dat', 'Select traits',
                                                 c('Thermal Affinity' = 'ThermalAffinity',
                                                   'Body Size'= 'BodySize',
                                                    'PLD' = 'PLD',
                                                   'Diet' = 'Diet',
                                                   'Aggregation' = 'Aggregation',
                                                   'Position' = 'Position',
                                                   'Parental Mode' = 'ParentalMode')),
                        
                             sliderInput("clusters","No. of clusters",min=1,max=10,value=1)),
                
                mainPanel("Hierarchical clustering ring",
                          fluidRow(
                            plotOutput("plotgraph1", height='550px'), plotOutput("plotgraph2", height = '250px')
                            ))
               )
  ))

