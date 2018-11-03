
library(shiny)

# Define UI for slider demo application
shinyUI(pageWithSidebar(

  #  Application title
  headerPanel("Single Country Summaries"),

  # Sidebar with sliders that demonstrate various available options
  sidebarPanel(
  
  selectInput("country", "Choose an country:", 
                choices = isos)
	),

  
  # plot output
  mainPanel(
  	tabsetPanel(
		tabPanel("Estimates",plotOutput("Est_Plot")),
		tabPanel("Transmission Parameter",plotOutput("transmission_Plot")),
		tabPanel("Reporting Parameter1",plotOutput("reporting1_Plot")),
		tabPanel("Reporting Parameter2",plotOutput("reporting2_Plot"))
		
		)
	)	
))