library(shiny)
library(plotly)
library(shinycssloaders)
anim.data = read.csv("anim_data.csv")
interp.resolution = 10
window.length = 10
regions = c( "EMR",  "EUR",  "AFR" , "AMR",  "WPR",  "SEAR")
# Define UI for application that plots random distributions
shinyUI(pageWithSidebar(

    # Application title
    headerPanel("Measles data 1990 - 2017"),

    # Sidebar with a slider input for number of observations
    sidebarPanel(
        radioButtons(
            inputId = "Scale",
            label   = "Plot Scale:",
            choices =  c("Non-scaled" = 'ns_',"Log Incidence" = 'li_'),
            selected="ns_",
            inline=TRUE
        ),
        checkboxGroupInput(
            inputId = "display_options",
            label   = "Trajectory Displays",
            # choices = c("Afterimage","2 Year Tail"),
            choices = c("Afterimage"),
            selected= c(""),
            inline=T
        ),

        checkboxGroupInput(
            inputId = "who_regions",
            label   = "Display Countries from which WHO region(s):",
            choices = regions,
            selected= regions
        ),
        tags$div(class = "header", checked = NA,
          tags$p(tags$b("Select Countries:")),
          actionButton("selectall", label="Select/Deselect all"),
          withSpinner(uiOutput('countries_selector'))
        )
    ),

    # Show a plot of the generated distribution
    mainPanel(
        # dataTableOutput('table1'),
        withSpinner(plotlyOutput("plot1"))
    )
))
