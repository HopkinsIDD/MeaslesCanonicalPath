library(shiny)
library(plotly)
library(RColorBrewer)
library(dplyr)
library(scales)
library(tidyr)
anim.data = read.csv("anim_data.csv")
regions = c("EMR","EUR","AFR","AMR","WPR","SEAR")
min=0.000001
scl <- function(x,y=x){ (x - min(y,na.rm = T))/(max(y,na.rm = T) - min(y,na.rm = T)) }

### Commenting out the code used to create data frames, but leaving it in the file in case they need to be recreated for some reason
# scl <- function(x){ (x - min(x,na.rm = T))/(max(x,na.rm = T) - min(x,na.rm = T)) }
# 
# prep.anim.data <- function(d){
#   ##' filter the data to only include the regions, and years required, along with only
#   ##' entries which have incidence greater than 0.
#   ##' d = d %>% filter(., Year %in% years, WHO_REGION %in% regions, Incidence > 0,
#   ##'                 Country != "")
#   d = d %>% filter(., Year >=1990, WHO_REGION %in% regions, Incidence > 0, Country != "")
#   
#   ##' remove years where the coefficient of variation is 0, and the year is before 1995,
#   ##' as before this time, no countries had reached elimination, so the only way to 
#   ##' have a 0 here, is if there is no data before this time (unless a country reported
#   ##' the exact same non-zero number of cases for every year, which didn't happen)
#   j = which(d$Year < 1995 & d$Coefficient.of.Variation == 0)
#   if(length(j) > 0){
#     d = d[-(j), ]
#   }
# }
# 
# rescale <- function(
#   d,
#   regions = c("EMR","EUR","AFR","AMR","WPR","SEAR"),
#   make.inc.cv.scale.same = TRUE,
#   sqrt.inc = F,
#   sqrt.cv=F,
#   log.incidence = F){
#   
#   ##' if we want to log the incidence before assigning countries to their location
#   ##' on the path, then we do that here. Along with adding a small non-zero value to each of the
#   ##' incidences, so that there are no -Inf values created.
#   if(log.incidence == T){
#     d$Incidence = log(d$Incidence + 0.000001)
#   }
#   
#   ##' if sqrt.inc=T then take sqrt of incidence
#   if(sqrt.inc == T){
#     d$Incidence <-  sqrt(d$Incidence)
#   }
#   
#   ##' if sqrt.cv=T then take sqrt of cv
#   if(sqrt.cv == T){
#     d$Coefficient.of.Variation = sqrt(d$Coefficient.of.Variation)
#   }
#   
#   ##' if wanted scale the incidence and coefficient of variation so that they are both in the 0-1 range
#   if(make.inc.cv.scale.same == T){
#     d$Incidence = scl(d$Incidence)
#     d$Coefficient.of.Variation = scl(d$Coefficient.of.Variation)
#   }
#   return(d)
#   
# }
# 
# 
# rescale_shiny_app <- function(d, type){
#   if (type=='li_'){
#     d = rescale(d=d,
#                 make.inc.cv.scale.same = TRUE,
#                 sqrt.inc=F,
#                 sqrt.cv=F,
#                 log.incidence=T)
#   } else if (type=='ns_'){
#     d = rescale(d=d,
#                 make.inc.cv.scale.same = F,
#                 sqrt.inc=F,
#                 sqrt.cv=F,
#                 log.incidence=F)
#   } else {
#   }
#   return(d)
# }
# 
# anim.data.orig <- read.csv("anim_data.csv")[,(-1)]
# anim.data <- prep.anim.data(d=anim.data.orig)
# d_li <- rescale_shiny_app(d=anim.data, type="li_")
# anim.data$li_Incidence <- d_li$Incidence
# anim.data$li_Coefficient.of.Variation <- d_li$Coefficient.of.Variation
# d_ns <- rescale_shiny_app(d=anim.data, type="ns_")
# anim.data$ns_Incidence <- sqrt(d_ns$Incidence) #sqrt of non-scaled for the sake of display
# anim.data$ns_Coefficient.of.Variation <- d_ns$Coefficient.of.Variation
# 
# anim.data = arrange(anim.data,Year,Country)
# 
# ## Trailers
# ## 2 Year Line
# anim.data2 = anim.data %>% filter(Country != '') %>%
#   gather('Variable','value',Coefficient.of.Variation,Incidence,ns_Coefficient.of.Variation,ns_Incidence,li_Coefficient.of.Variation,li_Incidence) %>%
#   group_by(Country,Variable) %>%
#   arrange(Year) %>%
#   mutate(
#     `0` = value,
#     `1` = lag(value,1),
#     `2` = lag(value,2),
#     `3` = lag(value,3),
#     `20` = lag(value,20),
#     `40` = lag(value,40)
#   ) %>%
#   gather('Lag','value',`0`,`1`,`2`,`3`,`20`,`40`) %>%
#   filter(!is.na(value)) %>%
#   spread('Variable','value')
# 
# yearly.anim.data = anim.data %>% filter(Year == round(Year))
# yearly.anim.data2 = anim.data2 %>% filter(Year == round(Year))
# 
# write.csv(file='yearly.anim.data2.csv',yearly.anim.data2)
# write.csv(file='yearly.anim.data.csv',yearly.anim.data)
# 
# #Ran the first few chunks of RMarkdown to get canonical.path and canonical.path2
# cc <- data.frame(li_avg.x = canonical.path$x,
#                  li_avg.y = canonical.path$y,
#                  ns_avg.x = canonical.path2$x, 
#                  ns_avg.y = sqrt(canonical.path2$y)) #sqrt of non-scaled for the sake of display
# write.csv(cc, "cc.csv")

yearly.anim.data = read.csv('yearly.anim.data.csv')
yearly.anim.data2 = group_by(read.csv('yearly.anim.data2.csv'),Country,Year)
cc = read.csv('cc.csv')

## One last data processing step:
li_min_inc = .1
# li_Incidence
yearly.anim.data = yearly.anim.data %>% mutate(li_Incidence = ifelse(li_Incidence > li_min_inc,li_Incidence,li_min_inc))
yearly.anim.data2 = yearly.anim.data2 %>% mutate(li_Incidence = ifelse(li_Incidence > li_min_inc,li_Incidence,li_min_inc))

#yearly.anim = anim.data %>% filter(Year %in% seq(min(anim.data$Year), max(anim.data$Year),.5))
#
# gg <- ggplot(yearly.anim, aes(Coefficient.of.Variation, Incidence, color = Mean.vaccination, frame = Year,label = Country)) +
#     geom_point(aes(size = size)) + xlim(0,1) +
#     scale_colour_gradientn(colors = c('firebrick', 'red', 'yellow', 'skyblue', 'blue'))
colfunc <- colorRampPalette(c('firebrick', 'red', 'yellow', 'skyblue', 'blue'))
# Define server logic required to generate and plot a random distribution
shinyServer(function(input, output, session) {
  output$country_select_header<- renderText({
    HTML(paste0("<p>","<b>","Select Countries:","</b>","</p>"))
  })
  selected_countries = reactiveVal(NULL)
  output$countries_selector <- renderUI({
    checkboxGroupInput(
      inputId = 'countries', # ID to be used in server.R
      # label="Select Countries:",
      label='',
      choices= sort(unique(subset(anim.data, anim.data$WHO_REGION %in% input$who_regions)$Country))[-1],
      # selected = intersect(input$countries, sort(unique(subset(anim.data, anim.data$WHO_REGION %in% input$who_regions)$Country)[-1]))
      selected = sort(unique(subset(anim.data, anim.data$WHO_REGION %in% input$who_regions)$Country))[-1]
    )
  })

  observe({
    if (input$selectall > 0){
      if (input$selectall %% 2 == 0){
        updateCheckboxGroupInput(
          session  = session,
          inputId  = "countries",
          choices= sort(unique(subset(anim.data, anim.data$WHO_REGION %in% input$who_regions)$Country))[-1],
          selected = sort(unique(subset(anim.data, anim.data$WHO_REGION %in% input$who_regions)$Country))[-1]
        )
      } else {
        updateCheckboxGroupInput(
          session  = session,
          inputId  = "countries",
          choices= sort(unique(subset(anim.data, anim.data$WHO_REGION %in% input$who_regions)$Country))[-1],
          selected = c() 
        )
      }
    }
  })

    # Expression that generates a plot of the distribution. The expression
    # is wrapped in a call to renderPlot to indicate that:
    #
    #  1- It is "reactive" and therefore should be automatically
    #     re-executed when inputs change
    #  2- Its output type is a plot
    #

    # output$info <- renderPrint({
    #     val <- nearPoints(anim.data, input$plot_click, xvar = "Coefficient.of.Variation", yvar = "Incidence",
    #                      panelvar1 = "Country", threshold = 100, maxpoints = 1,addDist = TRUE)
    #
    #     paste(val)
    # })
    #
  altered.data3 = reactive({
    yearly.anim = yearly.anim.data2 %>%
      filter(
        WHO_REGION %in% input$who_regions,
        Country %in% input$countries,
        Year == round(Year)
      )
    if(!("2 Year Tail" %in% input$display_options)){
      yearly.anim = yearly.anim %>%
        filter(Lag == 0)
    } else {
      yearly.anim = yearly.anim %>%
        filter(Lag %in% c(0,20,40))
    }
    yearly.anim$id = as.factor(paste(yearly.anim$WHO_REGION,yearly.anim$Country,yearly.anim$Lag,sep='_'))
    yearly.anim = yearly.anim %>% filter(Country != '')
    colnames = c("Country",paste0(input$Scale,"Coefficient.of.Variation"),"Year",paste0(input$Scale,"Incidence"),"Mean.vaccination","Lag")
    yearly.anim = yearly.anim[,colnames]
    yearly.anim$Mean.vaccination = round(yearly.anim$Mean.vaccination)+1
    colnames(yearly.anim) = c("Country","Coefficient.of.Variation","Year","Incidence","Mean.vaccination","Lag")
    return(yearly.anim)
  })

  altered.data2 = reactive({
    yearly.anim = yearly.anim.data2 %>%
      filter(
        WHO_REGION %in% input$who_regions,
        Country %in% input$countries,
        Year == round(Year)
      )
    if(!("Afterimage" %in% input$display_options)){
      yearly.anim = yearly.anim %>%
        filter(Lag == 0)
    } else {
      yearly.anim = yearly.anim %>%
        filter(Lag %in% c(0,1,2,3))
    }
    yearly.anim$id = as.factor(paste(yearly.anim$WHO_REGION,yearly.anim$Country,yearly.anim$Lag,sep='_'))
    yearly.anim = yearly.anim %>% filter(Country != '')
    colnames = c("Country",paste0(input$Scale,"Coefficient.of.Variation"),"Year",paste0(input$Scale,"Incidence"),"Mean.vaccination","Lag")
    yearly.anim = yearly.anim[,colnames]
    yearly.anim$Mean.vaccination = round(yearly.anim$Mean.vaccination)+1
    colnames(yearly.anim) = c("Country","Coefficient.of.Variation","Year","Incidence","Mean.vaccination","Lag")
    return(yearly.anim)
  })

  altered.data = reactive({
    yearly.anim = yearly.anim.data %>%
      filter(
        WHO_REGION %in% input$who_regions,
        Country %in% input$countries,
        Year == round(Year)
      )
    yearly.anim$id = as.factor(paste(yearly.anim$WHO_REGION,yearly.anim$Country,sep='_'))
    yearly.anim = yearly.anim %>% filter(Country != '')
    colnames = c("Country",paste0(input$Scale,"Coefficient.of.Variation"),"Year",paste0(input$Scale,"Incidence"),"Mean.vaccination")
    yearly.anim = yearly.anim[,colnames]
    yearly.anim$Mean.vaccination = round(yearly.anim$Mean.vaccination)+1
    colnames(yearly.anim) = c("Country","Coefficient.of.Variation","Year","Incidence","Mean.vaccination")
    return(yearly.anim)
  })

  xbounds = reactive({
    if(input$Scale == 'ns_'){return(c(0,3.2))}
    if(input$Scale == 'li_'){return(c(0,1))}
    return(c(0,4))
  })

  ybounds = reactive({
    if(input$Scale == 'ns_'){return(c(0,sqrt(15)))}
    if(input$Scale == 'li_'){return(c(0.1,1))}
    return(c(0,4))
  })
  aspect_ratio = reactive({
    dy = ybounds()[1] - ybounds()[2]
    dx = xbounds()[1] - xbounds()[2]
    return(.8)
    # return(dy/dx)
  })
  output$table1 <- renderDataTable({
    altered.data2()
  })
  ytickvals = reactive({
    if(input$Scale == 'ns_'){return(sqrt(c(0, 0.1, 0.5, 1, 1.5, 3, 5, 7, 10, 15)))}
    if(input$Scale == 'li_'){return(c(0.1,.2,.4,.6,.8,1))}
    return(c())
  })
  yticktext = reactive({
    if(input$Scale == 'ns_'){return(c(0, 0.1,0.5, 1, 1.5, 3, 5, 7, 10, 15))}
    if(input$Scale == 'li_'){return(c('<.01',.2,.4,.6,.8,1))}
    return(c())
  })
  ylabel = reactive({
    if(input$Scale == 'ns_'){return("Mean Incidence per 1000")}
    if(input$Scale == 'li_'){return("Log Mean Incidence")}
    return(c())
  })

  ### xtickvals = reactive({
  ###   if(input$Scale == 'ns_'){return(sqrt(c(0, 0.1, 0.5, 1, 1.5, 3, 5, 7, 10, 15)))}
  ###   if(input$Scale == 'li_'){return(c(0,.2,.4,.6,.8,1))}
  ###   return(c())
  ### })
  ### xticktext = reactive({
  ###   if(input$Scale == 'ns_'){return(c(0, 0.1,0.5, 1, 1.5, 3, 5, 7, 10, 15))}
  ###   if(input$Scale == 'li_'){return(c(0,.2,.4,.6,.8,1))}
  ###   return(c())
  ### })

  output$plot1 <- renderPlotly({
    if(length(input$countries) > 0){
    tmp <- animation_opts(
      layout(
        colorbar(
          add_paths(
              add_trace(
                plot_ly(
                altered.data2(),
                x= ~Coefficient.of.Variation,
                y= ~Incidence,# xlim = c(0,1), ylim = c(0,1),
                type='scatter',
                mode='markers', #'markers'
                frame= ~Year,
                color=~I('Gray'),
                ids=~interaction(Country,Lag),
                opacity=.3,
                hoverinfo='none',
                colors = colfunc(101),
                  width = session$clientData$output_plot1_width,
                  height = session$clientData$output_plot1_width * aspect_ratio()
                ),
                  data=altered.data(),
                  x= ~Coefficient.of.Variation,
                  y= ~Incidence,# xlim = c(0,1), ylim = c(0,1),
                  type='scatter',
                  mode='markers', #'markers'
                  frame= ~Year,
                  color=~as.integer(Mean.vaccination),
                  text=~Country,
                  ids=~Country,
                  # color =colfunc(101)[round(altered.data()$Mean.vaccination)+1],
                  #col = ~colfunc(101)[round(Mean.vaccination)+1],
                  hoveron='points',
                  hoverinfo='text',
                  colors = colfunc(101),
                  opacity=.75,
                  showlegend=F,
                  #pch = 16
                  #xlab = "Coefficient of Variation",
                  #ylab = "Incidence", bty = 'n', cex = 3.5, cex.axis=2, cex.lab = 2
                inherit=FALSE
              ),
            x=cc[[paste0(input$Scale,"avg.x")]],
            y=cc[[paste0(input$Scale,"avg.y")]],
            # line= list(line-opacity=.65),
            color=I('black'),
            showlegend=F,
            inherit=FALSE
          ),
          title='Mean Vaccination',
          x=xbounds()[2],
          y=mean(.5)
        ),
        xaxis=list(range=xbounds(),title="Coefficient of Variation"),#tickvals = xtickvals(),ticktext=xticktext()),
        showlegend=F,
        yaxis=list(range=ybounds(),title=ylabel(),tickvals = ytickvals(),ticktext=yticktext())
      ),
      frame = 1000,
      transition = 0,
      redraw = TRUE
    )
    } else {
    tmp <- animation_opts(
      layout(
        colorbar(
          add_trace(
            plot_ly(
              x=cc[[paste0(input$Scale,"avg.x")]],
              y=cc[[paste0(input$Scale,"avg.y")]],
              type = 'scatter',
              mode = 'line',
              color=I('black'),
              width = session$clientData$output_plot1_width,
              height = session$clientData$output_plot1_width * aspect_ratio(),
              colors = colfunc(101),
              showlegend=F
            ),
            type='scatter',
            mode='markers', #'markers'
            opacity=0.0000001,
            x=-1,
            y=-1, 
            frame = sort(unique(round(anim.data$Year))),
            color=rep(1:101,10)[1:length(unique(round(anim.data$Year)))],
            colors=colfunc(101),
            showlegend=F,
            inherit=FALSE
          ),
          title='Mean Vaccination',
          x=xbounds()[2],
          y=mean(ybounds())
        ),
        xaxis=list(range=xbounds(),title="Coefficient of Variation"),
        showlegend=F,
        yaxis=list(range=ybounds(),title="Incidence")
      ),
      frame = 1000,
      transition = 0,
      redraw = TRUE
    )
    }
    return(tmp)
  })
})
