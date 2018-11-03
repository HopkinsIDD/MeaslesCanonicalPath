rm(list=ls())

#set the working directory
# from the dropdown menu or using 
setwd("/Users/matthewferrari/Dropbox/GAVI measles project/current burden estimation code")
wd<-getwd()
#read in the files you want to use

pop <- read.csv("input-populations.csv", row.names=1, header=T)
births <- read.csv("input-births.csv", row.names=1, header=T)
deaths_under5 <- read.csv("input-deaths_5.csv", row.names=1, header=T)

##########################################################################################
# Deaths over 5
##########################################################################################
# change in population from year t to t+1
dpop <- pop[,-1] - pop[,-dim(pop)[2]]
# total deaths -- remainder of births minus population change
tdeaths <- births[,-dim(pop)[2]] - dpop
# deaths > 5 years is the total deaths minus the deaths under 5
deaths_over5 <- tdeaths - deaths_under5[,-dim(pop)[2]]

cases <- read.csv("input-cases.csv", row.names=1, header=T)
mcv1 <- read.csv("input-mcv1.csv", row.names=1, header=T)
mcv2 <- read.csv("input-mcv2.csv", row.names=1, header=T)
sia <- read.csv("input-sia(disabled).csv", row.names=1, header=T)
first.surveillance <- read.csv("input-first.surveillance.csv", row.names=1, header=T)


country <- "RWA"  # good candidates are COG, RWA, ZWE

index<-which(rownames(cases)==country)

country.list<-list(
cases = cases[index,],
population = pop[index,],
births = births[index,],
deaths_under5 = deaths_under5[index,],
deaths_over5 = deaths_over5[index,],
mcv1 = mcv1[index,],
mcv2 = mcv2[index,],
sia = sia[index,]*.95,	# initially SIA is coverage per year for children <5y -- will expand later
first.surveillance = first.surveillance[index,1]   # only relevant for countries missing early data
)

filename<-paste(country,"data_list.Rdat",sep="")

setwd(paste(wd,"/example countries",sep=""))
save(country.list,file=filename)
setwd(wd)
