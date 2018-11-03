# Reading files into R

#clear the memory
rm(list=ls())

#set the working directory
# from the dropdown menu or using 
setwd("/Users/matthewferrari/Dropbox/GAVI measles project/current burden estimation code")

#read in the files you want to use

pop <- read.csv("input-populations.csv", row.names=1, header=T)
births <- read.csv("input-births.csv", row.names=1, header=T)
deaths_5 <- read.csv("input-deaths_5.csv", row.names=1, header=T)

# change in population from year t to t+1
dpop <- pop[,-1] - pop[,-dim(pop)[2]]

# total deaths -- remainder of births minus population change
tdeaths <- births[,-dim(pop)[2]] - dpop

# deaths > 5 years is the total deaths minus the deaths under 5
deaths_over5 <- tdeaths - deaths_5[,-dim(pop)[2]]


# some sanity checks
pdeaths <- deaths_over5 / pop[,-dim(pop)[2]]
pdeaths <- (deaths_5[,-dim(pop)[2]]+deaths_over5) / births[,-dim(pop)[2]]
rdeaths <- deaths_5[,-dim(pop)[2]]/deaths_over5












