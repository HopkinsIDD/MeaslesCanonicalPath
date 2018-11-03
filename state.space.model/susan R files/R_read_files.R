# Reading files into R

#set the working directory
# from the dropdown menu or using 

setwd("/Users/matthewferrari/Dropbox/GAVI measles project/current burden estimation code")

#read in the files you want to use

cases <- read.csv("input-cases.csv", row.names=1, header=T)

# get rows of data

cases[1,] 

##########################################################
# naive way
# set up a new storage item for "Angola"

angola <- data.frame(cases=unlist(cases[2,]))

# add other columns in a similar way

births <- read.csv("input-births.csv", row.names=1, header=T)

# create a new drawer called "births" -- put births in it
angola$births <- unlist(births[2,]) 

# and on and on 

###########################################################
###########################################################
###########################################################
# specify index by ISO

country <- "AGO"

cases <- read.csv("input-cases.csv", row.names=1, header=T)
births <- read.csv("input-births.csv", row.names=1, header=T)
mcv1 <- read.csv("input-mcv1.csv", row.names=1, header=T)
mcv2 <- read.csv("input-mcv2.csv", row.names=1, header=T)
# . . . .etc


################################################################
# Make a cabinet - starting with case
country.index <- which(rownames(cases)==country)
country.summary <- data.frame(cases=unlist(cases[country.index,]))

################################################################
# Do births
country.index <- which(rownames(births)==country)
country.summary$births <- unlist(births[country.index,]) 

################################################################
# Do mcv1
country.index <- which(rownames(mcv1)==country)
country.summary$mcv1 <- unlist(mcv1[country.index,]) 

################################################################
# Do mcv2
country.index <- which(rownames(mcv2)==country)
country.summary$mcv2 <- unlist(mcv2[country.index,]) 








