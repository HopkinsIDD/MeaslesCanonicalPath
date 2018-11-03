setwd("/Users/mgraham/Dropbox/MeaslesTaxonomy/current burden estimation code") #change this to whatever you have named the file containing all the input files and scripts
require(MASS)
require(shiny)
wd<-getwd()

###################################################
### load-functions
###################################################
####################################################
# call all necessary functions and packages
source("datainput.R") # functions to process the data
source("Sfx.R") # functions to run the EKF
source("filtersmooth.R") # an accessory function for the EKF
source("ekf.optim.R") # a function to find the maximum likelihood
source("filtered_trans.R") # accessory functions to construct confidence intervals
source("estimate_incidence.R") # the function to return incidence estimates
source("estimate_mortality.R") # the function to return mortality estimates
source("read_sia.R")	# generate age-corrected SIA coverage
source("read_outbreaks.R")	# generate matrix of "outbreaks" or high reporting years
################################################
###################################################
###load-data
###################################################
allcases<-read.csv("input-cases.csv",header=T,row.names=1)
allbirths<-read.csv("input-births.csv",header=T,row.names=1)
alldeaths<-read.csv("input-deaths_5.csv",header=T,row.names=1)
allpop<-read.csv("input-populations.csv",header=T,row.names=1)
allmcv1<-read.csv("input-mcv1.csv",header=T,row.names=1)
mcv1.12mo<-read.csv("input-mcv1vaxeff.csv",header=T,row.names=1)
allmcv2<-read.csv("input-mcv2.csv",header=T,row.names=1)
allsia<-read.csv("output-sia_age_corrected.csv",header=T,row.names=1)
first.surveillance<-unlist(read.csv("input-first.surveillance.csv",header=T,row.names=1)[,2])
outbreaks<-read.csv("output-outbreaks.csv",header=T,row.names=1) 
###################################################
### estimate-incidence
###################################################
ignore.codes<-read.csv("input-ignore.indices.csv",header=T)
ignore.indices<-match(ignore.codes$iso,rownames(allcases))
incidence_estimate(binary.file="output-incidence.Rdat", allcases, allbirths, alldeaths,
allpop, allmcv1, mcv1.12mo, allmcv2, allsia, first.surveillance, outbreaks, ignore.indices,
display=TRUE)
load("output-incidence.Rdat")
###################################################
### load-mortality-data
###################################################
agedist<-read.csv("input-agedist_pred.csv",header=T)
agedist.se<-read.csv("input-agedist_se.csv",header=T)
cfr1<-read.csv("input-cfr1-constant.csv",header=T,row.names=1) 
cfr2<-read.csv("input-cfr2-constant.csv",header=T,row.names=1)
cfr3<-read.csv("input-cfr3-constant.csv",header=T,row.names=1)
mcv_5<-read.csv("input-mcv1x5yrs.csv",header=T,row.names=1)
codes<-read.csv("input-gbd_Region.csv",header=T,row.names=1)
###################################################
### estimate-mortality
###################################################
estimate_mortality(binary.file="output-mortality.Rdat",incidence_workspace="output-incidence.Rdat",
                   agedist, agedist.se, cfr1,cfr2,cfr3, mcv_5, codes,ignore.indices, 
                   ITER=1000, display=TRUE)
load("output-mortality.Rdat")

###################################################
### check against all-cause childhood mortality and define new outbreak matrix
###################################################
igme<-read.csv("input-igme.csv",header=T,row.names=1)[,1:(dim(allcases)[2]-1)]
# here is where you would change the "outbreak cutoff" from .20 to .10, etc
outbreak.cutoff<-.35
outbreaks.new<-!cbind(outbreaks[,1],(is.na(outbreaks[,2:(dim(allcases)[2])]) | measles_deaths>(outbreak.cutoff*igme)))  #NOTE: this will produce a warning -- it is not a problem, it happens because there are "NA" entries in outbreak
outbreaks.new[outbreaks.new==F]<-NA
write.table(outbreaks.new,file="output-outbreaks-new.csv",  sep=",",row.name=rownames(allcases),col.name=colnames(allcases))

######################################################################################################
######################################################################################################
### Re-run with the new outbreak matrix
######################################################################################################
######################################################################################################
######################################################################################################
### load-functions
###################################################
####################################################
# call all necessary functions and packages
source("datainput.R") # functions to process the data
source("Sfx.R") # functions to run the EKF
source("filtersmooth.R") # an accessory function for the EKF
source("ekf.optim.R") # a function to find the maximum likelihood
source("filtered_trans.R") # accessory functions to construct confidence intervals
source("estimate_incidence.R") # the function to return incidence estimates
source("estimate_mortality.R") # the function to return mortality estimates
################################################
###################################################
###load-data
###################################################
allcases<-read.csv("input-cases.csv",header=T,row.names=1)
allbirths<-read.csv("input-births.csv",header=T,row.names=1)
alldeaths<-read.csv("input-deaths_5.csv",header=T,row.names=1)
allpop<-read.csv("input-populations.csv",header=T,row.names=1)
allmcv1<-read.csv("input-mcv1.csv",header=T,row.names=1)
mcv1.12mo<-read.csv("input-mcv1vaxeff.csv",header=T,row.names=1)
allmcv2<-read.csv("input-mcv2.csv",header=T,row.names=1)
allsia<-read.csv("output-sia_age_corrected.csv",header=T,row.names=1)
first.surveillance<-unlist(read.csv("input-first.surveillance.csv",header=T,row.names=1)[,2])
outbreaks<-read.csv("output-outbreaks-new.csv",header=T,row.names=1) # the new file -- updated from the step above
###################################################
### estimate-incidence
###################################################
ignore.codes<-read.csv("input-ignore.indices.csv",header=T)
ignore.indices<-match(ignore.codes$iso,rownames(allcases))
incidence_estimate(binary.file="output-incidence.Rdat", allcases, allbirths, alldeaths,
allpop, allmcv1, mcv1.12mo, allmcv2, allsia, first.surveillance, outbreaks, ignore.indices,
display=TRUE)
load("output-incidence.Rdat")
###################################################
### load-mortality-data
###################################################
agedist<-read.csv("input-agedist_pred.csv",header=T)
agedist.se<-read.csv("input-agedist_se.csv",header=T)
cfr1<-read.csv("input-cfr1-constant.csv",header=T,row.names=1) 
cfr2<-read.csv("input-cfr2-constant.csv",header=T,row.names=1)
cfr3<-read.csv("input-cfr3-constant.csv",header=T,row.names=1)
mcv_5<-read.csv("input-mcv1x5yrs.csv",header=T,row.names=1)
codes<-read.csv("input-gbd_Region.csv",header=T,row.names=1)
###################################################
### estimate-mortality
###################################################
estimate_mortality(binary.file="output-mortality.Rdat",incidence_workspace="output-incidence.Rdat",agedist, agedist.se, cfr1,cfr2,cfr3, mcv_5, codes,ignore.indices, ITER=10000, display=TRUE)
load("output-mortality.Rdat")

###################################################
### check against all-cause childhood mortality and define new outbreak matrix
###################################################
igme<-read.csv("input-igme.csv",header=T,row.names=1)[,1:(dim(allcases)[2]-1)]
# here is where you would change the "outbreak cutoff" from .20 to .10, etc
outbreak.cutoff<-.35
outbreaks.new<-!cbind(outbreaks[,1],(is.na(outbreaks[,2:(dim(allcases)[2])]) | measles_deaths>(outbreak.cutoff*igme)))  #NOTE: this will produce a warning -- it is not a problem, it happens because there are "NA" entries in outbreak
outbreaks.new[outbreaks.new==F]<-NA
write.table(outbreaks.new,file="output-outbreaks-new.csv",  sep=",",row.name=rownames(allcases),col.name=colnames(allcases))

######################################################################################################
######################################################################################################
### Re-run with the new outbreak matrix
######################################################################################################
######################################################################################################
######################################################################################################
### load-functions
###################################################
####################################################
# call all necessary functions and packages
source("datainput.R") # functions to process the data
source("Sfx.R") # functions to run the EKF
source("filtersmooth.R") # an accessory function for the EKF
source("ekf.optim.R") # a function to find the maximum likelihood
source("filtered_trans.R") # accessory functions to construct confidence intervals
source("estimate_incidence.R") # the function to return incidence estimates
source("estimate_mortality.R") # the function to return mortality estimates
################################################
###################################################
###load-data
###################################################
allcases<-read.csv("input-cases.csv",header=T,row.names=1)
allbirths<-read.csv("input-births.csv",header=T,row.names=1)
alldeaths<-read.csv("input-deaths_5.csv",header=T,row.names=1)
allpop<-read.csv("input-populations.csv",header=T,row.names=1)
allmcv1<-read.csv("input-mcv1.csv",header=T,row.names=1)
mcv1.12mo<-read.csv("input-mcv1vaxeff.csv",header=T,row.names=1)
allmcv2<-read.csv("input-mcv2.csv",header=T,row.names=1)
allsia<-read.csv("output-sia_age_corrected.csv",header=T,row.names=1)
first.surveillance<-unlist(read.csv("input-first.surveillance.csv",header=T,row.names=1)[,2])
outbreaks<-read.csv("output-outbreaks-new.csv",header=T,row.names=1) # the new file -- updated from the step above
###################################################
### estimate-incidence
###################################################
ignore.codes<-read.csv("input-ignore.indices.csv",header=T)
ignore.indices<-match(ignore.codes$iso,rownames(allcases))
incidence_estimate(binary.file="output-incidence.Rdat", allcases, allbirths, alldeaths,
allpop, allmcv1, mcv1.12mo, allmcv2, allsia, first.surveillance, outbreaks, ignore.indices,
display=TRUE)
load("output-incidence.Rdat")
###################################################
### load-mortality-data
###################################################
agedist<-read.csv("input-agedist_pred.csv",header=T)
agedist.se<-read.csv("input-agedist_se.csv",header=T)
cfr1<-read.csv("input-cfr1-constant.csv",header=T,row.names=1) 
cfr2<-read.csv("input-cfr2-constant.csv",header=T,row.names=1)
cfr3<-read.csv("input-cfr3-constant.csv",header=T,row.names=1)
mcv_5<-read.csv("input-mcv1x5yrs.csv",header=T,row.names=1)
codes<-read.csv("input-gbd_Region.csv",header=T,row.names=1)
###################################################
### estimate-mortality
###################################################
estimate_mortality(binary.file="output-mortality.Rdat",incidence_workspace="output-incidence.Rdat",agedist, agedist.se, cfr1,cfr2,cfr3, mcv_5, codes,ignore.indices, ITER=10000, display=TRUE)
load("output-mortality.Rdat")

###################################################
### END the fitting to data -- launch graphics
###################################################
#runApp("single country output")
###################################################
### 
###################################################



