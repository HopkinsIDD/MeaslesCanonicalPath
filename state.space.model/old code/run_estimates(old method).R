setwd("/Users/matthewferrari/Documents/Current Projects/WHO/Burden Estimation/Geneva2013/working directory July 2013/")
require(MASS)
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
source("read_sia.R") # generate SIA matrix from input sia files
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
allsia<-read.csv("output-sia_coverage.csv",header=T,row.names=1)
first.surveillance<-unlist(read.csv("input-first.surveillance.csv",header=T,row.names=1)[,2])
outbreaks<-read.csv("input-outbreaks.csv",header=T,row.names=1) 
###################################################
### estimate-incidence
###################################################
ignore.indices<-c(4,5,6,8,9,10,13,18,19,21,22,23,24,25,26,30,31,32,39,40,42,43,44,45,46,48,49,50,52,55,56,58,60,63,70,71,72,73,74,75,76,77,80,83,84,85,86,88,89,91,94,95,100,101,102,104,105,106,108,112,116,118,124,128,130,131,134,135,137,138,142,144,145,148,153,156,157,159,161,162,163,164,166,172,173,174,175,176,182,183,185,186) #modified to excludes countries in lowest quartile of child mortality rates, those who have eliminated measles, and those with >95% coverage with 2 doses for >=5 years
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
estimate_mortality(binary.file="output-mortality.Rdat",incidence_workspace="output-incidence.Rdat",agedist, agedist.se, cfr1,cfr2,cfr3, mcv_5, codes,ignore.indices, ITER=1000, display=TRUE)
load("output-mortality.Rdat")

###################################################
### check against all-cause childhood mortality and define new outbreak matrix
###################################################
igme<-read.csv("input-igme.csv",header=T,row.names=1)[,1:(dim(allcases)[2]-1)]
# here is where you would change the "outbreak cutoff" from .20 to .10, etc
outbreak.cutoff<-.2
outbreaks.new<-cbind(outbreaks[,1],!(is.na(outbreaks[,2:(dim(allcases)[2])]) | measles_deaths>(outbreak.cutoff*igme)))  #NOTE: this will produce a warning -- it is not a problem, it happens because there are "NA" entries in outbreak
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
allsia<-read.csv("output-sia_coverage.csv",header=T,row.names=1)
first.surveillance<-unlist(read.csv("input-first.surveillance.csv",header=T,row.names=1)[,2])
outbreaks<-read.csv("output-outbreaks-new.csv",header=T,row.names=1) # the new file -- updated from the step above
###################################################
### estimate-incidence
###################################################
ignore.indices<-c(4,5,6,8,9,10,13,18,19,21,22,23,24,25,26,30,31,32,39,40,42,43,44,45,46,48,49,50,52,55,56,58,60,63,70,71,72,73,74,75,76,77,80,83,84,85,86,88,89,91,94,95,100,101,102,104,105,106,108,112,116,118,124,128,130,131,134,135,137,138,142,144,145,148,153,156,157,159,161,162,163,164,166,172,173,174,175,176,182,183,185,186) #modified to excludes countries in lowest quartile of child mortality rates, those who have eliminated measles, and those with >95% coverage with 2 doses for >=5 years
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
outbreak.cutoff<-.2
outbreaks.new<-cbind(outbreaks[,1],!(is.na(outbreaks[,2:(dim(allcases)[2])]) | measles_deaths>(outbreak.cutoff*igme)))  #NOTE: this will produce a warning -- it is not a problem, it happens because there are "NA" entries in outbreak
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
allsia<-read.csv("output-sia_coverage.csv",header=T,row.names=1)
first.surveillance<-unlist(read.csv("input-first.surveillance.csv",header=T,row.names=1)[,2])
outbreaks<-read.csv("output-outbreaks-new.csv",header=T,row.names=1) # the new file -- updated from the step above
###################################################
### estimate-incidence
###################################################
ignore.indices<-c(4,5,6,8,9,10,13,18,19,21,22,23,24,25,26,30,31,32,39,40,42,43,44,45,46,48,49,50,52,55,56,58,60,63,70,71,72,73,74,75,76,77,80,83,84,85,86,88,89,91,94,95,100,101,102,104,105,106,108,112,116,118,124,128,130,131,134,135,137,138,142,144,145,148,153,156,157,159,161,162,163,164,166,172,173,174,175,176,182,183,185,186) #modified to excludes countries in lowest quartile of child mortality rates, those who have eliminated measles, and those with >95% coverage with 2 doses for >=5 years
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
### END the fitting to data
###################################################

FILENAMES NOT CHANGED FROM THIS POINT DOWNWARD


# ######################################################################################################
# ######################################################################################################
# ### Simulate cases without vaccination
# ######################################################################################################
# ######################################################################################################
# #
# # This section calculates the cases and deaths that would have occurred without vaccination
# # This utilizes 2 additional functions "simulate_incidence", and "estimate_mortality_null".  
# # NOTE: in the call for these functions, you can specify the name of the output files in 2 steps,  
# #       a base name called "binary.file" and a file extension called "file.ext", which allows you to 
# #       run multiple versions of this, outputting differently named files each time by changing either of these
# #       two.  e.g. make "file.ext" equal to "scn1", "scn2", "scn3" for three separate scenarios.  
# #		NOTE: that when you make the call to "estimate_mortality_null" you source the appropriate output from 
# #             "simulate_incidence" if you are running multiple scenarios
# #		
# # 
# ###################################################
# ### calculate-mortality-without-vaccination: for all countries using mean parameters for AMRO and EURO
# ###################################################
# ###################################################
# ### load-functions
# ###################################################
# ####################################################
# # call all necessary functions and packages
# source("datainput.R") # functions to process the data
# source("Sfx.R") # functions to run the EKF
# source("filtersmooth.R") # an accessory function for the EKF
# source("ekf.optim.R") # a function to find the maximum likelihood
# source("filtered_trans.R") # accessory functions to construct confidence intervals
# source("estimate_incidence.R") # the function to return incidence estimates
# source("estimate_mortality.R") # the function to return mortality estimates
# ################################################
# ###################################################
# ###load-data
# ###################################################
# allcases<-read.csv("cases.csv",header=T,row.names=1)
# allbirths<-read.csv("births.csv",header=T,row.names=1)
# alldeaths<-read.csv("deaths_5.csv",header=T,row.names=1)
# allpop<-read.csv("populations.csv",header=T,row.names=1)
# allmcv1<-read.csv("mcv1.csv",header=T,row.names=1)
# mcv1.12mo<-read.csv("mcv1vaxeff.csv",header=T,row.names=1)
# allmcv2<-read.csv("mcv2.csv",header=T,row.names=1)
# allsia<-read.csv("output-sia_coverage.csv",header=T,row.names=1)
# first.surveillance<-unlist(read.csv("first.surveillance.csv",header=T,row.names=1)[,2])
# outbreaks<-read.csv("outbreaks-new.csv",header=T,row.names=1) 
# load("output-incidence.Rdat")
# source("simulate_cases.R")
# source("simulate_incidence.R")
# ignore.indices<-c(4,5,6,8,9,10,13,18,19,21,22,23,24,25,26,30,31,32,39,40,42,43,44,45,46,48,49,50,52,55,56,58,60,63,70,71,72,73,74,75,76,77,80,83,84,85,86,88,89,91,94,95,100,101,102,104,105,106,108,112,116,118,124,128,130,131,134,135,137,138,142,144,145,148,153,156,157,159,161,162,163,164,166,172,173,174,175,176,182,183,185,186) #modified to excludes countries in lowest quartile of child mortality rates, those who have eliminated measles, and those with >95% coverage with 2 doses for >=5 years
# PAR.t<-PAR
# mean.par<-apply(PAR,2,mean,na.rm=T)
# for(i in ignore.indices){PAR.t[i,]<-mean.par}
# ###################################################
# ### Currently implemented with vaccination set to 0
# ### This can be changed by either input hypothetical vaccination matrices (in the allmcv1, allmcv2, allsia slots) -- the most general option
# ### or multiply those matrices by a some proportional reduction (as I have done multiplying them all by 0 below) -- which imposes the same level of reduction everywhere
# ### NOTE: that "ignore.indices" is set to NA -- so this will estiamte incidence for all countries
# ###################################################
# simulate_incidence(binary.file="output-incidence",file.ext="-null", PAR.t,allcases, allbirths, alldeaths,allpop, allmcv1=allmcv1*0, mcv1.12mo, allmcv2=allmcv2*0, allsia=allsia*0, first.surveillance, outbreaks, ignore.indices=NA,display=TRUE,ITER=1)
# ###################################################
# ### load-mortality-data
# ###################################################
# agedist<-read.csv("agedist_pred.csv",header=T)
# agedist.se<-read.csv("agedist_se.csv",header=T)
# cfr1<-read.csv("cfr1-constant.csv",header=T,row.names=1) 
# cfr2<-read.csv("cfr2-constant.csv",header=T,row.names=1)
# cfr3<-read.csv("cfr3-constant.csv",header=T,row.names=1)
# mcv_5<-read.csv("mcv1x5yrs.csv",header=T,row.names=1)
# codes<-read.csv("gbd_Region.csv",header=T,row.names=1)
# ###################################################
# ### calculate-mortality-without-vaccination
# ###################################################
# # This is done using the regional age distributions for <60% coverage from Abhijeet and Kathleen
# # NOTE: the way we've been doing this is to apply an age distribution for each region and coverage class.  The question is then, under the assumption of no
# #       measles vaccination, should we use A) the age distribution associated with no vaccination, or B) the age distribution as if there had been vaccination
# #       (i.e. assuming that the age distribution might be reflective of other changes in the health system).
# #		RIGHT NOW I have this implemented according to scenario B.  The age distribution is chosen based on the 5-year average coverage (that Abhijeet and Kathleen used).
# #       So this is conservative.  If you want to go with scenario A, you would have to input a different matrix for "mcv1x5yrs.csv" OR, if you assuming NO vaccination,
# #       simply multiply that matrix by 0 (e.g. mcv_5<-mcv_5*0)
# # 
# # I have currently implemented this using the cfr's indexed to under 5 mortality 
# # Alternate scenarios can be run by simply changing the input files for the regional age distributions and/or CFR's
# # e.g. to use a contant age distribution everywhere, the simplest would be to create a file with the same form as "agedist_pred.csv" but with the same 
# #      age distribution in all rows.  NOTE: that you must also specify the "agedist_se.csv" file -- the standard errors for the age distribution -- setting
# #	   these all to 0 (or at least a really small value) should work if you have an arbitrary age distribution you want to use, but no corresponding standard errors
# # 
# # NOTE: that "ignore.indices" is set to NA so that all countries are run
# source("estimate_mortality_null.R")
# estimate_mortality_null(binary.file="output-mortality","-null",incidence_workspace="output-incidence-null.Rdat",agedist, agedist.se, cfr1,cfr2,cfr3, mcv_5, codes,ignore.indices=NA, ITER=100, display=TRUE)
# load("output-mortality-null.Rdat")



