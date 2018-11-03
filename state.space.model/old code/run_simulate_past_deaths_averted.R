wd<-getwd()
######################################################################################################
######################################################################################################
### Simulate cases for the null case with no vaccination
######################################################################################################
######################################################################################################
#
# This section calculates the cases and deaths that would have occurred without vaccination
# This utilizes 2 additional functions "simulate_incidence", and "estimate_mortality_null".  
# NOTE: in the call for these functions, you can specify the name of the output files in 2 steps,  
#       a base name called "binary.file" and a file extension called "file.ext", which allows you to 
#       run multiple versions of this, outputting differently named files each time by changing either of these
#       two.  e.g. make "file.ext" equal to "scn1", "scn2", "scn3" for three separate scenarios.  
#		NOTE: that when you make the call to "estimate_mortality_null" you source the appropriate output from 
#             "simulate_incidence" if you are running multiple scenarios
#		
# 
###################################################
### calculate-mortality-without-vaccination: for all countries using mean parameters for AMRO and EURO
###################################################
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
#source("read_sia.R")	# generate age-corrected SIA coverage
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
outbreaks<-NA
load("output-incidence.Rdat")
source("simulate_cases.R")
source("simulate_incidence.R")
ignore.codes<-read.csv("input-ignore.indices.csv",header=T)
ignore.indices<-match(ignore.codes$iso,rownames(allcases))
PAR.t<-PAR
mean.par<-apply(PAR,2,mean,na.rm=T)
for(i in ignore.indices){PAR.t[i,]<-mean.par}
###################################################
### Currently implemented with vaccination set to 0 -- allmcv1=allmcv1*0, allmcv2=allmcv2*0, allsia=allsia*0 in function call
### This can be changed by either input hypothetical vaccination matrices (in the allmcv1, allmcv2, allsia slots) -- the most general option
### or multiply those matrices by a some proportional reduction (as I have done multiplying them all by 0 below) -- which imposes the same level of reduction everywhere
### NOTE: that "ignore.indices" is set to NA -- so this will estimate incidence for all countries
###################################################
setwd(paste(wd,"/past deaths averted",sep=""))
simulate_incidence(binary.file="output-incidence",file.ext="-null", PAR.t,allcases, allbirths, alldeaths,allpop, allmcv1=allmcv1*0, mcv1.12mo, allmcv2=allmcv2*0, allsia=allsia*0, ignore.indices=NA,display=TRUE,ITER=1)
setwd(wd)
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
### calculate-mortality-without-vaccination
###################################################
# This is done using the regional age distributions for <60% coverage from Abhijeet and Kathleen
# NOTE: the way we've been doing this is to apply an age distribution for each region and coverage class.  The question is then, under the assumption of no
#       measles vaccination, should we use A) the age distribution associated with no vaccination, or B) the age distribution as if there had been vaccination
#       (i.e. assuming that the age distribution might be reflective of other changes in the health system).
#		RIGHT NOW I have this implemented according to scenario A.  The age distribution is chosen based on the 5-year average coverage (that Abhijeet and Kathleen used).
#       So this is conservative.  If you want to go with scenario A, you would have to input a different matrix for "mcv1x5yrs.csv" OR, if you assuming NO vaccination,
#       simply multiply that matrix by 0 (e.g. mcv_5<-mcv_5*0)
# 
# I have currently implemented this using the cfr's indexed to under 5 mortality 
# Alternate scenarios can be run by simply changing the input files for the regional age distributions and/or CFR's
# e.g. to use a contant age distribution everywhere, the simplest would be to create a file with the same form as "agedist_pred.csv" but with the same 
#      age distribution in all rows.  NOTE: that you must also specify the "agedist_se.csv" file -- the standard errors for the age distribution -- setting
#	   these all to 0 (or at least a really small value) should work if you have an arbitrary age distribution you want to use, but no corresponding standard errors
# 
# NOTE: that "ignore.indices" is set to NA so that all countries are run
source("estimate_mortality_null.R")
setwd(paste(wd,"/past deaths averted",sep=""))
estimate_mortality_null(binary.file="output-mortality","-null",incidence_workspace=paste(wd,"/past deaths averted/output-incidence-null.Rdat",sep=""),agedist, agedist.se, cfr1,cfr2,cfr3, mcv_5*0, codes,ignore.indices=NA, ITER=100, display=TRUE)
load("output-mortality-null.Rdat")
write.csv(measles_deaths,file="output-mortality-null.csv",row.names=rownames(allcases))
setwd(wd)

######################################################################################################
######################################################################################################
### Simulate cases for the null case with only MCV1
######################################################################################################
######################################################################################################
# description as above
###################################################
###################################################
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
outbreaks<-NA
load("output-incidence.Rdat")
source("simulate_cases.R")
source("simulate_incidence.R")
ignore.codes<-read.csv("input-ignore.indices.csv",header=T)
ignore.indices<-match(ignore.codes$iso,rownames(allcases))
PAR.t<-PAR
mean.par<-apply(PAR,2,mean,na.rm=T)
for(i in ignore.indices){PAR.t[i,]<-mean.par}
###################################################
### Currently implemented with mcv2 and sia set to 0 --  allmcv2=allmcv2*0, allsia=allsia*0 in function call
### This can be changed by either input hypothetical vaccination matrices (in the allmcv1, allmcv2, allsia slots) -- the most general option
### or multiply those matrices by a some proportional reduction (as I have done multiplying them all by 0 below) -- which imposes the same level of reduction everywhere
### NOTE: that "ignore.indices" is set to NA -- so this will estimate incidence for all countries
###################################################
setwd(paste(wd,"/past deaths averted",sep=""))
simulate_incidence(binary.file="output-incidence",file.ext="-mcv1", PAR.t,allcases, allbirths, alldeaths,allpop, allmcv1=allmcv1, mcv1.12mo, allmcv2=allmcv2*0, allsia=allsia*0, ignore.indices=NA,display=TRUE,ITER=1)
setwd(wd)
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
### calculate-mortality-without-vaccination
###################################################
# Mortality is calculated using the age distribution corresponding to the region and MCV1 coverage
# 
# NOTE: that "ignore.indices" is set to NA so that all countries are run
source("estimate_mortality_null.R")
setwd(paste(wd,"/past deaths averted",sep=""))
estimate_mortality_null(binary.file="output-mortality","-mcv1",incidence_workspace=paste(wd,"/past deaths averted/output-incidence-mcv1.Rdat",sep=""),agedist, agedist.se, cfr1,cfr2,cfr3, mcv_5, codes,ignore.indices=NA, ITER=100, display=TRUE)
load("output-mortality-mcv1.Rdat")
write.csv(measles_deaths,file="output-mortality-mcv1.csv",row.names=rownames(allcases))
setwd(wd)

######################################################################################################
######################################################################################################
### Simulate cases for the observed vaccination levels 
######################################################################################################
######################################################################################################
# description as above
###################################################
###################################################
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
outbreaks<-NA
load("output-incidence.Rdat")
source("simulate_cases.R")
source("simulate_incidence.R")
ignore.codes<-read.csv("input-ignore.indices.csv",header=T)
ignore.indices<-match(ignore.codes$iso,rownames(allcases))
PAR.t<-PAR
mean.par<-apply(PAR,2,mean,na.rm=T)
for(i in ignore.indices){PAR.t[i,]<-mean.par}
###################################################
### Currently implemented with vaccination set to 0
### This can be changed by either input hypothetical vaccination matrices (in the allmcv1, allmcv2, allsia slots) -- the most general option
### or multiply those matrices by a some proportional reduction (as I have done multiplying them all by 0 below) -- which imposes the same level of reduction everywhere
### NOTE: that "ignore.indices" is set to NA -- so this will estimate incidence for all countries
###################################################
setwd(paste(wd,"/past deaths averted",sep=""))
simulate_incidence(binary.file="output-incidence",file.ext="-full", PAR.t,allcases, allbirths, alldeaths,allpop, allmcv1=allmcv1, mcv1.12mo, allmcv2=allmcv2, allsia=allsia, ignore.indices=NA,display=TRUE,ITER=1)
setwd(wd)
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
### calculate-mortality-without-vaccination
###################################################
# Mortality is calculated using the age distribution corresponding to the region and MCV1 coverage
# 
# NOTE: that "ignore.indices" is set to NA so that all countries are run
source("estimate_mortality_null.R")
setwd(paste(wd,"/past deaths averted",sep=""))
estimate_mortality_null(binary.file="output-mortality","-full",incidence_workspace=paste(wd,"/past deaths averted/output-incidence-full.Rdat",sep=""),agedist, agedist.se, cfr1,cfr2,cfr3, mcv_5, codes,ignore.indices=NA, ITER=100, display=TRUE)
load("output-mortality-full.Rdat")
write.csv(measles_deaths,file="output-mortality-full.csv",row.names=rownames(allcases))
setwd(wd)




##########################################################################################
postscript(file=paste(wd,"/figures/past_deaths_averted.eps",sep=""),height=8,width=8,horizontal=F)
load(paste(wd,"/past deaths averted/output-mortality-null.Rdat",sep="")) 
plot(2000:2012,apply(measles_deaths,2,sum,na.rm=T)[20:32],type="l",ylim=c(0,2e6),col=2,xlab="year",ylab="predicted measles deaths")
polygon(c(2000:2012,2012:2000),c(rep(0,13),apply(measles_deaths,2,sum,na.rm=T)[32:20]),col=grey(.75))
load(paste(wd,"/past deaths averted/output-mortality-mcv1.Rdat",sep=""))
polygon(c(2000:2012,2012:2000),c(rep(0,13),apply(measles_deaths,2,sum,na.rm=T)[32:20]),col=4)
lines(2000:2012,apply(measles_deaths,2,sum,na.rm=T)[20:32],col=4)
load(paste(wd,"/past deaths averted/output-mortality-full.Rdat",sep=""))
polygon(c(2000:2012,2012:2000),c(rep(0,13),apply(measles_deaths,2,sum,na.rm=T)[32:20]),col=2)
lines(2000:2012,apply(measles_deaths,2,sum,na.rm=T)[20:32],col=2)
dev.off()

