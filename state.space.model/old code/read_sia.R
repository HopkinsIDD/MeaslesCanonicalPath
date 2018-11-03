sia.data<-read.csv("input-SIA data with population numbers.csv",header=T)

cases<-read.csv("input-cases.csv",header=T,row.names=1)

mcv_5<-read.csv("input-mcv1x5yrs.csv",header=T,row.names=1)
codes<-read.csv("input-gbd_Region.csv",header=T,row.names=1)
agedist<-read.csv("input-agedist_pred.csv",header=T)
codes_cov<-agedist[,1:4]

agedist.sum<-apply(agedist[,5:8],1,cumsum)
age_u1<-agedist.sum[1,]
age_u5<-agedist.sum[2,]
age_u10<-agedist.sum[3,]
age_u15<-agedist.sum[4,]

country.codes<-rownames(cases)
cty_gbd<-codes[,2]

sia.coverage<-matrix(0,nr=dim(cases)[1], nc=dim(cases)[2])
sia.effective<-matrix(0,nr=dim(cases)[1], nc=dim(cases)[2])
for(i in 1:dim(sia.data)[1]){
#browser()
	rownumber<-which(country.codes==sia.data$iso[i])
	colnumber<-sia.data$year[i]-1980+1

	if(mcv_5[rownumber,colnumber] <	.60){t<-1}
	if(mcv_5[rownumber,colnumber] >= .60 & mcv_5[rownumber,colnumber] < .85){t<-2}
	if(mcv_5[rownumber,colnumber] >= .85){t<-3}
	
	tmp<-paste(cty_gbd[rownumber],as.character(t),sep="")
	ind<-which(codes_cov[,2]==tmp)
	
	sia.coverage[rownumber,colnumber]<-(min(1,sia.data$coverage[i]/100)) * .95
	
	if(sia.data$agemax[i] <= 2) {sia.effective[rownumber,colnumber]<-.5*(age_u1[ind]+age_u5[ind])}
	if(sia.data$agemax[i] > 2 && sia.data$agemax[i] <=8) {sia.effective[rownumber,colnumber]<-age_u5[ind]}
	if(sia.data$agemax[i] > 8 && sia.data$agemax[i] <= 13) {sia.effective[rownumber,colnumber]<-age_u10[ind]}
	if(sia.data$agemax[i] > 13) {sia.effective[rownumber,colnumber]<-age_u15[ind]}
}
##########################################################################################
# Post-hoc changes to SIA coverage
##########################################################################################
sia.coverage[which(country.codes=="MHL"),(2003-1980+1)]<-sia.coverage[which(country.codes=="MHL"),(2002-1980+1)] # switch Marshall Islands SIA 2002 to 2003
sia.coverage[which(country.codes=="MHL"),(2002-1980+1)]<-0 # set Marshall islands 2002 to 0% in 2002
sia.coverage[which(country.codes=="NAM"),(2009-1980+1)]<-.5 # set Namibia to 50% in 2009
sia.coverage[which(country.codes=="ZWE"),(2009-1980+1)]<-.5 # set Zimbabwe to 50% in 2009
sia.coverage[which(country.codes=="LSO"),(2011-1980+1)]<-sia.coverage[which(country.codes=="LSO"),(2010-1980+1)] # switch SIA 2010 to 2011
sia.coverage[which(country.codes=="LSO"),(2010-1980+1)]<-0 # set Lesotho to 0% in 2010
sia.coverage[which(country.codes=="MWI"),(2011-1980+1)]<-sia.coverage[which(country.codes=="MWI"),(2010-1980+1)] # switch SIA 2010 to 2011
sia.coverage[which(country.codes=="MWI"),(2010-1980+1)]<-0 # set Malawi to 0% in 2010
sia.coverage[which(country.codes=="SWZ"),(2011-1980+1)]<-sia.coverage[which(country.codes=="SWZ"),(2010-1980+1)] # switch SIA 2010 to 2011
sia.coverage[which(country.codes=="SWZ"),(2010-1980+1)]<-0 # set Swaziland to 0% in 2010
sia.coverage[which(country.codes=="ZAF"),(2011-1980+1)]<-sia.coverage[which(country.codes=="ZAF"),(2010-1980+1)] # switch SIA 2010 to 2011
sia.coverage[which(country.codes=="ZAF"),(2010-1980+1)]<-0 # set South Africa to 0% in 2010
sia.coverage[which(country.codes=="ZMB"),(2011-1980+1)]<-sia.coverage[which(country.codes=="ZMB"),(2010-1980+1)] # switch SIA 2010 to 2011
sia.coverage[which(country.codes=="ZMB"),(2010-1980+1)]<-0 # set Zambia to 0% in 2010
sia.coverage[which(country.codes=="ZWE"),(2011-1980+1)]<-sia.coverage[which(country.codes=="ZWE"),(2010-1980+1)] # switch SIA 2010 to 2011
sia.coverage[which(country.codes=="ZWE"),(2010-1980+1)]<-0 # set Zimbabwe to 0% in 2010
sia.coverage[which(country.codes=="COD"),(2011-1980+1)]<-0 # set DR Congo to 0% in 2011


sia.new<-sia.coverage*sia.effective
write.table(sia.coverage,file="output-sia_coverage.csv",sep=",",row.names=rownames(cases),col.names=1980:(1979+dim(sia.coverage)[2]))
write.table(sia.effective,file="output-sia_effective.csv",sep=",",row.names=rownames(cases),col.names=1980:(1979+dim(sia.coverage)[2]))
write.table(sia.new,file="output-sia_age_corrected.csv",sep=",",row.names=rownames(cases),col.names=1980:(1979+dim(sia.coverage)[2]))



