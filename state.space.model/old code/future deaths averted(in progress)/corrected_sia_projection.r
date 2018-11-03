#placeholder for future sia file
sia.data<-read.csv("input-SIA_schedule_for_simulation.csv",header=T)

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
	rownumber<-which(country.codes==sia.data$iso[i])
	colnumber<-sia.data$year[i]-1980+1

	if(mcv_5[rownumber,colnumber] <	.60){t<-1}
	if(mcv_5[rownumber,colnumber] >= .60 & mcv_5[rownumber,colnumber] < .85){t<-2}
	if(mcv_5[rownumber,colnumber] >= .85){t<-3}
	
	tmp<-paste(cty_gbd[rownumber],as.character(t),sep="")
	ind<-which(codes_cov[,2]==tmp)
	
	sia.coverage[rownumber,colnumber]<-(min(1,sia.data$coverage[i]))/100 * .95
	
	if(sia.data$agemax[i] <= 2) {sia.effective[rownumber,colnumber]<-.5*(age_u1[ind]+age_u5[ind])}
	if(sia.data$agemax[i] > 2 && sia.data$agemax[i] <=8) {sia.effective[rownumber,colnumber]<-age_u5[ind]}
	if(sia.data$agemax[i] > 8 && sia.data$agemax[i] <= 13) {sia.effective[rownumber,colnumber]<-age_u10[ind]}
	if(sia.data$agemax[i] > 13) {sia.effective[rownumber,colnumber]<-age_u15[ind]}
}
##########################################################################################
# Post-hoc changes to SIA coverage
##########################################################################################


sia.new<-sia.mat.old*sia.mat
write.table(sia.coverage,file="output-sia_coverage.csv",sep=",",row.names=rownames(sia.mat.old),col.names=1980:(1981+dim(sia.coverage)[2]))
write.table(sia.effective,file="output-sia_effective.csv",sep=",",row.names=rownames(sia.mat.old),col.names=1980:(1981+dim(sia.coverage)[2]))
write.table(sia.new,file="output-sia_age_corrected.csv",sep=",",row.names=rownames(sia.mat.old),col.names=1980:(1981+dim(sia.coverage)[2]))


