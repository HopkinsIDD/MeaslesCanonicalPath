sia.data<-read.csv("/Users/matthewferrari/Documents/Current Projects/WHO/Burden Estimation/Geneva2013/new files/SIAs cleaned with single country year aggregation 2000 to 2011.v2.csv",header=T)
sia.mat.old<-read.csv("/Users/matthewferrari/Documents/Current Projects/WHO/Burden Estimation/Geneva2013/run2/MDB1 -12Dec_Somreportedcases adjusted/sia.csv",header=F,row.names=1)


mcv_5<-read.csv("/Users/matthewferrari/Documents/Current Projects/WHO/Burden Estimation/Geneva2013/run2/MDB1 -12Dec_Somreportedcases adjusted/mcv1x5yrs.csv",header=F,row.names=1)
codes<-read.csv("/Users/matthewferrari/Documents/Current Projects/WHO/Burden Estimation/Geneva2013/run2/MDB1 -12Dec_Somreportedcases adjusted/gbd_Region.csv",header=T,row.names=1)
agedist<-read.csv("/Users/matthewferrari/Documents/Current Projects/WHO/Burden Estimation/Geneva2013/run2/MDB1 -12Dec_Somreportedcases adjusted/agedist_pred.csv",header=T)
codes_cov<-agedist[,1:4]

agedist.sum<-apply(agedist[,5:8],1,cumsum)
age_u1<-agedist.sum[1,]
age_u5<-agedist.sum[2,]
age_u10<-agedist.sum[3,]
age_u15<-agedist.sum[4,]

country.codes<-rownames(sia.mat.old)
cty_gbd<-codes[,2]

sia.mat<-matrix(0,nr=dim(sia.mat.old)[1], nc=dim(sia.mat.old)[2])
for(i in 1:dim(sia.data)[1]){
	rownumber<-which(country.codes==sia.data$iso[i])
	colnumber<-sia.data$year[i]-1980+1
#browser()
	if(mcv_5[rownumber,colnumber] <	.60){t<-1}
	if(mcv_5[rownumber,colnumber] >= .60 & mcv_5[rownumber,colnumber] < .85){t<-2}
	if(mcv_5[rownumber,colnumber] >= .85){t<-3}
	
	tmp<-paste(cty_gbd[rownumber],as.character(t),sep="")
	ind<-which(codes_cov[,2]==tmp)
	
	if(sia.data$agemax[i] <= 2) {sia.mat[rownumber,colnumber]<-.5*(age_u1[ind]+age_u5[ind])}
	if(sia.data$agemax[i] > 2 && sia.data$agemax[i] <=8) {sia.mat[rownumber,colnumber]<-age_u5[ind]}
	if(sia.data$agemax[i] > 8 && sia.data$agemax[i] <= 13) {sia.mat[rownumber,colnumber]<-age_u10[ind]}
	if(sia.data$agemax[i] > 13) {sia.mat[rownumber,colnumber]<-age_u15[ind]}
}
sia.new<-sia.mat.old*sia.mat
write.table(sia.new,file="sia_age_corrected.csv",sep=",",row.names=rownames(sia.mat.old),col.names=F)
###RECHECK THE SIA COVERAGES FOR THE AGGREGATED FILE??


