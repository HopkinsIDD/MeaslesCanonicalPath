#read.outbreak.R

outbreak.data<-read.csv("input-outbreak_list.csv",header=T)
cases<-read.csv("input-cases.csv",header=T,row.names=1)

country.codes<-rownames(cases)
outbreaks<-matrix(1,nr=dim(cases)[1], nc=dim(cases)[2])
for(i in 1:dim(outbreak.data)[1]){
#browser()
	rownumber<-which(country.codes==outbreak.data$iso[i])
	colnumber<-outbreak.data$year[i]-1980+1
	outbreaks[rownumber,colnumber]<-NA
}
##########################################################################################
# Post-hoc changes to outbreaks coverage
##########################################################################################

write.table(outbreaks,file="output-outbreaks.csv",sep=",",row.names=rownames(cases),col.names=1980:(1979+dim(sia.coverage)[2]))
