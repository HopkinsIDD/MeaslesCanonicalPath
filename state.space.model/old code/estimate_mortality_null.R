estimate_mortality_null<-function(binary.file="output-mortality",file.ext="-null",incidence_workspace="output-incidence-null.Rdat",agedist,agedist.se,cfr1,cfr2,cfr3,mcv_5,codes,ignore.indices,ITER=10000,lb=0.025,ub=0.975,display=T){
#browser()
# should do some error checking here on the size of the inputs
require(MASS)
load(incidence_workspace)
n.countries<-dim(CASES)[1]
n.years<-dim(CASES)[2]
regions<-codes[,1]
codes_cov<-agedist[,1:4]
cty_gbd<-codes[,2]
mcv_5<-mcv_5
mcv_5[is.na(mcv_5)]<-0
#browser()
cfr<-array(0,c(n.countries,n.years,4))
cfr[,,1]<-as.matrix(cfr1[,1:n.years],n.countries,.n.years)
cfr[,,2]<-as.matrix(cfr2[,1:n.years],n.countries,.n.years)
cfr[,,3]<-as.matrix(cfr3[,1:n.years],n.countries,.n.years)

agedist2<-agedist[,5:8]
agedist.se<-agedist.se[,5:8]

measles_deaths<-matrix(NA,n.countries,n.years)
lb_measles_deaths<-matrix(NA,n.countries,n.years)
ub_measles_deaths<-matrix(NA,n.countries,n.years)
ageclass1.mortality<-matrix(NA,n.countries,n.years)
ageclass2.mortality<-matrix(NA,n.countries,n.years)
ageclass3.mortality<-matrix(NA,n.countries,n.years)
ageclass4.mortality<-matrix(NA,n.countries,n.years)
for(i in 1:n.countries){
if(is.na(match(i,ignore.indices))){

	for(j in 1:n.years){
		if(mcv_5[i,(j+1)] <	.60){t<-1}
		if(mcv_5[i,(j+1)] >= .60 & mcv_5[i,(j+1)] < .85){t<-2}
		if(mcv_5[i,(j+1)] >= .85){t<-3}
		tmp<-paste(cty_gbd[i],as.character(t),sep="")
		ind<-which(codes_cov[,2]==tmp)
		#make random draws
		
		p.mean<-unlist(agedist2[ind,])
		p.var<-(unlist(agedist.se[ind,]))^2
		
		age.dist.t<-rep(-1,4)
		age.dist.t<-mvrnorm(ITER,p.mean,diag(p.var))
		age.dist.t<-pmax(age.dist.t,1e-3)
		age.dist.t<-age.dist.t/matrix(apply(age.dist.t,1,sum),nr=ITER,nc=4,byrow=F)

		boot.point<-apply((age.dist.t*matrix(CASES[i,j],ITER,4)*matrix(unlist(cfr[i,j,]),nr=ITER,nc=4,byrow=T)),1,sum)
		ageclass1.mortality[i,j]<-(p.mean[1]*(CASES[i,j])*cfr[i,j,1])
		ageclass2.mortality[i,j]<-(p.mean[2]*(CASES[i,j])*cfr[i,j,2])
		ageclass3.mortality[i,j]<-(p.mean[3]*(CASES[i,j])*cfr[i,j,3])
		ageclass4.mortality[i,j]<-(p.mean[4]*(CASES[i,j])*cfr[i,j,4])
		boot.lower<-apply(age.dist.t*matrix(CIL[i,j],ITER,4)*matrix(unlist(cfr[i,j,]),nr=ITER,nc=4,byrow=T),1,sum)
		boot.upper<-apply(age.dist.t*matrix(CIU[i,j],ITER,4)*matrix(unlist(cfr[i,j,]),nr=ITER,nc=4,byrow=T),1,sum)	
		
		#if(!is.na(median(sort(boot.point))) && median(sort(boot.point))>1e10){browser()}
		
		measles_deaths[i,j]<-median(sort(boot.point))
		lb_measles_deaths[i,j]<-quantile(sort(boot.lower),probs=lb)
		ub_measles_deaths[i,j]<-quantile(sort(boot.upper),probs=ub)
	}
if(display==T){cat("country-",i,"-",rownames(codes)[i],".\n")}                                                 # display progress
}
}

#set AMRO countries to 0 mortality
#AMRO.ind<-which(regions=="AMR")
#measles_deaths[AMRO.ind,]<-0
#lb_measles_deaths[AMRO.ind,]<-0
#ub_measles_deaths[AMRO.ind,]<-0

#set EURO countries to 0 mortality
#EURO.ind<-which(regions=="EUR")
#measles_deaths[EURO.ind,]<-0
#lb_measles_deaths[EURO.ind,]<-0
#ub_measles_deaths[EURO.ind,]<-0


save(measles_deaths,lb_measles_deaths,ub_measles_deaths,ageclass1.mortality,ageclass2.mortality,ageclass3.mortality,ageclass4.mortality,file=paste(binary.file,file.ext,".Rdat",sep=""))

names(ageclass1.mortality)<-1981:(1980+dim(ageclass1.mortality)[2])
names(ageclass2.mortality)<-1981:(1980+dim(ageclass1.mortality)[2])
names(ageclass3.mortality)<-1981:(1980+dim(ageclass1.mortality)[2])
names(ageclass4.mortality)<-1981:(1980+dim(ageclass1.mortality)[2])

names(measles_deaths)<-1981:(1980+dim(ageclass1.mortality)[2])
names(lb_measles_deaths)<-1981:(1980+dim(ageclass1.mortality)[2])
names(ub_measles_deaths)<-1981:(1980+dim(ageclass1.mortality)[2])

ageclass1.mortality<-data.frame(ageclass1.mortality)
ageclass2.mortality<-data.frame(ageclass2.mortality)
ageclass3.mortality<-data.frame(ageclass3.mortality)
ageclass4.mortality<-data.frame(ageclass4.mortality)
measles_deaths<-data.frame(measles_deaths)
lb_measles_deaths<-data.frame(lb_measles_deaths)
ub_measles_deaths<-data.frame(ub_measles_deaths)

write.csv(ageclass1.mortality,file="output-ageclass1-mortality.csv",row.name=rownames(codes))
write.csv(ageclass2.mortality,file="output-ageclass2-mortality.csv",row.name=rownames(codes))
write.csv(ageclass3.mortality,file="output-ageclass3-mortality.csv",row.name=rownames(codes))
write.csv(ageclass4.mortality,file="output-ageclass4-mortality.csv",row.name=rownames(codes))

write.csv(measles_deaths,file="output-estimated_mortality.csv",row.name=rownames(codes))
write.csv(lb_measles_deaths,file="output-lb_estimated_mortality.csv",row.name=rownames(codes))
write.csv(ub_measles_deaths,file="output-ub_estimated_mortality.csv",row.name=rownames(codes))

}