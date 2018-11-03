incidence_estimate<-function(binary.file="output-incidence.Rdat",allcases,allbirths,
                             alldeaths,allpop,allmcv1,mcv1.12mo,allmcv2,allsia,
                             first.surveillance,outbreaks,ignore.indices,display=FALSE){
# test that all the input files are the same size
tst<-cbind(dim(allcases),dim(allbirths),dim(alldeaths),dim(allpop),dim(allmcv1),dim(allmcv2),dim(allsia),dim(outbreaks))
test<-tst!=tst[,1]
ifelse(sum(test)==0, print("OK: input files are the same size "),stop("Error: dimensions of input files not the same"))

n.countries<-dim(allcases)[1]
n.years<-dim(allcases)[2]

# test that all the input files are in the same order
tst<-cbind(rownames(allcases),rownames(allbirths),rownames(alldeaths),rownames(allpop),rownames(allmcv1),rownames(allmcv2),
           rownames(allsia),rownames(outbreaks))
test<-apply(tst,1,unique)
ifelse(length(unlist(test))>dim(tst)[1], stop("Error: ISO codes do not match in all input files"),print("OK: ISO codes match in all input files"))
if(length(unlist(test))==dim(tst)[1]){ISOcodes<-rownames(allcases)}


PAR<-matrix(NA,nr=n.countries,nc=6)                                                             # storage for parameters    
CASES<-matrix(NA,nr=n.countries,nc=n.years)                                                      # storage for incidence
SUS <- matrix(NA,nr=n.countries,nc=n.years)                                                      # storage for susceptibles
CIL<-matrix(NA,nr=n.countries,nc=n.years)                                                             # storage for lower bounds
CIU<-matrix(NA,nr=n.countries,nc=n.years)                                                             # storage for upper bounds
LIKE<-numeric(n.countries)                                                           		    # storage for likelihoods

 for(w in 1:n.countries){
#######################################################
#Calculate for Method 2 countries
#######################################################
 if(is.na(match(w,ignore.indices))){
allcases.w<-allcases[w,first.surveillance[w]:n.years]
allbirths.w<-allbirths[w,first.surveillance[w]:n.years]
alldeaths.w<-alldeaths[w,first.surveillance[w]:n.years]
allpop.w<-allpop[w,first.surveillance[w]:n.years]
allmcv1.w<-allmcv1[w,first.surveillance[w]:n.years]
mcv1.12mo.w<-mcv1.12mo[w,1]
allmcv2.w<-allmcv2[w,first.surveillance[w]:n.years]
allsia.w<-allsia[w,first.surveillance[w]:n.years]
outbreaks.w<-outbreaks[w,first.surveillance[w]:n.years]

T1<-first.surveillance[w]+1	
T2<-dim(CASES)[2]
vecs<-datainput(allbirths.w,
                allcases.w,
                alldeaths.w,
                allpop.w,
                allsia.w,
                allmcv1.w,
                mcv1.12mo.w,
                allmcv2.w,
                outbreaks.w)     # prepare input data file for a country

lower<-c(-3,-5,-10,-10,0,0)                                              # lower bound for parameter 1-5
upper<-c(5,-3,0,0,60,60)                                               # uper bound for parameter 1-5

init<-c(.5,-2,-2,-2,30,30)                                              # initial value for parameter 1-5

vecs$parameters<-c(NA,NA,NA,NA,5,5)                                     # hold parameter 4 and 5 constant

fit.par<-ekf.optim(vecs=vecs,init=init,lower=lower,upper=upper)      # estimate parameter 1-3
par.est<-fit.par$par                                                 # tentative estimates

v1<-seq(5,30,by=.1)                 				     # fit parameter4 using a search between 5 and 30
lik<-numeric(length(v1))                                             # length of search
for(k in 1:length(v1)){                                              # start the search
vecs$parameters<-c(par.est,v1[k],5)				     # hold parameters 1-3 constant and estimate 4
pred1<-S.fx(c(par.est,5,5),vecs,1,1)                                 # estimate true cases
lik[k]<-pred1$like                                                   # update likelihood
}
fitp<-v1[which.min(lik)]                                             # estimate for process variance (parameter 5)

vecs$parameters<-c(par.est,fitp,5)				     # hold parameters 1-3 constant and estimate 5
fitm<-log(S.fx(c(par.est,fitp,5),vecs,1,1)$varreturns)               # estimate for process variance (parameter 5)

vecs$parameters<-c(NA,NA,NA,NA,fitp,fitm)                               # hold parameters 4-5 constant and estimate 1-3
init<-c(par.est[1:4],fitp,fitm)                                      # use estimate for process/measurement variance as initial value
fit.par<-ekf.optim(vecs=vecs,init=init,lower=lower,upper=upper)      # search for all five parameters
 
par.est<-c(fit.par$par,fitp,fitm)                                    # estimates for all five parameters for one country
# End of the main program
##################################
# Here starts post-processing
PAR[w,]<-par.est                                                  	# put estimates into parameter list
pred1<-S.fx(par.est,vecs,1,1)                                        # estimate true cases
cases<-filt_cases(as.vector(pred1$xf),par.est,vecs=vecs)             # perform Kalman smoothing
CASES[w,T1:T2]<-cases	                                               	# put true cases into cases list
S.sd<-sqrt(as.vector(pred1$Pf))                                      # compute true cases standard deviation
cases.lb<- filt_cases(pmax(1,as.vector(pred1$xf)-2*S.sd),par.est,vecs=vecs) 
                                                                     # compute lower bounds of cases
cases.ub<- filt_cases(as.vector(pred1$xf)+2*S.sd,par.est,vecs=vecs)
                                                                     # compute upper bounds of cases
SUS[w,T1:T2] <- pred1$xf            # store mean estimate for susceptibles by year
CIL[w,T1:T2]<-cases.lb				# store the lower bound
CIU[w,T1:T2]<-cases.ub			# store the upper bound

like<-pred1$like                                                     # likelihood
LIKE[w]<-like                                                   # store likelihood

#if estimate is less then observed, replace with observed
if(any(allcases.w>CASES[w,],na.rm=T)){CASES[w,(which(allcases.w>CASES[w,]))]<-unlist(allcases.w[which(allcases.w>CASES[w,])])} 
#if upper bound is greater than population size, replace with population size
if(any(allpop.w<CIU[w,],na.rm=T)){CIU[w,(which(allpop.w>CIU[w,]))]<-unlist(allpop.w[which(allpop.w>CIU[w,])])} 
if(display==T){cat("country-",w,"-",ISOcodes[w],".\n")}                                                 # display progress
}
#######################################################
#Fill in for Method 1 countries
#######################################################
if(!is.na(match(w,ignore.indices))){
allcases.w<-unlist(allcases[w,])
CASES[w,]<-allcases.w/.2		# Assume that reporting rate is 20%
CIU[w,]<-allcases.w/.05			# Calculate upper bound by assuming that reporting rate is 5%
CIL[w,]<-allcases.w/.4			# Calculate upper bound by assuming that reporting rate is 40% 

}
}
#CASES[which(allcases>CASES,arr.ind=T)]<-cases[which(allcases>CASES,arr.ind=T)]

CASES<-CASES[,2:n.years]
CIL<-CIL[,2:n.years]
CIU<-CIU[,2:n.years]
SUS<-SUS[,2:n.years]

#untransform the parameters
PAR[,c(1,5,6)]<-exp(PAR[,c(1,5,6)])
PAR[,2]<-1/(exp(-PAR[,2])+1)
PAR[,3]<-1/(exp(-PAR[,3])+1)
PAR[,4]<-PAR[,3]+1/(exp(-PAR[,4])+1)

PAR<-data.frame(PAR)
CASES<-data.frame(CASES)
CIL<-data.frame(CIL)
CIU<-data.frame(CIU)
LIKE<-data.frame(LIKE)
SUS <- data.frame(SUS)

names(PAR)<-c("transmission","initial susceptibles","reporting","outbreak reporting","transmission variance","reporting variance")
names(CASES)<-1981:(1980+dim(CASES)[2])
names(CIL)<-1981:(1980+dim(CASES)[2])
names(CIU)<-1981:(1980+dim(CASES)[2])
names(SUS)<-1981:(1980+dim(SUS)[2])
names(LIKE)<-"likelihood"

write.csv(PAR,file="output-pars.csv",row.name=ISOcodes)
write.csv(CASES,file="output-estimated_incidence.csv",row.name=ISOcodes)
write.csv(CIL,file="output-lb_incidence.csv",row.names=ISOcodes)
write.csv(CIU,file="output-ub_incidence.csv",row.names=ISOcodes)
write.csv(LIKE,file="output-like.csv",row.name=ISOcodes)
write.csv(SUS,file="output-sus",row.name=ISOcodes)


save(PAR,CASES,CIL,CIU,LIKE,SUS,file=binary.file)                                   # save outputs in a R workspace for future loading and analysis
}