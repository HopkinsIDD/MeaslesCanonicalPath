simulate_incidence<-function(binary.file="output-incidence",file.ext="-null",PAR,allcases,allbirths,alldeaths,allpop,allmcv1,mcv1.12mo,allmcv2,allsia,ignore.indices,display=FALSE,ITER=1){
#browser()
# test that all the input files are the same size
tst<-cbind(dim(allcases),dim(allbirths),dim(alldeaths),dim(allpop),dim(allmcv1),dim(allmcv2),dim(allsia))
test<-tst!=tst[,1]
ifelse(sum(test)==0, print("OK: input files are the same size "),stop("Error: dimensions of input files not the same"))

n.countries<-dim(allcases)[1]
n.years<-dim(allcases)[2]

# test that all the input files are in the same order
tst<-cbind(rownames(allcases),rownames(allbirths),rownames(alldeaths),rownames(allpop),rownames(allmcv1),rownames(allmcv2),rownames(allsia))
test<-apply(tst,1,unique)
ifelse(length(unlist(test))>dim(tst)[1], stop("Error: ISO codes do not match in all input files"),print("OK: ISO codes match in all input files"))
if(length(unlist(test))==dim(tst)[1]){ISOcodes<-rownames(allcases)}

CASES<-matrix(NA,nr=n.countries,nc=n.years)                                                      # storage for incidence
CIL<-matrix(NA,nr=n.countries,nc=n.years)                                                      # storage for incidence
CIU<-matrix(NA,nr=n.countries,nc=n.years)                                                      # storage for incidence

for(w in 1:n.countries){
 if(is.na(match(w,ignore.indices))){
allcases.w<-allcases[w,]
allbirths.w<-allbirths[w,]
alldeaths.w<-alldeaths[w,]
allpop.w<-allpop[w,]
allmcv1.w<-allmcv1[w,]
mcv1.12mo.w<-mcv1.12mo[w,1]
allmcv2.w<-allmcv2[w,]
allsia.w<-allsia[w,]
par.w<-PAR[w,]
outbreaks.w<-NA

#T1<-first.surveillance[w]+1	
#T2<-dim(CASES.sim)[2]
vecs<-datainput(allbirths.w,allcases.w,alldeaths.w,allpop.w,allsia.w,allmcv1.w,mcv1.12mo.w,allmcv2.w,outbreaks.w)     # prepare input data file for a country
vecs$parameters<-par.w
cases<-simulate_cases(par.w,vecs,ITER=ITER)

CASES[w,]<-cases[2,]
CIL[w,]<-cases[1,]
CIU[w,]<-cases[3,]
if(display==T){cat("country-",w,"-",ISOcodes[w],".\n")}                                                 # display progress
}
}
CASES<-CASES[,2:n.years]
CIL<-CIL[,2:n.years]
CIU<-CIU[,2:n.years]
#browser()
write.table(CASES,file=paste("output-estimated-incidence",file.ext,".csv",sep=""),sep=",",row.name=ISOcodes,col.name=F)
file=paste("ageclass1-mortality",file.ext,".csv",sep="")
save(CASES,CIL,CIU,file=paste(binary.file,file.ext,".Rdat",sep=""))                                   # save outputs in a R workspace for future loading and analysis
}