simulate_cases<-function(par,vecs,ITER=1){
	#browser()
	#need num, y, mu0, Sigma0 and input in the right format
	k<-1 # dummy variable                              
	for(i in 1:length(vecs$parameters)){ # from theta1-5
		if(is.na(vecs$parameters[i])){assign(paste("beta",i,sep=""),par[k])
			k<-k+1}                                #
		else{assign(paste("beta",i,sep=""),vecs$parameters[i])}
		}
		
	births<-vecs$births # births
	pop<-vecs$pop # population size
	mcv1<-vecs$mcv1
	mcv2<-vecs$mcv2
  	sia<-vecs$sia
	X<-as.numeric(births - 0.84*(mcv1)*(1-mcv2)*births - 0.99*(mcv1)*mcv2*births - 0.925*(1-mcv1)*(mcv2)*births)
  	
  	mu0=pop[1]*beta2 # initial susceptible population size
  	Sigma0=1 # 
  	pop=pop # remove first population size data
  	X=X # remove first birth data
        # the data are removed because we need an initial condition to start the model 
	input=matrix(c(X,pop,sia),length(X),3) # prepare input data
	
#browser()  	 	
	Q=beta5 # process variance
 	R=beta6 # observation variance

 # definition of function f - the process sub-model (see Report)
 	f=function(s,x){                                                        
 		retval=(s+x[1]-((1-exp(-beta1*(s/x[2]-1e-2))))*s)*(1-x[3]) # system equation
 		
 		retval # return value
 	}  
 	g=function(s,x){    
 		retval=((1-exp(-as.numeric(beta1)*(s/x[2]-1e-2)))*s)
 		retval # return value
 	} 
 	#browser()
 	T<-dim(input)[1]
 	iter<-ITER
 	I<-matrix(NA,iter,T)
 	S<-matrix(NA,iter,T)
 	S[,1]<-as.numeric(input[1,2]*beta2)
 	for(t in 2:T){
 		S[,t]<- pmax(0,rnorm(iter, mean=as.numeric(f( S[,t-1], input[t,])), sd=as.numeric(ifelse(iter==1,0,sqrt(Q))) ))
 		I[,t]<- g( S[,t-1], as.numeric(input[t,]))
 		#if(median(I[,t])>1e10){browser()}
 	}
	I[!is.finite(I)]<-NA
	return(apply(I,2,quantile,prob=c(.25,.5,.75),na.rm=T))
}