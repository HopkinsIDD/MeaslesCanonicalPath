S.fx<-function(par,vecs,likesig,smoother=0){
	# main function of EKF
	#need num, y, mu0, Sigma0 and input in the right format
	k<-1 # dummy variable                              
	for(i in 1:length(vecs$parameters)){ # from theta1-5
		if(is.na(vecs$parameters[i])){assign(paste("theta",i,sep=""),par[k])
			k<-k+1}                                #
		else{assign(paste("theta",i,sep=""),vecs$parameters[i])}
		}

	beta1<-exp(theta1) # transmission rate 
	beta2<-1/(exp(-theta2)+1) # initial susceptibles
	beta3<-1/(exp(-theta3)+1) # Baseline underreporting rate
	beta4<-1/(exp(-theta4)+1) # underreporting rate
	beta5<-exp(theta5) # system/process variance	
	beta6<-exp(theta6) # observation/measurement variance

	
	X<-vecs$X # births - < 5 deaths - vaccination, potential susceptible pool
	sia<-vecs$sia # pulsed vaccination coverage
	Ir<-vecs$Ir # reported cases
	births<-vecs$births # births
	pop<-vecs$pop # population size
	outbreaks<-vecs$outbreaks
  
  	mu0=pop[1]*beta2 # initial susceptible population size
  	Sigma0=1 # 
  	Ir=Ir[-1] # remove first reported case (yr 1979)
  	outbreaks<-outbreaks[-1]
  	sia=sia[-length(sia)] # remove last population
  	pop=pop[-1] # remove first population size data
  	X=X[-1] # remove first birth-vaccine data
        # the data are removed because we need an initial condition to start the model 
  	num=length(Ir) # length of the data
  	y=Ir # transfer reported cases                                        
  	reporting<-numeric(length(Ir))
  	reporting[which(outbreaks==1)]<-beta3
	reporting[which(is.na(outbreaks))]<-beta3+beta4
	input=matrix(c(X,pop,sia,reporting),length(X),4) # prepare input data
	
#browser()  	 	
	Q=beta5 # process variance
 	R=beta6 # observation variance

 # definition of function f - the process sub-model (see Report)
 	f=function(s,x){                                                        
 		retval=(s+x[1]-((1-exp(-beta1*(s/x[2]-1e-2))))*s)*(1-x[3]) # system equation
 		retval # return value
 	}  
  # derivative of f with respect to s
  	divf=function(s,x){
                # these are for the Hessian Matrix	
		retval = (1-x[3])*(1- (( 1-exp(-beta1*s/x[2]) ) +  
		          s * (beta1/x[2])*exp(-beta1*s/x[2]))) #derivative of f with respect to s
		retval # return the derivative
	}
 
 # definition of function h - the observation sub-model (see Report)
 	h=function(s,x){                                                        
    	retval=x[4]*((1-exp(-beta1*s/x[2])))*s # observation equation
    	retval # return value
 	}
 
  # derivative of h with respect to s
 	divh=function(s,x){     
                # these are for the Hessian Matrix                                               
		retval=(x[4]*( 1-exp(-beta1*s/x[2]) ) + 
		        x[4]* s * (beta1/x[2])*exp(-beta1*s/x[2]))
                #derivative of h with respect to s            
    	retval # return value                   
 	}
 
 pdim=1 # nrow(Phi)                
 
 if (is.na(vecs$parameters[5])){ # smooth for theta4
 	for (j in 1:3){ # 3 iterations  
 	temp=filtersmooth(y,input,num,pdim,Q,R,f,
 	      divf,h,divh,mu0,Sigma0,1,1)    # call for the smooth function   
 	R=temp$varreturns   # return from the smooth function
 	}
 }
 
 returns=filtersmooth(y,input,num,pdim,Q,R,f,
         divf,h,divh,mu0,Sigma0,likesig,smoother) # call for the smoother
 returns   # return from the smoother

}
