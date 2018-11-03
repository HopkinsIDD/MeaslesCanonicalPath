#-- filtered_trans.R
#-- 
#
#--  Created by Matthew Ferrari on 8/24/09.
#--  Copyright 2009 __MyCompanyName__. All rights reserved.


filt_cases<-function(xf,par,vecs){

	theta2<-par[1]             # original theta2 is the 1st parameter here
	beta2<-exp(theta2)         # project theta2 onto [0,inf)
	
	pop<-vecs$pop              # extract population data
	
	S<-as.vector(xf)           # convert xf from array to vector
	cases<-((1-exp(-beta2*S/pop[-1])))*S
                               #compute cases
	return(cases)              # return cases
}


filt_cases_var<-function(xf,par,vecs,Pf){

	theta2<-par[1]             # original theta2 is the 1st parameter here
	beta2<-exp(theta2)         # project theta2 onto [0,inf)
	pop<-vecs$pop              # extract population data
	S<-as.vector(xf)           # convert xf from array to vector
	d_cases<-( 1 - exp(-beta2*S/pop[-1]) + S*beta2/pop[-1]*exp(-beta2*S/pop[-1]) )             
	                           # compute deviation of cases

	d_cases2<-d_cases^2        # variance of cases
	cases_var<-d_cases2*Pf     # adjust variance with future errors

	return(cases_var)          # return the variance of cases
}


filt_obs_var<-function(xf,par,vecs,Pf){

	theta2<-par[1]            # original theta2 is the 1st parameter here
	theta4<-par[3]            # original theta4 is the 3rd parameter here   
	beta2<-exp(theta2)        # project theta2 onto [0,inf)
	beta4<-1/(exp(-theta4)+1) # project theta4 onto [0,1)
	pop<-vecs$pop             # extract population data
	
	S<-as.vector(xf)          # compute susceptible
	d_cases<-beta4*( 1 - exp(-beta2*S/pop[-1]) + 
	S*beta2/pop[-1]*exp(-beta2*S/pop[-1]) )
	                          # compute actual observed cases deviation
	d_cases2<-d_cases^2       # variance of observed cases

	cases_var<-d_cases2*Pf    # adjust variance with future errors

	return(cases_var)         # return the variance of cases
}
