ekf.optim<-function(vecs,init,lower,upper, ...){  # Extended Kalman filter optimzation function
	ind<-which(is.na(vecs$parameters)) # find out which parameters to optimize
	fit<-optim(init[ind],S.fx,method="L-BFGS-B",
	lower=lower[ind],upper=upper[ind],vecs=vecs,likesig=1, 
	hessian=T)
# use L-BFGS-B method (bounded quasi-Newton optimization) with Hessian matrix
	return(fit)  # return fit value
}

ekf.optim.NM<-function(vecs,init,lower,upper, ...){ # Extended Kalman filter optimzation function 2
	ind<-which(is.na(vecs$parameters))          # find out which parameter to optimize
	fit<-optim(init[ind],S.fx,method="Nelder-Mead",vecs=vecs,likesig=1)
                                                    # use heuristic Nelder-Mead method 
	return(fit)                                 # return fit value
}



ekf.optim_ln<-function(vecs,init,lower,upper, ...){ # Extended Kalman filter optimzation function 3 
	ind<-which(is.na(vecs$parameters))          # find out which parameter to optimize
	fit<-optim(init[ind],S.fx_ln,method="L-BFGS-B",lower=lower[ind],upper=upper[ind],vecs=vecs,likesig=1)
                                                    # use L-BFGS-B method (bounded quasi-Newton optimization) without Hessian matrix
	return(fit)                                 # return fit value
}

