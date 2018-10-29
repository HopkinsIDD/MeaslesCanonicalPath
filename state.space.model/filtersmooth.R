filtersmooth = function(y,input,num,pdim,Q,R,f,divf,h,divh,mu0,Sigma0,likesig,smoother)  {
	 y = as.matrix(y)                            # y in matrix form
   qdim = ncol(y)                                  # dimension of q is the column length of y matrix
   rdim = ncol(as.matrix(input))                   # dimension of r is column length of input matrix
   ut = matrix(input,num,rdim)                     # combine into a new matrix
   xp = array(NA, dim = c(pdim,1,num))               # xp = x_t^{t-1}, xp is present estimate       
   Pp = array(NA, dim = c(pdim,pdim,num))            # Pp = P_t^{t-1}, Pp is present error covariance
   xf = array(NA, dim = c(pdim,1,num))               # xf = x_t^t, xf is future estimate
   Pf = array(NA, dim = c(pdim,pdim,num))            # Pf = x_t^t, Pf is future error covariance
   innov = array(NA, dim = c(qdim,1,num))            # initialize innovation
   sig = array(NA, dim = c(qdim,qdim,num))           # innovation covriance matrix
                                               # initialize (because R can't count from zero) 
   x00 = as.matrix(mu0, nrow = pdim, ncol = 1)       # initial estimate
   P00 = as.matrix(Sigma0, nrow = pdim, ncol = pdim) # initial error covariance matrix
   xp[,,1] = f(x00, ut[1,])                      # compute present estimate
   Phi = as.matrix(divf(x00,ut[1,]))              # extrinsic input
   Pp[,,1] = Phi %*% P00 %*% t(Phi) + Q                # compute present error covariance matrix
   B  =  matrix(divh(xp[,,1],ut[1,]), nrow = qdim, ncol = pdim) 
                                               # derivative of current estimates
   sigtemp = B %*% Pp[,,1] %*% t(B) + R
   sig[,,1] = (t(sigtemp) + sigtemp)/2             # innovation matrix - make sure it's symmetric
   siginv = solve(matrix(sig[,,1]))              # inverse of innovation matrix            
   K = Pp[,,1] %*% t(B) %*% siginv                   # compute Kalman filter gain
   innov[,,1] = y[1,]-h(xp[,,1],ut[1,])           # update innovation
   xf[,,1] = xp[,,1] + K %*% innov[,,1]              # update future estimate
   Pf[,,1] = Pp[,,1]-K %*% B %*% Pp[,,1]             # update future error covariance
   sigmat = as.matrix(sig[,,1], nrow = qdim, ncol = qdim)                # rearrange innovation
   like  =  log(det(sigmat))  +  t(innov[,,1]) %*% siginv %*% innov[,,1]   # -log(likelihood)
   hreturn = h(xp[,,1],ut[1,])                   # compute observation
   #browser()
############################# 
# start filter iterations
#############################
 for (i in 2:num){                                          # start iteration loop                       
	 yval = y[i,]                                         # transfer y value
	 xp[,,i] = f(xf[,,i-1], ut[i,])                       # compute current estimate
	 Phi = as.matrix(divf(xf[,,i-1],ut[i,]))              # extrinsic input
	 Pp[,,i] = Phi %*% Pf[,,i-1] %*% t(Phi) + Q                 # compute current covariance
	 B  =  matrix(divh(xp[,,i],ut[i,]), nrow = qdim, ncol = pdim) # derivative of current estimates
	 if (is.na(yval)) B = matrix(0, nrow = qdim, ncol = pdim) # deal with NA data in y
	 sigtemp = B %*% Pp[,,i] %*% t(B) + R                       # compute alternative innovation matrix
	 if (is.na(yval)) sigtemp[1,1] = 1                    # scale innovation covariance
	 sig[,,i] = (t(sigtemp) + sigtemp)/2                   # make sure innovation is symmetric
	 siginv = solve(matrix(sig[,,i]))                     # now siginv is sig[[i]]^{-1} （innovation inverse）
	 K = Pp[,,i] %*% t(B) %*% siginv                          # Kalman gain                                              
	 innov[,,i] = (yval-h(xp[,,i],ut[i,]))                # compute innovation
	 if (is.na(yval)) innov[,,i] = 0                      # adjust innovation if NA is present
	 xf[,,i] = xp[,,i] + K %*% innov[,,i]                     # compute future estimate
	 hreturn = c(hreturn,h(xp[,,i],ut[i,]))               # return observation
	 Pf[,,i] = Pp[,,i]-K %*% B %*% Pp[,,i]                    # compute future error covariance
	 sigmat = matrix(sig[,,i], nrow = qdim, ncol = qdim)      # innovation matrix
	 like =  like  +  log(det(sigmat))  +  t(innov[,,i]) %*% siginv %*% innov[,,i]   
                                                        # update likelihood	 
  }
#browser()
    like = 0.5*like                         # adjust likelihood by 50% because of log scale
    if (likesig == 0)                       # if likesig is 0
    { valreturn = list(xp = xp,Pp = Pp,xf = xf,Pf = Pf,like = like,innov = innov,sig = sig,Kn = K,ut = ut,hreturn = hreturn) }
                                          #  return list of outputs
    else 
    { valreturn = like }                    # return likelihood
    
    if (smoother == 1)                       # perform smoother
    { 
      xs = array(NA, dim = c(pdim,1,num))      # xs = x_t^n, smoothed estimates
      Ps = array(NA, dim = c(pdim,pdim,num))   # Ps = P_t^n, smoothed error covariance
      J = array(NA, dim = c(pdim,pdim,num))    # J = J_t, joint posterior distribution covariance
      xs[,,num] = xf[,,num]                  # initialize estimates
      Ps[,,num] = Pf[,,num]                  # initialize erro covariance
      
      for(k in num:2)  {
        # smoother loop
        Phi = divf(xf[,,k],ut[k,])           # change here
        J[,,k-1] = (Pf[,,k-1] %*% t(Phi)) %*% solve(Pp[,,k])  
        # smooth joint post. dist. cov.
        xs[,,k-1] = xf[,,k-1] + J[,,k-1] %*% (xs[,,k]-xp[,,k])
        # smooth estimates
        Ps[,,k-1] = Pf[,,k-1] + J[,,k-1] %*% (Ps[,,k]-Pp[,,k]) %*% t(J[,,k-1])
        # smooth error covariance
		}
       # and now for the initial values because R can't count backward to zero
       x00 = mu0                             # first estimate       
       P00 = Sigma0                          # first error covariance
       Phi = divf(xf[,,1],ut[1,])            # derivative of observation 
       J0 = as.matrix((P00 %*% t(Phi)) %*% solve(Pp[,,1]), nrow = pdim, ncol = pdim)
                                           # compute joint posterior distribution covariance
       x0n = as.matrix(x00 + J0 %*% (xs[,,1]-xp[,,1]), nrow = pdim, ncol = 1)
                                           # compute first estimate
       P0n =  P00  +  J0 %*% (Ps[,,k]-Pp[,,k]) %*% t(J0)
                                           # compute first error covariance
       hreturns = c()                        # list returns
		for (i in 1:num){                  
			hreturns = c(hreturns,h(xs[,,i],ut[i,]))} 
			                               # smoothed observations
		  varreturns = 0                       # return 0
       for (i in 1:num){
         
         if (!is.na(y[i,])) {           # deal with missing data
           B  =  matrix(divh(xs[,,i],ut[i,]), nrow = qdim, ncol = pdim) 
           # derivative of current estimates
           varreturns = varreturns + (y[i,]-hreturns[i])^2 + B %*% Ps[,,i] %*% t(B)}
       }
       		                               # compute return value
		  numtemp = length(y[!is.na(y[,1]),1]) 
		  # calculate length of NAs
		  varreturns = varreturns/numtemp  # adjust return values
		  valreturn = list(xs = xs,Ps = Ps,x0n = x0n,P0n = P0n,J0 = J0,J = J,xp = xp,
		                   Pp = Pp,xf = xf,Pf = Pf,like = like,innov = innov,sig = sig,
		                   Kn = K,ut = ut,hreturn = hreturn,hreturns = hreturns,varreturns = varreturns)
		  # list of return
    }                                       #endsmoother
    valreturn                               # function return
}
