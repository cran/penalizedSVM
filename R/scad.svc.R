scadsvc <- function(lambda = 0.01, x, y, a = 3.7, tol = 10^(-4), class.weights= NULL){
	# SCAD SVM classification
	#
	# Input:
	#   xtrain : n-by-d data matrix to train (n chips/patients, d clones/genes, d>>n )
	#   ytrain : column vector of target {-1, 1}'s (for n chips/patiens )
	#   a : tuning parameter in scad function (default: 3.7 or whatever the paper
	#uses)
	#		lambda : tuning parameter in scad function (default : 2 or whatever the
	#paper uses)
	#		tol: the cut-off value to be taken as 0
	#
	# By Axel Benner (15.03.2007)
	# 
	
	# SVM mit variable selection (clone selection) using scad penalty.
	

	# checks
	if(nrow(x) != length(y)) stop("Wrong dimensions: nrow(x) should be equal to length(y) !")
	if (nlevels(as.factor(y)) !=2 ) stop(paste("We need 2 classes, currently have:",  paste(levels(as.factor(y)), collapse=", ")) ) 
	
	xtrain <- x
	# change class labels to -1 and 1
	ytrain <- 2*as.numeric(as.factor(y))-3 

# start with linear svm:
	require(e1071, quietly=TRUE)
	require(corpcor) # for pseudoinverse
	require(statmod) # for matrixmultplications 
	
	linsvc <- svm(xtrain, factor(ytrain), kernel="linear", class.weights=class.weights, fitted =FALSE)
			
	# index of the support vectors
	index <- linsvc$index
	# type of svm
	type <- linsvc$type
	# w: coefficients times the training labels * SV
	w = apply(as.vector(linsvc$coefs) * linsvc$SV, 2, sum)  
	
	# rho = Bias: A scalar value representing the bias or threshold of the SVM
	#classifier,  which is the negative intercept.
	b = linsvc$rho
	diff = 1
	ntrain = nrow(xtrain)
	d = ncol(xtrain)   
	xind = 1:d
	i<-1
	print("start iterations:")
	while (diff > tol) {
		# cat(paste(diff, " ",sep=""))
		
		x1 = cbind(rep(1,nrow(xtrain)), xtrain)
		sgnres = as.vector(ytrain - x1 %*% c(b, w))
		# important!!!! : sometimes a point is lying exactly on the hyperline --> sgnres = 0
		# produce errors --> move this randomly at the one or the other size.
		sgnres[sgnres== 0] <- sample(c(1,-1),1) *  10^(-100)  
		res = abs(sgnres)
		y0 = ytrain / res
		
		# ###
		#D = 1/(2*ntrain) * diag(1/res)
		# ###
		# save as a vector D_vec
		D_vec = 1/(2*nrow(xtrain)) * (1/res)
		aw = abs(w)
		dp = lambda*(aw<=lambda)+(a*lambda-aw)/(a-1) * (lambda<aw &  aw <=a*lambda)
		Q1_vec<-c(0, dp/aw)
		
		P = 0.5 * t(ytrain + y0) %*% x1 / ntrain
		
		# ###
		#Q1 = diag( c(0, dp/aw))
		#Q = t(x1) %*% D %*% x1 + Q1
		# ps_Q<-pseudoinverse(Q)
		# ###
		
		# inv_Q is sometimes too large 
		#inv_Q<-.find.inverse(U=t(x1),D_vec=D_vec,A_vec=Q1_vec)
		# nwb = inv_Q %*% t(P)
		
		#nwb = inv( t(x1) %*% D %*% x1 + Q1) %*% t(P)
		# --> 
		nwb<-.calc.mult.inv_Q_mat2(U=t(x1),D_vec=D_vec, A_vec=Q1_vec, mat2=t(P) )
		#summary(nwb)
		
		## test # ok
		
		#Q = U %*% diag(D_vec) %*% t(U) + diag(Q1_vec)
		#inv_Q.0<-pseudoinverse(Q)
		#nbw.0<-inv_Q.0 %*% t(P)
		#summary(nbw.0)
		
		#nwb<-.calc.mult.inv_Q_mat2(U=t(x1),D_vec=D_vec, A_vec=Q1_vec, mat2=t(P), n.thr=1000)
		#summary(nwb)
		
		#inv_Q<-.find.inverse(U=t(x1),D_vec=D_vec,A_vec=Q1_vec )
		# nwb.inv= inv_Q %*% t(P)
		#summary(nwb.inv)
			
		## end of test 
		
		
		nw = nwb[-1]
		nb = nwb[1]
		diff = sqrt( sum( (nwb - c(b,w))^2 ) )
		ind = abs(nw) > 0.001  # oracle threshold
		#table(ind)
		#print(diff)
		if (sum(!is.na(ind))>0  & sum(ind)>0) {
		w = nw[ind]
		xtrain = xtrain[, ind, drop=FALSE]
		xind = xind[ind]
		b = nb
		
		}	 else {
			diff=tol/2
		}
		i<-i+1
		
	}
	print(paste("scad converged in", i, "iterations" )) 
	#

	ind = abs(nw) > 0.001 # oracle threshold
	#print("ind  ,   if (sum(ind)>0)")
	#print(ind)
	
	if (sum(!is.na(ind))>0  & sum(ind)>0){
			w = nw[ind]
			names(w)=colnames(xtrain)[ind]
			b = nb
			f = as.vector(xtrain %*% w + b)
			
			# xqx = 0.5 * x1 %*% inv_Q %*% t(x1)  - we don't have a huge quadratic matrix inv_Q, use the same trick as for nwb 
			qx<- .calc.mult.inv_Q_mat2(U=t(x1),D_vec=D_vec, A_vec=Q1_vec, mat2=t(x1) )
			xqx =  0.5 * x1 %*% qx
						
			# Output:
			#   w : direction vector of separating hyperplane
			#   b : the 'bias'
			#   xind : Indices of remained variables
			#   fitted : Fit of regularized SVM (for all patients with reduced set of genes )
			#ret <- list(w=w, b=b, xind=xind, index=index, xqx=xqx, fitted=f, type=type, x=x, y=y)
			ret <- list(w=w, b=b, xind=xind, index=index, xqx=xqx, fitted=f, type=type)
			class(ret) <- "scadsvm"
			return(ret)
		
	} else {
	return("No variable selected.")
	}
}
