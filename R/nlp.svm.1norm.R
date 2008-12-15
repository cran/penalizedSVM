lpsvm<-function (A, d, k=5, nu=0,output=1, delta=10^-3, epsi=10^-4,  my.seed=123){
#function [w,gamma,trainCorr, testCorr,  nu]=lpsvm(A,d,k,nu,output,delta)
## this implements L1 SVM classification, actually a fast Newton algorithm NLPSVM from Fung and Mangasarian(2004)
# lpsvm = nlp svm = L1 svm = 1-norm svm
#
# version 1.3
# last revision: 07/07/03
#===========================================================================
# Usage: lpsvm(A,d,k,nu,output,delta);
#
# A and d are both required, everything else has a default
# An example: [w gamma train test time nu] = lpsvm(A,d,10);
#
# Input parameters:
#    A: Data points = inputmatrix in R (m,n ); A row= pat, col=feature
#    d: 1's or -1's length= m 
#    k: k-fold for cv 
#     way to divide the data set into test and training set
#       if k = 0: simply run the algorithm without any correctness
#         calculation, this is the default
#       if k = 1: run the algorithm and calculate correctness on
#         the whole data set
#       if k = any value less than the number of rows in the data set:
#         divide up the data set into test and training
#         using k-fold method
#       if k = number of rows in the data set: use the 'leave 1' method
#		epsi: tuning parameter	
#
#    output: 0 - no output, 1 - produce output, default is 0
#    nu:             weighted parameter
#                    1 - easy estimation
#                    0  - hard estimation 
#                    any other value - used as nu by the algorithm
#                    default - 0
#    delta:  default is 10^-3
#===================================================================
# Output parameters:
#
#       w:              the normal vector of the classifier
#       gamma:          the threshold
#       trainCorr:      training set correctness
#       testCorr:       test set correctness
#       nu:             estimated value (or specified value) of nu
#==========================================================================
	
	require(statmod) # for matrixmultplications 
	require(corpcor)
	
	# be shure that d is a vector of 1 and -1
	tmp<-matrix(as.numeric(as.character(d)))
	rownames(tmp)<-names(d)
	d<-tmp
	#str(d)
	# d is a 1-row matrix: 
	
	if (nu==0){nu = .EstNuLong(A,d)}       # default is hard estimation
	if (nu==1){nu = .EstNuShort(A,d)} 
	# else the use defined nu
	
	set.seed(my.seed)
	r<- sample(1:length(d))
	d<-d[r]
	A <- A[r, ]
	
	trainCorr=0;
	testCorr=0;
	
	# flag: non empty model
	non.empty<-rep(NA, k)
	
	if ( k==0) { 
		#[w, gamma,iter] = .core(A,d,nu,delta);
		# no cv 
		List<-  .core(A=A,d=d,nu=nu,delta=delta, epsi=epsi)#
		if (!is.null(List)){
			non.empty[k]<- TRUE
			w<-  List$w
			gamma <- List$gamma
			iter <- List$iter
			xind<-  List$xind
			
			if (output==1) 	print(paste("Number of Iterations:" ,iter))
		}else {
			non.empty[k]<- FALSE
			if (output==1) print("empty model")
		} # end k=0
	} else{
		#if k==1 only training set correctness is calculated
		if (k==1){
			List = .core(A=A,d=d,nu=nu,delta=delta, epsi=epsi);
			if (!is.null(List)){
				non.empty[k]<- TRUE
				w<-  List$w
				gamma <- List$gamma
				iter <- List$iter
				xind<-  List$xind
				
				trainCorr = round(.correctness(AA=A[,xind, drop=FALSE],dd=d,w=w,gamma=gamma), 1);
				if (output == 1){
					print(paste("Training set correctness: ",trainCorr, "%"));
					print(paste("Number of Iterations:",iter));
				}
			}else {
				non.empty[k]<- FALSE
				if (output==1) print("empty model")
			} # end k==1
		}else{
			# do k-fold CV	
			sm<-nrow(A)
			sn<- ncol(A)
			# accumulative iterations
			accuIter = 0;
			
			num.w<- 0 # average number of genes in the final model
			
			indx = c(0:k);
			indx = floor(sm*indx/k);    #last row numbers for all 'segments'
			# split trainining set from test set
			for  (i in  1:k){
				#for  (i in  c(1,2,3,5) ){
				
				# ist ok to take the first nrow/k rows , because we have already permuted them.
				Ctest = A[(indx[i]+1):indx[i+1],];
				dtest = d[(indx[i]+1):indx[i+1]];
				
				# everything ecxept Ctrain 
				Ctrain = A[- ((indx[i]+1):indx[i+1]),];
				dtrain = d[- ((indx[i]+1):indx[i+1])];
				
				List = .core(Ctrain,dtrain,nu,delta, epsi);
				
				# if the list is not empty
				if (!is.null(List)){
					non.empty[i]<- TRUE
					w<-  List$w
					gamma <- List$gamma
					iter <- List$iter
					xind<-  List$xind
					
					tmp.num.w<-length(w)
					tmpTrainCorr = round(.correctness(Ctrain[,xind, drop=FALSE],dtrain,w=w,gamma=gamma),1)
					tmpTestCorr = round(.correctness(Ctest[,xind, drop=FALSE],dtest,w=w,gamma=gamma),1);
					
					if (output == 1){
						print("________________________________________________");
						print(paste("Fold",i));
						print(paste("Training set correctness:",tmpTrainCorr));
						print(paste("Testing set correctness:",tmpTestCorr));
						print(paste("Number of iterations:",iter));
						print("selected genes:")
						print(w)
					
					}
					
					trainCorr = trainCorr + tmpTrainCorr
					testCorr = testCorr + tmpTestCorr
					accuIter = accuIter + iter # accumulative iterations
					num.w<- num.w + tmp.num.w
				}else{
					non.empty[i]<- FALSE
					if (output == 1){
						print("________________________________________________");
						print(paste("Fold",i));
						print("empty model");
					}
				}
			} # end of for (looping through test sets)
				
				# take in accout the number of non-empty Lists (models)!
				trainCorr = trainCorr/sum(non.empty)
				testCorr = testCorr/sum(non.empty)
				
				if (output == 1){
					print("==============================================");
					print(paste("Training set correctness:",trainCorr))
					print(paste("Testing set correctness:" ,testCorr))
					print(paste("Average number of iterations: ",accuIter/sum(non.empty)))
					print(paste("Average number of genes: ",num.w/sum(non.empty)))
				}
			}
	
	}
	ret<- list(w=w, b=gamma, xind=xind, epsi=epsi, iter =iter, trainCorr=trainCorr,testCorr=testCorr, nu=nu )
	class(ret) <- "1norm"
	return(ret)
}
