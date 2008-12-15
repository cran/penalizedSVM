run.1norm<-function(x,y,k=5,nu=0, lambda1.set=NULL, output=1, seed=seed){  
	## find the best epsi via k-fold cv,  get the finla model with optimal epsi 
	#epsi - tuning parameter !
	# fine grid at the begin of the interval !!!
	ff.list<-list()
	
	# try with cv
	# fine grid 10^-5, 2*10^-5, 3*10^-5, ..., 9*10^-5, 10^-4, ..., 10^-1, 0.2, ..., 0.9
	if (is.null(lambda1.set)) {
		epsi.set<-vector(); for (num in (1:9)) epsi.set<-sort(c(epsi.set, c(num*10^seq(-5, -1, 1 ))) )
		lambda1.set<- epsi.set
	}
	# if the internal error occures, skip it, set the model to the empty model and go forther
	# error text: "Error in La.svd(x, nu, nv) : error code 1 from Lapack routine 'dgesdd' "
	catch.error<-function()ifelse ( length(grep("  from Lapack routine",geterrmessage() ) )>0, print("internal error of  Lapack routine 'dgesdd'"),print("something else")  )
	options(show.error.messages=TRUE,  error=catch.error)
	
	for(lam1 in lambda1.set) {
		cat(paste(lam1, "" ))
		#ff <- lpsvm(A=x, d=y, k=0, nu=0,output=1, delta=10^-3, epsi=epsi, my.seed=seed)
		# 5-fold cv
		try(ff <- lpsvm(A=x, d=y, k=5, nu=0,output=0, delta=10^-3, epsi=lam1, my.seed=seed))
		ff.list[[as.character(lam1)]]<-ff
		# important to empty ff, otherweise if the error occupes in the next step i+1 take the model from the step i
		ff<- NULL
	}
	
	# choose the lam1(=epsi) with max testing set correctness:  testCorr
	# revert to default
	options(error = NULL)
	
	print("cv done ")
	# all train Correctness values, find optimal epsi (max correctness)
	CV.testCorr<- sapply(ff.list, function(tmp.list) tmp.list$testCorr )
	# plot(names(CV.testCorr), CV.testCorr, main="CV correctness", xlab="lam1", type="b",ylim=c(0,100))
	opt.lam1<- as.numeric(names(CV.testCorr)[which.max(CV.testCorr)])
	
	# run L1 norm svm on the whole data with optimal lam1!
	f.final<- lpsvm(A=x, d=y, k=0, nu=0,output=1, delta=10^-3, epsi=opt.lam1, my.seed=seed)
	print(paste(" coef w by optimal lam1", opt.lam1))
	#print(f.final$w)
	
	return(f.final)
}
