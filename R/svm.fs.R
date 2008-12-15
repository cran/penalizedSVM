svm.fs <- function (x, ...)
  UseMethod ("svm.fs")



`svm.fs.default` <-
function(x,y, fs.method="1norm", cross.outer= 0, lambda1.set,  lambda2.set=NULL, calc.class.weights=FALSE, seed=240907, maxIter=NULL,...){

##  Input:
#			x: n-by-d data matrix to train (n chips/patients, d clones/genes, d>>n )
#			y: column vector (or factor vector) of target {-1, 1}'s (for n chips/patiens )
#			lam2.range : lambda2 range for elastic net (DrHSVM)

# feature selection - L1 or 1norm svm

	#require(MCRestimate)

	possible.fs<-c("1norm", "scad")
	nn<-length(y) # number of cases (patients)
	nlevels.class<- nlevels(as.factor(y))
	levels.class <- levels(as.factor(y))
	#print("!")
	
	# for scad
	if (calc.class.weights){
		class.weights = 100/ table(y)
	}else class.weights =  NULL

	#checks
	if (! (fs.method %in% possible.fs ))  stop(paste("You have to use one of following fecture selection methods", possible.fs))

	# if cross.outer>0 use cv !!!!!
	if(cross.outer>0){
		print(paste(cross.outer, "-fold cross validation ", sep="" ))
		res.cv<- .run.cv(x=x,y=y, fs.method=fs.method, lambda1.set=lambda1.set, cross.outer=cross.outer, class.weights=class.weights, seed=seed)
		print("cross validation done...")
	}   # end of if cross.outer>0 ###

  ########################################################################################
	# create final model

	# set seed again
	if (!is.null(seed)) set.seed(seed)
	
	if (fs.method %in% c("scad") ){
	#  take all posible lamda values and apply to the whole data set
	# 1.do scv scad  + gacv --> optimal lambda
		model<-run.scad(x=x,y=y, lambda1.set=lambda1.set, class.weights=class.weights)
	}
	
	
	if (fs.method=="1norm"){
	#epsi - tuning parameter !
	# find the best epsi via k-fold cv,  get the finla model with optimal epsi
		model<-run.1norm(x=x,y=y,k=5,nu=0, lambda1.set=lambda1.set, output=1, seed=seed)
	}
	

	########################################################################################
	
	# if no cv is done, no correct.prediction results ;-)
	if (cross.outer == 0) {
		rv <- list(classes=as.factor(y),
							sample.names = names(y),
							class.method=paste("svm",fs.method ),
							cross.outer=cross.outer,
							seed = seed,
							model =model )
	} else {
	## all data are collected
		rv <- list(votes=res.cv$vote.table,
							classes=as.factor(y),
							table=res.cv$confusion,
							#sample.names = rownames(vote.matrix),
							sample.names = names(y),
							gene.names = colnames(x),
							correct.prediction=res.cv$res$correct.prediction,
							correct.class.vote=res.cv$res$correct.class.vote,
							class.method=paste("svm",fs.method ),
							cross.outer=cross.outer,
							seed = seed,
							#cross.inner=cross.inner,
							#cross.repeat=cross.repeat,
							#sample.names=sample.names,
							model =model,
							lambda1.set=lambda1.set,
							lambda2.set=lambda2.set
							,cv.info=res.cv$model.info.list 
							)
	}

	#class(rv) <- "MCRestimate"
	class(rv) <- "penSVM"
	return(rv)
	# plot.MCRestimate
	# plot.MCRestimate ( rv, rownames.from.object=TRUE )
	
}

