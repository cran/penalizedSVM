run.scad<- function(x,y,  lambda1.set=NULL, class.weights){
	# for scad 
	#  take all posible lamda values and apply to the whole data set
	# 1.do scv scad  + gacv --> optimal lambda1
	
	ff.list<-list()
	
	# try with gacv
	# fine grid at the begin of the interval !!!
	##  skip lambda1=0 !!!!!
	if (is.null(lambda1.set)) lambda1.set <- c (0.001, 0.01,  seq(0.1 ,1, 0.1),5)
	
	for(lam1 in lambda1.set) {
		cat(lam1)
		try(ff <- scadsvc(as.matrix(x), y, lambda=lam1, class.weights=class.weights))
		ff.list[[as.character(lam1)]]<-ff
		# important to empty ff, otherweise if the error occupes in the next step i+1 take the model from the step i
		ff<- NULL
	}
	
	#if we have at least one model:
	if (any (sapply(ff.list, length)>1) ){
		#  gacv
		# important : y should be -1 and 1 !!!!!
		y.gacv<- factor(y)
		levels(y.gacv)=c(-1,1)
		y.gacv<-as.numeric(as.character(y.gacv))
		
		gacv<-rep(NA, length(ff.list))
		names(gacv) <- names(ff.list)
		# sort out the empty models
		for (i in 1:length(ff.list) )  if (length(ff.list[[i]])>1) gacv[i]<- findgacv.scad(y.gacv, model=ff.list[[i]])
		
		gacv<-na.omit(gacv)
		
		min.gacv<-min(gacv, na.rm=TRUE)
		
		# if we have several lamdas mit minimal gacv value, chose those with min number of features (genes)
		set.min.gacv<- names(gacv)[ gacv == min.gacv ]
		if (length(set.min.gacv)>1){
			ll<-set.min.gacv
			# number of features in end model
			ll.w.length<-unlist(lapply(ff.list[as.character(ll)], function(d){length(d$w)} ) )
			lam.opt.final <- ll[which.min(ll.w.length)]
		}else lam.opt.final <- lambda1.set[which(gacv==min.gacv)]
		
		# define the right model
		f.final <- ff.list[[as.character(lam.opt.final)]]
		f.final$lam.opt <- lam.opt.final
		f.final$gacv<- min.gacv
		f.final$class.weights<- class.weights
		# skip xqx slot: don't need any more
		f.final<-f.final[- which(names(f.final) == "xqx")]
	} else{
		f.final<- NULL
	}
	
	print("gacv is done")
	return(f.final)
}
