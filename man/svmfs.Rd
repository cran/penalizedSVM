\name{svmfs}

\alias{svmfs}
\alias{svmfs.default}
\alias{DrHSVM}
\alias{scad_L2.svc}


\title{  Fits  SVM with variable selection using penalties.}
\description{
  Fits  SVM with variable selection (clone selection) using penalties SCAD,  L1 norm, Elastic Net (L1 + L2 norms) and ELastic SCAD (SCAD + L1 norm). 
  Additionally tuning parameter search is presented by two approcaches: fixed grid or interval search.
  NOTE: The name of the function has been changed: svmfs instead of svm.fs! 

}
\usage{
\method{svmfs}{default}(x,y,
				fs.method = c("scad", "1norm", "scad+L2", "DrHSVM"),
				grid.search=c("interval","discrete"),
				lambda1.set=NULL,  
				lambda2.set=NULL,
				bounds=NULL, 
				parms.coding= c("log2","none"),
				maxevals=500, 
				inner.val.method = c("cv", "gacv"),
				cross.inner= 5,
				show= c("none", "final"),
				calc.class.weights=FALSE,
				class.weights=NULL, 
				seed=123, 
				maxIter=700, 
				verbose=TRUE,
				\dots) 
}
\arguments{
  \item{x}{ input matrix with genes in columns and samples in rows!  }
  \item{y}{ numerical vector of class labels, -1 , 1 }
  \item{fs.method}{ feature selection method. Availible 'scad',  '1norm' for 1-norm,  "DrHSVM" for Elastic Net 
  and "scad+L2" for Elastic SCAD }
  \item{grid.search}{ chose the search method for tuning lambda1,2: 'interval' or 'discrete', default: 'interval'}
  \item{lambda1.set}{ for fixed grid search: fixed grid for lambda1, default: NULL  }
  \item{lambda2.set}{ for fixed grid search: fixed grid for lambda2, default: NULL  }
  \item{bounds}{ for interval grid search: fixed grid for lambda2, default: NULL  }
  \item{parms.coding}{ for interval grid search: parms.coding: none or log2  , default: log2  }
  \item{maxevals}{ the maximum number of DIRECT function evaluations, default: 500.   }
  %\item{cross.outer}{ fold of outer cross validation, default is 0, no cv.     }
  \item{calc.class.weights}{ calculate class.weights for SVM, default: FALSE }
  \item{class.weights}{ a named vector of weights for the different
    classes, used for asymetric class sizes. Not all factor levels have
    to be supplied (default weight: 1). All components have to be named.  }
  \item{inner.val.method}{ method for the inner validation: cross validation, gacv , default cv  }
  \item{cross.inner}{ 'cross.inner'-fold cv, default: 5 }
  \item{show}{ for interval search: show plots of  DIRECT algorithm:    none, final iteration, all iterations. Default: none  }
  \item{seed}{ seed }
  \item{maxIter}{ maximal iteration, default: 700 }
  \item{verbose}{ verbose?, default: TRUE }
  \item{\dots}{ additional argument(s) }
}

\details{
  
  The goodness of the model is highly correlated with the choice of tuning parameter lambda.
  Therefore the model is trained with different lambdas and the best model with optimal tuning parameter 
  is used in futher analysises.
  For very small lamdas is recomended to use maxIter, otherweise the algorithms is slow or might not converge.
  
  The Feature Selection methods are using different techniques for finding optimal tunung parameters
  By SCAD SVM Generalized approximate cross validation (gacv) error is calculated for each pre-defined tuning parameter. 
  
  By L1-norm SVM the cross validation (default 5-fold) missclassification error is calculated for each lambda.  
   After training and cross validation, the optimal lambda with minimal missclassification error is choosen,
  and a final model with optimal lambda is created for the whole data set.    
}




\value{
  \item{classes }{ vector of class labels as input 'y' }
  \item{sample.names }{sample names}
  \item{class.method }{feature selection method}
  %\item{cross.outer }{ outer cv  }
  \item{seed }{ seed}
  \item{model }{ final model 
  	\itemize{
	  	\item w - coefficients of the hyperplane 
		  \item b - intercept of the hyperplane
		  \item xind -   the index of the selected features (genes) in the data matrix.
		  \item index - the index of the resulting support vectors in the data matrix. 
		  \item type - type of svm, from svm function 
		  \item lam.opt  - optimal lambda  
		  \item gacv - corresponding gacv 
	  }
	}
} 

\references{ 
Becker, N., Werft, W., Toedt, G., Lichter, P. and Benner, A.(2009) PenalizedSVM: a R-package for feature selection SVM classification, 
Bioinformatics, 25(13),p 1711-1712

}

\author{Natalia Becker\cr
\email{natalie_becker@gmx.de}
}

\seealso{ \code{\link{predict.penSVM}},   \code{\link[e1071]{svm}} (in package \pkg{e1071}) }

\examples{


\donttest{
		
		seed<- 123
		
		train<-sim.data(n = 200, ng = 100, nsg = 10, corr=FALSE, seed=seed )
		print(str(train)) 
		
		
		### Fixed grid ####
		
		# train SCAD SVM ####################
		# define set values of tuning parameter lambda1 for SCAD 
		lambda1.scad <- c (seq(0.01 ,0.05, .01),  seq(0.1,0.5, 0.2), 1 ) 
		# for presentation don't check  all lambdas : time consuming! 
		lambda1.scad<-lambda1.scad[2:3]
		# 
		# train SCAD SVM
		
		# computation intensive; for demostration reasons only for the first 100 features 
		# and only for 10 Iterations maxIter=10, default maxIter=700
		system.time(scad.fix<- svmfs(t(train$x)[,1:100], y=train$y, fs.method="scad", 
  		cross.outer= 0, grid.search = "discrete",  
  		lambda1.set=lambda1.scad,
  		parms.coding = "none", show="none",
  		maxIter = 10, inner.val.method = "cv", cross.inner= 5,
  		seed=seed, verbose=FALSE) 	)
			
		print(scad.fix)
			
		# train 1NORM SVM 	################	
		# define set values of tuning parameter lambda1 for 1norm
		#epsi.set<-vector(); for (num in (1:9)) epsi.set<-sort(c(epsi.set,
		#    c(num*10^seq(-5, -1, 1 ))) )
		## for presentation don't check  all lambdas : time consuming! 
		#lambda1.1norm <- 	epsi.set[c(3,5)] # 2 params
		#
		### train 1norm SVM
		## time consuming: for presentation only for the first 100 features    
		#norm1.fix<- svmfs(t(train$x)[,1:100], y=train$y, fs.method="1norm", 
		#			cross.outer= 0, grid.search = "discrete",  
		#			lambda1.set=lambda1.1norm,
		#			parms.coding = "none", show="none",
		#			maxIter = 700, inner.val.method = "cv", cross.inner= 5,
		#			seed=seed, verbose=FALSE ) 	
		#	
		#	print(norm1.fix)   
		
		### Interval  search  ####
		
		
		seed <- 123
		
		train<-sim.data(n = 200, ng = 100, nsg = 10, corr=FALSE, seed=seed )
		print(str(train)) 
		
		
		test<-sim.data(n = 200, ng = 100, nsg = 10, corr=FALSE, seed=seed+1 )
		print(str(test)) 
		
				
		bounds=t(data.frame(log2lambda1=c(-10, 10)))
						colnames(bounds)<-c("lower", "upper")	
		
		# computation intensive; for demostration reasons only for the first 100 features 
		# and only for 10 Iterations maxIter=10, default maxIter=700
		print("start interval search")
			system.time( scad<- svmfs(t(train$x)[,1:100], y=train$y,
			 fs.method="scad", bounds=bounds, 
			 cross.outer= 0, grid.search = "interval",  maxIter = 10, 
			 inner.val.method = "cv", cross.inner= 5, maxevals=500,
			 seed=seed, parms.coding = "log2", show="none", verbose=FALSE ) )
		print("scad final model")
		print(str(scad$model))
				
		(scad.5cv.test<-predict.penSVM(scad, t(test$x)[,1:100], newdata.labels=test$y)   )
		
		
		print(paste("minimal 5-fold cv error:", scad$model$fit.info$fmin, 
		"by log2(lambda1)=", scad$model$fit.info$xmin))
		
		print(" all lambdas with the same minimum? ")
		print(scad$model$fit.info$ points.fmin) 
		
		print(paste(scad$model$fit.info$neval, "visited points"))
		
		
		print(" overview: over all visitied points in tuning parameter space 
		with corresponding cv errors")
		print(data.frame(Xtrain=scad$model$fit.info$Xtrain, 
					cv.error=scad$model$fit.info$Ytrain))
		# 						 
		
		# create  3 plots on one screen: 
		# 1st plot: distribution of initial points in tuning parameter space
		# 2nd plot: visited lambda points vs. cv errors
		# 3rd plot: the same as the 2nd plot, Ytrain.exclude points are excluded. 
		# The value cv.error = 10^16 stays for the cv error for an empty model ! 
		.plot.EPSGO.parms (scad$model$fit.info$Xtrain, scad$model$fit.info$Ytrain,
				bound=bounds, Ytrain.exclude=10^16, plot.name=NULL )
		
} # end of \donttest

}

\keyword{ models }
\keyword{ multivariate }
\keyword{ optimize }

