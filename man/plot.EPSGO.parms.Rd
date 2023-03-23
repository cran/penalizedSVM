\name{.plot.EPSGO.parms}

\alias{.plot.EPSGO.parms}

\title{Plot Interval Search Plot Visited Points and the Q Values. }
\description{
For interval search plot visited points and the Q values (=Ytrain) exclude: for D=1 make  an additional plot: skip values for empty model, for example: Ytrain.exclude=10^16.
}

\usage{
.plot.EPSGO.parms(Xtrain, Ytrain,bounds, Ytrain.exclude=10^16, plot.name=NULL )
}

\arguments{
  \item{Xtrain}{X points to train}
  \item{Ytrain}{Y points to train}
  \item{bounds}{ bounds for parameters, see examples}
  \item{Ytrain.exclude}{If exclude for Ytrain exists, skip those  points. Defaults to 10^16.}
  \item{plot.name}{Defaults to \code{NULL}}
}

\value{None, only graphs are created.
} 


\details{
	The code source was adopted from MATLAB originals, special thanks to Holger Froehlich.
}

\author{Natalia Becker\cr
\email{natalie_becker@gmx.de}
}

\references{
Froehlich, H. and Zell, A. (2005) "Effcient parameter selection for support vector
machines in classification and regression via model-based global optimization"
\emph{In Proc. Int. Joint Conf. Neural Networks,  1431-1438 }.
}

\seealso{ \code{\link{svmfs}}   }

\examples{

	\donttest{
	seed <- 123
			
	train<-sim.data(n = 200, ng = 100, nsg = 10, corr=FALSE, seed=seed )
	print(str(train)) 
			
	Q.func<- ".calc.scad"
	
	bounds=t(data.frame(log2lambda1=c(-10, 10)))
							colnames(bounds)<-c("lower", "upper")	
			
	print("start interval search")
	# computation intensive; 
	# for demostration reasons only for the first 100 features 
	# and only for 10 iterations maxIter=10, default maxIter=700
	system.time(fit<-EPSGO(Q.func, bounds=bounds, parms.coding="log2", fminlower=0, 
		 show='none', N=21,  maxevals=500, 
		 pdf.name=NULL,  seed=seed,  
		 verbose=FALSE,
		 # Q.func specific parameters:
		 x.svm=t(train$x)[,1:100], y.svm=train$y,
		 inner.val.method="cv",
		 cross.inner=5, maxIter=10 ))
									 
	print(paste("minimal 5-fold cv error:", fit$fmin, "by log2(lambda1)=", fit$xmin))
		
	print(" all lambdas with the same minimum? ")
	print(fit$ points.fmin) 
			
	print(paste(fit$neval, "visited points"))
			
			
	print(" overview: over all visitied points in tuning parameter space 
				with corresponding cv errors")
	print(data.frame(Xtrain=fit$Xtrain, cv.error=fit$Ytrain))

	# create  3 plots om one screen: 
	# 1st plot: distribution of initial points in tuning parameter space
	# 2nd plot: visited lambda points vs. cv errors
	# 3rd plot: the same as the 2nd plot, Ytrain.exclude points are excluded. 
	#    The value cv.error = 10^16 stays for the cv error for an empty model ! 
	.plot.EPSGO.parms (fit$Xtrain, fit$Ytrain,bound=bounds, 
				Ytrain.exclude=10^16, plot.name=NULL )
	} # end of \donttest
}

\keyword{ models }
\keyword{ multivariate }
\keyword{ graphs }
\keyword{ iteration }
\keyword{ optimize }





