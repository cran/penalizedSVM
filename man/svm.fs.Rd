\name{svm.fs}

\alias{svm.fs}
\alias{svm.fs.default}
\alias{run.scad}
\alias{run.1norm}


\title{  Fits  SVM mit variable selection using penalties.}
\description{
  Fits  SVM with variable selection (clone selection) using penalties SCAD and L1 norm.
}
\usage{
\method{svm.fs}{default}(x, y, fs.method = "1norm", cross.outer = 0, lambda1.set, lambda2.set = NULL, calc.class.weights = FALSE, seed = 240907, maxIter = NULL,...)
run.scad(x,y,  lambda1.set=NULL, class.weights)
run.1norm(x,y,k=5,nu=0, lambda1.set=NULL, output=1, seed=seed)


}
\arguments{
  \item{x}{ input matrix with genes in columns and samples in rows!  }
  \item{y}{ vector of class labels }
  \item{fs.method}{ feature selection method. Availible 'scad' and '1norm' }
  \item{cross.outer}{ fold of outer cross validation, default is 0, no cv.     }
  \item{lambda1.set}{ set of tuning parameters lambda1 }
  \item{lambda2.set}{ set of tuning parameters lambda2, not yet in use }
  \item{calc.class.weights}{ calculate class.weights for SVM, default: FALSE }
  \item{class.weights}{ a named vector of weights for the different
    classes, used for asymetric class sizes. Not all factor levels have
    to be supplied (default weight: 1). All components have to be named.  }
  \item{k}{ k-fold cross validation, default: 5 }
 \item{nu}{ nu:  weighted parameter
  								\itemize{
                  \item  1 - easy estimation, 
                  \item  0 - hard estimation, 
                  \item any other value  - used as nu by the algorithm, default - 0. 
                  }
            }   
  \item{output}{ 0 - no output, 1 - produce output, default is 0 }
  \item{seed}{ seed }
  \item{maxIter}{ maximal iteration, default: not used yet  }
  \item{...}{ additional argument(s) }
}

\details{
  The goodness of the model is highly correlated with the choice of tuning parameter lambda.
  Therefore the model is trained with different lambdas and the best model with optimal tuning parameter 
  is used in futher analysises.
  
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
  \item{cross.outer }{ outer cv  }
  \item{seed }{ seed}
  \item{model }{ final model 
  	\itemize{
	  	\item w - coefficients of the hyperplane 
		  \item b - intercept of the hyperplane
		  \item xind -   the index of the selected features (genes) in the data matrix.
		  \item index - the index of the resulting support vectors in the data matrix. 
		  \item type - type of svm, from svm function 
		  \item lam.opt  - optimal lambda  
		  \item  gacv - corresponding gacv 
	  }
	}
} 

\references{ 

Zhang, H. H., Ahn, J., Lin, X. and Park, C. (2006). \emph{ Gene selection using
support vector machines with nonconvex penalty.} Bioinformatics, \bold{22}, pp. 88-95. 

Fung, G. and Mangasarian, O. L. (2004).  \emph{ A feature selection newton method for 
support vector machine classification.} Computational Optimization and Applications Journal, \bold{28(2)}, pp. 185-202.   

}

\author{ Natalia Becker natalia.becker at dkfz.de }


\seealso{ \code{\link{predict.penSVM}},   \code{\link[e1071]{svm}} (in package \pkg{e1071}) }

\examples{

my.seed<- 123

train<-sim.data(n = 200, ng = 100, nsg = 10, corr=FALSE, seed=my.seed )
print(str(train)) 


# train SCAD SVM ####################
# define set values of tuning parameter lambda1 for SCAD 
lambda1.scad <- c (seq(0.01 ,0.05, .01),  seq(0.1,0.5, 0.2), 1 ) 
# for presentation don't check  all lambdas : time consuming! 
lambda1.scad<-lambda1.scad[2:3]
# 
# train SCAD SVM
fit.scad<- svm.fs(x=t(train$x),y=train$y, fs.method="scad", cross.outer= 0, lambda1.set=lambda1.scad, seed=my.seed)
	
	
# train 1NORM SVM 	################	
# define set values of tuning parameter lambda1 for 1norm
epsi.set<-vector(); for (num in (1:9)) epsi.set<-sort(c(epsi.set, c(num*10^seq(-5, -1, 1 ))) )
# for presentation don't check  all lambdas : time consuming! 
lambda1.1norm <- 	epsi.set[c(3,5)] # 2 params

# train 1norm SVM
# time consuming: for presentation only for the first 100 samples    
## DON'T RUN : fit.1norm<- svm.fs(x=t(train$x),y=train$y, fs.method="1norm", cross.outer= 0, lambda1.set=lambda1.1norm, seed=my.seed)
fit.1norm<- svm.fs(x=t(train$x)[1:100,],y=train$y[1:100], fs.method="1norm", cross.outer= 0, lambda1.set=lambda1.1norm, seed=my.seed)
	
    

}

\keyword{ models }
\keyword{ multivariate }
\keyword{ optimize }

