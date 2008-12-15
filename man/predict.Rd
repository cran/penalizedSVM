\name{predict}
\alias{predict.penSVM}


\title{ Predict Method for Feature Selection SVM      }
\description{
  This function predicts values based upon a model trained by svm. 
	If class assigment is provided, confusion table, missclassification table,
	 sensitivity and specificity are calculated.
}

\usage{
\method{predict}{penSVM}(object, newdata, newdata.labels = NULL, \dots)
}

\arguments{
  \item{object}{ Object of class "penSVM", created by 'svm.fs'}
  \item{newdata}{ A matrix containing the new input data, samples in rows, features in columns }
  \item{newdata.labels}{ optional, new data class labels  }
	\item{...}{ additional argument(s) }
}

\value{
	returns a list of prediction values for classes
	\item{pred.class }{predicted class}
	\item{tab }{confusion table}
	\item{error }{missclassification error}
	\item{sensitivity }{sensitivity}
	\item{specificity }{specificity}
}

\author{ Natalia Becker }


\seealso{ \link[e1071]{svm}, \code{\link{svm.fs}} } 

\examples{

my.seed<- 123
train<-sim.data(n = 200, ng = 100, nsg = 10, corr=FALSE, seed=my.seed )
print(str(train)) 

#train standard svm
my.svm<-svm(x=t(train$x), y=train$y, kernel="linear")

# test with other data
test<- sim.data(n = 200, ng = 100, nsg = 10, seed=(my.seed+1) )

# Check accuracy standard SVM
my.pred <-ifelse( predict(my.svm, t(test$x)) >0,1,-1)
# Check accuracy:
table(my.pred, test$y)

# define set values of tuning parameter lambda1 for SCAD 
lambda1.scad <- c (seq(0.01 ,0.05, .01),  seq(0.1,0.5, 0.2), 1 ) 
# for presentation don't check  all lambdas : time consuming! 
fit.scad<- svm.fs(x=t(train$x),y=train$y, fs.method="scad", cross.outer= 0, lambda1.set=lambda1.scad[2:3], seed=my.seed)


# SCAD 
test.error.scad<-predict(fit.scad, newdata=t(test$x),newdata.labels=test$y )
 # Check accuracy SCAD SVM  	
print(test.error.scad$tab)	

#########################################
# analog for 1-norm SVM
epsi.set<-vector(); for (num in (1:9)) epsi.set<-sort(c(epsi.set, c(num*10^seq(-5, -1, 1 ))) )
lambda1.1norm <- 	epsi.set[c(3,5)] # 2 params

# train 1norm SVM
# time consuming: for presentation only for the first 100 samples    
## DON'T RUN : fit.1norm<- svm.fs(x=t(train$x),y=train$y, fs.method="1norm", cross.outer= 0, lambda1.set=lambda1.1norm, seed=my.seed)
fit.1norm<- svm.fs(x=t(train$x)[1:100,],y=train$y[1:100], fs.method="1norm", cross.outer= 0, lambda1.set=lambda1.1norm, seed=my.seed)

# L1-norm SVM
test.error.1norm<-predict(fit.1norm, newdata=t(test$x),newdata.labels=test$y )
 # Check accuracy L1-norm SVM  	
print(test.error.1norm$tab)	
}
