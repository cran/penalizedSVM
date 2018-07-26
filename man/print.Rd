\name{print}

\alias{print.penSVM}

\title{ Print Function for FS SVM}
\description{
 Print Function for FS SVM
}
\usage{
\method{print}{penSVM}(x,\dots)
}

\arguments{
  \item{x}{ model trained by scad or 1norm svm of class PenSVM }
  \item{\dots}{ additional argument(s) }
}



\author{ Natalia Becker }

\seealso{ \link[e1071]{svm}, \code{\link{svmfs}} } 

\examples{

seed<- 123
train<-sim.data(n = 20, ng = 100, nsg = 10, corr=FALSE, seed=seed )
print(str(train)) 

# for presentation don't check  all lambdas : time consuming! 
model <- scadsvc(as.matrix(t(train$x)), y=train$y, lambda=0.05)
print(str(model))

print(model)
}
