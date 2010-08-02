\name{print}
\alias{print.scad.svm}
\alias{print.1norm.svm}


\title{ Print Function for FS SVM}
\description{
 Print Function for FS SVM
}
\usage{
print.scad.svm(x,...)
print.1norm.svm(x,...)

}

\arguments{
  \item{x}{ model trained by scad or 1norm svm }
  \item{\dots}{ additional argument(s) }
}



\author{ Natalia Becker }

\seealso{ \link[e1071]{svm}, \code{\link{svm.fs}} } 

\examples{

seed<- 123
train<-sim.data(n = 20, ng = 100, nsg = 10, corr=FALSE, seed=seed )
print(str(train)) 

# for presentation don't check  all lambdas : time consuming! 
model <- scadsvc(as.matrix(t(train$x)), y=train$y, lambda=0.05)
print(str(model))

print.scad.svm(model)
}
