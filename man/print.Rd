\name{print}
\alias{print.scadsvm}
\alias{print.1norm.svm}


\title{ Print Function for FS SVM}
\description{
 Print Function for FS SVM
}
\usage{

\method{print}{scadsvm}(x,\dots)
\method{print}{1norm.svm}(x,\dots)

}

\arguments{
  \item{x}{ model trained by scad or 1norm svm }
  \item{...}{ additional argument(s) }
}



\author{ Natalia Becker }

\seealso{ \link[e1071]{svm}, \code{\link{svm.fs}} } 

\examples{

my.seed<- 123
train<-sim.data(n = 20, ng = 100, nsg = 10, corr=FALSE, seed=my.seed )
print(str(train)) 

# for presentation don't check  all lambdas : time consuming! 
model <- scadsvc(as.matrix(t(train$x)), y=train$y, lambda=0.05)
print(model)

}
