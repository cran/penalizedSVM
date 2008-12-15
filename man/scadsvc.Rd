\name{scadsvc}
\alias{scadsvc}

\title{ Fit SCAD SVM model  }
\description{
  SVM with variable selection (clone selection) using SCAD penalty.
}
\usage{
scadsvc(lambda = 0.01, x, y, a = 3.7, tol = 10^(-4), class.weights = NULL)
}

\arguments{
  \item{lambda}{ tuning parameter in SCAD function (default : 0.01)  }
  \item{x}{ n-by-d data matrix to train (n chips/patients, d clones/genes) }
  \item{y}{ vector of class labels  -1 or 1's (for n chips/patiens ) }
  \item{a}{ tuning parameter in scad function (default: 3.7) }
  \item{tol}{ the cut-off value to be taken as 0 }
  \item{class.weights}{ a named vector of weights for the different classes, 
  											used for asymetric class sizes. 
  											Not all factor levels have to be supplied (default weight: 1).
  										  All components have to be named. (default: NULL)}
}


\value{
  a list of 
 		\item{w }{ coefficients of the hyperplane }
	 	\item{b }{ intercept of the hyperplane}
	  \item{xind }{ the index of the selected features (genes) in the data matrix. }
	  \item{index }{ the index of the resulting support vectors in the data matrix.} 
	  \item{type }{type of svm, from svm function  }
	  \item{lam.opt }{ optimal lambda  }
	  \item{ gacv }{ corresponding gacv }
}
\details{
  Adopted from Matlab code: http://www4.stat.ncsu.edu/~hzhang/software.html
}
\references{
Zhang, H. H., Ahn, J., Lin, X. and Park, C. (2006). \emph{ Gene selection using
support vector machines with nonconvex penalty.} Bioinformatics, \bold{22}, pp. 88-95. 

 }

\author{ Axel Benner }


\seealso{  \code{\link{findgacv.scad}},  \code{\link{predict.penSVM}}, \code{\link{sim.data}}  }
\examples{


# simulate data
train<-sim.data(n = 200, ng = 100, nsg = 10, corr=FALSE, seed=12)
print(str(train)) 
	
# train data	
model <- scadsvc(as.matrix(t(train$x)), y=train$y, lambda=0.01)
print(model)


}
