\name{scadsvc}
\alias{scadsvc}



\title{ Fit SCAD SVM model  }
\description{
  SVM with variable selection (clone selection) using SCAD penalty.
}
\usage{
scadsvc(lambda1 = 0.01, x, y, a = 3.7, tol= 10^(-4), class.weights= NULL,
 seed=123, maxIter=700, verbose=TRUE)
}

\arguments{
  \item{lambda1}{ tuning parameter in SCAD function (default : 0.01)  }
  \item{x}{ n-by-d data matrix to train (n chips/patients, d clones/genes) }
  \item{y}{ vector of class labels  -1 or 1\'s (for n chips/patiens ) }
  \item{a}{ tuning parameter in scad function (default: 3.7) }
  \item{tol}{ the cut-off value to be taken as 0 }
  \item{class.weights}{ a named vector of weights for the different classes, 
  											used for asymetric class sizes. 
  											Not all factor levels have to be supplied (default weight: 1).
  										  All components have to be named. (default: NULL)}
  \item{seed}{ seed }
	\item{maxIter}{ maximal iteration, default: 700 }
	\item{verbose}{ verbose, default: TRUE }

}


\value{
 		\item{w}{coefficients of the hyperplane.}
	 	\item{b}{intercept of the hyperplane.}
	  \item{xind}{the index of the selected features (genes) in the data matrix. }
	  \item{xqx}{internal calculations product \eqn{xqx = 0.5 * x1 * inv_Q * t(x1)}, see code for more details.  }
	  \item{fitted}{fit of hyperplane f(x) for all _training_ samples with reduced set of features. }
	  \item{index}{the index of the resulting support vectors in the data matrix.} 
	  \item{type}{type of svm, from svm function.  }
	  \item{lambda1}{optimal lambda1.  }
	  \item{gacv}{corresponding gacv.}
	  \item{iter}{ nuber of iterations. }
}



\details{
  Adopted from Matlab code: http://www4.stat.ncsu.edu/~hzhang/software.html
}

\references{
Zhang, H. H., Ahn, J., Lin, X. and Park, C. (2006). \emph{ Gene selection using
support vector machines with nonconvex penalty. Bioinformatics, \bold{22}, pp. 88-95}. 
}

\author{ Axel Benner }


\seealso{  \code{\link{findgacv.scad}},  \code{\link{predict.penSVM}}, \code{\link{sim.data}}  }
\examples{


# simulate data
train<-sim.data(n = 200, ng = 100, nsg = 10, corr=FALSE, seed=12)
print(str(train)) 
	
# train data	
model <- scadsvc(as.matrix(t(train$x)), y=train$y, lambda=0.01)
print(str(model))

print.scad.svm(model)

}
