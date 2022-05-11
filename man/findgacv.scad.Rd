\name{findgacv.scad}
\alias{findgacv.scad}

\title{ Calculate Generalized Approximate Cross Validation Error Estimation for SCAD SVM model}

\description{
 calculate generalized approximate cross validation error (GACV) estimation 
 for SCAD SVM model}

\usage{
findgacv.scad(y, model)
}
\arguments{
  \item{y}{ vector of class labels (only for 2 classes) }
  \item{model}{ list, describing  SCAD SVM model, produced by function scadsvc }
}

\value{
  returns the GACV value
}



\references{

Zhang, H. H., Ahn, J., Lin, X. and Park, C. (2006). \emph{ Gene selection using
support vector machines with nonconvex penalty.} Bioinformatics, \bold{22}, pp. 88-95. 

Wahba G., Lin, Y. and Zhang, H. (2000). \emph{GACV for support vector machines, or, another way
to look at margin-like quantities, in A. J. Smola, P. Bartlett, B. Schoelkopf and D. Schurmans
(eds)}, Advances in Large Margin Classifiers, MIT Press, pp. 297-309.
}

\author{Natalia Becker\cr
\email{natalie_becker@gmx.de}
}

\seealso{  \code{\link{scadsvc}},  \code{\link{predict.penSVM}}, \code{\link{sim.data}}     }
\examples{

# simulate data
train<-sim.data(n = 200, ng = 100, nsg = 10, corr=FALSE, seed=12)
print(str(train)) 
	
# train data	
ff <- scadsvc(as.matrix(t(train$x)), y=train$y, lambda=0.01)
print(str(ff))

# estimate gacv error
(gacv<- findgacv.scad(train$y, model=ff))

}

\keyword{ models}
