\name{lpsvm}

\alias{lpsvm}



\title{Fit L1-norm SVM }
\description{
  SVM mit variable selection (clone selection) using L1-norm penalty. 
  ( a fast Newton algorithm NLPSVM from Fung and Mangasarian ) 
}
\usage{
lpsvm(A, d, k = 5, nu = 0, output = 1, delta = 10^-3, epsi = 10^-4, 
seed = 123, maxIter=700)
}

\arguments{
  \item{A}{ n-by-d data matrix to train (n chips/patients, d clones/genes).  }
  \item{d}{ vector of class labels  -1 or 1's (for n chips/patiens ). }
  \item{k}{  k-fold for cv, default k=5.  }
  \item{nu}{ weighted parameter, 1 - easy estimation,
                    0  - hard estimation, any other value - used as nu by the algorithm.
                    Default : 0. }
  \item{output}{ 0 - no output, 1 - produce output, default is 0. }
  \item{delta}{ some small value, default: \eqn{10^-3}. }
  \item{epsi}{ tuning parameter.  }
  \item{seed}{ seed. }
  \item{maxIter}{ maximal iterations, default: 700. }
}

\details{
    k: k-fold for cv, is a way to divide the data set into test and training set.\cr
       if k = 0: simply run the algorithm without any correctness
         calculation, this is the default. \cr
       if k = 1: run the algorithm and calculate correctness on
         the whole data set. \cr
       if k = any value less than the number of rows in the data set:
         divide up the data set into test and training
         using k-fold method. \cr
       if k = number of rows in the data set: use the 'leave one out' (loo) method
}
\value{
  a list of 
 		\item{w }{ coefficients of the hyperplane }
	 	\item{b }{ intercept of the hyperplane}
	  \item{xind }{ the index of the selected features (genes) in the data matrix. }
	  \item{epsi }{ optimal tuning parameter epsilon  }
	  \item{ iter }{ number of iterations }
	  \item{ k }{  k-fold for cv}
	  \item{ trainCorr }{ for cv: average train correctness }
	  \item{ testCorr }{ for cv: average test correctness }
	  \item{ nu }{ weighted parameter }
	  
	 
}

\references{ 
Fung, G. and Mangasarian, O. L. (2004). A feature selection newton method for 
support vector machine classification. \emph{Computational Optimization and Applications Journal 28(2)  pp. 185-202}. 
}


\author{ Natalia Becker }
\note{ Adapted from MATLAB code  http://www.cs.wisc.edu/dmi/svm/lpsvm/
}
\seealso{  \code{\link{sim.data}}  }

\examples{
set.seed(123)
train<-sim.data(n = 20, ng = 100, nsg = 10, corr=FALSE, seed=12)
print(str(train)) 
	
# train data	
model <- lpsvm(A=t(train$x), d=train$y, k=5, nu=0,output=0, delta=10^-3, epsi=0.001, seed=12)
print(model)


}
