\name{sim.data}
\alias{sim.data}

\title{ Simulation of microarray data }
\description{
 Simulation of 'n' samples. Each sample has 'sg' genes,  only 'nsg' of them are called significant and 
 have influence on class labels. All other '(ng - nsg)' genes are called ballanced.
 All gene ratios are drawn from a multivariate normal distribution. 
 There is a posibility to create blocks of highly correlated genes.     
}
\usage{
sim.data(n = 256, ng = 1000, nsg = 100,
		 p.n.ratio = 0.5, 
		 sg.pos.factor= 1, sg.neg.factor= -1,
		 # correlation info:
		 corr = FALSE, corr.factor = 0.8,
		 # block info:
		 blocks = FALSE, n.blocks = 6, nsg.block = 1, ng.block = 5, 
		 seed = 123, \dots)
}
\arguments{
  \item{n}{ number of samples, logistic regression works well if \eqn{n>200}!}
  \item{ng}{ number of genes}
  \item{nsg}{ number of significant genes}
  \item{p.n.ratio}{ ratio between positive and negative significant genes (default 0.5)} 

  \item{sg.pos.factor}{ impact factor of \_positive\_ significant genes on the classifaction, default:  1  } 
  \item{sg.neg.factor}{ impact factor of \_negative\_ significant genes on the classifaction,default:  -1  } 
  \item{corr}{ are the genes correalted to each other? (default FALSE). see Details}
  \item{corr.factor}{ correlation factorfor genes, between 0 and 1 (default 0.8)   }

 	\item{blocks}{ are blocks of highly correlated genes are allowed? (default FALSE)}
  \item{n.blocks}{ number of blocks }
  \item{nsg.block}{ number of significant genes per block}
  \item{ng.block}{ number of genes per block}
  \item{seed}{ seed}
  \item{...}{ additional argument(s) }
}


\details{
   If no blockes (n.blocks=0 or blocks=FALSE) are defined and corr=TRUE
   create covarance matrix for all genes! with decrease of correlation :  \eqn{cov(i,j)=cov(j,i)= corr.factor^(i-j)}
}

\value{
  \item{ x }{matrix of simulated data. Genes in rows and samples in columns}
  \item{ y }{named vector of class labels}
  \item{ seed }{seed}
}
\author{ Wiebke Werft, Natalia Becker }

 
\seealso{ \code{\link[MASS]{mvrnorm}} }

\examples{

my.seed<-123

# 1. simulate 20 samples, with 100 genes in each. Only the first two genes
# have an impact on the class labels.
# All genes are assumed to be i.i.d. 
train<-sim.data(n = 20, ng = 100, nsg = 3, corr=FALSE, seed=my.seed )
print(str(train)) 

# 2. change the proportion between positive and negative significant genes 
#(from 0.5 to 0.8)
train<-sim.data(n = 20, ng = 100, nsg = 10, p.n.ratio = 0.8,  seed=my.seed )
rownames(train$x)[1:15]
# [1] "pos1" "pos2" "pos3" "pos4" "pos5" "pos6" "pos7" "pos8" 
# [2] "neg1" "neg2" "bal1" "bal2" "bal3" "bal4" "bal5"

# 3. assume to have correlation for positive significant genes, 
# negative significant genes and 'balanced' genes separatly. 
train<-sim.data(n = 20, ng = 100, nsg = 10, corr=TRUE, seed=my.seed )
#cor(t(train$x[1:15,]))

# 4. add 6 blocks of 5 genes each and only one significant gene per block.
# all genes in the block are correlated with constant correlation factor
#  corr.factor=0.8 		
train<-sim.data(n = 20, ng = 100, nsg = 6, corr=TRUE, corr.factor=0.8,
			 blocks=TRUE, n.blocks=6, nsg.block=1, ng.block=5, seed=my.seed )
print(str(train)) 
# first block
#cor(t(train$x[1:5,]))
# second block
#cor(t(train$x[6:10,]))

}


\keyword{ distribution }
\keyword{ multivariate }
