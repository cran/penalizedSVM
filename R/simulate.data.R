
sim.data <- function(n = 256, ng = 1000, nsg = 100, p.n.ratio = 0.5, 
												sg.pos.factor= 1, sg.neg.factor= -1,
												# correlation info:
												corr=FALSE, corr.factor=0.8,
												# block info:
												blocks=FALSE, n.blocks=6, nsg.block=1, ng.block=5, 
												seed=12345, ...){ 
####################################################################################
# simulate microarray data with 
# n - number of samples, logistic regression works well if n>200!
# ng - number of genes
# nsg - number of significant genes 
# p.n.ratio - ratio between positive and negative significant genes (default 0.5)  
# sg.pos.factor - impact factor of _positive_ significant genes on the classifaction (default 1)  
# sg.neg.factor - impact factor of _negative_ significant genes on the classifaction (default -1)  
# all other non-significant genes have in both classes balanced ratios

# for correlated blocks of genes
# n.blocks -  number of correlated blocks of genes
# nsg.block - number of significant genes per block
# ng.block - number of genes per block

# if no blockes (n.blocks=0) are defined and corr=TRUE
# create covarance matrix for all genes! with decrease of correlation 

####################################################################################

	require(MASS)
	set.seed(seed)
	# assumption I: intercept = 0, beta0=0
	b0<- 0
	
	# effect for positive genes = +log2(2) = +1
	#            negative genes =  log2(1/2)= -1 
	m.pos<-m.neg<-m.bal<- 0
	sd.pos<-sd.neg<-sd.bal<- 1
	
	if (!( (p.n.ratio>=0) & (p.n.ratio<=1))) stop("rato between positive and negative significant genes should be in [0;1]")
	pos.nsg<- floor (nsg * p.n.ratio)
	neg.nsg<- nsg - pos.nsg
	
	# covariance matrix
	sigma<- .create.covariance.matrix(sd.pos, sd.neg, sd.bal,
																		ng,nsg, pos.nsg, neg.nsg,
																		corr,corr.factor, 
																		blocks, n.blocks, nsg.block, ng.block)
	
	# better to see
	# sigma.see<-sigma; sigma.see[ sigma.see==0]<- ""; rm(sigma.see) 
			
	#  1 if pos, -1 by neg
	bX<-rep(0,ng); bX[grep("pos",rownames(sigma))] <- sg.pos.factor 
	bX[grep("neg",rownames(sigma))] <- sg.neg.factor
	bX = matrix(bX, ncol=1)
	
	# means 
	means<- rep(0,ng); means[grep("pos",rownames(sigma))] <- m.pos; means[grep("neg",rownames(sigma))] <- m.neg
	
	X <-t(mvrnorm(n, means, sigma))
	colnames(X)<-c(1:n) 
	
	# Outcome 
	L <-  t(X) %*% bX
	# error of the model 
	Y <- ifelse(runif(n) < plogis(L), 1, -1)
		
	return(list("x"=X,"y"=Y[,1], "seed"=seed))
}
