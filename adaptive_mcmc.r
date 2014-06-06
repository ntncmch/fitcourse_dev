#
#
# program code for section 3.2, multivariate example
# here we update with a random walk using a multivariate normal distribution,
# with variance set to sigma^2*A, and tune the sigma^2 parameter, and also adapt
# A using the second half of the chain, the example comes from Roberts and Rosenthal (examples
# of adaptive MCMC 2006)
#
# 	written by
# 	Y Fan
#	7th April 2012
#
#################################################################
#################################################################
# The following code could be adjusted for any multivariate example. The only
# necessary change is to change the function that gives the target distribution.
# You may well want to change the number of iterations in the Markov chain (niter)
# and the starting value for theta (the log of standard deviation of the proposal distribution).
# Less frequently you might want to change the target acceptance proability to something other
# than 0.234.
#
# The last lines in the program are optional and draw diagnostic plots that monitor
# the Robins Monro search.
#
# For more complicated problems that include multivariate RWMH, you might want to copy and paste
# parts of this code into your program so as to estimate the variance of a proposal distribution,
# rather than having to specify that variance. (Usually, little of your code need be deleted.)
################################################################################################


library(MASS)
library(mvtnorm)

#function which adjusts the value of theta at the ith iteration
update.sigma <- function(sigma2, acc, p=p, i) {
	c=((1-1/d)*sqrt(2*pi)*exp(alpha^2/2)/(2*alpha) + 1/(d*p*(1-p)))
	Theta=log(sqrt(sigma2))
	Theta=Theta+c*(acc-p)/max(200, i/d)
	return(exp(Theta))
}


update.cov<-function(sigMat, i, thetaM, theta){
#function to recursively update covariance matrix, as part of adaptive MCMC sampling, updating the covariance matrix
	epsilon=1/i
	thetaM2=((thetaM*i)+theta)/(i+1)
	sigMat=(i-1)/i*sigMat + thetaM%*%t(thetaM)-(i+1)/i*thetaM2%*%t(thetaM2)+1/i*theta%*%t(theta) + epsilon*diag(d)
	return(list(sigMat=sigMat, thetaM=thetaM2))
}




#function that evaluations the target distribution
calc.target <- function(par) {
# Set the distribution of your target distribution here.
# That is, you must change this function to give your target distribution.
# (currenly this is Multivariate normal with mean 0, and covariance cov.true) 

#This part can be replaced when calculation is needed for a different target
	if (length(par)!=d) {
		cat("stop! wrong dimension of parameter!")
	}
	else {
		return(dmvnorm(par, mean=rep(0, d), sigma=cov.true))
	}
}




#initialise the MCMC program
niter=5000      # niter is the number of iterations in the Markov chain.
                # Change it to the value you want.
d=50 #dimension of parameters to be updated
output<-rep(0, d) #output parameter values
output.mat<-output #records all MCMC output for the parameters

sigMat<-diag(d)  #initial proposal covarince matrix, can be changed if a better one is available
acc.vec<-rep(NA, niter) #records the acceptance probabilities from each MCMC iteration

#############################################
# Begin
#initialise the RM section
#############################################
#the optimal acceptance probability for the multivariate case
pstar=0.234
alpha=-qnorm(pstar/2)
n0=round(5/(pstar*(1-pstar)))
#iMax, is the max number of iterations before the last restart
iMax=100
Numbig=0
Numsmall=0
sigma=1        #an arbitrary starting value for sigma, equivalent to theta=ln(sigma)=0
sigma2=sigma^2 #variance
sigma.start<- sigma
sigma.vec<-sigma
num.restart=0
i=1
#############################################
# End 
#initialise the RM section
#############################################

for (j in c(2:niter)) {

#propose a new value of theta
	output.prop<-mvrnorm(1, mu=output, Sigma=sigma2*sigMat)
	pi.old<- calc.target(output)
	pi.new<-calc.target(output.prop)
	u<-runif(1)
	acc=min(1, pi.new/pi.old)
	acc.vec=c(acc.vec,acc)
	i=i+1

	if (u < acc) {
		output<-output.prop
	}
	output.mat<-c(output.mat, output)
	

#################################################
# Begin
# update covariance matrix with adaptive MCMC
#################################################
	if (i > 100) {
		if (i==101) {
			sigMat=cov(output.mat)
			thetaM=apply(output.mat, 2, mean)
		} else
		{
			tmp=update.cov(sigMat, i, thetaM, output, sigma)
			sigMat=tmp$sigMat
			thetaM=tmp$thetaM
		}
	}

#################################################
# End 
# update covariance matrix with adaptive MCMC
#################################################



###############################################
# Begin 
# update sigma using RM
###############################################
	if (j>n0) {
		sigma<-update.sigma(sigma2, is.ACC, pstar, j)
		sigma2=sigma^2
		i=i+1
		sigma.vec<-c(sigma.vec, sigma)
		if ((j <= (iMax+n0)) && (Numbig<5 || Numsmall<5)) {
			Toobig<- (sigma > (3*sigma.start))
			Toosmall<-(sigma < (sigma.start/3))

			if (Toobig || Toosmall) {
						#restart the algorithm
				cat("restart the program at", i, "th iteration", "\n")
				sigma.restart<-c(sigma.restart, sigma)
				Numbig<- Numbig + Toobig
				Numsmall <- Numsmall + Toosmall
				j<-n0
				sigma.start<-sigma
			}
		} #end iMax
	}
###############################################
# Begin
# update sigma using RM
###############################################


} #end niter
#end of MCMC

##############################################################################################
## Diagnostic Plotting - The remaining lines are optional and the program will run without them
###############################################################################################
#
#
#

accept<-accept[-1]
#calculate running mean acceptance rate
meanacc<-rep(NA, length(acc.vec))
for (i in c(1:length(acc.vec))) {
	meanacc[i]=     mean(acc.vec[round(i/2) :i])
}


#begin plot
############
par(mfrow=c(2,2))
plot(output.mat[,1], type="l")
title("trace plot of MCMC for parameter 1")
hist(output.mat[,1], probability=T)
plot(sigma.vec, type="l", col=3, ylim=c(0, max(sigma.vec)), ylab="sigma")
plot(meanacc, type="l", col=2, ylim=c(0,1), ylab="acceptance probability")
abline(h=0.234)





