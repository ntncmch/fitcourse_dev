set_dir <- function(analysis) {

	dir_results <<- file.path(Sys.getenv("HOME"),"fitcourse",analysis)
	if(!file.exists(dir_results)){
		dir.create(dir_results)
	}

	dir_rds <<- file.path(dir_results,"rds")
	if(!file.exists(dir_rds)){
		dir.create(dir_rds)
	}

	dir_fig <<- file.path(dir_results,"figures")
	if(!file.exists(dir_fig)){
		dir.create(dir_fig)
	}

}

my_SEITL_createModelTdC <- function(deterministic=TRUE, verbose=TRUE) {

	library(plyr)
	library(truncnorm)

	# define theta using fitparam
	R0 <- fitparam(name="R0",value=10,support=c(0,Inf),sd.proposal=1,prior=list(distribution="dunif",parameters=c(min=1,max=100))) 

	LatentPeriod <- fitparam(name="LP",value=1.63,support=c(0,Inf),sd.proposal=0.2,prior=list(distribution="dtruncnorm",parameters=c(mean=1.63,sd=0.26,a=0,b=Inf))) 

	InfectiousPeriod <- fitparam(name="IP",value=1,support=c(0,Inf),sd.proposal=0.2,prior=list(distribution="dtruncnorm",parameters=c(mean=0.99,sd=0.96,a=0,b=Inf))) 

	TemporaryImmunePeriod  <- fitparam(name="TIP",value=10,support=c(0,Inf),sd.proposal=2,prior=list(distribution="dunif",parameters=c(min=0,max=50))) 

	ProbLongTermImmunity <- fitparam(name="alpha",value=0.5,support=c(0,1),sd.proposal=0.1,prior=list(distribution="dunif",parameters=c(min=0,max=1))) 

	ReportingRate <- fitparam(name="rho",value=0.7,support=c(0,Inf),sd.proposal=0.1,prior=list(distribution="dunif",parameters=c(min=0,max=2))) 

	proportionI0 <- fitparam(name="pI0",value=2/284,support=c(0,1),sd.proposal=1/284,prior=list(distribution="dunif",parameters=c(min=1/284,max=5/284))) 
	
	proportionL0 <- fitparam(name="pL0",value=0.04,support=c(0,1),sd.proposal=0.01,prior=list(distribution="dunif",parameters=c(min=0.0,max=0.1))) 
	
	PopSize <- fitparam(name="N",value=284) 

	# load and rename data
	data("FluTdC1971")
	data <- rename(FluTdC1971,c("day"="time","incidence"="Inc"))[c("time","Inc")]

	# simulator
	if(deterministic){
		simulate.model <- SEITL_simulateDeterministic
	}else{
		simulate.model <- SEITL_simulateStochastic
	}

	# create fitmodel
	SEITL <- fitmodel(
		verbose=verbose,
		name="SEITL",
		state.variables=c("S","E","I","T","L","Inc"),
		list.fitparam=list(R0,LatentPeriod,InfectiousPeriod,TemporaryImmunePeriod,ProbLongTermImmunity,ReportingRate,proportionI0,proportionL0,PopSize), 
		initialise.state=SEITL_initialiseState,
		log.prior.fitparam=SEITL_logPrior,
		simulate.model=simulate.model, 
		generate.observation=SEITL_generateObservation, 
		data=data, 
		log.likelihood=SEITL_logLikelihood,
		distance.ABC=SEITL_distanceOscillation
		) 

	return(SEITL)
}


run_MCMC_deterministic <- function() {

	SEITL <- my_SEITL_createModelTdC(deterministic=TRUE, verbose=TRUE) 

	theta.init <- SEITL$theta

	n_iteration <- 1000000
	adapt_size_start <- 50 
	adapt_size_cooling <- 0.99
	adapt_shape_start <- 100
	print_info_every <- n_iteration/10000


	# get env for replicate variable
	i_process <- as.numeric(Sys.getenv("ARG1")) + 1

	# set seed multiplicator (time difference in second since 01-12-2012) so that simulations at different time can be combined (different parameter)
	seed_mult <- as.numeric(Sys.time() - ISOdate(2012, 12, 1)) * 24 * 3600
	
	# set seed
	set.seed(i_process * seed_mult)
	
	analysis <- paste0("SEITL_deter_infoPrior_n=",n_iteration,"_size=",adapt_size_start,"_cool=",adapt_size_cooling,"_shape=",adapt_shape_start,"_run=",i_process)
	set_dir(analysis)

	# save fitmodel
	saveRDS(SEITL,file.path(dir_rds,"SEITL.rds"))

	ans <- mcmcMH(target=targetPosterior, target.args=list(log.prior=SEITL$log.prior, marginal.log.likelihood= marginalLogLikelihoodDeterministic, marginal.log.likelihood.args=list(fitmodel=SEITL)), theta.init=theta.init, gaussian.proposal=SEITL$gaussian.proposal, n.iterations=n_iteration, adapt.size.start=adapt_size_start, adapt.size.cooling=adapt_size_cooling, adapt.shape.start=adapt_shape_start, print.info.every=print_info_every)
	# save raw trace
	saveRDS(ans,file.path(dir_rds,"ans_mcmcMH.rds"))

	burn <- 0.1
	trim <- 1e4
	trace_trimed <- burnAndTrim(ans$trace,burn=burn,trim=trim)
	# save trimmed trace
	saveRDS(trace_trimed,file.path(dir_rds,paste0("trace_b=",burn,"_trim=",trim,".rds")))


	fit <- plotPosteriorFit(trace_trimed,SEITL,sample.size=300,plot=FALSE)
	# save fit
	saveRDS(fit,file.path(dir_rds,paste0("fit.rds")))

}





main <- function() {

	library(fitcourseR)

	run_MCMC_deterministic()

}

main()
