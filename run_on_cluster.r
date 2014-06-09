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

my_SEIT2L_createModelTdC <- function(deterministic=TRUE, verbose=TRUE) {

	library(plyr)
	library(truncnorm)

	# define theta using fitparam
	R0 <- fitparam(name="R0",value=10,support=c(0,Inf),sd.proposal=1,prior=list(distribution="dunif",parameters=c(min=1,max=100))) 
	
	LatentPeriod <- fitparam(name="LP",value=1.63,support=c(0,Inf),sd.proposal=0.2,prior=list(distribution="dtruncnorm",parameters=c(mean=1.63,sd=0.26,a=0,b=Inf))) 

	InfectiousPeriod <- fitparam(name="IP",value=1,support=c(0,Inf),sd.proposal=0.2,prior=list(distribution="dtruncnorm",parameters=c(mean=0.99,sd=0.96,a=0,b=Inf))) 

	TemporaryImmunePeriod  <- fitparam(name="TIP",value=10,support=c(0,Inf),sd.proposal=2,prior=list(distribution="dunif",parameters=c(min=0,max=50))) 

	ProbLongTermImmunity <- fitparam(name="alpha",value=0.5,support=c(0,1),sd.proposal=0.1,prior=list(distribution="dunif",parameters=c(min=0,max=1))) 

	ReportingRate <- fitparam(name="rho",value=0.7,support=c(0,2),sd.proposal=0.1,prior=list(distribution="dunif",parameters=c(min=0,max=2))) 

	proportionI0 <- fitparam(name="pI0",value=2/284,support=c(0,1),sd.proposal=1/284,prior=list(distribution="dunif",parameters=c(min=1/284,max=5/284))) 
	
	proportionL0 <- fitparam(name="pL0",value=0.1,support=c(0,1),sd.proposal=0.01,prior=list(distribution="dunif",parameters=c(min=0.0,max=0.5))) 
	
	PopSize <- fitparam(name="N",value=284) 

	# load and rename data
	data("FluTdC1971",envir = environment())
	data <- rename(FluTdC1971,c("day"="time","incidence"="Inc"))[c("time","Inc")]

	# simulator
	if(deterministic){
		simulate.model <- SEIT2L_simulateDeterministic
	}else{
		simulate.model <- SEIT2L_simulateStochastic
	}

	# create fitmodel
	SEIT2L <- fitmodel(
		verbose=verbose,
		name="SEIT2L",
		state.variables=c("S","E","I","T1","T2","L","Inc"),
		list.fitparam=list(R0,LatentPeriod,InfectiousPeriod,TemporaryImmunePeriod,ProbLongTermImmunity,ReportingRate,proportionI0,proportionL0,PopSize), 
		initialise.state=SEIT2L_initialiseState,
		log.prior.fitparam=SEIT2L_logPrior,
		simulate.model=simulate.model, 
		generate.observation=SEIT2L_generateObservation, 
		data=data, 
		log.likelihood=SEIT2L_logLikelihood,
		distance.ABC=SEIT2L_distanceOscillation
		) 

	return(SEIT2L)
}



run_MCMC_deterministic <- function() {

	SEIT2L <- my_SEIT2L_createModelTdC(deterministic=TRUE, verbose=TRUE) 

	theta.init <- SEIT2L$theta

	n_iteration <- 500000
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
	
	analysis <- paste0("SEIT2L_deter_infoPrior_n=",n_iteration,"_size=",adapt_size_start,"_cool=",adapt_size_cooling,"_shape=",adapt_shape_start,"_run=",i_process)
	set_dir(analysis)

	# save fitmodel
	saveRDS(SEIT2L,file.path(dir_rds,"SEIT2L.rds"))

	ans <- mcmcMH(target=targetPosterior, target.args=list(log.prior=SEIT2L$log.prior, marginal.log.likelihood= marginalLogLikelihoodDeterministic, marginal.log.likelihood.args=list(fitmodel=SEIT2L)), theta.init=theta.init, gaussian.proposal=SEIT2L$gaussian.proposal, n.iterations=n_iteration, adapt.size.start=adapt_size_start, adapt.size.cooling=adapt_size_cooling, adapt.shape.start=adapt_shape_start, print.info.every=print_info_every)
	# save raw trace
	saveRDS(ans,file.path(dir_rds,"ans_mcmcMH.rds"))

	burn <- 0.1
	trim <- 1e4
	trace_trimed <- burnAndTrim(ans$trace,burn=burn,trim=trim)
	# save trimmed trace
	saveRDS(trace_trimed,file.path(dir_rds,paste0("trace_b=",burn,"_trim=",trim,".rds")))


	fit <- plotPosteriorFit(trace_trimed,SEIT2L,sample.size=300,plot=FALSE)
	# save fit
	saveRDS(fit,file.path(dir_rds,paste0("fit.rds")))

}





main <- function() {

	library(fitcourseR)

	run_MCMC_deterministic()

}

main()
