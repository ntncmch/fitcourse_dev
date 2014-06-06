set_dir <- function(analysis) {

	dir_results <<- file.path(Sys.getenv("HOME"),"fitcourse",analysis)
	if(!file.exists(dir_results)){
		file.create(dir_results)
	}

	dir_rds <<- file.path(dir_results,"rds")
	if(!file.exists(dir_rds)){
		file.create(dir_rds)
	}

	dir_fig <<- file.path(dir_results,"figures")
	if(!file.exists(dir_fig)){
		file.create(dir_fig)
	}

}


run_MCMC_deterministic <- function() {

	SEITL <- SEITL_createModelTdC(deterministic=TRUE, verbose=TRUE) 

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
	
	analysis <- paste0("SEITL_deter_n=",n_iteration,"_size=",adapt_size_start,"_cool=",adapt_size_cooling,"_shape=",adapt_shape_start,"_run=",i_process)
	set_dir(analysis)

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
