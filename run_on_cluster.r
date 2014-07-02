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

build_posterior_unifPrior <- function(stochastic=FALSE, SEIT2L=FALSE) {

	library(fitR)

	# load and rename data
	data("FluTdC1971")

	# simulator
	if(SEIT2L){
		state.init <- c("S"=279,"E"=0,"I"=2,"T1"=3,"T2"=0,"L"=0,"Inc"=0)
		if(stochastic){
			example(SEIT2L_sto)
			my_fitmodel <- SEIT2L_sto
		}else{
			example(SEIT2L_deter)
			my_fitmodel <- SEIT2L_deter
		}	
	} else {
		state.init <- c("S"=279,"E"=0,"I"=2,"T"=3,"L"=0,"Inc"=0)
		if(stochastic){
			example(SEITL_sto)
			my_fitmodel <- SEITL_sto
		}else{
			example(SEITL_deter)
			my_fitmodel <- SEITL_deter
		}	
	}
	
# trajLogLike(fitmodel=my_fitmodel, theta=theta.init, state.init=state.init, data=FluTdC1971)
# my_fitmodel$logPrior(theta.init)

	my_posterior <- function(theta){

		return(logPosterior(fitmodel=my_fitmodel, theta=theta, state.init=state.init, data=FluTdC1971, margLogLike = trajLogLike))

	}

	return(my_posterior)
}



run_MCMC_deterministic <- function() {



	n_iteration <- 100000
	adapt_size_start <- 100 
	adapt_size_cooling <- 0.99
	adapt_shape_start <- 100
	print_info_every <- n_iteration/1000


	# get env for replicate variable
	i_process <- as.numeric(Sys.getenv("ARG1")) + 1

	if(i_process%in%c(1,3)){
		theta.init <- c("R0"=2, "D.lat"=2 , "D.inf"=2, "alpha"=0.8, "D.imm"=16, "rho"=0.85)
	} else {
		theta.init <- c("R0"=20, "D.lat"=2 , "D.inf"=2, "alpha"=0.1, "D.imm"=8, "rho"=0.3)
	}

	SEIT2L <- (i_process%in%c(1,2))


	targetPosterior <- build_posterior_unifPrior(sto=FALSE,SEIT2L=SEIT2L)


	proposal.sd <- c("R0"=1, "D.lat"=0.5 , "D.inf"=0.5, "alpha"=0.1, "D.imm"=2, "rho"=0.1)
	lower <- c("R0"=0, "D.lat"=0 , "D.inf"=0, "alpha"=0, "D.imm"=0, "rho"=0)
	upper <- c("R0"=Inf, "D.lat"=Inf , "D.inf"=Inf, "alpha"=1, "D.imm"=Inf, "rho"=1)

	# set seed multiplicator (time difference in second since 01-12-2012) so that simulations at different time can be combined (different parameter)
	seed_mult <- as.numeric(Sys.time() - ISOdate(2012, 12, 1)) * 24 * 3600

	# set seed
	set.seed(i_process * seed_mult)

	analysis <- paste0(ifelse(SEIT2L,"SEIT2L","SEITL"),"_deter_unifPrior_n=",n_iteration,"_size=",adapt_size_start,"_cool=",adapt_size_cooling,"_shape=",adapt_shape_start,"_run=",i_process)
	
	set_dir(analysis)

	ans <- mcmcMH(target=targetPosterior, theta.init=theta.init, proposal.sd=proposal.sd, limits=list(lower=lower,upper=upper), n.iterations=n_iteration, adapt.size.start=adapt_size_start, adapt.size.cooling=adapt_size_cooling, adapt.shape.start=adapt_shape_start, print.info.every=print_info_every)

	# save raw trace
	saveRDS(ans,file.path(dir_rds,"ans_mcmcMH.rds"))

	# trace <- mcmc(ans$trace)
	# xyplot(x=trace)
	# plotESSBurn(ans$trace,longest.burn.in=1200)

	# burn_until <- 800
	# xyplot(x=trace, start=burn_until)
	# acfplot(x=trace, start=burn_until, lag.max=100)
	# keep_every <- 10
	# xyplot(x=trace, start=burn_until, thin=keep_every)
	
	# trace <- mcmc(burnAndThin(ans$trace, burn=burn_until, thin=keep_every))

	# densityplot(x=trace)
	# summary(trace)

	# levelplot(x=trace)
	# crosscorr.plot(x=trace)

	# effectiveSize(trace)
	# HPDinterval(trace)
	# rejectionRate(trace)

}





main <- function() {

	library(fitR)

	run_MCMC_deterministic()

}

main()
