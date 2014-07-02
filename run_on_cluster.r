set_dir <- function(dir_name) {

	dir_results <<- file.path(Sys.getenv("HOME"),"fitcourse",dir_name)
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

build_posterior <- function(stochastic=FALSE, SEIT2L=FALSE, priorInfo=FALSE) {

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
	
	if(priorInfo){

		my_fitmodel$logPrior <- function(theta) {

			require(truncnorm)

			log.prior.R0 <- dunif(theta[["R0"]], min = 1, max = 50, log = TRUE)
			log.prior.latent.period <- log(dtruncnorm(theta[["D.lat"]], a = 0, b = Inf, mean = 2, sd = 1))
			log.prior.infectious.period <- log(dtruncnorm(theta[["D.inf"]], a = 0, b = Inf, mean = 2, sd = 1))
			log.prior.temporary.immune.period <- dunif(theta[["D.imm"]], min = 0, max = 50, log = TRUE)
			log.prior.probability.long.term.immunity <- dunif(theta[["alpha"]], min = 0, max = 1, log = TRUE)
			log.prior.reporting.rate <- dunif(theta[["rho"]], min = 0, max = 1, log = TRUE)

			return(log.prior.R0 + log.prior.latent.period + log.prior.infectious.period + log.prior.temporary.immune.period + log.prior.probability.long.term.immunity + log.prior.reporting.rate)

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



	adapt_size_start <- 100 
	adapt_size_cooling <- 0.999
	adapt_shape_start <- 200


	# get env for replicate variable
	i_process <- as.numeric(Sys.getenv("ARG1")) + 1

	df_set <- expand.grid(theta=1:2,priorInfo=c(FALSE,TRUE),SEIT2L=c(FALSE,TRUE),n_iteration=c(5000,100000))

	df_set <- df_set[i_process,]

	n_iteration <- df_set$n_iteration
	print_info_every <- n_iteration/1000


	theta <- list(
		theta1=c("R0"=2, "D.lat"=2 , "D.inf"=2, "alpha"=0.8, "D.imm"=16, "rho"=0.85),
		theta2=c("R0"=20, "D.lat"=2 , "D.inf"=2, "alpha"=0.1, "D.imm"=8, "rho"=0.3)
		)

	theta.init <- theta[[df_set$theta]]
	
	targetPosterior <- build_posterior(sto=FALSE,SEIT2L=df_set$SEIT2L,priorInfo=df_set$priorInfo)

	proposal.sd <- c("R0"=1, "D.lat"=0.5 , "D.inf"=0.5, "alpha"=0.1, "D.imm"=2, "rho"=0.1)
	lower <- c("R0"=0, "D.lat"=0 , "D.inf"=0, "alpha"=0, "D.imm"=0, "rho"=0)
	upper <- c("R0"=Inf, "D.lat"=Inf , "D.inf"=Inf, "alpha"=1, "D.imm"=Inf, "rho"=1)

	# set seed multiplicator (time difference in second since 01-12-2012) so that simulations at different time can be combined (different parameter)
	seed_mult <- as.numeric(Sys.time() - ISOdate(2012, 12, 1)) * 24 * 3600

	# set seed
	set.seed(i_process * seed_mult)

	analysis <- paste0(ifelse(df_set$SEIT2L,"SEIT2L","SEITL"),"_deter_",ifelse(df_set$priorInfo,"info","unif"),"Prior_n=",n_iteration,"_size=",adapt_size_start,"_cool=",adapt_size_cooling,"_shape=",adapt_shape_start,"_set=",i_process)
	dir_name <- "mcmc_deter"
	set_dir(dir_name)

	ans <- mcmcMH(target=targetPosterior, theta.init=theta.init, proposal.sd=proposal.sd, limits=list(lower=lower,upper=upper), n.iterations=n_iteration, adapt.size.start=adapt_size_start, adapt.size.cooling=adapt_size_cooling, adapt.shape.start=adapt_shape_start, print.info.every=print_info_every)

	# saveRDS(ans,paste0("/Users/Tonton/edu/Fit_course/dev/dataset/",analysis,".rds"))

	# save raw trace
	saveRDS(ans,file.path(dir_rds,paste0("mcmc_",analysis,".rds")))

	# ans <- readRDS(paste0("/Users/Tonton/edu/Fit_course/dev/dataset/SEITL_deter_unifPrior_n=5000_size=100_cool=0.99_shape=100_run=4.rds"))

	# trace <- mcmc(ans$trace)
	# xyplot(x=trace)
	# plotESSBurn(ans$trace,longest.burn.in=2000)

	# burn_until <- 1000
	# xyplot(x=trace, start=burn_until)
	# acfplot(x=trace, start=burn_until, lag.max=100)
	# keep_every <- 30
	# xyplot(x=trace, start=burn_until, thin=keep_every)
	
	# trace <- mcmc(burnAndThin(ans$trace, burn=burn_until, thin=keep_every))

	# densityplot(x=trace)
	# summary(trace)

	# levelplot(x=trace)
	# crosscorr.plot(x=trace)

	# effectiveSize(trace)
	# HPDinterval(trace)
	# rejectionRate(mcmc(ans$trace))

}

main <- function() {

	library(fitR)

	run_MCMC_deterministic()

}

main()
