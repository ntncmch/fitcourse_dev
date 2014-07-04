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

build_posterior <- function(stochastic=FALSE, SEIT2L=FALSE, priorInfo=FALSE, n_particles=48) {

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

	if(stochastic){
		my_posterior <- function(theta){

			return(logPosterior(fitmodel=my_fitmodel, theta=theta, state.init=state.init, data=FluTdC1971, margLogLike = margLogLikeSto, n.particles=n_particles, n.cores=NULL))

		}
	} else {
		my_posterior <- function(theta){

			return(logPosterior(fitmodel=my_fitmodel, theta=theta, state.init=state.init, data=FluTdC1971, margLogLike = trajLogLike))

		}		
	}

	return(my_posterior)
}

test_smc <- function(n_iter=10, n_particles=c(120,240)) {

	n.cores <- detectCores()
	cat("SMC runs on ",n.cores," cores\n")
	flush.console()
	
	library(fitR)

	theta <- c("R0"=6.4, "D.lat"=1.3 , "D.inf"=2.8, "alpha"=0.49, "D.imm"=10, "rho"=0.67)

	i_process <- as.numeric(Sys.getenv("ARG1")) + 1

	n_particles <- n_particles[i_process]

	example(SEIT2L_sto)

	state.init <- c("S"=279,"E"=0,"I"=2,"T1"=3,"T2"=0,"L"=0,"Inc"=0)
	data(FluTdC1971)

	start_smc  <- Sys.time()
	sample_ll <- vector("numeric",length=n_iter)
	for(i in seq_along(sample_ll)){
		sample_ll[i] <- margLogLikeSto(fitmodel=SEIT2L_sto, theta=theta, state.init=state.init, data=FluTdC1971, n.particles=n_particles, n.cores=NULL)
	}

	end_smc  <- Sys.time()
	suppressMessages(time.estimation <- round(as.period((end_smc-start_smc)*10000/length(sample_ll))))		

	sample_finite_ll <- sample_ll[is.finite(sample_ll)]

	ans <- list(ll=sample_ll, stat=c(mean=mean(sample_finite_ll), sd=sd(sample_finite_ll), prop_finite=length(sample_finite_ll)/length(sample_ll)),time=time.estimation)

	print(ans)

	set_dir("test_smc")	
	name <- paste0(n_particles,"_particles.rds")
	saveRDS(ans,file.path(dir_rds,name))

}

analyse_smc <- function(n_particles) {

	name <- paste0(n_particles,"_particles.rds")
	dir_rds <- "/Users/Tonton/edu/Fit_course/dev/dataset/test_smc/rds"
	list_ans <- lapply(name,function(x) {readRDS(file.path(dir_rds,x))})
	names(list_ans) <- n_particles
	# stat <- ldply(list_ans,function(x) {as.data.frame(x[c("mean","sd")])})
	stat <- ldply(list_ans,function(x) {t(data.frame(x$stat))},.id="n_particles")
	time <- ldply(list_ans,function(x) {data.frame(time10000=as.numeric(as.duration(x[["time"]])))},.id="n_particles")
	df_bench <- join(stat,time)
	df_bench <- mutate(df_bench,n_particles=as.numeric(as.character(n_particles)),time10000_day=time10000/3600/24)

	df_plot <- melt(df_bench,id.vars="n_particles")
	ggplot(df_plot, aes(x=n_particles,y=value))+facet_wrap(~variable,scales="free_y")+geom_line()+geom_vline(xintercept=408,col="red")

}

run_MCMC <- function(stochastic=FALSE) {


	if(stochastic){
		adapt_size_start <- 50 
		adapt_size_cooling <- 0.99
		adapt_shape_start <- 100
	} else {
		adapt_size_start <- 100 
		adapt_size_cooling <- 0.999
		adapt_shape_start <- 200	
	}
	


	# get env for replicate variable
	i_process <- as.numeric(Sys.getenv("ARG1")) + 1
	# set seed multiplicator (time difference in second since 01-12-2012) so that simulations at different time can be combined (different parameter)
	seed_mult <- as.numeric(Sys.time() - ISOdate(2012, 12, 1)) * 24 * 3600

	# set seed
	set.seed(i_process * seed_mult)


	if(stochastic){
		df_set <- expand.grid(theta=1,priorInfo=TRUE,SEIT2L=TRUE,n_iteration=3000)
	} else {
		df_set <- expand.grid(theta=1:2,priorInfo=c(FALSE,TRUE),SEIT2L=c(FALSE,TRUE),n_iteration=c(5000,100000))		
		df_set <- df_set[i_process,]
	}	

	n_iteration <- df_set$n_iteration
	print_info_every <- 1 #n_iteration/1000


	if(stochastic){

		theta.init <- c(
			"R0"=rnorm(1,mean=6.4,sd=1.52579/10),
			"D.lat"=rnorm(1,mean=1.3,sd=0.27554/10),
			"D.inf"=rnorm(1,mean=2.8,sd=0.80270/10), 
			"alpha"=rnorm(1,mean=0.49,sd=0.04039/10), 
			"D.imm"=rnorm(1,mean=10,sd=1.11616/10), 
			"rho"=rnorm(1,mean=0.67,sd=0.04671/10)
			)


	} else {

		theta <- list(
			theta1=c("R0"=2, "D.lat"=2 , "D.inf"=2, "alpha"=0.8, "D.imm"=16, "rho"=0.85),
			theta2=c("R0"=20, "D.lat"=2 , "D.inf"=2, "alpha"=0.1, "D.imm"=8, "rho"=0.3)
			)

		theta.init <- theta[[df_set$theta]]

	}

	targetPosterior <- build_posterior(stochastic=stochastic,SEIT2L=df_set$SEIT2L,priorInfo=df_set$priorInfo, n_particles=408)

	if(stochastic){
		data(mcmc_TdC_deter_longRun)
		covmat <- mcmc_SEIT2L_infoPrior_theta1$covmat.empirical		
		proposal.sd <- NULL
	} else {
		proposal.sd <- c("R0"=1, "D.lat"=0.5 , "D.inf"=0.5, "alpha"=0.1, "D.imm"=2, "rho"=0.1)	
		covmat <- NULL	
	}

	lower <- c("R0"=0, "D.lat"=0 , "D.inf"=0, "alpha"=0, "D.imm"=0, "rho"=0)
	upper <- c("R0"=Inf, "D.lat"=Inf , "D.inf"=Inf, "alpha"=1, "D.imm"=Inf, "rho"=1)

	
	analysis <- paste0(ifelse(df_set$SEIT2L,"SEIT2L","SEITL"),"_",ifelse(stochastic,"sto","deter"),"_",ifelse(df_set$priorInfo,"info","unif"),"Prior_n=",n_iteration,"_size=",adapt_size_start,"_cool=",adapt_size_cooling,"_shape=",adapt_shape_start,"_set=",i_process)
	dir_name <- ifelse(stochastic,"mcmc_sto_2","mcmc_deter")
	set_dir(dir_name)


	n.cores <- detectCores()
	cat("SMC runs on ",n.cores," cores\n")
	flush.console()
	print(theta.init)


	ans <- mcmcMH(target=targetPosterior, theta.init=theta.init, proposal.sd=proposal.sd, covmat=covmat, limits=list(lower=lower,upper=upper), n.iterations=n_iteration, adapt.size.start=adapt_size_start, adapt.size.cooling=adapt_size_cooling, adapt.shape.start=adapt_shape_start, print.info.every=print_info_every)

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

analyse_mcmc <- function() {


	df_set <- expand.grid(theta=1:2,priorInfo=c(FALSE,TRUE),SEIT2L=c(FALSE,TRUE),n_iteration=c(5000,100000))
	dir_name <- "/Users/Tonton/edu/Fit_course/dev/dataset/mcmc_deter"
	dir_rds <- file.path(dir_name,"rds")
	dir_fig <- file.path(dir_name,"figures")

	adapt_size_start <- 100 
	adapt_size_cooling <- 0.999
	adapt_shape_start <- 200


	i <- 13:14
	df <- df_set[i,]
	analysis <- paste0("mcmc_",ifelse(df$SEIT2L,"SEIT2L","SEITL"),"_deter_",ifelse(df$priorInfo,"info","unif"),"Prior_n=",df$n_iteration,"_size=",adapt_size_start,"_cool=",adapt_size_cooling,"_shape=",adapt_shape_start,"_set=",i,".rds")

	if(length(analysis)==1){

		ans <- readRDS(file.path(dir_rds,analysis))
		trace <- mcmc(ans$trace)

	} else {

		ans <- lapply(analysis,function(x) {readRDS(file.path(dir_rds,x))})
		trace <- lapply(ans,function(x) {mcmc(x$trace)})

	}

	library(fitR)

	# trace <- mcmc(ans[[1]]$trace)

	trace <- mcmc.list(trace)
	effectiveSize(trace)
	rejectionRate(trace)
	xyplot(x=trace)
	plotESSBurn(ans[[2]]$trace[1:5000,],longest.burn.in=1000)

	burn_until <- 1000
	xyplot(x=trace, start=burn_until)
	acfplot(x=trace, start=burn_until, lag.max=100)
	keep_every <- 50
	xyplot(x=trace, start=burn_until, thin=keep_every)
	xyplot(x=trace,thin=keep_every)
	
	trace2 <- burnAndThin(trace, burn=burn_until, thin=keep_every)
	effectiveSize(trace2)
	densityplot(x=trace2)
	summary(trace2)

	levelplot(x=trace2)
	crosscorr.plot(x=trace)

	densityplot(x=trace2)

	# 
	# HPDinterval(trace)
	# 


	

	

}



main <- function() {

	library(fitR)

	# n_particles <- 12*c(seq(4,30,4),seq(34,88,8))
	# test_smc(n_iter=100,n_particles=n_particles)

	run_MCMC(stochastic=TRUE)
	
}

main()
