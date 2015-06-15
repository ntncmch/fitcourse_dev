
start_me <- function() {

	dir_dev <<- "/Users/Tonton/edu/Fit_course/dev"
	dir_pkg <<- "/Users/Tonton/edu/Fit_course/fitR"
	dir_knitr <<- "/Users/Tonton/edu/Fit_course/mfiidd/knitr"
	dir_md <<- "/Users/Tonton/edu/Fit_course/mfiidd"

	# setwd(dir_dev)

	library(jsonlite)
	library(testthat)
	library(plyr)
	library(dplyr)
	library(lubridate)
	library(ggplot2)
	library(stringr)
	library(reshape2)
	library(devtools)
	library(profr)
	library(deSolve)
	library(adaptivetau)
	library(compiler)
	library(rlecuyer)
	library(mvtnorm)
	library(tmvtnorm)


}

install_fitR <- function() {

	install_github("sbfnk/fitR")

}


create_data <- function() {

	FluTdC1971 <- readRDS("/Users/Tonton/edu/Fit_course/dev/dataset/TdC_1971.rds")
	FluTdC1971 <- rename(FluTdC1971,c("day"="time","incidence"="obs"))#[c("time","Inc")]
# str(FluTdC1971)
	save(FluTdC1971,file=file.path(dir_pkg,"data","FluTdC1971.rdata"))

}

create_SMC_benchmark <- function() {

	n_particles <- 12*c(seq(4,30,4),seq(34,88,8))
	name <- paste0(n_particles,"_particles.rds")
	dir_rds <- "/Users/Tonton/edu/Fit_course/dev/dataset/test_smc/rds"
	list_ans <- lapply(name,function(x) {readRDS(file.path(dir_rds,x))})
	names(list_ans) <- n_particles
	# stat <- ldply(list_ans,function(x) {as.data.frame(x[c("mean","sd")])})
	stat <- ldply(list_ans,function(x) {t(data.frame(x$stat))},.id="n_particles")
	time <- ldply(list_ans,function(x) {data.frame(time10000=as.numeric(as.duration(x[["time"]])))},.id="n_particles")
	df_bench <- join(stat,time)
	df_bench <- mutate(df_bench,n_particles=as.numeric(as.character(n_particles)),time_10000iter_day=time10000/3600/24, prop_depleted=1-prop_finite)
	df_bench$time10000 <- NULL
	df_bench$prop_finite <- NULL

	calibrateSMC <- df_bench
	save(calibrateSMC,file=file.path(dir_pkg,"data","calibrateSMC.rdata"))
	

}

create_trace_mcmc_deter <- function() {



	df_set <- 	expand.grid(theta=1:2,priorInfo=c(FALSE,TRUE),model=c("SEITL","SEIT2L"),n_iteration=c(5000,100000)) 
	df_set$set <- 1:nrow(df_set)
	# Fix for SEIT4L as I resimulated only SEIT4L and not SEITL so "set" is not appropriate
	# df_set <- df_set %>% filter(!SEIT4L)

	df_set_SEITL4 <- expand.grid(theta=1:2,priorInfo=c(FALSE,TRUE),model="SEIT4L",n_iteration=c(5000,100000))
	df_set_SEITL4$set <- 1:nrow(df_set_SEITL4)
	df_set <- rbind(df_set, df_set_SEITL4)

	dir_name <- "/Users/Tonton/edu/Fit_course/dev/dataset/mcmc_deter"
	dir_rds <- file.path(dir_name,"rds")
	dir_fig <- file.path(dir_name,"figures")

	adapt_size_start <- 100 
	adapt_size_cooling <- 0.999
	adapt_shape_start <- 200

	# short trace
	df_prel <- subset(df_set, n_iteration== 5000 & !priorInfo & theta==1)

	list_trace <- dlply(df_prel,c("theta","priorInfo","model","n_iteration"),function(df) {

		analysis <- paste0("mcmc_",df$model,"_deter_",ifelse(df$priorInfo,"info","unif"),"Prior_n=",df$n_iteration,"_size=",adapt_size_start,"_cool=",adapt_size_cooling,"_shape=",adapt_shape_start,"_set=",df$set,".rds")

		ans <- readRDS(file.path(dir_rds,analysis))

		return(ans)
	})
	
	name <- sapply(df_prel$model, function(x) {paste0("mcmc_",x)})
	
	names(list_trace) <- name

	# change names of SEITL and SEIT2L (since 2015 we use _ in param names instead of .)
	change_names <- name[!str_detect(name, "SEIT4L")]
	keep_name <- setdiff(name, change_names)
	for(x in change_names){
		names(list_trace[[x]]$trace) <- names(list_trace[[keep_name]]$trace)
		rownames(list_trace[[x]]$covmat.empirical) <- rownames(list_trace[[keep_name]]$covmat.empirical)
		colnames(list_trace[[x]]$covmat.empirical) <- colnames(list_trace[[keep_name]]$covmat.empirical)
	}


	attach(list_trace)

	save(list=name,file=file.path(dir_pkg,"data","mcmc_TdC_deter_shortRun.rdata"))

	# data("mcmc_TdC_deter_longRun",envir = environment())
# 	ls()
	# names(mcmc_SEITL_theta1)

	# long trace
	df_long <- subset(df_set,n_iteration>5000)

	list_trace_name <- dlply(df_long,c("theta","priorInfo","model","n_iteration"),function(df) {

		analysis <- paste0("mcmc_",df$model,"_deter_",ifelse(df$priorInfo,"info","unif"),"Prior_n=",df$n_iteration,"_size=",adapt_size_start,"_cool=",adapt_size_cooling,"_shape=",adapt_shape_start,"_set=",df$set,".rds")

		ans <- readRDS(file.path(dir_rds,analysis))

		name <- paste0("mcmc_",df$model,ifelse(df$priorInfo,"_infoPrior_","_"),"theta",df$theta)


		return(list(ans,name))
	})

	list_trace <- sapply(list_trace_name,function(x) {x[1]})
	list_name <- sapply(list_trace_name,function(x) {x[[2]]})

	names(list_trace) <- list_name

	change_names <- list_name[!str_detect(list_name, "SEIT4L")]
	keep_name <- setdiff(list_name, change_names) %>% first
	
	for(x in change_names){
		names(list_trace[[x]]$trace) <- names(list_trace[[keep_name]]$trace)
		rownames(list_trace[[x]]$covmat.empirical) <- rownames(list_trace[[keep_name]]$covmat.empirical)
		colnames(list_trace[[x]]$covmat.empirical) <- colnames(list_trace[[keep_name]]$covmat.empirical)
	}

	attach(list_trace)

	save(list=list_name,file=file.path(dir_pkg,"data","mcmc_TdC_deter_longRun.rdata"))


}

create_trace_mcmc_sto <- function() {

	dir_name <- "/Users/Tonton/edu/Fit_course/dev/dataset/mcmc_sto_SEIT2L_n408"
	dir_rds <- file.path(dir_name,"rds")
	dir_fig <- file.path(dir_name,"figures")

	analysis <- paste0("mcmc_SEIT2L_sto_infoPrior_n=3000_size=50_cool=0.99_shape=100_set=",1:5,".rds")

	pmcmc_SEIT2L_infoPrior_n400 <- lapply(analysis,function(x) {readRDS(file.path(dir_rds,x))})
	names(pmcmc_SEIT2L_infoPrior_n400) <- paste0("chain",seq_along(pmcmc_SEIT2L_infoPrior_n400))
	save(pmcmc_SEIT2L_infoPrior_n400,file=file.path(dir_pkg,"data","pmcmc_SEIT2L_infoPrior_n400.rdata"))
	


	dir_name <- "/Users/Tonton/edu/Fit_course/dev/dataset/mcmc_sto_SEIT2L_n48"
	dir_rds <- file.path(dir_name,"rds")
	dir_fig <- file.path(dir_name,"figures")

	analysis <- paste0("mcmc_SEIT2L_sto_infoPrior_n=3000_size=50_cool=0.99_shape=100_set=",1:5,".rds")

	pmcmc_SEIT2L_infoPrior_n50 <- lapply(analysis,function(x) {readRDS(file.path(dir_rds,x))})

	x <- pmcmc_SEIT2L_infoPrior_n50[[5]]
	x$trace[1,] <- x$trace[2,]
	pmcmc_SEIT2L_infoPrior_n50[[5]] <- x

	names(pmcmc_SEIT2L_infoPrior_n50) <- paste0("chain",seq_along(pmcmc_SEIT2L_infoPrior_n50))

	save(pmcmc_SEIT2L_infoPrior_n50,file=file.path(dir_pkg,"data","pmcmc_SEIT2L_infoPrior_n50.rdata"))
	


	# trace <- lapply(ans,function(x) {mcmc(x$trace)})

	# trace <- mcmc.list(trace)

	# 1-rejectionRate(trace)

	# effectiveSize(trace)

	# xyplot(trace)

	# plotESSBurn(trace[[1]])

	# burn_until <- 250
	# trace.burn <- burnAndThin(trace, burn=burn_until)
	# effectiveSize(trace.burn)

	# xyplot(x=trace.burn)

	# acfplot(x=trace, start=burn_until, lag.max=50)
	# keep_every <- 10

	# trace.burn.thin <- burnAndThin(trace.burn, thin=keep_every)
	# xyplot(x=trace.burn.thin)

	# effectiveSize(trace.burn.thin)
	# densityplot(x=trace.burn)
	# densityplot(x=trace.burn.thin)
	# summary(trace.burn.thin)



	# data(SEIT2L_stoch)
	# state.init <- c("S"=279,"E"=0,"I"=2,"T1"=3,"T2"=0,"L"=0,"Inc"=0)

	# plotPosteriorFit(trace.burn.thin,SEIT2L_stoch,state.init,data=FluTdC1971,posterior.summary="mean")
	# quartz()
	# plotPosteriorFit(trace.burn.thin,SEIT2L_stoch,state.init,data=FluTdC1971,posterior.summary="sample")

	# levelplot(x=trace2)
	# crosscorr.plot(x=trace)


	# # compare with deter

	# # create mcmc object
	# trace.info1 <- mcmc(mcmc_SEIT2L_infoPrior_theta1$trace)
	# trace.info2 <- mcmc(mcmc_SEIT2L_infoPrior_theta2$trace)

	# # combine in a mcmc.list
	# trace.info <- mcmc.list(trace.info1,trace.info2)

	# # burn and thin as the chain with uniform prior (see above sections)
	# trace.info.burn.thin <-  burnAndThin(trace.info, burn=5000, thin=40)

	# plotPosteriorDensity(list(deter=trace.info.burn.thin,sto=trace.burn.thin))


	# x <- summary(trace.info.burn.thin)
	# theta.deter <- x$statistics[SEIT2L_stoch$theta.names,"Mean"]


	# x <- summary(trace.burn.thin)
	# theta.sto <- x$statistics[SEIT2L_stoch$theta.names,"Mean"]

	# plotFit(SEIT2L_stoch,theta.sto,init.state,data=FluTdC1971,n=100,p.extinction=TRUE)
	# quartz()
	# plotFit(SEIT2L_stoch,theta.deter,init.state,data=FluTdC1971,n=1000,p.extinction=TRUE)

}

create_data_fitmodel_SIR <- function() {


# create a simple deterministic SIR model with constant population size

	SIR_name <- "SIR with constant population size"
	SIR_state.names <- c("S","I","R")
	SIR_theta.names <- c("R0","D.inf")

	SIR_simulateDeterministic <- function(theta,init.state,times) {

		SIR_ode <- function(time, state, parameters) {

            ## parameters
			beta <- parameters[["R0"]] / parameters[["D.inf"]]
			nu <- 1 / parameters[["D.inf"]]

            ## states
			S <- state[["S"]]
			I <- state[["I"]]
			R <- state[["R"]]

			N <- S + I + R

			dS <- -beta * S * I/N
			dI <- beta * S * I/N - nu * I
			dR <- nu * I

			return(list(c(dS, dI, dR)))
		}

		trajectory <- data.frame(ode(y=init.state,times=times,func=SIR_ode,parms=theta, method = "ode45"))

		return(trajectory)
	}

	## function to compute log-prior
	SIR_logPrior <- function(theta) {

        ## uniform prior on R0: U[1,100]
		log.prior.R0 <- dunif(theta[["R0"]], min = 1, max = 100, log = TRUE)
        ## uniform prior on infectious period: U[0,30]
		log.prior.D <- dunif(theta[["D.inf"]], min = 0, max = 30, log = TRUE)

		return(log.prior.R0 + log.prior.D)
	}

	## function to compute the log-likelihood of one data point
	SIR_pointLogLike <- function(data.point, model.point, theta){

        ## the prevalence is observed through a Poisson process
		return(dpois(x=data.point[["obs"]], lambda=model.point[["I"]], log=TRUE))
	}

	## function to generate observation from a model simulation
	SIR_genObsPoint <- function(model.point, theta){

        ## the prevalence is observed through a Poisson process
		obs.point <- rpois(n=1, lambda=model.point[["I"]])

		return(c(obs=obs.point))
	}

	## create deterministic SIR fitmodel
	SIR <- fitmodel(
		name=SIR_name,
		state.names=SIR_state.names,
		theta.names=SIR_theta.names,
		simulate=SIR_simulateDeterministic,
		genObsPoint=SIR_genObsPoint,
		logPrior=SIR_logPrior,
		pointLogLike=SIR_pointLogLike)

	save(SIR,file=file.path(dir_pkg,"data","SIR.rdata"))


	SIR_stochastic_name <- "stochastic SIR with constant population size"

	SIR_simulateStochastic <- function(theta,init.state,times) {

        ## transitions
		SIR_transitions <- list(
			c(S = -1, I = 1), # infection
			c(I = -1, R = 1) # recovery
			)

        ## rates
		SIR_rateFunc <- function(x, parameters, t) {

			beta <- parameters[["R0"]]/parameters[["D"]]
			nu <- 1/parameters[["D"]]

			S <- x[["S"]]
			I <- x[["I"]]
			R <- x[["R"]]

			N <- S + I + R

			return(c(
				beta * S * I / N, # infection
				nu * I # recovery
				))
		}

        # make use of the function simulateModelStochastic that
        # returns trajectories in the correct format
		return(simulateModelStochastic(theta,init.state,times,SIR_transitions,SIR_rateFunc))

	}

	# create stochastic SIR fitmodel
	SIR_stoch <- fitmodel(
		name=SIR_stochastic_name,
		state.names=SIR_state.names,
		theta.names=SIR_theta.names,
		simulate=SIR_simulateStochastic,
		genObsPoint=SIR_genObsPoint,
		logPrior=SIR_logPrior,
		pointLogLike=SIR_pointLogLike)

	save(SIR_stoch,file=file.path(dir_pkg,"data","SIR_stoch.rdata"))


	SIR_reporting_name <- "SIR with constant population size and incomplete reporting"
	SIR_reporting_theta.names <- SIR_theta.names <- c("R0","D.inf", "RR")

	## function to compute log-prior
	SIR_logPrior <- function(theta) {

        ## uniform prior on R0: U[1,100]
		log.prior.R0 <- dunif(theta[["R0"]], min = 1, max = 100, log = TRUE)
        ## uniform prior on infectious period: U[0,30]
		log.prior.D.inf <- dunif(theta[["D.inf"]], min = 0, max = 30, log = TRUE)
        ## uniform prior on the reporting rate: U[0,1]
		log.prior.RR <- dunif(theta[["RR"]], min = 0, max = 1, log = TRUE)

		return(log.prior.R0 + log.prior.D.inf + log.prior.RR)
	}

	## function to compute the log-likelihood of one data point
	SIR_reporting_pointLogLike <- function(data.point, model.point, theta){

        ## the prevalence is observed through a Poisson process with a reporting rate
		return(dpois(x=data.point[["obs"]], lambda=model.point[["I"]]*theta[["RR"]], log=TRUE))
	}

	## function to generate observation from a model simulation
	SIR_reporting_genObsPoint <- function(model.point, theta){

       ## the prevalence is observed through a Poisson process
		obs.point <- rpois(n=1, lambda=model.point[["I"]]*theta[["RR"]])

		return(c(obs=obs.point))
	}

	## create deterministic SIR fitmodel
	SIR_reporting <- fitmodel(
		name=SIR_reporting_name,
		state.names=SIR_state.names,
		theta.names=SIR_theta.names,
		simulate=SIR_simulateDeterministic,
		genObsPoint=SIR_reporting_genObsPoint,
		logPrior=SIR_logPrior,
		pointLogLike=SIR_reporting_pointLogLike)

	save(SIR_reporting,file=file.path(dir_pkg,"data","SIR_reporting.rdata"))


	SIR_exp_name <- "SIR with constant population size, parameters transformed to the exponential scale"

	SIR_exp_simulateDeterministic <- function(theta,init.state,times) {

		SIR_ode <- function(time, state, parameters) {

                ## parameters
			beta <- exp(parameters[["R0"]]) / exp(parameters[["D.inf"]])
			nu <- 1 / exp(parameters[["D.inf"]])

                ## states
			S <- state[["S"]]
			I <- state[["I"]]
			R <- state[["R"]]

			N <- S + I + R

			dS <- -beta * S * I/N
			dI <- beta * S * I/N - nu * I
			dR <- nu * I

			return(list(c(dS, dI, dR)))
		}

		trajectory <- data.frame(ode(y=init.state,times=times,func=SIR_ode,parms=theta, method = "ode45"))

		return(trajectory)
	}

	## function to compute log-prior
	SIR_exp_logPrior <- function(theta) {

        ## uniform prior on R0: U[1,100]
		log.prior.R0 <- dunif(exp(theta[["R0"]]), min = 1, max = 100, log = TRUE)
        ## uniform prior on infectious period: U[0,30]
		log.prior.D <- dunif(exp(theta[["D.inf"]]), min = 0, max = 30, log = TRUE)

		return(log.prior.R0 + log.prior.D)
	}

	## create deterministic SIR fitmodel
	SIR_exp <- fitmodel(
		name=SIR_name,
		state.names=SIR_state.names,
		theta.names=SIR_theta.names,
		simulate=SIR_exp_simulateDeterministic,
		genObsPoint=SIR_genObsPoint,
		logPrior=SIR_exp_logPrior,
		pointLogLike=SIR_pointLogLike)

	save(SIR_exp,file=file.path(dir_pkg,"data","SIR_exp.rdata"))


}


create_data_fitmodel_SEITL <- function() {


	SEITL_deter_name <- "deterministic SEITL model with daily incidence and constant population size"
	SEITL_state.names <- c("S","E","I","T","L","Inc")
	SEITL_theta.names <- c("R0", "D.lat", "D.inf", "alpha", "D.imm", "rho")


	SEITL_simulateDeterministic <- function(theta,init.state,times) {

	# Solves the system of ordinary differential equations for the SEITL model.
	# note the new state Inc for the daily incidence

		SEITL_ode <- function(time, state, theta) {

			# param
			beta <- theta[["R0"]]/theta[["D.inf"]]
			epsilon <- 1/theta[["D.lat"]]
			nu <- 1/theta[["D.inf"]]
			alpha <- theta[["alpha"]]
			tau <- 1/theta[["D.imm"]]

			# states
			S <- state[["S"]]
			E <- state[["E"]]
			I <- state[["I"]]
			T <- state[["T"]]
			L <- state[["L"]]
			Inc <- state[["Inc"]]

			N <- S + E +I + T + L

			dS <- -beta*S*I/N + (1-alpha)*tau*T
			dE <- beta*S*I/N - epsilon*E
			dI <- epsilon*E - nu*I
			dT <- nu*I - tau*T
			dL <- alpha*tau*T
			dInc <- epsilon*E

			return(list(c(dS,dE,dI,dT,dL,dInc)))
		}


		# put incidence at 0 in init.state
		init.state["Inc"] <- 0

		traj <- as.data.frame(ode(init.state, times, SEITL_ode, theta, method = "ode45"))

		# compute incidence of each time interval
		traj <- mutate(traj,Inc=c(0,diff(Inc)))

		return(traj)

	}


	SEITL_genObsPoint <- function(model.point, theta){
		# Generate an observed incidence under a Poisson observation process.  

		obs.point <- rpois(n=1, lambda=theta[["rho"]]*model.point[["Inc"]])

		return(c(obs=obs.point))
	}

	SEITL_logPrior <- function(theta) {
		# Evaluate the log of the prior density distribution of the parameter values.

		log.prior.R0 <- dunif(theta[["R0"]], min = 1, max = 50, log = TRUE)
		log.prior.latent.period <- dunif(theta[["D.lat"]], min = 0, max = 10, log = TRUE)
		log.prior.infectious.period <- dunif(theta[["D.inf"]], min = 0, max = 15, log = TRUE)
		log.prior.temporary.immune.period <- dunif(theta[["D.imm"]], min = 0, max = 50, log = TRUE)
		log.prior.probability.long.term.immunity <- dunif(theta[["alpha"]], min = 0, max = 1, log = TRUE)
		log.prior.reporting.rate <- dunif(theta[["rho"]], min = 0, max = 1, log = TRUE)

		return(log.prior.R0 + log.prior.latent.period + log.prior.infectious.period + log.prior.temporary.immune.period + log.prior.probability.long.term.immunity + log.prior.reporting.rate)

	}


	SEITL_pointLogLike <- function(data.point, model.point, theta){
		# Computes the log-likelihood of a data point given the state of the model and under a poisson observation process.
		return(dpois(x=data.point[["obs"]],lambda=theta[["rho"]]*model.point[["Inc"]],log=TRUE))

	}


	# create fitmodel
	SEITL_deter <- fitmodel(
		name=SEITL_deter_name,
		state.names=SEITL_state.names,
		theta.names=SEITL_theta.names,
		simulate=SEITL_simulateDeterministic,
		genObsPoint=SEITL_genObsPoint,
		logPrior=SEITL_logPrior,
		pointLogLike=SEITL_pointLogLike)

	save(SEITL_deter,file=file.path(dir_pkg,"data","SEITL_deter.rdata"))


	SEITL_sto_name <- "stochastic SEITL model with daily incidence and constant population size"

	SEITL_simulateStochastic <- function(theta,init.state,times) {

		# Simulate realisation of the stochastic version of the SEITL model.

		SEITL_transitions <- list(
			c(S=-1,E=1),# infection
			c(E=-1,I=1,Inc=1),# infectiousness + incidence
			c(I=-1,T=1),# recovery + short term protection
			c(T=-1,L=1),# efficient long term protection
			c(T=-1,S=1)# deficient long term protection
			)

		SEITL_rateFunc <- function(state,theta,t) {

			# param
			beta <- theta[["R0"]]/theta[["D.inf"]]
			epsilon <- 1/theta[["D.lat"]]
			nu <- 1/theta[["D.inf"]]
			alpha <- theta[["alpha"]]
			tau <- 1/theta[["D.imm"]]

			# states
			S <- state[["S"]]
			E <- state[["E"]]
			I <- state[["I"]]
			T <- state[["T"]]
			L <- state[["L"]]
			Inc <- state[["Inc"]]

			N <- S + E +I + T + L

			return(c(
				beta*S*I/N, # infection
				epsilon*E, # infectiousness + incidence
				nu*I, # recovery + short term protection
				alpha*tau*T, # efficient long term protection
				(1-alpha)*tau*T # deficient long term protection
				)
			)
		}

		# put incidence at 0 in init.state
		init.state["Inc"] <- 0

		traj <- simulateModelStochastic(theta,init.state,times,SEITL_transitions,SEITL_rateFunc) 

		# compute incidence of each time interval
		traj <- mutate(traj,Inc=c(0,diff(Inc)))

		return(traj)

	}


	SEITL_stoch <- fitmodel(
		name=SEITL_sto_name,
		state.names=SEITL_state.names,
		theta.names=SEITL_theta.names,
		simulate=SEITL_simulateStochastic,
		genObsPoint=SEITL_genObsPoint,
		logPrior=SEITL_logPrior,
		pointLogLike=SEITL_pointLogLike)

	save(SEITL_stoch,file=file.path(dir_pkg,"data","SEITL_stoch.rdata"))


	SEIT2L_deter_name <- "deterministic SEIT2L model with daily incidence and constant population size"
	SEIT2L_state.names <- c("S","E","I","T1", "T2","L","Inc")

	SEIT2L_simulateDeterministic <- function(theta,init.state,times) {

		SEIT2L_ode <- function(time, state, theta) {

			# param
			beta <- theta[["R0"]]/theta[["D.inf"]]
			epsilon <- 1/theta[["D.lat"]]
			nu <- 1/theta[["D.inf"]]
			alpha <- theta[["alpha"]]
			tau <- 1/theta[["D.imm"]]

			# states
			S <- state[["S"]]
			E <- state[["E"]]
			I <- state[["I"]]
			T1 <- state[["T1"]]
			T2 <- state[["T2"]]
			L <- state[["L"]]
			Inc <- state[["Inc"]]

			N <- S + E +I + T1 + T2 + L

			dS <- -beta*S*I/N + (1-alpha)*2*tau*T2
			dE <- beta*S*I/N - epsilon*E
			dI <- epsilon*E - nu*I
			dT1 <- nu*I - 2*tau*T1
			dT2 <- 2*tau*T1 - 2*tau*T2
			dL <- alpha*2*tau*T2
			dInc <- epsilon*E

			return(list(c(dS,dE,dI,dT1,dT2,dL,dInc)))
		}


		# put incidence at 0 in init.state
		init.state["Inc"] <- 0

		traj <- as.data.frame(ode(init.state, times, SEIT2L_ode, theta, method = "ode45"))

		# compute incidence of each time interval
		traj <- mutate(traj,Inc=c(0,diff(Inc)))

		return(traj)

	}


	SEIT2L_deter <- fitmodel(
		name=SEIT2L_deter_name,
		state.names=SEIT2L_state.names,
		theta.names=SEITL_theta.names,
		simulate=SEIT2L_simulateDeterministic,
		genObsPoint=SEITL_genObsPoint,
		logPrior=SEITL_logPrior,
		pointLogLike=SEITL_pointLogLike)

	save(SEIT2L_deter,file=file.path(dir_pkg,"data","SEIT2L_deter.rdata"))



	SEIT2L_sto_name <- "stochastic SEIT2L model with daily incidence and constant population size"
	SEIT2L_state.names <- c("S","E","I","T1", "T2","L","Inc")

	SEIT2L_simulateStochastic <- function(theta,init.state,times) {

	# Simulate realisation of the stochastic version of the SEIT2L model.

		SEIT2L_transitions <- list(
			c(S=-1,E=1),# infection
			c(E=-1,I=1,Inc=1),# infectiousness + incidence
			c(I=-1,T1=1),# recovery + temporary protection
			c(T1=-1,T2=1),# progression of temporary protection
			c(T2=-1,L=1),# efficient long term protection
			c(T2=-1,S=1)# deficient long term protection
			)

		SEIT2L_rateFunc <- function(state,theta,t) {

			# param
			beta <- theta[["R0"]]/theta[["D.inf"]]
			epsilon <- 1/theta[["D.lat"]]
			nu <- 1/theta[["D.inf"]]
			alpha <- theta[["alpha"]]
			tau <- 1/theta[["D.imm"]]

			# states
			S <- state[["S"]]
			E <- state[["E"]]
			I <- state[["I"]]
			T1 <- state[["T1"]]
			T2 <- state[["T2"]]
			L <- state[["L"]]
			Inc <- state[["Inc"]]

			N <- S + E +I + T1 + T2 + L

			return(c(
				beta*S*I/N, # infection (S -> E)
				epsilon*E, # infectiousness + incidence (E -> I)
				nu*I, # recovery + short term protection (I -> T1)
				2*tau*T1, # progression of temporary protection (T1 -> T2)
				alpha*2*tau*T2, # efficient long term protection (T2 -> L)
				(1-alpha)*2*tau*T2 # deficient long term protection (T2 -> S)
				)
			)
		}

		# put incidence at 0 in init.state
		init.state["Inc"] <- 0

		traj <- simulateModelStochastic(theta,init.state,times,SEIT2L_transitions,SEIT2L_rateFunc) 

		# compute incidence of each time interval
		traj <- mutate(traj,Inc=c(0,diff(Inc)))

		return(traj)

	}


	SEIT2L_stoch <- fitmodel(
		name=SEIT2L_sto_name,
		state.names=SEIT2L_state.names,
		theta.names=SEITL_theta.names,
		simulate=SEIT2L_simulateStochastic,
		genObsPoint=SEITL_genObsPoint,
		logPrior=SEITL_logPrior,
		pointLogLike=SEITL_pointLogLike)


	save(SEIT2L_stoch,file=file.path(dir_pkg,"data","SEIT2L_stoch.rdata"))


}

simulate_SEITL <- function(SEITL) {
	# create
	SEITL <-  createSEITLmodelTDC(deter=FALSE,FALSE)

	# simulate model
	tmp <- simulateModelReplicates(SEITL,0:60,50)
	plotTraj(SEITL,tmp,state=c("I"),alpha=0.25,plot=TRUE)

}


SEITL_smc <- function() {

	enableJIT(1)
	SEITL <- create_SEITL()
	# test_par <- c(8.5524028575751,1.13799444869168,3.16879027702902,9.78149740477471, 0.873570387495651, 0.663430402689141, 0.069814645950058, 0.0212610261591055, 284)
	# names(test_par) <- names(getParameterValue(SEITL$parameters))
	# SEITL$parameters <- revalueParameters(SEITL$parameters,revalue=test_par)

	system.time(smc <- particleFilter(SEITL,10))


	smc <- particleFilter(SEITL,10)
	print(smc$logLikelihood)
	plotSMC(smc,SEITL,0.1)


	np <- as.list(c(10,25,50,75,100))
	names(np) <- np
	tmp <- ldply(np,function(nParticles) {
		ll <- llply(1:10,function(x) {
			smc <- particleFilter(fitmodel=SEITL,nParticles=nParticles,n.cores=NULL)
			return(smc$logLikePoint)
		},.progress="text")
		ll <- unlist(ll)
		return(data.frame(mean=mean(ll),sd=sd(ll)))
	},.id="nParticles")

	

	# profiling

	pf_ex <- profr(particleFilter(SEITL,100),0.02)
	summary(pf_ex)
	plot(pf_ex)
	head(pf_ex)
	df_time <- ddply(pf_ex,c("f","level"),function(df){return(data.frame(time=sum(df$time)))})
	df_time <- arrange(df_time,level,time)

}

test_update <- function() {

	n <- 100
	x <- rmvnorm(n,mean=c(1,1),sigma=matrix(c(1,0.5,0.5,1),nrow=2,byrow=T))
	x <- as.data.frame(x)
	names(x) <- c("a","b")

	covmat <- cov(x)
	covmat[,] <- 0
	theta.mean <- colMeans(x)
	n <- n + 1

	theta <- rmvnorm(1,mean=c(0,0),sigma=matrix(c(1,0.5,0.5,1),nrow=2,byrow=T))
	theta <- as.vector(theta)
	names(theta) <- c("a","b")
	# test update.covmat
	tmp <- update.covmat(covmat,theta.mean,theta,n)
	print(tmp)

	# test update.cov
	# tmp <- update.cov(covmat,n,theta.mean,theta)
	# print(tmp)

	# true
	x <- rbind(x,theta)
	theta.mean <- colMeans(x)

	print(list(covmat=cov(x),theta.mean=colMeans(x)))

	# naive computation
	for(i in 1:length(theta)){
		for(j in 1:length(theta)){
			covmat[i,j] <- sum((x[,i]-theta.mean[i])*(x[,j]-theta.mean[j]))/n
		}
	}
	print(list(covmat=covmat,theta.mean=theta.mean))

	#
	print(two_pass_cov(x$a,x$b))

	print(var(x$a))
	print(var(x$b))

}

test_smc <- function() {

	data(SEITL_stoch)

	theta <- c("R0"=10, "D.lat"=2 , "D.inf"=3, "alpha"=0.5, "D.imm"=15, "rho"=0.7)
	state.init <- c("S"=280,"E"=0,"I"=2,"T"=0,"L"=4,"Inc"=0)
	data("FluTdC1971",envir = environment())
	data <- FluTdC1971[1:5,]

	x <- margLogLikeSto(fitmodel=SEITL, theta=theta, state.init=state.init, data=data, n.particles=10)
	expect_true(is.numeric(x))


}


test_parallel <- function() {


	library(parallel)
	library(doParallel)

	SEITL <- SEITL_createFitmodel(deterministic=FALSE, verbose=TRUE) 

	system.time(x <- particleFilter(fitmodel= SEITL, n.particles=1000, progress = TRUE, n.cores = 1))
	system.time(xp <- particleFilter(fitmodel= SEITL, n.particles=1000, progress = TRUE, n.cores = NULL))
	system.time(xp2 <- particleFilter(fitmodel= SEITL, n.particles=1000, progress = TRUE, n.cores = 4))

}

test_mcmcMH <- function() {


	SEITL <- SEITL_createFitmodel(deterministic=TRUE, verbose=TRUE) 

	theta.init <- SEITL$theta

	ans <- mcmcMH(target=posteriorDensity, target.args=list(logPrior=SEITL$logPrior, marginal.logLikePoint= margLogLikeDeter, marginal.logLikePoint.args=list(fitmodel=SEITL)), theta.init=theta.init, gaussian.proposal=SEITL$gaussian.proposal, n.iterations=100, adapt.size.start=10, adapt.size.cooling=0.99, adapt.shape.start=10, print.info.every=200)
	
	trace_thined <- burnAndThin(ans$trace)


	plotTrace(trace_thined)
	plotPosteriorDistribution(trace_thined)

	SEITL_stoch <- SEITL_createFitmodel(deterministic=FALSE, verbose=TRUE) 

	fit <- plotPosteriorFit(trace_thined,SEITL_stoch,sample.size=300)


}


test_plot <- function() {

	SEIT2L_deter <- SEIT2L_createFitmodel("deterministic")
	SEIT2L_stoch <- SEIT2L_createFitmodel("stochastic")
	list_model <- list(SEITL_deter,SEITL_stoch)

	theta <- c("R0"=10, "D.lat"=2 , "D.inf"=3, "alpha"=0.5, "D.imm"=15, "rho"=0.7)
	state.init <- c("S"=280,"E"=0,"I"=2,"T"=0,"L"=4,"Inc"=0)
	data("FluTdC1971",envir = environment())
	times <- c(0,FluTdC1971$time)

	traj <- genObsTraj(SEITL_deter, theta, state.init, times)



	traj <- simulateModelReplicates(fitmodel=SEITL_stoch,theta=theta, state.init=state.init, times=times, n=200, observation=TRUE)

	plotTraj(traj,data=FluTdC1971,summary=TRUE,alpha=0.2)

	theta.guess3 <- c("R0"=10, "D.lat"=2 , "D.inf"=2, "alpha"=0.5, "D.imm"=15, "rho"=0.70)
	state.init.guess3 <- c("S"=282,"E"=0,"I"=2,"T"=0,"L"=0,"Inc"=0)

	plotFit(SEITL_stoch,theta.guess3,state.init.guess3,data=FluTdC1971,summary=TRUE, n=100, p.extinct=TRUE)

	theta.guess3 <- c("R0"=10, "D.lat"=2 , "D.inf"=2, "alpha"=0.5, "D.imm"=15, "rho"=0.70)
	state.init.guess3 <- c("S"=282,"E"=0,"I"=2,"T1"=0,"T2"=0,"L"=0,"Inc"=0)

	plotFit(SEIT2L_stoch,theta.guess3,state.init.guess3,data=FluTdC1971,summary=TRUE, n=100, p.extinct=TRUE)

}

test_ABC <- function() {

	SEITL <- SEITL_createFitmodel(deterministic=TRUE, verbose=TRUE) 

	theta.init <- SEITL$theta

	# choose epsilon by looking at the distance for theta.init
	epsilon <- computeDistanceABC(theta.init,SEITL)
	# 0.9
	# look at the fit
	plotThetaFit(theta.init,SEITL,100)
	# looks good
	# choose epsilon = 1 to be less stringent

	# non adaptive
	ans <- mcmcMH(target=targetPosteriorABC, target.args=list(fitmodel=SEITL,epsilon=1), theta.init=theta.init, gaussian.proposal=SEITL$gaussian.proposal, n.iterations=1000)

	# adaptative
	ans <- mcmcMH(target=targetPosteriorABC, target.args=list(fitmodel=SEITL,epsilon=1), theta.init=theta.init, gaussian.proposal=SEITL$gaussian.proposal, n.iterations=1000, adapt.size.start=20, adapt.size.cooling=0.99, adapt.shape.start=50)

	trace_thined <- burnAndThin(ans$trace,0.1,500)

	plotTrace(trace_thined,TRUE)
	plotPosteriorDistribution(trace_thined,TRUE)
	plotPosteriorFit(trace_thined,SEITL)

}

step_by_step <- function() {

	# define parameters using the fitparam class
	R0 <- fitparam(name="R0",value=3)

	InfectiousPeriod <- fitparam(name="D.inf",value=3)

	ReportingRate <- fitparam(name="rho",value=0.7)

	proportionI0 <- fitparam(name="pI0",value=1/300)
	proportionR0 <- fitparam(name="pR0",value=0)

	PopSize <- fitparam(name="N",value=300)


	# create NAMED list
	parameters <- list(R0,InfectiousPeriod,ReportingRate,proportionI0,proportionR0,PopSize)
	
	x <- getParameterValues(parameters)
	y <- setParameterValues(parameters,c("rho"=1))

	# state variables
	state.names <- c("S","I","R")

	# function to initialise the model
	SIR_initialiseState <- function(parameters) {

		# constant pop size
		N <- parameters$N$value

		# number of infected and immune
		I <- round(parameters$pI0$value*N)
		R <- round(parameters$pR0$value*N)

		if(I+R>N){
			stop("Initial conditions not valid")
		}

		return(c(S=N-I-R,I=I,R=R,Inc=0))
	}

}


play_guess <- function() {

	data(SEITL_deter)
	data(SEIT2L_deter)
	data(SEIT2L_stoch)
	data(FluTdC1971)
	# example(SEITL_guess_values)
	theta.guess3 <- c("R0"=10, "D.lat"=2 , "D.inf"=2, "alpha"=0.5, "D.imm"=13, "rho"=0.7)
	init.state.guess3 <- c("S"=279,"E"=0,"I"=2,"T"=3,"L"=0,"Inc"=0)
	init.state.guess5 <- c("S"=277,"E"=0,"I"=1,"T"=0,"L"=6,"Inc"=0)

	theta.guess4 <- c("R0"=11.20, "D.lat"=2.07 , "D.inf"=2.44, "alpha"=0.53, "D.imm"=11.57, "rho"=0.71)

	theta.guess5 <- c("R0"=6.3721, "D.lat"=1.3424 , "D.inf"=2.7904, "alpha"=0.4867, "D.imm"=10.0413, "rho"=0.6745)
	theta.guess6 <- c("R0"=6.4120, "D.lat"=1.3541 , "D.inf"=2.7832, "alpha"=0.4859, "D.imm"=10.0247, "rho"=0.6732)

	plotFit(SEITL_deter,theta.guess5,init.state.guess3,FluTdC1971,all.vars=TRUE)
	trajLogLike(SEITL_deter,theta.guess5,init.state.guess3,FluTdC1971)
	
	quartz()
	init.state.guess4 <- c("S"=279,"E"=0,"I"=2,"T1"=3,"T2"=0,"L"=0,"Inc"=0)
	init.state.guess6 <- c("S"=277,"E"=0,"I"=1,"T1"=0,"T2"=0,"L"=6,"Inc"=0)
	quartz()
	plotFit(SEIT2L_deter,theta.guess5,init.state.guess4,FluTdC1971,all.vars=TRUE)
	trajLogLike(SEIT2L_deter,theta.guess6,init.state.guess4,FluTdC1971)


}

analyse_mcmc <- function() {

	library(fitR)
	data(mcmc_TdC_deter_longRun)
	data(SEITL_deter)
	data(SEIT2L_deter)
	data(SEIT2L_stoch)
	# combine traces
	trace1 <- mcmc(mcmc_SEITL_theta1$trace)
	trace2 <- mcmc(mcmc_SEITL_theta2$trace)
	trace <- mcmc.list(list(trace1,trace2))

	1-rejectionRate(trace)
	#

	xyplot(trace)

	# visual assessment cut at 5000
	trace.burn <- burnAndThin(trace, burn=5000)

	effectiveSize(trace.burn)
	acfplot(trace.burn)

	# big ESS cut at 40
	trace.burn.thin <- burnAndThin(trace.burn, thin=40)
	xyplot(trace.burn.thin)

	# check
	effectiveSize(trace.burn.thin)
	acfplot(trace.burn.thin)

	plotPosteriorDensity()

	# 
	quartz()
	densityplot(trace.burn.thin)


	init.state1 <- c("S"=279,"E"=0,"I"=2,"T"=3,"L"=0,"Inc"=0) 
	# init.state2 <- c("S"=279,"E"=0,"I"=2,"T1"=3,"T2"=0,"L"=0,"Inc"=0) 
	data <- FluTdC1971

	plotPosteriorFit(trace,SEITL_deter,init.state1,data,posterior.summary="max",all.vars=FALSE)
	

	summary(trace1)
	summary(trace2)

	quartz()
	init.state2 <- c("S"=279,"E"=0,"I"=2,"T1"=3,"T2"=0,"L"=0,"Inc"=0) 
	plotPosteriorFit(trace1,SEIT2L_deter,init.state2,data,posterior.summary="mean",all.vars=TRUE)

	summary(trace.burn.thin)

	if(0){
		trace <- trace1
		fitmodel <- SEITL_deter
		posterior.summary <- "max"
	}

}

test_SEIT2L <- function() {

	trace.info1 <- mcmc(mcmc_SEITL_infoPrior_theta1$trace)
	trace.info2 <- mcmc(mcmc_SEITL_infoPrior_theta2$trace)

# combine in a mcmc.list
	trace.info <- mcmc.list(trace.info1,trace.info2)

# burn and thin as the chain with uniform prior (see above sections)
	trace.info.burn.thin <-  burnAndThin(trace.info, burn=5000, thin=40)

	summary(trace.info.burn.thin)
	trace.combined <- ldply(trace.info.burn.thin)

	theta.bar <- colMeans(trace.combined[SEIT2L_deter$theta.names])

	theta.mean.posterior <- c("R0"=7.63, "D.lat"=1.29 , "D.inf"=3.67, "alpha"=0.48, "D.imm"=9.14, "rho"=0.65)
	init.state <- c("S"=279,"E"=0,"I"=2,"T1"=3,"T2"=0,"L"=0,"Inc"=0) 
	plotFit(SEIT2L_stoch,theta.max.loglike,state.init,data=FluTdC1971,n=100)



}

test_pmcmc <- function() {

	theta.mean.posterior <- c("R0"=7.63, "D.lat"=1.29 , "D.inf"=3.67, "alpha"=0.48, "D.imm"=9.14, "rho"=0.65)
	init.state <- c("S"=279,"E"=0,"I"=2,"T1"=3,"T2"=0,"L"=0,"Inc"=0) 
	
	theta.mean.abc <- c(R0 = 6.83480393, D.lat = 2.06938072, D.inf = 0.53114347, alpha = 0.47594650, D.imm = 12.54261436, rho = 0.78203204)
	init.state.abc <- c(S = 250, E = 0, I = 4, T1 = 0, T2 = 0, L = 30, Inc = 0)
	init.state.abc <- c(S = 279, E = 0, I = 2, T1 = 3, T2 = 0, L = 0, Inc = 0)

	theta <- theta.mean.abc
	init.state <- init.state.abc

	my_particleFilter(SEIT2L_stoch,theta,init.state, data=FluTdC1971, n.particles=100)

	plotFit(SEIT2L_stoch,theta,init.state,data=FluTdC1971,n=100)

}

generate_knitr <- function() {

	require(knitr)
	require(plyr)
	
	wd <- getwd()
	setwd(dir_knitr)

	all_files <- grep(".Rmd",list.files(dir_knitr),value=TRUE)

	# files <- c("README","fitcourseR","first_fitmodel","data_likelihood","mcmc","play_with_seitl","smc","smc_solution","smc_example")
	# files <- c("play_with_seitl.Rmd")
	# files <- c("pmcmc","pmcmc_solution","smc_example","smc_example_solution")
	# files <- c("play")
	# files <- c("ABC","ABC_example")
	# files <- paste0(c("README","introduction","mcmc","mcmc_diagnostics","data_likelihood","generate_samples"),".Rmd")

	# files <- paste0(c("introduction","posterior_example","posterior_example_solution"),".Rmd")
	# files <- paste0(c("mcmc","mcmc_example","mcmc_example_solution","generate_samples"),".Rmd")
	# files <- paste0(c("mcmc_diagnostics"),".Rmd")
	# files <- paste0(c("play_with_seitl","mcmc_and_model_comparison"),".Rmd") #"play_with_seitl_example"
	# files <- paste0(c("mcmc_and_model_comparison","example_mcmc_seitl","play_with_seitl_example","pmcmc","pmcmc_solution"),".Rmd")
	# files <- paste0(c("mcmc_and_model_comparison","our_ppc_insert","our_ppc"),".Rmd")
	# files <- paste0(c("play_with_seitl","play_with_seitl_example"),".Rmd")
	# files <- paste0(c("pmcmc","pmcmc_solution","smc_example","smc_example_solution"),".Rmd")
	files <- paste0(c("mcmc_and_model_comparison"),".Rmd")

	all_input <- file.path(dir_knitr,files)
	l_ply(all_input,knit)
	setwd(wd)
}

dev <- function(){

	# create R package
	# create(dir_pkg)
	# start_me()
	document(dir_pkg)
	# load_all(dir_pkg)
	# test(dir_pkg)
	# test(dir_pkg,"classe")
	# test(dir_pkg,"simu")
	# test(dir_pkg,"logLike")
	# test(dir_pkg,"mcmc")
	# check(dir_pkg, check_dir=dir_dev, cleanup =FALSE)		

	# dev_help("fitmodel")
	# dev_help("FluTdC1971")
	# dev_help("particleFilter")
	# dev_help("logLikelihood_stochastic")
	# dev_help("SEITL_deter")
	# dev_help("setParameterValues")
	# dev_help("plotThetaFit")
	# dev_help("plotPosteriorFit")
	# dev_help("SEITL_createFitmodel")

}

main <- function() {

	start_me()
	# dev_mode()
	# install_fitR()

	# create_trace_mcmc_deter()
	# create_SMC_benchmark()
	# create_trace_mcmc_sto()
	# create_data_fitmodel_SEITL()
	# create_data_fitmodel_SIR()
	
	dev()
	# test_bootstrap()
	# generate_knitr()
	# analyse_mcmc_SEITL_deter()
	# create_data()
	# document(dir_pkg)
	# load_all(dir_pkg)

	# test_update()


}

main()
