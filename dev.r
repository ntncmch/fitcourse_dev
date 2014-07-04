
start_me <- function() {

	dir_dev <<- "/Users/Tonton/edu/Fit_course/dev"
	dir_pkg <<- "/Users/Tonton/edu/Fit_course/fitR"
	dir_knitr <<- "/Users/Tonton/edu/Fit_course/mfiidd/knitr"
	dir_md <<- "/Users/Tonton/edu/Fit_course/mfiidd"

	# setwd(dir_dev)

	library(jsonlite)
	library(testthat)
	library(plyr)
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

create_trace_mcmc_deter <- function() {

	df_set <- 	expand.grid(theta=1:2,priorInfo=c(FALSE,TRUE),SEIT2L=c(FALSE,TRUE),n_iteration=c(5000,100000))
	df_set$set <- 1:nrow(df_set)
	dir_name <- "/Users/Tonton/edu/Fit_course/dev/dataset/mcmc_deter"
	dir_rds <- file.path(dir_name,"rds")
	dir_fig <- file.path(dir_name,"figures")

	adapt_size_start <- 100 
	adapt_size_cooling <- 0.999
	adapt_shape_start <- 200

	# long trace
	df_prel <- subset(df_set,n_iteration== 5000 & !priorInfo & theta==1)

	list_trace <- dlply(df_long,c("theta","priorInfo","SEIT2L","n_iteration"),function(df) {

		analysis <- paste0("mcmc_",ifelse(df$SEIT2L,"SEIT2L","SEITL"),"_deter_",ifelse(df$priorInfo,"info","unif"),"Prior_n=",df$n_iteration,"_size=",adapt_size_start,"_cool=",adapt_size_cooling,"_shape=",adapt_shape_start,"_set=",df$set,".rds")

		ans <- readRDS(file.path(dir_rds,analysis))

		return(ans)
	})
	
	name <- unique(paste0("mcmc_",ifelse(df_long$SEIT2L,"SEIT2L","SEITL")))
	
	names(list_trace) <- name

	attach(list_trace)

	save(list=name,file=file.path(dir_pkg,"data","mcmc_TdC_deter_shortRun.rdata"))

	
#	data("mcmc_TdC_deter_longRun",envir = environment())
# 	ls()
	names(mcmc_SEITL_theta1)

	# long trace
	df_long <- subset(df_set,n_iteration>5000)

	list_trace <- dlply(df_long,c("theta","priorInfo","SEIT2L","n_iteration"),function(df) {

		analysis <- paste0("mcmc_",ifelse(df$SEIT2L,"SEIT2L","SEITL"),"_deter_",ifelse(df$priorInfo,"info","unif"),"Prior_n=",df$n_iteration,"_size=",adapt_size_start,"_cool=",adapt_size_cooling,"_shape=",adapt_shape_start,"_set=",df$set,".rds")

		ans <- readRDS(file.path(dir_rds,analysis))

		return(ans)
	})
	
	name <- unique(paste0("mcmc_",ifelse(df_long$SEIT2L,"SEIT2L","SEITL"),ifelse(df_long$priorInfo,"_infoPrior_","_"),"theta",df_long$theta))
	
	names(list_trace) <- name

	attach(list_trace)

	save(list=name,file=file.path(dir_pkg,"data","mcmc_TdC_deter_longRun.rdata"))


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

	example(SEITL_sto)

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

	SEITL_sto <- SEITL_createFitmodel(deterministic=FALSE, verbose=TRUE) 

	fit <- plotPosteriorFit(trace_thined,SEITL_sto,sample.size=300)


}


test_plot <- function() {

	SEIT2L_deter <- SEIT2L_createFitmodel("deterministic")
	SEIT2L_sto <- SEIT2L_createFitmodel("stochastic")
	list_model <- list(SEITL_deter,SEITL_sto)

	theta <- c("R0"=10, "D.lat"=2 , "D.inf"=3, "alpha"=0.5, "D.imm"=15, "rho"=0.7)
	state.init <- c("S"=280,"E"=0,"I"=2,"T"=0,"L"=4,"Inc"=0)
	data("FluTdC1971",envir = environment())
	times <- c(0,FluTdC1971$time)

	traj <- genObsTraj(SEITL_deter, theta, state.init, times)



	traj <- simulateModelReplicates(fitmodel=SEITL_sto,theta=theta, state.init=state.init, times=times, n=200, observation=TRUE)

	plotTraj(traj,data=FluTdC1971,summary=TRUE,alpha=0.2)

	theta.guess3 <- c("R0"=10, "D.lat"=2 , "D.inf"=2, "alpha"=0.5, "D.imm"=15, "rho"=0.70)
	state.init.guess3 <- c("S"=282,"E"=0,"I"=2,"T"=0,"L"=0,"Inc"=0)

	plotFit(SEITL_sto,theta.guess3,state.init.guess3,data=FluTdC1971,summary=TRUE, n=100, p.extinct=TRUE)

	theta.guess3 <- c("R0"=10, "D.lat"=2 , "D.inf"=2, "alpha"=0.5, "D.imm"=15, "rho"=0.70)
	state.init.guess3 <- c("S"=282,"E"=0,"I"=2,"T1"=0,"T2"=0,"L"=0,"Inc"=0)

	plotFit(SEIT2L_sto,theta.guess3,state.init.guess3,data=FluTdC1971,summary=TRUE, n=100, p.extinct=TRUE)

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

analyse_mcmc_SEITL_deter <- function() {

	dir_MCMC <- "/Users/Tonton/edu/Fit_course/MCMC"
	# analyse <- "SEITL_deter_n=1e+06_size=50_cool=0.99_shape=100"
	analyse <- "SEITL_deter_infoPrior_n=1e+06_size=50_cool=0.99_shape=100"
	n_run <- 4

	# for(i in 1:n_run){
	# 	run <- paste0(analyse,"_run=",i)
	# 	file <- file.path(dir_MCMC,run,"tracer.txt")
	# 	trace <- readRDS(file.path(dir_MCMC,run,"rds","trace_b=0.1_thin=10000.rds"))
	# 	export2Tracer(trace,file)		
	# }

	# assess fit
	i <- 1
	run <- paste0(analyse,"_run=",i)
	trace <- readRDS(file.path(dir_MCMC,run,"rds","trace_b=0.1_thin=10000.rds"))
	fitmodel <- readRDS(file.path(dir_MCMC,run,"rds","SEITL.rds"))

	SEIT2L <- SEIT2L_createFitmodel(deterministic=TRUE, verbose=TRUE) 

	plotPosteriorFit(trace, fitmodel=SEIT2L,posterior.median=TRUE ,summary=TRUE, sample.size = 100, plot=TRUE)
	quartz()
	plotPosteriorFit(trace, fitmodel,posterior.median=FALSE ,summary=TRUE, sample.size = 100, plot=TRUE)

	# take one sample, plot autocorr, thin and burn, 

}

generate_knitr <- function() {

	require(knitr)
	require(plyr)
	
	wd <- getwd()
	setwd(dir_knitr)
	# files <- c("README","fitcourseR","first_fitmodel","data_likelihood","mcmc","play_with_seitl","smc","smc_solution","smc_example")
	files <- c("play_with_seitl")
	# files <- c("smc","smc_solution","smc_example")
	# files <- c("play")
	# files <- c("ABC","ABC_example")
	# files <- c("README")
	all_input <- file.path(dir_knitr,paste0(files,".Rmd"))
	
	l_ply(all_input,knit)
	
	setwd(wd)
}

dev <- function(){

	# create R package
	# create(dir_pkg)
	# start_me()
	document(dir_pkg,clean=FALSE)
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
	# dev_help("mcmcMH")
	# dev_help("setParameterValues")
	# dev_help("plotThetaFit")
	# dev_help("plotPosteriorFit")
	# dev_help("SEITL_createFitmodel")

}

main <- function() {

	# start_me()
	# dev_mode()
	# install_fitR()

	dev()
	# test_bootstrap()
	# generate_knitr()
	# analyse_mcmc_SEITL_deter()
	# create_data()
	# document(dir_pkg)
	# load_all(dir_pkg)

	# test_update()


# TODO:
	# ABC TdC add second distance on second wave (heigh of the peak? oscillation around second wave?)
	# mutate model
	# use soda and mcmc class for diagnostic
	# add plots for prior/likelihood/posterior
	# model selection


}

main()
