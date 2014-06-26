library(fitcourseR)

my_SEITL <- SEITL_createModelTdC(deterministic=FALSE)

distanceOscillation

# my summary functions

time.data <- my_SEITL$data$time
inc.data <- my_SEITL$data$Inc

# code
my_summaryStatPeak <- function(time,inc){

	# highest peak
	return(max(inc))
}


my_summaryStatPeakTime <- function(time,inc){
	#
	return(time[max(inc)])
}

my_summaryStatFinalSize <- function(time,inc){

	# final size
	return(sum(inc))
}

my_summaryStatPeak0to20 <- function(time,inc){

	# max inc between day 0 and 20
	return(max(inc[time<20]))
}

my_summaryStatPeakAfter20 <- function(time,inc){

	# max inc between day 0 and 20
	return(max(inc[time>20]))
}

my_summaryStatInc0to20 <- function(time,inc){

	# inc between day 0 and 20
	return(inc[time >= 0 & time <= 20])
}


my_summaryStatInc25to45 <- function(time,inc){

	# inc between day 25 and 45
	return(inc[time >= 25 & time <= 45])
}

# test
ss1.data <- my_summaryStatPeak(time=time.data,inc=inc.data)
ss2.data <- my_summaryStatPeakTime(time=time.data,inc=inc.data)
ss3.data <- my_summaryStatFinalSize(time=time.data,inc=inc.data)
ss4.data <- my_summaryStatPeak0to20(time=time.data,inc=inc.data)
ss5.data <- my_summaryStatPeakAfter20(time=time.data,inc=inc.data)
ss6.data <- my_summaryStatInc0to20(time=time.data,inc=inc.data)
ss7.data <- my_summaryStatInc25to45(time=time.data,inc=inc.data)

ss.data <- list(ss1.data,ss2.data,ss3.data,ss4.data,ss5.data,ss6.data,ss7.data)
print(ss.data)



simu <- my_SEITL$simulate.model(theta=my_SEITL$theta, state.init=my_SEITL$initialise.state(my_SEITL$theta),times=c(0,my_SEITL$data$time))
simu.obs <- my_SEITL$generate.observation(model.traj=simu, theta=my_SEITL$theta)
time.simu <- simu.obs$time[-1] # remove initial time
inc.simu <- simu.obs$observation[-1] # remove initial state

# test
ss1.simu <- my_summaryStatPeak(time=time.simu,inc=inc.simu)
ss2.simu <- my_summaryStatPeakTime(time=time.simu,inc=inc.simu)
ss3.simu <- my_summaryStatFinalSize(time=time.simu,inc=inc.simu)
ss4.simu <- my_summaryStatPeak0to20(time=time.simu,inc=inc.simu)
ss5.simu <- my_summaryStatPeakAfter20(time=time.simu,inc=inc.simu)
ss6.simu <- my_summaryStatInc0to20(time=time.simu,inc=inc.simu)
ss7.simu <- my_summaryStatInc25to45(time=time.simu,inc=inc.simu)

ss.simu <- list(ss1.simu,ss2.simu,ss3.simu,ss4.simu,ss5.simu,ss6.simu,ss7.simu)
print(ss.simu)

# distances

my_distance1 <- function(data,simu) {
	# highest peak
	return((data-simu)^2/data)
}


my_distance2 <- function(data,simu) {
	# peak time
	return(abs(data-simu))
}


my_distance3 <- function(data,simu) {
	# final size
	return(abs(data-simu))
}

my_distance4 <- function(data,simu) {
	# peak first wave
	return((data-simu)^2/data)
}

my_distance5 <- function(data,simu) {
	# peak second wave
	return((data-simu)^2/data)
}

my_distance6 <- function(data,simu) {
	# oscillation distance first wave
	return(distanceOscillation(data,simu))
}

my_distance7 <- function(data,simu) {
	# oscillation distance second wave
	return(distanceOscillation(data,simu))
}

# test
d1 <- my_distance1(data=ss1.data,simu=ss1.simu)
d2 <- my_distance2(data=ss2.data,simu=ss2.simu)
d3 <- my_distance3(data=ss3.data,simu=ss3.simu)
d4 <- my_distance4(data=ss4.data,simu=ss4.simu)
d5 <- my_distance5(data=ss5.data,simu=ss5.simu)
d6 <- my_distance6(data=ss6.data,simu=ss6.simu)
d7 <- my_distance7(data=ss7.data,simu=ss7.simu)

print(c(d1,d2,d3,d4,d5,d6,d7))


# acceptance ABC

ABCacceptance <- function(theta,fitmodel,tol) {

    # simulate model
	simu <- fitmodel$simulate.model(theta=theta, state.init=fitmodel$initialise.state(theta),times=c(0,fitmodel$data$time))

    # generate observation
	simu.obs <- fitmodel$generate.observation(model.traj=simu, theta=theta)
	time.simu <- simu.obs$time[-1] # remove initial time
	inc.simu <- simu.obs$observation[-1] # remove initial state

    # compute data and simulated summary statistics
	# data
	time.data <- fitmodel$data$time
	inc.data <- fitmodel$data$Inc

	ss6.data <- my_summaryStatInc0to20(time=time.data,inc=inc.data)
	ss7.data <- my_summaryStatInc25to45(time=time.data,inc=inc.data)
	# simu
	ss6.simu <- my_summaryStatInc0to20(time=time.simu,inc=inc.simu)
	ss7.simu <- my_summaryStatInc25to45(time=time.simu,inc=inc.simu)


    # compute distances between summary statistics
	d6 <- my_distance6(data=ss6.data,simu=ss6.simu)
	d7 <- my_distance7(data=ss7.data,simu=ss7.simu)

    # check that distances are within the tolerances
    # return 0/1
	if(d6<tol[6] && d7 < tol[7]){
		return(1)
	} else {
		return(0)
	}

}



for(i in 1:10){
	print(ABCacceptance(my_SEITL$theta,my_SEITL,tol=c(0,0,0,0,0,10,10)))
}










