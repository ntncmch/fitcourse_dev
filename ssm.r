# create SEIT4L model with SSMinR

# install_github("ntncmch/SSMinR")
library(fitR)
library(SSMinR)
library(dplyr)


SEIT4L_inputs <- list(
	input(name="N",description="population size", value=284, tag="pop_size"),
	input(name="S", description="initial number of susceptible", tag="remainder"),
	input(name="E", description="initial number of exposed", value=0),
	input(name="I", description="initial number of infectious", value=2), 
	input(name="T1", description="initial number of temporary immune", value=3),
	input(name="T2", description="initial number of temporary immune", value=0),
	input(name="T3", description="initial number of temporary immune", value=0),
	input(name="T4", description="initial number of temporary immune", value=0),
	input(name="L", description="initial number of recovered", value=0),
	input(name="R0", description="basic reproduction number", prior=unif(1,50), value=2), 
	input(name="D_lat", description="latent period", prior=truncnorm(mean=2,sd=1,a=0), value=2),
	input(name="D_inf", description="infectious period", prior=truncnorm(mean=2,sd=1,a=0), value=2),
	# input(name="D_lat", description="latent period", prior=unif(0,10), value=2),
	# input(name="D_inf", description="infectious period", prior=unif(0,15), value=2),
	input(name="alpha", description="infectious period", prior=unif(0,1), value=0.8),
	input(name="D_imm", description="infectious period", prior=unif(0,50), value=16),
	input(name="rho", description="reporting rate", prior=unif(0,1), value=0.85),
	input(name="beta", description="effective contact rate",transformation="R0/D_inf")
	)

SEIT4L_reactions <- list(
	reaction(from="S", to="E", description="infection", rate="beta*I/N"),
	reaction(from="E", to="I", description="onset of infectiosity", rate="1/D_lat", accumulators="Inc"),
	reaction(from="I", to="T1", description="recovery + temporary immunity", rate="1/D_inf"),
	reaction(from="T1", to="T2", description="progression of temporary immunity", rate="4/D_imm"),
	reaction(from="T2", to="T3", description="progression of temporary immunity", rate="4/D_imm"),
	reaction(from="T3", to="T4", description="progression of temporary immunity", rate="4/D_imm"),
	reaction(from="T4", to="L", description="long-term immunity", rate="4*alpha/D_imm"),
	reaction(from="T4", to="S", description="no long-term immunity", rate="4*(1-alpha)/D_imm")
	)

SEIT4L_observations <- list(
	poisson_obs(state="Inc", reporting="rho")
	)

data(FluTdC1971)

dataset <- FluTdC1971 %>% select(-time) %>% rename(Inc_obs = obs)

# the model will be created in the default temporary directory. Change the path to "wherever/you/want".
# dir_model <- path.expand("~/edu/Fit_course/dev/ssm")
dir_model <- path.expand("~/fitcourse/dev/ssm")

my_ssm <- new_ssm(
	model_path=file.path(dir_model,"SEIT4L"),
	pop_name="TdC",
	data=dataset,
	start_date=min(dataset$date), # start model integration 1 day before the first observation
	inputs=SEIT4L_inputs,
	reactions=SEIT4L_reactions,
	observations=SEIT4L_observations
	)


# my_ssm_lhs <- my_ssm %>% do_lhs(n=20, do="simplex", trace=FALSE, iter=100, prior=TRUE) %>% get_max_lhs %>% print

# my_ssm_fit_ode <- my_ssm %>% SSMinR::pmcmc(iter=10000, cooling=0.99, eps_switch=100, switch=500) %>% 
# SSMinR::pmcmc(iter=10000, traj=TRUE, trace=TRUE, n_traj=500, cooling=0.99, eps_switch=100, switch=500) %>% print %>% to_tracer %>% plot_X

# my_ssm_fit_ode %>% plot_theta

# my_ssm_fit_ode$covmat



