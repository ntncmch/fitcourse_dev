#!/bin/bash

usage() { echo "Usage: $0 [-m <string>: path to model] [-a <string>: path to analysis] [-r <numeric>: replicate id]" 1>&2; exit 1; }

while getopts ":m:a:r:" o; do
    case "${o}" in
        m)
            m=${OPTARG}
            # echo "model: ${m}"
            ;;
        a)
            a=${OPTARG}
            # echo "analysis: ${a}"
            ;;
        r)
            r=${OPTARG}
            # echo "replicate: ${r}"
            ;;
        *)
    
            usage
            ;;
    esac
done
shift $((OPTIND-1))

if [ -z "${m}" ] || [ -z "${a}" ] || [ -z "${r}" ]; then
    usage
fi

cd ${m}/bin

## LHS SIMPLEX ------------------------------------------------------------------------------------------------
# cat ${m}/lhs/theta_${r}.json | ./simplex --iter 10000 --prior --id ${r} --root ${a} > ${a}/map_${r}.json
# rm blocked jobs with
# condor_rm -constraint 'ClusterId == 842805  && ProcId < 197' 

## LHS KSIMPLEX ------------------------------------------------------------------------------------------------
# cat ${m}/lhs/theta_${r}.json | ./ksimplex --iter 10000 --prior --id ${r} --root ${a} > ${a}/map_${r}.json
# rm blocked jobs with
# condor_rm -constraint 'ClusterId == 843431  && ProcId < 450' 

## LHS + kMCMC
# cat ${m}/theta.json | ./simplex --iter 10000 --prior --trace | ./kmcmc --no_dem_sto --iter 10000 --eps_switch 20 --cooling 0.99 --switch 500 --seed_time | ./kmcmc --iter 100000 --eps_switch 20 --trace --acc --traj --n_traj 1000 --seed_time --id 1 --root ${a} > ${a}/theta_mean_kmcmc.json

## Kalman-MCMC after simplex LHS ------------------------------------------------------------------------------------------------
# mkdir ${a}/burning
# cat ${m}/theta_map_simplex.json | ./kmcmc --iter 10000 --eps_switch 20 --switch 500 --cooling 0.99 --seed_time --root ${a}/burning --trace --acc --id ${r} | ./kmcmc --iter 100000 --eps_switch 20 --trace --acc --traj --n_traj 1000 --seed_time --id ${r} --root ${a} > ${a}/theta_mean_kmcmc_${r}.json

# cat ${m}/theta.json | ./kmcmc --iter 10000 --eps_switch 20 --switch 500 --cooling 0.99 --seed_time --root ${a}/burning --trace --acc --id ${r} | ./kmcmc --iter 100000 --eps_switch 20 --trace --acc --traj --n_traj 1000 --seed_time --id ${r} --root ${a} > ${a}/theta_mean_kmcmc_${r}.json
# cat ${m}/map_*.json | ./kmcmc --iter 10000 --eps_switch 20 --switch 500 --cooling 0.99 --seed_time --root ${a}/burning --trace --acc --id ${r} | ./kmcmc --iter 100000 --eps_switch 20 --trace --acc --traj --n_traj 1000 --seed_time --id ${r} --root ${a} > ${a}/theta_mean_kmcmc_${r}.json
# cat ${m}/theta_map_simplex.json | ./kmcmc --iter 10000 --eps_switch 20 --cooling 0.99 --switch 500 --seed_time | ./kmcmc --iter 10000 --eps_switch 20 --switch 500 --cooling 0.99 --seed_time --root ${a}/burning --trace --acc --id ${r} | ./kmcmc --iter 100000 --eps_switch 20 --trace --acc --traj --n_traj 1000 --seed_time --id ${r} --root ${a} > ${a}/theta_mean_kmcmc_${r}.json

# cat ${m}/theta_map_simplex.json | ./kmcmc --iter 5000 --eps_switch 20 --cooling 0.99 --switch 500 --seed_time | ./kmcmc --iter 10000 --eps_switch 20 --cooling 0.99 --switch 500 --seed_time --root ${a}/burning --trace --acc --id ${r} | ./kmcmc --iter 20000 --eps_switch 20 --trace --acc --traj --n_traj 500 --seed_time --id ${r} --root ${a} > ${a}/theta_mean_kmcmc_${r}.json
# cat ${m}/theta_map_simplex.json | ./kmcmc --iter 10000 --eps_switch 20 --cooling 0.99 --switch 500 --seed_time --root ${a}/burning --trace --acc --id ${r} | ./kmcmc --iter 20000 --eps_switch 20 --trace --acc --traj --n_traj 500 --seed_time --id ${r} --root ${a} > ${a}/theta_mean_kmcmc_${r}.json
# cat ${m}/theta_map_simplex.json | ./kmcmc --iter 5000 --eps_switch 20 --cooling 0.99 --switch 500 --seed_time --root ${a}/burning --trace --acc --id ${r} | ./kmcmc --iter 5000 --eps_switch 20 --cooling 0.99 --switch 500 --seed_time --root ${a}/burning --trace --acc --id ${r} > ${a}/theta_mean_kmcmc_short_${r}.json
# cat ${m}/theta_mean_kmcmc_short.json | ./kmcmc --iter 100000 --eps_switch 20 --trace --acc --traj --n_traj 2000 --seed_time --id ${r} --root ${a} > ${a}/theta_mean_kmcmc_${r}.json
# cat ${m}/theta_combined.json | ./kmcmc --iter 10000 --eps_switch 20 --cooling 0.99 --seed_time --root ${a}/burning --trace --acc --id ${r} | ./kmcmc --iter 100000 --eps_switch 20 --trace --acc --traj --n_traj 3000 --seed_time --id ${r} --root ${a} > ${a}/theta_mean_kmcmc_${r}.json
# rm blocked jobs with
# condor_rm -constraint 'ClusterId == 842777  && ProcId < 780'

## SMC ------------------------------------------------------------------------------------------------
# cat ../theta_map_simplex.json | ./smc -J 500 -N 5 -v --seed_time
# cat ../theta_map_simplex.json | ./pmcmc psr -D 0.1 -J 500 --iter 10000 --eps_switch 20 --cooling 0.99 --switch 500 --trace --acc --traj --n_traj 1000 --id 3 -v > ../theta_mean_mcmc.json

## MCMC ------------------------------------------------------------------------------------------------
# mkdir ${a}/burning
# cat ${m}/theta_map_simplex.json | ./pmcmc --iter 10000 --eps_switch 20 --cooling 0.99 --switch 500 --seed_time | ./pmcmc --iter 10000 --eps_switch 20 --switch 500 --cooling 0.99 --seed_time --root ${a}/burning --trace --acc --id ${r} | ./pmcmc --iter 100000 --eps_switch 20 --trace --acc --traj --n_traj 2000 --seed_time --id ${r} --root ${a} > ${a}/theta_mean_mcmc_${r}.json
# cat ${m}/theta_map_simplex.json | ./pmcmc --iter 10000 --eps_switch 20 --cooling 0.99 --switch 500 --seed_time | ./pmcmc --iter 10000 --eps_switch 20 --cooling 0.99 --switch 500 --seed_time | ./pmcmc --iter 200000 --eps_switch 20 --trace --acc --traj --n_traj 2000 --seed_time --id ${r} --root ${a} > ${a}/theta_mean_mcmc_${r}.json
# cat ${m}/theta_combined.json | ./pmcmc --iter 10000 --eps_switch 20 --cooling 0.99 --switch 500 --seed_time | ./pmcmc --iter 10000 --eps_switch 20 --cooling 0.99 --switch 500 --seed_time | ./pmcmc --iter 200000 --eps_switch 20 --trace --acc --traj --n_traj 1000 --seed_time --id ${r} --root ${a} > ${a}/theta_mean_mcmc_${r}.json
# cat ${m}/theta_map_simplex.json | ./pmcmc --no_dem_sto --no_diff --iter 10000 --eps_switch 20 --cooling 0.99 --switch 500 --seed_time | ./pmcmc --no_dem_sto --no_diff --iter 10000 --eps_switch 20 --cooling 0.99 --switch 500 --seed_time | ./pmcmc --no_dem_sto --no_diff --iter 100000 --eps_switch 20 --cooling 0.99 --trace --acc --traj --n_traj 1000 --id ${r} --root ${a} > ${a}/theta_mean_mcmc_${r}.json
# rm blocked jobs with
# condor_rm -constraint 'ClusterId == 842777  && ProcId < 780'

## pMCMC ------------------------------------------------------------------------------------------------
mkdir ${a}/burning
cat ${m}/theta.json | ./pmcmc --iter 10000 --eps_switch 100 --cooling 0.99 --switch 500 --seed_time | ./pmcmc --iter 10000 --eps_switch 100 --cooling 0.99 --switch 500 --seed_time | ./pmcmc psr -J 20 -D 0.05 --iter 3000 --eps_switch 50 --cooling 0.99 --switch 100 --seed_time --root ${a}/burning --trace --traj --n_traj 500 --acc --id ${r} 

#| ./pmcmc psr -J 396 -N 12 -D 0.05 --iter 3000 --eps_switch 50 --cooling 0.99 --switch 100 --trace --acc --traj --n_traj 500 --seed_time --id ${r} --root ${a} > ${a}/theta_out_${r}.json

# cat ../theta_mean_pmcmc.json | ./pmcmc --iter 5000 --eps_switch 20 --cooling 0.99 --switch 500 --seed_time --trace --traj --n_traj 500 --acc --n_obs 10

# cat ${m}/theta_map_simplex.json | ./pmcmc psr -J 2000 -N 12 -D 1 --iter 10000 --eps_switch 20 --cooling 0.99 --switch 500 --seed_time --root ${a}/burning --trace --traj --n_traj 500 --acc --id ${r} | ./pmcmc psr -J 2000 -N 12 -D 1 --iter 100000 --eps_switch 20 --trace --acc --traj --n_traj 3000 --seed_time --id ${r} --root ${a} > ${a}/theta_mean_pmcmc_${r}.json
# cat ${m}/theta_mean_kmcmc.json | ./pmcmc psr -J 2000 -N 12 -D 1 --iter 10000 --eps_switch 20 --cooling 0.99 --switch 500 --seed_time --root ${a}/burning --trace --traj --n_traj 500 --acc --id ${r} | ./pmcmc psr -J 2000 -N 12 -D 1 --iter 100000 --eps_switch 20 --trace --acc --traj --n_traj 3000 --seed_time --id ${r} --root ${a} > ${a}/theta_mean_pmcmc_${r}.json

# cat ${m}/theta_map_simplex.json | ./pmcmc psr -J 2000 -N 12 -D 1 --iter 10000 --eps_switch 20 --cooling 0.99 --switch 500 --seed_time --root ${a}/burning --trace --traj --n_traj 500 --acc --id ${r} | ./pmcmc psr -J 2000 -N 12 -D 1 --iter 100000 --eps_switch 20 --trace --acc --traj --n_traj 3000 --seed_time --id ${r} --root ${a} > ${a}/theta_mean_pmcmc_${r}.json
# cat ${m}/theta_map_simplex.json | ./pmcmc psr -J 2000 -N 12 -D 1 --iter 10000 --eps_switch 20 --cooling 0.99 --switch 500 --trace --acc --traj --n_traj 1000 --seed_time --id ${r} --root ${a} > ${a}/theta_mean_pmcmc_${r}.json
# cat ${m}/map_*.json | ./pmcmc psr -J 1000 -N 10 -D 1 --iter 1000 --eps_switch 20 --cooling 0.99 --switch 500 --seed_time --root ${a}/burning --trace --acc --id ${r} | ./pmcmc psr -J 1000 -N 10 -D 1 --iter 10000 --eps_switch 20 --trace --acc --traj --n_traj 2000 --seed_time --id ${r} --root ${a} > ${a}/theta_mean_pmcmc_${r}.json
# cat ${m}/theta_mean_kmcmc.json | ./pmcmc psr -J 2000 -N 10 -D 1 --iter 10000 --eps_switch 20 --cooling 0.99 --switch 500 --seed_time --root ${a}/burning --trace --traj --n_traj 500 --acc --id ${r} | ./pmcmc psr -J 2000 -N 10 -D 1 --iter 100000 --eps_switch 20 --trace --acc --traj --n_traj 2000 --seed_time --id ${r} --root ${a} > ${a}/theta_mean_pmcmc_${r}.json
# cat ${m}/theta_mean_pmcmc.json | ./pmcmc psr -J 2000 -N 10 -D 1 --iter 10000 --eps_switch 20 --cooling 0.99 --switch 500 --seed_time --root ${a}/burning --trace --traj --n_traj 500 --acc --id ${r} | ./pmcmc psr -J 2000 -N 10 -D 1 --iter 100000 --eps_switch 20 --trace --acc --traj --n_traj 2000 --seed_time --id ${r} --root ${a} > ${a}/theta_mean_pmcmc_${r}.json
# rm blocked jobs with
# condor_rm -constraint 'ClusterId == 842777  && ProcId < 780'


## assess forecast ------------------------------------------------------------------------------------------------
# mkdir ${a}/burning
# cat ${m}/theta_mean_pmcmc.json | ./pmcmc psr -J 2000 -N 12 -D 1 --iter 5000 --eps_switch 20 --cooling 0.99 --switch 500 --seed_time --root ${a}/burning --trace --traj --n_traj 500 --acc --id ${r} --n_obs $((${r}+3)) | ./pmcmc psr -J 2000 -N 12 -D 1 --iter 50000 --eps_switch 20 --trace --acc --traj --n_traj 5000 --seed_time --id ${r} --n_obs $((${r}+3)) --root ${a} > ${a}/theta_mean_pmcmc_${r}.json
# cat ${m}/theta_mean_pmcmc.json | ./pmcmc --iter 5000 --eps_switch 20 --cooling 0.99 --switch 500 --seed_time --root ${a}/burning --trace --traj --n_traj 500 --acc --id ${r} --n_obs $((${r}+3)) | ./pmcmc --iter 50000 --eps_switch 20 --trace --acc --traj --n_traj 5000 --seed_time --id ${r} --n_obs $((${r}+3)) --root ${a} > ${a}/theta_mean_pmcmc_${r}.json


## FORECAST ------------------------------------------------------------------------------------------------
# cd ${m}/forecast/models/strategy_${r}/bin
# ssm-predict ${m}/mcmc/mcmc_0.json ${m}/mcmc/X_0.csv ${m}/mcmc/trace_0.csv 2014-09-14 | ./simul --traj --eps_abs_integ 1e-6 --eps_rel_integ 1e-6 --interpolator linear --start 2014-09-14 --end 2015-12-28 --freq 7 --verbose --id ${r} --n_thread 12 --root ${a}







