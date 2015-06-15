#!/usr/bin/perl -w

#use strict;

#######################################################
# run SSM on the cluster                              #
#######################################################
my $HOME=$ENV{HOME};

my $Condor_job_identifier="SSM";
my $model="SEIT4L";
my $analysis="pmcmc_48";

my $Rscript="$HOME/fitcourse/dev/ssm.r";

my $exe= "$HOME/fitcourse/dev/run_ssm.sh";
my $uni= "vanilla";
my $mem= "1*1024";
my $replicate= "1";

# for LHS generate as many points as replicate
$ENV{SSM_N_LHS} = "$replicate";

# for FORECAST generate as many scenarios as replicate
$ENV{SSM_N_FORECAST} = "$replicate";

my @node_list= ("08","09","10","11","12");

#request reserved machines
my $my_requirements="";
foreach $node (@node_list)
{
	$my_requirements .= "(machine == \"ecoevo$node.ecoevo.biologie.ens.fr\")+";
}
#remove the + at the end of the string
$my_requirements=substr($my_requirements,0,-1);


my $job_dir="$HOME/fitcourse/dev/ssm/$model/$analysis";

print "$analysis $model\n";

# name of analysis
$ENV{SSM_ANALYSIS} = "$analysis";

# for build use model name
$ENV{SSM_MODEL} = "$model";

# run R script first
system("Rscript $Rscript");

# rm previous analysis, create the folder
system("rm -rf $job_dir");
system("mkdir -p $job_dir");

my $in= "$Rscript";
my $out= "$job_dir/out_R\$(Process).txt";

my $Condor_jobname="SSM_$analysis";
print "generate $Condor_jobname\n";
open(Condor_SCRIPT,">$Condor_jobname.txt");
print Condor_SCRIPT "executable = $exe\n";
print Condor_SCRIPT "universe = $uni\n";
print Condor_SCRIPT "request_memory = $mem\n";
# print Condor_SCRIPT "+RequiresWholeMachine = False\n";
print Condor_SCRIPT "+RequiresWholeMachine = True\n";
print Condor_SCRIPT "Requirements = $my_requirements\n";
#print Condor_SCRIPT "Requirements = machine == \"ecoevo13.ecoevo.biologie.ens.fr\"\n";
print Condor_SCRIPT "getenv = True\n";

print Condor_SCRIPT "arguments = -m $HOME/fitcourse/dev/ssm/$model -a $job_dir -r \$(Process)\n";
#print Condor_SCRIPT "transfer_files = ALWAYS\n";
# print Condor_SCRIPT "transfer_input_files = $in\n";
print Condor_SCRIPT "notification = NEVER\n";
# print Condor_SCRIPT "input = $in\n";
print Condor_SCRIPT "output = $out\n";
print Condor_SCRIPT "error = $job_dir/error_R\$(Process).txt\n";
print Condor_SCRIPT "log = $job_dir/log_R\$(Process).txt\n";
# print Condor_SCRIPT "environment = ARG1=\$(Process);analysis=$analysis\n";
# print Condor_SCRIPT "environment = analysis=$analysis\n";

print Condor_SCRIPT "queue $replicate\n";

print "submit $Condor_jobname with $replicate replicate\n";

system("condor_submit $Condor_jobname.txt")==0
or die "failed submitting Condor job \n";

#### GIVE THE JOB A SECOND TO GET SUBMITTED ####
sleep(1);
print "submitted.\n";



