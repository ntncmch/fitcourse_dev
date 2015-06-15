#!/usr/bin/perl -w
#######################################################
# run R script on the cluster                         #
#######################################################
my $HOME="/users/ecologie/camacho";

my $Condor_job_identifier="pdRscript";
# my $jobDir="$HOME/fitcourse/test_smc_SEIT2L";
# my $jobDir="$HOME/fitcourse/mcmc_sto_SEIT2L_bad";
# my $jobDir="$HOME/fitcourse/mcmc_deter_SEIT4L";
my $jobDir="$HOME/fitcourse/mcmc_sto_SEIT4L_48";
my $scriptDir="$HOME/fitcourse/dev";

my @Rscript_list= ("run_on_cluster");

# my @node_list= ("02","03","05","06","08","09","10","11","12","13","14","15","16");
my @node_list= ("07","08","09","11","12");

my $exe= "$HOME/bin/R";
my $uni= "vanilla";
my $mem= "1*1024";
my $replicate= "5";

#request reserved machines
my $my_requirements="";
foreach $node (@node_list)
{
    $my_requirements .= "(machine == \"ecoevo$node.ecoevo.biologie.ens.fr\")+";
}
#remove the + at the end of the string
$my_requirements=substr($my_requirements,0,-1);

# remove the .txt (submission, out, log and error files)
system("rm -r $jobDir/*.txt");


foreach $Rscript (@Rscript_list)
{
    #my $scriptDir= "$jobDir";
    #system("rm -r $jobDir/Rsave $jobDir/pdf");
    system("mkdir -p $jobDir");
    my $in= "$scriptDir/$Rscript.r";
    my $out= "$jobDir/out_R\$(Process).txt";

    
    my $Condor_jobname="$Condor_job_identifier.$Rscript";
    print "generate $Condor_jobname\n";
    open(Condor_SCRIPT,">$Condor_jobname.txt");
    print Condor_SCRIPT "executable = $exe\n";
    print Condor_SCRIPT "universe = $uni\n";
    #print Condor_SCRIPT "request_memory = $mem\n";
    print Condor_SCRIPT "+RequiresWholeMachine = True\n";
    print Condor_SCRIPT "Requirements = $my_requirements\n";
    #print Condor_SCRIPT "Requirements = machine == \"ecoevo13.ecoevo.biologie.ens.fr\"\n";
    print Condor_SCRIPT "getenv = True\n";

    print Condor_SCRIPT "arguments = --vanilla\n";
	#print Condor_SCRIPT "transfer_files = ALWAYS\n";
    #print Condor_SCRIPT "transfer_input_files = $in\n";
    print Condor_SCRIPT "notification = NEVER\n";
    print Condor_SCRIPT "input = $in\n";
    print Condor_SCRIPT "output = $out\n";
    print Condor_SCRIPT "error = $jobDir/error_R\$(Process).txt\n";
    print Condor_SCRIPT "log = $jobDir/log_R\$(Process).txt\n";
    print Condor_SCRIPT "environment = ARG1=\$(Process)\n";

    print Condor_SCRIPT "queue $replicate\n";
    
    print "submit $Condor_jobname with $replicate replicate\n";
    
    system("condor_submit $Condor_jobname.txt")==0
    or die "failed submitting Condor job \n";
    
    #### GIVE THE JOB A SECOND TO GET SUBMITTED ####
    sleep(1);
    print "submitted.\n";
}

