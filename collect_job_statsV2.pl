#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;
use File::Glob ':glob';
use File::Basename;
use Getopt::Long;

my $version = "0.1";

#perl collect_job_statsV2.pl -jobs 2030748-2046748 -outputFile logOutput.txt

#Commandline variables
my ($help, $list, $output);

#### get options
GetOptions(
                "h"             => \$help,
                "jobs=s"        => \$list, #list of jobs to generate for
                "outputFile=s"  => \$output #output file
          );
usage() and exit(1) if $help;
#Obligatory args
usage() and exit(1) unless $list;
usage() and exit(1) unless $output;

my @indices;
my $idx;
my @JobBatchArray;
my @listOfJobs;
my $Job;
my $JobBatch;

chomp($list);
chomp($output);

my $outputHeader = "AllocCPUS\tAccount\tAssocID\tAveCPU\tAvePages\tAveRSS\tAveVMSize\tBlockID\tCluster\tComment\tCPUTime\tCPUTimeRAW\tDerivedExitCode\tElapsed\tEligible\tEnd\tExitCode\tGID\tGroup\tJobID\tJobName\tLayout\tMaxPages\tMaxPagesNode\tMaxPagesTask\tMaxRSS\tMaxRSSNode\tMaxRSSTask\tMaxVMSize\tMaxVMSizeNode\tMaxVMSizeTask\tMinCPU\tMinCPUNode\tMinCPUTask\tNCPUS\tNNodes\tNodeList\tNTasks\tPriority\tPartition\tQOS\tQOSRAW\tReqCPUS\tReserved\tResvCPU\tResvCPURAW\tStart\tState\tSubmit\tSuspended\tSystemCPU\tTimelimit\tTotalCPU\tUID\tUser\tUserCPU\tWCKey\tWCKeyID\n";
open (OUTPUT, ">$output") or die "Cannot open outputfile: $output\n";
print OUTPUT $outputHeader;

#Check if input was a single job or range of jobs
if ($list =~ m/([0-9]{1,})-([0-9]{1,})/gs){ #Range of jobs
    my $first = $1;
    my $last = $2;
    for (my $k = $first; $k <= $last; $k++){
        push(@listOfJobs, $k);
    }
}elsif ($list =~ m/([0-9]{1,})/gs){ #Check if single jobID was provided
    my $singleJob = $1;
    push (@listOfJobs, $singleJob);
}else{
    #Format not recognized, throw error
    die("Format of jobIDs not recognized!\n");
}

#foreach job retrieve stats by using a single query command
    #Retrieve stats for jobID
    #collectStats($jobID);
    
    my $jobID = join(",", @listOfJobs);
    my $statsToRetrieve = "AllocCPUS,Account,AssocID,AveCPU,AvePages,AveRSS,AveVMSize,BlockID,Cluster,Comment,CPUTime,CPUTimeRAW,DerivedExitCode,Elapsed,Eligible,End,ExitCode,GID,Group,JobID,JobName,Layout,MaxPages,MaxPagesNode,MaxPagesTask,MaxRSS,MaxRSSNode,MaxRSSTask,MaxVMSize,MaxVMSizeNode,MaxVMSizeTask,MinCPU,MinCPUNode,MinCPUTask,NCPUS,NNodes,NodeList,NTasks,Priority,Partition,QOS,QOSRAW,ReqCPUS,Reserved,ResvCPU,ResvCPURAW,Start,State,Submit,Suspended,SystemCPU,Timelimit,TotalCPU,UID,User,UserCPU,WCKey,WCKeyID";

    my $cmd = "sacct --parsable2 --format $statsToRetrieve -j $jobID"; #Print output with "|" as seperator
    
    #Retrieve output from SLURM database
    my $exec = `$cmd`;
    print OUTPUT "$exec";
    #my @array = split("\n", $exec);
    #
    ##Push header into array
    #my $header = $array[0];
    #chomp $header;
    #
    ##Push normal job line into array
    #my $Job = $array[1];
    #chomp $Job;
    #my @JobArray = split(quotemeta("|"), $Job);
    #
    ##Push batch job line into array
    #my $JobBatch = $array[2];
    #if (defined $JobBatch) {
    #
    #chomp $JobBatch;
    #
    #@JobBatchArray = split(quotemeta("|"), $JobBatch);
    #
    ##Retrieve JobName idx from header
    #my @colNames = qw(JobName);
    #getColumnIdx($header, \@colNames);
    #
    #my $JobNameIdx = $indices[0];
    #my $correctJobName = $JobArray[$JobNameIdx];
    ##Splice new correct job name into Job Batch Array
    #splice(@JobBatchArray, $JobNameIdx, 1, $correctJobName);
    #}else{
    #    next;
    #}
    #
    #
    ##Print new output into file
    #my $lastIdx = $#JobBatchArray;
    #for (my $i = 0; $i < $lastIdx; $i++){
    #    print OUTPUT $JobBatchArray[$i] . "\t";
    #}
    #print OUTPUT $JobBatchArray[$lastIdx] . "\n";

    undef(@JobBatchArray);
    undef($Job);
    undef($JobBatch);


#Close output file
close(OUTPUT);


##########SUBS#########SUBS###########SUBS##########

#Collect stats per job
#sub collectStats {
#    my $jobID = shift;
#    my $statsToRetrieve = "AllocCPUS,Account,AssocID,AveCPU,AvePages,AveRSS,AveVMSize,BlockID,Cluster,Comment,CPUTime,CPUTimeRAW,DerivedExitCode,Elapsed,Eligible,End,ExitCode,GID,Group,JobID,JobName,Layout,MaxPages,MaxPagesNode,MaxPagesTask,MaxRSS,MaxRSSNode,MaxRSSTask,MaxVMSize,MaxVMSizeNode,MaxVMSizeTask,MinCPU,MinCPUNode,MinCPUTask,NCPUS,NNodes,NodeList,NTasks,Priority,Partition,QOS,QOSRAW,ReqCPUS,Reserved,ResvCPU,ResvCPURAW,Start,State,Submit,Suspended,SystemCPU,Timelimit,TotalCPU,UID,User,UserCPU,WCKey,WCKeyID";
#
#    my $cmd = "sacct --parsable2 --format $statsToRetrieve -j $jobID"; #Print output with "|" as seperator
#    
#    #Retrieve output from SLURM database
#    my $exec = `$cmd`;
#    my @array = split("\n", $exec);
#    
#    #Push header into array
#    my $header = $array[0];
#    chomp $header;
#    
#    #Push normal job line into array
#    my $Job = $array[1];
#    chomp $Job;
#    my @JobArray = split(quotemeta("|"), $Job);
#    
#    #Push batch job line into array
#    my $JobBatch = $array[2];
#    if (defined $JobBatch) {
#    
#    chomp $JobBatch;
#    
#    @JobBatchArray = split(quotemeta("|"), $JobBatch);
#    
#    #Retrieve JobName idx from header
#    my @colNames = qw(JobName);
#    getColumnIdx($header, \@colNames);
#    
#    my $JobNameIdx = $indices[0];
#    my $correctJobName = $JobArray[$JobNameIdx];
#    #Splice new correct job name into Job Batch Array
#    splice(@JobBatchArray, $JobNameIdx, 1, $correctJobName);
#    }else{
#        next;
#    }
#    return(@JobBatchArray);
#    
#}

#Grep column indices and push into array
sub getColumnIdx {
    my $header = shift;
    my ($colnames) = @_;
    my @headerarray = split(quotemeta("|"), $header); #Split header on space!
    undef(@indices);
    foreach my $columnName (@$colnames){
        my( $idx )= grep { $headerarray[$_] eq $columnName } 0..$#headerarray;
        push(@indices, $idx);
    }
    return(@indices);
}


#Usage of software
sub usage {
        print <<EOF;
#########################################################################################################

#########################################################################################################

Usage: ./CoNVaDING-$version.pl
-h\t\t\tThis manual.
-jobs\t\t\tThe (range of) jobIDs to generate list for, can either be:
\t\t\t\tSingle jobID: <jobID>\t\t\t\texample: 34567
\t\t\t\tRange of jobIDs: <first jobID>-<last jobID>\texample: 34567-34590
-outputFile\t\tOutput file to write results to

#########################################################################################################
EOF
 
}