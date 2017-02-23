#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;
use File::Glob ':glob';
use File::Basename;
use Getopt::Long;

my $sleep_duration=900; #Number of seconds interval to check if new jobs need to be submitted
my $maxJobsQueued=9500; #Number of maximum jobs
my $jobIntervalPrinter=10000;
my $currentWorkingDir = `pwd`;
chomp($currentWorkingDir);
print "\nWorking directory: $currentWorkingDir\n";
print "\nReading job directory..\n";

my @jobsToSubmit;
#Read file
open(INPUT, "< ./jobsToResubmit.txt") or die "Failed to open file: $!\n";
while (my $lin=<INPUT>){
    chomp($lin);
    push(@jobsToSubmit, $lin);
}
close(INPUT);
#my @jobsToSubmit=glob("./testJobs/chr*/*.sh");
print "Done reading job file!\n\n";

#for i in {1..20}; do echo -e -n '#!/bin/sh\n\nsleep 30' > test.chr1.chunk$i.sh; done

#foreach my $shellScript (@jobsToSubmit){
#    chomp($shellScript);
#    print "Script: $shellScript\n";
#}

my $lastShellScriptId = $#jobsToSubmit;
print "\nNumber of jobs detected: " . ($lastShellScriptId+1) . "\n";
print "Queue checking interval time: $sleep_duration seconds\n";
print "Maximum total number of jobs in queue: $maxJobsQueued\n\n";


my $waketime;
my $counter=0; #Counter which outputs update in counts every n-th job
while (1) {
    $waketime = time + $sleep_duration;
    &checkQueue;
    #print "blaat\n";
    sleep($waketime-time);
}


sub checkQueue{
    my $sqQueued = `squeue -u umcg-fvandijk | wc -l`;
    chomp $sqQueued;
    my $queued = ($sqQueued - 1); #total jobs in queue
    my $toSubmit = ($maxJobsQueued-$queued); #number of jobs which can be submitted
    for (my $i=1; $i<=$toSubmit; $i++){
        if (@jobsToSubmit) { #If not empty continue
            my $jobToSubmit = $jobsToSubmit[0]; #Submit first job from array
            #Submit this job now, afterwards delete from array
            my ($file,$dir,$ext) = fileparse($jobToSubmit, qr/\.[^.]*/); #extract filename, directory and extension from full path name
            chdir($dir);
            `sbatch "$file$ext"`;
            chdir($currentWorkingDir);
            #`sbatch --qos=dev --cpus-per-task=1 --time=00:01:00 $jobToSubmit`;
            shift(@jobsToSubmit); #Remove first element/job from array
        }else{
            if ($#jobsToSubmit < $toSubmit) { #If remainder of jobs in array is less than can be submitted we're done
                print "Job array empty, all jobs submitted!\nPlease be aware that some jobs might still be queued or running!!\n\n";
                exit;
            }
        }
        print "Processed $counter jobs..\n" if $counter % $jobIntervalPrinter == 0;
        $counter++;
    }
}
