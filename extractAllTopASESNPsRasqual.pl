#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;
use File::Glob ':glob';
use File::Basename;
use Getopt::Long;
use Term::ProgressBar;

#module load Term-ProgressBar/2.17-foss-2015b

#Retrieve all files
my @files = glob("/groups/umcg-bios/tmp04/projects/ASE_GoNL/rasqual/GoNL_WGS_meta-exons/results/Rasqual/metaExon/chr*/*rasqual.output.txt");

my $progress_bar = Term::ProgressBar->new($#files);

#Open output file handle
open(OUTPUT, "> overlapRasqualClinVarExac.txt") or die("Unable to open comparison file");

my $count=0;
foreach my $file (@files){
    #open file
    open(FILE, "< $file") or die("Unable to open chunk file: $file"); #Read file
    my %rank;
    my %fullLine;
    while (my $line = <FILE>){
        chomp $line;
        my @arr=split("\t", $line);
        my $rsID = $arr[1];
        my $chr = $arr[2];
        my $pos = $arr[3];
        my $loc = "$chr\t$pos";
        my $chiSq = $arr[10];
        #push location and chiSq in hash to sort later
        $rank{ $loc } = $chiSq;
        #push location and complete line in hash for later use, after ranking
        $fullLine{ $loc } = $line;
    }
    close(FILE);
    #Reverse sort the list to put highest chiSquare first
    my @locs = reverse sort { $rank{$a} <=> $rank{$b} } keys(%rank); #Reverse sort the values, making highest chiSq first
    #Print corresponding line to highest ranked event
    my $toPrint = $fullLine{ $locs[0] };
    print OUTPUT "$toPrint\n";
    $progress_bar->update($count);
    $count++;
}

close(OUTPUT);
