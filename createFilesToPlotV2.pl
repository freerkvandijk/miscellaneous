#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;
use List::Util qw(min max);
use List::Util qw(sum);
use File::Glob ':glob';
use File::Basename;
use Getopt::Long;

my $cov=10;


for (my $i=1; $i<=22; $i++){
    my $input="output.chr$i.cov$cov.txt";
    my $output="output.chr$i.cov$cov.toPlot.txt";


    my $minSamples = 100;


    #Open file
    print "Processing chr$i ..\n";
    open(INPUT, "< $input") || die "Can't open file: $input\n";
    open(OUTPUT, "> $output") || die "Can't open outputfile: $output\n";
    print OUTPUT "Sample\tEnsemblGeneID\tRatio\n";
    while (my $line=<INPUT>) {
        chomp($line);
        my @array = split("\t", $line);
        my $lastIdx = $#array;
        my $geneID = $array[0];
        my $headerLine = `head -1 $input`;
        chomp $headerLine;
        my @header = split("\t", $headerLine);
        if ($line =~ m/Ensem.+/gs) { #Header line
            #;
        }else{ #Gene line, create and print output
            #EnsemblGeneID   AC1C40ACXX-1-18_BD2D5MACXX-6-18
            my $count = $array[$lastIdx];
            if ($count > $minSamples-1) { #If count higher than user specified, write to output, else ignore this line
                for (my $j=1; $j<=$lastIdx-1; $j++){ #Loop over all sample ratios
                    my $sam=$header[$j];
                    my $ratio=$array[$j];
                    print OUTPUT "$sam\t$geneID\t$ratio\n";
                }
            }
        }
    }
    close(INPUT);
    close(OUTPUT);
}

