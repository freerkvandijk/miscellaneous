#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;
use List::Util qw(min max);
use File::Glob ':glob';
use File::Basename;
use Getopt::Long;


for(my $i=1; $i<=22; $i++){
my $inputFile="/groups/umcg-bios/tmp03/projects/rasqualMaskedBAMs/inputVCFs/unpacked/BIOS_freeze2.chr$i.concatenated.shapeit.phased.sorted.5samplesRemoved.vcf";
my $outputFile="/groups/umcg-bios/tmp03/projects/rasqualMaskedBAMs/inputVCFs/unpacked/BIOS_freeze2.chr$i.concatenated.shapeit.phased.sorted.5samplesRemoved.dupEntriesRemoved.vcf";
##
print "Deduplicating file: $inputFile..\n";
open(VCF, "< $inputFile") || die "can't open input file!\n";
open(OUTPUT, "> $outputFile") || die "can't open output file!\n";
my %vcf;
while (my $lin = <VCF>){
    chomp($lin);
    if ($lin !~ m/^#.+/gs){
        my @array = split("\t", $lin);
        my $chr = $array[0];
        my $pos = $array[1];
        $vcf{ $pos } = $lin;
    }else{#Header line, just write away as is
        print OUTPUT "$lin\n";
    }
}
close(VCF);
print "Done reading input file\n";

print "Sorting and writing output file..\n";
#Print sorted output VCF
for my $key (sort {$a<=>$b} keys %vcf){
    print OUTPUT $vcf{ $key } . "\n";
}
close(OUTPUT);
print "Done deduplicating file: $inputFile, results can be found in: $outputFile\n";
}