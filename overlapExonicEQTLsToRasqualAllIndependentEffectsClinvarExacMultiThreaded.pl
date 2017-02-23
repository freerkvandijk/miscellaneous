#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;
use File::Glob ':glob';
use File::Basename;
use Getopt::Long;
#use Term::ProgressBar;
use lib "/home/umcg-fvandijk/perl_modules/threads-2.09/lib/";
use threads;
use lib "/home/umcg-fvandijk/perl_modules/threads-shared-1.52/lib/threads/";
use threads::shared;

####LOAD THIS MODULE FIRST####
#module load Perl/5.22.0-foss-2015b-bare
my $num_of_threads = 24;

#module load tabix
#module load Term-ProgressBar/2.17-foss-2015b

#Read ClinVar
#my $clinvar = "/apps/data/ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar_20161101.vcf.gz";
#open(CV, "gunzip -c $clinvar |") || die "can't open pipe to VCF file $clinvar";
#open(CV, "<clinvar_20161101.vcf") || die "can't open pipe to VCF file clinvar_20161101.vcf";

my %clnsig;
my %exac;
my %exacAF;
share(%clnsig);
share(%exac);
share(%exacAF);

#my $totalLines = `zcat $clinvar | wc -l`;
#my $progress_bar = Term::ProgressBar->new($totalLines);

print "Reading Clinvar file..\n";

my $c=0;

# use the initThreads subroutine to create an array of threads.
my @threads = initThreads();

# Loop through the array:
my $chunkNum=0;

foreach(@threads){
    my $twoDig = sprintf ("%02d",$chunkNum);
    # Tell each thread to perform our 'doOperation()' subroutine.
    $_ = threads->create(\&doOperation);
    $chunkNum++;
}

# This tells the main program to keep running until all threads have finished.
foreach(@threads){
    $_->join();
}


print "Done reading ClinVar.\n";

#Loop over Rasqual top gene effects
my $Rr="/groups/umcg-bios/tmp04/projects/ASE_GoNL/rasqual/GoNL_WGS_meta-exons_allSNPs/comparisonExonicEQTLsToRasqualAllIndependentEffects.1e-8.txt";
open(RAS, "< $Rr") || die "can't open file $Rr";

#Open output file
open(OUTPUT, "> /groups/umcg-bios/tmp04/projects/ASE_GoNL/rasqual/GoNL_WGS_meta-exons_allSNPs/overlapExonicEQTLsToRasqualAllIndependentEffects.1e-8.ClinvarExac.txt") or die("Unable to open output overlap file");

my $count=0;
my $disCount=0;
my $tot=0;
my $exacCount=0;
my $exacDisCount=0;

print "Processing RASQUAL top SNPs..\n";
while (my $lin = <RAS>){
    $tot++;
    chomp $lin;
    my @arr= split("\t", $lin);
    my $chr = $arr[0];
    my $pos = $arr[1];
    my $loc = "$chr\t$pos";
    my $chiSq = $arr[11]; #Column 12 in this case
    #if ($chiSq >= 5) {
    print OUTPUT "$lin\t";
    if (exists $clnsig{ $loc }) {
        $count++;
        print OUTPUT $clnsig{ $loc } . "\t";
    }else{
        $disCount++;
        print OUTPUT "NA\t";
    }
    
    if (exists $exac{ $loc }) {
        my $af = $exacAF{ $loc };
        print OUTPUT "YES\t$af\n";
        $exacCount++;
    }else{
        print OUTPUT "NO\tNA\n";
        $exacDisCount++;
    }
    #}
}
close(RAS);
close(OUTPUT);
print "Done processing RASQUAL top snps file.\n";

print "Total number of events: $tot\n";
print "Overlap: $count\n";
print "Non-overlap: $disCount\n";
print "ExAC_Overlap: $exacCount\n";
print "ExAC_Non-overlap: $exacDisCount\n";



#####SUBS#####
sub initThreads{
    # An array to place our threads in
    my @initThreads;
    for(my $i = 1;$i<=$num_of_threads;$i++){
        push(@initThreads,$i);
    }
    return @initThreads;
}

sub doOperation{
    #my $chunkNum = shift;
    # Get the thread id. Allows each thread to be identified.
    my $id = threads->tid();
    my $n = ($id-1);
    my $chunkNum = sprintf ("%02d",$n);
    print "Processing file: clinvar_20161101.chunk$chunkNum..\n";
    open(CV, "< clinvar_20161101/clinvar_20161101.chunk$chunkNum") || die "can't open file clinvar_20161101.chunk$chunkNum";
    while (my $line = <CV>) {
        chomp($line);
        if ($line !~ m/^#/gs) {
            my @array = split("\t", $line);
            my $chr = $array[0];
            my $pos = $array[1];
            my $rsID = $array[2];
            my $ref = $array[3];
            my $alt = $array[4];
            my $info = $array[7];
            my $loc = "$chr\t$pos";
            #extract clinical significance from ClinVar
            if ($info =~ m/.+CLNSIG=([0-9]){1,};.+/gs){
                my $clnsig = $1;
                #print "$loc\t$clnsig\n";
                $clnsig{ $loc } = $clnsig;
            }
            #Extract relevant info from ExAC VCF
            my $grepExAC = `tabix /apps/data/ExAC/release0.3.1/ExAC.r0.3.1.sites.vep.vcf.gz $chr:$pos-$pos`;
            if (defined $grepExAC && $grepExAC ne "") {
                $exac{ $loc } = $line;
                my @array = split("\t", $grepExAC);
                my @altArr = split(",", $array[4]);
                my( $index )= grep { $altArr[$_] eq $alt } 0..$#altArr;
                if ($grepExAC =~ m/.+;AF=(.+);AN=.+/gs) {
                    my $AF = $1;
                    my @alt = split(",", $AF);
                    #print "AF: $AF\n";
                    if (defined $index && $index ne "") {
                        $exacAF{ $loc } = $alt[$index];
                    }else{
                        $exacAF{ $loc } = "NA";
                    }
                }
            }
        }
        #$progress_bar->update($c);
        $c++;
    }
    close(CV);
    
    # Inform us that the thread is done and exit the thread.
    print "Thread $id done processing file: clinvar_20161101.chunk$chunkNum!\n";
    threads->exit();
}
