#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;
use List::Util qw(min max);
use List::Util qw(sum);
use File::Glob ':glob';
use File::Basename;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);
use lib '/home/umcg-fvandijk/perl_modules';

#time perl calculateASEv3.pl --chr 1 --gtfFile /apps/data/ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf --phasedVCF TEST1chr1_ASVCF.vcf --annotatedVCF chr1_ASVCF.snpEff.caddsnv.exac.gonl.1kg.clinvar.snpEff.cgd.vcf --outputFile output.txt


####Variables
my $help;
my $chrIN;
my $gtfFile;
my $phasedVCF;
my $annotatedVCF;
my $coverage;
my $output;
my $counter=0;
my $lcounter=0;
my $countIntervalPrinter=100;
my $tabixPath="/apps/software/tabix/0.2.6-foss-2015b/bin/";

####Get options
GetOptions(
                "h"                     => \$help,
                "chr=s"                 => \$chrIN,
                "gtfFile=s"             => \$gtfFile,
                "phasedVCF=s"           => \$phasedVCF,
                "annotatedVCF=s"        => \$annotatedVCF,
                "coverage=s"            => \$coverage,
                "outputFile=s"          => \$output
);
usage() and exit(1) if $help;
####Obligatory args
usage() and exit(1) unless $gtfFile;
usage() and exit(1) unless $phasedVCF;
usage() and exit(1) unless $annotatedVCF;
usage() and exit(1) unless $coverage;
usage() and exit(1) unless $output;

###Print given variables from user input
print "\n\nParameters used for this analysis:\n";
print "chr: $chrIN\n";
print "gtfFile: $gtfFile\n";
print "phasedVCF: $phasedVCF\n";
print "annotatedVCF: $annotatedVCF\n";
print "coverage: $coverage\n";
print "outputFile: $output\n\n\n";

###Read GTF file
print "Processing GTF file..\n";
open(GTF, "< $gtfFile") || die "Can't open file: $gtfFile\n";
my %genes;
my $geneCount=0;
my %exons;
my %regionGenes;
my %regionExons;
while (my $line=<GTF>){ #Read GTF file line by line
    chomp($line);
    if ($line !~ m/^#.+/gs) { #If not header line process further
        my @array = split("\t", $line);
        my $chr = $array[0];
        $chr =~ s/X/23/gs;
        $chr =~ s/Y/24/gs;
        $chr =~ s/MT/25/gs;
        my $feature = $array[2];
        my $start = $array[3];
        my $stop = $array[4];
        my $info = $array[8];
        my $geneID = "NA";
        if (looks_like_number($chr)) { #Check if chromosome is a number
            if ($chrIN == $chr) { #If chromosome from input matches the one in GTF file
                if ($info =~ m/gene_id "(ENSG[0-9]{1,})";.+/gs) { #Extract Ensembl gene ID
                    $geneID = $1;
                }
                #print "$chr\t$start\t$stop\t$feature\t$geneID\n";
                if ($feature eq "gene") { #If feature is a gene, add it to hash with corresponding start-stop coordinates
                    $geneCount++; #Count number of unique genes
                    if (not exists $genes{ $geneID }) {
                        $genes{ $geneID } = "$start-$stop";
                        $regionGenes{ "$start-$stop" } = $geneID;
                    }else{
                        exit("Gene ID: $geneID already exists!\n");
                    }
                }
                if ($feature eq "exon") { #If feature is exon, append the regions to corresponding Ensembl gene ID
                    if (not exists $exons{ $geneID }) {
                        $exons{ $geneID } = "$start-$stop";
                    }else{ #Append this exon to current value of existing Ensembl GeneID
                        my $val = $exons{ $geneID };
                        $val .= ",$start-$stop";
                        $exons{ $geneID } = $val;
                    }
                    $regionExons{ "$start-$stop" } = $geneID;
                }
            }
        }
    }
}
close(GTF);
print "#Genes detected: $geneCount\n";
print "Done processing GTF file.\n\n";

# Foreach gene
## Per exon
### het SNPs
#### check if within ratio 0.05-0.95  .. calc: ref/(ref+alt)
##### Per haplotype sum reads .. output


###Proces phased VCF file
print "\nProcessing phased VCF file..\n";
open(OUTPUT, "> $output") || die "Can't open file: $output\n";
my %SNPs;
my $SNPcount=0;
my $SNPexonic=0;
my $headerLine = `zcat $phasedVCF | head -500 | grep '^#CHR'`; #Extract header line from VCF
chomp($headerLine);
my @header = split("\t", $headerLine);
my $lastSampleIdx = $#header;
print OUTPUT "EnsemblGeneID";
for (my $j=9; $j <= $lastSampleIdx; $j++){
    print OUTPUT "\t" . $header[$j];
}
print OUTPUT "\tCOUNT\n";
my $curGene;
#Check if position is within a specific gene
foreach my $key (sort keys %regionGenes){
    my $count=0;
    $curGene = $regionGenes{ $key };
    my $exonString = $exons{ $regionGenes{ $key } }; #Extract exon regions by using the key (gene ID)
    my @ar = split(",", $exonString);
    my %HoI;
    my %posits;
    my $lastIdx;
    foreach my $exon (@ar){ #Foreach exon
        my @arr = split("-", $exon);
        my $start = $arr[0];
        my $stop = $arr[1];
        #Foreach exon query the VCF file using tabix
        my $tabixCMD="$tabixPath/tabix $phasedVCF $chrIN:$exon";
        my $exeCMD=`$tabixCMD`;
        my @vars = split("\n", $exeCMD); #Push tabix results in array
        foreach my $lin (@vars){
            chomp($lin);
            my @array=split("\t", $lin); #Split every line using tab character
            $lastIdx = $#array;
            my $pos = $array[1];
            $posits{$pos} = 1; #Ensure every SNP position is put in there once!!
            #push(@positions, $pos); #Push all exonic SNPs in array
            $SNPexonic++;
            for (my $k=9; $k<=$lastIdx; $k++){ #Loop over all individuals, store position, individuals index and genotype/coverage in HoH
                $HoI{$pos}{$k} = $array[$k];
            }
        }
    }
    
    my @positions = keys %posits;
    if (@positions) { #If at least one exonic SNP found in this gene
        print OUTPUT "$curGene"; #Print gene name
        for (my $i=9; $i <= $lastIdx; $i++){ #Iterate over all samples
            my @haplotype1;
            my @haplotype2;
            foreach my $snp (@positions){ #Foreach SNP select genotype from this individual
                my $indGT = $HoI{$snp}{$i};
                chomp $indGT;
                my @arr = split(":", $indGT);
                my $gt = $arr[0]; #Genotype of this individual
                my $cov = $arr[1]; #Coverage of this individual
                if ($gt eq "0|1" || $gt eq "1|0") { #Check if heterozygous SNPs
                    my @covArr = split(",", $cov);
                    my $ref = $covArr[0]; #Reference allele counts
                    my $alt = $covArr[1]; #Alternative allele counts
                    my $tot = ($ref+$alt); #Total depth
                    if ($tot > 9) { #If both ref and alt are above threshold, do calculation of ratio and push in correct haplotype array for this gene
                        my $ratio = ($ref/$tot); #Ratio
                        if ($ratio > 0.05 && $ratio < 0.95) { #Ratio passes threshold, determine haplotype and add to count
                            if ($gt eq "0|1") { #Push ref and alt in hap1 and hap2 arrays respectively
                                push(@haplotype1, $ref);
                                push(@haplotype2, $alt);
                            }else{ #Else phased other way, so swap the ref and alt counts into correct array
                                push(@haplotype1, $alt);
                                push(@haplotype2, $ref);
                            }
                        }
                    }
                }
            }
            
            #Sum reads on each haplotype and calculate haplotype ratio
            my $hp1 = "NA";
            my $hp2 = "NA";
            my $div = "NA";
            if (@haplotype1) { #If haplotype1 array counts values, we can calculate
                $hp1 = sum(@haplotype1);
                $hp2 = sum(@haplotype2);
                my $hpCov = ($hp1 + $hp2);
                if ($hpCov > $coverage) {
                    $div = ($hp1/$hpCov);
                    $count++;
                }
            }
            #print "$k:$i:$hp1:$hp2:$div\n";
            #print OUTPUT "\t$hp1:$hp2:$div";
            print OUTPUT "\t$div";
        }
        print OUTPUT "\t$count\n";
    }
    print "Processed $lcounter genes..\n" if $lcounter % $countIntervalPrinter == 0;
    $lcounter++;
}

close(OUTPUT);

print "Total exonic SNPs: $SNPexonic\n";
print "Done processing phased VCF file.\n\n";



####Usage
sub usage {
        print <<EOF;
        
        
###################################################
This script ....
###################################################
Usage: calculateASEv4.pl
Author: f.van.dijk @ UMCG
############## Obligatory parameters ##############
--help\t\tThis manual
--chr\t\tChromosome to process.
--gtfFile\tInput GTF file to merge genes from.
--phasedVCF\tbla
--annotatedVCF\tbla
--outputFile\tOutput file containing chunks.
###################################################
EOF
 
}