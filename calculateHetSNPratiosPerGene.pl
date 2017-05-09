#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;
use List::Util qw(min max);
use List::Util qw(sum);
use List::Util qw(first);
use File::Glob ':glob';
use File::Basename;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);
use lib '/home/umcg-fvandijk/perl_modules';

my $sampleOfInterest = "20160941_510719_DNA089416_429943_5GPM1603";
my $geneOfInterest = "CTBP1";
my $coverage=100;

my $tabixPath="/apps/software/tabix/0.2.6-foss-2015b/bin/";


###Read GTF file
print "Processing GTF file..\n";
open(GTF, "< /apps/data/ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf") || die "Can't open GTF file\n";
my %genes;
my %geneNames;
my %geneIDs;
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
        my $geneName = "NA";
        if (looks_like_number($chr)) { #Check if chromosome is a number
                if ($info =~ m/gene_id "(ENSG[0-9]{1,})";.+gene_name "(.+)"; gene_sour.+/gs) { #Extract Ensembl gene ID
                    $geneID = $1;
                    $geneName = $2;
                }
                #print "$chr\t$start\t$stop\t$feature\t$geneID\n";
                if ($feature eq "gene") { #If feature is a gene, add it to hash with corresponding start-stop coordinates
                    $geneCount++; #Count number of unique genes
                    if (not exists $genes{ $geneID }) {
                        $genes{ $geneID } = "$chr:$start-$stop";
                        #if ($geneID eq "ENSG00000166340") {
                        #    print "\nENSG00000166340: $start-$stop\n\n";
                        #}
                        
                        $geneNames{ $geneID } = "$geneName";
                        $geneIDs{ $geneName } = $geneID;
                        $regionGenes{ "$chr:$start-$stop" } = $geneID;
                    }#else{
                    #   exit("Gene ID: $geneID already exists!\n");
                    #}
                }
                if ($feature eq "exon") { #If feature is exon, append the regions to corresponding Ensembl gene ID
                    if (not exists $exons{ $geneID }) {
                        $exons{ $geneID } = "$start-$stop";
                    }else{ #Append this exon to current value of existing Ensembl GeneID
                        my $val = $exons{ $geneID };
                        $val .= ",$start-$stop";
                        $exons{ $geneID } = $val;
                    }
                    $regionExons{ "$chr:$start-$stop" } = $geneID;
                }
        }
    }
}
close(GTF);
print "#Genes detected: $geneCount\n";
print "Done processing GTF file.\n\n";



my $VCFheaderLine = `zcat /groups/umcg-gdio/tmp04/projects/BioSB_5GPM/ASE/results/bins/chr4_ASVCF.vcf.gz | head -500 | grep '^#CHR'`; #Extract header line from VCF
chomp($VCFheaderLine);
my @VCFheader = split("\t", $VCFheaderLine);
my $lastSampleIdx = $#VCFheader;
my %VCFindex;
@VCFindex{@VCFheader} = (0..$#VCFheader);

#Extract sampleOfInterest index from header;
my $search = $sampleOfInterest;
my $sampleIdx = first { $VCFheader[$_] eq $search } 0..$#VCFheader;


my $geneSymbol = $geneOfInterest;
my $exonString = $exons{ $geneIDs{ $geneSymbol } }; #GeneID is retrieved from geneIDs hash, use it to retrieve exons from this gene

my @ar = split(",", $exonString); #split all exons in string
my %HoI;
my %posits;
my $lastIdx;
foreach my $exon (@ar){ #Foreach exon
    my @arr = split("-", $exon);
    my $start = $arr[0];
    my $stop = $arr[1];
    #Foreach exon query the VCF file using tabix
    my $tabixCMD="$tabixPath/tabix /groups/umcg-gdio/tmp04/projects/BioSB_5GPM/ASE/results/bins/chr4_ASVCF.vcf.gz 4:$exon";
    my $exeCMD=`$tabixCMD`;
    my @vars = split("\n", $exeCMD); #Push tabix results in array
    foreach my $lin (@vars){
        chomp($lin);
        my @array=split("\t", $lin); #Split every line using tab character
        $lastIdx = $#array;
        my $pos = $array[1];
        print $array[$sampleIdx] . "\n";
        $posits{$pos} = 1; #Ensure every SNP position is put in there once!!
        #push(@positions, $pos); #Push all exonic SNPs in array
        for (my $k=9; $k<=$lastIdx; $k++){ #Loop over all individuals, store position, individuals index and genotype/coverage in HoH
            $HoI{$pos}{$k} = $array[$k];
        }
    }
}

print "\n\n\n\n\n\n\n";

my @positions = keys %posits;
if (@positions) { #If at least one exonic SNP found in this gene
    print "$geneSymbol\tPOS\tRATIO\n";
    my @haplotype1;
    my @haplotype2;
    foreach my $snp (@positions){ #Foreach SNP select genotype from this individual
        my $indGT = $HoI{$snp}{$sampleIdx};
        print "$indGT\n";
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
                    print "$geneSymbol\t$snp\t$ratio\n";
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
        if ($hpCov >= $coverage) {
            $div = ($hp1/$hpCov);
        }
    }
    #print OUTPUT "\t$div"; #print geneRatio
}
