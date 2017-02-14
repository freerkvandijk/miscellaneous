#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;
use List::Util qw(min max);
use File::Glob ':glob';
use File::Basename;
use Getopt::Long;

my $chr = $ARGV[0];

my %incl;
open(IN, "< all_positions/all_positions_chr$chr.txt") || die "can't open file all_positions/all_positions_chr$chr.txt";
while (my $line=<IN>){
    chomp($line);
    $incl{ $line } = 1;
}
close(IN);


my %vcf;
open(VCF, "< VCFs/unpacked/TESTchr$chr\_ASVCF.vcf") || die "can't open file VCFs/unpacked/chr$chr\_ASVCF.vcf";
while (my $lin = <VCF>){
    chomp($lin);
    if ($lin !~ m/^#.+/gs){
        my @array = split("\t", $lin);
        my $chr = $array[0];
        $chr =~ s/X/23/gs;
        $chr =~ s/Y/24/gs;
        my $pos = $array[1];
        if (exists $incl{ "$chr\t$pos" }){
            $vcf{ $pos } = $lin;
        }
    }
}
close(VCF);


my %VCFanno;
open(VCFANNO, "< VCFs/unpacked/chr$chr\_ASVCF.snpEff.caddsnv.exac.gonl.1kg.clinvar.snpEff.cgd.vcf") || die "can't open file VCFs/unpacked/chr$chr\_ASVCF.snpEff.caddsnv.exac.gonl.1kg.clinvar.snpEff.cgd.vcf";
while (my $li = <VCFANNO>){
    chomp($li);
    if ($li !~ m/^#.+/gs){
        my @array = split("\t", $li);
        my $chr = $array[0];
        $chr =~ s/X/23/gs;
        $chr =~ s/Y/24/gs;
        my $pos = $array[1];
        my $info = $array[7];
        if (exists $incl{ "$chr\t$pos" }){
            $VCFanno{ $pos } = $info;
        }
    }
}
close(VCF);


my @infoToPrint = qw(AC AF AN DP NUMALT set ANN LOF NMD CADD CADD_SCALED EXAC_AF EXAC_AC_HOM EXAC_AC_HET GoNL_GTC GoNL_AF Thousand_Genomes_AF CLINVAR_CLNSIG CLINVAR_CLNALLE HGNC_ID ENTREZ_GENE_ID CGDCOND CGDINH CGDGIN CGDAGE ALLELIC_CONDITIONS MANIFESTATION_CATEGORIES INTERVENTION_CATEGORIES COMMENTS INTERVENTION_RATIONALE REFS);
my $headerInfo = join("\t", @infoToPrint);
my $header = "CHR\tPOS\tRSID\tREF\tALT\ttotHomRef\ttotHomAlt\ttotHet\ttotHetImbalance\tcovHet10\tcovHet20\tcovHet30\t$headerInfo";

#Extract sample names from VCF header
my $headerLineVCF=`grep '^#CH' VCFs/unpacked/chr$chr\_ASVCF.vcf`;
chomp($headerLineVCF);
my @headerVCFarray = split("\t", $headerLineVCF);
for (my $n=9; $n<=$#headerVCFarray; $n++){
    $header .= "\t" . $headerVCFarray[$n];
}

open(OUTPUT, "> TESToutput.annotated.chr$chr.txt") || die "can't open output file!\n";
print OUTPUT "$header\n";

for my $key (sort {$a<=>$b} keys %vcf){
    my @arr=split("\t", $vcf{ $key });
    my $lastIdx = $#arr;
    my @result;
    my $chr = $arr[0];
    my $pos = $arr[1];
    for (my $j=0;$j<=4;$j++){
        push(@result, $arr[$j]);
    }

    my $numHet10=0;
    my $numHet20=0;
    my $numHet30=0;
    my $numHomRef=0;
    my $numHomAlt=0;
    my $numHet=0;
    my $numHetImB=0;

    #loop through all individuals
    my @individuals;
    for (my $k=9; $k<=$lastIdx; $k++){
        my $ind= $arr[$k];
        my @indArray=split(":",$ind);
        my $gt = $indArray[0];
        my $counts = $indArray[1];
        my @countsArray = split(",", $counts);
        my $ref = $countsArray[0];
        my $alt = $countsArray[1];
        my $total = ($ref + $alt);
        my $ratio;
        my $absRatio;
        if ($total > 0) {
            $ratio = ($ref/$total);
            $absRatio = abs($ratio);
        }else{
            $ratio = "NA";
            $absRatio=1;
        }
        my $indToPrint = "$gt:$counts:$ratio";
        push(@individuals, $indToPrint);
        #Check genotype, if het and ratio between 0.6 and 0.95 there is imballance, start counting the coverage thresholds
        if ($gt eq "0|0") {
            $numHomRef++;
        }elsif($gt eq "1|1"){
            $numHomAlt++;
        }elsif($gt eq "0|1" || $gt eq "1|0"){
            $numHet++;
            if ($absRatio >= 0.60 && $absRatio <= 0.95) {
                $numHetImB++;
                if ($total > 29) {
                    $numHet30++;
                }elsif($total > 19){
                    $numHet20++;
                }elsif($total > 9){
                    $numHet10++;
                }
            }
        }else{
            print "Unknown genotype at position: $chr:$pos-$pos\n";
            exit(1);
        }
    }

    ##Search for ClinVar status in clinvar VCF
    #my $grepClinVar = `tabix /groups/umcg-bios/tmp03/projects/rasqualMaskedBAMs/overlapClinVarEtc/clinvar/clinvar_20161101.chr$chr.vcf.gz $chr:$pos-$pos`;
    #my $clnsig = "NA";
    #if (defined $grepClinVar && $grepClinVar ne "") {
    #    chomp($grepClinVar);
    #    #Extract ClinVar disease status
    #    if ($grepClinVar =~ m/.+CLNSIG=([0-9]){1,};.+/gs){
    #        $clnsig = $1;
    #        #print "$loc\t$clnsig\n";
    #    }
    #}
    
    
    ##Search annotation from CMDlineAnnotator files
    my @inf = split(";", $VCFanno{ $pos });
    my $infLastIdx = $#inf;
    
    my %infos;
    foreach my $ele (@inf){
        chomp $ele;
        my @ar=split("=", $ele);
        if ($ar[0] ne "DB") {
            my $key = $ar[0];
            my $val = $ar[1];
            if (defined $key && $key ne "") {
                if (defined $val && $val ne ""){
                    #if ($key eq "CLINVAR_CLNSIG") {
                    #    my @arr = split("|", $val);
                    #    my %vals = map { $_ => 1 } @arr;
                    #    if (exists ($vals{ 5 })) {
                    #        $val = 5;
                    #    }else{
                    #        $val = $arr[0]; #Value is first status from array
                    #    }
                    #}
                    $infos{ $key } = $val; #Push all info annotation with corresponding value in hash
                }else{
                    $infos{ $key } = "NA";
                }
            }
        }
    }
    
    #Make string from all annotation
    my $annoToPrint;
    foreach my $el (@infoToPrint){
        my $toPrint = "NA";
        if (exists $infos{ $el }) {
            $toPrint = $infos{ $el };
            #print "El: $el\t" . $infos{ $el } . "\t";
            #print "$toPrint\t";
        }
        $annoToPrint .= "\t$toPrint";
    }
            
    #Push all other results in result array
    my $merge= "$numHomRef\t$numHomAlt\t$numHet\t$numHetImB\t$numHet10\t$numHet20\t$numHet30" . "$annoToPrint";
    push(@result, $merge);
    my $allIndvsToPrint = join("\t", @individuals);
    push(@result, $allIndvsToPrint);
    my $resultToString = join("\t", @result);
    print OUTPUT "$resultToString\n";

}

close(OUTPUT);

