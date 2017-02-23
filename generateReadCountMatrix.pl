#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;
use File::Glob ':glob';
use File::Basename;
use List::Util qw(sum);
use Getopt::Long;


q#Read eQTL file
open(EQTL, "</groups/umcg-bios/tmp02/umcg-fvandijk/projects/ASEconcordance/readCounts/exon_level_eQTLs_Top_effects.exonicSNPsOnly.txt") or die("Unable to open input file"); #Read file
my @EQTLfile=<EQTL>;
#my $string = <FILE>;
close(EQTL);

my $lastEQTLIdx = $#EQTLfile;
my %eQTLs;
for (my $w = 1; $w <= $lastEQTLIdx; $w++){
    my $line = $EQTLfile[$w];
    chomp $line;
    my @array = split("\t", $line);
    #PValue  SNPName SNPChr  SNPChrPos       GeneId  GeneChr GeneCenterChrPos        SNPType AlleleAssessed  OverallZScore   HGNCName        FDR
    my $chr = $array[2];
    my $pos = $array[3];
    my $type = $array[7];
    my $alAs = $array[8];
    my $zScore = $array[9];
    my $all1;
    my $all2;
    my $key = "$chr\_$pos";
    if ($type =~ m/(.+)\/(.+)/gs){
        $all1 = $1;
        $all2 = $2;
    }
    my $updatedZscore = $zScore;
    #print "$chr\t$pos\t$type\t$alAs\t$zScore\t$all1\t$all2\n";
    #Check if SNPType alleles match ref and alt used, otherwise swap
    if (exists $wgs{$key} ) {
        if ($all2 eq $alAs) { #If allele assessed matches alternative, swap effect direction
            #swap effect direction
            if ($zScore =~ m/^-(.+)/gs){
                $updatedZscore = $1;
            }elsif ($zScore =~ m/^[0-9]{1,}.+/gs){
                $updatedZscore = "-$zScore";
            }
        }
    }
    my $val = "$updatedZscore";
    $eQTLs{ $key } = $val; #Push zScore in eQTL hash
}


#Iterate over counts VCF and generate matrix/table
open(COUNTS, "<All.allelicDepths.vcf") or die("Unable to open input file"); #Read file
my @COUNTS=<COUNTS>;
close(COUNTS);

my @header;
my @wgsSamples;
my @output;
my $lastIdxCounts = $#COUNTS;
for (my $n=0; $n<= $lastIdxCounts; $n++){
    my $string = $COUNTS[$n];
    chomp $string;
    if ($string =~ m/^#CHR.+/gs) { #Headerline
        my @array = split("\t", $string);
        my $lastIdx = $#array;
        for (my $b=9; $b <= $lastIdx; $b++){
            push(@wgsSamples, $array[$b]);
            my $gonlID = $samples{$array[$b]};
            push(@header, "$gonlID\_REF\t$gonlID\_ALT");
        }
        my $headerSampleNames = join("\t", @header);
        push(@output, "CHR\tPOS\tID\tREF\tALT\t$headerSampleNames\tTOTAL_REF\tTOTAL_ALT\tRATIO\tZSCORE");
    }
    if ($string !~ m/^#.+/gs) {
        my $count=0;
        my @array = split("\t", $string);
        my $lastIdx = $#array;
        my $chr = $array[0];
        my $pos = $array[1];
        my $id = $array[2];
        my $ref = $array[3];
        my $alt = $array[4];
        my $key = "$chr\_$pos";
        my @altAll = split(",", $alt);
        
        if ($key ne "#CHROM_POS") {

        
        my $searchAlt = $wgsAltAll{$key};
        my %index;
        @index{@altAll} = (0..$#altAll); #Retrieve index of alternative allele
        my $index;
        if (exists $index{$searchAlt}) {
            $index = $index{$searchAlt};
        }else{
            $index = 0;
        }
        
        #if (exists $index{$searchAlt}) {
        #    #code
        #}else{
        #    print "KEY: $key\n";
        #}
        
        
        #print "$chr $pos $ref $alt $key $searchAlt Index:$index \n";
        
        #Iterate over samples, retrieve counts of alt allele
        my @co;
        my @ref;
        my @alt;
        push(@co, "$chr\t$pos\t$id\t$ref\t$searchAlt");
        for (my $a=9; $a<=$lastIdx; $a++){
            my @sampleInfo=split(":", $array[$a]); #Second element is allelic depth
            my @adArray=split(",", $sampleInfo[1]);
            my $refCount = $adArray[0];
            my $altCount = $adArray[$index+1]; #First alternative, while reference is always displayed first. So index + 1.
            
            my $sample = $wgsSamples[$a-9];
            #print "SAMPLE: $sample\n";
            if (exists $samples{$sample}) {

                my $gonlID = $samples{$sample};
                my $gonlVCFindex = $indexHash{$gonlID};
                my $line = $wgs{$key};
                my @wgsArr = split("\t", $line);
    #            if (defined $gonlVCFindex && length $gonlVCFindex > 0) {
    #                #code
    #            }else{
    #                print $string . "\n$gonlID\n";
    #            }
                
                #DEBUGGING PURPOSES
                #if ($key eq "1_6640116" && $sample eq "BC48DMACXX-4-2") { #gonl-145a
                #    print "$key\t$sample\t$gonlID\t$refCount\t$altCount\n$line\n";
                #}
        
                my $sampleGT = $wgsArr[$gonlVCFindex];
                
                if ($sampleGT eq "1|0" || $sampleGT eq "0|1") { #Check if genotype is heterozygous;
                    #code
                    push(@ref, $refCount);
                    push(@alt, $altCount);
                }else{
                    $refCount = "NA";
                    $altCount = "NA";
                }
            }else{
                $refCount = "NA";
                $altCount = "NA";
            }
            push(@co, "$refCount\t$altCount");
            #print "R:$refCount A:$altCount S:$sample G:$gonlID GI:$gonlVCFindex GT:$sampleGT\n";
        }
        my $sumRef;
        my $sumAlt;
        my $total;
        if (@ref && @alt) {
            $sumRef = sum(@ref);
            $sumAlt = sum(@alt);
            $total = ($sumRef + $sumAlt);
        }else{
            $sumRef = "0";
            $sumAlt = "0";
            $total = "0";
        }
        
        my $ratio = "NA";
        if ($total ne 0 && $sumAlt ne 0) {
            $ratio = ($sumRef/$total);
        }
        my $zScore = "NA";
        if (exists $eQTLs{$key}) {
            $zScore = $eQTLs{$key}
        }

        push(@co, "$sumRef\t$sumAlt\t$ratio\t$zScore");
        my $toAdd = join("\t", @co);
        push(@output, $toAdd);
        }
    }
    #print "\n";

}

open(OUTPUT, ">countsTable.txt") or die("Unable to open output file"); #

foreach my $ele(@output){
    print OUTPUT "$ele\n";
}

close(OUTPUT);
