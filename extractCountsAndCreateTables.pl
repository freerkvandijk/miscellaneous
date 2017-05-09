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

my $minSamples = 100;
my $cov=30;
my $chrIN=11;
#my $chrIN = 1;
my $tabixPath="/apps/software/tabix/0.2.6-foss-2015b/bin/";




#/groups/umcg-bios/tmp04/umcg-fvandijk/projects/GoNLpathogenicVariants/CGD.20161110.txt
print "Processing CGD file..\n";
open(CGD, "< /groups/umcg-bios/tmp04/umcg-fvandijk/projects/GoNLpathogenicVariants/CGD.20161110.txt") || die "Can't open CGD file\n";
my @cgd=<CGD>;
close(CGD);
my %CGDgenes;
for (my $k=1; $k<=$#cgd; $k++){
    my $line = $cgd[$k];
    chomp($line);
    my @array=split("\t", $line);
    my $gene = $array[0];
    my $inheritance = $array[4];
    my $manifestation = $array[7];
    my $intervention = $array[9];
    my $val = "$gene\t$inheritance\t$manifestation\t$intervention";
    $CGDgenes{ $gene } = $val;
}
print "Done processing CGD file.\n\n";


my $countTable = "../output.chr$chrIN.cov$cov.txt";

print "\n\n\nProcessing chromosome: $chrIN\n";

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
            if ($chrIN == $chr) { #If chromosome from input matches the one in GTF file
                if ($info =~ m/gene_id "(ENSG[0-9]{1,})";.+gene_name "(.+)"; gene_sour.+/gs) { #Extract Ensembl gene ID
                    $geneID = $1;
                    $geneName = $2;
                }
                #print "$chr\t$start\t$stop\t$feature\t$geneID\n";
                if ($feature eq "gene") { #If feature is a gene, add it to hash with corresponding start-stop coordinates
                    $geneCount++; #Count number of unique genes
                    if (not exists $genes{ $geneID }) {
                        $genes{ $geneID } = "$start-$stop";
                        #if ($geneID eq "ENSG00000166340") {
                        #    print "\nENSG00000166340: $start-$stop\n\n";
                        #}
                        
                        $geneNames{ $geneID } = "$geneName";
                        $geneIDs{ $geneName } = $geneID;
                        $regionGenes{ "$start-$stop" } = $geneID;
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
                    $regionExons{ "$start-$stop" } = $geneID;
                }
            }
        }
    }
}
close(GTF);
print "#Genes detected: $geneCount\n";
print "Done processing GTF file.\n\n";


print "Processing ClinVar file..\n";
open(CLNVR, "< /groups/umcg-bios/tmp04/umcg-fvandijk/ASE/clinvarPathogenic.txt") || die "Can't open ClinVar file\n";
my @clnvr=<CLNVR>;
close(CLNVR);
my %clinvar;
for (my $h=1; $h<=$#clnvr; $h++){
    my $line = $clnvr[$h];
    chomp($line);
    my @array=split("\t", $line);
    my $geneSymbol = $array[4];
    my $chr = $array[18];
    my $pos = $array[19];
    my $clnsig = $array[6];
    print "$geneSymbol\t$chr\t$pos\t$clnsig\n";
    if (looks_like_number($chr)) { #Check if chromosome is a number
        if ($chrIN == $chr) {
            $clinvar{ $pos } = $geneSymbol; #Push position as key, HUGO geneSymbol as value in hash
        }
    }
}
print "Done processing ClinVar file.\n\n";


print "Processing ClinVar annotation file";
open(CLNVRANN, "< clinvar.1dec2016.all.fix.vcf") || die "Can't open ClinVar annotation file\n";
my @clnvrANN=<CLNVRANN>;
close(CLNVR);
my %clinvarANN;
for (my $p=28; $p<=$#clnvrANN; $p++){
    my $line = $clnvrANN[$p];
    chomp($line);
    my @array=split("\t", $line);
    my $chr = $array[0];
    my $pos = $array[1];
    my $info = $array[7];
    my @annot = split(/\|/, $info);
}

print "Done processing ClinVar annotation file";


print "Processing count table\n";
my %table;
my %tableSamples;
my $TABheaderLine = `head -1 $countTable`; #Extract header from output table
chomp $TABheaderLine;
my @TABheader = split("\t", $TABheaderLine);
pop(@TABheader); #Remove first and last element, since they are not sampleName
shift(@TABheader);
open(INPUT, "< $countTable") || die "Can't open file: $countTable\n"; #Read count table
while (my $line=<INPUT>) {
    chomp($line);
    my @array = split("\t", $line);
    my $lastIdx = $#array;
    my $geneID = $array[0];
    my $geneName = $geneNames{ $geneID }; #Extract HUGO symbol based on Ensembl Gene ID
    if ($line =~ m/Ensem.+/gs) { #Header line
        #;
    }else{ #Gene line, create and print output
        #EnsemblGeneID   AC1C40ACXX-1-18_BD2D5MACXX-6-18
        my $count = $array[$lastIdx];
        if (exists $geneNames{ $geneID }) {
            $table{ $geneName } = $count; #Push counts per geneName in hash
        }
        pop(@array); #Remove first and last element from array
        shift(@array);
        $tableSamples{ $geneName } = join("\t",@array); #Push all sample values in hash, using geneName as key
    }
}
print "Done processing count table\n";


my $VCFheaderLine = `zcat ../../VCFs/chr$chrIN\_ASVCF.vcf.gz | head -500 | grep '^#CHR'`; #Extract header line from VCF
chomp($VCFheaderLine);
my @VCFheader = split("\t", $VCFheaderLine);
my $lastSampleIdx = $#VCFheader;
my %VCFindex;
@VCFindex{@VCFheader} = (0..$#VCFheader);


    my %posits;
    my $geneSymbol = "TPP1";
    my $geneID = $geneIDs{ "$geneSymbol" };
    print "GeneID: $geneID\n";
    my $curGene = $regionGenes{ "6634000-6640692" };
    print "CurGene: $curGene\n";
    my $exonString = $exons{ $regionGenes{ "6634000-6640692" } }; #Extract exon regions by using the key (gene ID)
    my @ar = split(",", $exonString);
    print "\n\n\n";
    
    foreach my $exon (@ar){ #Foreach exon
        my @arr = split("-", $exon);
        my $start = $arr[0];
        my $stop = $arr[1];
        #Foreach exon query the VCF file using tabix
        my $tabixCMD="$tabixPath/tabix ../../VCFs/chr$chrIN\_ASVCF.vcf.gz $chrIN:$exon";
        my $exeCMD=`$tabixCMD`;
        my @vars = split("\n", $exeCMD); #Push tabix results in array
        foreach my $lin (@vars){
            chomp($lin);
            my @array=split("\t", $lin); #Split every line using tab character
            my $lastIdx = $#array;
            my $pos = $array[1];
            $posits{$pos} = 1; #Ensure every SNP position is put in there once!!
            #push(@positions, $pos); #Push all exonic SNPs in array

        }
    }
    
    
    
    
    
    my @positions = keys %posits;
    open(GENEOUTPUT, "> ./chr$chrIN.$geneID.$geneSymbol.gene.txt") || die "Can't open file: chr$chrIN.$geneID.$geneSymbol.gene.txt\n"; #Read count table
    open(HETOUTPUT, "> ./chr$chrIN.$geneID.$geneSymbol.het.txt") || die "Can't open file: chr$chrIN.$geneID.$geneSymbol.het.txt\n"; #Read count table
    open(GTOUTPUT, "> ./chr$chrIN.$geneID.$geneSymbol.gt.txt") || die "Can't open file: chr$chrIN.$geneID.$geneSymbol.gt.txt\n"; #Read count table
    open(ANNOUTPUT, "> ./chr$chrIN.$geneID.$geneSymbol.ann.txt") || die "Can't open file: chr$chrIN.$geneID.$geneSymbol.gt.txt\n"; #Read count table
    
    my $sampleHead = join("\t", @TABheader);
    
    print GENEOUTPUT "SNP\t$sampleHead\n";
    print HETOUTPUT "SNP\t$sampleHead\n";
    print GTOUTPUT "SNP\t$sampleHead\n";
    print ANNOUTPUT "SNP\tANNOTATION\n";
    
    foreach my $pos (@positions){
    #print "Pos: $pos\n";
    my $tabixCMD="$tabixPath/tabix ../../VCFs/chr$chrIN\_ASVCF.vcf.gz $chrIN:$pos-$pos";
    my $exeCMD=`$tabixCMD`;
    
    
    my $tabixANNCMD="$tabixPath/tabix ../../VCFs/chr$chrIN.ASVCF.snpEffAnnotated.vcf.gz $chrIN:$pos-$pos";
    my $exeANNCMD=`$tabixANNCMD`;
    chomp($exeANNCMD);


    if (defined $exeCMD && $exeCMD ne "") {
        my @ann=split("\t", $exeANNCMD);
        my @se = split(";", $ann[7]);
        my @sl = split(/\|/, $se[2]);
        my $anno = $sl[1];
        my $cgdDes = $CGDgenes{ $geneSymbol };
        if (exists $table{ $geneSymbol }) {
            my @hetsDP10;
            my @hets;
            my @homs;
            my $tableCount = $table{ $geneSymbol };
            #if ($tableCount >= $minSamples) {

            print GENEOUTPUT "$chrIN\_$pos";
            print HETOUTPUT "$chrIN\_$pos";
            print GTOUTPUT "$chrIN\_$pos";
            print ANNOUTPUT "$chrIN\_$pos";
            
            my $geneID = $geneIDs{ $geneSymbol };
            #print OUTPUT "$chrIN\t$pos\t$geneSymbol\t$geneID\t$cgdDes\t$tableCount\t";
            chomp($exeCMD);
            my @array = split("\t",$exeCMD); #format per sample  0|0:0,0
            my $ref = $array[3];
            my $alt = $array[4];
            my @tableLine = split("\t",$tableSamples{ $geneSymbol });
            for (my $l=0; $l<=$#tableLine; $l++){
                my $sample = $tableLine[$l];
                #if ($sample ne "NA") {
                    my $search = $TABheader[$l]; #$search is sampleName to search
                    #print GENEOUTPUT "$search\t$chrIN\_$pos\_$ref\_$alt";
                    my $idx = first { $VCFheader[$_] eq $search } 0..$#VCFheader;
                    my $sampleVCF = $array[$idx];
                    #print "\t$sample\t$sampleVCF";
                    my @samVCF = split(":", $sampleVCF);
                    my $gt = $samVCF[0];
                    my @dp = split(",", $samVCF[1]);
                    my $cov = ($dp[0] + $dp[1]);
                    my $d = $dp[0];
                    my $e = $dp[1];
                    my $geneRatioABS="NA";
                    my $hetRatioABS="NA";
                    
                    my $geneRatio = $sample;
                    my $samRatio="NA";
                    if (looks_like_number($sample)) {
                        if ($geneRatio > 0.5) {
                            $geneRatio=(1-$geneRatio);
                        }
                        
                    }
                    if ($cov > 9) {
                        if ($gt eq "0|1") {
                            $samRatio=($dp[0]/$cov); #reference/total
                            #$hetRatio = $sample;
                        }elsif($gt eq "1|0"){
                            $samRatio=($dp[1]/$cov); #Alternative/total
                            #if ($sample ne "NA") {
                            #    $hetRatio = (1-$sample);
                            #}
                        }
                        if (looks_like_number($samRatio)) {
                            if ($samRatio > 0.5) {
                                $samRatio =(1-$samRatio);
                            }
                        }
                    }
                    
                    print GENEOUTPUT "\t$geneRatio";
                    print HETOUTPUT "\t$samRatio";
                    print GTOUTPUT "\t$gt";
            }
            print GENEOUTPUT "\n";
            print HETOUTPUT "\n";
            print GTOUTPUT "\n";
            print ANNOUTPUT "\t$anno\n";
            #print OUTPUT "$het10\t$het\t$hom\t$hets10\t$hets\t$homs\n";
            
            
            #}
        }
    }
    #print "\n";
    }
    
    

print "Done processing chromosome: $chrIN\n";


