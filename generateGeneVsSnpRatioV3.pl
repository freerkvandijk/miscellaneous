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


#open(OUTPUT, "> overview.testOutput.txt") || die "Can't open output file\n";
#1	45797228	MUTYH	ENSG00000132781	MUTYH	AR	Dermatologic; Oncologic	Mutations may also be inolved in susceptibility to other types of malignancies	1232	4	0	1228
#print OUTPUT "CHR\tPOS\tGENENAME\tENSEMBLGENEID\tGENENAME\tINHERITANCE\tCOMM\tCOMM2\tTOTCOUNTS100COVGENE\tHETS10COV\tHETSLESS10COV\tHOM\tHETSAM10COV\tHETSAMLESS10COV\tHOMSAM\n";
`rm ./geneVsSnpRatioV3/chr*.txt`;
for (my $chrIN=1; $chrIN<=22; $chrIN++){
    
    
    my $countTable = "output.chr$chrIN.cov$cov.txt";
    
    
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
                        #if ($geneID eq "ENSG00000126522") {
                        #    print "\nENSG00000126522: $start-$stop\n\n";
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
    #print "$geneSymbol\t$chr\t$pos\n";
    if (looks_like_number($chr)) { #Check if chromosome is a number
        if ($chrIN == $chr) {
            $clinvar{ $pos } = $geneSymbol; #Push position as key, HUGO geneSymbol as value in hash
        }
    }
}
print "Done processing ClinVar file.\n\n";


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


my $VCFheaderLine = `zcat ../VCFs/chr$chrIN\_ASVCF.vcf.gz | head -500 | grep '^#CHR'`; #Extract header line from VCF
chomp($VCFheaderLine);
my @VCFheader = split("\t", $VCFheaderLine);
my $lastSampleIdx = $#VCFheader;
my %VCFindex;
@VCFindex{@VCFheader} = (0..$#VCFheader);

#for (my $j=9; $j <= $lastSampleIdx; $j++){
#    print OUTPUT "\t" . $header[$j];
#}


foreach my $key (sort {$a<=>$b} keys %clinvar){
    my $geneSymbol = $clinvar{ $key };
    
    
    ##if ($geneSymbol eq "ASL") {
    ##my $geneID = $geneIDs{ $geneSymbol };
    ##print "GeneID: $geneID\n";
    ##my $curGene = $regionGenes{ "65540785-65558545" };
    ##print "CurGene: $curGene\n";
    ##my $exonString = $exons{ $regionGenes{ "65540785-65558545" } }; #Extract exon regions by using the key (gene ID)
    ##my @ar = split(",", $exonString);
    #print "\n\n\n";
    #foreach my $el (@ar){
    #    print "$el\t";
    #    my $bla=`$tabixPath/tabix ./VCFs/chr$chrIN\_ASVCF.vcf.gz $chrIN:$el | wc -l`;
    #    chomp($bla);
    #    print "$bla\n";
    #}
    #print "\n\n\n";
    
    
    my $pos = $key;
    my $tabixCMD="$tabixPath/tabix ../VCFs/chr$chrIN\_ASVCF.vcf.gz $chrIN:$pos-$pos";
    my $exeCMD=`$tabixCMD`;
    if (defined $exeCMD && $exeCMD ne "") {
        my $cgdDes = $CGDgenes{ $geneSymbol };
        if (exists $table{ $geneSymbol }) {
            my @hetsDP10;
            my @hets;
            my @homs;
            my $tableCount = $table{ $geneSymbol };
            if ($tableCount >= $minSamples) {

            my $geneID = $geneIDs{ $geneSymbol };
            #print OUTPUT "$chrIN\t$pos\t$geneSymbol\t$geneID\t$cgdDes\t$tableCount\t";
            chomp($exeCMD);
            my @array = split("\t",$exeCMD); #format per sample  0|0:0,0
            my $ref = $array[3];
            my $alt = $array[4];
            if ( -e "./geneVsSnpRatioV3/chr$chrIN.$geneID.$geneSymbol.txt" ) {
                open(GENEOUTPUT, ">> ./geneVsSnpRatioV3/chr$chrIN.$geneID.$geneSymbol.txt") || die "Can't open file: chr$chrIN.$geneID.$geneSymbol.txt\n"; #Read count table
            }else{
                open(GENEOUTPUT, "> ./geneVsSnpRatioV3/chr$chrIN.$geneID.$geneSymbol.txt") || die "Can't open file: chr$chrIN.$geneID.$geneSymbol.txt\n"; #Read count table
                print GENEOUTPUT "SAMPLE\tSNP\tGENERATIO\tHETRATIO\tHETRATIOSAMPLE\n";
            }
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
                    if ($gt eq "1|0" || $gt eq "0|1") { #heterozygous
                        if ($cov > 9){ #homozygous and enough coverage
                            #print "SUFFICIENT COVERAGE Sample: $search\t\tHET id table: $l\t\tHet id VCF: $idx\t\tTAB: $sample\t\tVCF: $gt:$d/$e\n";
                            my $hetRatio = $sample;
                            if ($gt eq "0|1") {
                                $hetRatio=($dp[0]/$cov); #reference/total
                                #$hetRatio = $sample;
                            }elsif($gt eq "1|0"){
                                $hetRatio=($dp[1]/$cov); #Alternative/total
                                #if ($sample ne "NA") {
                                #    $hetRatio = (1-$sample);
                                #}
                            }else{
                                exit("Something horrible happened!!\n");
                            }
                            $hetRatioABS = $hetRatio;
                            if ($hetRatio > 0.5) {
                                $hetRatio=(1-$hetRatio);
                            }
                            
                            print GENEOUTPUT "$search\t$chrIN\_$pos\_$ref\_$alt\tNA\t$hetRatio\t$search\n";
                            push(@hetsDP10, $search);
                        }else{
                            #print GENEOUTPUT "\t$sample\tNA\tNA";
                            push(@hets, $search); #Heterozygous and not enough coverage
                        }
                    }else{ #homozygous
                        #print GENEOUTPUT "\t$sample\tNA\tNA";
                        push(@homs, $search);
                        if ($cov > 9) {
                            if (looks_like_number($sample)) {
                                if ($sample > 0.5) {
                                    $sample=(1-$sample);
                                }
                            }
                            print GENEOUTPUT "$search\t$chrIN\_$pos\_$ref\_$alt\t$sample\tNA\tNA\n";
                        }
                    }
                    #print GENEOUTPUT "\n";
                #}
            }
            close(GENEOUTPUT);
            my $hets10 = join(",", @hetsDP10);
            my $hets = join(",", @hets);
            my $homs = join(",", @homs);
            my $het10 = scalar(@hetsDP10);
            my $het = scalar(@hets);
            my $hom = scalar(@homs);
            #print OUTPUT "$het10\t$het\t$hom\t$hets10\t$hets\t$homs\n";
            
            
            }
        }
        #for (my $s=9; $s <= $lastSampleIdx; $s++){
        #    my @sample=split(":",$array[$s]);
        #    my $gt = $sample[0];
        #    my @dp = split(",", $sample[1]);
        #    my $cov = ($dp[0] + $dp[1]);
        #}
    }
    #print "\n";
    
    
    ##}
    
    
}
print "Done processing chromosome: $chrIN\n";

}
#close(OUTPUT);


#for (my $i=1; $i<=1; $i++){
    #Open file
#    print "Processing chr$i ..\n";
#    open(INPUT, "< output.chr$i.txt") || die "Can't open file: output.chr$i.txt\n";
#
#    print OUTPUT "Sample\tEnsemblGeneID\tRatio\n";
    #while (my $line=<INPUT>) {
    #    chomp($line);
    #    my @array = split("\t", $line);
    #    my $lastIdx = $#array;
    #    my $geneID = $array[0];
    #    my $geneName = $geneNames{ $geneID };
    #    my $headerLine = `head -1 output.chr$i.txt`;
    #    chomp $headerLine;
    #    my @header = split("\t", $headerLine);
    #    if ($line =~ m/Ensem.+/gs) { #Header line
#            #;
#        }else{ #Gene line, create and print output
#            #EnsemblGeneID   AC1C40ACXX-1-18_BD2D5MACXX-6-18
#            my $count = $array[$lastIdx];
#            if ($count > $minSamples-1) { #If count higher than user specified, write to output, else ignore this line
#                for (my $j=1; $j<=$lastIdx-1; $j++){ #Loop over all sample ratios
#                    my $sam=$header[$j];
#                    my $ratio=$array[$j];
#                    print "$sam\t$geneID\t$ratio\n";
#                }
#            }
#        }
#    }
#    close(INPUT);
#    close(OUTPUT);
#}



