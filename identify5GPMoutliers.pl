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
use Statistics::Descriptive;


#module load PerlPlus/5.22.0-foss-2015b-v17.01.1
#module load List-MoreUtils/0.416-foss-2015b-Perl-5.20.2-bare

my $cov = 100;
my $minSamples = 100;
#my $chrIN = 1;
my $tabixPath="/apps/software/tabix/0.2.6-foss-2015b/bin/";


print "Processing CGD file..\n";
open(CGD, "< /groups/umcg-bios/tmp04/umcg-fvandijk/projects/GoNLpathogenicVariants/CGD.20161110.txt") || die "Can't open CGD file\n"; #Open CGD file to extract inheritance information
my @cgd=<CGD>;
close(CGD);
my %CGDgenes;
my %CGDinheritance;
for (my $k=1; $k<=$#cgd; $k++){ #Loop over lines in file
    my $line = $cgd[$k];
    chomp($line);
    my @array=split("\t", $line); #Put in array and extract necessary information
    my $geneName = $array[0];
    my $inheritance = $array[4];
    my $manifestation = $array[7];
    my $intervention = $array[9];
    my $val = "$geneName\t$inheritance\t$manifestation\t$intervention";
    $CGDgenes{ $geneName } = $val; 
    $CGDinheritance{ $geneName } = $inheritance; #geneName and inheritance mode
}
print "Done processing CGD file.\n\n";


print "Processing GTF file..\n";
open(GTF, "< /apps/data/ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf") || die "Can't open GTF file\n";
my %genes;
my %genesChr;
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
                    $genes{ $geneID } = "$chr:$start-$stop"; #geneName chr:start-stop
                    $genesChr{ $geneID } = "$chr"; #geneName chr
                    $geneNames{ $geneID } = "$geneName"; #geneID geneName
                    $geneIDs{ $geneName } = $geneID; #geneName geneID
                    $regionGenes{ "$chr:$start-$stop" } = $geneID; #chr:start-stop geneID
                }
            }
        }
    }
}
close(GTF);
print "#Genes detected: $geneCount\n";
print "Done processing GTF file.\n\n";


print "Processing BIOS count files..\n";
my %BIOScounts;
my @BIOScountFiles = glob("/groups/umcg-bios/tmp04/umcg-fvandijk/ASE/final/output.chr*.cov$cov.txt"); #Retrieve gene ratio files
foreach my $countFile (@BIOScountFiles){
    print "Processing $countFile ..\n";
    open(COUNT, "< $countFile") || die "Can't open file: $countFile\n";
    my @counts=<COUNT>;
    close(COUNT);
    for (my $k=1; $k<=$#counts; $k++){ #Start at 1 to skip headerline 
        my $line = $counts[$k];
        chomp($line);
        my @array = split("\t", $line);
        my $geneID = $array[0]; #geneID
        shift(@array); #remove first element from array
        pop(@array); #remove last element as well, only gene ratios left now
        my $vals = join("\t", @array);
        my $gn = $geneNames { $geneID };
        $BIOScounts{ $geneNames { $geneID } } = $vals; #Retrieve geneName by using geneID and use it as key in hash
    }
}
print "Done processing BIOS count files\n";


print "Processing 5GPM count files..\n";
my %GPMcounts;
my @GPMcountFiles = glob("/groups/umcg-gdio/tmp04/projects/BioSB_5GPM/ASE/outputCounts/output.chr*.cov$cov.txt"); #Retrieve gene ratio files
my %GPMsampleIdx;
my %GPMsampleMatch;
foreach my $count5GPMFile (@GPMcountFiles){
    print "Processing $count5GPMFile ..\n";
    open(GPM, "< $count5GPMFile") || die "Can't open file: $count5GPMFile\n";
    my @counts5GPM=<GPM>;
    close(GPM);
    my $header = $counts5GPM[0];
    chomp($header);
    my @head = split("\t", $header);
    for (my $m=1; $m <= $#head; $m++){
        my $sample = $head[$m];
        chomp $sample;
        $GPMsampleIdx{ $sample } = $m;
        my @sampleToSplit = split("_", $sample); #Split sampleName by underscore, this because the last part of sampleName can be matched with naming scheme of the *.txt.genes files
        #foreach my $samplePart (@sampleToSplit){
        #    if ($samplePart =~ m/^5GPM.+/gs) { #Match on 5GPM tag in front of sampleNumber
        #        #print "$sample\t$samplePart\n";
        #        $GPMsampleMatch{ $samplePart } = $sample;
        #    }
        #}
	if ($sample =~ m/.+\_(5GPM.+)/gs) {
	    chomp($1);
	    $GPMsampleMatch{ $1 } = $sample;
	}
	
    }
    ##
    for (my $l=1; $l<=$#counts5GPM; $l++){ #Start at 1 to skip headerline 
        my $line = $counts5GPM[$l];
        chomp($line);
        my @array = split("\t", $line);
        my $geneID = $array[0]; #geneID
        shift(@array); #remove first element from array
        pop(@array); #remove last element as well, only gene ratios left now
        my $vals = join("\t", @array);
        my $gn = $geneNames { $geneID };
        $GPMcounts{ $geneNames { $geneID } } = $vals; #Retrieve geneName by using geneID and use it as key in hash
    }
}
print "Done processing 5GPM count files\n";


#5GPM counts: /groups/umcg-gdio/tmp04/projects/BioSB_5GPM/ASE/outputCounts/output.chr*.cov*.txt
#BIOS counts: /groups/umcg-bios/tmp04/umcg-fvandijk/ASE/final/output.chr*.cov*.txt
#Genes per sample to assess: /groups/umcg-gdio/prm02/projects/5gpmRna/gavin/1512.outputGavin.txt.genes

#Retrieve genes per sample
#my @sampleGeneFiles = glob("/groups/umcg-gdio/prm02/projects/5gpmRna/gavin/*.outputGavin.txt.genes"); #Loop over genes of interest per sample
my @sampleGeneFiles = glob("/groups/umcg-bios/tmp04/umcg-fvandijk/projects/compare5GPMtoBIOS/gavin/*.outputGavin.txt.genes"); #Loop over genes of interest per sample

#Process all samples and check if genes are present in BIOS gene ratios
print "Processing sample files..\n\n";
foreach my $file (@sampleGeneFiles){
    open(GENES, "< $file") || die "Can't open file: $file\n";
    my @sampleGenes=<GENES>;
    close(GENES);
    #Process every gene and check if it occured
    my $sampleGeneCount=0;
    my $sampleGeneCountExist=0;
    my ($name, $path, $suffix) = fileparse($file, '\.[^\.]*'); #Extract name, path and suffix from file variable
    foreach my $gene (@sampleGenes){ #Loop over genes for this sample
        chomp($gene);
        $sampleGeneCount++;
        #$gene = "TPP1"; #Debugging purpose
        if (exists $BIOScounts{ $gene }) {
            $sampleGeneCountExist++;
            if (exists $GPMcounts{ $gene }) {

            my $mergedVals = $BIOScounts{ $gene } . "\t" . $GPMcounts{ $gene }; #Merge BIOS counts value line with 5GPM counts value line
            my @arrayMergedVals = split("\t", $mergedVals);
            my @ratios;
            foreach my $ratio (@arrayMergedVals){
                if (looks_like_number($ratio)){
                    if ($ratio > 0.5) {
                        $ratio = (1-$ratio);
                    }
                    push(@ratios, $ratio);
                }
            }
            
            my $stat = Statistics::Descriptive::Full->new();
            $stat->add_data(@ratios);
            #print $stat->quantile(0) . "\n"; # => 1       0 => zero quartile (Q0) : minimal value
            #print $stat->quantile(1) . "\n"; # => 3.25    1 => first quartile (Q1) : lower quartile = lowest cut off (25%) of data = 25th percentile
            #print $stat->quantile(2) . "\n"; # => 5.5     2 => second quartile (Q2) : median = it cuts data set in half = 50th percentile
            #print $stat->quantile(3) . "\n"; # => 7.75    3 => third quartile (Q3) : upper quartile = highest cut off (25%) of data, or lowest 75% = 75th percentile
            #print $stat->quantile(4) . "\n"; # => 10      4 => fourth quartile (Q4) : maximal value
            my $firstQ = $stat->quantile(1);
            my $thirdQ = $stat->quantile(3);
            
            #Check if array has at least user specified (100) values, so 100 samples with gene ratio
            my $arrayLength = scalar(@ratios);
            #foreach my $testE (@ratios){
            #    if ($testE <= $firstQ || $testE >= $thirdQ) {
            #        print "Test: $testE\n";
            #        
            #    }
            #}
            if ($arrayLength >= $minSamples) {
                my @GPMarray = split("\t", $GPMcounts{ $gene });
                $name =~ s/.outputGavin.txt//gs;
                my $sampleIdxToRetrieve = "5GPM" . $name;
                if (exists $GPMsampleMatch{ $sampleIdxToRetrieve }) {
		    #if ($sampleIdxToRetrieve eq "5GPM1512") {

                    my $sampleIdx = (($GPMsampleIdx{ $GPMsampleMatch{ $sampleIdxToRetrieve } }) - 1); #Substract 1 because this is based on the header from output counts file, which starts with an EnsemblGeneID column
                    #print "5GPM$name\t$sampleIdxToRetrieve\t$sampleIdx\n";
                    my $sampleRatio = $GPMarray[$sampleIdx];
                    if (looks_like_number($sampleRatio)) {
                        #print "SMRbefore: $gene\t$sampleRatio\n";
			if ($sampleRatio > 0.5){
			#print "SMRafterCorrection: $gene\t$sampleRatio\n";
                            $sampleRatio=(1-$sampleRatio);
                        }
                        if ($sampleRatio <= $firstQ || $sampleRatio >= $thirdQ) { #Possible outlier, write full array of ratios, plus sample away
                            my $CGDstatus = "NA";
                            if (exists $CGDinheritance{ $gene }) {
                                $CGDstatus = $CGDinheritance{ $gene };
                            }
                            $CGDstatus =~ s|/||gs;
                            open(OUTPUT, "> /groups/umcg-bios/tmp04/umcg-fvandijk/projects/compare5GPMtoBIOS/outputs/$sampleIdxToRetrieve.$gene.cov$cov.status.$CGDstatus.txt") || die "Can't open file: /groups/umcg-bios/tmp04/umcg-fvandijk/projects/compare5GPMtoBIOS/outputs/$sampleIdxToRetrieve.$gene.cov$cov.status.$CGDstatus.txt\n";
                            #print "$sampleIdxToRetrieve: $sampleRatio\n"
                            print OUTPUT "SAMPLE\tGENE\tGENERATIO\t5GPMSAMPLE\n";
                            for (my $n=0; $n<=$#ratios; $n++){
                                print OUTPUT "$n\t$gene\t" . $ratios[$n] . "\tNA\n";
                            }
			    #print "SMRprinted: $gene\t$sampleRatio\n\n";
                            print OUTPUT "$sampleIdxToRetrieve\t$gene\t" . $sampleRatio . "\t$sampleIdxToRetrieve\n";
                        }
                        close(OUTPUT);
                    }
		    
		    #}
		    
                }
            }
            }
        }
    }
    print "Sample file: $file\n";
    print "#Genes detected: $sampleGeneCount\n";
    print "#Genes detected in BIOS counts: $sampleGeneCountExist\n";
    print "\n\n";
}


