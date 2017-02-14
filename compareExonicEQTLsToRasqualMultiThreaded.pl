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
my $num_of_threads = 20;

##Input dir
my $dir="./exon_level_eQTLs_independent_effects_chunks_1e-5/";

#Shared hash to store all results in
my %results;
share(%results);

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



#Loop overal original file to output all results
print "\n\nProcessing results and generating output file..\n\n";
open(OUTPUT, "> comparisonExonicEQTLsToRasqual.txt") or die("Unable to open comparison file");

open(E, "< exon_level_eQTLs_independent_effects.10e-5.txt") or die("Unable to open eQTL file"); #Read file
while (my $lin=<E>) {
    chomp $lin;
    if (exists $results{ $lin }) {
        print OUTPUT $results{ $lin };
    }
}
close(E);
close(OUTPUT);
print "\nDone processing all results!\n\n";


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

    my $id = threads->tid();
    my $n = ($id-1);
    my $chunkNum = sprintf ("%02d",$n);
    print "Processing file: $dir/exon_level_eQTLs_independent_effects.chunk$chunkNum..\n";
    open(EQTL, "< $dir/exon_level_eQTLs_independent_effects.chunk$chunkNum") or die("Unable to open eQTL file"); #Read file
    #open(EQTL, "< exon_level_eQTLs_Top_effects.exonicSNPsOnly_10E-8.txt") or die("Unable to open eQTL file"); #Read file

    my $conc=0;
    my $disconc=0;
    my $total=0;
    my $concThres=0;
    my $disconcThres=0;
    my $totalThres=0;
    my $counter=0;
    while (my $line = <EQTL>) {
        chomp $line;
        my $key = $line;
        my $outputLine;
        my @array = split("\t", $line);
        #PValue  SNPName SNPChr  SNPChrPos       ExonId  ExonChr ExonCenterChrPos        SNPType AlleleAssessed  OverallZScore  OverallZScore   GeneId  FDR     GeneName
        #3.27167E-310    rs72885163      1       46807117        chr_1_46806859_46809966 1       46808412        C/A     A       48.4722542      48.4722542      ENSG00000117481 0       NSUN4
        my $rsID = $array[1];
        my $chr = $array[2];
        my $pos = $array[3];
        my $loc = "$chr\t$pos";
        my $exonID = $array[4];
        $exonID =~ s/chr_//gs;
        #Split exonID further because sometimes end positions of exons don't match between ensemebl v71 and v75
        my @exonArr = split("_", $exonID);
        my $exonIDToSearch = $exonArr[0] . "_" . $exonArr[1] . "_";
        my $snpType = $array[7];
        my $refAll;
        my $altAll;
        if ($snpType =~ m/(.+)\/(.+)/gs) {
            $refAll = $1;
            $altAll = $2;
        }
        my $assessedAll = $array[8];
        my $zScore = $array[9];
        my $geneIDstring = $array[11];
        $geneIDstring =~ s/"//gs;
        my @geneIDs = split(",", $geneIDstring);
        #print "$rsID\t$chr\t$pos\t$refAll\t$altAll\t$assessedAll\t$Zscore\t$geneID\n";
        #Grep Ensembl geneID in rasqual result directory to retrieve correct file.
        #Order SNPs in file by Rsquare value to retrieve rank of SNP in gene
        #Also extract alleles and effect size to compare to eQTL results
        if ($#geneIDs == 0) {
        $counter++;
        #foreach my $geneID (@geneIDs){
            my $grep = `grep -n -P "\t$exonIDToSearch" /apps/data/UMCG/Ensembl.GRCh37.75-Exon_And_GeneList/meta-exonlist_sorted.chr$chr.txt | awk '{print \$1}' FS=":"`;
            chomp $grep;
            #print "$exonID\n";
            #print "blaat.$grep.blaat\n";
            if (defined $grep && $grep ne "") {
                my $chunkNum = ($grep-1);
                my $dig = sprintf ('%05s', $chunkNum);
                #print "$chunkNum\t$dig\n";
                #results/Rasqual/gene/chr1/genelist_sorted.chr1.chunk0859.rasqual.output.txt 
                #1 26856252 26901521 + ENSG00000117676 RPS6KA1 ensembl_havana protein_coding 0.464866 0.535134 11031 0 1 26856252 26901521 860 Output:   rs2736831       1       26357656        C       A       0.712687        0.068426        1.000000        -0.0725920838   0.1147856803   0.501321        0.004090        0.503934        27.188050       0       15      151     6       6       26758773        -7761.805364    0       0.989698        0.990427
                my $chunkFile = "results/Rasqual/metaExon/chr$chr/meta-exonlist_sorted.chr$chr.chunk$dig.rasqual.output.txt";
                open(CHUNK, "< $chunkFile") or die("Unable to open chunk file: $chunkFile"); #Read file
                my %rank;
                my $lCount=0;
                my $updatedZscore = $zScore;
                while (my $line = <CHUNK>){
                    chomp $line;
                    my @arr=split("\t", $line);
                    my $crsID = $arr[1];
                    my $cchr = $arr[2];
                    my $cpos = $arr[3];
                    my $cloc = "$cchr\t$cpos";
                    my $crefAll = $arr[4];
                    my $caltAll = $arr[5];
                    my $chiSq = $arr[10];
                    my $effSize = $arr[11];
                    my $NoFSNPs = $arr[16];
                    my $NoSNPs = $arr[17];
                    my $SqCor = $arr[24];
                    my $AF = $arr[6];
                    #If locations match, do comparison and extract effectSize
                    if ($loc eq $cloc) {
                                        my $effSize = $arr[11];
                        #check if alleles match or need to be swapped
                        #print "$cloc\t$crefAll\t$caltAll\t$chiSq\t$effSize\t";
                        my $correctedEffSize = $effSize;
                        $outputLine .= "$loc\t$refAll\t$altAll\t";
                        if ($refAll eq $crefAll && $altAll eq $caltAll) { #No swaps needed
                            
                        }elsif($refAll eq $caltAll && $altAll eq $crefAll){ #Alleles swapped, so swap effect direction
                            #$correctedEffSize = (1-$effSize);
                            $altAll = $refAll;
                        }
                        #Check which allele was assessed in eQTL data, if alternative allele was assessed swap Zscore of eQTL
                        if ($altAll ne $assessedAll) { #Alt was assessed, swap Zscore direction
                            if ($zScore =~ m/^-(.+)/gs){
                                $updatedZscore = $1;
                            }elsif ($zScore =~ m/^[0-9]{1,}.+/gs){
                                $updatedZscore = "-$zScore";
                            }
                        }
                        my $zScoreDirection;
                        if ($updatedZscore =~ m/^-(.+)/gs) {
                            $zScoreDirection = "NEG";
                        }else{
                            $zScoreDirection = "POS";
                        }
                        my $effSizeDirection = "NEG";;
                        
                        #print "\n\nBLAAT: $correctedEffSize\n\n";
                        
                        if ($correctedEffSize < 0.5 ) {
                            $effSizeDirection = "NEG";
                        }else{
                            $effSizeDirection = "POS";
                        }
                        
                        if ($zScoreDirection eq $effSizeDirection) {
                            $conc++;
                            if (abs(log($chiSq)) >= 3) {
                                $concThres++;
                                $totalThres++;
                            }
                        }else{
                            $disconc++;
                            if (abs(log($chiSq)) >= 3) {
                                $disconcThres++;
                                $totalThres++;
                            }
                        }
                        $total++;
                        
                        $outputLine .= "$assessedAll\t$zScore\t$cloc\t$crefAll\t$caltAll\t$AF\t$chiSq\t$effSize\t\t$updatedZscore\t$zScoreDirection\t$correctedEffSize\t$effSizeDirection\t$NoFSNPs\t$NoSNPs\t$SqCor\t";
                    }
                    #push rsID and chiSq in hash to sort later
                    $rank{ $cloc } = $chiSq;
                    $lCount++;  
                }
                close(CHUNK);
                
                my @clocs = reverse sort { $rank{$a} <=> $rank{$b} } keys(%rank); #Reverse sort the values, making highest chiSq first
                my ($index) = grep { $clocs[$_] ~~ $loc } 0 .. $#clocs;
                if (defined $index && $index ne "") {
                
                my $rank = ($index+1);
                my $normRank = 1-($rank/$lCount);
                my $rounded = sprintf "%.2f", $normRank;
                $outputLine .= "$rank\t$normRank\t$rounded\n";
                $results{ $line } = $outputLine;
                }
            }else{
                #`echo $line >> tmp.txt`;
            }
        }
        #}
    }
    
    # Inform us that the thread is done and exit the thread.
    print "Thread $id done processing file: $dir/exon_level_eQTLs_independent_effects.chunk$chunkNum!\n";
    threads->exit();
}
