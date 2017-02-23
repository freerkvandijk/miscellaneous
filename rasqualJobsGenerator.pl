#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;
use File::Glob ':glob';
use File::Basename;
use Getopt::Long;
use Term::ProgressBar;
use lib "/home/umcg-fvandijk/perl_modules/threads-2.09/lib/";
use threads;
use lib "/home/umcg-fvandijk/perl_modules/threads-shared-1.52/lib/threads/";
use threads::shared;

#module load Term-ProgressBar/2.17-foss-2015b


#Read param list
print "Reading param csv..\n";
#open(PARAM, "< test.csv") or die("Unable to open param file"); #Read file
open(PARAM, "< ./ASE/chromosomeMetaExonChunks.csv") or die("Unable to open param file"); #Read file
my @params=<PARAM>;
close(PARAM);
print "Done reading param csv!\n";

my $lastIdx = $#params;
#Loop over array and generate jobs

print "\nGenerating jobs..\n\n";
my $progress_bar = Term::ProgressBar->new($lastIdx);

for (my $i=1; $i<=$lastIdx; $i++){
my $line=$params[$i];
chomp $line;
my @array=split(",", $line);
my $CHR=$array[0];
my $featureChunkFile=$array[1];
my @naming=split(/\./, $featureChunkFile); #split on dot
my $jobName = $naming[1] . "." . $naming[2] . ".runRasqual";

#Open output file
#`mkdir -p ./testJobs/chr$CHR/`;
#open(OUTPUT, "> ./testJobs/chr$CHR/$jobName.sh") or die("Unable to open param file"); #Read file
`mkdir -p ./jobs/runRasqual/chr$CHR/`;
open(OUTPUT, "> ./jobs/runRasqual/chr$CHR/$jobName.sh") or die("Unable to open param file"); #Read file

my $project="GoNL";
my $stage="module load";
my $checkStage="module list";
my $projectDir="/groups/umcg-bios/tmp04/projects/ASE_GoNL/rasqual/GoNL_WGS_meta-exons_allSNPs/results/";
my $binDir="$projectDir/bins/";
my $ASVCF="$projectDir/bins/chr$CHR\_ASVCF.vcf.gz";
my $RASQUALDIR="/groups/umcg-wijmenga/tmp04/tools/rasqual/";
my $cisWindow="1000000";
my $kfilebinMetaExon="$binDir/KmetaExon.chr$CHR.bin";
my $yfilebinMetaExon="$binDir/YmetaExon.chr$CHR.bin";
my $yfiletxtMetaExon="$binDir/YmetaExon.chr$CHR.txt";
my $featureType="metaExon";
my $rasqualOutDir="$projectDir/Rasqual";
my $regionsFile="/apps/data/UMCG/Ensembl.GRCh37.75-Exon_And_GeneList/genomeWideRegion.chr$CHR.bed";
my $GSLVersion="2.1-foss-2015b";
my $tabixVersion="0.2.6-goolf-1.7.20";
my $minCoveragePerFeature="35";
my $insertSize="250";
my $featureChunkDir="/apps/data/UMCG/Ensembl.GRCh37.75-Exon_And_GeneList/";
my $rasqualFeatureChunkOutput="$rasqualOutDir/$featureType/chr$CHR/$featureChunkFile.rasqual.output.txt";
my $rasqualFeatureChunkPermutationOutput="$rasqualOutDir/$featureType/chr$CHR/$featureChunkFile.rasqual.output.permutation.txt";

print OUTPUT <<EOF;
#!/bin/bash
#SBATCH --job-name=$jobName
#SBATCH --output=$jobName.out
#SBATCH --error=$jobName.err
#SBATCH --qos=leftover
#SBATCH --time=05:59:00
#SBATCH --cpus-per-task 1
#SBATCH --mem 2gb
#SBATCH --nodes 1
  
echo \"## \"\$(date)\" ##  \$0 Start \"

getFile $featureChunkFile
getFile $ASVCF
getFile $regionsFile


$stage GSL/$GSLVersion
$stage tabix/$tabixVersion
$checkStage

touch $jobName.sh.started

mkdir -p $rasqualOutDir/$featureType/chr$CHR/

echo LAST START >> $rasqualFeatureChunkOutput

kfilebin=$kfilebinMetaExon
yfilebin=$yfilebinMetaExon
yfiletxt=$yfiletxtMetaExon
featureDir=$featureChunkDir/meta-exonlistChunksPerFeature/chr$CHR/


rm -f $rasqualFeatureChunkOutput
rm -f $rasqualFeatureChunkPermutationOutput

#Copy feature file, bgzip and tabix it (this to reduce number of files on shared storage space)
cp \$featureDir/$featureChunkFile \$TMPDIR/$featureChunkFile
bgzip -c \$TMPDIR/$featureChunkFile > \$TMPDIR/$featureChunkFile.gz
tabix -s 1 -b 2 -e 3 \$TMPDIR/$featureChunkFile.gz


#Run analysis
window=\$(($cisWindow/2)) # 1Mb
samples_num=\$(awk -F'\\t' '{print NF-1; exit}' \$yfiletxt)
cutoff=\$((\$samples_num * $minCoveragePerFeature))
Top=\$(tabix $ASVCF $CHR: | tail -n 1 | cut -f 2)
while read region;do
while read line;do
		echo \"Analyzing feature: \$line\";
        #INIT#########################
        array=(\$line)
        id=\$(echo \$line| cut -f1-5,7-9,13,16)
        chr=\"\${array[0]}\"
        start=\"\${array[1]}\"
        end=\"\${array[2]}\"
        featureStarts=\"\${array[13]}\"
        featureEnds=\"\${array[14]}\"
        seq_len=\"\${array[12]}\"
        line_number=\"\${array[15]}\"
        Cutoff_for_this=\$((\$cutoff * (\$seq_len/$insertSize)))
        Coverage=\$(sed "\${line_number}q;d" \$yfiletxt | awk '{ for(i=2; i<=NF;i++) j+=\$i; print \$j; j=0 }') 
        if (( \$start < \$window )); then L=1
                else L=\$((\$start - \$window)); fi
        if (( (\$start + \$window) > \$Top )); then R=\$Top
                else R=\$((\$start + \$window)); fi
        Totalsnps="\$(tabix $ASVCF \$chr:\$L-\$R | wc -l)"
        echo \"Search space: chr\$chr:\$L-\$R\"
        echo \"TotalSNPs: \$Totalsnps\"
        echo \"\"
        #PREFILTERS###########################################
#       if (( \$Coverage < \$Cutoff_for_this)); then continue; fi
        ######################################################
        tabix $ASVCF \$chr:\$L-\$R | $RASQUALDIR/bin/rasqual --force -y \$yfilebin -k \$kfilebin -n \$samples_num -j \$line_number -l \$Totalsnps -m \$Totalsnps -s \$featureStarts -e \$featureEnds -f \"\$id Output:\" --n-threads 1 >> $rasqualFeatureChunkOutput
        tabix $ASVCF \$chr:\$L-\$R | $RASQUALDIR/bin/rasqual --force -y \$yfilebin -k \$kfilebin -n \$samples_num -j \$line_number -l \$Totalsnps -m \$Totalsnps -s \$featureStarts -e \$featureEnds -f \"\$id Output:\" --n-threads 1 -r >> $rasqualFeatureChunkPermutationOutput
done < <(tabix \$TMPDIR/$featureChunkFile.gz \"\$region\" )
done < <(awk 'F\"\\t\" \$(\$1 == $CHR) {printf (\"%s:%s-%s\\n\", \$1, \$2, \$3)}' $regionsFile)

#Putfile the results
if [ -f \"$rasqualFeatureChunkOutput\" ];
then
 echo \"returncode: \$?\"; 
 putFile $rasqualFeatureChunkOutput
 putFile $rasqualFeatureChunkPermutationOutput
 echo \"succes moving files\";
else
 echo \"returncode: \$?\";
 echo \"fail\";
 exit 1;
fi


echo \"## \"\$(date)" ##  \$0 Done \"

touch $jobName.sh.finished


EOF

close(OUTPUT);

$progress_bar->update($i);

}

print "\n\nDone generating jobs..\n";
