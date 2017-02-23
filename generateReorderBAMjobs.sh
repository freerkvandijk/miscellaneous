
JOBDIR="/groups/umcg-bios/tmp04/projects/masked_BAMs/reorderJobs/"

for BAM in $( ls /groups/umcg-bios/tmp04/projects/masked_BAMs/BAMsFromGrid/*.bam )
do
    
PREFIX=$(basename $BAM .bam)
SAMPLE=$(basename $BAM .mdup.bam)
OUTPUTDIR="/groups/umcg-bios/tmp04/projects/masked_BAMs/BAMsSorted"
    

echo "#!/bin/bash
#SBATCH --job-name=reorderBAM.$SAMPLE
#SBATCH --output=reorderBAM.$SAMPLE.out
#SBATCH --error=reorderBAM.$SAMPLE.err
#SBATCH --time=01:59:00
#SBATCH --cpus-per-task 2
#SBATCH --mem 4gb
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=30L
#SBATCH --qos=leftover

set -e
set -u

ENVIRONMENT_DIR='.'

module load tabix/0.2.6-foss-2015b
module load picard/2.7.2-foss-2015b-Java-1.8.0_74
module list

echo \"## \"\$(date)\" Start \$0\"

mkdir -p $OUTPUTDIR


if java -jar -Xmx4g -XX:ParallelGCThreads=2 -Djava.io.tmpdir=\$TMPDIR \${EBROOTPICARD}/picard.jar \\
ReorderSam \\
INPUT=$BAM \\
OUTPUT=$OUTPUTDIR/$PREFIX.sorted.bam \\
REFERENCE=/apps/data/ftp.broadinstitute.org/bundle/2.8/b37/human_g1k_v37.fasta \\
CREATE_INDEX=TRUE
then
 echo \"returncode: \$?\";
 cd $OUTPUTDIR
 bname=\$(basename $OUTPUTDIR/$PREFIX.sorted.bam)
 md5sum \${bname} > \${bname}.md5
 bname=\$(basename $OUTPUTDIR/$PREFIX.sorted.bai)
 md5sum \${bname} > \${bname}.md5
 cd -
 echo \"succes moving files\";
fi

touch reorderBAM.$SAMPLE.sh.finished

echo \"## \"\$(date)\" ##  \$0 Done \"

">$JOBDIR/reorderBAM.$SAMPLE.sh

done
