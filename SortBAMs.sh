





module load picard/2.7.2-foss-2015b-Java-1.8.0_74


java -jar -Xmx4g -XX:ParallelGCThreads=2 ${EBROOTPICARD}/picard.jar \
ReorderSam \
INPUT=./BAMsFromGrid/AC43GJACXX-6-7.mdup.bam \
OUTPUT=test.output.bam \
REFERENCE=/apps/data/ftp.broadinstitute.org/bundle/2.8/b37/human_g1k_v37.fasta \
CREATE_INDEX=TRUE




