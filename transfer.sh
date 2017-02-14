
DIR="/groups/umcg-bios/tmp04/projects/masked_BAMs/transfer/"
TARGETDIR="/groups/umcg-bios/tmp04/projects/masked_BAMs/"

###Create proxy certificate
###On ui.grid.sara.nl type:

#voms-proxy-info --all
#voms-proxy-init -voms bbmri.nl:/bbmri.nl
#voms-proxy-init -voms bbmri.nl:/bbmri.nl/RP3
#voms-proxy-info --all

###Tar all certificates
#tar -cvzf cert.tgz /etc/grid-security/certificates/
###Afterwards retrieve them from grid ui into this directory

###Create file listing to transfer
#uberftp -ls gsiftp://gridftp.grid.sara.nl:2811/pnfs/grid.sara.nl/data/bbmri.nl/RP3/RNASeq/v2.1.3/*.mdup.bam >> allBiosBamsToTransfer.txt
###Afterwards copy file list to this directory
#scp fvandijk@ui.grid.sara.nl:allBiosBamsToTransferFileNames.txt .
###Modify paths
#cp allBiosBamsToTransferFileNames.txt filesToCopy.txt
#perl -pi -e 's|/pnfs/grid.sara.nl/data/bbmri.nl|srm://srm.grid.sara.nl/pnfs/grid.sara.nl/data/bbmri.nl|gs' filesToCopy.txt
#


module load ngs-utils/16.12.1
module load SRM-client/2.16.12-Java-1.7.0_80

perl $EBROOTNGSMINUTILS/copy2grid.pl \
-t g2c \
-f $DIR/filesToCopy2.txt \
-d $TARGETDIR/ \
-p $DIR/x509up_u37197 \
-c $DIR/etc/grid-security/certificates/ \
-l TRACE



