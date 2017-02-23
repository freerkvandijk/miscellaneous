#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;
use File::Glob ':glob';
use File::Basename;
use Getopt::Long;
use POSIX 'floor';
use GD::Graph::bars;
use GD::Graph::Data;


#Commandline variables
my ($help, $cramDir, $bamDir, $samplesheet, $output);

#### get options
GetOptions(
                "h"             => \$help,
                "cramDir=s"     => \$cramDir, #
                "bamDir=s"      => \$bamDir, #
                "samplesheet:s" => \$samplesheet,
                "outputFilePrefix=s"  => \$output #output file
          );
usage() and exit(1) if $help;
#Obligatory args
usage() and exit(1) unless $cramDir;
usage() and exit(1) unless $bamDir;
usage() and exit(1) unless $output;

my @bamFiles;
my @indices;
my $getColumnIdx;
#Read samplesheet.csv if exists
if (defined $samplesheet) { #If samplesheet specified, use that one, otherwise do glob in BAM directory
    open(FILE, "<$samplesheet") or die("Unable to open samplesheet: $!"); #Read file
    my @file=<FILE>;
    close(FILE);
    my $header = $file[0];
    my @colNames = qw(internalId project sampleName reads1FqGz reads2FqGz sortedBamFile);
    getColumnIdx($header, \@colNames);
    my $iIdIdx = $indices[0];
    my $projectIdx = $indices[1];
    my $sampleNameIdx = $indices[2];
    my $r1Idx = $indices[3];
    my $r2Idx = $indices[4];
    my $numSamples = $#file; #Number of samples is total - headerline, is the same as returning the last Idx of the array (because we have 0 off-set)
    for (my $i=1; $i <= $numSamples; $i++){
        my @array = split(",", $file[$i]);
        my $internalId = $array[$iIdIdx];
        $internalId =~ s/_extr/_extra/g;
        my $sampleName = $array[$sampleNameIdx];
        push(@bamFiles, "$bamDir/$sampleName\_$internalId.bam");
    }
}else{
    #retrieve BAM files
    @bamFiles = glob("$bamDir/*.bam");
    #my @cramFiles = glob("$cramDir/*.cram");
}


#open output file
open (OUTPUTALL, ">$output.all.txt") or die "Cannot open outputfile: $output.all.txt\n";
open (OUTPUTCHECK, ">$output.check.txt") or die "Cannot open outputfile: $output.check.txt\n";
open (OUTPUTREMOVE, ">$output.remove.txt") or die "Cannot open outputfile: $output.remove.txt\n";
open (OUTPUTNOCRAMFILE, ">$output.missingCRAMfile.txt") or die "Cannot open outputfile: $output.missingCRAMfile.txt\n";

my $header = "SampleName\tBAM_file\tCRAM_file\tBAM_file_size\tCRAM_file_size\tCompression_rate";
print OUTPUTALL "$header\n";
print OUTPUTCHECK "$header\n";
print OUTPUTREMOVE "$header\n";
print OUTPUTNOCRAMFILE "$header\n";

my $toPlot;
my %counts;
my $count=0;
foreach my $bamFile (@bamFiles){
    print "Processing BAM file: $bamFile\n";
    #Extract samplename
    my($file, $dir, $ext) = fileparse($bamFile, qr/\.[^.]*/);
    #print $file . "\n";
    my $bamFileSize = `ls -la $bamFile | awk '{print \$5}' FS=" "`;
    chomp $bamFileSize;
    print OUTPUTALL "$file\t$bamFile\t";
    #Check if accompanying CRAM file exists
    my $cramFileOutput = "$cramDir/$file.cram";
    my $cramFileSizeOutput = "N/A";
    my $compressionRatioOutput = "N/A";
    if (-e "$cramDir/$file.cram") {
        print "$file exists\n";
        #Retrieve cram file size
        my $cramFileSize = `ls -la $cramDir/$file.cram | awk '{print \$5}' FS=" "`;
        $cramFileSizeOutput = $cramFileSize;
        chomp $cramFileSizeOutput;
	if ($bamFileSize == 0){
	    $bamFileSize = 1;
	}
        my $compressionRatio = ($cramFileSize/$bamFileSize);
        $compressionRatioOutput = $compressionRatio;
        #chomp $compressionRatioOutput;
        $toPlot = $compressionRatioOutput;
        $toPlot = round(($compressionRatioOutput*10));
        $counts{$toPlot}++;
        $count++;
        if ($compressionRatioOutput <= 0.30 ) { #|| $compressionRatioOutput >= 0.70
            print OUTPUTCHECK "$file\t$bamFile\t$cramFileOutput\t$bamFileSize\t$cramFileSizeOutput\t$compressionRatioOutput\n";
        }else{
            print OUTPUTREMOVE "$file\t$bamFile\t$cramFileOutput\t$bamFileSize\t$cramFileSizeOutput\t$compressionRatioOutput\n";
        }
    }else {
        print "Accompanying CRAM file does not exist, please check again.";
        $cramFileOutput = "N/A";
        chomp $cramFileSizeOutput;
        chomp $compressionRatioOutput;
        print OUTPUTNOCRAMFILE "$file\t$bamFile\t$cramFileOutput\t$bamFileSize\t$cramFileSizeOutput\t$compressionRatioOutput\n";
        #$cramFileSizeOutput = "N/A";
        #$compressionRatioOutput = "N/A";
    }
    chomp $cramFileSizeOutput;
    chomp $compressionRatioOutput;
    print OUTPUTALL "$cramFileOutput\t$bamFileSize\t$cramFileSizeOutput\t$compressionRatioOutput\n";
    print "\n";
}

close(OUTPUTALL);
close(OUTPUTCHECK);
close(OUTPUTREMOVE);
close(OUTPUTNOCRAMFILE);

my @xAxis;
my @yAxis;
foreach my $key (sort {$a<=>$b} keys %counts){
    #print $key . "\t" . $counts{ $key } . "\n";
    my $xAxisFrom = ($key*10);
    my $xAxisTo = ($xAxisFrom+10);
    push(@xAxis, "$xAxisFrom-$xAxisTo");
    push(@yAxis, $counts{ $key });
}

#Plot graph
print "\nCreating graph..\n";
my $data = GD::Graph::Data->new([
    [@xAxis],
    [@yAxis],
]) or die GD::Graph::Data->error;

my $graph = GD::Graph::bars->new(1200, 800);

my $plotTitle = "Compression rates (N=$count)";
$graph->set( 
    x_label         => '% Compression',
    y_label         => '#Samples',
    title           => $plotTitle,
 
    # Sepearte the bars with 4 pixels
    bar_spacing => 4,
    # Show the grid
    #long_ticks  => 1,
    
    y_long_ticks => 1,
    # Show values on top of each bar
    show_values => 1,
    
    #y_max_value     => 7,
    #y_tick_number   => 8,
    #y_label_skip    => 3,
 
    #x_labels_vertical => 1,
 
    #bar_spacing     => 10,
    #shadow_depth    => 4,
    #shadowclr       => 'dred',
 
    #transparent     => 0,
) or die $graph->error;

#$graph->set( 'x_number_format' => \&x_format );
$graph->plot($data) or die $graph->error;

open(my $outGraph, '>', "$output.graph.png") or die "Cannot open '$output.graph.png' for write: $!";
binmode $outGraph;
print $outGraph $graph->gd->png;
close $outGraph;




###SUBS###SUBS###SUBS###

sub x_format
    {
        my $value = shift;
        my $ret;

        if ($value >= 0)
        {
            $ret = sprintf("\$%d", $value + 0.05);
        }
        else
        {
            $ret = sprintf("-\$%d", abs($value) + 0.05);
        }

        return $ret;
    }

sub round {
  my $x = shift;
  $toPlot=floor($x); #To round use (floor($x + 0.5))
  return($toPlot);
}

#Grep column indices and push into array
sub getColumnIdx {
    my $header = shift;
    my ($colnames) = @_;
    my @headerarray = split(quotemeta(","), $header); #Split header on space!
    undef(@indices);
    foreach my $columnName (@$colnames){
        my( $idx )= grep { $headerarray[$_] eq $columnName } 0..$#headerarray;
        push(@indices, $idx);
    }
    return(@indices);
}

#Usage of software
sub usage {
        print <<EOF;
#########################################################################################################

#########################################################################################################

Usage: ./CheckCramFiles.pl
-h\t\t\tThis manual.
-bamDir\t\t\tDirectory containing BAM files to check
-cramDir\t\t\tInput directory containing CRAM files to compare with
-samplesheet\t\t\tOPITONAL: Samplesheet to check. The script only checks the samples
\t\t\t\tspecified in the samplesheet in the cram and bam directory. 
-outputFilePrefix\t\t\tOutput file to write results to

#########################################################################################################
EOF
 
}