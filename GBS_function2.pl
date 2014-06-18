#!/usr/bin/perl -w

##### STEP 2 : TRIM DEMULTIPLEXED READS #####
##### Usage: f2 trimmomatic_path trim_file
##### Required input:
#####   trimmomatic_path : The full pathname to the trimmomatic jar file
#####   trim_file : A list of sequences to filter out the reads (ie. Illumina adaptors)
##### Output:
#####   trim/$index_$sample_R1-s.fastq
#####   trim/$index_$sample_R1-p.fastq
#####   trim/$index_$sample_R2-s.fastq
#####   trim/$index_$sample_R2-p.fastq
#####   trim/$index_$sample_output.log

use strict;
use warnings;

sub f2
{
    #############################################
    ##### DEFAULT VARIABLES FOR TRIMMOMATIC #####
    #############################################
    # Version
        my $version = '0.32';
    # ILLUMINACLIP:
        my $seed_mismatches = '2';
        my $palindrome_clip_threshold = '30';
        my $simple_clip_threshold = '10';
    # SLIDINGWINDOW:
        my $window_size = '4';
        my $required_quality = '15';
    # LEADING:
        my $leading = '3';
    # TRAILING:
        my $trailing = '3';
    # MINLEN:
        my $minlen = '36';
    #############################################

    # Collect function-specific parameters
    my $trimmomatic_path = $_[0];
    my $trim_file = $_[1];
    my $sample = $_[2];
    my $index_file = $_[3];

    # Ensure index_file exists in the cwd
    unless ( -f $index_file && -r $index_file )
    {   die "ERROR: $index_file does not exist or is unreadable.\n";    }

    # Check for trimmomatic
    unless ( -s "$trimmomatic_path" )
    {
        print "WARNING: Can't locate Trimmomatic at $trimmomatic_path.\n",
              "Would you like to install v$version now in the cwd? (yes/no) ";
        chomp (my $usr_input = <STDIN>);
        if ($usr_input =~ /yes/i)
        {
            system("curl -O http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-$version.zip");
            system("unzip Trimmomatic-$version.zip; rm -f Trimmomatic-$version.zip;");

            $trimmomatic_path = "Trimmomatic-$version/trimmomatic-$version.jar";
            my $path = `pwd`;
            if (-e "$trimmomatic_path")
            {
                print "$trimmomatic_path now located in $path\n";
            } else {
                die "ERROR: $trimmomatic_path was not successfully installed in $path\n";
            }
        }
        else {   die "Exiting\n";    }
    }

    # Check for trim_file
    unless ( -r "$trim_file" && -f "$trim_file" )
    {   die "ERROR: $trim_file does not exist or is not readable.\n";   }

    system("mkdir -p trim/");

    my @indices = `cat $index_file`;
    foreach my $index (@indices)
    {
        chomp($index);

        my $R1_reads = "$index\_$sample\_R1-clip.fastq";
        my $R2_reads = "$index\_$sample\_R2.fastq";

        unless ( -f "$R1_reads" && -r "$R1_reads" )
        {   die "ERROR: $R1_reads does not exist or is unreadable.\n"; }
        unless ( -f "$R2_reads" && -r "$R2_reads" )
        {   die "ERROR: $R2_reads does not exist or is unreadable.\n"; }

        # Run Trimmomatic
        my $cmd = "java -classpath $trimmomatic_path org.usadellab.trimmomatic.TrimmomaticPE ";
        $cmd .= "-phred33 ";
        $cmd .= "-trimlog trim/$index.trim.log ";
        $cmd .= "$R1_reads $R2_reads ";
        $cmd .= "trim/$index\_$sample\_R1-p.fastq trim/$index\_$sample\_R1-s.fastq ";
        $cmd .= "trim/$index\_$sample\_R2-p.fastq trim/$index\_$sample\_R2-s.fastq ";
        $cmd .= "ILLUMINACLIP:$trim_file:$seed_mismatches:";
        $cmd .= "$palindrome_clip_threshold:$simple_clip_threshold ";
        $cmd .= "LEADING:$leading ";
        $cmd .= "TRAILING:$trailing ";
        $cmd .= "SLIDINGWINDOW:$window_size:$required_quality ";
        $cmd .= "MINLEN:$minlen ";

        my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
            run( command => $cmd, verbose => 0 );
        if ($success)
        {
            print "Trimmomatic successfully finished trimming $index indexed reads.\n";
            my $trimlog = "trim/$index\_$sample\_output.log";
            open TRIMLOG, ">$trimlog";
            print TRIMLOG join " ", @$full_buf;
            my $summary = `grep 'Input Read Pairs' $trimlog`;
            print "$summary\n";
        } else
        {
            print "ERROR: Trimmomatic did not complete successfully:\n";
            die "$error_message\n@$stderr_buf\n";
        }
    }
}

sub create_summary
{

}

1;
