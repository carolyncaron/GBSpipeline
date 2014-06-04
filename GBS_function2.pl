#!/usr/bin/perl -w

use strict;
use warnings;
use IPC::Cmd qw[run];

##### STEP 2 : TRIM DEMULTIPLEXED READS #####
##### Usage: f2 trimmomatic_path trim_file
##### Required input:
#####   trimmomatic_path - The full pathname to the trimmomatic jar file
#####   trim_file - A list of sequences to filter out the reads (ie. Illumina adaptors)
#####
#####

sub f2
{
    ##### DEFAULT VARIABLES FOR TRIMMOMATIC #####
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

    # Collect function-specific parameters
    my $trimmomatic_path = $_[0];
    my $trim_file = $_[1];
    my $sample = $_[2];
    my $index_file = $_[3];

    # Check for trimmomatic
    unless ( -e "$trimmomatic_path" )
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
    {
        die "ERROR: $trim_file does not exist or is not readable.\n";
    }

    system("mkdir -p trim/");

    # Run trimmomatic
    my @indices = `cat $index_file`;
    foreach my $index (@indices)
    {
        chomp($index);
        my $cmd = "java -classpath $trimmomatic_path org.usadellab.trimmomatic.TrimmomaticPE ";
        $cmd .= "-phred33 ";
        $cmd .= "-trimlog trim/$index.trim.log ";
        $cmd .= "$index\_$sample\_R1-clip.fastq $index\_$sample\_R2.fastq ";
        $cmd .= "trim/$index\_$sample\_R1-p.fastq trim/$index\_$sample\_R1-s.fastq ";
        $cmd .= "trim/$index\_$sample\_R2-p.fastq trim/$index\_$sample\_R2-s.fastq ";
        $cmd .= "ILLUMINACLIP:$trim_file:$seed_mismatches:$palindrome_clip_threshold:$simple_clip_threshold ";
        $cmd .= "LEADING:$leading ";
        $cmd .= "TRAILING:$trailing ";
        $cmd .= "SLIDINGWINDOW:$window_size:$required_quality ";
        $cmd .= "MINLEN:$minlen ";

        my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
            run( command => $cmd, verbose => 0 );
        if ($success)
        {
            print "Trimmomatic successfully finished trimming $index indexed reads.\n";
            open TRIMLOG, ">trim/$index\_$sample\_output.log";
            print TRIMLOG join " ", @$full_buf;
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
