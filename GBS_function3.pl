#!/usr/bin/perl -w

use strict;
use warnings;
use IPC::Cmd qw[run];

##### STEP 3 : ALIGNMENT OF READS TO REFERENCE GENOME #####
##### Usage: f3 bowtie2_path reference_genome
##### Required input:
#####   bowtie2_path - location of bowtie2 installation
#####   reference_genome - basename of the reference genome index file
#####
#####

sub f3
{
    ##### DEFAULT VARIABLES FOR BOWTIE2 #####
    my $max_valid_alignments = '1';     # k
    my $num_threads = '24';             # p
    my $max_fragment_length = '11000';  # X
    my $max_reseed_rate = '5';          # R

    my $version = '2.2.3';

    # Collect function-specific parameters
    my $bowtie2_path = $_[0];
    my $reference_genome = $_[1];
    my $sample = $_[2];
    my $index_file = $_[3];

    # Ensure index_file exists in the cwd
    unless ( -f $index_file && -r $index_file )
    {
        die "ERROR: $index_file does not exist or is unreadable.\n";
    }

    # Check for bowtie2
    unless ( -e "$bowtie2_path" )
    {
        print "WARNING: Can't locate bowtie2 at $bowtie2_path.\n",
              "Would you like to install v$version now in the cwd? (yes/no) ";
        chomp (my $usr_input = <STDIN>);
        if ($usr_input =~ /yes/i)
        {
            system("curl -O http://sourceforge.net/projects/bowtie-bio/files/latest/download?source=files");
            system("unzip bowtie2-*.zip; rm -f Trimmomatic-$version.zip;");

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

    system("mkdir -p align/");

    my @indices = `cat $index_file`;
    foreach my $index (@indices)
    {
        chomp($index);

        my $R1_trimmed = "trim/$index\_$sample\_R1-p.fastq";
        my $R2_trimmed = "trim/$index\_$sample\_R2-p.fastq";

        unless ( -f $R1_trimmed && -r $R1_trimmed )
        {   die "ERROR: $R1_trimmed does not exist or is unreadable.\n"; }
        unless ( -f $R2_trimmed && -r $R2_trimmed )
        {   die "ERROR: $R2_trimmed does not exist or is unreadable.\n"; }

        # Run bowtie2
        my $cmd = "bowtie2 --end-to-end --no-mixed --no-discordant --no-sq --no-head ";
        $cmd .= "-k $max_valid_alignments -X $max_fragment_length -R $max_reseed_rate -p $num_threads ";
        $cmd .= "-x $reference_genome -1 $R1_trimmed -2 $R2_trimmed -S align/$index_$sample.sam";

        my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
            run( command => $cmd, verbose => 0 );
        if ($success)
        {
            print ".\n";
            my $alignlog = "align/$index\_$sample\_output.log";
            open ALIGNLOG, ">$alignlog";
            print ALIGNLOG join " ", @$full_buf;
            #my $summary = `grep 'Input Read Pairs' $alignlog`;
            #print "$summary\n";
        } else
        {
            print "ERROR: bowtie2 did not complete successfully:\n";
            die "$error_message\n@$stderr_buf\n";
        }
    }
}

1;
