#!/usr/bin/perl -w

##### STEP 3 : ALIGNMENT OF READS TO REFERENCE GENOME #####
##### Usage: f3 bowtie2_dir reference_genome
##### Required input:
#####   bowtie2_dir - location of bowtie2 installation
#####   reference_genome - the pathname of the reference genome sequence
##### Output:
#####   align/$index_$sample.sam
#####   align/$index_$sample_align.log
#####   *If not already present, reference genome index files*

use strict;
use warnings;

sub f3
{
    #########################################
    ##### DEFAULT VARIABLES FOR BOWTIE2 #####
    #########################################
    my $max_valid_alignments = '1';      # k
    my $num_threads = '4';               # p
    my $max_fragment_length = '11000';   # X
    my $max_reseed_rate = '5';           # R
    #########################################

    # Collect function-specific parameters
    my $bowtie2_dir = $_[0];
    my $reference_genome = $_[1];
    my $sample = $_[2];
    my $index_file = $_[3];

    my @bowtie2_options;
    my $num_options = $#_ + 1;
    if ($num_options > 4)
    {
        @bowtie2_options = @{$_[4]};
    }

    # Ensure index_file exists in the cwd
    unless ( -f $index_file && -r $index_file )
    {   die "ERROR: $index_file does not exist or is unreadable.\n";    }

    # Check for bowtie2
    unless ( -s "$bowtie2_dir/bowtie2" )
    {   die "ERROR: Can't locate bowtie2 in $bowtie2_dir\n";   }

    # Check for the reference genome
    unless (-f $reference_genome && -r $reference_genome)
    {   die "ERROR: $reference_genome does not exist or is unreadable.\n";  }

    system("mkdir -p align/");

    # Remove file extension to get the basename
    my $reference_basename = $reference_genome;
    $reference_basename =~ s/\.f\w*$//;

    # Check for the bowtie2 index files in the same directory as the reference genome,
    # otherwise build the bowtie2 index from the reference genome sequence
    if ( -f "$reference_basename.1.bt2" && -f "$reference_basename.2.bt2" && \
         -f "$reference_basename.3.bt2" && -f "$reference_basename.4.bt2" && \
         -f "$reference_basename.rev.1.bt2" && -f "$reference_basename.rev.2.bt2" )
    {
        print "Reference genome index files present, proceeding to make alignments.\n";
    }
    else
    {
        unless ( -s "$bowtie2_dir/bowtie2-build")
        {   die "ERROR: Could not locate bowtie2-build in $bowtie2_dir\n";   }

        my $cmd = "$bowtie2_dir/bowtie2-build $reference_genome $reference_basename";
        my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
            run( command => $cmd, verbose => 0 );
        if ($success)
        {
            print "Successfully built bowtie index files for $reference_genome\n";
        }
        else {
            print "ERROR: bowtie2_build did not complete successfully:\n";
            die "$error_message\n@$stderr_buf";
        }
    }

    # Print the parameters being given to bowtie2
    print "Running bowtie2 with the following parameters: \n";
    if (@bowtie2_options)
    {
        print join(" ",@bowtie2_options);
        print " --no-sq --no-head\n";
    }
    else
    {
        print "\t--end-to-end --no-mixed --no-discordant --no-sq --no-hd ",
              "-k $max_valid_alignments -X $max_fragment_length -R $max_reseed_rate ",
              "-p $num_threads\n";
    }

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

        # Run bowtie2 using either:
        #  -Default parameters if no options given
        #  -Using custom parameters
        my $cmd;
        if(@bowtie2_options)    # Run with custom paramaters
        {
            $cmd = "$bowtie2_dir/bowtie2 --no-sq --no-head @bowtie2_options ";
            $cmd .= "-x $reference_basename -1 $R1_trimmed -2 $R2_trimmed -S align/$index\_$sample.sam";
        } else                  # Run with default parameters
        {
            $cmd = "$bowtie2_dir/bowtie2 --end-to-end --no-mixed --no-discordant --no-sq --no-head ";
            $cmd .= "-k $max_valid_alignments -X $max_fragment_length -R $max_reseed_rate -p $num_threads ";
            $cmd .= "-x $reference_basename -1 $R1_trimmed -2 $R2_trimmed -S align/$index\_$sample.sam";
        }

        my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
            run( command => $cmd, verbose => 0 );
        if ($success)
        {
            print "Completed alignments for index $index\n";
            # Create a log documenting alignment stats
            my $alignlog = "align/$index\_$sample\_align.log";
            open ALIGNLOG, ">$alignlog" or die "ERROR: Unable to create log file $alignlog\n";
            print ALIGNLOG join " ", @$full_buf;
            close ALIGNLOG or die "ERROR: Unable to close log file $alignlog\n";
        } else
        {
            print "ERROR: bowtie2 did not complete successfully:\n";
            die "$error_message\n@$stderr_buf\n";
        }
    }
}

1;
