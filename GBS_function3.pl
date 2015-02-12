#!/usr/bin/perl -w

##### STEP 3 : ALIGNMENT OF READS TO REFERENCE GENOME #####
##### Usage: f3 bowtie2_dir reference_genome
##### Required input:
#####   bowtie2_dir - location of bowtie2 installation
#####   reference_genome - the pathname of the reference genome sequence
##### Output:
#####   $output_dir/align/$index_$sample.sam
#####   $output_dir/align/$index_$sample_align.log
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
    my $output_dir = $_[4];

    # Save parameters if provided for bowtie2
    my @bowtie2_options;
    my $num_options = $#_ + 1;
    if ($num_options > 5)
    {
        @bowtie2_options = @{$_[5]};
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

    system("mkdir -p $output_dir/align/");
    system("mkdir -p logs/");

    # Remove file extension to get the basename
    my $reference_basename = $reference_genome;
    $reference_basename =~ s/\.f\w*$//;

    # Check for the bowtie2 index files in the same directory as the reference genome,
    # otherwise build the bowtie2 index from the reference genome sequence
    if ( -s "$reference_basename.1.bt2" && -s "$reference_basename.2.bt2" && \
         -s "$reference_basename.3.bt2" && -s "$reference_basename.4.bt2" && \
         -s "$reference_basename.rev.1.bt2" && -s "$reference_basename.rev.2.bt2" )
    {
        print " Reference genome index files present, proceeding to make alignments.\n";
    }
    else
    {
        unless ( -s "$bowtie2_dir/bowtie2-build")
        {   die "ERROR: Could not locate bowtie2-build in $bowtie2_dir\n";   }

        print " Creating reference index files. Please wait... \n";
        my $cmd = "$bowtie2_dir/bowtie2-build $reference_genome $reference_basename";
        my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
            run( command => $cmd, verbose => 0 );
        if ($success)
        {
            print " Successfully built bowtie index files for $reference_genome\n";
        }
        else {
            print "ERROR: bowtie2_build did not complete successfully:\n";
            die "$error_message\n@$stderr_buf";
        }
    }

    # Print the parameters being given to bowtie2
    print " Running bowtie2 with the following parameters: \n";
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

    # Create summary file
    system("mkdir -p summary_files/");
    my $summary_file = "summary_files/$sample\_align\_summary.txt";
    open SUMMARY, ">$summary_file" or die "ERROR: Could not open $summary_file\n";
    print SUMMARY "Index\tInput Reads\tReads Paired\t% Reads Paired\tOverall Alignment Rate\n";
    close SUMMARY or die "ERROR: Could not close $summary_file\n";

    my @indices = `cat $index_file`;
    my $num_indices = $#indices + 1;
    if ($num_indices == 0){    die "ERROR: $index_file exists but appears to be empty.\n";    }
    my $index_count = 0;
    #print_progress($index_count, $num_indices);

    # BEGIN aligning by index
    foreach my $index (@indices)
    {
        chomp($index);

        my $R1_trimmed = "$output_dir/trim/$index\_$sample\_R1-p.fastq";
        my $R2_trimmed = "$output_dir/trim/$index\_$sample\_R2-p.fastq";

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
            $cmd .= "-x $reference_basename -1 $R1_trimmed -2 $R2_trimmed -S $output_dir/align/$index\_$sample.sam";
        } else                  # Run with default parameters
        {
            $cmd = "$bowtie2_dir/bowtie2 --end-to-end --no-mixed --no-discordant --no-sq --no-head ";
            $cmd .= "-k $max_valid_alignments -X $max_fragment_length -R $max_reseed_rate -p $num_threads ";
            $cmd .= "-x $reference_basename -1 $R1_trimmed -2 $R2_trimmed -S $output_dir/align/$index\_$sample.sam";
        }

        my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
            run( command => $cmd, verbose => 0 );
        if ($success)
        {
            # Create a log documenting alignment stats
            my $alignlog = "logs/$index\_$sample\_bowtie2_output.log";
            open ALIGNLOG, ">$alignlog" or die "ERROR: Unable to create log file $alignlog\n";
            print ALIGNLOG join " ", @$full_buf;
            close ALIGNLOG or die "ERROR: Unable to close log file $alignlog\n";
        } else
        {
            print "ERROR: bowtie2 did not complete successfully:\n";
            die "$error_message\n@$stderr_buf\n";
        }

        # Summary of progress
        $index_count++;
        print_progress($index_count, $num_indices, "Current index: $index");
        summarize_align($index, $sample, \@$stderr_buf);
    }
    print "\n",
    " Processed reads located in:\n  $output_dir/align/ \n",
    " Summary (open in Excel or use command more): $summary_file\n";
}

##############################
##### Summarize step 3
##### Index     # of reads aligned  % of reads aligned
##############################
sub summarize_align
{
    my $index = $_[0];
    my $sample = $_[1];
    my @bowtie2_output = @{$_[2]};

    # All the info is in the first array element, so just store it as a string
    # such that it can be split up into pieces
    my $sample_summary = $bowtie2_output[0];
    # Eliminate indenting to clean it up
    $sample_summary =~ s/^\s+//;
    # Store each line as its own array element
    my @align_info = split(/\n/, $sample_summary);

    ## Can go on a line-by-line basis, which appears to be consistent:
    # 1st line: # of input reads
    $align_info[0] =~ /(\d+) reads;/;
    my $input_reads = $1;
    # 2nd line: # and % of paired reads
    $align_info[1] =~ /(\d+)\s+\(\s?(\d+\.\d+\s?%).+were paired;/;
    my $paired_reads = $1;
    my $percent_paired = $2;
    # Last line: % overall alignment rate
    $align_info[$#align_info] =~ /(\d+\.\d+\s?%)/;
    my $overall_alignment = $1;

    my $summary_file = "summary_files/$sample\_align\_summary.txt";
    open SUMMARY, ">>$summary_file" or die "ERROR: Could not open $summary_file\n";
    print SUMMARY "$index\t$input_reads\t$paired_reads\t$percent_paired\t$overall_alignment\n";
    close SUMMARY or die "ERROR: Could not close $summary_file\n";
}

1;
