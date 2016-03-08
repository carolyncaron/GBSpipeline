#!/usr/bin/perl -w

##### STEP 3 : ALIGNMENT OF READS TO REFERENCE GENOME #####
##### Usage: f3 bowtie2_dir reference_genome
##### Required input:
#####   bowtie2_dir - location of bowtie2 installation
#####   reference_genome - the pathname of the reference genome sequence
##### Output:
#####   $output_dir/align/$sample_$population.sam
#####   $output_dir/align/$sample_$population_align.log
#####   *If not already present, reference genome index files*

use strict;
use warnings;

sub f3
{
    #########################################
    ##### DEFAULT VARIABLES FOR BOWTIE2 #####
    #########################################
    my %align_options = (
    # p
        'ALIGN_THREADS' => '1',
    # k
    ### Default is set to 3 to allow a single read to map up to 3 times. Note that setting
    ### this value higher will greatly increase the amount of time bowtie2 needs to run.
        'MAX_VALID_ALIGNMENTS' => '3',
    # X
        'MAX_FRAGMENT_LENGTH' => '11000',
    # R
        'MAX_RESEED_RATE' => '5',
    );
    #########################################

    # Collect function-specific parameters
    my $bowtie2_dir = $_[0];
    my $reference_genome = $_[1];
    my $population = $_[2];
    my $index_file = $_[3];
    my $output_dir = $_[4];

    # Collect bowtie2-specific options
    if ($_[5])
    {
        my %options = %{$_[5]};
        while ( my ($key, $value) = each %options )
        {
            if ($options{ $key })
            {
                $align_options{ $key } = $value;
            }
        }
    }
    print Dumper(\%align_options);

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

    # Begin recording progress since building reference index files can take a while
    my @samples = `cut -f1 $index_file | sort -u`;
    my $num_samples = $#samples + 1;
    if ($num_samples == 0){    die "ERROR: $index_file exists but appears to be empty.\n";    }
    my $sample_count = 0;

    # Remove directory and file extension to get the basename
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

    # Create summary file
    system("mkdir -p summary_files/");
    my $summary_file = "summary_files/$population\_align\_summary.txt";
    open SUMMARY, ">$summary_file" or die "ERROR: Could not open $summary_file\n";
    print SUMMARY "Sample\tInput Reads\tUnique Reads\t% Unique\tOverall Alignment Rate\n";
    close SUMMARY or die "ERROR: Could not close $summary_file\n";

    # BEGIN aligning by sample
    foreach my $sample (@samples)
    {
        chomp($sample);
        $sample =~ s/ //g;
        $sample_count++;
        print_progress($sample_count, $num_samples, "Current sample: $sample       ");

        my $R1_trimmed = "$output_dir/trim/$sample\_$population\_R1-p.fastq";
        my $R2_trimmed = "$output_dir/trim/$sample\_$population\_R2-p.fastq";

        unless ( -f $R1_trimmed && -r $R1_trimmed )
        {   die "ERROR: $R1_trimmed does not exist or is unreadable.\n"; }
        unless ( -f $R2_trimmed && -r $R2_trimmed )
        {   die "ERROR: $R2_trimmed does not exist or is unreadable.\n"; }

        my $cmd = "$bowtie2_dir/bowtie2 --end-to-end --no-mixed --no-discordant --no-sq --no-head ";
        $cmd .= "-k $align_options{'MAX_VALID_ALIGNMENTS'} -X $align_options{'MAX_FRAGMENT_LENGTH'} ";
        $cmd .= "-R $align_options{'MAX_RESEED_RATE'} -p $align_options{'ALIGN_THREADS'} ";
        $cmd .= "-x $reference_basename -1 $R1_trimmed -2 $R2_trimmed -S $output_dir/align/$sample\_$population.sam";

        ############## SEQUENTIAL RUN ###############
        # Run bowtie2 for each sample sequentially. If desired to run them concurrently,
        #   comment out the following and uncomment the section below
        my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
            run( command => $cmd, verbose => 0 );
        if ($success)
        {
            # Create a log documenting alignment stats
            my $alignlog = "logs/$sample\_$population\_bowtie2_output.log";
            open ALIGNLOG, ">$alignlog" or die "ERROR: Unable to create log file $alignlog\n";
            print ALIGNLOG join " ", @$full_buf;
            close ALIGNLOG or die "ERROR: Unable to close log file $alignlog\n";
        } else
        {
            print "ERROR: bowtie2 did not complete successfully:\n";
            die "$error_message\n@$stderr_buf\n";
        }

        # Summary of progress
        summarize_align($sample, $population, \@$stderr_buf);
        print_progress($sample_count, $num_samples, "                             ");
        ##############################################

        ############## CONCURRENT RUN ################
        # The following code allows the individual bowtie2 commands to be run in the background,
        # effectively aligning reads from all samples concurrently.
        # CAUTION: Ensure you have at least as many threads as samples!
        #my $alignlog = "logs/$sample\_$population\_bowtie2_output.log";
        #system("$cmd 2>$alignlog &");
        ##############################################
    }

    print "\n",
    " Processed reads will be located in:  $output_dir/align/ \n",
    " Summary file:  $summary_file\n";
}

##############################
##### Summarize step 3
##### sample     # of input reads   # of unique reads   % unique reads    % overall alignment
##############################
sub summarize_align
{
    my $sample = $_[0];
    my $population = $_[1];
    my @bowtie2_output = @{$_[2]};
    #print join " ", @bowtie2_output;

    # All the info is in the first array element, so just store it as a string
    # such that it can be split up into pieces
    my $sample_summary = $bowtie2_output[0];
    # Eliminate indenting to clean it up
    $sample_summary =~ s/^\s+//;
    # Store each line as its own array element
    my @align_info = split(/\n/, $sample_summary);

    ## Can go on a line-by-line basis, which appears to be consistent:
    # 1st line: # of input reads
    $align_info[0] =~ /(\d+)\s+reads;/;
    my $input_reads = $1;
    # Fourth line: # of uniquely mapping reads
    $align_info[3] =~ /(\d+)\s+\(\s?(\d+\.\d+\s?%).+ aligned concordantly exactly 1 time/;
    my $unique_reads = $1;
    my $percent_unique = $2;
    # Last line: % overall alignment rate
    $align_info[$#align_info] =~ /(\d+\.\d+\s?%)\s+overall alignment rate/;
    my $overall_alignment = $1;

    my $summary_file = "summary_files/$population\_align\_summary.txt";
    open SUMMARY, ">>$summary_file" or die "ERROR: Could not open $summary_file\n";
    print SUMMARY "$sample\t$input_reads\t$unique_reads\t$percent_unique\t$overall_alignment\n";
    close SUMMARY or die "ERROR: Could not close $summary_file\n";
}

1;
