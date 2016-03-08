#!/usr/bin/perl -w

##### STEP 2 : TRIM DEMULTIPLEXED READS #####
##### Usage: trim_reads trimmomatic_path trim_file
##### Required input:
#####   trimmomatic_path : The full pathname to the trimmomatic jar file
#####   trim_file : A list of sequences to filter out the reads (ie. Illumina adaptors)
##### Output:
#####   $output_dir/trim/$sample_$population_R1-s.fastq
#####   $output_dir/trim/$sample_$population_R1-p.fastq
#####   $output_dir/trim/$sample_$population_R2-s.fastq
#####   $output_dir/trim/$sample_$population_R2-p.fastq
#####   $output_dir/trim/$sample_$population_output.log

use strict;
use warnings;

sub f2
{
    #############################################
    ##### DEFAULT VARIABLES FOR TRIMMOMATIC #####
    #############################################
    my %trim_options = (
    # Version
        'VERSION' => '0.33',
    # THREADS
        'TRIM_THREADS' => '1',
    # ILLUMINACLIP:
        'SEED_MISMATCHES' => '2',
        'PALINDROME_CLIP_THRESHOLD' => '30',
        'SIMPLE_CLIP_THRESHOLD' => '10',
    # SLIDINGWINDOW:
        'WINDOW_SIZE' => '4',
        'REQUIRED_QUALITY' => '15',
    # LEADING:
        'LEADING_QUALITY' => '3',
    # TRAILING:
        'TRAILING_QUALITY' => '3',
    # MINLEN:
        'MINLEN' => '36',
    );
    #############################################

    # Collect function-specific parameters
    my $trimmomatic_path = $_[0];
    my $trim_file = $_[1];
    my $population = $_[2];
    my $index_file = $_[3];
    my $output_dir = $_[4];

    # Collect trimmomatic-specific options
    if ($_[5])
    {
        my %options = %{$_[5]};
        while ( my ($key, $value) = each %options )
        {
            if ($options{ $key })
            {
                $trim_options{ $key } = $value;
            }
        }
    }

    # Ensure index_file exists in the cwd
    unless ( -f $index_file && -r $index_file )
    {   die "ERROR: $index_file does not exist or is unreadable.\n";    }

    # @TODO: Lookup the version of the user's trimmomatic
    unless ( -s "$trimmomatic_path" )
    {

    #}
    #else
    #{
        print "WARNING: Can't locate Trimmomatic at $trimmomatic_path.\n",
              "Would you like to attempt to install v$trim_options{'VERSION'} there now? (yes/no) ";
        chomp (my $usr_input = <STDIN>);
        if ($usr_input =~ /yes/i)
        {
            system("curl -O http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-$trim_options{'VERSION'}.zip");
            system("unzip Trimmomatic-$trim_options{'VERSION'}.zip; rm -f Trimmomatic-$trim_options{'VERSION'}.zip; mv Trimmomatic-$trim_options{'VERSION'}/ $trimmomatic_path/");
#
            $trimmomatic_path = "Trimmomatic-$trim_options{'VERSION'}/trimmomatic-$trim_options{'VERSION'}.jar";
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

    system("mkdir -p $output_dir/trim/");
    system("mkdir -p logs/");

    # Create summary file
    system("mkdir -p summary_files/");
    my $summary_file = "summary_files/$population\_trim\_summary.txt";
    open SUMMARY, ">$summary_file" or die "ERROR: Could not create $summary_file\n";
    print SUMMARY "Sample\tInput Read Pairs\tSurviving Read Pairs\t% Both Surviving\t",
        "Only Forward Surviving\t% Forward Surviving\tOnly Reverse Surviving\t% Reverse Surviving\t",
        "Dropped Reads\t% Dropped\n";
    close SUMMARY or die "Unable to close $summary_file\n";

    # The sample names are in the first column if the index_file is in 2-column format
    # Otherwise, we only have a list of indices so the following cut command should still work
    my @samples = `cut -f1 $index_file | sort -u`;
    my $num_samples = $#samples + 1;
    if ($num_samples == 0){    die "ERROR: $index_file exists but appears to be empty.\n";    }
    my $sample_count = 0;

    # BEGIN trimming by sample
    foreach my $sample (@samples)
    {
        chomp($sample);
        $sample =~ s/ //g;
        print_progress($sample_count, $num_samples, "Current sample: $sample       ");

        my $R1_reads = "$output_dir/demultiplex/$sample\_$population\_R1.fastq";
        my $R2_reads = "$output_dir/demultiplex/$sample\_$population\_R2.fastq";

        unless ( -f "$R1_reads" && -r "$R1_reads" )
        {   die "ERROR: $R1_reads does not exist or is unreadable.\n"; }
        unless ( -f "$R2_reads" && -r "$R2_reads" )
        {   die "ERROR: $R2_reads does not exist or is unreadable.\n"; }

        # Run Trimmomatic
         my $cmd = "java -classpath $trimmomatic_path org.usadellab.trimmomatic.TrimmomaticPE ";
         $cmd .= "-threads $trim_options{'TRIM_THREADS'} -phred33 ";
         $cmd .= "-trimlog logs/$sample.trim.log ";
         $cmd .= "$R1_reads $R2_reads ";
         $cmd .= "$output_dir/trim/$sample\_$population\_R1-p.fastq ";
         $cmd .= "$output_dir/trim/$sample\_$population\_R1-s.fastq ";
         $cmd .= "$output_dir/trim/$sample\_$population\_R2-p.fastq ";
         $cmd .= "$output_dir/trim/$sample\_$population\_R2-s.fastq ";
         $cmd .= "ILLUMINACLIP:$trim_file:$trim_options{'SEED_MISMATCHES'}:";
         $cmd .= "$trim_options{'PALINDROME_CLIP_THRESHOLD'}:$trim_options{'SIMPLE_CLIP_THRESHOLD'} ";
         $cmd .= "LEADING:$trim_options{'LEADING_QUALITY'} ";
         $cmd .= "TRAILING:$trim_options{'TRAILING_QUALITY'} ";
         $cmd .= "SLIDINGWINDOW:$trim_options{'WINDOW_SIZE'}:$trim_options{'REQUIRED_QUALITY'} ";
         $cmd .= "MINLEN:$trim_options{'MINLEN'} ";

         my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
             run( command => $cmd, verbose => 0 );
         if ($success)
         {
             my $trimlog = "logs/$sample\_$population\_trimmomatic\_output.log";
             open TRIMLOG, ">$trimlog";
             print TRIMLOG join " ", @$full_buf;
         } else
         {
             print "ERROR: Trimmomatic did not complete successfully:\n";
             die "$error_message\n@$stderr_buf\n";
         }

         $sample_count++;
         summarize_trim($sample, $population, $output_dir, \@$stderr_buf);
         print_progress($sample_count, $num_samples, "                              ");
    }
    print "\n",
        " Processed reads located in:  $output_dir/trim/ \n",
        " Summary file:  $summary_file\n";
}

#################################
##### Summarize step 2 using trimmomatic output:
##### Sample    Input read pairs    Both surviving   Both surviving %    Forward only surviving  Forward surviving %
#####   Reverse only surviving  Reverse surviving %     Dropped    Dropped %
#################################
sub summarize_trim
{
    my $sample = $_[0];
    my $population = $_[1];
    my $output_dir = $_[2];
    my @trimmomatic_output = @{$_[3]};

    my $summary_file = "summary_files/$population\_trim\_summary.txt";
    open SUMMARY, ">>$summary_file" or die "ERROR: Could not open $summary_file\n";

    # Grep for the line containing stats about the trimming process
    # Save values of interest into variables
    my ($line_of_interest) = grep /Input Read Pairs/, @trimmomatic_output;
    my @trim_info = split(/ /, $line_of_interest);
    my $input_read_pairs = $trim_info[3];
    my $both_surviving = $trim_info[6];
    my $both_surviving_percent = $trim_info[7];
    my $forward_surviving = $trim_info[11];
    my $forward_surviving_percent = $trim_info[12];
    my $reverse_surviving = $trim_info[16];
    my $reverse_surviving_percent = $trim_info[17];
    my $dropped = $trim_info[19];
    chomp( my $dropped_percent = $trim_info[20] );

    print SUMMARY "$sample\t$input_read_pairs\t$both_surviving\t$both_surviving_percent\t",
        "$forward_surviving\t$forward_surviving_percent\t$reverse_surviving\t",
        "$reverse_surviving_percent\t$dropped\t$dropped_percent\n";

    close SUMMARY or die "ERROR: Could not close $summary_file\n";
}

1;
