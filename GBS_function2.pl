#!/usr/bin/perl -w

##### STEP 2 : TRIM DEMULTIPLEXED READS #####
##### Usage: f2 trimmomatic_path trim_file
##### Required input:
#####   trimmomatic_path : The full pathname to the trimmomatic jar file
#####   trim_file : A list of sequences to filter out the reads (ie. Illumina adaptors)
##### Output:
#####   $output_dir/trim/$index_$sample_R1-s.fastq
#####   $output_dir/trim/$index_$sample_R1-p.fastq
#####   $output_dir/trim/$index_$sample_R2-s.fastq
#####   $output_dir/trim/$index_$sample_R2-p.fastq
#####   $output_dir/trim/$index_$sample_output.log

use strict;
use warnings;

sub f2
{
    #############################################
    ##### DEFAULT VARIABLES FOR TRIMMOMATIC #####
    #############################################
    # Version
        my $version = '0.33';
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
    my $output_dir = $_[4];

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

    system("mkdir -p $output_dir/trim/");
    system("mkdir -p logs/");

    # Create summary file
    system("mkdir -p summary_files/");
    my $summary_file = "summary_files/$sample\_trim\_summary.txt";
    open SUMMARY, ">$summary_file" or die "ERROR: Could not create $summary_file\n";
    print SUMMARY "Barcode\tInput Read Pairs\tSurviving Read Pairs\t% Both Surviving\t",
        "Only Forward Surviving\t% Forward Surviving\tOnly Reverse Surviving\t% Reverse Surviving\t",
        "Dropped Reads\t% Dropped\n";
    close SUMMARY or die "Unable to close $summary_file\n";

    my @indices = `cat $index_file`;
    my $num_indices = $#indices + 1;
    if ($num_indices == 0){    die "ERROR: $index_file exists but appears to be empty.\n";    }
    my $index_count = 0;
    print_progress($index_count, $num_indices);

    # BEGIN trimming by index
    foreach my $index (@indices)
    {
        chomp($index);

        my $R1_reads = "$output_dir/demultiplex/$index\_$sample\_R1-clip.fastq";
        my $R2_reads = "$output_dir/demultiplex/$index\_$sample\_R2.fastq";

        unless ( -f "$R1_reads" && -r "$R1_reads" )
        {   die "ERROR: $R1_reads does not exist or is unreadable.\n"; }
        unless ( -f "$R2_reads" && -r "$R2_reads" )
        {   die "ERROR: $R2_reads does not exist or is unreadable.\n"; }

        # Run Trimmomatic
        my $cmd = "java -classpath $trimmomatic_path org.usadellab.trimmomatic.TrimmomaticPE ";
        $cmd .= "-phred33 ";
        $cmd .= "-trimlog logs/$index.trim.log ";
        $cmd .= "$R1_reads $R2_reads ";
        $cmd .= "$output_dir/trim/$index\_$sample\_R1-p.fastq $output_dir/trim/$index\_$sample\_R1-s.fastq ";
        $cmd .= "$output_dir/trim/$index\_$sample\_R2-p.fastq $output_dir/trim/$index\_$sample\_R2-s.fastq ";
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
            #print "Trimmomatic successfully finished trimming $index indexed reads.\n";
            my $trimlog = "logs/$index\_$sample\_trimmomatic\_output.log";
            open TRIMLOG, ">$trimlog";
            print TRIMLOG join " ", @$full_buf;
        } else
        {
            print "ERROR: Trimmomatic did not complete successfully:\n";
            die "$error_message\n@$stderr_buf\n";
        }

        # Summary of progress
        $index_count++;
        print_progress($index_count, $num_indices, "Current index: $index");
        summarize_trim($index, $sample, $output_dir, \@$stderr_buf);
    }
    print "\n",
        " Processed reads located in:\n  $output_dir/trim/ \n",
        " Summary (open in Excel or use command more): $summary_file\n";
}

#################################
##### Summarize step 2 using trimmomatic output:
##### Index    Input read pairs    Both surviving   Both surviving %    Forward only surviving  Forward surviving %
#####   Reverse only surviving  Reverse surviving %     Dropped    Dropped %
#################################
sub summarize_trim
{
    my $index = $_[0];
    my $sample = $_[1];
    my $output_dir = $_[2];
    my @trimmomatic_output = @{$_[3]};

    my $summary_file = "summary_files/$sample\_trim\_summary.txt";
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
    my $dropped_percent = $trim_info[20];

    print SUMMARY "$index\t$input_read_pairs\t$both_surviving\t$both_surviving_percent\t",
        "$forward_surviving\t$forward_surviving_percent\t$reverse_surviving\t",
        "$reverse_surviving_percent\t$dropped\t$dropped_percent\n";

    close SUMMARY or die "ERROR: Could not close $summary_file\n";
}

1;
