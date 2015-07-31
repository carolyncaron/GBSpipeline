#!/usr/bin/perl -w

##### STEP 1 : DEMULTIPLEX READS #####
##### Usage: demultiplex population index_file RE_site R1_file R2_file output_dir
##### Required user input:
#####   population : A name for the population in the GBS run (no whitespace allowed)
#####   index_file : A file containing the indices (aka barcodes) for distinguishing samples
#####   RE_site: A short nucleotide string representing the rare-cutter restriction site used in extracting the reads
#####   R1_file : A FASTQ file (zipped or unzipped) containing raw Read 1 reads
#####   R2_file : A FASTQ file (zipped or unzipped) containing raw Read 2 reads
#####   output_dir : A directory in which processed reads will be placed
##### Output:
#####   $output_dir/demultiplex/$sample_population_R1.fastq for trimmed, filtered R1 reads belonging to $sample
#####   $output_dir/demultiplex/$sample_population_R2.fastq for R2 reads which pair with filtered R1 reads
#####   $output_dir/demultiplex/unindexed_population_R1.fastq for unindexed R1 reads
#####   $output_dir/demultiplex/unindexed_population_R2.fastq for unindexed R2 reads

use strict;
use warnings;

sub function1
{
    # Collect function-specific parameters and check their validity
    my $population = $_[0];
    my $index_file = $_[1];
    unless ( -f $index_file && -r $index_file )
    {   die "ERROR: $index_file does not exist or is unreadable.\n";    }

    my $RE_site = $_[2];
    unless ( $RE_site !~ /(ATGC)/i ) { die "ERROR: Restriction enzyme site must contain only nucleotides [ACGT].\n"; }

    my $R1_file = $_[3];
    my $R2_file = $_[4];
    unless ( -f $R1_file && -r $R1_file && -f $R2_file && -r $R2_file )
    {   die "ERROR: Files $R1_file and/or $R2_file do not exist or are unreadable.\n";  }

    # Determine in which directory to place processed reads
    my $output_dir = $_[5];
    unless ( -d $output_dir )
    {
        print " No directory $output_dir exists. Create it? (yes/no) ";
        chomp (my $usr_input = <STDIN>);
        if ($usr_input =~ /yes/i) { system("mkdir -p $output_dir"); }
        else {
            die "Exiting. Please create a directory for processed reads or use ./ for the current directory.\n";
        }
    }

    # Check the files containing reads that they are in FASTQ format, version Illumina 1.8+
    foreach my $R_file ($R1_file, $R2_file)
    {
        open (READS, $R_file) or die "ERROR: Cannot open $R_file.\n";
        my $id = <READS>;
        my $seq = <READS>;
        my $id2 = <READS>;
        my $qual = <READS>;
        if ($id !~ /^\@/) { die "ERROR: Header does not begin with \@; $R_file is not a valid FASTQ file.\n"; }
        if ($seq !~ /[ACTGN]/i) { die "ERROR: Sequence is not composed of only nucleotides [ACGTN] in $R_file.\n"; }
        if ($id2 !~ /^\+/) { die "ERROR: Second header does not begin with \+; $R_file is not a valid FASTQ file.\n"; }
        if (length($qual) != length($seq)) { die "ERROR: Quality is not the same length as sequence in $R_file.\n"; }
        # Check if the FASTQ format is Illumina version 1.3
        #if ($id =~ /\/[12]$/ ) { die "ERROR: Illumina 1.3 FASTQ file format detected. Please use version 1.8+.\n"; }
        close READS or die "ERROR: Could not close $R_file\n";
    }

    # Create a directory for demultiplexed files
    system("mkdir -p $output_dir/demultiplex/");
    # Create a directory for files that will only be kept temporarily
    system("mkdir -p $output_dir/tmp/");

    # Check that the index_file is valid. Either:
    # 1. It just contains the barcodes, each on its own line (1 column)
    # 2. Each line contains the name of the population followed by the barcode (tab delimited thus 2 columns)
    my %samples = ();
    my ( $sname, $barcode );
    my $num_indices = 0;

    open (INDICES, $index_file) or die "ERROR: Unable to open $index_file.\n";
    my $first_line = <INDICES>;
    chomp $first_line;
    my @columns = split(/\t/, $first_line);
    my $num_col = ($#columns)+1;
    if ($num_col > 2) { die "ERROR: Unexpected number of columns in $index_file: $num_col \n"; }
    # Since we already looked at the first line, add its contents to the array of indices and samples hash,
    # then iterate over the remainder of the file
    elsif ($num_col == 2)
    {
        unless ( $columns[1] !~ /(ATGCN)/i ) { die "ERROR: Indices must contain only nucleotides [ACGTN].\n"; }
        push @{ $samples{ $columns[0] }}, $columns[1];
        $num_indices++;
    }
    elsif ($num_col == 1)
    {
        unless ( $columns[0] !~ /(ATGCN)/i ) { die "ERROR: Indices must contain only nucleotides [ACGTN].\n"; }
        push @{ $samples{ $columns[0] }}, $columns[0];
        $num_indices++;
    }
    else {  die "ERROR: $index_file appears to be empty.\n";  }

    # The first index is inserted into the %samples hash, now enter the remainder.
    # Key -> sample, Value -> index (sample = index if single column file was provided)
    while (<INDICES>)
    {
        chomp;
        if ($num_col == 2)
        {   ($sname, $barcode) = split (/\t/, $_); }
        else
        {    $sname = $_; $barcode = $sname;  }

        # Remove whitespace
        $sname =~ s/ //g;

        #@TODO: If a barcode is duplicated for a sample (sadly I've run into this), ensure it isn't added to the hash twice
        push @{ $samples{ $sname }}, $barcode;

        $num_indices++;
    }
    close INDICES or die "ERROR: Unable to close $index_file\n";
    if ($num_indices == 0){    die "ERROR: $index_file exists but appears to be empty.\n";    }

    # Start the progress bar since wc command can take a while depending on the file.
    my $index_count = 0;
    print_progress($index_count, $num_indices, "Counting reads...");

    # Keep counts on the curent index, R2 reads for an index, and the total raw read pair counts
    my ( $R2_count, $R1_raw_count, $R2_raw_count ) = 0;
    my $wc_cmd = "wc -l $R1_file $R2_file";
    my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
        run( command => $wc_cmd, verbose => 0 );
    if ($success)
    {
        # Split the buffer to extract only the line count
        my @wc_buffer = split(' ',"@$stdout_buf");
        $R1_raw_count = ($wc_buffer[0] / 4);
        $R2_raw_count = ($wc_buffer[2] / 4);
    } else
    {
        die "Unable to do line counts of files $R1_file and $R2_file: $error_message\n@$stderr_buf\n";
    }

    # Create summary file
    system("mkdir -p summary_files/");
    my $summary_file = "summary_files/$population\_demultiplex\_summary.txt";
    open SUMMARY, ">$summary_file" or die "ERROR: Could not create $summary_file\n";
    print SUMMARY "Sample\tRead1 count\t% of Raw Read1\tRead2 count\t% of Raw Read2\n";
    close SUMMARY or die "ERROR: Unable to close $summary_file\n";

    # BEGIN demultiplexing
    my $sample;
    foreach $sample ( sort keys %samples )
    {
        # Declare files for each sample which may need to be appended-to if there are
        # multiple indices for this sample. Hence going through the trouble of making
        # sure they're empty if the pipeline gets re-run...
        my $R1_matches = "$output_dir/demultiplex/$sample\_$population\_R1.fastq";
        my $R2_matches = "$output_dir/demultiplex/$sample\_$population\_R2.fastq";

        # Make sure we are starting with empty output files
        my @output_files = ( "$R1_matches", "$R2_matches" );
        foreach my $outfile (@output_files) { system(">$outfile");  }

        # Reset R2 counter
        $R2_count = 0;

        # Iterate through each index for this sample
        foreach my $index (@{ $samples{$sample} })
        {
            chomp($index);
            $index =~ s/ //g;

            # Create the index-specific files
            my $R1_with_index = "$output_dir/tmp/$index\_$population\_R1-unclipped.fastq";
            my $R1_clipped = "$output_dir/tmp/$index\_$population\_R1-clipped.fastq";
            my $R1_headers = "$output_dir/tmp/$index\_$population.list";

            # Four steps per index. Break it down so the progress bar reports more often.
            my $step_count = $index_count*4;
            my $total_steps = $num_indices*4;
            print_progress($step_count++, $total_steps, " Current sample: $sample        ");

            ## 1. Search R1 reads for each barcode at the start of the sequence
            # (zgrep -A 2 -B 1 options: include 2 lines after, 1 line before match)
            # Then remove the "--" separator between each read and save as a FASTQ file
            system("zgrep -A 2 -B 1 ^$index $R1_file | grep -v ^- > $R1_with_index");
            print_progress($step_count++, $total_steps);

            ## 2.  Trim off barcodes from each read using the barcode length after accounting
            # for the restriction site within the barcode
            my $clip = $index;
            $clip =~ s/$RE_site$//i;
            my $len = length($clip);
            open CLIPPED, ">$R1_clipped" or die "ERROR: Could not create $R1_clipped\n";
            select(CLIPPED);
            fix_r1($R1_with_index, $len);
            #CLIPPED->flush();
            select(STDOUT);
            close CLIPPED or die "ERROR: Could not close $R1_clipped\n";
            system("cat $R1_clipped >> $R1_matches");
            print_progress($step_count++, $total_steps);

            ## 3. Place FASTQ headers for clipped R1 reads into a separate barcode-specific file
            open HEADERS, ">$R1_headers" or die "ERROR: Could not create $R1_headers\n";
            select(HEADERS);
            get_id($R1_clipped);
            #HEADERS->flush();
            select(STDOUT);
            close HEADERS or die "ERROR: Could not close $R1_headers\n";
            print_progress($step_count++, $total_steps);

            ## 4. Use the FASTQ headers from R1 reads to extract R2 reads into sample-specific files
            open R2_READS, ">>$R2_matches" or die "ERROR: Could not create $R2_matches\n";
            select(R2_READS);
            $R2_count = ($R2_count + get_r2($R2_file, "$R1_headers"));
            #R2_READS->flush();
            select(STDOUT);
            close R2_READS or die "ERROR: Could not close $R2_matches\n";
            print_progress($step_count++, $total_steps, "                                ");

            $index_count++;
        }
        summarize_demultiplex($sample, $population, $output_dir, $R2_count, $R1_raw_count, $R2_raw_count, $samples{$sample});
    }

    # Search for reads that are missing a barcode
    # First create a master file containing headers of R1 reads with found indices
    print "\n Handling reads without an index... ";
    my $indexed_list = "$output_dir/tmp/indexed.list";
    my $concat_cmd = "cat $output_dir/tmp/*\_$population.list > $indexed_list";
    ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
            run( command => $concat_cmd, verbose => 0 );
    if ($success)
    {
        foreach my $read_pair ('R1', 'R2')
        {
            my $unindexed = "$output_dir/demultiplex/unindexed\_$population\_$read_pair.fastq";
            my $read_file;
            if ($read_pair eq 'R1')
            {
                $read_file = $R1_file;
            } else
            {
                $read_file = $R2_file;
            }
            open NO_INDEX, ">$unindexed" or die "ERROR: Could not create $unindexed\n";
            select(NO_INDEX);
            $R2_count = get_unindexed($read_file, $indexed_list, $read_pair);
            select(STDOUT);
            close NO_INDEX or die "ERROR: Could not close $unindexed\n";
        }
        summarize_demultiplex("unindexed", $population, $output_dir, $R2_count, $R1_raw_count, $R2_raw_count);
    }
    else
    {
        print "ERROR: Could not concatenate R1 header files into one file:\n";
        die "$error_message\n@$stderr_buf\n";
    }
    print "\n",
          " Processed reads located in:\n  $output_dir/demultiplex/\n",
          " Summary: $summary_file\n";
}

#################################
##### Summarize step 1 by listing demultiplexed read counts
##### Index   Read 1 count  % Raw R1    Read 2 count   % Raw R2
#################################
sub summarize_demultiplex
{
    my $sample = $_[0];
    my $population = $_[1];
    my $output_dir = $_[2];
    my $R2_count = $_[3];
    my $R1_raw_count = $_[4];
    my $R2_raw_count = $_[5];
    my @indices;
    if ( defined $_[6] ){ @indices = @{$_[6]}; }

    my $summary_file = "summary_files/$population\_demultiplex\_summary.txt";
    open SUMMARY, ">>$summary_file" or die "ERROR: Could not open $summary_file\n";

    my $R1_count = 0;
    my $R1_demultiplexed;
    my $wc_cmd;
    if ($sample eq 'unindexed')
    {
        #$R1_demultiplexed = "$output_dir/tmp/indexed.list";
        $wc_cmd = "wc -l $output_dir/tmp/indexed.list";
    }
    else
    {
        #$R1_demultiplexed = "$output_dir/tmp/$sample\_$population.list";
        $wc_cmd = "cat ";
        foreach my $index ( @indices ) { $wc_cmd .= "$output_dir/tmp/$index\_$population.list "; }
        $wc_cmd .= "| wc -l";
    }
    #my $wc_cmd = "wc -l $R1_demultiplexed";
    my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
        run( command => $wc_cmd, verbose => 0 );
    if ($success)
    {
        # Split the buffer to extract only the line count
        my @wc_buffer = split(' ',"@$stdout_buf");
        # No need to divide by 4 as this file only contained FASTQ header lines
        $R1_count = $wc_buffer[0];
        #print "R1_count for $sample: $R1_count\n";
    }
    else
    {
        die "ERROR: Unable to do line count of file $R1_demultiplexed: $error_message\n@$stderr_buf\n";
    }

    # If summarizing unindexed reads, adjust $R1_count (currently counts # of indexed reads)
    if ($sample eq 'unindexed')
    {
        $R1_count = ( $R1_raw_count - $R1_count );
    }

    # Calculate the percentage of raw reads demultiplexed
    my $R1_percent = ( $R1_count / $R1_raw_count ) * 100;
    my $R2_percent = ( $R2_count / $R2_raw_count ) * 100;

    printf (SUMMARY "%s\t%d\t%.2f%%\t%d\t%.2f%%\n", $sample, $R1_count, $R1_percent, $R2_count, $R2_percent);
    close SUMMARY or die "ERROR: Could not close $summary_file\n";
}

##################################
##### ADDITIONAL SUBROUTINES #####
##################################

## Written by Larissa
# Input: A FASTQ file of filtered R1 reads, barcode length (bcl)
sub fix_r1
{
    my $seq_file = $_[0];
    my $bcl = $_[1];

    open(FH, $seq_file) or die "ERROR: Could not open $seq_file\n";

    while(<FH>) {
        my $q1    = $_;
        my $seq1  = <FH>;
        $seq1 =~ s/^.{$bcl}//;
        my $qq1   = <FH>;
        my $qual1 = <FH>;
        $qual1 =~ s/^.{$bcl}//;
        print $q1.$seq1.$qq1.$qual1;
    }
    close FH or die "ERROR: Could not close $seq_file\n";
}

## Written by Larissa
# Input: A FASTQ file of raw R2 reads, a list of R1 read headers
sub get_r2
{
    my $file = $_[0];
    my $headers = $_[1];
    my $href;
    my $r2_count = 0;

    open(HEAD, $headers) or die "ERROR: Could not open $headers\n";
    while(<HEAD>) {
        chomp;
        # Depending on the read format, here are two cases:
        $_ =~ s/1:N:0:/2:N:0:/;
        $_ =~ s/(\w+)\/1$/$1\/2/;
        $href->{$_}->{'found'} = 1;
    }
    close HEAD or die "ERROR: Could not close $headers\n";

    open(FH, $file) or die "ERROR: Could not open $file\n";
    while(<FH>)
    {
        my $q1    = $_;
        my $search = $q1;
        chomp $search;
        my $seq1  = <FH>;
        my $qq1   = <FH>;
        my $qual1 = <FH>;
        if (defined $href->{$search}->{'found'}){
            print $q1.$seq1.$qq1.$qual1;
            $r2_count++;
        }
    }
    close FH or die "ERROR: Could not close $file\n";
    return $r2_count;
}

# Extracts the header from each read in a provided FASTQ file
sub get_id
{
    my $file = $_[0];

    open FILE, $file or die "ERROR: Could not open $file\n";
    do
    {
        my $id1     = <FILE>;
        my $seq     = <FILE>;
        my $id2     = <FILE>;
        my $qual    = <FILE>;

        print $id1;
    }
    while (!eof(FILE));
    close FILE or die "ERROR: Could not close $file\n";
}

## Written by Larissa
# Input: A file of raw reads (reads), a file of read headers that have been demultiplexed (indexed_reads)
sub get_unindexed
{
    my $reads = $_[0];
    my $indexed_reads = $_[1];
    my $read_pair = $_[2];
    my $href;
    my $read_count = 0;

    # Iterate through all headers of R1 reads that have a barcode, save in a hash
    open(IDX, $indexed_reads) or die "ERROR: Could not open $indexed_reads\n";
    while(<IDX>)
    {
        chomp;
        if ($read_pair eq 'R2')
        {
            $_ =~ s/1:N:0:/2:N:0:/;
            $_ =~ s/(\w+)\/1$/$1\/2/;
        }
        $href->{$_}->{'found'} = 1;
    }
    close IDX or die "ERROR: Could not close $indexed_reads\n";

    # Iterate through the raw reads and if not present in the hash, print to a new file
    open(FH, $reads) or die "ERROR: Could not open $reads\n";
    while(<FH>)
    {
        my $q1    = $_;
        my $search = $q1;
        chomp $search;
        my $seq1  = <FH>;
        my $qq1   = <FH>;
        my $qual1 = <FH>;
        if (!defined $href->{$search}->{'found'}){
            print $q1.$seq1.$qq1.$qual1;
            $read_count++;
        }
    }
    close FH or die "ERROR: Could not close $reads\n";
    return $read_count;
}

1;
