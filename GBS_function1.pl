#!/usr/bin/perl -w

##### STEP 1 : DEMULTIPLEX READS #####
##### Usage: f1 sample_name index_file R1_file R2_file
##### Required user input:
#####   sample_name : A sample name for naming files after processing reads (no spaces allowed)
#####   index_file : A file containing the indices (aka barcodes) for distinguishing samples, provided by Illumina
#####   R1_file : A FASTQ file (zipped or unzipped) containing raw Read 1 reads
#####   R2_file : A FASTQ file (zipped or unzipped) containing raw Read 2 reads
##### Output:
#####   $index_sample_R1.fastq for R1 reads beginning with barcode $index
#####   $index_sample_R2.fastq for R2 reads which pair with filtered R1 reads
#####   /tmp/$index_sample_R1.tmp for filtered R1 reads starting with $index (not yet in FASTQ format)
#####   /tmp/$index_sample_R1-clip.fastq for R1 reads trimmed of barcode $index
#####   /tmp/$index.list A list of R1 FASTQ headers starting with $index

use strict;
use warnings;

sub function1
{
    # Collect function-specific parameters
    my $sample = $_[0];
    my $index_file = $_[1];
    unless ( -f $index_file && -r $index_file )
    {   die "ERROR: $index_file does not exist or is unreadable.\n";    }

    my $R1_file = $_[2];
    my $R2_file = $_[3];
    unless ( -f $R1_file && -r $R1_file && -f $R2_file && -r $R2_file )
    {   die "ERROR: Files $R1_file and/or $R2_file do not exist or are unreadable.\n";  }

    # Check the files containing reads that they are in FASTQ format, version Illumina 1.8+
    foreach my $R_file ($R1_file, $R2_file)
    {
        open (READS, $R_file) or die "ERROR: Cannot open $R_file.\n";
        my $id = <READS>;
        my $seq = <READS>;
        my $id2 = <READS>;
        my $qual = <READS>;
        if ($id !~ /^\@/) { die "ERROR: Header does not begin with \@; $R_file is not a valid FASTQ file.\n"; }
        if ($seq !~ /[ACTGN]/i) { die "ERROR: Sequence is not composed of only nucleotides [ACTGN] in $R_file.\n"; }
        if ($id2 !~ /^\+/) { die "ERROR: Second header does not begin with \+; $R_file is not a valid FASTQ file.\n"; }
        if (length($qual) != length($seq)) { die "ERROR: Quality is not the same length as sequence in $R_file.\n"; }
        # Check if the FASTQ format is Illumina version 1.3
        #if ($id =~ /\/[12]$/ ) { die "ERROR: Illumina 1.3 FASTQ file format detected. Please use version 1.8+.\n"; }
        close READS or die "ERROR: Could not close $R_file\n";
    }

    # Place indices in an array
    my @indices = `cat $index_file`;
    my $num_indices = $#indices + 1;
    if ($num_indices == 0){    die "ERROR: $index_file exists but appears to be empty.\n";    }

    my $count = 0;
    system("mkdir -p tmp/");

    foreach my $index (@indices)
    {
        chomp($index);

        ## a.
        # Search R1 reads for each barcode at the start of the sequence,
        # save to a temp file (-A 2 -B 1 options: include 2 lines after, 1 line before match)
        system("zgrep -A 2 -B 1 ^$index $R1_file > tmp/$index\_$sample\_R1.tmp");

        ## b.
        # Remove the "--" separator between each read in each temp file and save as a FASTQ file
        system("grep -v ^- tmp/$index\_$sample\_R1.tmp > $index\_$sample\_R1.fastq");

        ## c.
        # Trim off barcodes from each read using the barcode length after accounting for the
        # restriction site (TGCA) within the barcode
        my $clip = $index;
## TODO # Likely have to request RE site from user to clip from indices
        $clip =~ s/$//;
        my $len = length($clip);
        open CLIPPED, ">$index\_$sample\_R1-clip.fastq" or die "ERROR: Could not create $index\_$sample\_R1-clip\n";
        select(CLIPPED);
        fix_r1("$index\_$sample\_R1.fastq", $len);
        #CLIPPED->flush();
        select(STDOUT);
        close CLIPPED or die "ERROR: Could not close $index\_$sample\_R1-clip\n";

        ## d.
        # Place FASTQ headers for clipped R1 reads into a separate barcode-specific file
        open HEADERS, ">tmp/$index\_$sample.list" or die "ERROR: Could not create tmp/$index\_$sample.list\n";
        select(HEADERS);
        get_id("$index\_$sample\_R1-clip.fastq");
        #HEADERS->flush();
        select(STDOUT);
        close HEADERS or die "ERROR: Could not close tmp/$index.list\n";

        ## e.
        # Use the FASTQ headers from R1 reads to extract R2 reads into index-specific files
        open R2_READS, ">$index\_$sample\_R2.fastq" or die "ERROR: Could not create $index\_$sample\_R2.fastq";
        select(R2_READS);
        get_r2($R2_file, "tmp/$index\_$sample.list");
        #R2_READS->flush();
        select(STDOUT);
        close R2_READS or die "ERROR: Could not close $index\_$sample\_R2.fastq\n";

        # Summary of progress
        $count++;
        print "Demultiplexed barcode $index [$count/$num_indices]\n";
    }

    #create_summary($index_file, $sample);
}

##### Summarize step 1 by listing raw read counts, demultiplexed read counts
sub create_summary
{
    my $index_file = $_[0];
    my $sample = $_[1];
    my $summary_file = "$sample\_demultiplex\_summary.txt";
    open SUMMARY, ">$summary_file" or die "ERROR: Could not create $summary_file\n";
    print SUMMARY "Barcode\tRead1 count\tRead2 count\n";

    my @indices = `cat $index_file`;
    foreach my $index (@indices)
    {
        my ( $read1_count, $read2_count ) = 0;
        my $R1_file = "$index\_$sample\_R1-clip.fastq";
        my $line_count = system("wc -l $R1_file");
        print "LC: $line_count\n";
        $line_count =~ s/^\s+(\d+)\s.*$/$1/;

        # Divide line_count by 4 to determine read count in a fastq file
        $read1_count = ( $line_count / 4 );

        my $R2_file = "$index\_$sample\_R2.fastq";
        $line_count = system("wc -l $R2_file");
        $line_count =~ s/^\s+(\d+)\s/$1/;
        $read2_count = ( $line_count / 4 );

        print SUMMARY "$index\t$read1_count\t$read2_count\n";
    }
    print "Demultiplexing completed - see summary file:\n",
          "\t$summary_file\n";
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

    open(HEAD, $headers) or die "ERROR: Could not open $headers\n";
    while(<HEAD>) {
        chomp;
        # Depending on the read format, here are two cases:
        $_ =~ s/1:N:0:/2:N:0:/;   # ***** Another hardcoded assumption. How to distinguish R1 from R2?? *****
        $_ =~ s/(\w+)\/1$/$1\/2/;   # This works with the test data set
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
        }
    }
    close FH or die "ERROR: Could not close $file\n";
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

1;
