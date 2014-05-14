#!/usr/bin/perl -w

# GBS PIPELINE
# MAR 26, 2014
#
# This script contains functions that can be called individually to perform a GBS pipeline
# It is designed to be placed in the cwd (ideally a directory created specifically for
# running the pipeline) along with an indices.list file.
#
#
#

use strict;
use warnings;
use File::Basename;
use IO::File;
use Getopt::Std;
# Only applicable for Perl 5.10+, so hesitant to use it, but better than use Switch ...
use feature qw/switch/;
use Cwd qw(cwd);
use Time::HiRes;

######################
##### USAGE MENU #####
######################

sub usage_menu
{
    use Term::ANSIColor qw(:constants);
    $Term::ANSIColor::AUTORESET = 1;

    print BOLD "GBS_pipeline.pl: ";
    print "A GBS pipeline which aligns raw paired read data to a reference genome and performs raw variant calls.\n\n",
        "To access this usage menu in the future, type Perl GBS_pipeline.pl -u\n\n";
    print BOLD "USAGE:\n";
    print "Perl GBS_pipeline.pl [function] [arg1] [arg2] ... \n\n";
    print BOLD "FILES NEEDED:\n";
    print "indices.list\tA list of indices that were used in sequencing to uniquely tag reads.\n",
        "\n";
    print BOLD "FUNCTIONS:\n";
    print "1. function1 sample_name input1 input2\n",
        "2.\n";
    exit();
}

################
##### MAIN #####
################

# First ensure that "indices.list" is present in the cwd

#print cwd, "\n";

my $index_list = "indices.list";
unless ( -f $index_list && -r $index_list )
{
    print "ERROR: Missing or unable to read indices.list file in the cwd!\n\n";
    usage_menu();
}

# Check for user-specified functions to run- if none specified, display usage menu #

# Check if command line flags were provided: this overrides any non-flag options
my %options;
getopts('u', \%options);

if ($options{u})    { usage_menu(); }

# Check that at least one argument is given
if ( exists ( $ARGV[0] ) )
{
    my $FUNCTION = $ARGV[0];
    shift @ARGV;
    my @args = @ARGV;
    print "User specified function $FUNCTION with args: @args\n";

    # Use a switch statement for accessing functions
    given($FUNCTION)
    {
        when (/f1/)     {
                            my $start = Time::HiRes::gettimeofday();
                            print "Calling $FUNCTION ...\n";
                            f1($index_list);
                            #summarize($FUNCTION, $start);
                        }
        when (/f2/)       { print "Calling $FUNCTION ...\n"; }
        default           { print "$FUNCTION does not exist.\n"; usage_menu(); }
    }
}
else
{
    # Looks at the cwd for files/directories present and estimates at what stage in the
    # pipeline it left off




    # If no expected files present, print the usage menu
    usage_menu();
}

#########################
##### PRINT SUMMARY #####
#########################

sub summarize
{
    my $function = $_[0];
    print "Completed $function.\n";
    my $start = $_[1];
    my $end = Time::HiRes::gettimeofday();
    printf("Time elapsed: %.2f s\n", $end - $start);
    print("The following files were created:\n");
}


##########################
##### USER FUNCTIONS #####
##########################

##### STEP 1 : DEMULTIPLEX READS #####
##### Usage: f1 input1
##### Required input:
#####   A FASTQ file (zipped or unzipped) containing raw Read 1 reads
#####   A FASTQ file (zipped or unzipped) containing raw Read 2 reads
##### Output:
#####   $index_sample_R1.fastq for R1 reads beginning with barcode $index
#####   $index_sample_R2.fastq for R2 reads which pair with filtered R1 reads
#####   /tmp/$index_sample_R1.tmp for filtered R1 reads starting with $index (not yet in FASTQ format)
#####   /tmp/$index_sample_R1-clip.fastq for R1 reads trimmed of barcode $index
#####   /tmp/$index.list A list of R1 FASTQ headers starting with $index

sub f1
{
    # Collect function-specific parameters
    my $index_file = $_[0];
    my ($sample, $R1_file, $R2_file);
    my $num_args = $#ARGV + 1;
    print "Index file: $index_file\n ARGS: @ARGV \n# of ARGS: $num_args\n";

    if (($num_args) == 3)
    {
        $sample = $ARGV[0];
        $R1_file = $ARGV[1];
        $R2_file = $ARGV[2];
        unless ( -f $R1_file && -r $R1_file && -f $R2_file && -r $R2_file )
        {
            die "ERROR: Files $R1_file and/or $R2_file do not exist or are unreadable.\n";
        }

        # Check the files containing reads that they are in FASTQ format, version Illumina 1.8+
        foreach my $R_file ($R1_file, $R2_file)
        {
            open (READS, $R_file) or die "Cannot open $R_file.\n";
            my $id = <READS>;
            my $seq = <READS>;
            my $id2 = <READS>;
            my $qual = <READS>;
            if ($id !~ /^\@/) { die "Header does not begin with \@; $R_file is not a valid FASTQ file.\n"; }
            if ($seq !~ /[ACTGN]/i) { die "Sequence is not composed of only nucleotides [ACTGN] in $R_file.\n"; }
            if ($id2 !~ /^\+/) { die "Second header does not begin with \+; $R_file is not a valid FASTQ file.\n"; }
            if (length($qual) != length($seq)) { die "Quality is not the same length as sequence in $R_file.\n"; }
            # Check if the FASTQ format is Illumina version 1.3
            #if ($id =~ /\/[12]$/ ) { die "Illumina 1.3 FASTQ file format detected. Please use version 1.8+.\n"; }
            close READS;
        }

    }
    else    # User did not provide sufficient arguments; prompt
    {
        print "Unexpected number of parameters given ($num_args). Program will exit.\n";
        die "Try: Perl GBS_pipeline.pl f1 sample_name /path/to/file1/filename1.fastq /path/to/file2/filename2.fastq\n";
    }

    open (INDICES, $index_file) or die "indices.list exists but cannot be opened.\n";
    while ( my $index = <INDICES> )
    {
        chomp($index);
        system("mkdir -p tmp/");
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
# ***** Likely have to request RE site from user to clip from indices
        $clip =~ s/$//;
        my $len =length($clip);
        open CLIPPED, ">tmp/$index\_$sample\_R1-clip.fastq";
        select CLIPPED;
        fix_r1("$index\_$sample\_R1.fastq", $len);
        select STDOUT;
        ## d.
        # Place FASTQ headers for filtered R1 reads into a separate barcode-specific file
        # ****** Since @ is also a quality symbol, this is not a reliable method for all test cases!
        # ****** system("grep ^@ tmp/$index\_$sample\_R1-clip.fastq >$index.list");
        open HEADERS, ">tmp/$index.list";
        select HEADERS;
        get_id("tmp/$index\_$sample\_R1-clip.fastq");
        select STDOUT;
        ## e.
        # Use the FASTQ headers from R1 reads to extract R2 reads into temporary files
        open (R2_READS, ">".$R2_file);
        select R2_READS;
        get_r2($R2_file, "$index.list");
        select STDOUT;
        #Clean up
        close CLIPPED; close HEADERS; close R2_READS;
    }
    close INDICES;
}

##################################
##### ADDITIONAL SUBROUTINES #####
##################################

## Written by Larissa
# Input: A FASTQ file of filtered R1 reads, barcode length
sub fix_r1
{
    my $seq_file = $_[0];
    my $bcl = $_[1];

    open(FH, $seq_file);

    while(<FH>) {
        my $q1    = $_;
        my $seq1  = <FH>;
        $seq1 =~ s/^.{$bcl}//;
        my $qq1   = <FH>;
        my $qual1 = <FH>;
        $qual1 =~ s/^.{$bcl}//;
        print $q1.$seq1.$qq1.$qual1;
    }
    close FH;
}

## Written by Larissa
# Input: A FASTQ file of raw R2 reads, a list of R1 read headers
sub get_r2
{
    my $file = $_[0];
    my $headers = $_[1];
    my $href;

    open(HEAD, $headers);
    while(<HEAD>) {
        chomp;
        # Depending on the read format, here are two cases:
        $_ =~ s/1:N:0:/2:N:0:/;   # ***** Another hardcoded assumption. How to distinguish R1 from R2?? *****
        $_ =~ s/(\w+)\/1$/$1\/2/;   # This works with the test data set
        $href->{$_}->{'found'} = 1;
    }

    open(FH, $file);
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
    close FH;
}

# Extracts the header from each read in a provided FASTQ file
sub get_id
{
    my $file = $_[0];

    open FILE, $file;
    while(<FILE>)
    {
        my $id1      = $_;
        my $seq     = <FILE>;
        my $id2     = <FILE>;
        my $qual    = <FILE>;
        print $id1;
    }
}
