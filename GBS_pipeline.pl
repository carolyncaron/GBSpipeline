#!/usr/bin/perl -w

####################################################################
##### GBS_pipeline.pl - a multi-step pipeline for GBS analysis #####
####################################################################
##### Usage: ./GBS_pipeline.pl [options] function
##### Functions:
#####   - demultiplex (also: f1, function1)
#####   - trim_reads (also: trim, f2, function2)
#####   - align_reads (also: align, f3, function3)
#####   - call_SNPs (also: SNP_calling, f4, function4)
##### Options:
#####   - c (config): Specify name of the config file (Default: GBS.conf)
#####   - help: Read the perldoc
##### Requirements:
#####   - Trimmomatic (v0.17-0.33), Bowtie2 (v2.x.x), SAMtools (v1.x), BCFtools (v1.x)
#####   - A file listing barcodes (indices) for demultiplexing (optionally, each barcode associated with a sample name)
#####   - Paired read data in separate FASTQ files (R1 reads and R2 reads)
#####   - A file containing adaptor and other sequences for trimming purposes
#####   - A FASTA file containing the reference genome
##### See the perldoc for full documention: perldoc GBS_Pipeline.pl
#####   or refer to: github.com/carolyncaron/GBSpipeline/wiki/Using-GBSpipeline-for-the-first-time
#####   for a tutorial.
#######################

use strict;
use warnings;
use File::Basename;
use IO::File;
use feature qw/switch/;
use Time::HiRes;
use FindBin qw($Bin);
use IPC::Cmd qw[run];
use Pod::Usage;
use Data::Dumper;

#######################
##### GET OPTIONS #####
#######################

# Check if command line flags were provided: this overrides any non-flag options
use Getopt::Long;

my ( $config_file, $help );
GetOptions ('c|config=s' => \$config_file,
            'help' => \$help,
           );

pod2usage(-verbose => 2) if $help;

#######################
##### CONFIG FILE #####
#######################

my $CONFIG_FILE;

# Check for the presence of a config file, otherwise assume it is named "GBS.conf"
unless ( $config_file ) {   $config_file = "./GBS.conf";    }

# Parse configuration file into a hash
my %config_hash;
open CONFIG, "<$config_file" or die "ERROR: Unable to open config file $CONFIG_FILE\n";
while (<CONFIG>)
{
    chomp;
    s/#.*//; # Remove comments
    s/^\s+//; # Remove opening whitespace
    s/\s+$//;  # Remove closing whitespace
    next unless length;
    my ($key, $value) = split(/\s*=\s*/, $_, 2);
    $config_hash{$key} = $value;
}

################
##### MAIN #####
################
# Check for user-specified functions to run- if none specified, display POD #

# Check that at least one argument is given
if ( exists ( $ARGV[0] ) )
{
    my $FUNCTION = $ARGV[0];
    shift @ARGV;
    my @args = @ARGV;
    my $num_args = $#args + 1;

    my $start = Time::HiRes::gettimeofday();

    # Regardless of the function specified, each function requires access to these parameters:
    # population, index_file and output_dir
    my ( $population, $index_file, $output_dir, $reference_genome );

    if ($num_args == 0)
    {
        if ($config_hash{'POPULATION'}) { $population = $config_hash{'POPULATION'}; }
            else { report_missing('POPULATION'); }
        if ($config_hash{'INDEX_FILE'}) { $index_file = $config_hash{'INDEX_FILE'}; }
            else { report_missing('INDEX_FILE'); }
        if ($config_hash{'READS_DIR'})  { $output_dir = $config_hash{'READS_DIR'};  }
            else { report_missing('READS_DIR'); }
    }
    # Grab them from the recently created config file from step 1 (or step 3 for reference genome)
    # Note that we don't need these to be non-empty, since this may be the first step and they will
    # be provided at that time.
    else
    {
        chomp($population = `grep 'POPULATION' $config_file | cut -d'=' -f2`);
        chomp($index_file = `grep 'INDEX_FILE' $config_file | cut -d'=' -f2`);
        chomp($output_dir = `grep 'READS_DIR' $config_file | cut -d'=' -f2`);
        chomp($reference_genome = `grep 'REFERENCE' $config_file | cut -d'=' -f2`);
    }

    # Call the appropriate function
    given($FUNCTION)
    {
        ##### DEMULTIPLEX #####
        when ( /function1/ || /f1/ || /demultiplex/ )
        {
            print "Calling $FUNCTION ...\n";
            my ( $RE_site, $R1_file, $R2_file );

            # Check config file for parameters
            if ($num_args == 0)
            {
                if ($config_hash{'RE_SITE'})    { $RE_site    = $config_hash{'RE_SITE'};    }
                    else { $RE_site = ''; }
                if ($config_hash{'R1_FILE'})    { $R1_file    = $config_hash{'R1_FILE'};    }
                    else { report_missing('R1_FILE'); }
                if ($config_hash{'R2_FILE'})    { $R2_file    = $config_hash{'R2_FILE'};    }
                    else { report_missing('R2_FILE'); }
            }
            # Else, parameters are provided via command line
            elsif ($num_args == 6)
            {
                $population = $args[0];
                $index_file = $args[1];
                $RE_site    = $args[2];
                $R1_file    = $args[3];
                $R2_file    = $args[4];
                $output_dir = $args[5];

                # Save the sample and index file into the configuration file
                add_to_config("POPULATION", $population, "The generic name used for this population");
                add_to_config("INDEX_FILE", $index_file, "The filename of the list of indices (aka barcodes)");
                add_to_config("READS_DIR", $output_dir, "The location where output of processed reads are placed");
            }
            # No config file or parameters!
            else
            {
                print "ERROR: Unexpected number of parameters given ($num_args). Program will exit.\n";
                die "Try ./GBS_pipeline.pl -help for options.\n";
            }

            require "$Bin/GBS_function1.pl";
            function1($population, $index_file, $RE_site, $R1_file, $R2_file, $output_dir);

            summarize($start);
        }

        ##### TRIM READS #####
        when ( /function2/ || /f2/ || /trim_reads/ || /trim/ )
        {
            print "Calling $FUNCTION ...\n";
            my ( $trimmomatic_path, $trim_file, %trim_options );

            # Check config file for parameters
            if ($num_args == 0)
            {
                if ($config_hash{'TRIMMOMATIC_PATH'}) { $trimmomatic_path = $config_hash{'TRIMMOMATIC_PATH'}; }
                    else { report_missing('TRIMMOMATIC_PATH'); }
                if ($config_hash{'TRIM_FILE'}) { $trim_file = $config_hash{'TRIM_FILE'}; }
                    else { report_missing('TRIM_FILE'); }

                # Save options for Trimmomatic in a hash
                %trim_options = map { $_ => $config_hash{$_} } qw/TRIM_THREADS SEED_MISMATCHES PALINDROME_CLIP_THRESHOLD SIMPLE_CLIP_THRESHOLD WINDOW_SIZE REQUIRED_QUALITY LEADING_QUALITY TRAILING_QUALITY MINLEN/;
            }
            # Check for parameters on command line
            elsif ($num_args == 2)
            {
                $trimmomatic_path = $args[0];
                $trim_file = $args[1];
            }
            else
            {
                print "ERROR: Unexpected number of parameters given ($num_args). Program will exit.\n";
                die "Try ./GBS_pipeline.pl -help for options.\n";
            }

            require "$Bin/GBS_function2.pl";
            f2($trimmomatic_path, $trim_file, $population, $index_file, $output_dir, \%trim_options);

            summarize($start);
        }

        ##### ALIGN READS #####
        when ( /function3/ || /f3/ || /align_reads/ || /align/ )
        {
            print "Calling $FUNCTION ...\n";
            my ( $bowtie2_dir, %align_options );

            # Check config file for parameters
            if ($num_args == 0)
            {
                if ($config_hash{'BOWTIE2_PATH'}) { $bowtie2_dir = $config_hash{'BOWTIE2_PATH'}; }
                    else { report_missing('BOWTIE2_PATH'); }
                if ($config_hash{'REFERENCE'}) { $reference_genome = $config_hash{'REFERENCE'}; }
                    else { report_missing('REFERENCE'); }

                # Save options for Bowtie2 in a hash
                %align_options = map { $_ => $config_hash{$_} } qw/ALIGN_THREADS MAX_VALID_ALIGNMENTS MAX_FRAGMENT_LENGTH MAX_RESEED_RATE/;
            }
            # Check for parameters on command line
            elsif ($num_args == 2)
            {
                $bowtie2_dir = $args[0];
                $reference_genome = $args[1];
                add_to_config("REFERENCE",$reference_genome,"The pathname of the reference genome sequence.");
            }
            else
            {
                print "ERROR: Unexpected number of parameters given ($num_args). Program will exit.\n";
                die "Try ./GBS_pipeline.pl -help for options.\n";
            }

            require "$Bin/GBS_function3.pl";
            f3($bowtie2_dir, $reference_genome, $population, $index_file, $output_dir, \%align_options);

            summarize($start);
        }

        ##### SNP CALLING #####
        when ( /function4/ || /f4/ || /SNP_calling/ || /call_SNPs/ )
        {
            print "Calling $FUNCTION ...\n";
            my ( $samtools_dir, $bcftools_dir, $call_mode, $samtools_threads );

            # Check config file for parameters
            if ($num_args == 0)
            {
                if ($config_hash{'SAMTOOLS_PATH'}) { $samtools_dir = $config_hash{'SAMTOOLS_PATH'}; }
                    else { report_missing('SAMTOOLS_PATH'); }
                if ($config_hash{'BCFTOOLS_PATH'}) { $bcftools_dir = $config_hash{'BCFTOOLS_PATH'}; }
                    else { $bcftools_dir = $samtools_dir; }
                if ($config_hash{'MODE'}) { $call_mode = $config_hash{'MODE'}; }
                    else { $call_mode = "single"; }
                if ($config_hash{'SAMTOOLS_THREADS'}) { $samtools_threads = $config_hash{'SAMTOOLS_THREADS'}; }
                    else { $samtools_threads = "1"; }
                if ($config_hash{'REFERENCE'}) { $reference_genome = $config_hash{'REFERENCE'}; }
                    else { report_missing('REFERENCE'); }
            }
            # Check for parameters on command line
            elsif ($num_args == 2)
            {
                $samtools_dir = $args[0];
                $bcftools_dir = $args[1];
                ## Here the call mode is set to "single" since it is assumed that parameters were provided on the
                ## command line to repeat a run from a time before the pipeline supported a config file, hence it
                ## was also before the pipeline handled SNP calling on multiple samples.
                $call_mode = "single";
            }
            else
            {
                print "ERROR: Unexpected number of parameters given ($num_args). Program will exit.\n";
                die "Try ./GBS_pipeline.pl -help for options.\n";
            }

            require "$Bin/GBS_function4.pl";
            f4($samtools_dir, $bcftools_dir, $call_mode, $population, $index_file, $output_dir, $reference_genome, $samtools_threads);

            summarize($start);
        }

        default
        {
            # Invoke usage information
            pod2usage(-verbose => 1, -msg => "ERROR: $FUNCTION does not exist. Try ./GBS_pipeline.pl -help for options.");
        }
    }
}
else
{
    # If no function is provided, invoke usage information
    pod2usage(-verbose => 1, -msg => "ERROR: No function specified. See ./GBS_pipeline.pl -help for options.");
}

#########################
##### PRINT SUMMARY #####
#########################
# Input: The start time

sub summarize
{
    print "Complete! ";
    my $start = $_[0];
    my $end = Time::HiRes::gettimeofday();
    printf("Time elapsed: %.2f s\n", $end - $start);
}

#########################
##### ADD TO CONFIG #####
#########################
# Input: Parameter name, value and a comment describing the parameter.

sub add_to_config
{
    my $variable_name = $_[0];
    my $value = $_[1];
    my $comment = $_[2];

    open CONFIG, ">>$config_file" or die "ERROR: Unable to open config file $config_file\n";

    # If the variable exists but has been altered, alter the file accordingly
    if ( `grep $variable_name $config_file` )
    {
        # Finds the variable name and replaces the whole line with the variable name + the new value
        # Yes, it's a perl one-liner within a perl script... but it does the trick
        # First prevent any special characters in $value from being used in the regex
        my $cmd = "perl -pi -e s{^$variable_name=.*}{$variable_name=$value}g $config_file";
        my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
            run( command => $cmd, verbose => 0 );
        unless ($success)
        {   print "ERROR: Unable to replace $variable_name in $CONFIG_FILE:\n$error_message\n@$stderr_buf";   }
    } else # The variable does not yet exist in the file, so append it.
    {
        print CONFIG "#$comment\n$variable_name=$value\n\n";
    }
    close CONFIG or die "ERROR: Unable to close config file $config_file\n";
}

####################################
##### REPORT MISSING PARAMETER #####
####################################
# Input: Name of parameter as it appears in the config file.

sub report_missing
{
    my $param = $_[0];
    die "ERROR: Could not find a value for required parameter \"$param\" in $config_file. Please ensure you've filled in all required fields for this step.\n";
}

###########################
##### PROGRESS REPORT #####
###########################
# Input: A current count of the number of steps completed, the total number of steps

sub print_progress
{
    my $step_count = $_[0];
    my $num_steps = $_[1];
    my $message = $_[2];

    select(STDOUT);

    # Check if a message was given, otherwise give it the empty string
    unless ( defined $message and length $message ) { $message = ""; }
    # Calculate percentage of steps completed
    my $percent_complete = ($step_count/$num_steps)*100;
    my $steps;
    # Remove decimal places for simplicity
    $percent_complete = int $percent_complete;

    print " [";
    # Determine the length of the bar
    for ($steps=0; $steps<$percent_complete; $steps=$steps+4)
    {
        print "=";
    }
    # Determine the space remaining after the bar
    my $steps_remaining = (100 - $steps);
    for (my $i = 0; $i < $steps_remaining; $i = $i+4)
    {
        print " ";
    }
    # Output percentage then shift cursor to beginning of the line
    print "] $percent_complete %  $message\r";
}

###############
##### POD #####
###############

=pod

=head1 NAME

GBS_Pipeline.pl

=head1 SYNOPSIS

B<./GBS_Pipeline.pl> [options] B<function> [arg1] [arg2] ...

=head1 DESCRIPTION

A complete set of commands to demultiplex, trim, align and call raw
variants on paired-end reads for the purpose of analyzing reads generated by genotyping
by sequencing (GBS).

The pipeline contains a set of steps (functions) that can be called individually to
complete a GBS analysis given multiplexed paired-end read data (Illumina 1.8+).
Each function requires a specific set of arguments which can be parsed directly
from the provided configuration file (recommended) or directly on the command line
(see below for expected arguments and order).

It is recommended that the GBS_pipeline files be placed in a new directory created
specifically for the GBS analysis to be performed. Running all steps of the pipeline will
create the following directories:

F<demultiplex/ trim/ align/ variants/>

where outputs from each step are placed. Summary files are placed in a separate
F<summary_files> directory.

=head2 OPTIONS

=over 4

=item B<-c> or B<-config>

Specifies the name of the config file (default: F<GBS.conf>)

=item B<-help>

Prints this perldoc

=back

=head2 FUNCTIONS

=over 4

=item B<demultiplex> [POPULATION] [F<INDEX_FILE>] [RE_SITE] [F<R1_FILE>] [F<R2_FILE>] [READS_DIR]

Purpose: Demultiplex reads based on a provided index file.

B<POPULATION> can be any generic name to represent this GBS run. It is only used in the
naming of output files. Avoid the use of whitespace! (Ex: lens culinaris => lens_culinaris)

The B<INDEX_FILE> is a simple text file and can have one of two formats. Either it is a single
column list of the indices used for this GBS run, or it consists of 2 tab-delimited columns:

=over 1

=item 1. The first column contains the sample names

=item 2. The second column consists of the index that is associated with the sample name in the
  first column.

=back

Either format provides enough information for the pipeline to demultiplex the read files,
however it is recommended to use the second format as the output files will be named using
the sample names rather than the indices, which saves a step for the user.

B<RE_SITE> is the overhang sequence (ex. TGCA in the case of ApaL1) from the rare-cutter
restriction enzyme site used in the GBS protocol (2-enzyme GBS only). This parameter is I<only
necessary if the restriction site is appended to each index>. If this is not the case for your
data set, it is strongly suggested that you use the configuration file and leave this blank.

B<R1_FILE> and B<R2_FILE> should be provided in FASTQ format version Illumina 1.8+.
Ensure the full pathname is provided.

B<OUTPUT_DIR> is a user-specified directory for placement of processed reads. This can be
beneficial when running analysis on a machine or server where space is limited, since
very large files can be directed to a separate storage unit or location. You may
specify "." to indicate the current working directory.

=item B<trim_reads> [F<TRIMMOMATIC_PATH>] [F<TRIM_FILE>]

Purpose: Trim low quality bases and remove Illumina adaptors from raw FASTQ reads using a
command-line tool called B<Trimmomatic>  I<(Bolger, A. M., Lohse, M., & Usadel, B. (2014).
Trimmomatic: A flexible trimmer for Illumina Sequence Data. Bioinformatics, btu170.)>.

B<TRIMMOMATIC_PATH> is the full pathname to the user's copy of Trimmomatic. Ensure that
the executable .jar is included in your path, and not just the directory. If Trimmomatic
is not found at TRIMMOMATIC_PATH, you will be prompted to install the latest version
to the current working directory.

B<TRIM_FILE> is the filename for a list of Illumina adaptor sequences or other sequences in
FASTA format which are desired to be trimmed from the raw reads.

Keep in mind that some output from the previous step (demultiplexing reads) is also needed
as input into the trimming step. Thus, the input files are expected to named using a particular
format (F<sample_population_R1.fastq> and F<sample_population_R2.fastq>) and located in
OUTPUT_DIR.

The following options are given to trimmomatic as defaults, but may be altered through the
configuration file (recommended) or by editing the code (see: F<GBS_function2.pl>). Refer to
the manual for trimmomatic for full option descriptions:

 TRIM_THREADS: 1
    The maximum mismatch count to allow a full match
 SEED_MISMATCHES: 2
    Specifies accuracy needed for a match between two
    adapter-ligated reads
 PALINDROME_CLIP_THRESHOLD: 30
    Specifies accuracy needed for a match between any adapter
    and read
 SIMPLE_CLIP_THRESHOLD: 20
    The number of bases to average across in sliding window
    trimming
 WINDOW_SIZE: 4
    The average quality required in a sliding window trimming
    procedure
 REQUIRED_QUALITY: 15
    Threshold for minimum quality of leading bases
 LEADING_QUALITY: 3
    Threshold for minimum quality of trailing bases
 TRAILING_QUALITY: 3
    The minimum length for trimmed reads to be kept
 MINLEN: 36

=item B<align_reads> [F<BOWTIE2_PATH>] [F<REFERENCE>]

Purpose: Aligns FASTQ reads to a reference genome using a command-line tool called
B<bowtie2> I<(Langmead B,
Salzberg S. Fast gapped-read alignment with Bowtie 2. Nature Methods. 2012, 9:357-359.)>.

B<BOWTIE2_PATH> is the full path of the I<directory> (not the executable) of the your copy of bowtie2, and is required.

B<REFERENCE> is the pathname of the reference genome sequence, normally in FASTA format.
bowtie2 uses index files created from the genome, and not the FASTA sequence directly. If index files
have not been pre-built, then the pipeline will create them for you using bowtie2.
If you already created index files and want the pipeline to use them, ensure they are in the
same directory as the reference genome with the same basename as the FASTA sequence file
(ie. reference_genome_file.1.bt2 etc... in addition to reference_genome_file.FASTA).

A small subset of bowtie2 options have been optimized for GBS (particularly with plant genomes).
You can change these in the configuration file or by altering the code (see: F<GBS_function3.pl>).
Please refer to the latest bowtie2 manual for full descriptions of these options.
The options and their defaults (as provided by the pipeline) are:

 ALIGN_THREADS: 1
    the number of parallel search threads to run
 MAX_VALID_ALIGNMENTS: 3
    The max number of distinct, valid alignments for each
    read
 MAX_FRAGMENT_LENGTH: 11000
    The max fragment length for valid paired-end
    alignment
 MAX_RESEED_RATE: 5
    The max number of times reads with repetitive seeds
    will be re-seeded

=item B<SNP_calling> [F<SAMTOOLS_PATH>] [F<BCFTOOLS_PATH>]

Purpose: Processes alignments output by bowtie2 in the previous step, calls raw variants
using B<SAMtools> I<(Li H.*, Handsaker B.*, Wysoker A., Fennell T., Ruan J., Homer N., Marth G.,
Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup (2009) The Sequence
alignment/map (SAM) format and SAMtools. Bioinformatics, 25, 2078-9.)> and calls SNPs using
B<bcftools>.

NOTE: It is not recommended to use different versions of SAMtools and bcftools together.

B<SAMTOOLS_PATH> is the location of the user's copy of SAMtools (v1.0+). Warning: previous versions
will not work!

B<BCFTOOLS_PATH> is the location of the user's copy of bcftools (v1.0+). Warning: previous versions
will not work! If you choose to have both SAMtools and bcftools executables in a common directory
(as was the default prior to v1.0+), then this argument can be left blank in the configuration file
(by default, the pipeline will assume the same path as SAMTOOLS_PATH).

B<MODE> can be set to one of: single, multi or both (case insensitive). Prior to a recent update, this
step would call SNPs for each sample against the reference individually. You can still achieve this by
specifying "single" as the mode for SNP calling in the config file. If you want to call SNPs on the
entire population together, you can specify "multi". If you are interested in both of these, you can
specify "both".

Lastly, you can specify the number of threads for SAMtools to use where it is possible (samtools view
and samtools sort in this case) through the config file:

    B<SAMTOOLS_THREADS>: 1

=back

=head1 SUMMARIES

=over 4

=item F<population_demultiplex_summary.txt>

Provides an overview of demultiplexing the raw reads.

 Sample: The sample name OR the index used to distinguish this sample
 Read1 count: The number of R1 reads that contained the index
 % of Raw Read1: The % of raw R1 reads that contained the index
 Read2 count: The number of R2 reads that contained the index
 $ of Raw Read2: The % of raw R2 reads that contained the index

=item F<population_trim_summary.txt>

Provides an overview of trimming the demultiplexed reads.

 Sample: The sample name OR the index used to distinguish this sample
 Input Read Pairs: The number of paired-end reads prior to trimming
 Surviving Read Pairs: The number of reads where both pairs survived
 % Both Surviving: The % of both paired reads surviving
 Only Forward Surviving: The number of forward-orientated reads that
    survived but the other pair did not
 % Forward Surviving: The % of forward-orientated reads that survived
 Only Reverse Surviving: The number of reverse-oriented reads that
    survived but the other pair did not
 % Reverse Surviving: The % of reverse-orientated reads that survived
 Dropped Reads: The number of reads dropped due to contamination or
    low quality
 % Dropped: The % of reads dropped

=item F<population_align_summary.txt>

Provides an overview of the alignment of reads to a reference genome.

 Sample: The sample name OR the index used to distinguish this sample
 Input Reads: Number of reads when the alignment began
 Total Hits: Total number of successful alignments made
 Overall Alignment Rate: The % of input reads that successfully aligned
 Unique Hits: Number of reads that mapped uniquely to the genome, including the "best" hit for multi-mapped reads
    (NOTE: this is different from the "Unique Reads" column in previous versions of the pipeline, as the best hits
    were not determined until the SNP-calling step. Please now refer to the bowtie2 logs for the previous value.
    Unique Hits is more informative as it indicates the number of reads that are used in the subsequent step.)
 Percent Unique: The percent of unique hits

=back

=head1 AUTHOR

Carolyn Caron - <carolyn[dot]caron[at]usask[dot]ca>

=cut
