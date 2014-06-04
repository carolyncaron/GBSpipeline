#!/usr/bin/perl -w

=pod

=head1 NAME

GBS_Pipeline.pl - A complete set of commands to demultiplex, trim, align and call raw variants on paired-end reads for the purpose of genotyping by sequencing (GBS).

=head1 SYNOPSIS

./GBS_Pipeline.pl [function] [arg1] [arg2] ...

=head1 DESCRIPTION

The pipeline contains a set of steps (functions) that can be called individually to
complete a GBS analysis given multiplexed paired-end read data.

The GBS_Pipeline directory should be placed in the same directory containing all the
necessary files (ideally, create copies of these files and place them in a specific
directory for the GBS analysis).

=head2 Functions

=over 6

=item function1 sample_name F<barcode_file read1_file read2_file>

Demultiplex reads based on a barcode file (provided by Illumina to distinguish samples
used in sequencing).

sample_name can be any name assigned to the reads, to be used in naming output files.
Avoid use of whitespace (Ex: lens culinaris => lens_culinaris)

read1_file and read2_file should be provided in FASTQ format version Illumina 1.8+

=item function2 F<trimmomatic_path trim_file>

=back

=cut

use strict;
use warnings;
use File::Basename;
use IO::File;
use Getopt::Long;
use feature qw/switch/;
use Time::HiRes;
use FindBin qw($Bin);

# use Cwd 'abs_path';
# my $SCRIPT_PATH = abs_path($0);
# print "$SCRIPT_PATH\n";
# $SCRIPT_PATH =~ s/\/$0$/\//;
# print "$SCRIPT_PATH\n";

######################
##### USAGE MENU #####
######################

sub usage_menu
{
    print "For help, try:  perldoc ./GBS_Pipeline.pl\n";
#     use Term::ANSIColor qw(:constants);
#     $Term::ANSIColor::AUTORESET = 1;
#
#     print BOLD "\nGBS_pipeline.pl: ";
#     print "A GBS pipeline which aligns raw paired read data to a reference genome \n",
#         "and performs raw variant calls.\n\n";
#     print BOLD "USAGE:\n";
#     print "Perl GBS_pipeline.pl [function] [arg1] [arg2] ... \n\n";
#     print BOLD "FILES NEEDED:\n";
#     print "indices.list\tA list of indices that were used in sequencing to uniquely tag reads.\n",
#         "\n";
#     print BOLD "FUNCTIONS:\n";
#     print "(Avoid use of whitespace within parameters. Ex: Lens culinaris => Lens_culinaris)\n";
#     print "1. function1 sample_name barcode_file input1 input2\n",
#         "2.\n";
#     print "\n";
#     exit();
}

######################
#### CONFIG FILE #####
######################
my $CONFIG_FILE = "GBS.conf";

sub add_to_config
{
    open CONFIG, ">>$CONFIG_FILE" or die "ERROR: Unable to open config file $CONFIG_FILE\n";
    my $variable_name = $_[0];
    my $value = $_[1];
    my $comment = $_[2];
# TODO #
    # What if the variable exists but the value has changed??
    unless ( grep { /^\$variable_name/ } (<CONFIG>) )
    {
        print CONFIG "#$comment\n$variable_name = $value\n\n";
    }
    close CONFIG or die "ERROR: Unable to close config file $CONFIG_FILE\n";
}

################
##### MAIN #####
################

# Check for user-specified functions to run- if none specified, display usage menu #

########### TODO ##########
# Check if command line flags were provided: this overrides any non-flag options
# Getopt::Long::Configure ('bundling');
#
# my ( $SAMPLE_NAME, $R1, $R2, $TRIM_DIR );
# GetOptions ('s|sample=s' => \$SAMPLE_NAME,
#             'r1=s' => \$R1, 'r2=s' => \$R2,
#             't|trimdir=s' => \$TRIM_DIR
#             );
###########################

# Check that at least one argument is given
if ( exists ( $ARGV[0] ) )
{
    my $FUNCTION = $ARGV[0];
    shift @ARGV;
    my @args = @ARGV;
    my $num_args = $#args + 1;

    print "User specified function $FUNCTION with $num_args args: @args\n";

    # Use a switch statement for accessing functions
    given($FUNCTION)
    {
        when ( /f1/ || /demultiplex/ )
        {
            unless ($num_args == 4)
            {
                print "ERROR: Unexpected number of parameters given ($num_args). Program will exit.\n";
                die "--Try: Perl GBS_pipeline.pl function1 sample_name index_file /path/to/file1/filename1.fastq /path/to/file2/filename2.fastq\n";
            }

            my $sample = $args[0];
            my $index_file = $args[1];
            my $R1_file = $args[2];
            my $R2_file = $args[3];

            #my $start = Time::HiRes::gettimeofday();

            print "Calling $FUNCTION ...\n";

            require "$Bin/GBS_function1.pl";
            function1($sample, $index_file, $R1_file, $R2_file);

            # Save the sample and index file into the configuration file
            if (-e $CONFIG_FILE) {   system("rm $CONFIG_FILE");  } #Start a fresh config file if one exists already
            add_to_config("SAMPLE",$sample,"The generic sample name used in naming files during processing");
            add_to_config("INDEX_FILE",$index_file,"The filename of the list of indices (aka barcodes) provided by Illumina");

            #summarize($FUNCTION, $start);
        }
        when ( /f2/ || /trim_reads/ )
        {
            unless ($num_args >= 2)
            {
                print "ERROR: Unexpected number of parameters given ($num_args). Program will exit.\n";
                die "--Try: ./GBS_pipeline.pl function2 /path/to/trimmomatic path/to/trim_file.fasta\n";
            }

            my $trimmomatic_path = $args[0];
            my $trim_file = $args[1];

            # Extract sample name and indices from the config file
            chomp(my $sample = `grep 'SAMPLE' $CONFIG_FILE | cut -d' ' -f3`);
            chomp(my $index_file = `grep 'INDEX_FILE' $CONFIG_FILE | cut -d' ' -f3`);
            print "index file = $index_file\n";

            print "Calling $FUNCTION ...\n";
            require "$Bin/GBS_function2.pl";
            f2($trimmomatic_path, $trim_file, $sample, $index_file);
        }
        default
        {
            print "ERROR: $FUNCTION does not exist.\n";
            usage_menu();
        }
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
    print "The following files were created:\n";

    print "Continue onto the next step? (Yes/No)\n",
        "Hint: You can read your summary file then call the function for the next step,\n",
        "by typing: Perl ./GBS_pipeline.pl $function [arg1] [arg2] ... \n";

}


##########################
##### USER FUNCTIONS #####
##########################

