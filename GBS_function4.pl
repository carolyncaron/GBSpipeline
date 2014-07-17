#!/usr/bin/perl -w

##### STEP 4 : SAM file processing and SNP calling #####
##### Usage: f4 samtools_dir
##### Required input:
#####   samtools_dir - location of user's copy of samtools
##### Output:
#####   align/$index_$sample_mapped.sam
#####   align/$index_$sample_mapped.sorted.sam
#####   align/$index_$sample_mapped.bam
#####   align/$index_$sample_mapped.sorted.bam
#####   variants/$index_$sample_mapped.bcf
#####   variants/$index_$sample_mapped.vcf

use strict;
use warnings;

sub f4
{
    # Collect function-specific parameters
    my $samtools_dir = $_[0];
    my $sample = $_[1];
    my $index_file = $_[2];
    my $reference_genome = $_[3];

    # Check for samtools/bcftools
    unless ( -s "$samtools_dir/samtools")
    {   die "ERROR: Can't locate samtools in $samtools_dir\n";  }
    unless ( -s "$samtools_dir/bcftools")
    {   die "ERROR: Can't locate bcftools in $samtools_dir\n";  }

    unless ( -f $index_file && -r $index_file )
    {   die "ERROR: $index_file does not exist or is unreadable.\n";    }

    # Check for the reference genome
    unless (-f $reference_genome && -r $reference_genome)
    {   die "ERROR: $reference_genome does not exist or is unreadable.\n";  }

    system("mkdir -p variants/");

    my @indices = `cat $index_file`;
    foreach my $index (@indices)
    {
        chomp($index);
        my $sam_file = "align/$index\_$sample.sam";

        # 4a. Filter samfiles of sequences which did not align to save space
        my $cmd = "$samtools_dir/samtools view -ShF 4 -T $reference_genome $sam_file > align/$index\_$sample\_mapped.sam";
        my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
            run( command => $cmd, verbose => 0 );
        unless ($success)
        {   die "ERROR: Could not filter $sam_file for unmapped reads:\n$error_message\n@$stderr_buf";   }

        # 4b. Sort sam files so they can be indexed when converted to BAM
        $cmd = "sort -k3,3 -k4,4n align/$index\_$sample\_mapped.sam > align/$index\_$sample\_mapped.sorted.sam";
        ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
            run( command => $cmd, verbose => 0 );
        unless ($success)
        {   die "ERROR: Could not sort align/$index\_$sample\_mapped.sam: $error_message\n@$stderr_buf";   }

        # 4c. Convert to BAM format (Creates .fai files from the reference genome then BAM)
        $cmd = "$samtools_dir/samtools view -bT $reference_genome align/$index\_$sample\_mapped.sorted.sam ";
        $cmd .= "> align/$index\_$sample\_mapped.bam";
        ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
            run( command => $cmd, verbose => 0 );
        unless ($success)
        {   die "ERROR: Failed to convert $index\_$sample\_mapped.sorted.sam to BAM format: $error_message\n@$stderr_buf";   }

        print "Index $index files converted from .sam to .bam\n";

        # 4c. Sort BAM
        $cmd = "$samtools_dir/samtools sort align/$index\_$sample\_mapped.bam align/$index\_$sample\_mapped.sorted";
        ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
            run( command => $cmd, verbose => 0 );
        unless ($success)
        {   die "ERROR: Failed to sort $index\_$sample\_mapped.bam: $error_message\n@$stderr_buf";   }

        # 4d. Index BAM
        $cmd = "$samtools_dir/samtools index align/$index\_$sample\_mapped.sorted.bam";
        ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
            run( command => $cmd, verbose => 0 );
        unless ($success)
        {   die "ERROR: Failed to index $index\_$sample\_mapped.sorted.bam: $error_message\n@$stderr_buf";   }

        #4e. Identify genomic variants using mpileup
        $cmd = "$samtools_dir/samtools mpileup -f $reference_genome -g -I -B -D ";
        $cmd .= "align/$index\_$sample\_mapped.sorted.bam > variants/$index\_$sample\_mapped.bcf";
        ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
            run( command => $cmd, verbose => 0 );
        unless ($success)
        {   die "ERROR: Failed to identify genomic variants with mpileup: $error_message\n@$stderr_buf";   }

        $cmd = "bcftools view -c -g -I -v variants/$index\_$sample\_mapped.bcf > variants/$index\_$sample\_mapped.vcf";
        ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
            run( command => $cmd, verbose => 0 );
        unless ($success)
        {   die "ERROR: Failed to call SNPs/indels using bcftools:\n$error_message\n@$stderr_buf";   }

        print "Completed variant calls for index $index\n";

    }
}

1;
