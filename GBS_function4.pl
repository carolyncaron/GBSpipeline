#!/usr/bin/perl -w

##### STEP 4 : SAM file processing and SNP calling #####
##### Usage: f4 samtools_dir bcftools_dir
##### Required input:
#####   samtools_dir - location of user's copy of samtools
##### Output:
#####   $output_dir/align/$index_$sample_mapped.sam
#####   $output_dir/align/$index_$sample_mapped.sorted.sam
#####   $output_dir/align/$index_$sample_mapped.bam
#####   $output_dir/align/$index_$sample_mapped.sorted.bam
#####   $output_dir/variants/$index_$sample_mapped.bcf
#####   $output_dir/variants/$index_$sample_mapped.vcf

use strict;
use warnings;

sub f4
{
    # Collect function-specific parameters
    my $samtools_dir = $_[0];
    my $bcftools_dir = $_[1];
    my $sample = $_[2];
    my $index_file = $_[3];
    my $output_dir = $_[4];
    my $reference_genome = $_[5];

    # Check for samtools/bcftools
    unless ( -s "$samtools_dir/samtools")
    {   die "ERROR: Can't locate samtools in $samtools_dir\n";  }
    unless ( -s "$bcftools_dir/bcftools")
    {   die "ERROR: Can't locate bcftools in $bcftools_dir\n";  }

    unless ( -f $index_file && -r $index_file )
    {   die "ERROR: $index_file does not exist or is unreadable.\n";    }

    # Check for the reference genome
    unless (-f $reference_genome && -r $reference_genome)
    {   die "ERROR: $reference_genome does not exist or is unreadable.\n";  }

    system("mkdir -p $output_dir/variants/");

    my @indices = `cat $index_file`;
    foreach my $index (@indices)
    {
        chomp($index);
        my $sam_file = "$output_dir/align/$index\_$sample.sam";

        # 4a. Filter samfiles of sequences which did not align to save space
        ### Version 1.0: -S not necessary, SAM format auto detected
        my $cmd = "$samtools_dir/samtools view -ShF 4 -T $reference_genome $sam_file > ";
        $cmd .= "$output_dir/align/$index\_$sample\_mapped.sam";
        my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) = run( command => $cmd, verbose => 0 );
        unless ($success)
        {   die "ERROR: Could not filter $sam_file for unmapped reads:\n$error_message\n@$stderr_buf";   }

        # 4b. Sort sam files so they can be indexed when converted to BAM
        $cmd = "sort -k3,3 -k4,4n $output_dir/align/$index\_$sample\_mapped.sam > ";
        $cmd .= "$output_dir/align/$index\_$sample\_mapped.sorted.sam";
        ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) = run( command => $cmd, verbose => 0 );
        unless ($success)
        {   die "ERROR: Could not sort $index\_$sample\_mapped.sam: $error_message\n@$stderr_buf";   }

        # 4c. Convert to BAM format (Creates .fai files from the reference genome then BAM)
        $cmd = "$samtools_dir/samtools view -bT $reference_genome $output_dir/align/$index\_$sample\_mapped.sorted.sam ";
        $cmd .= "> $output_dir/align/$index\_$sample\_mapped.bam";
        ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) = run( command => $cmd, verbose => 0 );
        unless ($success)
        {   die "ERROR: Failed to convert $index\_$sample\_mapped.sorted.sam to BAM format: $error_message\n@$stderr_buf";   }

        print "Index $index files converted from .sam to .bam\n";

        # 4c. Sort BAM
        $cmd = "$samtools_dir/samtools sort $output_dir/align/$index\_$sample\_mapped.bam ";
        $cmd .= "$output_dir/align/$index\_$sample\_mapped.sorted";
        ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) = run( command => $cmd, verbose => 0 );
        unless ($success)
        {   die "ERROR: Failed to sort $index\_$sample\_mapped.bam: $error_message\n@$stderr_buf";   }

        # 4d. Index BAM
        $cmd = "$samtools_dir/samtools index $output_dir/align/$index\_$sample\_mapped.sorted.bam";
        ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) = run( command => $cmd, verbose => 0 );
        unless ($success)
        {   die "ERROR: Failed to index $index\_$sample\_mapped.sorted.bam: $error_message\n@$stderr_buf";   }

        #4e. Identify genomic variants using mpileup
        ### Version 1.0: -D deprecated, use -t DP instead
        $cmd = "$samtools_dir/samtools mpileup -f $reference_genome -g -I -B -D ";
        $cmd .= "$output_dir/align/$index\_$sample\_mapped.sorted.bam > $output_dir/variants/$index\_$sample\_mapped.bcf";
        ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) = run( command => $cmd, verbose => 0 );
        unless ($success)
        {   die "ERROR: Failed to identify genomic variants with mpileup: $error_message\n@$stderr_buf";   }

        $cmd = "$bcftools_dir/bcftools view -c -g -I -v $output_dir/variants/$index\_$sample\_mapped.bcf > ";
        $cmd .= "$output_dir/variants/$index\_$sample\_mapped.vcf";
        ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) = run( command => $cmd, verbose => 0 );
        unless ($success)
        {   die "ERROR: Failed to call SNPs/indels using bcftools:\n$error_message\n@$stderr_buf";   }

        print "Completed variant calls for index $index\n";
    }
}

1;
