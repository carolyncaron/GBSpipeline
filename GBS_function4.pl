#!/usr/bin/perl -w

##### STEP 4 : SAM file processing and SNP calling #####
##### Usage: call_SNPs [SAMTOOLS_PATH] [BCFTOOLS_PATH]
##### Required input:
#####   SAMTOOLS_PATH : location of user's copy of samtools
#####   BCFTOOLS_PATH : location of user's copy of bcftools
##### Output:
#####   [OUTPUT_DIR]/align/[SAMPLE]_[POPULATION]_mapped.bam
#####   [OUTPUT_DIR]/align/[SAMPLE]_[POPULATION]_mapped.sorted.bam
#####   [OUTPUT_DIR]/variants/[SAMPLE]_[POPULATION].bcf and/or [OUTPUT_DIR]/variants/[POPULATION].bcf
#####   [OUTPUT_DIR]/variants/[SAMPLE]_[POPULATION].vcf and/or [OUTPUT_DIR]/variants/[POPULATION].vcf

use strict;
use warnings;

sub f4
{
    # Collect function-specific parameters
    my $samtools_dir = $_[0];
    my $bcftools_dir = $_[1];
    my $call_mode = $_[2];
    my $population = $_[3];
    my $index_file = $_[4];
    my $output_dir = $_[5];
    my $reference_genome = $_[6];
    my $samtools_threads = $_[7];

    # Check for samtools/bcftools
    unless ( -s "$samtools_dir/samtools")
    {   die "ERROR: Can't locate samtools in $samtools_dir\n";  }
    unless ( -s "$bcftools_dir/bcftools")
    {   die "ERROR: Can't locate bcftools in $bcftools_dir\n";  }

    # Check samtools + bcftools versions
    my $sam_version = `$samtools_dir/samtools --version-only` || die "ERROR: Please ensure that samtools version is 1.0 or greater.\n";
    my $bcf_version = `$bcftools_dir/bcftools --version-only` || die "ERROR: Please ensure that bcftools version is 1.0 or greater.\n";

    # Check that call mode is valid
    unless ($call_mode =~ /single|multi|both/i)
    {   die "ERROR: \"$call_mode\" is not a valid mode for calling SNPs. Please fix this in the configuration file.\n";  }

    unless ( -f $index_file && -r $index_file )
    {   die "ERROR: $index_file does not exist or is unreadable.\n";    }

    my @samples = `cut -f1 $index_file | sort -u`;
    my $num_samples = $#samples + 1;
    if ($num_samples == 0){    die "ERROR: $index_file exists but appears to be empty.\n";    }
    my ($sample_count, $step_count, $total_steps) = 0;
    print_progress($sample_count, $num_samples, " Building genome index");

    # Check for the reference genome
    unless (-f $reference_genome && -r $reference_genome)
    {   die "ERROR: $reference_genome does not exist or is unreadable.\n";  }

    # Create the index from the reference genome for samtools to use
    my $cmd = "$samtools_dir/samtools faidx $reference_genome";
    my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) = run( command => $cmd, verbose => 0 );
        unless ($success)
    {   die "ERROR: Failed to index the reference genome with samtools faidx: $error_message\n@$stderr_buf";   }

    system("mkdir -p $output_dir/variants/");

    # These are for multi-mode
    my $mpileup_input_files = "";
    my $sample_list = "samples.tmp";

    # BEGIN SNP calling
    foreach my $sample (@samples)
    {
        # Number of steps per sample depends on single or multi SNP calling mode
        if ( $call_mode =~ /single/i ) {
            $step_count = $sample_count*5;
            $total_steps = $num_samples*5;
        } elsif ( $call_mode =~ /multi/i ) {
            $step_count = $sample_count*3;
            $total_steps = ($num_samples*3)+2;
        } else { # both
            $step_count = $sample_count*5;
            $total_steps = ($num_samples*5)+2;
        }
        chomp($sample);
        $sample =~ s/ //g;

        # Filenames
        my $sam_file = "$output_dir/align/$sample\_$population\_besthits.sam";
        my $bam_file = "$output_dir/align/$sample\_$population\_mapped.bam";
        my $mpileup_input_file = "$output_dir/align/$sample\_$population\_mapped.sorted.bam";
        my $bcf_file = "$output_dir/variants/$sample\_$population\.bcf";
        my $vcf_file = "$output_dir/variants/$sample\_$population\.vcf";

        ## 1. Filter out sequences which did not align, convert to BAM format
        print_progress($step_count++, $total_steps, " $sample: Converting to BAM ");
        $cmd = "$samtools_dir/samtools view -b -@ $samtools_threads -F 4 -T $reference_genome ";
        $cmd .= "$sam_file > $bam_file";
        ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) = run( command => $cmd, verbose => 0 );
        unless ($success)
        {   die "ERROR: Failed to convert $output_dir/align/$sample\_$population\_besthits.sam to BAM: \
            $error_message\n@$stderr_buf";   }

        ## 2. Sort BAM
        print_progress($step_count++, $total_steps, " $sample: Sorting BAM        ");
        $cmd = "$samtools_dir/samtools sort -@ $samtools_threads -O bam -T $bam_file\_tmp ";
        $cmd .= "$bam_file > $mpileup_input_file";
        ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) = run( command => $cmd, verbose => 0 );
        unless ($success)
        {   die "ERROR: Failed to sort $sample\_$population\_mapped.bam: $error_message\n@$stderr_buf";   }

        ## 3. Index BAM
        print_progress($step_count++, $total_steps, " $sample: Indexing BAM ");
        $cmd = "$samtools_dir/samtools index $mpileup_input_file";
        ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) = run( command => $cmd, verbose => 0 );
        unless ($success)
        {   die "ERROR: Failed to index $mpileup_input_file: $error_message\n@$stderr_buf";   }

        if ($call_mode =~ /single|both/i)
        {
            ## 4. Identify genomic variants using mpileup
            ### Since version 1.0: -D deprecated, use -t DP instead
            print_progress($step_count++, $total_steps, " $sample: Identifying variants");
            $cmd = "$samtools_dir/samtools mpileup -f $reference_genome -g -I -B -t DP,AD ";
            $cmd .= "$mpileup_input_file > $bcf_file";
            ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) = run( command => $cmd, verbose => 0 );
            unless ($success)
            {   die "ERROR: Failed to identify genomic variants with mpileup: $error_message\n@$stderr_buf";   }

            ## 5. Use bcftools to call SNPs using bcf files
            ### Parameters prior to 1.0: -c -g -I -v
            print_progress($step_count++, $total_steps, " $sample: Calling SNPs        ");
            if ($bcf_version =~ /1.0|1.1|1.2/)
            {
                $cmd = "$bcftools_dir/bcftools call -c -v -V indels -A $bcf_file > $vcf_file";
            } else
            {
                $cmd = "$bcftools_dir/bcftools call -m -v -V indels -A --threads $samtools_threads $bcf_file > $vcf_file";
            }
            #$cmd = "$bcftools_dir/bcftools call -c -v -V indels -A $bcf_file > $vcf_file";
            ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) = run( command => $cmd, verbose => 0 );
            unless ($success)
            {   die "ERROR: Failed to call SNPs/indels using bcftools:\n$error_message\n@$stderr_buf";   }
        }
        if ($call_mode =~ /multi|both/i)
        {
            # If we are using mpileup on all samples, add the input file to a list
            $mpileup_input_files .= "$mpileup_input_file ";
            open SAMPLES, ">>$sample_list" or die "ERROR: Couldn't open temporary file $sample_list\n";
            print SAMPLES "$mpileup_input_file\t$sample\n";
            close SAMPLES or die "ERROR: Couldn't close temporary file $sample_list\n";
        }

        print_progress($step_count++, $total_steps, "");
        $sample_count++;
        #summarize_SNPcall($sample, $population);
    }

    if ($call_mode =~ /multi|both/i)
    {
        #$step_count = $total_steps - 2;
        my $bcf_file = "$output_dir/variants/$population\.bcf";
        my $rehead_bcf = "$output_dir/variants/$population\.rehead.bcf";
        my $vcf_file = "$output_dir/variants/$population\.vcf";

        # Call mpileup on ALL samples at once, then call SNPs
        print_progress($step_count++, $total_steps, " Identifying variants for all samples ");
        my $cmd = "$samtools_dir/samtools mpileup -f $reference_genome -g -I -B -t DP ";
        $cmd .= "$mpileup_input_files> $bcf_file";
        my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) = run( command => $cmd, verbose => 0 );
        unless ($success)
        {   die "ERROR: Failed to identify genomic variants with mpileup: $error_message\n@$stderr_buf";   }

        # Change the header in the BCF file to be the sample names rather than the full path to each input file
        $cmd = "$bcftools_dir/bcftools reheader -s $sample_list $bcf_file > $rehead_bcf";
        ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) = run( command => $cmd, verbose => 0 );
        unless ($success)
        {   die "ERROR: Failed to change the header in $bcf_file: $error_message\n@$stderr_buf";   }

        # Depending on our version of bcftools, use the old caller (v1.0-1.2) or the new caller (v1.3+)
        print_progress($step_count++, $total_steps, " Calling SNPs for all samples        ");
        if ($bcf_version =~ /1.0|1.1|1.2/)
        {
            $cmd = "$bcftools_dir/bcftools call -c -v -V indels -A $rehead_bcf > $vcf_file";
        } else
        {
            $cmd = "$bcftools_dir/bcftools call -m -v -V indels -A $rehead_bcf --threads $samtools_threads > $vcf_file";
        }
        ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) = run( command => $cmd, verbose => 0 );
        unless ($success)
        {   die "ERROR: Failed to call SNPs/indels using bcftools:\n$error_message\n@$stderr_buf";   }

        system("rm $sample_list");
    }

    print "\n Processed alignment and variant call format files located in:\n  $output_dir/align/\n  $output_dir/variants/ \n";
}


##############################
##### Summarize
##### sample     # of unique hits
##############################
sub summarize_SNPcall
{
     my $sample = $_[0];
     my $population = $_[1];

}


1;
