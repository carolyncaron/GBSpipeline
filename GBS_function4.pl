#!/usr/bin/perl -w

##### STEP 4 : SAM file processing and SNP calling #####
##### Usage: f4 samtools_dir bcftools_dir
##### Required input:
#####   samtools_dir - location of user's copy of samtools
#####   bcftools_dir - location of user's copy of bcftools
##### Output:
#####   $output_dir/align/$sample_$population_mapped.sam
#####   $output_dir/align/$sample_$population_mapped.sorted.sam
#####   $output_dir/align/$sample_$population_mapped.bam
#####   $output_dir/align/$sample_$population_mapped.sorted.bam
#####   $output_dir/variants/$sample_$population_mapped.bcf
#####   $output_dir/variants/$sample_$population_mapped.vcf

use strict;
use warnings;

sub f4
{
    # Collect function-specific parameters
    my $samtools_dir = $_[0];
    my $bcftools_dir = $_[1];
    my $population = $_[2];
    my $index_file = $_[3];
    my $output_dir = $_[4];
    my $reference_genome = $_[5];

    # Check for samtools/bcftools
    unless ( -s "$samtools_dir/samtools")
    {   die "ERROR: Can't locate samtools in $samtools_dir\n";  }
    unless ( -s "$bcftools_dir/bcftools")
    {   die "ERROR: Can't locate bcftools in $bcftools_dir\n";  }

    # Check samtools + bcftools versions
    my $sam_version = `$samtools_dir/samtools --version-only` || die "ERROR: Please ensure that samtools version is 1.0 or greater.\n";
    my $bcf_version = `$bcftools_dir/bcftools --version-only` || die "ERROR: Please ensure that bcftools version is 1.0 or greater.\n";

    unless ( -f $index_file && -r $index_file )
    {   die "ERROR: $index_file does not exist or is unreadable.\n";    }

    # Check for the reference genome
    unless (-f $reference_genome && -r $reference_genome)
    {   die "ERROR: $reference_genome does not exist or is unreadable.\n";  }

    system("mkdir -p $output_dir/variants/");

    #print "Default action is to combine \n",
    #          "Would you like to attempt to install v$version there now? (yes/no) ";
    #chomp (my $usr_input = <STDIN>);
    #if ($usr_input =~ /yes/i)
    #{}

    my @samples = `cut -f1 $index_file | sort -u`;
    my $num_samples = $#samples + 1;
    if ($num_samples == 0){    die "ERROR: $index_file exists but appears to be empty.\n";    }
    my $sample_count = 0;
    print_progress($sample_count, $num_samples, " Filtering multi-mapped reads");

    # BEGIN SNP calling
    foreach my $sample (@samples)
    {
        # Eight steps per sample... break it down so progress bar reports more often
        my $step_count = $sample_count*8;
        my $total_steps = $num_samples*8;
        chomp($sample);
        $sample =~ s/ //g;
        my $sam_file = "$output_dir/align/$sample\_$population.sam";

        ## 1. Filter alignments in which reads multi-mapped for the best hits
        print_progress($step_count++, $total_steps, " Filtering multi-mapped reads   ");
        my $besthits_log = "logs/$sample\_$population\_multimap_processing.log";
        bwt2_besthits($sam_file, "$output_dir/align/$sample\_$population\_besthits.sam", $besthits_log);

        ## 2. Filter out sequences which did not align to save space
        print_progress($step_count++, $total_steps, " Filtering non-aligned reads   ");
        my $cmd = "$samtools_dir/samtools view -F 4 -T $reference_genome ";
        $cmd .= "$output_dir/align/$sample\_$population\_besthits.sam > ";
        $cmd .= "$output_dir/align/$sample\_$population\_mapped.sam";
        my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) = run( command => $cmd, verbose => 0 );
        unless ($success)
        {   die "ERROR: Could not filter $output_dir/align/$sample\_$population\_besthits.sam \
                 for unmapped reads:\n$error_message\n@$stderr_buf";   }

        ## 3. Sort sam files so they can be indexed when converted to BAM
        print_progress($step_count++, $total_steps, " Sorting SAM                 ");
        $cmd = "sort -k3,3 -k4,4n $output_dir/align/$sample\_$population\_mapped.sam > ";
        $cmd .= "$output_dir/align/$sample\_$population\_mapped.sorted.sam";
        ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) = run( command => $cmd, verbose => 0 );
        unless ($success)
        {   die "ERROR: Could not sort $sample\_$population\_mapped.sam: $error_message\n@$stderr_buf";   }

        ## 4. Convert to BAM format (Creates .fai files from the reference genome then BAM)
        print_progress($step_count++, $total_steps, " Converting to BAM ");
        $cmd = "$samtools_dir/samtools view -bT $reference_genome ";
        $cmd .= "$output_dir/align/$sample\_$population\_mapped.sorted.sam ";
        $cmd .= "> $output_dir/align/$sample\_$population\_mapped.bam";
        ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) = run( command => $cmd, verbose => 0 );
        unless ($success)
        {   die "ERROR: Failed to convert $sample\_$population\_mapped.sorted.sam to BAM format: \
                 $error_message\n@$stderr_buf";   }

        ## 5. Sort BAM
        print_progress($step_count++, $total_steps, " Sorting BAM        ");
        $cmd = "$samtools_dir/samtools sort $output_dir/align/$sample\_$population\_mapped.bam ";
        $cmd .= "$output_dir/align/$sample\_$population\_mapped.sorted";
        ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) = run( command => $cmd, verbose => 0 );
        unless ($success)
        {   die "ERROR: Failed to sort $sample\_$population\_mapped.bam: $error_message\n@$stderr_buf";   }

        ## 6. Index BAM
        print_progress($step_count++, $total_steps, "Indexing BAM");
        $cmd = "$samtools_dir/samtools index $output_dir/align/$sample\_$population\_mapped.sorted.bam";
        ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) = run( command => $cmd, verbose => 0 );
        unless ($success)
        {   die "ERROR: Failed to index $sample\_$population\_mapped.sorted.bam: $error_message\n@$stderr_buf";   }

        ### @TODO: Give the user the option to call mpileup on all samples rather than
        ###        each one individually
        ## 7. Identify genomic variants using mpileup
        ### Since version 1.0: -D deprecated, use -t DP instead
        print_progress($step_count++, $total_steps, " Identifying variants");
        $cmd = "$samtools_dir/samtools mpileup -f $reference_genome -g -I -B -t DP ";
        $cmd .= "$output_dir/align/$sample\_$population\_mapped.sorted.bam > ";
        $cmd .= "$output_dir/variants/$sample\_$population\_mapped.bcf";
        ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) = run( command => $cmd, verbose => 0 );
        unless ($success)
        {   die "ERROR: Failed to identify genomic variants with mpileup: $error_message\n@$stderr_buf";   }

        ## 8. Use bcftools to call SNPs using bcf files
        ### Parameters prior to 1.0: -c -g -I -v
        print_progress($step_count++, $total_steps, " Calling SNPs        ");
        $cmd = "$bcftools_dir/bcftools call -c -v -V indels -A $output_dir/variants/$sample\_$population\_mapped.bcf > ";
        $cmd .= "$output_dir/variants/$sample\_$population\_mapped.vcf";
        ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) = run( command => $cmd, verbose => 0 );
        unless ($success)
        {   die "ERROR: Failed to call SNPs/indels using bcftools:\n$error_message\n@$stderr_buf";   }

        # Summary of progress
        $sample_count++;
        print_progress($step_count++, $total_steps, "                ");
        #summarize_SNPcall($sample, $population);
    }

    print "\n Processed reads located in:\n  $output_dir/align/\n  $output_dir/variants/ \n";
}

# Filters bowtie2 results to choose the "best" alignments for reads which multimapped.
sub bwt2_besthits
{
    my $sam_file = $_[0];
    my $out_file = $_[1];
    my $hitslog = $_[2];

    open(SAM,"<".$sam_file) || die "ERROR: Could not open $sam_file\n";
    open(OUT,">".$out_file) || die "ERROR: Could not write out to $out_file\n";

    my %data;
    # Don't know the significance of this variable
    my $cons = 18;
    while(<SAM>) {
        chomp;
        if ($_ =~ /^@/)
        {   print OUT $_."\n";  }
        else
        {   my ($id1, $bw_flag1, $s1, $s_s1, $mapq1, $maps1, $ja1, $p_s1, $size1, $seq1, $jb1, $aln_score1, @junks1) = (split("\t",$_));
            $aln_score1=~s/AS:i://g;
            next if $s1 eq '*';
            my $l1 = $_."\n";
            my $l2 = <SAM>;

            my ($id2, $bw_flag2, $s2, $s_s2, $mapq2, $maps2, $ja2, $p_s2, $size2, $seq2, $jb2, $aln_score2, @junks2) = (split("\t",$l2));
            $aln_score2=~s/AS:i://g;

            unless ($id1 eq $id2 && $s1 eq $s2 && $s_s1 eq $p_s2 && $p_s1 eq $s_s2)
            {
                die "ERROR: Unexpected file format\n";
            }

            if (exists($data{$id1}))
            {
                if( $data{$id1}->{score} <= ($aln_score1 + $aln_score2+$cons))
                {
                    $data{$id1}->{num_best}+=1;
                }
            }
            else
            {
                $data{$id1} = {l1=>$l1, l2=>$l2,num_best=>1,score=>$aln_score1+$aln_score2};
            }
        }
    }
    my $num_uniq_hits = 0;
    my $num_multi_hits = 0;

    foreach my $id (keys %data)
    {
        if ($data{$id}->{num_best} eq 1) {
            $num_uniq_hits++;
            print OUT $data{$id}->{l1}.$data{$id}->{l2};
        }
        else {   $num_multi_hits++; }
    }

    open HITSLOG, ">$hitslog";
    print HITSLOG "Found $num_uniq_hits number of uniquely mapped reads.\n";
    print HITSLOG "Found $num_multi_hits number of multi-mapped reads.\n";
    print HITSLOG "There is a total of ".($num_uniq_hits+$num_multi_hits)." hits.\n";

}

1;
