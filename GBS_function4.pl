#!/usr/bin/perl -w

##### STEP 4 : SAM file processing and SNP calling #####
##### Usage: f4 samtools_dir bcftools_dir
##### Required input:
#####   samtools_dir - location of user's copy of samtools
#####   bcftools_dir - location of user's copy of bcftools
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

    # Check samtools + bcftools versions
    my $sam_version = `$samtools_dir/samtools --version-only` || die "ERROR: Ensure that samtools version is 1.0 or greater.\n";
    my $bcf_version = `$bcftools_dir/bcftools --version-only` || die "ERROR: Ensure that bcftools version is 1.0 or greater.\n";

    unless ( -f $index_file && -r $index_file )
    {   die "ERROR: $index_file does not exist or is unreadable.\n";    }

    # Check for the reference genome
    unless (-f $reference_genome && -r $reference_genome)
    {   die "ERROR: $reference_genome does not exist or is unreadable.\n";  }

    system("mkdir -p $output_dir/variants/");

    my @indices = `cat $index_file`;
    my $num_indices = $#indices + 1;
    if ($num_indices == 0){    die "ERROR: $index_file exists but appears to be empty.\n";    }
    my $index_count = 0;
    print_progress($index_count, $num_indices);

    # BEGIN SNP calling
    foreach my $index (@indices)
    {
        # Seven steps per index... break it down so progress bar reports more often
        my $num_steps = $index_count*8;
        chomp($index);
        my $sam_file = "$output_dir/align/$index\_$sample.sam";

        # Multi-mapping filtering
        bwt2_besthits($sam_file, "$output_dir/variants/$index\_$sample\_besthits.sam");
        $num_steps++;
        print_progress($num_steps, $num_indices*8, " Filtering multi-mapped reads");

        # 4a. Filter samfiles of sequences which did not align to save space
        ### Version 1.0: -S not necessary, SAM format auto detected
        my $cmd = "$samtools_dir/samtools view -ShF 4 -T $reference_genome ";
        $cmd .= "$output_dir/variants/$index\_$sample\_besthits.sam > ";
        $cmd .= "$output_dir/align/$index\_$sample\_mapped.sam";
        my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) = run( command => $cmd, verbose => 0 );
        unless ($success)
        {   die "ERROR: Could not filter $output_dir/variants/$index\_$sample\_besthits.sam for unmapped reads:\n$error_message\n@$stderr_buf";   }
        $num_steps++;
        print_progress($num_steps, $num_indices*8, " Filtering non-aligned reads ");

        # 4b. Sort sam files so they can be indexed when converted to BAM
        $cmd = "sort -k3,3 -k4,4n $output_dir/align/$index\_$sample\_mapped.sam > ";
        $cmd .= "$output_dir/align/$index\_$sample\_mapped.sorted.sam";
        ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) = run( command => $cmd, verbose => 0 );
        unless ($success)
        {   die "ERROR: Could not sort $index\_$sample\_mapped.sam: $error_message\n@$stderr_buf";   }
        $num_steps++;
        print_progress($num_steps, $num_indices*8, " Sorting SAM                 ");

        # 4c. Convert to BAM format (Creates .fai files from the reference genome then BAM)
        $cmd = "$samtools_dir/samtools view -bT $reference_genome $output_dir/align/$index\_$sample\_mapped.sorted.sam ";
        $cmd .= "> $output_dir/align/$index\_$sample\_mapped.bam";
        ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) = run( command => $cmd, verbose => 0 );
        unless ($success)
        {   die "ERROR: Failed to convert $index\_$sample\_mapped.sorted.sam to BAM format: $error_message\n@$stderr_buf";   }
        $num_steps++;
        print_progress($num_steps, $num_indices*8, " Converting to BAM");

        # 4d. Sort BAM
        $cmd = "$samtools_dir/samtools sort $output_dir/align/$index\_$sample\_mapped.bam ";
        $cmd .= "$output_dir/align/$index\_$sample\_mapped.sorted";
        ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) = run( command => $cmd, verbose => 0 );
        unless ($success)
        {   die "ERROR: Failed to sort $index\_$sample\_mapped.bam: $error_message\n@$stderr_buf";   }
        $num_steps++;
        print_progress($num_steps, $num_indices*8, " Sorting BAM      ");

        # 4e. Index BAM
        $cmd = "$samtools_dir/samtools index $output_dir/align/$index\_$sample\_mapped.sorted.bam";
        ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) = run( command => $cmd, verbose => 0 );
        unless ($success)
        {   die "ERROR: Failed to index $index\_$sample\_mapped.sorted.bam: $error_message\n@$stderr_buf";   }
        $num_steps++;
        print_progress($num_steps, $num_indices*8, "Indexing BAM");

        # 4f. Identify genomic variants using mpileup
        ### Version 1.0: -D deprecated, use -t DP instead
        $cmd = "$samtools_dir/samtools mpileup -f $reference_genome -g -I -B -D ";
        $cmd .= "$output_dir/align/$index\_$sample\_mapped.sorted.bam > $output_dir/variants/$index\_$sample\_mapped.bcf";
        ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) = run( command => $cmd, verbose => 0 );
        unless ($success)
        {   die "ERROR: Failed to identify genomic variants with mpileup: $error_message\n@$stderr_buf";   }
        $num_steps++;
        print_progress($num_steps, $num_indices*8, " Identifying variants");

        # 4g. Use bcftools to call SNPs using bcf files
        ### -c -g -I -v
        $cmd = "$bcftools_dir/bcftools call -c -v $output_dir/variants/$index\_$sample\_mapped.bcf > ";
        $cmd .= "$output_dir/variants/$index\_$sample\_mapped.vcf";
        ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) = run( command => $cmd, verbose => 0 );
        unless ($success)
        {   die "ERROR: Failed to call SNPs/indels using bcftools:\n$error_message\n@$stderr_buf";   }
        $num_steps++;

        # Summary of progress
        $index_count++;
        print_progress($num_steps, $num_indices*8, "Completed index $index");
    }

    # Merge mpileup results into a single file


    print "\n Processed reads located in:\n  $output_dir/align/\n  $output_dir/variants/ \n";
}

# Filters bowtie2 results to choose the "best" alignments for reads which multimapped.
sub bwt2_besthits
{
    my $sam_file = $_[0];
    my $out_file = $_[1];

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
}

1;
