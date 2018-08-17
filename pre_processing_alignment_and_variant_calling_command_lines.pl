=head
DESCRIPTION: taking as input one or more directories of paired-end fastq files (one directory per sample, of format, e.g. "sample1/sample1.1.fq.gz" and "sample1/sample1.2.fq.gz"), this Perl script will output a shell script containing all commands necessary to pre-process these reads, align them to a reference genome, call variants per sample and, finally, to merge all variant call sets into one gVCF (genomic variant call file).

DISCLAIMER: the following should be considered "lab-ready" code and is provided for reproducibility purposes only. No attempt has been made to optimise the code or to test it upon other datasets. The software versions used are detailed below; note the associated command lines are subject to change should the software be upgraded.

PREREQUISITES: prior to use, the following prequisites must be met:
1. Two tools must be accessible from the command line: (a) bwa v0.7.8 (http://bio-bwa.sourceforge.net/), and (b) samtools v1.3 (http://samtools.sourceforge.net/).
2. Paths to the JAR files for Trimmomatic v0.35 (http://www.usadellab.org/cms/?page=trimmomatic), GATK v3.5 (https://software.broadinstitute.org/gatk/download/), and Picard Tools v2.1.1 (https://broadinstitute.github.io/picard/) must be hard-coded (see the section of script marked 'requirements', below).
3. Both an indexed reference genome (Amel_4.5) and a set of known SNPs (from dbSNP build 140) must be provided for use with GATK. These paths are hard-coded as the variables $ref_genome and $known_bee_snps, respectively (see the section of script marked 'requirements', below).

## obtain the A. mellifera reference genome. This must then be indexed both with SAMtools faidx (necessary for compatibility with GATK) and BWA (necessary for alignment). Also necessary for compatibility with GATK is a 'dictionary' file, created using Picard Tools. Note the output files from the indexing steps must remain in the same directory, and have the same root name, as the reference genome .fa file.
mkdir am_genome
cd am_genome
curl -O -L ftp://ftp.ensemblgenomes.org/pub/metazoa/release-30/fasta/apis_mellifera/dna/Apis_mellifera.GCA_000002195.1.30.dna.toplevel.fa.gz
gunzip Apis_mellifera.GCA_000002195.1.30.dna.toplevel.fa.gz
bwa index Apis_mellifera.GCA_000002195.1.30.dna.toplevel.fa
samtools faidx Apis_mellifera.GCA_000002195.1.30.dna.toplevel.fa
java -jar picard.jar CreateSequenceDictionary REFERENCE=Apis_mellifera.GCA_000002195.1.30.dna.toplevel.fa OUTPUT=Apis_mellifera.GCA_000002195.1.30.dna.toplevel.dict

## obtain the set of known bee SNPs from NCBI dbSNP. These comprise a set of VCFs (one per chromosome/contig/scaffold), which we shall merge into a single file - "bee_dbsnp_snps.vcf" - using GATK CombineVariants
mkdir bee_snps
cd bee_snps
wget -r --no-parent ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/archive/bee_7460/VCF/
cd ftp.ncbi.nlm.nih.gov/snp/organisms/archive/bee_7460/VCF/
gunzip *.gz # this error will happen if you don't unzip these files: http://gatkforums.broadinstitute.org/gatk/discussion/5349/genotypegvcfs-warn-track-variant-doesnt-have-a-sequence-dictionary-built-in
# at this point, we must pay attention to what the chromosome/contig names in each VCF are considered to be - they are 'LG1', 'LG2' etc., inconsistent with those given in the reference genome .fa (ignoring the many contig/scaffold IDs in the .fa, the chromosome names here are purely numeric: 1, 2, 3, etc. To see this for yourself, use "grep '>', am_genome/Apis_mellifera.GCA_000002195.1.30.dna.toplevel.fa"). This discrepancy can be rectified manually by renaming the chromosome IDs in each VCF, 'LG1' to '1', and so on. Delete the now-obsolete indices for each of the original VCFs before proceeding (GATK will automatically create an index should one be missing anyway).
java -jar GenomeAnalysisTK.jar -T CombineVariants -R am_genome/Apis_mellifera.GCA_000002195.1.30.dna.toplevel.fa -o bee_dbsnp_snps.vcf --variant vcf_chr_LG1.vcf --variant vcf_chr_LG2.vcf --variant vcf_chr_LG3.vcf --variant vcf_chr_LG4.vcf --variant vcf_chr_LG5.vcf --variant vcf_chr_LG6.vcf --variant vcf_chr_LG7.vcf --variant vcf_chr_LG8.vcf --variant vcf_chr_LG9.vcf --variant vcf_chr_LG10.vcf --variant vcf_chr_LG11.vcf --variant vcf_chr_LG12.vcf --variant vcf_chr_LG13.vcf --variant vcf_chr_LG14.vcf --variant vcf_chr_LG15.vcf --variant vcf_chr_LG16.vcf
=cut

use strict;
use warnings;

# REQUIREMENTS
my $gatk_path 		 = 'GenomeAnalysisTK.jar'; # available from https://software.broadinstitute.org/gatk/download/auth?package=GATK-archive&version=3.5-0-g36282e4
my $picard_path		 = 'picard.jar'; # available from https://github.com/broadinstitute/picard/releases/download/2.1.1/picard-tools-2.1.1.zip
my $trimmomatic_path = 'trimmomatic-0.35.jar'; # available from http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.35.zip
my $ref_genome		 = 'am_genome/Apis_mellifera.GCA_000002195.1.30.dna.toplevel.fa'; # available from ftp://ftp.ensemblgenomes.org/pub/metazoa/release-30/fasta/apis_mellifera/dna/Apis_mellifera.GCA_000002195.1.30.dna.toplevel.fa.gz
my $known_bee_snps	 = 'bee_dbsnp_snps.vcf'; # available by combining the content of ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/archive/bee_7460/VCF/ (see above)
my $raw_reads_dir 	 = 'raw_reads'; # this directory should contain numerous subdirectories, one per sample, that each contain paired-end .fq.gz files

# OUTPUT
my $shell_script = 'script.sh';
open(OUT,'>',$shell_script) or die $!;
print OUT "#!/bin/bash\n";
my $cleaned_reads_dir = 'cleaned_reads';
my $alignments_dir 	  = 'aligned_cleaned_reads';
my $assembly_dir	  = 'assembled_reads';
my $var_call_dir	  = 'variant_calling';
print OUT "mkdir $cleaned_reads_dir\n" unless (-d($cleaned_reads_dir));
print OUT "mkdir $alignments_dir\n" unless (-d($alignments_dir));
print OUT "mkdir $var_call_dir\n" unless (-d($var_call_dir));

opendir(DIR,$raw_reads_dir) or die $!;
my @subdirs = readdir(DIR);
closedir(DIR) or die $!;
my $sample_num = 0;
my @bam_files_for_variant_calling = ();
foreach my $sample_id (@subdirs)
	{ next if (($sample_id eq '.') or ($sample_id eq '..'));
	  $sample_num++;
	  opendir(DIR,"$raw_reads_dir/$sample_id") or die $!;
	  my @fastq_files = readdir(DIR);
	  closedir(DIR) or die $!;
	  
	  # 1. CLEAN ALL READS
	  my %ids_seen = ();
	  foreach my $fastq_file (@fastq_files)
		{ next if (($fastq_file eq '.') or ($fastq_file eq '..'));
		  next if ($fastq_file !~ /^(.*?)\_\d+\.fq\.gz$/);
		  my $fastq_id = $1;
		  next if (exists($ids_seen{$fastq_id}));
		  my $raw_read_id_file1 = "$raw_reads_dir/$sample_id/$fastq_id.1.fq.gz";
		  my $raw_read_id_file2 = "$raw_reads_dir/$sample_id/$fastq_id.2.fq.gz";
		  my $processed_read_id_file1 = "$cleaned_reads_dir/$fastq_id.forward_paired.fq.gz";
		  my $processed_read_id_file2 = "$cleaned_reads_dir/$fastq_id.forward_unpaired.fq.gz";
		  my $processed_read_id_file3 = "$cleaned_reads_dir/$fastq_id.reverse_paired.fq.gz";
		  my $processed_read_id_file4 = "$cleaned_reads_dir/$fastq_id.reverse_unpaired.fq.gz";
		  print OUT "java -jar $trimmomatic_path PE -phred33 $raw_read_id_file1 $raw_read_id_file2 $processed_read_id_file1 $processed_read_id_file2 $processed_read_id_file3 $processed_read_id_file4 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:100\n";
		  $ids_seen{$fastq_id}++;
		}
	  
	  # 2. ALIGN ALL READS TO THE REFERENCE GENOME, CREATING BAM FILES WHICH WE'LL LATER POST-PROCESS (PARTITIONING OFF UNMAPPED READS; PREPARING THE REMAINDER FOR GATK, ETC.)
	  my $merge_line = '';
	  while((my $file_root,my $irrel)=each(%ids_seen))
		{ my $fwd_paired   = "$cleaned_reads_dir/$file_root.forward_paired.fq.gz";
		  my $rev_paired   = "$cleaned_reads_dir/$file_root.reverse_paired.fq.gz";
		  my $fwd_unpaired = "$cleaned_reads_dir/$file_root.forward_unpaired.fq.gz";
		  my $rev_unpaired = "$cleaned_reads_dir/$file_root.reverse_unpaired.fq.gz";
		  my $unpaired     = "$cleaned_reads_dir/$file_root.unpaired.fq";
		  print OUT "zcat $fwd_unpaired $rev_unpaired > $unpaired\n";
		  print OUT "rm $fwd_unpaired $rev_unpaired\n";
		  my $bam_paired   = "$alignments_dir/$file_root.paired.bam";
		  my $bam_unpaired = "$alignments_dir/$file_root.unpaired.bam";
		  print OUT "bwa mem -R '@RG\tID:group$sample_num\tSM:sample$sample_num\tPL:Illumina\tLIB:lib$sample_num\tPU:unit$sample_num' -M $ref_genome $fwd_paired $rev_paired | samtools view -Sbh - > $bam_paired\n";
		  print OUT "bwa mem -R '@RG\tID:group$sample_num\tSM:sample$sample_num\tPL:Illumina\tLIB:lib$sample_num\tPU:unit$sample_num' -M $ref_genome $unpaired | samtools view -Sbh - > $bam_unpaired\n";
		  $merge_line .= "INPUT=$bam_paired INPUT=$bam_unpaired ";
		  print OUT "rm $fwd_paired $rev_paired\n";
		}
	  my $bam_merged   = "$alignments_dir/$sample_num.bam";
	  print OUT "java -jar $picard_path MergeSamFiles $merge_line OUTPUT=$bam_merged TMP_DIR=$alignments_dir\n";
	  my $bam_mapped   = "$alignments_dir/$sample_num.mapped.bam";
	  my $bam_unmapped = "$alignments_dir/$sample_num.unmapped.bam";
	  print OUT "samtools view -b -F 12 -q 20 $bam_merged > $bam_mapped\n";
	  my $bam_mapped_sorted = "$alignments_dir/$sample_num.mapped.sorted.bam";
	  print OUT "java -jar $picard_path SortSam INPUT=$bam_mapped OUTPUT=$bam_mapped_sorted TMP_DIR=$alignments_dir\n";
	  print OUT "rm $bam_mapped\n";
	  my $bam_mapped_sorted_deduped = "$alignments_dir/$sample_num.mapped.sorted.deduped.bam";
	  my $metrics_file = "$alignments_dir/$sample_num.mapped.sorted.deduped.metrics";
	  print OUT "java -jar $picard_path MarkDuplicates INPUT=$bam_mapped_sorted OUTPUT=$bam_mapped_sorted_deduped METRICS_FILE=$metrics_file ASSUME_SORTED=true\n";
	  print OUT "rm $bam_mapped_sorted\n";
	  print OUT "rm $metrics_file\n";
	  print OUT "java -jar $picard_path BuildBamIndex INPUT=$bam_mapped_sorted_deduped\n";
	  push(@bam_files_for_variant_calling,$bam_mapped_sorted_deduped);
	}
	
# 3. VARIANT CALLING, PART 1: PER-SAMPLE COMMANDS
my $vcfs_for_merging = '';
foreach my $bam (@bam_files_for_variant_calling)
	{ my $sample_num;
	  if ($bam =~ /^.*?\/(\d+)\.mapped\.sorted\.deduped\.bam$/) { $sample_num = $1; }
	  my $intervals_list = "$var_call_dir/$sample_num.intervals.list";
	  print OUT "java -jar $gatk_path -T RealignerTargetCreator -R $ref_genome -I $bam -o $intervals_list\n"; # DETERMINE 'SUSPICIOUS' INTERVALS IN NEED OF LOCAL REALIGNMENT AROUND INDELS
	  my $realigned_bam = "$var_call_dir/$sample_num.realigned.bam";
	  print OUT "java -jar $gatk_path -T IndelRealigner -R $ref_genome -I $bam -targetIntervals $intervals_list -o $realigned_bam\n"; # LOCALLY REALIGN BAM FILES TO CORRECT FOR MISMATCHING BASES THAT COULD OTHERWISE BE MISTAKEN FOR SNPs (CAUSED BY INDELS IN THE GNEOME THAT ARE ABSENT FROM THE REFERENCE)
	  print OUT "rm $intervals_list\n";
	  my $recal_table = "$var_call_dir/$sample_num.recal.data.table";
	  print OUT "java -jar $gatk_path -T BaseRecalibrator -R $ref_genome -I $realigned_bam -knownSites $known_bee_snps -o $recal_table\n"; # RECALIBRATE BASE QUALITY SCORES USING THE SET OF dbSNP SNPs (SO AS TO ENSURE THE RECALIBRATED SCORES ARE MORE ACCURATE IN THAT THEY'RE CLOSER TO THE ACTUAL PROBABILITY OF MISMATCHING THE REFERENCE GENOME)
	  my $realigned_recalibrated_bam = "$var_call_dir/$sample_num.realigned.recalibrated.bam";
	  print OUT "java -jar $gatk_path -T PrintReads -R $ref_genome -I $realigned_bam -BQSR $recal_table -o $realigned_recalibrated_bam\n"; # APPLY THE RECALIBRATED BASE QUALITY SCORES TO THE SEQUENCE DATA
	  print OUT "rm $realigned_bam\n";
	  print OUT "rm $recal_table\n";
	  my $realigned_recalibrated_vcf = "$var_call_dir/$sample_num.realigned.recalibrated.vcf";
	  print OUT "java -jar $gatk_path -T HaplotypeCaller -R $ref_genome -I $realigned_recalibrated_bam --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 -o $realigned_recalibrated_vcf\n"; # CALL VARIANTS FOR EACH SAMPLE INDIVIDUALLY, PRODUCING A VCF (VARIANT CALL FILE)
	  print OUT "rm $realigned_recalibrated_bam\n";
	  $vcfs_for_merging .= "--variant $realigned_recalibrated_vcf ";
	}

# 4. VARIANT CALLING, PART 2: ALL-SAMPLE COMMANDS
my $variants_vcf = "$var_call_dir/variants.vcf";
print OUT "java -jar $gatk_path -T GenotypeGVCFs -R $ref_genome $vcfs_for_merging --max_alternate_alleles 12 -o $variants_vcf\n"; # MERGE THE VCFs OF ALL SAMPLES, SO AS TO CALL VARIANTS ACROSS ALL SAMPLES SIMULTANEOUSLY
my $variants_vcf2 = "$var_call_dir/variants2.vcf";
print OUT "java -jar $gatk_path -T VariantAnnotator -R $ref_genome --dbsnp $known_bee_snps --variant $variants_vcf -o $variants_vcf2 -A QualByDepth -A MappingQualityRankSumTest -A MappingQualityZero -A ReadPosRankSumTest -A FisherStrand -A ClippingRankSumTest -A DepthPerAlleleBySample\n"; # INCORPORATE NECESSARY ANNOTATIONS FOR A VARIANT RECALIBRATION MODELS (AND ADD IN KNOWN BEE SNP IDs, i.e. rsIDs, TOO)
print OUT "rm $variants_vcf\n";
my $recal_file = "$var_call_dir/variants.recal";
my $tranches_file = "$var_call_dir/variants.tranches";
print OUT "java -jar $gatk_path -T VariantRecalibrator -R $ref_genome -input $variants_vcf2 -resource:ensembl,known=false,training=true,truth=false,prior=12.0 $known_bee_snps -resource:dbsnp,known=false,training=true,truth=true,prior=15.0 $known_bee_snps -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an DP -mode SNP -badLodCutoff -3 -recalFile $recal_file -tranchesFiles $tranches_file\n"; # RECALIBRATE VARIANT QUALITY SCORES FOR EACH SNP, WITH dbSNP AS BOTH THE 'TRUTH' SET AND 'TRAINING' SET (THESE ARE LISTS OF SITES FOR WHICH A PRIOR PROBABILITY OF BEING CORRECT HAS BEEN APPLIED, BUT WHICH AREN'T USED BY THE GATK ALGORITHM ITSELF - 'TRAINING' SITES ARE USED TO BUILD A RECALIBRATION MODEL USING INPUT VARIANTS THAT OVERLAP WITH THEM; 'TRUTH' SITES ARE USED TO DETERMINE WHERE TO SET A SENSITIVITY CUTOFF)
my $variants_vqsr_vcf = "$var_call_dir/variants.vqsr.vcf";
print OUT "java -jar $gatk_path -T ApplyRecalibration -R $ref_genome -input $variants_vcf2 --ts_filter_level 99.0 -mode SNP -recalFile $recal_file -tranchesFile $tranches_file -o $variants_vqsr_vcf\n"; # APPLY THE RECALIBRATED VARIANT QUALITY SCORES TO THE DATA, CREATING A FILTERED VCF (HERE NAMED 'VQSR', i.e. VARIANT QUALITY SCORE RECALIBRATION)
print OUT "rm $recal_file $tranches_file $variants_vcf2\n";

close(OUT) or die $!;
exit 1;