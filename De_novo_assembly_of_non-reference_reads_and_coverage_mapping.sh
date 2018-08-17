#!/bin/bash

#########################################################################################################

#Begin in this directory: <home>
#Requirement: Must have this python script in this directory: "Contig_filter_len.cov.py"
#Leftover reads following host ref assembly = Reads form each sample are concatenated keeping 1st, 2nd pairs and unpaired together (3 .fq files from 23 samples). Here called combined.1.fq combined.2.fq combined.unpaired.fq 
#Original leftover reads located here: "<sample_name>.<1_or_2_or_unpaired>.fq"
#Cutoff for contigs are set at length 1000 bp and coverage of 2.

#########################################################################################################
##Add the following modules:
# module add apps/gcc/bwa/0.7.8
# module add apps/gcc/samtools/1.2
# module add compilers/java/1.8.0_60 # for Picard
# module add compilers/python/2.7.9 # for Blobtools
# module add apps/gcc/SPAdes/3.8.1
# module load apps/gcc/samtools/1.4

##5 rounds of assembly

#iteration1
cd <home>
mkdir iteration1
spades.py --meta -t 96 -1 combined.1.fq -2 combined.2.fq -s combined.unpaired.fq -o <home>/iteration1
python Contig_filter_len.cov.py 1000 2 iteration1/scaffolds.fasta

bwa index -a bwtsw iteration1/scaffolds.fltr.fasta
bwa mem -M -t 96 iteration1/scaffolds.fltr.fasta combined.unpaired.fq | samtools view -bh - | samtools sort - iteration1/it1.unpaired
bwa mem -M -t 96 iteration1/scaffolds.fltr.fasta combined.1.fq combined.2.fq| samtools view -bh - | samtools sort - iteration1/it1.paired
java -jar picard-tools-1.138/picard.jar MergeSamFiles SORT_ORDER='coordinate' OUTPUT= iteration1/it1.bam INPUT=iteration1/it1.paired.bam INPUT=iteration1/it1.unpaired.bam
rm iteration1/it1.unpaired.bam iteration1/it1.paired.bam 
samtools view -h -b -f 4 iteration1/it1.bam > iteration1/it1.unmapped.bam
java -jar picard-tools-1.138/picard.jar SamToFastq INPUT=iteration1/it1.unmapped.bam FASTQ=iteration1/it1.unmapped.1.fq SECOND_END_FASTQ=iteration1/it1.unmapped.2.fq UNPAIRED_FASTQ=iteration1/it1.unmapped.unpaired.fq
rm iteration1/it1.bam iteration1/it1.unmapped.bam
gzip *.fq 

#iteration2
mkdir iteration2
spades.py --meta -t 96 -1 iteration1/it1.unmapped.1.fq -2 iteration1/it1.unmapped.2.fq -s iteration1/it1.unmapped.unpaired.fq -o <home>/iteration2
python Contig_filter_len.cov.py 1000 2 iteration2/scaffolds.fasta

bwa index -a bwtsw iteration2/scaffolds.fltr.fasta
bwa mem -M -t 96 iteration2/scaffolds.fltr.fasta iteration1/it1.unmapped.unpaired.fq | samtools view -bh - | samtools sort - iteration2/it2.unpaired
bwa mem -M -t 96 iteration2/scaffolds.fltr.fasta iteration1/it1.unmapped.1.fq iteration1/it1.unmapped.2.fq| samtools view -bh - | samtools sort - iteration2/it2.paired
java -jar picard-tools-1.138/picard.jar MergeSamFiles SORT_ORDER='coordinate' OUTPUT= iteration2/it2.bam INPUT=iteration2/it2.paired.bam INPUT=iteration2/it2.unpaired.bam
rm iteration2/it2.unpaired.bam iteration2/it2.paired.bam 
samtools view -h -b -f 4 iteration2/it2.bam > iteration2/it2.unmapped.bam
java -jar picard-tools-1.138/picard.jar SamToFastq INPUT=iteration2/it2.unmapped.bam FASTQ=iteration2/it2.unmapped.1.fq SECOND_END_FASTQ=iteration2/it2.unmapped.2.fq UNPAIRED_FASTQ=iteration2/it2.unmapped.unpaired.fq
rm iteration2/it2.bam iteration2/it2.unmapped.bam
gzip iteration1/*.fq

#iteration3
mkdir iteration3
spades.py --meta -t 96 -1 iteration2/it2.unmapped.1.fq -2 iteration2/it2.unmapped.2.fq -s iteration2/it2.unmapped.unpaired.fq -o <home>/iteration3
python Contig_filter_len.cov.py 1000 2 iteration3/scaffolds.fasta

bwa index -a bwtsw iteration3/scaffolds.fltr.fasta
bwa mem -M -t 96 iteration3/scaffolds.fltr.fasta iteration2/it2.unmapped.unpaired.fq | samtools view -bh - | samtools sort - iteration3/it3.unpaired
bwa mem -M -t 96 iteration3/scaffolds.fltr.fasta iteration2/it2.unmapped.1.fq iteration2/it2.unmapped.2.fq| samtools view -bh - | samtools sort - iteration3/it3.paired
java -jar picard-tools-1.138/picard.jar MergeSamFiles SORT_ORDER='coordinate' OUTPUT= iteration3/it3.bam INPUT=iteration3/it3.paired.bam INPUT=iteration3/it3.unpaired.bam
rm iteration3/it3.unpaired.bam iteration3/it3.paired.bam 
samtools view -h -b -f 4 iteration3/it3.bam > iteration3/it3.unmapped.bam
java -jar picard-tools-1.138/picard.jar SamToFastq INPUT=iteration3/it3.unmapped.bam FASTQ=iteration3/it3.unmapped.1.fq SECOND_END_FASTQ=iteration3/it3.unmapped.2.fq UNPAIRED_FASTQ=iteration3/it3.unmapped.unpaired.fq
rm iteration3/it3.bam iteration3/it3.unmapped.bam
gzip iteration2/*.fq 

#iteration4
mkdir iteration4
cd iteration3
spades.py --careful -t 40 -s it3.unmapped.unpaired.fq -o <home>/iteration4
cd ..
python Contig_filter_len.cov.py 1000 2 iteration4/scaffolds.fasta

bwa index -a bwtsw iteration4/scaffolds.fltr.fasta
bwa mem -M -t 40 iteration4/scaffolds.fltr.fasta iteration3/it3.unmapped.unpaired.fq | samtools view -bh - | samtools sort - iteration4/it4.unpaired
java -jar picard-tools-1.138/picard.jar MergeSamFiles SORT_ORDER='coordinate' OUTPUT=iteration4/it4.bam INPUT=iteration4/it4.unpaired.bam
rm iteration4/unpaired.bam iteration4/paired.bam 
samtools view -h -b -f 4 iteration4/it4.bam > iteration4/it4.unmapped.bam
java -jar picard-tools-1.138/picard.jar SamToFastq INPUT=iteration4/it4.unmapped.bam FASTQ=iteration4/it4.unmapped.1.fq SECOND_END_FASTQ=iteration4/it4.unmapped.2.fq UNPAIRED_FASTQ=iteration4/it4.unmapped.unpaired.fq
rm iteration4/it4.bam iteration4/it4.unmapped.bam
gzip iteration3/*.fq

#iteration5
mkdir iteration5
cd iteration4
spades.py --careful -t 40 -s it4.unmapped.unpaired.fq -o <home>/iteration5
cd ..
python Contig_filter_len.cov.py 1000 2 iteration5/scaffolds.fasta

bwa index -a bwtsw iteration5/scaffolds.fltr.fasta
bwa mem -M -t 40 iteration5/scaffolds.fltr.fasta iteration4/it4.unmapped.unpaired.fq | samtools view -bh - | samtools sort - iteration5/it5.unpaired
java -jar picard-tools-1.138/picard.jar MergeSamFiles SORT_ORDER='coordinate' OUTPUT=iteration5/it5.bam INPUT=iteration5/it5.unpaired.bam
rm iteration5/unpaired.bam iteration5/paired.bam 
samtools view -h -b -f 4 iteration5/it5.bam > iteration5/it5.unmapped.bam
java -jar picard-tools-1.138/picard.jar SamToFastq INPUT=iteration5/it5.unmapped.bam FASTQ=iteration5/it5.unmapped.1.fq SECOND_END_FASTQ=iteration5/it5.unmapped.2.fq UNPAIRED_FASTQ=iteration5/it5.unmapped.unpaired.fq
rm iteration5/it5.bam iteration5/it5.unmapped.bam
gzip iteration4/*.fq
gzip iteration5/*.fq

###########################################################################################################################
#After 5 iterations the assembly is complete. Now we combine assemblies, and align reads from each sample to these contigs.
###########################################################################################################################

mkdir Samples_aligned
cat */*fltr.fasta > Samples_aligned/All_contigs.fasta
bwa index -a bwtsw Samples_aligned/All_contigs.fasta
cd Samples_aligned

#Each sample name is located in: <Sample_names_list>
#We now align original leftover reads to the assembled contigs
#idxstats from samtools gives mapped reads/contig in 3rd column (unmapped in 4th, 1st is name, 2nd is length of contig)

for i in `cat <Sample_names_list> | awk {'print $1'}`; do
	mkdir $i
	bwa mem -M -t 96 All_contigs.fasta ../../../Unmapped_post_alignment/$i.unpaired.fq | samtools view -bh - | samtools sort - $i/unpaired
	bwa mem -M -t 96 All_contigs.fasta ../../../Unmapped_post_alignment/$i.1.fq ../../../Unmapped_post_alignment/$i.2.fq | samtools view -bh - | samtools sort - $i/paired
	java -jar picard-tools-1.138/picard.jar MergeSamFiles SORT_ORDER='coordinate' OUTPUT= $i/$i.bam INPUT=$i/paired.bam INPUT=$i/unpaired.bam
	samtools index $i/$i.bam
	samtools idxstats $i/$i.bam | sort -k 1 > $i/$i.reads.txt
done

####################################################
##Open Blobtools##
####################################################
cd blobtools/
source bin/activate

#For each sample, get a bam2cov for base coverage etc.
#Needs Samtools

for i in `cat <Sample_names_list> | awk {'print $1'}`; do
	./blobtools map2cov -i <home>/Samples_aligned/All_contigs.fasta -b <home>/Samples_aligned/$i/$i.bam -o <home>/Samples_aligned/$i/
done
#done for individual samples, now for all samples together
cd <home>/Samples_aligned
gunzip ../*.gz

bwa mem -M -t 96 All_contigs.fasta ../combined.unpaired.fq | samtools view -bh - | samtools sort - ../unpaired
bwa mem -M -t 96 All_contigs.fasta ../combined.1.fq ../combined.2.fq | samtools view -bh - | samtools sort - ../paired
java -jar picard-tools-1.138/picard.jar MergeSamFiles SORT_ORDER='coordinate' OUTPUT= ../combined.bam INPUT=../paired.bam INPUT=../unpaired.bam
rm ../paired.bam ../unpaired.bam

cd blobtools/
source bin/activate

./blobtools map2cov -i <home>/Samples_aligned/All_contigs.fasta -b <home>/Samples_aligned/combined.bam -o <home>/

cd <home>/Samples_aligned

#Now take base cov from each sample/contig and paste together in a table
#Required python file Base_cov.py
for i in `cat <Sample_names_list> | awk {'print $1'}`; do sort $i/*.cov -k 1 > $i/sorted.cov; done
for i in `cat <Sample_names_list> | awk {'print $1'}`; do python Base_cov.py $i $i/sorted.cov $i/$i.base_score.txt; done

paste Contigs.txt */*base_score.txt > Basecov_per_contig.txt
