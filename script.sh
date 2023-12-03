#!bin/bash

conda activate ngs

# Creation index
idxFile=chr16.fa.gz
bwa index -a bwtsw index/"$idxFile"
#
# --------- Pipeline --------------------
if [ ! -d fastqc ]; then
	mkdir fastqc
fi

# Exécution du fastqc sur tous les fichiers fastq
fastqc -o fastqc patient7.exome/*.fastq.gz

if [ ! -d trimmomatic ]; then
	mkdir trimmomatic
fi

# Trimmomatic sur toutes les paires de fastq
for read1_file in patient7.exome/*_r1F.fastq.gz; do
	# Extract the name without the extension
	name="${read1_file%_r1F.fastq.gz}"

	# Define the output file name
	output=$(basename "${name}")
	# Run Trimmomatic to trim and filter reads
	trimmomatic PE "${name}_r1F.fastq.gz" "${name}_r2F.fastq.gz" \
		trimmomatic/"${output}.R1_paired.fastq" trimmomatic/"${output}.R1_unpaired.fastq" \
		trimmomatic/"${output}.R2_paired.fastq" trimmomatic/"${output}.R2_unpaired.fastq" \
		LEADING:20 TRAILING:20 MINLEN:50
done

if [ ! -d sam ]; then
	mkdir sam
fi

#
for read1_file in trimmomatic/*.R1_paired.fastq; do
	# Extract the name without the extension
	name=$(basename "${read1_file%.R1_paired.fastq}")

	bwa mem -M -t 2 -A 2 -E 1 index/"$idxFile" trimmomatic/"$name".R1_paired.fastq trimmomatic/"$name".R2_paired.fastq > \
		sam/"$name"_map.sam
done

gunzip index/"$idxFile"

for sam in sam/*.sam; do
	name="${sam%.sam}"
	#Sam 2 Bam
	samtools view -S -b "$sam" >"$name".bam

	# flagstats
	samtools flagstat "$sam"

	#Sort Bam
	samtools sort "$name".bam >"$name"_sorted.bam

	#Index bam file
	samtools index "$name"_sorted.bam
	#Convert to Mpileup
	samtools mpileup -B -A -f index/chr16.fa "$name"_sorted.bam > \
		"$name".msf
done

if [ ! -d varscan ]; then
	mkdir varscan
fi

if [ ! -d results ]; then
	mkdir results
fi

varscan somatic sam/TCRBOA7-N-WEX-chr16_map.msf sam/TCRBOA7-T-WEX-chr16_map.msf varscan/varscan.vcf

grep -i 'somatic' varscan/varscan.vcf.indel >results/filtered.indel
awk '{OFS="\t"; if (!/^#/){print $1,$2-1,$2,$4"/"$5,"+"}}' \
	results/filtered.indel >results/indel.bed
grep -i 'somatic' varscan/varscan.vcf.snp >results/filtered.snp
awk '{OFS="\t"; if (!/^#/){print $1,$2-1,$2,$4"/"$5,"+"}}' \
	results/filtered.snp >results/snp.bed

gunzip annotation/gencode.v24lift37.basic.annotation.gtf.gz

bedtools intersect -a annotation/gencode.v24lift37.basic.annotation.gtf -b results/indel.bed >results/intersecIndel
grep '\sgene\s' results/intersecIndel | awk '{print " " $1 " " $4 " " $5 " " $16}' >results/resultatIndel

bedtools intersect -a annotation/gencode.v24lift37.basic.annotation.gtf -b results/snp.bed >results/intersecSnp
grep '\sgene\s' results/intersecSnp | awk '{print " " $1 " " $4 " " $5 " " $16}' >results/resultatSnp
