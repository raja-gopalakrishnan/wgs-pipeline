#!/usr/bin/env python

configfile: "config.yaml"
SAMPLES=config["samples"]

#level_one is a list of all sample names
#level_two is a list of the dictionaries associated with each sample name
level_one=[k for k,v in SAMPLES.items()]
level_two=[v for k,v in SAMPLES.items()]

#List of all parent strains
PARENTS=[]
for i in range(0,len(level_two)):
    if level_two[i]["parent"]=="parent" :
        PARENTS.append(level_one[i])

#List of all mutant strains
MUTANTS=[]
for i in range(0,len(level_two)):
    if level_two[i]["parent"]!="parent" :
        MUTANTS.append(level_one[i])

rule all:
	input:
		#make_barcode_file
		"barcodes.tsv",
		#fastqc_raw
		"qual_ctrl/raw",
		#demultiplex
		expand("fastq/{sample}.fastq.gz", sample=SAMPLES),
		#samtools_sort
		expand("alignment/{sample}.bam", sample=SAMPLES),
		#bowtie_summary
		"alignment/bowtie_read_number_summary.txt",
		#gzip_loose_fastq
		expand("fastq/trimmed/{sample}.trimmed.fastq.gz",sample=SAMPLES),
		#index_bam
		expand("alignment/{sample}_cleaned.bam.bai", sample =SAMPLES),
		#call_variants
		expand("raw_variants/{sample}_raw_variants.vcf", sample = SAMPLES),
		#combine_filtered_variants
		expand("raw_variants/{sample}_raw_variants_all_filtered.vcf", sample = SAMPLES),
		#merge_vcf
		expand("mergedVCFs/{mutant}_parent_merged.vcf", mutant=MUTANTS),
		#get_mutant_snps
		expand("unique_variants/{mutant}_only_passfilter.vcf", mutant=MUTANTS),
		#vcf_to_table
		#expand("unique_variants/{mutant}_only_passfilter.table", mutant=MUTANTS),
		#get_genes
		#expand("unique_variants/{mutant}_only_genes.txt", mutant=MUTANTS),
		#join_table_genes
		expand("unique_variants/{mutant}_only_passfilter_annotated.table", mutant=MUTANTS),
		#select_orfs
		expand("unique_variants/{mutant}_only_passfilter_ORFs.table", mutant=MUTANTS)
#Make a tab separated file with column 1 containing name of lib and column 2 containing the barcode sequence
rule make_barcode_file:
	output:
		"barcodes.tsv"
	params:
		bc = config["samples"]
	run:
		with open(output[0],"w") as f:
			for x in params.bc:
				f.write('\t'.join((x,params.bc[x]["barcode"])) + '\n')
#fastqc analysis
rule fastqc_raw:
	input:
		config["fastq"]
	output:
		"qual_ctrl/raw"
	threads: config["threads"]
	log: "logs/fastqc/fastqc-raw.log"
	shell: """
	mkdir {output}
	(fastqc -o {output} --noextract -t {threads} {input}) &> {log}
	"""
##Demultiplex fastq file
rule demultiplex:
	input:
		fastq = config["fastq"],
		barcode = "barcodes.tsv"
	output:
		expand("fastq/{sample}.fastq.gz", sample=SAMPLES)
	log: "logs/demultiplex.log"
	params:
		mismatch = config["demultiplex"]["mismatch"]
	shell:"""
	(fastq-multx -B {input.barcode} -b {input.fastq} -m {params.mismatch} -o fastq/%.fastq.gz) &> {log}
	"""
#Cutadapt - trim 1 base from 5' end of every file
#Also trims low quality reads from the 3' end and removes reads that are less than 'n' nucleotides in length after trimming
rule cutadapt:
	input:
		"fastq/{sample}.fastq.gz"
	output:
		temp("fastq/trimmed/{sample}.trimmed.fastq")
	params:
		qual_cutoff = config["cutadapt"]["qual_cutoff"],
		min_length = config["cutadapt"]["min_length"]
	log: "logs/cutadapt/cutadapt-{sample}.log"
	shell:"""
	(cutadapt -u 1 --nextseq-trim={params.qual_cutoff} -m {params.min_length} -o {output} {input}) &> {log}
	"""
#Alignment to the genome using bowtie
rule bowtie:
	input:
		fastq = "fastq/trimmed/{sample}.trimmed.fastq"
	output:
		sam = temp("alignment/{sample}.sam"),
		unaligned = temp("alignment/fastq/unaligned-{sample}.fastq"),
		aligned = temp("alignment/fastq/aligned-{sample}.fastq")
	params:
		index = "genome/" + config["species"] + "/" + config["basename"]
	threads: config["threads"]
	log: "logs/bowtie/bowtie-align-{sample}.log"
	shell:"""
	(bowtie2 -p {threads} -x {params.index} -U {input.fastq} --al {output.aligned} --un {output.unaligned} -S {output.sam}) &> {log} 
	"""
#Convert SAM file to BAM and sort the file
rule samtools_sort:
	input:
		"alignment/{sample}.sam"
	output:
		"alignment/{sample}.bam"
	log: "logs/samtools_sort/samtools-sort-{sample}.log"
	threads: config["threads"]
	shell:"""
	(samtools view -buh -q 3 {input} | samtools sort -T {wildcards.sample} -@ {threads} -o {output} -) &> {log}
	"""
#Gather summary statistics for alignment
rule bowtie_summary:
	input:
		expand("alignment/{sample}.bam", sample=SAMPLES)
	output:
		number = "alignment/bowtie_read_number_summary.txt",
		percentage = "alignment/bowtie_read_percentage_summary.txt"
	script: "scripts/extract_bowtie_log.R"
#Zip the trimmed fastq and unaligned/aligned fastq files
rule gzip_loose_fastq:
	input:
		"alignment/{sample}.sam",
		trim = "fastq/trimmed/{sample}.trimmed.fastq",
		aligned = "alignment/fastq/unaligned-{sample}.fastq",
		unaligned = "alignment/fastq/aligned-{sample}.fastq"
	output:
		"fastq/trimmed/{sample}.trimmed.fastq.gz",
		"alignment/fastq/unaligned-{sample}.fastq.gz",
		"alignment/fastq/aligned-{sample}.fastq.gz"
	shell:"""
	pigz -fk {input.trim}
	pigz -fk {input.aligned}
	pigz -fk {input.unaligned}
	"""
#Add read groups to BAM files
rule add_read_group:
	input:
		"alignment/{sample}.bam"
	output:
		"alignment/{sample}_cleaned.bam"
	params:
		rgid = lambda wildcards: config["fastq"] + "_" + wildcards.sample,
		rglb = lambda wildcards: "library_" + wildcards.sample,
		rgpl = "illumina",
		rgsm = lambda wildcards: wildcards.sample,
		rgpu = lambda wildcards: config["samples"][wildcards.sample]["barcode"]
	log: "logs/add_read_group/add_read_group-{sample}.log"
	shell:"""
	(picard AddOrReplaceReadGroups I={input} O={output} RGLB={params.rglb} RGPL={params.rgpl} RGSM={params.rgsm} RGPU={params.rgpu}) &> {log}
	"""
#Index the bam files (required for GATK functions)
rule index_bam:
	input:
		"alignment/{sample}_cleaned.bam"
	output:
		"alignment/{sample}_cleaned.bam.bai"
	log: "logs/index_bam/index_bam-{sample}.log"
	shell:"""
	(samtools index {input}) &> {log}
	"""
# Call low quality variants using HaplotypeCaller in GATK
rule call_variants:
	input:
		"alignment/{sample}_cleaned.bam.bai",
		bam = "alignment/{sample}_cleaned.bam"
	output:
		vcf = "raw_variants/{sample}_raw_variants.vcf"
	params:
		genome = "genome/" + config["genome"] + ".fa",
		conf = config["call_variants"]["stand-call-conf"]
	log: "logs/call_variants/call_variants-{sample}.log"
	shell:"""
	(gatk-launch HaplotypeCaller -R {params.genome} -I {input.bam} --genotyping-mode DISCOVERY -stand-call-conf {params.conf} -O {output.vcf}) &> {log}
	"""
# Select snps from variants
rule select_snps:
	input:
		"raw_variants/{sample}_raw_variants.vcf"
	output:
		"raw_variants/{sample}_raw_variants_snps.vcf"
	params:
		genome = "genome/" + config["genome"] + ".fa"
	log: "logs/select_snps/select_snps-{sample}.log"
	shell:"""
	(gatk-launch SelectVariants -R {params.genome} -V {input} -select-type SNP -O {output}) &> {log}
	"""
# Select indels from variants
rule select_indels:
	input:
		"raw_variants/{sample}_raw_variants.vcf"
	output:
		"raw_variants/{sample}_raw_variants_indels.vcf"
	params:
		genome = "genome/" + config["genome"] + ".fa"
	log: "logs/select_snps/select_snps-{sample}.log"
	shell:"""
	(gatk-launch SelectVariants -R {params.genome} -V {input} -select-type INDEL -O {output}) &> {log}
	"""
#Hard-filtering of SNPs
rule filter_snps:
	input:
		"raw_variants/{sample}_raw_variants_snps.vcf"
	output:
		"raw_variants/{sample}_raw_variants_snps_filtered.vcf"
	params:
		genome = "genome/" + config["genome"] + ".fa",
		qd = config["filter_snps"]["qd"],
		fs = config["filter_snps"]["fs"],
		mq = config["filter_snps"]["mq"],
		sor = config["filter_snps"]["sor"],
		mqrs = config["filter_snps"]["mqrs"],
		rprs = config["filter_snps"]["rprs"]
	log: "logs/filter_snps/filter_snps-{sample}.log"
	shell:"""
	(gatk-launch VariantFiltration -R {params.genome} -V {input} --filter-expression 'QD < {params.qd} || FS > {params.fs} || MQ < {params.mq} || SOR > {params.sor} || MQRankSum < {params.mqrs} || ReadPosRankSum < {params.rprs}' --filter-name 'Fails GATK filter' -O {output}) &> {log}
	"""
#Hard-filtering of INDELs
rule filter_indels:
	input:
		"raw_variants/{sample}_raw_variants_indels.vcf"
	output:
		"raw_variants/{sample}_raw_variants_indels_filtered.vcf"
	params:
		genome = "genome/" + config["genome"] + ".fa",
		qd = config["filter_indels"]["qd"],
		fs = config["filter_indels"]["fs"],
		rprs = config["filter_indels"]["rprs"]
	log: "logs/filter_snps/filter_snps-{sample}.log"
	shell:"""
	(gatk-launch VariantFiltration -R {params.genome} -V {input} --filter-expression 'QD < {params.qd} || FS > {params.fs} || ReadPosRankSum < {params.rprs}' --filter-name 'Fails GATK filter' -O {output}) &> {log}
	"""
#Combine filtered indels and snps into single file
#This rule uses GATK3.8
rule combine_filtered_variants:
	input:
		snps = "raw_variants/{sample}_raw_variants_snps_filtered.vcf",
		indels = "raw_variants/{sample}_raw_variants_indels_filtered.vcf"
	output:
		"raw_variants/{sample}_raw_variants_all_filtered.vcf"
	params:
		genome = "genome/" + config["genome"] + ".fa"
	log: "logs/combine_filtered_variants/combine_filtered_variants-{sample}.log"
	shell:"""
	(gatk -T CombineVariants -R {params.genome} -V:SNPs {input.snps} -V:INDELs {input.indels} -o {output} -genotypeMergeOptions PRIORITIZE -priority SNPs,INDELs) &> {log}
	"""
#Combine parent and mutant VCFs
#This rule uses GATK3.8
rule merge_vcf:
	input:
		mutant = "raw_variants/{mutant}_raw_variants_all_filtered.vcf",
		parent = lambda wildcards: "raw_variants/" + config["samples"][wildcards.mutant]["parent"] + "_raw_variants_all_filtered.vcf"
	output:
		"mergedVCFs/{mutant}_parent_merged.vcf"
	params:
		genome = "genome/" + config["genome"] + ".fa"
	log: "logs/merge_vcf/merge_vcf-{mutant}.log"
	shell:"""
	(gatk -T CombineVariants -R {params.genome} -V:mutant {input.mutant} -V:parent {input.parent} -o {output}) &> {log}
	"""
#Pull out variants that pass filter and are only there in mutant
#This rule uses GATK3.8
rule get_mutant_snps:
	input:
		"mergedVCFs/{mutant}_parent_merged.vcf"
	output:
		"unique_variants/{mutant}_only_passfilter.vcf"
	params:
		genome = "genome/" + config["genome"] + ".fa"
	log: "logs/get_mutant_snps/get_mutant_snps-{mutant}.log"
	shell:"""
	(gatk -T SelectVariants -R {params.genome} -V {input} -select 'set=="mutant"' -o {output}) &> {log}
	"""
#Convert VCF file to a table by picking out necessary fields only
rule vcf_to_table:
	input:
		"unique_variants/{mutant}_only_passfilter.vcf"
	output:
		temp("unique_variants/{mutant}_only_passfilter.table")
	params:
		genome = "genome/" + config["genome"] + ".fa"
	log: "logs/vcf_to_table/vcf_to_table-{mutant}.log"
	shell:"""
	(gatk-launch VariantsToTable -R {params.genome} -V {input} -F CHROM -F POS -F REF -F ALT -F QUAL -GF AD -GF DP -O {output}) &> {log}
	"""
#Convert the vcf file to a bed file and map the snps to genes
rule get_genes:
	input:
		"unique_variants/{mutant}_only_passfilter.vcf"
	output:
		temp("unique_variants/{mutant}_only_passfilter_genes.txt")
	params:
		anno = "genome/annotations/" + config["species"] + "_ORFs.bed"
	log: "logs/set_genes/set_genes-{mutant}.log"
	shell:"""
	(intersectBed -a <(vcf2bed <{input}) -b {params.anno} -wao | sort -k1,1 -k2,2n - |awk '{{FS=OFS="\t"}}{{print $16}}' > {output}) &> {log}
	"""
#Append gene name to table file
rule join_table_genes:
	input:
		genes = "unique_variants/{mutant}_only_passfilter_genes.txt",
		table = "unique_variants/{mutant}_only_passfilter.table"
	output:
			"unique_variants/{mutant}_only_passfilter_annotated.table"
	log:"logs/join_table_genes/join_table_genes-{mutant}.log"
	shell:"""
	(paste -d '\t' <(cat <(head -1 {input.table}) <(sort -k1,1 -k2,2n <(sed '1d' {input.table}))) <(sed '1i GENE_ID' {input.genes}) > {output}) &> {log}
	"""
#Select ORFs only
rule select_orfs:
	input:
		"unique_variants/{mutant}_only_passfilter_annotated.table"
	output:
		"unique_variants/{mutant}_only_passfilter_ORFs.table"
	log:"logs/select_orfs/select_orfs-{mutant}.log"
	shell:"""
	(cat <(head -1 {input}) <(awk '{{OFS=FS="\t"}}{{if($10!="."){{print $0}}}}' {input} | sed '1d' | sort -k5,5nr) > {output}) &> {log}
	"""