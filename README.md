# Whole genome sequencing pipeline
This pipeline can be used for identifying SNPs and indels in yeast mutants (*S. cerevisiae* and *S. pombe*) obtained from a genetic screen or selection.

## Getting started

### Required software
- Unix operating system
- git
- conda

### Input files
- FASTQ file (samples multiplexed with in-line barcodes). Sequening runs must contain at least one parent strain - the original strain that was used for the genetic screen or selection, which will be used to filter out non-causal variants.


## Instructions

### 1. Clone the repository
```
git clone https://github.com/winston-lab/wgs_snp_analysis_2018.git
```

### 2. Edit the config file
Copy the template_config.yaml file to make the config.yaml file
```
cp template_config.yaml config.yaml
```
Open the config.yaml file and edit the following items:
- ```samples``` : Edit sample names, barcode and parent sample name. If doing bulk segregant analysis, the library made from pooled spores that do not have the mutant phenotype should be treated as the parent.
- species specific parameters: Specify the ```species```, ```basename``` and ```genome``` - *S. cerevisiae* or *S. pombe*
- ```fastq```: Provide location of fastq file
Change other parameters as required. Feel free to play around with the parameters for hard filtering of SNPs.

### 3. Create the environment to run the pipeline
```
# Create an environment where all packages required for running the pipeline are installed
conda env create -f envs/snakemake_default.yaml

# Activate the environment
source activate snakemake_default
```

### 4. Install GATK 3.8
This pipeline primarily uses GATK 4, which is installed by conda. But for a couple of steps, it also uses GATK 3.8. Installing GATK 3.8 requires downloading the jar file for GATK from [here](https://software.broadinstitute.org/gatk/download/auth?package=GATK-archive&version=3.8-0-ge9d806836). Copy the file into the home directory and run the following command:
```
gatk-update GenomeAnalysisTK-3.8-0-ge9d806836.tar.bz2
```

### 5. Run the pipeline
```
# Do a dry run of the pipeline which will show you the number of jobs that will be submitted
snakemake -np
```

If running the pipeline locally on your computer use the following command to run it:
```
snakemake -p
```

If running the pipeline on an HPC cluster (recommended) that uses the SLURM job scheduler, use the following command:
```
sbatch slurm.sh
```
If using a HPC cluster with a different job scheduler, modify the ```slurm.sh``` and ```cluster.yaml``` files accordingly

## Interpreting the results
The ```unique_variants/``` directory will have the results of interest. A variant is defined as a SNP or indel that is present in the sample as compared to the reference genome. Every mutant sample will have six files (plus one indexed VCF file) in this folder:
#### 1. ```sample_only_nofilter.vcf```
VCF file that contains a list of variants that are present only in the mutant and not in the parent.
#### 2. ```sample_only_nofilter.table```
Contains same information as the ```sample_only_nofilter.vcf``` file, except in a user-friendly table format while displaying only certain fields of interest.
#### 3. ```sample_only_nofilter_genes.txt```
Map the variant in the ```sample_only_nofilter.vcf``` to genes and print a list of the gene names. The presence of '.' indicates that a variant was mapped to an intergenic region.
#### 4. ```sample_only_nofilter_annotated.table```
Append the gene names to the ```sample_only_nofilter.table``` file.
#### 5. ```sample_only_nofilter_ORFs.table```
Contains a list of variants that are found in the coding regions of genes.
#### 6. ```sample_only_passfilter_ORFs.table```
Contains a list of variant that are found in the coding regions of genes and passed the hard filtering of GATK.

The pipeline was designed primarily for identifying variants present within coding regions. So, when the pipeline is done, I usually take a look at ```sample_only_passfilter_ORFs.table```. But also check the ```sample_only_nofilter_ORFs.table``` and ```sample_only_nofilter_annotated.table``` files for variants that did not pass the GATK filtering parameters and for variants that mapped to the intergenic regions respectively.

## Troubleshooting
The ```call_variants``` step is the most memory intensive process. For samples that I have used, providing 16 GB of memory has proved to be sufficient. However, depending on the number of reads in the sample, the memory requirement might be higher. You can edit the amount of memory provided for this step in the ```cluster.yaml``` file.

## Acknowledgements
- Thank you James Chuang, for providing the *S. cerevisiae* and *S. pombe* annotation files!
- Thanks to Dan Pagano, whose script helped me understand how to use GATK!
