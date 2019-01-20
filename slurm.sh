#!/bin/bash

#SBATCH -p priority
#SBATCH -t 12:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH -c 1
#SBATCH -e superjob.err
#SBATCH -o superjob.log
#SBATCH -J WGS-snakemake

snakemake -p --latency-wait 300 --rerun-incomplete --cluster-config cluster.yaml --use-conda --jobs 999 --cluster "sbatch -p {cluster.queue} -c {cluster.n} -t {cluster.time} --mem-per-cpu={cluster.mem} -J {cluster.name} -e {cluster.err} -o {cluster.log}"
