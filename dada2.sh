#!/bin/bash

#SBATCH -J qiime2_dada2
#SBATCH -A b1042
#SBATCH --mail-type=ALL
#SBATCH --mail-user=elizabeth.mallott@northwestern.edu
#SBATCH -N 1
#SBATCH	-n 8
#SBATCH -t 48:00:00
#SBATCH --output=/home/ekm9460/output_dada2.out
#SBATCH --error=/home/ekm9460/output_dada2.err
#SBATCH -p genomicsguest

module purge all

module load singularity

singularity exec -B /projects/b1057 -B /projects/b1057/liz -B /projects/b1057/liz/leaf_monkeys /projects/b1057/qiime2-core-2019-4.simg qiime dada2 denoise-paired --i-demultiplexed-seqs /projects/b1057/liz/leaf_monkeys/paired-end-demux.qza --p-trunc-len-f 290 --p-trunc-len-r 290 --p-trim-left-f 20 --p-trim-left-r 20 --p-max-ee 5 --p-n-threads 8 --o-table /projects/b1057/liz/leaf_monkeys/table.qza --o-representative-sequences /projects/b1057/liz/leaf_monkeys/rep-seqs.qza --o-denoising-stats /projects/b1057/liz/leaf_monkeys/dada2-stats.qza



