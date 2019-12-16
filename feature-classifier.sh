#!/bin/bash

#SBATCH -J qiime2_classifier
#SBATCH -A b1042
#SBATCH --mail-type=ALL
#SBATCH --mail-user=elizabeth.mallott@northwestern.edu
#SBATCH -N 1
#SBATCH	-n 8
#SBATCH -t 24:00:00
#SBATCH --output=/home/ekm9460/output_classifier.out
#SBATCH --error=/home/ekm9460/output_classifier.err
#SBATCH -p genomicsguest

module purge all

module load singularity

singularity exec -B /projects/b1057 -B /projects/b1057/liz /projects/b1057/qiime2-core-2019-4.simg qiime feature-classifier classify-sklearn --i-classifier /projects/b1057/gg-13-8-99-ng-classifier.qza --i-reads /projects/b1057/liz/rep-seqs.qza --o-classification /projects/b1057/liz/taxonomy.qza



