#!/bin/bash

cd /home/creggian/dataset/bam
ls ENC*bam | cut -d\. -f1 | xargs -I {} sh -c 'bedtools genomecov -bga -split -ibam {}.bam -g ${HOME}/software/genomes/hg19.chrom.sizes > {}.bdga'
