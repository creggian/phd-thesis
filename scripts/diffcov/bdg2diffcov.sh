#!/bin/bash

cd /home/creggian/dataset/bam

files=$(ls *bdga)
for f in $files;
do
     bedtools closest -io -iu -a ${f} -b ${f} -D a | awk '{printf("%s\t%s\t%s\t%s\n", $1, $3, $4-$8, $9)}' > ${f}.dc
done
