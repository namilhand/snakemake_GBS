#/bin/bash

# This script convert [BC.complete.tiger.txt], which is shared by Ian, to [BC.complete.vcf] which will be used as argument of [bcftools call -T] option.
# Note that the file for [bcftools call -T] option must be compressed and indexed.
# Use [tabix], one of the samtools, to index tab-delimited genome position files.

awk 'BEGIN{IFS="\t"; OFS="\t"}{print "Chr"$2, $3, substr($4, 2, 1), $5, substr($6, 2, 1), $7}' ./from_ian/BC.complete.tiger.txt | tail -n +2 > BC.complete.tsv

