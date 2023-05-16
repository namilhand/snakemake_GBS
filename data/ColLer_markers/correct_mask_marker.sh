#!/bin/bash

origin="/home/nison/work/pipelines/GBS/fromIan/bc.hqsnv.mask.txt"
outdir="/home/nison/work/pipelines/GBS/ColLer_markers"

awk 'BEGIN{IFS="\t"; OFS="\t"}{print substr($2, 2, 4), $3, substr($4, 2, 1), 0, substr($5, 2, 1), 0}' $origin | tail -n +2 > $outdir/BC.mask.tsv 




#awk 'BEGIN{IFS="\t"; OFS="\t"}{print "Chr"$2, $3, substr($4, 2, 1), $5, substr($6, 2, 1), $7}' ./from_ian/BC.complete.tiger.txt | tail -n +2 > BC.complete.tsv
