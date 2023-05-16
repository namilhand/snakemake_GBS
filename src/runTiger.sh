#!/bin/bash

# Title : runTiger.sh
# Maintainer : Namil Son
# First edition : 12. 24. 2020
# Last modified : 12. 24. 2020
# Description :
#  This file takes *mask.txt and *complete.txt, which are tiger input files, as input and runs tiger scripts.
# Usage :
# runTiger.sh <*.complete.txt> <*.mask.txt> <dirout>

#---- Setting working environnment --------
# Use JAVA7
#export JAVA_HOME="/usr/local/java/java-se-7u75-ri"
#export PATH="$JAVA_HOME/bin:$PATH"

dirtiger="/usr/local/src/GBS_suite/TIGER_Scripts-for-distribution"

completecall=$1
maskcall=$2
dirout=$3
#maskcall="/home/nison/work/GBS/sandbox/toy/tigerin/toy.mask.txt"
#completecall="/home/nison/work/GBS/sandbox/toy/tigerin/toy.complete.txt"
#dirout="/home/nison/work/GBS/sandbox/toy/tigerout"
prefix=$(basename $completecall)
prefix=${prefix%.complete.txt}

#cd $dirout
##---- Run the base caller on the corrected input file
java -jar $dirtiger/base_caller.jar -r $maskcall -o ${dirout}/${prefix}.mask.basecalls -n bi

#---- Run the allele frequency estimator
java -jar $dirtiger/allele_freq_estimator.jar -r $maskcall -o ${dirout}/${prefix}.mask.allele.frequencies.for.bmm -n bi -w 1000

#---- Apply the beta mixture model
R --slave --vanilla --args ${dirout}/${prefix}.mask.allele.frequencies.for.bmm ${dirout}/${prefix}.mask.bmm.intersections < $dirtiger/beta_mixture_model.R

#---- Prepare files for HMM probabilitiy estimation using the BASECALLER output and the output of the beta mixture model
# Usage :
# perl prep_probl.pl -s <LABEL> -m <corrected input file> -b <base_call_output> -c <chrsizes> -o <output>
# -s specifies sample label
perl $dirtiger/prep_prob.pl -s $prefix -m $maskcall -b ${dirout}/${prefix}.mask.basecalls -c $dirtiger/TAIR10_chrSize.txt -o ${dirout}/${prefix}.mask.probabilities

#---- Calculate transmission and emission probabilities for the HMM
# Usage :
# perl hmm_prob.pl -s <allele_frequencies_for_bmm> -p <probabilities> -o <prefix for the output> -a <output from beta mixture model> -c <chromosome sizes file>
perl $dirtiger/hmm_prob.pl -s ${dirout}/${prefix}.mask.allele.frequencies.for.bmm -p ${dirout}/${prefix}.mask.probabilities -o ${dirout}/${prefix} -a ${dirout}/${prefix}.mask.bmm.intersections -c $dirtiger/TAIR10_chrSize.txt

#---- Run the HMM
# Usage :
# java -jar hmm_play.jar -r <base_call_output> -o <hmm_output> -t bi -z <sample_hmm_model>
java -jar $dirtiger/hmm_play.jar -r ${dirout}/${prefix}.mask.basecalls -o ${dirout}/${prefix}.mask.hmm.out -t bi -z ${dirout}/${prefix}_hmm_model

#---- Get rough estimate of recombination breakpoint positions
# Usage :
# perl prepare_break.pl -s <sample label> -m <corrected_input_file> -b <hmm.out> -c <chromosome size file> -o <output>
perl $dirtiger/prepare_break.pl -s ${dirout}/${prefix} -m $maskcall -b ${dirout}/${prefix}.mask.hmm.out -c $dirtiger/TAIR10_chrSize.txt -o ${dirout}/${prefix}.rough.co
# two output files : 
# 1. ROUGH_CO.txt
# 2. ROUGH_CO.breaks.txt

#---- Refine recombination breaks
# Usage :
# perl refine_recombination_break.pl <complete_input_file> <ROUGH_CO.breaks.txt>
# outputs :
# 1. *.rough.co.recomb.txt
# 2. *.rough.co.refined.breaks.txt
# 3. *.rough.co.refined.recomb.txt
perl $dirtiger/refine_recombination_break_revised.pl $completecall ${dirout}/${prefix}.rough.co.breaks.txt

#---- Smooth out breaks
# Usage : 
# perl breaks_smoother.pl -b <ROUGH_CO.refined.breaks.txt> -o <SMOOTH_CO>
perl $dirtiger/breaks_smoother.pl -b ${dirout}/${prefix}.rough.co.refined.breaks.txt -o ${dirout}/${prefix}.smooth.co.txt

#---- Visualize output
# Usage :

#R --slave --vanilla --args "toy" $dirout/toy.visual.pdf $dirout/toy.rough.co.breaks.txt $dirout/toy.rough.co.refined.breaks.txt $dirout/toy.smooth.co.txt $dirout/toy.mask.allele.frequencies.for.bmm $dirout/toy_sliding_window.breaks.txt < $dirtiger/plot_genotyping.R
