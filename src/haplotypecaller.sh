#!/bin/bash

## how to check alternatives of javan version?
##: use `alternatives`

# alternatives --config java

export JAVA_HOME="/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.222.b10-0.el7_6.x86_64/jre/bin/java"
export PATH="$JAVA_HOME:$PATH"
export dirbam="/datasets/data_4/nison/GBS/20201028_GBS_Col0-96/results/02_bowtie2/filtered"
export dirqtl="/datasets/data_5/nison2/GBS/GBS_WT_Rowan/results/snp_filtering"
export dirvcf="$dirqtl/vcf/wt_2020"
export ler_col_snp="/datasets/data_5/nison2/GBS/GBS_WT_Rowan/results/snp_filtering/BC.complete.bed"
export ref="/home/nison/work/refgenome/TAIR10/TAIR10.fasta"
export dirtmp="$dirvcf/tmp"
export dirlog="$dirvcf/log"
export picard="/home/nison/opt/picard/picard.jar"
export rgpu="HHKJ7CCX2"

mkdir -p $dirvcf
mkdir -p $dirtmp
mkdir -p $dirlog

function addrg {
    input=$1
    file=$(basename $input)
    sample=${file%_MappedOn_tair10_sort.md.bam}
    java -jar $picard AddOrReplaceReadGroups \
        I=$input \
        O=${input%.bam}.rg.bam \
        RGID=$sample \
        RGLB=$sample \
        RGPL=ILLUMINA \
        RGPU=$rgpu \
        RGSM=$sample
    }
export -f addrg

parallel --jobs 32 addrg {} ::: $dirbam/*.md.bam
wait;

parallel --jobs 32 samtools index {} ::: $dirbam/*.md.rg.bam
wait;

function hc {
		input=$1
		filename=$(basename $input)
		index=${filename%_MappedOn_tair10_sort.md.rg.bam}

		if [[ -f $dirvcf/${filename%.bam}.g.vcf.gz.tbi ]]; then
				echo ">>>>>>>> skip $index <<<<<<<<"
				exit;
		else
				echo "======== run $index ========"
		fi
		gatk --java-options "-Xmx1g -Djava.io.tmpdir=$dirtmp"\
				HaplotypeCaller \
				-R $ref \
				--intervals $ler_col_snp \
				-ERC GVCF \
				--interval-padding 100 \
				-I $input \
				-O $dirvcf/${filename%.bam}.g.vcf.gz &> $dirlog/${index}.g.vcf.log
		}
export -f hc

parallel --jobs 32 hc {} ::: $dirbam/*.md.rg.bam



# -Xmx4g option in java: set maximum Java heap size as 4g
# --interval-padding100: added 100 bp padding to interval so that the caller sees enough context to reevaluate the call appropriately
# --spark-master local[32]: run on the local machine using 32 cores
