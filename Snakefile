# Snakemake workflow for TIGER analysis

# Usage (snakemake --cores should reflect available cores):
# conda env create --file environment.yaml --name GBS 
# conda activate GBS
# snakemake -p --cores 30 
# conda deactivate

import pandas as pd
import os

# To make a shell script invoked using the "shell" directive,
# we need to determine the base path of Snakefile since we Expect the scripts directory to be there as well
SRCDIR = srcdir("")

# Specify config file parameters
configfile: "config.yaml"
## library info
SAMPLE_SIZE = config["LIBRARY_INFO"]["sample_size"]
SAMPLE = list(range(1, SAMPLE_SIZE + 1))
LIB_NAME = config["LIBRARY_INFO"]["lib_name"]
## reference genome fasta
REFGENOME = config["REFGENOME"]["fasta"]
REF_FAI = config["REFGENOME"]["fai"]
## cutadapt
P2_adapter_index = config["FILTER"]["P2_adapter_index"]
#P2_adapter_sequence=config["FILTER"]["P2_adapter_index"][P2_adapter_index]
TRIM_R1_3P = config["FILTER"]["P2_adapter_sequence"][3]
TRIM_R2_3P_FILE = config["FILTER"]["cutadapt"]["trim_R2_3prime_file"]
TRIM_R2_3P_DF = pd.read_csv(TRIM_R2_3P_FILE, header = None, sep = " ") # read the reverse complemented 96 adapter1 sequences

## bowtie2
REFERENCE = config["MAPPING"]["reference"]
REFBASE = os.path.basename(REFERENCE)
## markers to use
BCF_TARGET = config["MARKER"]["bcf_target"]
MARKCOMP = config["MARKER"]["complete"]
MARKMASK = config["MARKER"]["mask"]
## TSV
binSize = config["TSV"]["binSize"]
binName = config["TSV"]["binName"]
## sources
TIGERIN_src = config["SCRIPTS"]["tigercall2tigerin"]
TIGER_src = config["SCRIPTS"]["tiger"]
RUNTIGER_src = config["SCRIPTS"]["runtiger"]
COTABLE_src = config["SCRIPTS"]["cotable"]
TSV_src = config["SCRIPTS"]["toTSV"]


# to be deleted
# BARCODE = glob_wildcards("raw/{libnum}_{barcode}_{read_num}.fastq.gz").barcode
# =======


# Specify the desired end target file(s)
rule all:
    input:
       expand("results/01_trimmed/{sample}_R{read_num}.tr.fastq.gz", sample = SAMPLE, read_num = [1, 2]),
       expand("results/02_bowtie2/lib{sample}_MappedOn_{refbase}.bam",
               sample = SAMPLE,
               refbase = REFBASE),
       expand("results/02_bowtie2/filtered/lib{sample}_MappedOn_{refbase}_sort.bam",
              sample = SAMPLE,
              refbase = REFBASE),
       expand("results/02_bowtie2/filtered/lib{sample}_MappedOn_{refbase}_sort.md.bam",
           sample = SAMPLE,
           refbase = REFBASE),
       expand("results/03_vcf/lib{sample}_MappedOn_{refbase}.vcf.gz",
           sample = SAMPLE,
           refbase = REFBASE),
       expand("results/03_vcf/lib{sample}_MappedOn_{refbase}.vcf.gz.csi",
           sample = SAMPLE,
           refbase = REFBASE),
       expand("results/04_tigercall/lib{sample}_MappedOn_{refbase}.tigercalls.txt",
           sample = SAMPLE,
           refbase = REFBASE),
       expand("results/05_tigerin/lib{sample}_MappedOn_{refbase}.mask.txt",
           sample = SAMPLE,
           refbase = REFBASE),
       expand("results/06_tiger/lib{sample}_MappedOn_{refbase}.smooth.co.txt",
           sample = SAMPLE,
           refbase = REFBASE),
       expand("results/{lib_name}_cotable.txt",
           lib_name = LIB_NAME),
        expand("results/{lib_name}_genomeBin{binName}.tsv", lib_name = LIB_NAME, binName = binName)

# Run fastqc on single-end raw data
# Trim off adapters

## Get reverse complemented adapter sequence for trimming 3' end of Read 2
## https://snakemake.readthedocs.io/en/stable/tutorial/advanced.html#step-3-input-functions
def getadapt(wildcards):
    return(TRIM_R2_3P_DF.iloc[int(wildcards.sample)-1][1])

rule cutadapt:
    """Remove adapters"""
    output:
        tr_read1 = temp("results/01_trimmed/{sample}_R1.tr.fastq.gz"),
        tr_read2 = temp("results/01_trimmed/{sample}_R2.tr.fastq.gz"),
        qc    = "qc/cutadapt/{sample}_cutadapt.qc.txt"
    # input:
    #     read1 = "data/fastq/{sample}_{barcode}_1.fastq.gz",
    #     read2 = "data/fastq/{sample}_{barcode}_2.fastq.gz"
    resources:
        tmpdir = "tmp"
    params:
        trim_R1_3prime = TRIM_R1_3P,
        trim_R2_3prime = getadapt,
        cut_R1_5prime = config["FILTER"]["cutadapt"]["cut_R1_5prime"],
        cut_R2_5prime = config["FILTER"]["cutadapt"]["cut_R2_5prime"],
        quality_filter = config["FILTER"]["cutadapt"]["quality-filter"],
        minimum_overlap = config["FILTER"]["cutadapt"]["minimum-overlap"],
        P2_index = P2_adapter_index
    log:
        "logs/cutadapt/{sample}_trimmed.log"
    shell:
        "cutadapt -u {params.cut_R1_5prime}"
        " -U {params.cut_R2_5prime}"
        " -a {params.trim_R1_3prime}"
        " -A {params.trim_R2_3prime}"
        " -O {params.minimum_overlap}"
        " -q {params.quality_filter}"
        # " --cores=0"
        " raw/P1_{wildcards.sample}_P2_{params.P2_index}_1.fastq.gz raw/P1_{wildcards.sample}_P2_{params.P2_index}_2.fastq.gz"
        " -o {output.tr_read1}"
        " -p {output.tr_read2}"
        " > {output.qc} &> {log}"

# Align to reference genome
rule bowtie2:
    """Map reads using bowtie2 and filter alignments using samtools"""
    # output: temp("results/02_bowtie2/lib{sample}_MappedOn_{refbase}.bam")
    output: temp("results/02_bowtie2/lib{sample}_MappedOn_{refbase}.bam")
    input:
        tr_1 = "results/01_trimmed/{sample}_R1.tr.fastq.gz",
        tr_2 = "results/01_trimmed/{sample}_R2.tr.fastq.gz"
    params:
        MAPQmaxi = config["MAPPING"]["MAPQmaxi"]
    threads: config["THREADS"]
    log:
        "logs/bowtie2/lib{sample}_MappedOn_{refbase}.log"
    shell:
        # -F 2308 excludes unmapped reads,
        # as well as secondary and supplementary alignments
        "(bowtie2 --very-sensitive"
        " --threads {threads}"
        " -x {REFERENCE} -1 {input.tr_1} -2 {input.tr_2}"
        " | samtools view -bh -@ {threads} -F 2308 -o {output} - ) &> {log}"

# Filter out multireads by using MAPQscore
rule samtools:
    output: temp("results/02_bowtie2/filtered/lib{sample}_MappedOn_{refbase}_sort.bam")
    input: "results/02_bowtie2/lib{sample}_MappedOn_{refbase}.bam"
    params:
        sortMemory = config["MAPPING"]["sortMemory"],
        MAPQmaxi = config["MAPPING"]["MAPQmaxi"]
    threads: config["THREADS"]
    log: "logs/samtools/lib{sample}_MappedOn_{refbase}_sort.log"
    shell:
        "(samtools view -h {input} -q {params.MAPQmaxi} -u "
        "| samtools sort -@ {threads} -m {params.sortMemory} -o {output} -) &> {log};"
        # http://biofinysics.blogspot.com/2014/05/how-does-bowtie2-assign-mapq-scores.html
        # https://sequencing.qcfail.com/articles/mapq-values-are-really-useful-but-their-implementation-is-a-mess/

rule markdup:
    output: 
        bam="results/02_bowtie2/filtered/lib{sample}_MappedOn_{refbase}_sort.md.bam",
        metric="logs/markdup/lib{sample}_MappedOn_{refbase}_sort.md.txt",
        index = "results/02_bowtie2/filtered/lib{sample}_MappedOn_{refbase}_sort.md.bam.bai"
    resources:
        tmpdir = "tmp"
    input: "results/02_bowtie2/filtered/lib{sample}_MappedOn_{refbase}_sort.bam"
    # threads: config["THREADS"]
    shell:
        "picard MarkDuplicates -I {input}"
        " -O {output.bam}"
        " -M {output.metric}"
        " --TMP_DIR {resources.tmpdir}"
        " --REMOVE_DUPLICATES true;"
        "samtools index {output.bam} -o {output.index}"


# variant calling (VCF)
rule varcall:
    output: 
        vcf = "results/03_vcf/lib{sample}_MappedOn_{refbase}.vcf.gz",
    input: "results/02_bowtie2/filtered/lib{sample}_MappedOn_{refbase}_sort.md.bam"
    shell:
        r"""
        bcftools mpileup -T {BCF_TARGET} -C {BCF_TARGET} -f {REFGENOME} {input} -Ou |
        bcftools call -m -T {BCF_TARGET} -Oz -o {output.vcf} -
        """
rule index_vcf:
    output:
        index = "results/03_vcf/lib{sample}_MappedOn_{refbase}.vcf.gz.csi"
    input: "results/03_vcf/lib{sample}_MappedOn_{refbase}.vcf.gz"
    shell:
        "bcftools index {input} -o {output.index}"

rule vcf2tigercall:
    output: "results/04_tigercall/lib{sample}_MappedOn_{refbase}.tigercalls.txt"
    input: "results/03_vcf/lib{sample}_MappedOn_{refbase}.vcf.gz"
    shell:
        r"""
        bcftools query -f '%CHROM %POS %REF %ALT %QUAL [ %INDEL %DP %DP4]\n' {input} |
        awk 'BEGIN{{IFS="\t"; OFS="\t"}}\
            {{split($4, alt, ","); split($8, dp4, ",")}}\
            length(alt)==1 && $1 != "ChrC" && $1 != "ChrM" && $6 != 1\
            {{print $1, $2, $3, dp4[1]+dp4[2], $4, dp4[3]+dp4[4]}}' > {output}
        """
rule tigercall2tigerin:
    output: 
        complete="results/05_tigerin/lib{sample}_MappedOn_{refbase}.complete.txt",
        mask="results/05_tigerin/lib{sample}_MappedOn_{refbase}.mask.txt"
    input:
        "results/04_tigercall/lib{sample}_MappedOn_{refbase}.tigercalls.txt"
    params:
        markcomp = MARKCOMP,
        markmask = MARKMASK
    shell:
        r"""
        mkdir -p tmp/tigerin_temp;
        Rscript {TIGERIN_src} {input} tmp/tigerin_temp lib{wildcards.sample} {params.markcomp} {params.markmask}
        mv tmp/tigerin_temp/lib{wildcards.sample}.mask.txt {output.mask}
        mv tmp/tigerin_temp/lib{wildcards.sample}.complete.txt {output.complete}
        """

rule tiger:
    output: "results/06_tiger/lib{sample}_MappedOn_{refbase}.smooth.co.txt"
    input:
        complete = "results/05_tigerin/lib{sample}_MappedOn_{refbase}.complete.txt",
        mask= "results/05_tigerin/lib{sample}_MappedOn_{refbase}.mask.txt"
    log:
        "logs/tiger/{sample}_{refbase}.tiger.log"
    shell:
        r"""
        bash {RUNTIGER_src} {TIGER_src} {input.complete} {input.mask} results/06_tiger &> {log}
        """
        # mkdir -p results/06_tiger
rule cotable:
    output: "results/{lib_name}_cotable.txt"
    input:
        input_files = expand("results/06_tiger/lib{sample}_MappedOn_{refbase}.smooth.co.txt", sample=SAMPLE, refbase=REFBASE)
    shell:
        r"""
        Rscript {COTABLE_src} results/06_tiger {output} {SAMPLE_SIZE}
        """
rule toGenomeBin:
    output: "results/{lib_name}_genomeBin{binName}.tsv"
    input: "results/{lib_name}_cotable.txt"
    shell:
        r"""
        Rscript {TSV_src} {input} {binSize} {REF_FAI} {output}
        """
