# Snakemake workflow for TIGER analysis

# Usage (snakemake --cores should reflect available cores):
# conda env create --file environment.yaml --name GBS 
# conda activate GBS
# snakemake -p --cores 30 
# conda deactivate

import pandas as pd
import os

# To make the per_base_coverage rule work with a shell script invoked using the "shell" directive,
# we need to determine the base path of Snakefile since we Expect the scripts directory to be there as well
SRCDIR = srcdir("")

# Specify config file parameters
configfile: "config.yaml"
SAMPLE = list(range(1, 97))
REFERENCE = config["MAPPING"]["reference"]
REFBASE = os.path.basename(REFERENCE)
BARCODE = glob_wildcards("raw/{libnum}_{barcode}_{read_num}.fastq.gz").barcode
ADAPTER_R2_3P = config["FILTER"]["cutadapt"]["adapter_R2_3prime_file"]
MARKCOMP_GZ = config["MARKER"]["complete_gz"]
MARKCOMP = config["MARKER"]["complete"]
MARKMASK = config["MARKER"]["mask"]
REFGENOME = config["REFGENOME"]
TIGERIN_src = config["SCRIPTS"]["tigercall2tigerin"]
TIGER_src = config["SCRIPTS"]["tiger"]
COTABLE_src = config["SCRIPTS"]["cotable"]
N_LIB = config["COTABLE"]["num_lib"]
LIB_NAME = config["COTABLE"]["lib_name"]

adapter_R2_3prime_df = pd.read_csv(ADAPTER_R2_3P, header = None, sep = " ")

# Specify the desired end target file(s)
rule all:
    input:
        # expand("qc/fastqc/{sample}_{barcode}_{read_num}_fastqc.html", sample = [10,11,12], barcode = BARCODE[0:2], read_num = [1, 2]),
        expand("results/01_trimmed/{sample}_R{read_num}.tr.fastq.gz", sample = SAMPLE, read_num = [1, 2]),
        expand("results/02_bowtie2/lib{sample}_MappedOn_{refbase}.bam",
                 sample = SAMPLE,
                 refbase = REFBASE),
        expand("results/02_bowtie2/filtered/mapq_{mapq}/lib{sample}_MappedOn_{refbase}_mapq{mapq}_sort.bam",
               sample = SAMPLE,
               refbase = REFBASE,
               mapq = [5]),
        expand("results/02_bowtie2/filtered/mapq_{mapq}/lib{sample}_MappedOn_{refbase}_mapq{mapq}_sort.md.bam",
            sample = SAMPLE,
            refbase = REFBASE,
            mapq = [5]),
        expand("results/03_vcf/lib{sample}_MappedOn_{refbase}_mapq{mapq}.vcf.gz",
            sample = SAMPLE,
            refbase = REFBASE,
            mapq = [5]),
        expand("results/04_tigercall/lib{sample}_MappedOn_{refbase}_mapq{mapq}.tigercalls.txt",
            sample = SAMPLE,
            refbase = REFBASE,
            mapq = [5]),
        expand("results/05_tigerin/lib{sample}_MappedOn_{refbase}_mapq{mapq}.mask.txt",
            sample = SAMPLE,
            refbase = REFBASE,
            mapq = [5]),
        expand("results/06_tiger/mapq_{mapq}/lib{sample}_MappedOn_{refbase}_mapq{mapq}.smooth.co.txt",
            sample = SAMPLE,
            refbase = REFBASE,
            mapq = [5])
        # expand("results/{lib_name}_mapq{mapq}_cotable.txt",
        #     lib_name = LIB_NAME,
        #     mapq = [5])

# Run fastqc on single-end raw data
# rule fastqc_raw:
#     """Create fastqc report"""
#     output:
#         html = "qc/fastqc/{sample}_{barcode}_{read_num}_fastqc.html",
#         zip  = "qc/fastqc/{sample}_{barcode}_{read_num}_fastqc.zip"
#     input:
#         "raw/{sample}_{barcode}_{read_num}.fastq.gz"
#     params:
#         " --extract" +
#         " --adapters " + str(config["FILTER"]["fastqc"]["adapters"])
#     log:
#         "logs/fastqc/{sample}_{barcode}_{read_num}.log"
#     wrapper:
#         "v1.28.0/bio/fastqc"

# Trim off adapters

## Get Read2 3'end adapter sequence
## https://snakemake.readthedocs.io/en/stable/tutorial/advanced.html#step-3-input-functions
def getadapt(wildcards):
    return(adapter_R2_3prime_df.iloc[int(wildcards.sample)-1][1])

rule cutadapt:
    """Remove adapters"""
    output:
        tr_read1 = "results/01_trimmed/{sample}_R1.tr.fastq.gz",
        tr_read2 = "results/01_trimmed/{sample}_R2.tr.fastq.gz",
        qc    = "qc/cutadapt/{sample}_cutadapt.qc.txt"
    # input:
    #     read1 = "data/fastq/{sample}_{barcode}_1.fastq.gz",
    #     read2 = "data/fastq/{sample}_{barcode}_2.fastq.gz"
    params:
        adapter_R1_3prime = config["FILTER"]["cutadapt"]["adapter_R1_3prime"],
        adapter_R2_3prime = getadapt,
        cut_R1_5prime = config["FILTER"]["cutadapt"]["cut_R1_5prime"],
        quality_filter = config["FILTER"]["cutadapt"]["quality-filter"],
        minimum_overlap = config["FILTER"]["cutadapt"]["minimum-overlap"]
    log:
        "logs/cutadapt/{sample}_trimmed.log"
    shell:
        "cutadapt -u {params.cut_R1_5prime}"
        " -a {params.adapter_R1_3prime}"
        " -A {params.adapter_R2_3prime}"
        " -O {params.minimum_overlap}"
        " -q {params.quality_filter}"
        " --cores=0"
        " raw/{wildcards.sample}_*_1.fastq.gz raw/{wildcards.sample}_*_2.fastq.gz"
        " -o {output.tr_read1}"
        " -p {output.tr_read2}"
        " > {output.qc} 2> {log}"

# Align to reference genome
rule bowtie2:
    """Map reads using bowtie2 and filter alignments using samtools"""
    output: temp("results/02_bowtie2/lib{sample}_MappedOn_{refbase}.bam")
    input:
        # fastq = "data/dedup/trimmed/{sample}_dedup_trimmed.fastq.gz",
        tr_1 = "results/01_trimmed/{sample}_R1.tr.fastq.gz",
        tr_2 = "results/01_trimmed/{sample}_R2.tr.fastq.gz"
    params:
        # alignments = config["MAPPING"]["alignments"],
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
        " | samtools view -bh -@ {threads} -F 2308 -o {output} - ) 2> {log}"

# Filter alignments for mismatches and extract unique alignments
rule samtools:
    output: "results/02_bowtie2/filtered/mapq_{mapq}/lib{sample}_MappedOn_{refbase}_mapq{mapq}_sort.bam"
        # mapq10   = "results/bowtie2/filtered/{sample}_{barcode}_MappedOn_{refbase}_mapq10_sort.bam",
        # mapq5   = "results/bowtie2/filtered/{sample}_{barcode}_MappedOn_{refbase}_mapq5_sort.bam",
        # mapq0   = "results/bowtie2/filtered/{sample}_{barcode}_MappedOn_{refbase}_mapq0_sort.bam"
    input: "results/02_bowtie2/lib{sample}_MappedOn_{refbase}.bam"
    params:
        sortMemory = config["MAPPING"]["sortMemory"]
    threads: config["THREADS"]
    log: "logs/samtools/lib{sample}_MappedOn_{refbase}_mapq{mapq}_sort.log"
        # mapq10   = "logs/samtools/{sample}_{barcode}_MappedOn_{refbase}_mapq10_sort.log",
        # mapq5 = "logs/samtools/{sample}_{barcode}_MappedOn_{refbase}_mapq5_sort.log",
        # mapq0 = "logs/samtools/{sample}_{barcode}_MappedOn_{refbase}_mapq0_sort.log"
    shell:
        "(samtools view -h {input} -q {wildcards.mapq} -u "
        "| samtools sort -@ {threads} -m {params.sortMemory} -o {output} -) 2> {log};"
        # "(samtools view -h {input} -q 5 -u "
        # "| samtools sort -@ {threads} -m {params.sortMemory} -o {output.mapq5} -) 2> {log.mapq5};"
        # "(samtools view -h {input} -q 0 -u "
        # "| samtools sort -@ {threads} -m {params.sortMemory} -o {output.mapq0} -) 2> {log.mapq0}"
        # Extract unique alignments, excluding alignments with MAPQ scores < config["MAPPING"]["MAPQunique"]
        # http://biofinysics.blogspot.com/2014/05/how-does-bowtie2-assign-mapq-scores.html
        # https://sequencing.qcfail.com/articles/mapq-values-are-really-useful-but-their-implementation-is-a-mess/
        # "(samtools view -h -q {params.MAPQunique} {input} "
        # "| grep -e '^@' -e 'XM:i:[012][^0-9]' "
        # "| samtools view -u - "
        # "| samtools sort -@ {threads} -m {params.sortMemory} - "
        # "| samtools rmdup -s - {output.unique}) 2> {log.unique}"

rule markdup:
    output: 
        bam="results/02_bowtie2/filtered/mapq_{mapq}/lib{sample}_MappedOn_{refbase}_mapq{mapq}_sort.md.bam",
        metric="logs/markdup/lib{sample}_MappedOn_{refbase}_mapq{mapq}_sort.md.txt",
        index = "results/02_bowtie2/filtered/mapq_{mapq}/lib{sample}_MappedOn_{refbase}_mapq{mapq}_sort.md.bam.bai"
    input: "results/02_bowtie2/filtered/mapq_{mapq}/lib{sample}_MappedOn_{refbase}_mapq{mapq}_sort.bam"
    threads: config["THREADS"]
    shell:
        "picard MarkDuplicates -I {input}"
        " -O {output.bam}"
        " -M {output.metric}"
        " --REMOVE_DUPLICATES true;"
        "samtools index {output.bam} -o {output.index}"

# variant calling (VCF)
rule varcall:
    output: 
        vcf = "results/03_vcf/lib{sample}_MappedOn_{refbase}_mapq{mapq}.vcf.gz",
        index = "results/03_vcf/lib{sample}_MappedOn_{refbase}_mapq{mapq}.vcf.gz.csi"
    input: "results/02_bowtie2/filtered/mapq_{mapq}/lib{sample}_MappedOn_{refbase}_mapq{mapq}_sort.md.bam"
    shell:
        r"""
        bcftools mpileup -f {REFGENOME} {input} -Ou |
        bcftools call -m -T {MARKCOMP_GZ} -Oz -o {output.vcf} -

        bcftools index {output.vcf} -o {output.index}
        """
rule vcf2tigercall:
    output: "results/04_tigercall/lib{sample}_MappedOn_{refbase}_mapq{mapq}.tigercalls.txt"
    input: "results/03_vcf/lib{sample}_MappedOn_{refbase}_mapq{mapq}.vcf.gz"
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
        complete="results/05_tigerin/lib{sample}_MappedOn_{refbase}_mapq{mapq}.complete.txt",
        mask="results/05_tigerin/lib{sample}_MappedOn_{refbase}_mapq{mapq}.mask.txt"
    input:
        "results/04_tigercall/lib{sample}_MappedOn_{refbase}_mapq{mapq}.tigercalls.txt"
    shell:
        r"""
        mkdir -p tmp/tigerin_temp;
        Rscript {TIGERIN_src} {input} tmp/tigerin_temp lib{wildcards.sample}_mapq{wildcards.mapq}
        mv tmp/tigerin_temp/lib{wildcards.sample}_mapq{wildcards.mapq}.mask.txt {output.mask}
        mv tmp/tigerin_temp/lib{wildcards.sample}_mapq{wildcards.mapq}.complete.txt {output.complete}
        """

rule tiger:
    output: "results/06_tiger/mapq_{mapq}/lib{sample}_MappedOn_{refbase}_mapq{mapq}.smooth.co.txt"
    input:
        complete = "results/05_tigerin/lib{sample}_MappedOn_{refbase}_mapq{mapq}.complete.txt",
        mask= "results/05_tigerin/lib{sample}_MappedOn_{refbase}_mapq{mapq}.mask.txt"
    shell:
        r"""
        bash {TIGER_src} {input.complete} {input.mask} results/06_tiger/mapq_{wildcards.mapq}
        """
        # mkdir -p results/05_tiger
# rule cotable:
#     output: "results/{lib_name}_mapq{mapq}_cotable.txt"
#     shell:
#         r"""
#         Rscript {COTABLE_src} results/06_tiger/mapq_{wildcards.mapq} {output} {N_LIB}
#         """