# Title : tigercalls2tiger.R
# Maintainer : Namil Son
# First Edition : 12. 24. 2020
# Last modified : 12. 24. 2020
# Description :
#  This script takes *.tigercalls.txt file, which is a modified form of vcf, and outputs to two tiger input files. 1) tiger input with complete Ler marker; 2) tiger input with masked Ler marker.
# Usage :
# tigercalls2tiger.R <input> <outdir> <prefix>
# input and outdir should be the full path to the file.

suppressMessages({
        library(tidyverse)
})

#---- Prepare inputs --------
args <- commandArgs(trailingOnly=TRUE)

input <- args[1]
outdir <- args[2]
prefix <- args[3]
file.complete <- args[4]
file.mask <- args[5]

#input <- "/datasets/data_4/nison/GBS/GBS_marker_v2/2019MMDD_GBS_Col/results/04_tigercall/lib16_MappedOn_tair10.tigercalls.txt"
#file.complete <- "/datasets/data_4/nison/GBS/GBS_marker_v2/2019MMDD_GBS_Col/data/ColLer_markers/coller_complete_v2.tsv"
#file.mask <- "/datasets/data_4/nison/GBS/GBS_marker_v2/2019MMDD_GBS_Col/data/ColLer_markers/coller_gbs_marker_v2.tsv"

colname=c("Chr", "Pos", "Ref", "ref.count", "Alt", "alt.count")

# tigercalls.txt
tigercall <- read_delim(file=input, delim="\t", col_names=colname)
# complete marker set
complete <- read_delim(file=file.complete, delim="\t", col_names=c("Chr", "Pos", "Ref", "Alt"))
# corrected marker set (hqsnv)
mask <- read_delim(file=file.mask, delim="\t", col_names=c("Chr", "Pos", "Ref", "Alt"))

#---- make tiger input files --------
# {onlynum} : edit "Chr" column from "Chr_" to "_"
onlynum <- function(x){
        y=substr(x, 4, 4)
        return(y)
}
# Update complete marker set with detected SNPs
complete <- complete %>%
        left_join(tigercall[,-c(3,5)], by=c("Chr", "Pos")) %>%
        replace_na(list(ref.count=0, alt.count=0)) %>%
        relocate(ref.count, .after=Ref) %>%
        mutate_at(vars(Chr), funs(onlynum))

# Update corrected marker set with detected SNPs
mask <- mask %>%
        left_join(tigercall[, -c(3,5)], by=c("Chr", "Pos")) %>%
        replace_na(list(ref.count=0, alt.count=0)) %>%
        relocate(ref.count, .after=Ref) %>%
        mutate_at(vars(Chr), funs(onlynum))

write_tsv(complete, file=paste0(outdir, "/", prefix, ".complete.txt"), col_names=FALSE)
write_tsv(mask, file=paste0(outdir, "/", prefix, ".mask.txt"), col_names=FALSE)
