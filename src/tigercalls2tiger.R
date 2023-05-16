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

colname=c("Chr", "Pos", "Ref", "ref.count", "Alt", "alt.count")

# tigercalls.txt
tigercall <- read_delim(file=input, delim="\t", col_names=colname)
# complete marker set
complete <- read_delim(file="/usr/local/src/GBS_suite/ColLer_markers/BC.complete.tsv", delim="\t", col_names=colname)
# corrected marker set (hqsnv)
mask <- read_delim(file="/usr/local/src/GBS_suite/ColLer_markers/BC.mask.tsv", delim="\t", col_names=colname)

#---- make tiger input files --------
# {onlynum} : edit "Chr" column from "Chr_" to "_"
onlynum <- function(x){
        y=substr(x, 4, 4)
        return(y)
}
# Update complete marker set with detected SNPs
complete <- complete %>%
        dplyr::select(-c("ref.count", "alt.count")) %>%
        left_join(tigercall[,-c(3,5)], by=c("Chr", "Pos")) %>%
        replace_na(list(ref.count=0, alt.count=0)) %>%
        relocate(ref.count, .after=Ref) %>%
        mutate_at(vars(Chr), funs(onlynum))

# Update corrected marker set with detected SNPs
mask <- mask %>%
        dplyr::select(-c("ref.count", "alt.count")) %>%
        left_join(tigercall[,-c(3,5)], by=c("Chr", "Pos")) %>%
        replace_na(list(ref.count=0, alt.count=0)) %>%
        relocate(ref.count, .after=Ref) %>%
        mutate_at(vars(Chr), funs(onlynum))

write_tsv(complete, file=paste0(outdir, "/", prefix, ".complete.txt"), col_names=FALSE)
write_tsv(mask, file=paste0(outdir, "/", prefix, ".mask.txt"), col_names=FALSE)
