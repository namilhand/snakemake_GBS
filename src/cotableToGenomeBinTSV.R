library(tidyverse)
library(parallel)

args <- commandArgs(trailingOnly=T)

input <- args[1]
binSize <- as.numeric(args[2])
fai <- args[3]
output <- args[4]


#input <- "/datasets/data_4/nison/GBS/20230512_GBS_MIGS-H2AW-SUVH/results/20230512_MIGS-H2AW-SUVH_cotable.txt"
#binSize <- 10000
#fai <- "/home/nison/work/refgenome/TAIR10/TAIR10.fasta.fai"

cotable <- read_csv(input)
fai <- read_tsv(fai, col_names=F)
genomeSize <- fai$X2[1:5]
cumGenomeSize <- c(0, cumsum(genomeSize))



# makeWindow(): slice chromosome by window size and record coordinate and cumulative coordinate of each slice
makeWindow <- function(winsize){
        wg.bin <- tibble(chrs=character(), bin.start=numeric(), bin.end=numeric(), cum.start=numeric(), cum.end=numeric())
        for(i in 1:5){
                coord.start <- seq(1, genomeSize[i], by=winsize)
                coord.end <- c(coord.start[-1]-1, genomeSize[i])
                cum.start <- coord.start+cumGenomeSize[i]
                cum.end <- coord.end + cumGenomeSize[i] 
                nchr <- rep(i, length(coord.start))
                bins <- tibble(chrs=paste0("Chr",nchr), bin.start=coord.start, bin.end=coord.end, cum.start=cum.start, cum.end=cum.end)
                wg.bin <- add_row(wg.bin, bins)
        }
        return(wg.bin)
}

genomeBin <- makeWindow(binSize)

# count COs within each window
binCotable <- function(dat){
        binned <- NULL
        collect.bin <- NULL
        for(i in 1:5){
                chr.dat <- dat %>%
                        filter(chrs==paste0("Chr",i))
                chr.bin <- genomeBin %>%
                        filter(chrs==paste0("Chr",i))
                for(j in 1:nrow(chr.bin)){
                        dat.in.bin <- sum(chr.dat[["cos"]] >= chr.bin$bin.start[j] & chr.dat[["cos"]] <= chr.bin$bin.end[j])
                        collect.bin <- c(collect.bin, dat.in.bin)
                }
                
        }
        coBin <- tibble(genomeBin, coInBin=collect.bin, libSize=length(unique(dat$lib)))

        resultList <- mclapply(1:5, function(x){
            chrCoBin <- filter(coBin, chrs == paste0("Chr", x))
            if(with(chrCoBin[nrow(chrCoBin),], bin.end - bin.start) + 1 < binSize){
                chrCoBin[nrow(chrCoBin),]$coInBin <- chrCoBin[nrow(chrCoBin)-1,]$coInBin
            }
            return(chrCoBin)
        }, mc.cores = 5)

        result <- bind_rows(resultList)
        return(result)
}
genomeBinCotable <- binCotable(cotable)
genomeBinCotableTSV <- genomeBinCotable %>%
    dplyr::select(-c(bin.end, cum.end))
colnames(genomeBinCotableTSV) <- c("chr", "window", "cumwindow", "coInWindow", "libSize")
write_tsv(genomeBinCotableTSV, file=output)