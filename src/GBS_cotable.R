# Title : GBS_cotable.R
# Maintainer : Namil Son
# Date of written : Mar 1st, 2021
# Description :
# This script reads *.smooth.co.txt files within the "tiger_out" directory, then summarize into one table.
# Usage : Rscript GBS_cotable.R <full/path/to/tiger_out> <full/apath/to/output/directory> <group name(string between "lib number." and ".smooth.co.txt" in "*.smooth.co.txt" file> <#lib>
# e.g) Rscript GBS_cotable.R /home/work/GBS/Col-0/tiger_out /home/work/GBS/Col-0 Col-0 96

library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)

dirin <- args[1]
output_name <- args[2]
#dirout <- args[2]
nlib <- args[3]
#outname <- args[4]


#cotable() : read TIGER smooth results in bulk, then summarize into a crossover table
cotable <- function(lib.nums, dirin){
		lib.nums <- as.numeric(lib.nums)
        co.all <- NULL
		files=list.files(dirin)
		files.smooth=files[grep("smooth", files)]
        for(k in 1:lib.nums){
                smooth=paste0(dirin, "/", files.smooth[k])
                dat <- read.table(file=smooth, header=F)
                all <- NULL
                for(i in 1:5){
                        chr <- dat[which(dat[,2]==i),]
                        if((length(chr[,1])>1)==T){
                                start <- chr[,4]
                                stop <- chr[,3]
                                stop <- stop[-1]
                                start <- start[-length(start)]
                                diff <- stop-start
                                mid <- round(diff/2)
                                cos <- start + mid
                                chrs <- rep(i, length(cos))
                                lib <- rep(dat[1,1], length(chrs))
                                width <- stop-start
                                bind <- cbind(lib, chrs, start, stop, cos, width)
                                all <- rbind(all, bind)
                                print(dim(all))
                        }
                }
                co.all <- rbind(co.all, all)
        }
        co.all <- as_tibble(co.all)
        co.all <- co.all %>%
                mutate_at(vars(chrs), funs(paste0("Chr",.)))
        return(co.all)
}

output <- cotable(nlib, dirin)
write_csv(output, file=output_name)
