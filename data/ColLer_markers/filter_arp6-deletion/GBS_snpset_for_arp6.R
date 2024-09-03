library(tidyverse)

dirin <- "/home/nison/work/pipelines/GBS/ColLer_markers"
dirout <- "/home/nison/work/pipelines/GBS/ColLer_markers/filter_arp6-deletion"

full <- file.path(dirin, "BC.complete.tsv")
mask <- file.path(dirin, "BC.mask.tsv")
deletion <- file.path(dirin, "filter_arp6-deletion", "arp6_deletion.txt")

full <- read_tsv(full, col_names=F)
mask <- read_tsv(mask, col_names=F)

black <- read_tsv(deletion, col_names=F, skip=1)
# arp6 blacklist region where deleted by proton-induced mutagenesis

full.filt <- full %>%
		filter(!((X1 == "Chr3") & X2 >= black$X2[1] & X2 <= black$X3[1]))

mask.filt <- mask %>%
		filter(!((X1 == "Chr3") & X2 >= black$X2[1] & X2 <= black$X3[1]))

write_tsv(full.filt, file.path(dirout, "BC.complete.arp6_filter.tsv"), col_names=F)
write_tsv(mask.filt, file.path(dirout, "BC.mask.arp6_filter.tsv"), col_names=F)
