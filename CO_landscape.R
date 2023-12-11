# Description 
library(tidyverse)
library(scales)

# GENOME INFO SETTING ====================================
#' 0. chromosome info 
centromeres <- c(15086045,3607929,13587786,3956021,11725024)
chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")
chr.ends <- c(30427671,19698289,23459830,18585056,26975502)

tha.cum <- c(0, cumsum(chr.ends))
tha.tot <- tha.cum[length(tha.cum)]
centromeres.cum <- centromeres + tha.cum[1:5]
# pericentromeric + centromeric region
# (Ref.: Underwood, 2018, Genome Research)
north.start <- c(11420001, 910001, 10390001, 1070001, 8890001)
south.end <- c(18270000, 7320000, 16730000, 6630000, 15550000)

# INPUT DATA =============================================================

### assuming using snakemake
dir_gbs <- "data/tsv/gbs"
dir_dname <- "data/tsv/dname"

args <- commandArgs(trailingOnly=T)
dat_test <- args[1]
libname <- args[2]

###

## GBS landscape
col_1 <- read_tsv(file.path(dir_gbs, "WT_2019_genomeBin100kb.tsv"))
col_2 <- read_tsv(file.path(dir_gbs, "WT_2020_genomeBin100kb.tsv"))
col_3 <- read_tsv(file.path(dir_gbs, "WT_2021_genomeBin100kb.tsv"))
dat_test <- read_tsv(file.path("results", paste0(libname, "_genomeBin100kb.tsv")))

## dname landscape
mC_col <- read_tsv(file.path(dir_dname, "1-1_C_TAIR10_window100kb_step100kb.tsv"), col_names=T)
colnames(mC_col) <- c("chr", "window", "cumwindow", "mC.ratio")

## merge GBS results from the same genotype
col_tsv <- dplyr::select(col_1, -c(coInWindow, libSize)) %>%
    add_column(coInWindow = col_1$coInWindow + col_2$coInWindow + col_3$coInWindow) %>%
    add_column(libSize = col_1$libSize + col_2$libSize + col_3$libSize) %>%
    mutate(mean.coInWindow = coInWindow/libSize)


test_tsv <- dat_test %>% 
    mutate(mean.coInWindow = coInWindow/libSize)

# co_tsv.list <- list(wt=col_tsv, h2aw=h2aw_tsv, suvh=suvh_tsv, dmigs=dmigs_tsv, cmt3=cmt3_tsv)
co_tsv.list <- list(wt=col_tsv, test=test_tsv)
names(co_tsv.list) <- c("wt", libname)


# cumulative coordinates of pericentromere
coord_pericen <- tibble(chr=1:5, north=north.start + tha.cum[1:5], south=south.end + tha.cum[1:5])

# MA smoothing
# k = 5 is recommended for drawing ChIP landscape

ma_width.co <- 7
ma_width.chip <- 5

mafilter <- function(dat, valueColumn, k){
        # k = arm width for MA smoothing
        filt <- 1/(2*k+1)
        filt <- rep(filt, 2*k+1)
        filt.collect <- NULL

        for(i in 1:5){
                chr <- dat %>%
                        filter(chr==paste0("Chr",i))
                filt.chr <- stats::filter(chr[[valueColumn]], filt)
                filt.chr[1:k] <- filt.chr[k+1]
                len <- length(filt.chr)
                filt.chr[(len-k+1):len] <- filt.chr[(len-k)]

                filt.collect <- c(filt.collect, filt.chr)
        }
        fin <- tibble(dat, smooth=filt.collect)
        return(fin)
}

# moving average genomeBin100kb
co_ma.list <- lapply(co_tsv.list, mafilter, "mean.coInWindow", ma_width.co)
co_ma <- bind_rows(co_ma.list, .id="genotype")
co_ma <- mutate(co_ma, genotype = factor(genotype, levels=c("wt", libname)))  %>%
    mutate(cMMb = smooth *10 * 100)
    

# palette
# pal_GBS <- colour("bright")(5)
pal_GBS <- c("wt"="blue", "test"="red")
names(pal_GBS) <- c("wt", libname)

# co landscape

## y-axis scale
landscape_ylim <- c(0, max(co_ma$cMMb)*1.05) 

drawCoLandscape <- function(dat, pal, prefix, plot_width, plot_height){

    ## convert CO per F2 per 100 kb into cM/Mb
    ## *10 : CO per F2 per 100kb --> CO per F2 per Mb
    ## *100: CO per F2 per Mb (M/Mb) --> cM/Mb

    # dat_cMMb <- mutate(dat, cMMb = smooth * 10 * 100 )


#    pal <- pal_GBS
#    dat <- co_ma
    p <- ggplot() +
        annotate("rect", xmin=coord_pericen$north[1], xmax=coord_pericen$south[1], ymin=-Inf, ymax=Inf, fill="grey90") +
        annotate("rect", xmin=coord_pericen$north[2], xmax=coord_pericen$south[2], ymin=-Inf, ymax=Inf, fill="grey90") +
        annotate("rect", xmin=coord_pericen$north[3], xmax=coord_pericen$south[3], ymin=-Inf, ymax=Inf, fill="grey90") +
        annotate("rect", xmin=coord_pericen$north[4], xmax=coord_pericen$south[4], ymin=-Inf, ymax=Inf, fill="grey90") +
        annotate("rect", xmin=coord_pericen$north[5], xmax=coord_pericen$south[5], ymin=-Inf, ymax=Inf, fill="grey90") +
        geom_line(data=dat, aes(x=cumwindow, y=cMMb, colour=genotype), size=0.3) +
        scale_x_continuous(name = "Coordinates (Mb)", labels=scales::label_number(scale = 1/1000000), breaks=c(seq(1, max(mC_col$cumwindow), 20*10^6), 120000000)) +
        scale_y_continuous(limits=landscape_ylim) +
        scale_colour_manual(values = pal) +
        geom_vline(xintercept=tha.cum, colour="Black", size=0.2) +
        geom_vline(xintercept=centromeres.cum, colour="Black", linetype="dashed", size=0.2) +
        labs(y="cM/Mb") +
        theme_classic() +
        theme(legend.key.size=unit(0.2, "inches"),
        legend.title=element_text(size=7),
        legend.text=element_text(size=7)) +
        theme(legend.title=element_blank()) +
        theme(legend.position = c(1, 1)) +
        theme(legend.justification = c(1,0.8)) +
        theme(legend.box.margin = margin(0.5,0.5,0.5,0.5)) +
        theme(text=element_text(size=9, colour="black"),
        axis.text=element_text(colour="black", size=7))

        pdf(file=paste0(prefix, "_co_landscape_ma-co_", ma_width.co,".pdf"), width=plot_width, height=plot_height)
        print(p)
        dev.off()
        png(file=paste0(prefix, "_co_landscape_ma-co_", ma_width.co,".png"), width=plot_width, height=plot_height, unit="in", res=300)
        print(p)
        dev.off()
}

## all merged
drawCoLandscape(co_ma, pal_GBS, paste0("results/07_default_analysis/", libname), 5.7, 2.0)
## without cmt3
# drawCoLandscape(filter(co_ma, genotype != "cmt3"), pal_GBS, "plots/h2aw_dmigs", 5.7, 2.0)

drawCoDiffLandscape <- function(dat, gt_test, gt_ctrl, linecolour, y_label, prefix, plot_width, plot_height, y_min, y_max){
    # dat <- co_ma
    # linecolour <- pal_GBS[c("cmt3", "h2aw")]
    # dat_cMMb <- mutate(dat, cMMb = smooth * 10 * 100 )

    # gt_test <- c("h2aw", "dmigs")
    # gt_ctrl <- "wt"

    diff_landscape_ylim <- c(y_min, y_max)

    dat_ctrl <- filter(dat, genotype == gt_ctrl)
    dat_diff <- list()

    for(i in 1:length(gt_test)){
        dat_test <- filter(dat, genotype == gt_test[i])
        dat_diff[[i]] <- dplyr::select(dat_ctrl, !c("genotype", "smooth", "cMMb"))
        dat_diff[[i]] <- add_column(dat_diff[[i]], difference =  dat_test$cMMb - dat_ctrl$cMMb)
    }
    names(dat_diff) <- gt_test
    dat_diff_bind <- bind_rows(dat_diff, .id="genotype")

    y_label <- bquote(Delta ~ " cM/Mb (" ~ .(y_label) ~ ")")

    p <- ggplot() +
        annotate("rect", xmin=coord_pericen$north[1], xmax=coord_pericen$south[1], ymin=-Inf, ymax=Inf, fill="grey90") +
        annotate("rect", xmin=coord_pericen$north[2], xmax=coord_pericen$south[2], ymin=-Inf, ymax=Inf, fill="grey90") +
        annotate("rect", xmin=coord_pericen$north[3], xmax=coord_pericen$south[3], ymin=-Inf, ymax=Inf, fill="grey90") +
        annotate("rect", xmin=coord_pericen$north[4], xmax=coord_pericen$south[4], ymin=-Inf, ymax=Inf, fill="grey90") +
        annotate("rect", xmin=coord_pericen$north[5], xmax=coord_pericen$south[5], ymin=-Inf, ymax=Inf, fill="grey90") +
        geom_line(data=dat_diff_bind, aes(x=cumwindow, y=difference, colour=genotype), size=0.3) +
        scale_x_continuous(name = "Coordinates (Mb)", labels=scales::label_number(scale = 1/1000000), breaks=c(seq(1, max(mC_col$cumwindow), 20*10^6), 120000000)) +
        scale_y_continuous(limits=diff_landscape_ylim) +
        scale_colour_manual(values=linecolour) +
        geom_vline(xintercept=tha.cum, colour="Black", size=0.2) +
        geom_vline(xintercept=centromeres.cum, colour="Black", linetype="dashed", size=0.2) +
        geom_hline(yintercept=0, size=0.2, colour="Black") +
        # labs(y=expression(delta)) +
        labs(y=y_label) +
        theme_classic() +
        theme(legend.key.size=unit(0.2, "inches"),
        legend.title=element_text(size=7),
        legend.text=element_text(size=7)) +
        theme(legend.title=element_blank()) +
        theme(legend.position = c(1, 1)) +
        theme(legend.justification = c(1,0.8)) +
        theme(legend.box.margin = margin(0.5,0.5,0.5,0.5)) +
        theme(text=element_text(size=9, colour="black"),
        axis.text=element_text(colour="black", size=7))

        pdf(file=paste0(prefix, "_co-diff_landscape_ma-co_", ma_width.co,".pdf"), width=plot_width, height=plot_height)
        print(p)
        dev.off()
        png(file=paste0(prefix, "_co-diff_landscape_ma-co_", ma_width.co,".png"), width=plot_width, height=plot_height, unit="in", res=300)
        print(p)
        dev.off()
}

drawCoDiffLandscape(co_ma, libname, "wt", pal_GBS, "-WT", paste0("results/07_default_analysis/", libname, "-wt"), 5.7, 2, -10.5, 10)
# chromatin landscape

#' 7 TEL-CEN plotting
# distFromTel() : calculate the distance from telomere in two directions (north or south).
# It takes *.cotable as input

#-----Note-----
# This code is kind of trikcy at a first glance.
# The first idea of TEL-CEN scaling that pop in my mind is below:
#   1. scale the CO coordinate into proportion
#   2. tile the proportion into bin
#   However, this does not work because when you scale the coordinate before tiling,
# you will eventually mixup the CO coordinates completely. It hides the pattern.
#   Therefore, before converting into proportion, just tile first so that the CO sites fromthe same chromosome can be tied together a little. Maybe it's the point where the method needs to be modified. 
#---------------


# distFromTel() : convert coordinate into proportion of distance from cent to tel
distFromTel <- function(dat, binsize, spline_df, target_column){
#		dat <- cos.all.list.bin100$hcr2
#		binsize <- 0.01
#		mafilt.size <- 9
        # dat <- filter(co_ma, genotype == "wt")
        # dat <- mC_col
        # target_column <- "mC.ratio"
        # binsize <- 0.01
        # mafilt.size <- 7
        # n_knots <- 10

        dat <- arrange(dat, chr, window)
        dat$cent <- rep(centromeres, times=table(dat$chr))
        dat$end <- rep(chr.ends, times=table(dat$chr))
        
        left <- filter(dat, window < cent) %>%
                add_column(arm="north") %>%
                mutate(coord.prop=window/cent) %>%
				filter(chr %in% c("Chr1", "Chr3", "Chr5"))

		# Excluded north arm of Chr2 and Chr4 as these regions is too short for fair comparison with other chromosomes
        right <- filter(dat, window >= cent) %>%
                add_column(arm="south") %>%
                mutate(coord.prop=(end-window)/(end-cent))
        
        prop <- bind_rows(left, right) %>%
            arrange(coord.prop)

        # bin.spline <- smooth.spline(prop$coord.prop, prop$mean.coInWindow, nknots=10)

        # plot(prop$coord.prop, prop$mean.coInWindow)
        # lines(prop$coord.prop, c(bin.spline$y, 0), lty=2, col=2, lwd=2)
        
        collect_bin <- NULL
        wins <- seq(0, 1, by=binsize)

        for(j in 1:(length(wins)-1)){
                bin <- mean(prop[which(prop$coord.prop >= wins[j] & prop$coord.prop < wins[j+1]),][[target_column]])
                collect_bin <- c(collect_bin, bin)
        }

        # collect_bin.smooth <- smooth.spline(wins[2:length(wins)], collect_bin, nknots=n_knots)$y
        collect_bin.smooth <- smooth.spline(wins[2:length(wins)], collect_bin, df=spline_df)$y

        # k=mafilt.size
        # filt <- 1/(2*k+1)
        # filt <- rep(filt, 2*k+1)
        # filt.dat <- stats::filter(collect_bin, filt)
        # filt.dat <- loess(formula = paste(target_column, "coord.prop", sep="~"), span=0.1, degree=1, data=prop)
        # filt.dat <- filt.dat$fitted

        res <- tibble(prop=wins[2: length(wins)],
                        bin=collect_bin,
                        bin.smooth=collect_bin.smooth
                    #     bin=collect_bin,
                    #   bin.ma=filt.dat,
                      )
        
        return(res)
}

spline_df <- 11
co_telcen.list <- lapply(co_tsv.list, distFromTel, 0.01, spline_df, "mean.coInWindow")
co_telcen.bind <- bind_rows(co_telcen.list, .id="sample") %>%
    mutate(sample = factor(sample, levels=c("wt", libname)))
mC_col_telcen <- distFromTel(mC_col, 0.01, spline_df, "mC.ratio")


## define centromeric-pericentromeric region
## region where DNA methylation level is higher than the average
define_telcen_pericen <- mC_col_telcen %>%
    mutate(distFromAvg = abs(bin.smooth - mean(bin.smooth))) %>%
    arrange(distFromAvg)
telcen_pericen_border <- define_telcen_pericen$prop[1]

drawTelCenProfile <- function(dat_co, dat_chromatin, chrom_label, prefix, plot_width, plot_height){
    # dat_co <- filter(co_telcen.bind, sample != "cmt3")
    # dat_chromatin <- mC_col_telcen
    # dat_chromatin <- k9me2_telcen
    # chrom_label <- "mC ratio"
    # chrom_label <- "k9me2"

    # If smoothed mean CO is negative value, offset by adding the minimum smoothed mean CO
    if(min(dat_co$bin.smooth) < 0){
        dat_co$bin.smooth <- dat_co$bin.smooth - min(dat_co$bin.smooth)
    }

    pal <- pal_GBS

    # scale two different y-axis
    range_co <- summary(dat_co$bin.smooth)["Max."] - summary(dat_co$bin.smooth)["Min."]
    range_chromatin <- summary(dat_chromatin$bin.smooth)["Max."] - summary(dat_chromatin$bin.smooth)["Min."]
    scale_y <- range_co/range_chromatin
    shift_y <- (min(dat_co$bin.smooth, na.rm=T) - min(dat_chromatin$bin.smooth, na.rm=T)*scale_y)*0.9
    # min(dat_co$bin.smooth, na.rm=T) - min(dat_chromatin$bin.smooth, na.rm=T)*scale_y

    mean_cos <- dat_co %>%
            group_by(sample) %>%
            summarise(meanCO=mean(bin, na.rm=TRUE))
    
    yscale_lim <- c(min(dat_co$bin.smooth, na.rm=T)*0.9, max(dat_co$bin.smooth, na.rm=T)*1.1)



    p <- ggplot() +
        # annotate("rect", xmin=telcen_pericen_border, xmax=Inf, ymin=-Inf, ymax=Inf, fill="grey90") +
        geom_vline(xintercept=telcen_pericen_border, size=0.2, linetype = "dashed", colour="grey60") +
        geom_line(data=dat_co, aes(x=prop, y=bin.smooth, colour=sample), size=0.4) +
        # geom_smooth(data=dat_co, aes(x=prop, y=bin, colour=sample), size=0.4, method="lm", formula = y ~ splines::bs(x, 3), se=FALSE) +
        geom_hline(data=mean_cos, aes(yintercept=meanCO, colour=sample), linetype="dashed", size=0.3) +
        geom_area(data=dat_chromatin, aes(x=prop, y=(bin*scale_y) + shift_y), fill="grey80", alpha=0.5) +
        scale_y_continuous(name = "cM/Mb", sec.axis=sec_axis(~(. - shift_y)/scale_y, name=chrom_label)) +
        coord_cartesian(ylim=yscale_lim) +
        scale_x_continuous(name="Distance from the telomere (ratio)",
                           breaks=seq(0, 1, by=0.25),
                           label=c("TEL", 0.25, 0.5, 0.75, "CEN"))+
        scale_colour_manual(values=pal) +
        theme_classic() +
        theme(legend.key.size=unit(0.1, "inches"),
              legend.title=element_text(size=7),
              legend.text=element_text(size=7),
              legend.position="top") +
        theme(legend.title = element_blank()) +
        theme(text=element_text(size=9, colour="black"),
        axis.text=element_text(size=7, colour="black"))
    
    pdf(file=paste0(prefix, "_spline-df_", spline_df, ".pdf"), plot_width, plot_height)
    print(p)
    dev.off()

    png(file=paste0(prefix, "_spline-df_", spline_df, ".png"), plot_width, plot_height, unit="in", res=300)
    print(p)
    dev.off()
}
drawTelCenProfile(filter(co_telcen.bind), mC_col_telcen, "DNA methylation (ratio)", paste0("results/07_default_analysis/", libname, "-wt_telcen"), 2.75, 2.5)

drawTelCenDiffProfile <- function(dat_co, ctrl_group, prefix, plot_width, plot_height){
    # dat_co <- filter(co_telcen.bind, sample != "cmt3")
    # ctrl_group <- "wt"
    # dat_chromatin <- mC_col_telcen
    # dat_chromatin <- k9me2_telcen
    # chrom_label <- "mC ratio"
    # chrom_label <- "k9me2"

    dat_co.ctrl <- filter(dat_co, sample == ctrl_group)
    dat_co.diff <- filter(dat_co, sample != ctrl_group) %>%
        add_column(bin.smooth.wt = rep(dat_co.ctrl$bin.smooth, times=length(unique(.$sample)))) %>%
        mutate(bin.smooth.diff = bin.smooth - bin.smooth.wt)

    pal <- pal_GBS

    yscale_lim <- c(min(dat_co.diff$bin.smooth.diff, na.rm=T)*0.9, max(dat_co.diff$bin.smooth.diff, na.rm=T)*1.1)

    y_label <- paste0("-", ctrl_group)
    y_label <- bquote(Delta ~ " cM/Mb (" ~ .(y_label) ~ ")")



    p <- ggplot() +
        # annotate("rect", xmin=telcen_pericen_border, xmax=Inf, ymin=-Inf, ymax=Inf, fill="grey90") +
        geom_vline(xintercept=telcen_pericen_border, size=0.2, linetype = "dashed", colour="grey60") +
        geom_line(data=dat_co.diff, aes(x=prop, y=bin.smooth.diff, colour=sample), size=0.4) +
        geom_hline(yintercept=0, size=0.2, colour="Black") +
        # geom_smooth(data=dat_co, aes(x=prop, y=bin, colour=sample), size=0.4, method="lm", formula = y ~ splines::bs(x, 3), se=FALSE) +
        # geom_hline(data=mean_cos, aes(yintercept=meanCO, colour=sample), linetype="dashed", size=0.3) +
        # geom_area(data=dat_chromatin, aes(x=prop, y=(bin*scale_y) + shift_y), fill="grey90", alpha=0.5) +
        scale_y_continuous(name = y_label) +
        # coord_cartesian(ylim=yscale_lim) +
        scale_x_continuous(name="Distance from the telomere (ratio)",
                           breaks=seq(0, 1, by=0.25),
                           label=c("TEL", 0.25, 0.5, 0.75, "CEN"))+
        scale_colour_manual(values=pal) +
        theme_classic() +
        theme(legend.key.size=unit(0.1, "inches"),
              legend.title=element_text(size=7),
              legend.text=element_text(size=7),
              legend.position="top") +
        theme(legend.title = element_blank()) +
        theme(text=element_text(size=9, colour="black"),
        axis.text=element_text(size=7, colour="black"))
    
    pdf(file=paste0(prefix, "_spline-df_", spline_df, ".pdf"), plot_width, plot_height)
    print(p)
    dev.off()

    png(file=paste0(prefix, "_spline-df_", spline_df, ".png"), plot_width, plot_height, unit="in", res=300)
    print(p)
    dev.off()
}

drawTelCenDiffProfile(co_telcen.bind, "wt", paste0("results/07_default_analysis/", libname, "_telcen-diff-to-wt"), 2.3, 2.5)
