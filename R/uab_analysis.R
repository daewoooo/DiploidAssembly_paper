## Plot position of assembly gaps in the human reference genome ##
##################################################################
## Load required libraries
library(GenomicRanges)
library(ggnewscale)
library(biovizBase)
library(ggplot2)
library(ggbio)
library(primatR)
library(cowplot)
library(reshape2)
library(BSgenome.Hsapiens.UCSC.hg38)
library(scales)
library(RColorBrewer)
library(dplyr)
library(gsubfn)

inputfolder <- "/media/porubsky/19AD2DA754D7CFAF/Projects/HG002/CommonBreaksAsm/"
binsize <- 500000
chromosomes <- paste0('chr', c(1:22,'X'))
bsgenome <- BSgenome.Hsapiens.UCSC.hg38
source("/media/porubsky/19AD2DA754D7CFAF/Projects/HG002/CommonBreaksAsm/uab_functions.R")

## Load the StrandS phased assemblies ##
HG002.H1.gr <- bed2ranges(bedfile = "/media/porubsky/19AD2DA754D7CFAF/Projects/HG002/diploidAsm2GRCh38/HG002_H1_p_ctg_cns.bed",
                          index = 'HG002.H1.StrandS',
                          min.align=10000, min.ctg.size=500000, max.gap=100000)
HG002.H1.gr$ctg <- paste0("H1_", HG002.H1.gr$ctg)
HG002.H2.gr <- bed2ranges(bedfile = "/media/porubsky/19AD2DA754D7CFAF/Projects/HG002/diploidAsm2GRCh38/HG002_H2_p_ctg_cns.bed",
                          index = 'HG002.H2.StrandS',
                          min.align=10000, min.ctg.size=500000, max.gap=100000)
HG002.H2.gr$ctg <- paste0("H2_", HG002.H2.gr$ctg)

HG002.StrandS <- suppressWarnings( c(HG002.H1.gr, HG002.H2.gr) )
HG002.StrandS <- keepSeqlevels(HG002.StrandS, value = chromosomes, pruning.mode = 'coarse')
ranges2UCSC(gr = HG002.StrandS, outputDirectory = inputfolder, index = "strandS_ctg_new", colorRGB = '200,0,200', id.field = 'ctg')

## Load the HiC phased assemblies (Shilpa biorXiv) ##
HG002.H1.gr <- bed2ranges(bedfile = "/media/porubsky/19AD2DA754D7CFAF/Projects/HG002/diploidAsm2GRCh38/shilpaAsm/NA24385-denovo-H1.bed",
                          index = 'HG002.H1.HiC',
                          min.align=10000, min.ctg.size=500000, max.gap=100000)
HG002.H2.gr <- bed2ranges(bedfile = "/media/porubsky/19AD2DA754D7CFAF/Projects/HG002/diploidAsm2GRCh38/shilpaAsm/NA24385-denovo-H2.bed",
                          index = 'HG002.H2.HiC',
                          min.align=10000, min.ctg.size=500000, max.gap=100000)

HG002.hic <- suppressWarnings( c(HG002.H1.gr, HG002.H2.gr) )
HG002.hic <- keepSeqlevels(HG002.hic, value = chromosomes, pruning.mode = 'coarse')
ranges2UCSC(gr = HG002.hic, outputDirectory = inputfolder, index = "hic_ctg_new", colorRGB = '0,0,200', id.field = 'ctg')

## Load genome annotations ##
## Add centromere information
cent.df <- read.table(file.path(inputfolder, "centromeres_GRCh38.bed"))
cent.df <- cent.df[,c(2,3,4)]
colnames(cent.df) <- c('seqnames', 'start', 'end')
## Add gap information
gaps.df <- read.table(file.path(inputfolder, "gaps_GRCh38"))
gaps.df <- gaps.df[,c(2,3,4)]
colnames(gaps.df) <- c('seqnames', 'start', 'end')
gaps.df <- gaps.df[gaps.df$seqnames %in% chromosomes,]
## Add SD information
SDs.df <- read.table(file.path(inputfolder, "segDup_GRCh38.gz"))
SDs.df <- SDs.df[,c(2,3,4,27)]
colnames(SDs.df) <- c('seqnames', 'start', 'end', 'fracMatch')
SDs.df <- SDs.df[SDs.df$seqnames %in% chromosomes,]
SDs.df$seqnames <- factor(SDs.df$seqnames, levels = chromosomes)
## Convert annotation to Genomic ranges
SDs.gr <- makeGRangesFromDataFrame(SDs.df, keep.extra.columns = TRUE)
SDs.gr <- reduce(SDs.gr)
SDs.gr <- SDs.gr[width(SDs.gr) >= 10000]
cent.gr <- makeGRangesFromDataFrame(cent.df, keep.extra.columns = TRUE)
gaps.gr <- makeGRangesFromDataFrame(gaps.df, keep.extra.columns = TRUE)

## Load inversion calls for HG002
hg002.inv <- read.table("/media/porubsky/19AD2DA754D7CFAF/Projects/HG002/alignments2GRCh38/HG002/Selected_inversions/BreakpointR_results/HG002_INVcalls.bed", skip = 1, stringsAsFactors = FALSE)
hg002.inv  <- hg002.inv[,c(1,2,3,4,6)]
colnames(hg002.inv) <- c('seqnames','start','end','SVclass','gen')
hg002.inv.gr <- makeGRangesFromDataFrame(hg002.inv, keep.extra.columns = TRUE)
## Keep only simple inversions and inverted duplications
hg002.inv.gr <- hg002.inv.gr[hg002.inv.gr$SVclass %in% c('INV', 'invDup')]
hg002.inv.gr$gen <- recode(hg002.inv.gr$gen, '+' = 'HET', '-' = 'HOM')

## Merge contigs ranges
HG002.ctg.gr <- suppressWarnings( c(HG002.StrandS, HG002.hic) )
## Test collapses
#gr <- GRanges(seqnames='chr9', ranges=IRanges(start=123967027, end=124003122))

breaks.roi.ranges <- getUABs(gr = HG002.ctg.gr, chromosomes = chromosomes, binsize = binsize, stepsize = binsize/2, bsgenome = bsgenome, ID.col = 'ID')

## Export FASTA file from UAB regions
breaks.roi.ranges.fasta <- regions2FASTA(gr = breaks.roi.ranges, bsgenome = BSgenome.Hsapiens.UCSC.hg38)

## Find overlap with genome annotation within 10kb buffer
maxgap <- 10000
breaks.roi.ranges$annot <- ''
hits <- findOverlaps(breaks.roi.ranges, cent.gr, maxgap = maxgap)
breaks.roi.ranges$annot[unique(queryHits(hits))] <- paste0(breaks.roi.ranges$annot[unique(queryHits(hits))], "cent.")
hits <- findOverlaps(breaks.roi.ranges, SDs.gr, maxgap = maxgap)
breaks.roi.ranges$annot[unique(queryHits(hits))] <- paste0(breaks.roi.ranges$annot[unique(queryHits(hits))], "SDs.")
hits <- findOverlaps(breaks.roi.ranges, gaps.gr, maxgap = maxgap)
breaks.roi.ranges$annot[unique(queryHits(hits))] <- paste0(breaks.roi.ranges$annot[unique(queryHits(hits))], "Gap.")
hits <- findOverlaps(breaks.roi.ranges, hg002.inv.gr, maxgap = maxgap)
breaks.roi.ranges$annot[unique(queryHits(hits))] <- paste0(breaks.roi.ranges$annot[unique(queryHits(hits))], "INVs.")
breaks.roi.ranges$annot[breaks.roi.ranges$annot == ''] <- 'unannot'
## Merge categories
# breaks.roi.ranges$annot[grep(breaks.roi.ranges$annot, pattern = 'cent')] <- "Centromere"
# breaks.roi.ranges$annot[grep(breaks.roi.ranges$annot, pattern = 'SDs')] <- "SDs"
# breaks.roi.ranges$annot[grep(breaks.roi.ranges$annot, pattern = 'Gap')] <- "Gap"
# breaks.roi.ranges$annot[grep(breaks.roi.ranges$annot, pattern = 'INVs')] <- "INVs"
# breaks.roi.ranges.df <- as.data.frame(breaks.roi.ranges)

## Export common assembly breaks
destination <- "/media/porubsky/19AD2DA754D7CFAF/Projects/HG002/CommonBreaksAsm/sharedBreaks.HG002.txt"
breaks.roi.ranges.df <- as.data.frame(breaks.roi.ranges)
write.table(x = breaks.roi.ranges.df, file = destination, quote = FALSE, row.names = FALSE)

## Plot genome-wide ideograms of all breaks
colors <- c('cadetblue4', 'darkgoldenrod4','cadetblue2','darkgoldenrod2')
plt1 <- plotUABlollipop(gr = HG002.ctg.gr, breaks.gr = breaks.roi.ranges, ID.col = 'ID', chromosomes = chromosomes, genome.ideo = 'hg38', colors = colors, annot.SD = SDs.df)
plt2 <- plotUABgenome(breaks.gr = breaks.roi.ranges, chromosomes = chromosomes, bsgenome = bsgenome, centromeres.gr = cent.gr, gaps.gr = gaps.gr)
## Export final plots
destination <- "/media/porubsky/19AD2DA754D7CFAF/Projects/HG002/CommonBreaksAsm/sharedBreaks.HG002.lillipops.pdf"
ggsave(filename = destination, plot = plt1, width = 20, height = 7)
destination <- "/media/porubsky/19AD2DA754D7CFAF/Projects/HG002/CommonBreaksAsm/sharedBreaks.HG002.pdf"
ggsave(filename = destination, plot = plt2, width = 10, height = 5)


## Get average perc identity per region
breaks.roi.ranges.df <- as.data.frame(breaks.roi.ranges)
breaks.roi.ranges.df$annot <- factor(breaks.roi.ranges$annot, levels=unique(breaks.roi.ranges$annot))
plt2 <- breaks.roi.ranges.df %>% group_by(annot) %>% summarise(count = n()) %>%
  ggplot(aes(x=annot, y=count, fill=annot)) +
  geom_col() +
  scale_fill_manual(values = brewer.pal(n = 9, name = 'Set1'), name="") +
  geom_text(aes(x=annot, y = count, label=count), vjust=0, hjust=0) +
  theme_bw() +
  ylab("# of assembly breaks overlapping with genomic features") +
  xlab("") +
  coord_flip()
## Export final plot
destination <- "/media/porubsky/19AD2DA754D7CFAF/Projects/HG002/CommonBreaksAsm/sharedBreaks.HG002.annot.pdf"
ggsave(filename = destination, plot = plt2, width = 8, height = 4)

## Get overlap of UAB predited for HG00733
uab.hg00733 <- read.table("/media/porubsky/19AD2DA754D7CFAF/Projects/HG002/CommonBreaksAsm/supplementary_table_6.csv", sep = "," ,header = TRUE, stringsAsFactors = FALSE)
uab.hg00733.gr <- makeGRangesFromDataFrame(uab.hg00733, keep.extra.columns = TRUE)

hits <- findOverlaps(uab.hg00733.gr, breaks.roi.ranges, maxgap = 10000)

breaks.roi.ranges$uab <- "unseen"
breaks.roi.ranges$uab[subjectHits(hits)] <- "seen"

plt.df <- as.data.frame(breaks.roi.ranges)
plt3 <- ideo +
  geom_point(data=plt.df, aes(x = midpoint, y = 0.5, color=uab), size=3, shape=18) +
  facet_grid(seqnames ~ ., switch = 'y') +
  scale_x_continuous(expand = c(0,0)) +
  theme_void() +
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(strip.text.y = element_text(angle = 180))
## Export final plot
destination <- "/media/porubsky/19AD2DA754D7CFAF/Projects/HG002/CommonBreaksAsm/sharedBreaks.HG002.compared.pdf"
ggsave(filename = destination, plot = plt3, width = 10, height = 5)


## TO CLEAN !!!
## Plot additional statistics on contig breakage regions ##
###########################################################
## Expand ranges to the size of the largest range
max.width <- max(width(breaks.roi.ranges))
lookup <- max.width %/% 2
#lookup <- 2500000
lookup.ranges <- GRangesList()
for (i in seq_along(breaks.roi.ranges)) {
  roi <- breaks.roi.ranges[i]
  chr <- unique(seqnames(roi))
  chr.len <- seqlengths(roi)[chr]
  lookup.gr <- GRanges(seqnames=chr, ranges=IRanges(start=max(roi$midpoint - lookup, 1), 
                                                    end=min(roi$midpoint + lookup, chr.len)
  ))
  lookup.ranges[[i]] <- lookup.gr
}
breaks.roi.lookup <- unlist(lookup.ranges, use.names = FALSE)

## Get overlap with SDs
SDs.gr <- makeGRangesFromDataFrame(SDs.df, keep.extra.columns = TRUE)
seqlengths(SDs.gr) <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)[seqlevels(SDs.gr)]
SD.content <- list()
for (i in seq_along(breaks.roi.lookup)) {
  roi <- breaks.roi.lookup[i]
  SDs.roi <- subsetByOverlaps(SDs.gr, roi)
  if (length(SDs.roi) > 0) {
    SDs.roi$length <- width(SDs.roi)
    SDs.roi$categ <- as.numeric(findInterval(SDs.roi$fracMatch, c(0.9, 0.98, 0.99)))
    SDs.roi$categ <- dplyr::recode(SDs.roi$categ, `1` = '0.9-0.98', `2` = '0.98-0.99', `3` = '>0.99')
    SDs.roi$categ <- factor(SDs.roi$categ, levels=c('0.9-0.98','0.98-0.99','>0.99'))
    categ.widths <- split(SDs.roi$length, SDs.roi$categ, drop = FALSE)
    categ.widths <- lapply(categ.widths, sum)
    values <- unlist(categ.widths)
    df <- data.frame(ID=names(categ.widths), values=values, region.ID=i, SD.perc=values/width(roi), row.names = NULL)
    SD.content[[i]] <- df
  } else {
    df <- data.frame(ID=c('0.9-0.98','0.98-0.99','>0.99'), values=0, region.ID=i, SD.perc=0, row.names = NULL)
    SD.content[[i]] <- df
  }
}
SD.content.df <- do.call(rbind, SD.content)

plt1 <- ggplot() + 
  geom_col(data=SD.content.df, aes(y = SD.perc, x = region.ID, fill = ID, group = region.ID)) +
  scale_fill_manual(values = c('0.9-0.98'='gray53', '0.98-0.99'='goldenrod2', '>0.99'='darkorange1')) +
  coord_flip() + 
  theme_bw() +
  theme(legend.position="top") +
  xlab(paste0("Region 1-", length(breaks.roi.lookup))) +
  ylab("SD fraction per SD category")

## Get average perc identity per region
breaks.roi.ranges.df <- as.data.frame(breaks.roi.ranges)
breaks.roi.ranges.df$annot <- factor(breaks.roi.ranges$annot, levels=unique(breaks.roi.ranges$annot))
plt2 <- breaks.roi.ranges.df %>% group_by(annot) %>% summarise(count = n()) %>%
 ggplot(aes(x=annot, y=count, fill=annot)) +
  geom_col() +
  scale_fill_manual(values = brewer.pal(n = 8, name = 'Set1'), name="") +
  geom_text(aes(x=annot, y = count, label=count), vjust=0, hjust=0) +
  theme_bw() +
  ylab("# of assembly breaks overlapping with genomic features") +
  xlab("") +
  coord_flip()

## Get overlap with GC content
GC.cont <- breaks.roi.lookup
GC.cont$gc.perc <- GCcontent(BSgenome.Hsapiens.UCSC.hg38, GC.cont)
GC.cont.df <- as.data.frame(GC.cont)
GC.cont.df$region.ID <- 1:nrow(GC.cont.df)
plt3 <- ggplot() + 
  geom_point(data=GC.cont.df , aes(y = C.G, x = region.ID)) +
  coord_flip() + 
  theme_bw() +
  xlab(paste0("Region 1-", length(breaks.roi.lookup))) +
  ylab('GC percentage')

## Get overlap with Gaps
gaps.df <- read.table(file.path(inputfolder, "gaps_GRCh38"))
gaps.df <- gaps.df[,c(2,3,4)]
colnames(gaps.df) <- c('seqnames', 'start', 'end')
gaps.df <- gaps.df[gaps.df$seqnames %in% chromosomes,]
gaps.df$fracMatch <- 1
gaps.gr <- makeGRangesFromDataFrame(gaps.df, keep.extra.columns = TRUE)
hits <- findOverlaps(breaks.roi.lookup, gaps.gr)
overl.gaps <- length(unique(subjectHits(hits)))
overl.gaps.frac <- overl.gaps / length(breaks.roi.lookup)
plt.df <- data.frame(overlap.gaps = overl.gaps.frac, total.gaps = 1 - overl.gaps.frac)
plt.df <- melt(plt.df)
plt.df$variable <- factor(plt.df$variable, levels = c('total.gaps', 'overlap.gaps'))
plt4 <- ggplot() +
  geom_col(data=plt.df, aes(x=1, y=value, fill=variable)) +
  scale_fill_manual(values = c('gray49', 'gray74')) +
  theme_bw() +
  xlab("") +
  ylab("Overlap with known assembly gaps") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

## Plot un-annotated assembly breaks
## Plot the ideogram using ggplot2
ideo <- ggplot() + 
  geom_rect(ideo.df, aes(ymin=start, ymax=end, xmin=0, xmax=0.9, fill=gieStain), color='black', show.legend=FALSE) + 
  scale_fill_giemsa() +
  facet_grid(. ~ seqnames, switch = 'x') + 
  scale_y_continuous(breaks = breaks, labels = labels) +
  theme_vertical +
  xlab("") +
  ylab("") +
  theme(panel.spacing.x = unit(3, "lines"))

breaks.roi.ranges.noSD.noGap.noCent.df <- as.data.frame(breaks.roi.ranges.noSD.noGap.noCent)
#plt5 <- ideo + geom_point(data=breaks.roi.ranges.noSD.df, aes(x=1, y=midpoint, color=ID), shape=60, size=7) +
#  scale_color_manual(values = c('magenta', 'orange1'), name="")
plt5 <- ideo + geom_point(data=breaks.roi.ranges.noSD.noGap.noCent.df, aes(x=1, y=midpoint), color='magenta', shape=60, size=7) +
  ggtitle("Genomic positions of unannotated assembly breaks")

bottom <- plot_grid(plt1, plt2, plt3, plt4, nrow = 1, align = 'h', axis = 't', rel_widths = c(2,2,2,1))
final.plt <- plot_grid(plt, bottom, plt5, ncol = 1, rel_heights = c(2,1,1.5))  

ggsave(filename = file.path(inputfolder, 'assembly_comparison_v6.pdf'), plot = final.plt, width = 20, height = 16)

## Export ROI regions
destination <- file.path(inputfolder, 'assembly.breaks.csv')
write.table(breaks.roi.ranges.df, file = destination, quote = FALSE, row.names = FALSE, sep = ',')
destination <- file.path(inputfolder, 'unannotated.assembly.breaks.csv')
write.table(breaks.roi.ranges.df[breaks.roi.ranges.df$annot == 'unannot',], file = destination, quote = FALSE, row.names = FALSE, sep = ',')
