## Load required libraries
library(GenomicRanges)
library(RColorBrewer)
library(primatR)

inputfolder <- "/media/porubsky/Elements/HG002/CommonBreaksAsm/"
chromosomes <- paste0('chr', c(1:22, 'X'))

## Load segDup Track
SDs.df <- read.table(file.path(inputfolder, "segDup_GRCh38.gz"))
SDs.df <- SDs.df[,c(2,3,4,27)]
colnames(SDs.df) <- c('seqnames', 'start', 'end', 'fracMatch')
SDs.df <- SDs.df[SDs.df$seqnames %in% chromosomes,]
SDs.df$seqnames <- factor(SDs.df$seqnames, levels = chromosomes)
SDs.gr <- makeGRangesFromDataFrame(SDs.df, keep.extra.columns = TRUE)
SDs.gr <- reduce(SDs.gr)
#SDs.gr <- SDs.gr[width(SDs.gr) >= 10000]

## Load StrandS based diploid assemblies
hg002_h1 <- read.table("/media/porubsky/Elements/HG002/diploidAsm2GRCh38/HG002_H1_p_ctg_cns.bed")
colnames(hg002_h1) <- c('seqnames','start','end','ctg.name','mapq','strand')
hg002_h1.gr <- makeGRangesFromDataFrame(hg002_h1, keep.extra.columns = TRUE)
hg002_h2 <- read.table("/media/porubsky/Elements/HG002/diploidAsm2GRCh38/HG002_H2_p_ctg_cns.bed")
colnames(hg002_h2) <- c('seqnames','start','end','ctg.name','mapq','strand')
hg002_h2.gr <- makeGRangesFromDataFrame(hg002_h2, keep.extra.columns = TRUE)

## Export contigs as UCSC browser
ranges2UCSC(gr = hg002_h1.gr[seqnames(hg002_h1.gr) %in% chromosomes], outputDirectory = "/media/porubsky/Elements/HG002/eval_uab/", index = "HG002_H1_StrandS", id.field = 'ctg.name')
ranges2UCSC(gr = hg002_h2.gr[seqnames(hg002_h2.gr) %in% chromosomes], outputDirectory = "/media/porubsky/Elements/HG002/eval_uab/", index = "HG002_H2_StrandS", id.field = 'ctg.name')

## Load HG002 inversions
hg002_inv <- read.table("/home/porubsky/WORK/HG002/alignments2GRCh38/HG002/Selected_inversions/BreakpointR_results/HG002_INVcalls.bed", skip = 1)
hg002_inv <- hg002_inv[,c(1,2,3,4,6)]
colnames(hg002_inv) <- c('seqnames','start','end','SV.class','gen')
hg002_inv.gr <- makeGRangesFromDataFrame(hg002_inv, keep.extra.columns = TRUE)
## Keep only simple inversions and inverted duplications
hg002_inv.gr <- hg002_inv.gr[hg002_inv.gr$SV.class %in% c('INV', 'invDup')]

## Select contigs that overlaps with inversion calls
#hg002_h1_selected.gr <- subsetByOverlaps(hg002_h1.gr, hg002_inv.gr)
#hg002_h2_selected.gr <- subsetByOverlaps(hg002_h2.gr, hg002_inv.gr)

## Get inversions that overlaps with predicted UAB
hg002.uab <- read.table("/media/porubsky/Elements/HG002/CommonBreaksAsm/sharedBreaks.HG002.txt", header = TRUE, stringsAsFactors = FALSE)
hg002.uab.gr <- makeGRangesFromDataFrame(hg002.uab, keep.extra.columns = TRUE)
hg002_inv.gr <- subsetByOverlaps(hg002_inv.gr, hg002.uab.gr)

## Split contigs by inversions overlap
h1.hits <- findOverlaps(hg002_h1.gr, hg002_inv.gr, maxgap = 1000000)
hg002_h1_selected.grl <- split(hg002_h1.gr[queryHits(h1.hits)], subjectHits(h1.hits))
h2.hits <- findOverlaps(hg002_h2.gr, hg002_inv.gr, maxgap = 1000000)
hg002_h2_selected.grl <- split(hg002_h2.gr[queryHits(h2.hits)], subjectHits(h2.hits))

## Get Rcolorbrewer color schemes
qual.col.pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
color.vector <- col.vector <- unlist(mapply(RColorBrewer::brewer.pal, qual.col.pals$maxcolors, rownames(qual.col.pals)))

## Plot overlap of inversions with contigs
## H1 ##
all.plts <- list()
for (i in seq_along(hg002_h1_selected.grl)) {
  ctg.gr <- hg002_h1_selected.grl[[i]]
  ctg.gr$level <- disjointBins(ctg.gr)
  inv.gr <- hg002_inv.gr[i]
  inv.id <- as.character(inv.gr)
  
  sd.sub.gr <- subsetByOverlaps(SDs.gr, range(ctg.gr, ignore.strand=TRUE))
  sd.df <- as.data.frame(sd.sub.gr)
  
  ctg.df <- as.data.frame(ctg.gr)  
  inv.df <- as.data.frame(inv.gr)  
  
  plt <- ggplot(ctg.df) + 
    geom_rect(data=ctg.df, aes(xmin=start, xmax=end, ymin=level, ymax=level+1, fill='contig', color=ctg.name), size=1) +
    geom_rect(data=inv.df, aes(xmin=start, xmax=end, ymin=0, ymax=1, fill='inv_region'), color='red') +
    geom_rect(data=sd.df, aes(xmin=start, xmax=end, ymin=-1, ymax=0, fill='SDs')) +
    scale_fill_manual(values = c('black','red','orange')) +
    scale_color_manual(values = color.vector, guide='none') +
    theme_bw() +
    ggtitle(inv.id)
  all.plts[[length(all.plts) + 1]] <- plt
}
destination <- "/media/porubsky/Elements/HG002/eval_uab/HG002_H1_asmBreak_byINV_new.pdf"
grDevices::pdf(destination, width=10, height=3)
bquiet = lapply(all.plts, print)
d <- grDevices::dev.off()

## H2 ##
all.plts <- list()
for (i in seq_along(hg002_h2_selected.grl)) {
  ctg.gr <- hg002_h2_selected.grl[[i]]
  ctg.gr$level <- disjointBins(ctg.gr)
  inv.gr <- hg002_inv.gr[i]
  inv.id <- as.character(inv.gr)
  
  sd.sub.gr <- subsetByOverlaps(SDs.gr, range(ctg.gr, ignore.strand=TRUE))
  sd.df <- as.data.frame(sd.sub.gr)
  
  ctg.df <- as.data.frame(ctg.gr)  
  inv.df <- as.data.frame(inv.gr)  
  
  plt <- ggplot(ctg.df) + 
    geom_rect(data=ctg.df, aes(xmin=start, xmax=end, ymin=level, ymax=level+1, fill='contig', color=ctg.name), size=1) +
    geom_rect(data=inv.df, aes(xmin=start, xmax=end, ymin=0, ymax=1, fill='inv_region'), color='red') +
    geom_rect(data=sd.df, aes(xmin=start, xmax=end, ymin=-1, ymax=0, fill='SDs')) +
    scale_fill_manual(values = c('black','red','orange')) +
    scale_color_manual(values = color.vector, guide='none') +
    theme_bw() +
    ggtitle(inv.id)
  all.plts[[length(all.plts) + 1]] <- plt
}
destination <- "/media/porubsky/Elements/HG002/eval_uab/HG002_H2_asmBreak_byINV_new.pdf"
grDevices::pdf(destination, width=10, height=3)
bquiet = lapply(all.plts, print)
d <- grDevices::dev.off()

