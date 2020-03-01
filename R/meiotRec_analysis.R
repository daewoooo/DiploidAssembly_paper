## Load required libraries
library(SaaRclust)
library(data.table)
library(fastseg)
library(BSgenome.Hsapiens.UCSC.hg38)

chromosomes <- paste0("chr", c(1:22, "X"))
#inputdirectory <- "/home/porubsky/WORK/Diploid_assembly/eval_meiotRecs/"
inputdirectory <- "/media/porubsky/Elements/Diploid_assembly/eval_meiotRecs/"

## Load required functions
source("/home/porubsky/DiploidAssembly_paper/R/meiotRec_functions.R")

## Load SNV data
HG00733.snvs <- loadSNVtable(snv.tab = file.path(inputdirectory, 'HG00733/snv_snv_h12.bed.gz'), rmdup = TRUE, chromosomes = chromosomes, family = 'child', sampleName = 'HG00733')
HG00732.snvs <- loadSNVtable(snv.tab = file.path(inputdirectory, 'HG00732/snv_snv_h12.bed.gz'), rmdup = TRUE, chromosomes = chromosomes, family = 'parent', sampleName = 'HG00732')
HG00731.snvs <- loadSNVtable(snv.tab = file.path(inputdirectory, 'HG00731/snv_snv_h12.bed.gz'), rmdup = TRUE, chromosomes = chromosomes, family = 'parent', sampleName = 'HG00731')

## Find meiotic recombination breakpoints
hap.comparisons <- getRecombMap(child = HG00733.snvs, par1 = HG00732.snvs, par2 = HG00731.snvs)

## Plot meiotic recombination breakpoints
plt.df <- as.data.frame(hap.comparisons$parental.snvs)
plt1 <- ggplot() +
  geom_point(data=plt.df[plt.df$H1.comp != 'N',], aes(x=start, y=0.5, color=H1.comp), shape=108) +
  geom_point(data=plt.df[plt.df$H2.comp != 'N',], aes(x=start, y=1, color=H2.comp), shape=108) +
  scale_color_manual(values = c('cadetblue4','darkgoldenrod4','cadetblue3','darkgoldenrod2'), name='') +
  facet_grid(seqnames ~ ., switch = 'y') +
  theme_void() +
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(strip.text.y = element_text(angle = 180)) +
  theme(legend.position = 'top')

## Plot predicted haplotype segments
plt.df <- as.data.frame(hap.comparisons$parental.segments)

theme_vertical <- theme(legend.position ="top",
                   axis.line = element_blank(),
                   axis.text.x=element_blank(), 
                   axis.ticks.x=element_blank(),   
                   strip.text.y = element_text(angle = 180),
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank())

plt.df$ymin <- 0
plt.df$ymax <- 0
plt.df$ymin[plt.df$inherited == 'HG00731'] <- 0
plt.df$ymax[plt.df$inherited == 'HG00731'] <- 2
plt.df$ymin[plt.df$inherited == 'HG00732'] <- 3
plt.df$ymax[plt.df$inherited == 'HG00732'] <- 5

## Load GAPS and CENTROMERS
gaps <- read.table("/media/porubsky/Elements/Diploid_assembly/GRCh38_annot/gaps_GRCh38", stringsAsFactors = FALSE)
gaps.gr <- makeGRangesFromDataFrame(gaps, seqnames.field = 'V2', start.field = 'V3', end.field = 'V4')
gaps.gr <- keepSeqlevels(gaps.gr, value = chromosomes, pruning.mode = 'coarse')
gaps.df <- as.data.frame(gaps.gr)

cent <- read.table("/media/porubsky/Elements/Diploid_assembly/GRCh38_annot/centromeres_GRCh38.bed", stringsAsFactors = FALSE)
cent.gr <- makeGRangesFromDataFrame(cent, seqnames.field = 'V2', start.field = 'V3', end.field = 'V4')
cent.gr <- keepSeqlevels(cent.gr, value = chromosomes, pruning.mode = 'coarse')
cent.df <- as.data.frame(cent.gr)

#vertical
plt2 <- ggplot(plt.df) + 
  geom_rect(aes(xmin=ymin, xmax=ymax, ymin=start, ymax=end, fill=ID)) + 
  geom_rect(data=gaps.df, aes(xmin=0, xmax=2, ymin=start, ymax=end), fill='white') + 
  geom_rect(data=gaps.df, aes(xmin=3, xmax=5, ymin=start, ymax=end), fill='white') + 
  geom_rect(data=cent.df, aes(xmin=0, xmax=2, ymin=start, ymax=end), fill='red') + 
  geom_rect(data=cent.df, aes(xmin=3, xmax=5, ymin=start, ymax=end), fill='red') + 
  facet_grid(. ~ seqnames, switch = 'x') + 
  scale_fill_manual(values=c('cadetblue4','darkgoldenrod4','cadetblue3','darkgoldenrod2')) + 
  xlab("") +
  ylab("") +
  theme_vertical

## Compare predicted meiotic breakpoints with HGSVC data
hgsvc.pur <- read.table(file.path(inputdirectory, "PUR_refined_breakpoints_HGSVC.txt"), header = TRUE, stringsAsFactors = FALSE)
hgsvc.pur.gr <- makeGRangesFromDataFrame(hgsvc.pur, keep.extra.columns = TRUE)
hgsvc.pur.gr <- hgsvc.pur.gr[hgsvc.pur.gr$parent != 'WH']
hgsvc.pur.gr$parent[hgsvc.pur.gr$parent == 'F'] <- 'HG00731'
hgsvc.pur.gr$parent[hgsvc.pur.gr$parent == 'M'] <- 'HG00732'

hgsvc.pur.gr$level <- 0
hgsvc.pur.gr$level[hgsvc.pur.gr$parent == 'HG00732'] <- 4
hgsvc.pur.gr$level[hgsvc.pur.gr$parent == 'HG00731'] <- 1
hgsvc.pur.df <- as.data.frame(hgsvc.pur.gr)

plt2 <- plt2 + 
  geom_point(data=hgsvc.pur.df, aes(y=start, x=level, color=parent)) +
  facet_grid(. ~ seqnames, switch = 'x') +
  scale_color_manual(values = c('blue3', 'brown1'), name="HGSVC_recomb")

## Export final plot
final.plt <- plot_grid(plt1, plt2, ncol = 1, labels = c('A','B'))
destination <- file.path(inputdirectory, "HG00733_meiot_recs.pdf")
ggsave(filename = destination, plot = final.plt, width = 12, height = 12)

#recomb.break.gr <- unlist(recomb.break.grl, use.names = FALSE)
#findOverlaps(hgsvc.pur.gr, recomb.break.gr, maxgap = 100000)
#subsetByOverlaps(hgsvc.pur.gr, recomb.break.gr)

## Plot distance of male and female specific breakpoint from telomeres ##
chr.len <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)[chromosomes]
hg38.starts <- GRanges(seqnames=chromosomes, ranges=IRanges(start=1, end=1))
hg38.ends <- GRanges(seqnames=chromosomes, ranges=IRanges(start=chr.len, end=chr.len))
hg38.tels <- sort(c(hg38.starts, hg38.ends))

recomb.break.gr <- hap.comparisons$recombination.breakpoints

male.recomb <- recomb.break.gr[recomb.break.gr$inherited == 'HG00731']
female.recomb <- recomb.break.gr[recomb.break.gr$inherited == 'HG00732']

male.dists <- distanceToNearest(hg38.tels, male.recomb)
female.dists <- distanceToNearest(hg38.tels, female.recomb)

plt.df <- data.frame(dists = rbind(mcols(male.dists), mcols(female.dists)), ID = rep(c('male', 'female'), c(length(male.dists), length(female.dists))))

plt.df <- plt.df[order(plt.df$distance),]
plt.df$x <- 1:nrow(plt.df)
plt3 <- ggplot(plt.df) + 
  geom_point(aes(x=x, y = distance, color=ID), size=2) +
  coord_trans(y = 'log10') +
  geom_hline(yintercept = 5000000) +
  scale_y_continuous(breaks = c(2000000, 5000000, 10000000, 50000000, 100000000), labels = c('2000000', '5000000', '10000000','50000000','100000000')) +
  ylab("Distance to 1st break from TEL (bp)") +
  xlab("Sorted distances from TEL") + theme_bw()
  #geom_boxplot(aes(x=ID, y = distance, fill=ID))

final.plt <- plot_grid(plt2, plt3, ncol = 1, rel_heights = c(2,1.5))
destination <- file.path(inputdirectory, "HG00733_meiot_recs2.pdf")
ggsave(filename = destination, plot = final.plt, width = 12, height = 8)
