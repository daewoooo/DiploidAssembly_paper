#' Function to plot haplotype assignment of phased assemblies to the reference 
#' haplotypes (e.g. trio family based).
#' 
#' Input .tsv files are made by first splitting phased assemblies (in FASTA format) into a
#' 1Mb long chunks, aligning them to the reference genome of interest (GRCh38) and then 
#' haplotagging these chunks using whatshap haplotag function.
#' 
#'
#' @param H1.tsv Table of assignments of H1 assembly to the reference haplotypes.
#' @param H2.tsv Table of assignments of H2 assembly to the reference haplotypes.
#' @param genome A \pkg{biovizBase} reference genome id. default: 'hg38'
#' @param chromosomes User defined set of chromosomes to plot.
#' @param title A \code{character} to use as a title of the plot.
#' @inheritParams counts2ranges
#' @return A \code{ggplot} object.
#' @author David Porubsky
#' @export
#'

plotHapAssignments <- function(H1.tsv, H2.tsv, genome = 'hg38', chromosomes=NULL, title=NULL) {
  ## Load the data
  H1.tsv <- read.table(H1.tsv, col.names=c('contig','chromosome','position','flags','length','haplotype'))
  H2.tsv <- read.table(H2.tsv, col.names=c('contig','chromosome','position','flags','length','haplotype'))
  ## Add end positions
  H1.tsv$end <- with(H1.tsv, position + length - 1)
  H2.tsv$end <- with(H2.tsv, position + length - 1)
  ## Convert to Genomic ranges
  H1.gr1 <- makeGRangesFromDataFrame(H1.tsv[H1.tsv$haplotype=="HP:i:1",], seqnames.field='chromosome', start.field='position', end.field='end')
  H1.gr2 <- makeGRangesFromDataFrame(H1.tsv[H1.tsv$haplotype=="HP:i:2",], seqnames.field='chromosome', start.field='position', end.field='end')
  H2.gr1 <- makeGRangesFromDataFrame(H2.tsv[H2.tsv$haplotype=="HP:i:1",], seqnames.field='chromosome', start.field='position', end.field='end')
  H2.gr2 <- makeGRangesFromDataFrame(H2.tsv[H2.tsv$haplotype=="HP:i:2",], seqnames.field='chromosome', start.field='position', end.field='end')
  H1.gr1$hap <- 'H1'
  H1.gr2$hap <- 'H2'
  H2.gr1$hap <- 'H1'
  H2.gr2$hap <- 'H2'
  H1.gr <- c(H1.gr1, H1.gr2)
  H2.gr <- c(H2.gr1, H2.gr2)
  
  ## Get coressponding ideogram from the database
  suppressMessages( hg38IdeogramCyto <- biovizBase::getIdeogram(genome, cytobands = TRUE) )
   
  ## Remove chromosomes that are not defined in chromosomes
  if (!is.null(chromosomes)) {
    hg38IdeogramCyto <- keepSeqlevels(hg38IdeogramCyto, value = chromosomes, pruning.mode = 'coarse')
    seqlevels(hg38IdeogramCyto) <- chromosomes

    H1.gr <- keepSeqlevels(H1.gr, value = chromosomes, pruning.mode = 'coarse')
    H2.gr <- keepSeqlevels(H2.gr, value = chromosomes, pruning.mode = 'coarse')
  }  
   
  ## Set the ggplot theme for final chromosome ideogram
  theme_vertical <- theme(legend.position ="top",
                           axis.line = element_blank(),
                           axis.text.x=element_blank(), 
                           axis.ticks.x=element_blank(),   
                           strip.text.y = element_text(angle = 180),
                           strip.text.x = element_text(size = 7, face="bold"),
                           panel.grid.major = element_blank(), 
                           panel.grid.minor = element_blank(),
                           panel.background = element_blank(),
                           panel.spacing.x=unit(1.5, "lines"))
  
  
  ## Convert ideogram bands stored in GRanges object into the data.frame
  ideo.df <- as.data.frame(hg38IdeogramCyto)
  
  ## Get chromosome breaks and labels
  max.len <- signif(max(ideo.df $end), digits = 2)
  breaks <- seq(from = 0, to = max.len, length.out = 6)
  labels <- breaks / 1000000
  labels <- paste0(labels, 'Mb')
  chr.num <- length(unique(ideo.df $seqnames))
  
  ## Plot the ideogram using ggplot2
  ideo <- ggplot() + 
    geom_rect(data=ideo.df, aes(ymin=start, ymax=end, xmin=0, xmax=1, fill=gieStain), color='black', show.legend=FALSE) + 
    scale_fill_giemsa() +
    facet_grid(. ~ seqnames, switch = 'x') + 
    scale_y_continuous(breaks = breaks, labels = labels) +
    theme_vertical +
    xlab("") +
    ylab("")
  
  ## Prepare data for plotting
  H1.gr$xmin <- 1
  H1.gr$xmax <- 2
  H2.gr$xmin <- -1
  H2.gr$xmax <- 0
  plt.gr <- c(H1.gr, H2.gr)
  names(plt.gr) <- NULL
  
  ## Convert user supplied GRanges object into the data.frame
  plt.df <- as.data.frame(plt.gr)
  plt <- ideo + 
    new_scale("fill") + 
    geom_rect(data=plt.df, aes(ymin=start, ymax=end, xmin=xmin, xmax=xmax, fill=hap)) +
    scale_fill_manual(values = c('deepskyblue', 'gold'), name="")
  
  ## Decrese spacing between facets  
  plt <- plt + theme(panel.spacing.x=unit(0.5, "lines") , panel.spacing.y=unit(1,"lines"))
  
  ## Add title
  if (!is.null(title)) {
    plt <- plt + ggtitle(title)
  }
  return(plt)
}
