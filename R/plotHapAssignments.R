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
#' @param custom.ideo A \code{data.frame} 
#' @param chromosomes User defined set of chromosomes to plot.
#' @param title A \code{character} to use as a title of the plot.
#' @inheritParams counts2ranges
#' @return A \code{ggplot} object.
#' @author David Porubsky
#' @export
#'

plotHapAssignments <- function(H1.tsv, H2.tsv, genome = 'hg38', custom.ideo=NULL, chromosomes=NULL, title=NULL) {
  ## Helper function
  getPhasingAcc <- function(tsv.file, chromosomes=chromosomes) {
    df <- read.table(tsv.file)
    if (!is.null(chromosomes)) {
      ## Filter user defined chromsomes
      df <- df[df$V2 %in% chromosomes,]
    }  
    ## Make sure haplotype assignment is defined as factor
    df$V6 <- factor(df$V6, levels = unique(df$V6))
    ## Split by chromosomes
    dfl <- split(df, df$V2)
    dfl.counts <- lapply(dfl, function(x) table(x$V6))
    dfl.counts <- do.call(rbind, dfl.counts)
    correct <- pmax(dfl.counts[,1], dfl.counts[,2])
    incorrect <- pmin(dfl.counts[,1], dfl.counts[,2])
    result <- data.frame(chr = rownames(dfl.counts), correct = correct, incorrect = incorrect)
    return(result)
  }
  
  ## Load the data
  H1.tsv.df <- read.table(H1.tsv, col.names=c('contig','chromosome','position','flags','length','haplotype'))
  H2.tsv.df <- read.table(H2.tsv, col.names=c('contig','chromosome','position','flags','length','haplotype'))
  ## Add end positions
  H1.tsv.df$end <- with(H1.tsv.df, position + length - 1)
  H2.tsv.df$end <- with(H2.tsv.df, position + length - 1)
  ## Convert to Genomic ranges
  H1.gr1 <- makeGRangesFromDataFrame(H1.tsv.df[H1.tsv.df$haplotype=="HP:i:1",], seqnames.field='chromosome', start.field='position', end.field='end')
  H1.gr2 <- makeGRangesFromDataFrame(H1.tsv.df[H1.tsv.df$haplotype=="HP:i:2",], seqnames.field='chromosome', start.field='position', end.field='end')
  H2.gr1 <- makeGRangesFromDataFrame(H2.tsv.df[H2.tsv.df$haplotype=="HP:i:1",], seqnames.field='chromosome', start.field='position', end.field='end')
  H2.gr2 <- makeGRangesFromDataFrame(H2.tsv.df[H2.tsv.df$haplotype=="HP:i:2",], seqnames.field='chromosome', start.field='position', end.field='end')
  H1.gr1$hap <- 'H1'
  H1.gr2$hap <- 'H2'
  H2.gr1$hap <- 'H1'
  H2.gr2$hap <- 'H2'
  H1.gr <- c(H1.gr1, H1.gr2)
  H2.gr <- c(H2.gr1, H2.gr2)
  
  ## Get coressponding ideogram from the database
  if (genome %in% c('hg19', 'hg38') & is.null(custom.ideo)) {
    suppressMessages( ideogramCyto <- biovizBase::getIdeogram(genome, cytobands = TRUE) )
  } else if (!is.null(custom.ideo) & file.exists(custom.ideo)) {
    ideogramCyto <- utils::read.table(custom.ideo, stringsAsFactors = FALSE)
    colnames(ideogramCyto) <- c("seqnames","start","end", "name", "gieStain")
    ideogramCyto <- makeGRangesFromDataFrame(ideogramCyto, keep.extra.columns = TRUE)
  }
  
  ## Remove chromosomes that are not defined in chromosomes
  if (!is.null(chromosomes)) {
    ideogramCyto <- keepSeqlevels(ideogramCyto, value = chromosomes, pruning.mode = 'coarse')
    seqlevels(ideogramCyto) <- chromosomes

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
  ideo.df <- as.data.frame(ideogramCyto)
  
  ## Get chromosome breaks and labels
  max.len <- signif(max(ideo.df $end), digits = 2)
  breaks <- seq(from = 0, to = max.len, length.out = 6)
  labels <- breaks / 1000000
  labels <- paste0(labels, 'Mb')
  chr.num <- length(unique(ideo.df $seqnames))
  
  ## Plot the ideogram using ggplot2
  ideo <- ggplot() + 
    geom_rect(data=ideo.df, aes(ymin=start, ymax=end, xmin=0, xmax=1, fill=gieStain), color='black', show.legend=FALSE) + 
    ggbio::scale_fill_giemsa() +
    facet_grid(. ~ seqnames, switch = 'x') + 
    scale_y_continuous(breaks = breaks, labels = labels, expand = c(0,0)) +
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
    ggnewscale::new_scale("fill") + 
    geom_rect(data=plt.df, aes(ymin=start, ymax=end, xmin=xmin, xmax=xmax, fill=hap)) +
    scale_fill_manual(values = c('deepskyblue', 'gold'), name="")
  
  ## Decrese spacing between facets  
  plt <- plt + theme(panel.spacing.x=unit(0.5, "lines") , panel.spacing.y=unit(1,"lines"))
  
  ## Calculate concordance
  H1.acc <- getPhasingAcc(H1.tsv, chromosomes = chromosomes)
  H2.acc <- getPhasingAcc(H2.tsv, chromosomes = chromosomes)
  H1.err <- (sum(H1.acc$incorrect) / (sum(H1.acc$correct) + sum(H1.acc$incorrect))) * 100
  H2.err <- (sum(H2.acc$incorrect) / (sum(H2.acc$correct) + sum(H2.acc$incorrect))) * 100
  
  subtitle <- paste(paste0("Wrongly assigned segments to H1: ", round(H1.err, digits = 3), "%"),
                    paste0("Wrongly assigned segments to H2: ", round(H2.err, digits = 3), "%"),
                    sep = '\n')
  plt <- plt + labs(subtitle = subtitle)
  
  ## Add title
  if (!is.null(title)) {
    plt <- plt + ggtitle(title)
  }
  return(plt)
}
