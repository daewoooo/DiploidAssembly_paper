## Required libraries to load
library(diffHic)
library(ggplot2)
library(scales)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(cowplot)
library(InteractionSet)


#' Costruct a contact/interation matrix from a Hi-C BAM file.
#'
#' Function to construct a contact matrix using a BAM aligned reads resulting from proximity-ligation experiment (Hi-C).
#'
#' @param bamfile A name sorted BAM file containing aligned reads from an Hi-C experiment.
#' @param outputdirectory A path to a directory to store the results.
#' @param bsgenome A BSgenome object used to predict restriction enzyme cut sites.
#' @param chromosomes A selected chromosomes/scaffolds to be processed.
#' @param restriction.site A nucleotide sequence to define restriction sites in 'bsgenome'.
#' @param overhang A number of bases to that creates an overhang at the defined restriction sites.
#' @param min.mapq A minimum mapping quality for a read-pair to be included.
#' @param min.counts A minimum number of interactions (read-pairs) between two genomic bins.
#' @param min.seq.len A minimum chromosome/scaffold length to be considered (default: 10Mb).
#' @param dedup Set to \code{TRUE} if duplicate reads should be removed.
#' @param bin.size A size of the genomic bin to count interactions within.
#' @param use.existing.h5 Set to \code{TRUE} if the existing .h5 file, for a submitted BAM file, should be used.
#' @param n.pairs.chunk ...
#' @return A \code{\link{InteractionSet-class}} object with Hi-C interaction counted in specific genomic bin.
#' @author David Porubsky
#' 
bam2contactMatrix <- function(bamfile=NULL, outputdirectory=NULL, bsgenome=NULL, chromosomes=NULL, restriction.site=NULL, overhang=NULL, min.mapq=10, min.counts=1, min.seq.len=10000000, dedup=TRUE, bin.size=1e6, use.existing.h5=TRUE, n.pairs.chunk=1000L) {
  ## Check if output directory is defined
  if (!is.null(outputdirectory) & is.character(outputdirectory)) {
    if (!dir.exists(outputdirectory)) {
      dir.create(outputdirectory)
    }
  } else {
    outputdirectory <- dirname(bamfile)
    warning("Parameter 'outputdirectory' is not defined, will write into 'bamfile' location!!!")
  }
    
  ## Load BSgenome object
  if (class(bsgenome) != 'BSgenome') {
    if (is.character(bsgenome)) {
      bsgenome <- tryCatch({
        suppressPackageStartupMessages(library(bsgenome, character.only=TRUE))
        bsgenome <- eval(parse(text=bsgenome)) ## replacing string by object
      }, error = function(e) {return(NULL)})
    }
  }
  
  ## Get contigs/scaffolds/chromosomes names and sizes
  if (!is.null(bamfile) & is.character(bamfile)) {
    file.header <- Rsamtools::scanBamHeader(bamfile)[[1]]
    chrom.lengths <- file.header$target
  } else if (class(bsgenome) == 'BSgenome') {
    chrom.lengths <- seqlengths(bsgenome)
  }
  
  ## Order chromosome/scaffolds by their length
  chrom.lengths <- chrom.lengths[order(chrom.lengths, decreasing = TRUE)]
  
  ## Cut the reference genome at user defined restriction sites
  if (!is.null(restriction.site) & is.character(restriction.site)) {
    message("Cutting genome at restriction sites ...", appendLF=FALSE); ptm <- proc.time()
    if (nchar(restriction.site) & (is.numeric(overhang) & overhang > 0)) {
      res.frags <- diffHic::cutGenome(bsgenome, pattern = restriction.site, overhang = overhang)
    } else {
      res.frags <- diffHic::cutGenome(bsgenome, pattern = restriction.site)
    }
    time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
    ## Keep only chromosomes that are present in the submitted bamfile
    keep.chr <- GenomeInfoDb::seqlevels(res.frags)[GenomeInfoDb::seqlevels(res.frags) %in% names(chrom.lengths)]
    if (length(keep.chr) < length(GenomeInfoDb::seqlevels(res.frags))) {
      num.removed <- length(GenomeInfoDb::seqlevels(res.frags)) - length(keep.chr)
      warning(paste0("    ", num.removed, " chromosomes from restriction map of " ,attributes(bsgenome)$pkgname ," were not present in the 'bamfile' and have been removed!!!"))
      res.frags <- GenomeInfoDb::keepSeqlevels(res.frags, value = keep.chr, pruning.mode = 'coarse')
    }  
  } else {
    ## If no restriction site is submitted use empty GRanges object using the 
    ## contigs/scaffolds/chromosomes names and sizes from the submitted bamfile
    warning("No restriction site defined => assuming DNase HiC experiment")
    res.frags <- GenomicRanges::GRanges()
    GenomeInfoDb::seqlevels(res.frags) <- names(chrom.lengths)
    GenomeInfoDb::seqlengths(res.frags) <- chrom.lengths
  }
  
  ## Convert 'res.frags' to object that specify read pair loading parameters
  res.frags.param <- diffHic::pairParam(res.frags)
  
  ## Select only specific user-defined chromosomes
  if (!is.null(chromosomes) & is.character(chromosomes)) {
    if (all(chromosomes %in% GenomeInfoDb::seqlevels(res.frags))) {
      res.frags.param <- diffHic::reform(res.frags.param, restrict=chromosomes)
      ## Subset chrom.lengths in use to those defined in chromosomes parameter
      chrom.lengths <- chrom.lengths[names(chrom.lengths) %in% chromosomes]
    } else {
      warning("Not all chromosomes defined in submitted 'bamfile', skipping chromosomes selection!!!")
    } 
  }  
  
  ## Keep only chromosome/scaffold of a certain length
  if (!is.null(min.seq.len) & min.seq.len > 0) {
    chrom.lengths <- chrom.lengths[chrom.lengths >= min.seq.len]
    res.frags.param <- diffHic::reform(res.frags.param, restrict=names(chrom.lengths))
  }
  
  ## Get paired Hi-C data in h5 file
  outfile.h5 <- paste0(basename(bamfile), '.h5')
  destination.h5 <- file.path(outputdirectory, outfile.h5)
  if (use.existing.h5) {
    if (file.exists(destination.h5)) {
      message(paste0("Using previously exported read pairs in .h5 file: ", destination.h5))
    } else {
      message("Extracting read pairs [Time required: high] ...", appendLF=FALSE); ptm <- proc.time()
      diagnostic <- diffHic::preparePairs(bam = bamfile, param = res.frags.param, file = destination.h5, dedup = dedup, minq = min.mapq, storage = n.pairs.chunk)
      time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
      message("    Total number of processed pairs: ", diagnostic$pairs[1])
    }
  } else {
    message("Extracting read pairs [Time required: high] ...", appendLF=FALSE); ptm <- proc.time()
    diagnostic <- diffHic::preparePairs(bam = bamfile, param = res.frags.param, file = destination.h5, dedup = dedup, minq = min.mapq, storage = n.pairs.chunk)
    time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
    message("    Total number of processed pairs: ", diagnostic$pairs[1])
  }
  
  ## Count combinations for interactions between all pairs of bins
  message("Counting binned interactions [Time required: medium] ...", appendLF=FALSE); ptm <- proc.time()
  interaction.data <- diffHic::squareCounts(files = destination.h5, param = res.frags.param, width = bin.size, filter = min.counts)
  time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
  
  ## Save interaction set as RData object
  message("Exporting RData object  ...", appendLF=FALSE); ptm <- proc.time()
  outfile.rdata <- paste0(basename(bamfile), '.RData')
  destination <- file.path(outputdirectory, outfile.rdata)
  save(interaction.data, file = destination)
  time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
  
  return(interaction.data)
}


#' Plot contact/interation matrix constructed from a Hi-C BAM file.
#'
#' @param interaction.obj A \code{\link{InteractionSet-class}} object with Hi-C interaction counted in specific genomic bin.
#' @param genome.coord Set to \code{TRUE} if coordinates of each chromosome/scaffold should be transformed into the whole genome space.
#' @param chromosomes A selected chromosomes/scaffolds to be plotted.
#' @param title A charecter string to use as a title of the final plot.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object.
#' @author David Porubsky
#' 
plotContactMatrix <- function(interaction.obj=NULL, genome.coord=TRUE, chromosomes=NULL, title=NULL) {
  message("Plotting HiC interactions ...", appendLF=FALSE); ptm <- proc.time()
  ## Load interactionSet object
  interaction.data <- get(load(interaction.obj))
  
  ## Check if submited object is of 'InteractionSet' class
  if (class(interaction.data)[1] != 'InteractionSet') {
    stop("Aborting, submitted 'interaction.obj' is not of 'InteractionSet' class!!!")
  } 
  
  ## Get genomic bins for pairs of bins
  bins.first <- InteractionSet::anchors(interaction.data, type="first")
  bins.second <- InteractionSet::anchors(interaction.data, type="second")
  ## Report final HiC interaction object
  n.links <- SummarizedExperiment::assay(interaction.data)
  n.links <- n.links[,1]
  
  ## Order chromosomes/scaffolds by size (decreasing)
  GenomeInfoDb::seqlevels(bins.first) <- GenomeInfoDb::seqlevels(bins.first)[order(GenomeInfoDb::seqlengths(bins.first), decreasing = TRUE)]
  GenomeInfoDb::seqlevels(bins.second) <- GenomeInfoDb::seqlevels(bins.second)[order(GenomeInfoDb::seqlengths(bins.second), decreasing = TRUE)]
  
  if (!is.null(chromosomes) & is.character(chromosomes)) {
    mask <- GenomeInfoDb::seqnames(bins.first) %in% chromosomes & GenomeInfoDb::seqnames(bins.second) %in% chromosomes
    bins.first <- bins.first[mask]
    bins.second <- bins.second[mask]
    GenomeInfoDb::seqlevels(bins.first) <- chromosomes
    GenomeInfoDb::seqlevels(bins.second) <- chromosomes
    n.links <- n.links[as.vector(mask)]
  }
  
  ## Reorder interactions
  gi <- InteractionSet::GInteractions(bins.first, bins.second)
  gi <- InteractionSet::swapAnchors(gi)
  bins.first <- InteractionSet::anchors(gi, type="first")
  bins.second <- InteractionSet::anchors(gi, type="second")
  
  ## Convert chromosome/scaffold based coordinates into genome-wide coordinates
  if (genome.coord) {
    bins.first <- transCoord(bins.first)
    bins.second <- transCoord(bins.second)
  }  
  
  ## Contruct data.frame used for plotting
  if (genome.coord) {
    link.df <- data.frame('M1.seqnames'=as.character(GenomicRanges::seqnames(bins.first)),
                          'M1.start'=bins.first$start.genome,
                          'M1.end'=bins.first$end.genome,
                          'M2.seqnames'=as.character(GenomicRanges::seqnames(bins.second)),
                          'M2.start'=bins.second$start.genome,
                          'M2.end'=bins.second$end.genome,
                          value=n.links)
  } else {
    link.df <- data.frame('M1.seqnames'=as.character(GenomicRanges::seqnames(bins.first)),
                          'M1.start'=as.numeric(start(bins.first)),
                          'M1.end'=as.numeric(end(bins.first)),
                          'M2.seqnames'=as.character(GenomicRanges::seqnames(bins.second)),
                          'M2.start'=as.numeric(start(bins.second)),
                          'M2.end'=as.numeric(end(bins.second)),
                          value=n.links)
  }  
  
  ## Get chromosome/scaffold boundaries
  region.boundaries <- range(bins.first)
  region.boundaries <- keepSeqlevels(region.boundaries, value = as.character(seqnames(region.boundaries)), pruning.mode = 'coarse')
  region.boundaries <- transCoord(region.boundaries)
  
  ## Transform boundaries of each chromosome into genomewide coordinates (cumulative sum of bin coordintes)
  cum.seqlengths <- cumsum(as.numeric(GenomeInfoDb::seqlengths(region.boundaries)))
  cum.seqlengths.0 <- c(0, cum.seqlengths[-length(cum.seqlengths)])
  ## Get positions of ends of each chromosome to plot lines between the chromosomes
  if (length(cum.seqlengths) > 1) {
    chr.lines <- data.frame(y = c(0, cum.seqlengths))
  } else {
    chr.lines <- data.frame(y = 0)
  }    
  ## Get positions of each chromosomes names
  chr.label.pos <- round(cum.seqlengths.0 + (0.5 * GenomeInfoDb::seqlengths(region.boundaries)))
  
  ## Get y-labels of chromosome/scaffold sizes
  max.len <- signif(max(region.boundaries$end.genome), digits = 2)
  y.breaks <- seq(from = 0, to = max.len, length.out = 6)
  y.labels <- y.breaks / 1000000
  y.labels <- paste0(y.labels, 'Mbp')
  
  ## Remove outliers
  #outlier <- stats::quantile(link.df$value, 0.999)
  #set outlier bins to the limit
  #link.df$value[link.df$value >= outlier] <- outlier
  
  plt <- ggplot2::ggplot(link.df) +
    geom_rect(aes(xmin=M1.start, xmax=M1.end, ymin=M2.start, ymax=M2.end, fill=value)) +
    geom_vline(xintercept = chr.lines$y, linetype="dashed", color="red") +
    scale_x_continuous(breaks = chr.label.pos, labels = names(chr.label.pos)) +
    scale_y_continuous(breaks = y.breaks, labels = y.labels) +
    scale_fill_gradient(low = "white", high = "black", trans = 'log', name="# of Interactions (log)") +
    xlab("Chromosome/scaffold name") +
    ylab("Genomic position (Mbp)") +
    theme_bw() +
    theme(aspect.ratio=1) + #Make sure plot is always a square
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  
  ## Add title if defined
  if (is.character(title)) {
    plt <- plt + ggtitle(title)
  }
  time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
  
  ## Return final plot
  return(plt)
}


#' Plot statistics of read-pairs stored in Hi-C BAM file.
#'
#' @param read.pairs.stat.obj A \code{data.frame} object containing fragment, orientaton and insert size distribution.
#' @param h5.file ...
#' @param chromosomes ... 
#' @param bsgenome ...
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object.
#' @author David Porubsky
#' 
plotContactMatrixStat <- function(bamfile=NULL, h5.file=NULL, chromosomes=NULL, bsgenome=NULL) {
  ## Load BSgenome object
  if (class(bsgenome) != 'BSgenome') {
    if (is.character(bsgenome)) {
      suppressPackageStartupMessages(library(bsgenome, character.only=TRUE))
      bsgenome <- eval(parse(text=bsgenome)) ## replacing string by object
    }
  }
  
  ## Get contigs/scaffolds/chromosomes names and sizes
  if (!is.null(bamfile) & is.character(bamfile)) {
    file.header <- Rsamtools::scanBamHeader(bamfile)[[1]]
    chrom.lengths <- file.header$target
  } else if (class(bsgenome) == 'BSgenome') {
    chrom.lengths <- GenomeInfoDb::seqlengths(bsgenome)
  }
  
  res.frags <- GenomicRanges::GRanges()
  GenomeInfoDb::seqlevels(res.frags) <- names(chrom.lengths)
  GenomeInfoDb::seqlengths(res.frags) <- chrom.lengths

  ## Convert 'res.frags' to object that specify read pair loading parameters
  res.frags.param <- diffHic::pairParam(res.frags)

  ## Select only specific user-defined chromosomes
  if (!is.null(chromosomes) & is.character(chromosomes)) {
    if (all(chromosomes %in% GenomeInfoDb::seqlevels(res.frags))) {
      res.frags.param <- diffHic::reform(res.frags.param, restrict=chromosomes)
    } else {
      warning("Not all chromosomes defined in submitted 'bamfile', skipping chromosomes selection!!!")
    } 
  } 
  
  ## Calculate read pair statistics (fragment length, orientation and insert size)
  message("Calculating read pair statistics [Time required: medium] ...", appendLF=FALSE); ptm <- proc.time()
  stat.data <- diffHic::getPairData(file = h5.file, param = res.frags.param)
  ## Save read pair statistics as RData object
  #outfile.rdata <- paste0(basename(bamfile), '.hicPairsStat.RData')
  #destination <- file.path(outputdirectory, outfile.rdata)
  #save(read.stat, file = destination)
  time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
  
  message("Plotting HiC read-pair statistics ...", appendLF=FALSE); ptm <- proc.time()

  ## Plot proportion of intra- versus inter-chromosomal links
  inter <- length(stat.data$insert[is.na(stat.data$insert)])
  intra <- length(stat.data$insert[!is.na(stat.data$insert)])
  plt.df <- data.frame(counts=c(inter, intra), labels=c('inter.chr','intra.chr'))
  plt.df$link.frac <- (plt.df$counts / sum(plt.df$counts)) * 100
  ## Make the plot
  plt1 <- ggplot2::ggplot(plt.df) + 
    geom_col(aes(x='HiC interactions', y=link.frac, fill=labels)) +
    scale_fill_manual(values = c('chocolate2', 'chartreuse4'), name="") +
    xlab("") +
    ylab("% of intra vs inter links") +
    theme_bw() 
  
  ## Plot histogram of insert size distributions
  h.counts <- hist(stat.data$insert, plot = FALSE, breaks = 100)
  plt.df <- data.frame(counts=h.counts$counts, xpos=h.counts$mids)
  ## Prepare y-axis breaks
  max.count <- signif(max(plt.df$counts), digits = 2)
  y.breaks <- 10^(1:12)
  y.breaks <- y.breaks[y.breaks < max.count]
  ## Prepare x-axis breaks
  max.len <- signif(max(plt.df$xpos), digits = 2)
  x.breaks <- seq(from = 0, to = max.len, length.out = 6)
  x.labels <- x.breaks / 1000000
  x.labels <- paste0(x.labels, 'Mbp')
  ## Make the plot
  plt2 <- ggplot2::ggplot(plt.df) + 
    geom_col(aes(x=xpos, y=counts)) +
    geom_hline(yintercept = y.breaks, linetype='dashed', color='white') +
    scale_x_continuous(breaks = x.breaks, labels = x.labels, name = "Interaction distance (Mbp)") +
    scale_y_continuous(trans = 'log10', breaks = y.breaks, labels = comma, name = "Frequency (log10)") + 
    theme_bw()
  
  ## Categorize insert sizes by their length
  insert.categs <- findInterval(stat.data$insert, vec = c(1000, 10000, 100000, 1000000, 10000000))
  categ.counts <- as.numeric(table(insert.categs))
  categ.frac <- (categ.counts / sum(categ.counts)) * 100
  plt.df <- data.frame(counts=categ.counts, frac=categ.frac)
  labels <- c('<1kbp', '1kbp-10kbp','10kbp-100kbp','100kbp-1000kbp','1000kbp-10000kbp','>10000kbp')
  plt.df$labels <- factor(labels, levels = labels)
  ## Make the plot
  plt3 <- ggplot2::ggplot(plt.df) + 
    geom_col(aes(x=labels, y=frac)) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    xlab("Interaction size category (kbp)") +
    ylab("Fraction (%)")
  
  ## Construct final report
  final.plt <- plot_grid(plt1, plt2, plt3, nrow = 1, rel_widths = c(1.5,4,2), align = 'h')
  
  time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
  return(final.plt)
}


#' Transform genomic coordinates
#'
#' Add two columns with transformed genomic coordinates to the \code{\link{GRanges-class}} object. This is useful for making genomewide plots.
#'
#' @param gr A \code{\link{GRanges-class}} object.
#' @return The input \code{\link{GRanges-class}} with two additional metadata columns 'start.genome' and 'end.genome'.
#' @author David Porubsky
transCoord <- function(gr) {
  cum.seqlengths <- cumsum(as.numeric(GenomeInfoDb::seqlengths(gr)))
  cum.seqlengths.0 <- c(0,cum.seqlengths[-length(cum.seqlengths)])
  names(cum.seqlengths.0) <- GenomeInfoDb::seqlevels(gr)
  gr$start.genome <- start(gr) + cum.seqlengths.0[as.character(GenomicRanges::seqnames(gr))]
  gr$end.genome <- end(gr) + cum.seqlengths.0[as.character(GenomicRanges::seqnames(gr))]
  return(gr)
}
