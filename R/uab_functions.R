#' Load a BED file into a \code{\link{GRanges-class}} object.
#' 
#' This function will take a BED file of contig alignments to a reference genome and converts them 
#' into a \code{\link{GRanges-class}} object. This function also aims to collapse single unique contigs
#' tha are splitted in multipe pieces after the alignment to the reference.
#'
#' @param bedfile A BED file of contig alignments to a reference genome.
#' @param index A unique identifier to be added as an 'ID' field. 
#' @param min.mapq A minimum mapping quality of alignments reported in submitted BED file.
#' @param min.align A minimum length of an aligned sequence to a reference genome.
#' @param min.ctg.size A minimum length a final contig after gaps are collapsed.
#' @param max.gap A maximum length of a gap within a single contig alignments to be collapsed.
#' @return A \code{\link{GRanges-class}} object.
#' @author David Porubsky
#' 
bed2ranges <- function(bedfile=NULL, index=NULL, min.mapq=10, min.align=10000, min.ctg.size=500000, max.gap=100000) {
  bed.df <- utils::read.table(file = bedfile, header = FALSE, stringsAsFactors = FALSE)
  ## Add column names
  colnames(bed.df) <- c('seqnames','start','end','ctg','mapq','strand')[1:ncol(bed.df)]
  ## Add index if defined
  if (!is.null(index) & is.character(index)) {
    bed.df$ID <- index
  }
  ## Filter by mapping quality
  if (min.mapq > 0) {
    bed.df <- bed.df[bed.df$mapq >= min.mapq,]
  }
  if (nrow(bed.df) == 0) {
    stop("None of the BED alignments reach user defined mapping quality (min.mapq) !!!")
  }
  ## Convert data.frame to GRanges object
  bed.gr <- GenomicRanges::makeGRangesFromDataFrame(bed.df, keep.extra.columns = TRUE)
  ## Ignore strand
  strand(bed.gr) <- '*'
  ## Filter out small alignments
  if (min.align > 0) {
    bed.gr <- bed.gr[width(bed.gr) >= min.align]
  }
  ## Split ranges per contig
  bed.grl <- split(bed.gr, bed.gr$ctg)
  ## Make sure data frame is sorted
  #bed.gr <- GenomicRanges::sort(bed.gr, ignore.strand=TRUE)
  ## Collapse merge splitted continuous alignments of the same contig
  #bed.gr <- primatR::collapseBins(bed.gr, id.field = 1)
  #bed.grl <- endoapply(bed.grl, function(gr) primatR::collapseBins(gr, id.field = 1))
  bed.grl <- S4Vectors::endoapply(bed.grl, function(gr) fillGaps(gr=gr, max.gap=max.gap))
  bed.gr <- unlist(bed.grl, use.names = FALSE)
  ## Filter out small contigs after contig concatenation
  if (min.ctg.size > 0) {
    bed.gr <- bed.gr[width(bed.gr) >= min.ctg.size]
  }
  return(bed.gr)
}


#' This function will takes in a \code{\link{GRanges-class}} object of alignments of a single contigs to a reference
#' and collapses gaps in alignment based user specified maximum allowed gap.
#'
#' @param gr A \code{\link{GRanges-class}} object of regions of a single contigs aligned to a reference.
#' @param max.gap A maximum length of a gap within a single contig alignments to be collapsed.
#' @return A \code{\link{GRanges-class}} object..
#' @author David Porubsky
#' 
fillGaps <- function(gr, max.gap=100000) {
  gr <- GenomeInfoDb::keepSeqlevels(gr, value = as.character(unique(GenomeInfoDb::seqnames(gr))), pruning.mode = 'coarse')
  gr <- GenomicRanges::sort(gr)
  gap.gr <- GenomicRanges::gaps(gr, start = min(start(gr)))
  gap.gr <- gap.gr[width(gap.gr) <= max.gap]
  if (length(gap.gr) > 0) {
    red.gr <- GenomicRanges::reduce(c(gr[,0], gap.gr))
    mcols(red.gr) <- mcols(gr)[length(gr),]
  } else {
    red.gr <- GenomicRanges::reduce(gr)
    mcols(red.gr) <- mcols(gr)[length(gr),]
  }  
  return(red.gr)
}


#' Predict universal assembly breaks (UAB's) in a set of de novo assemblies
#'
#' This function will takes in a \code{\link{GRanges-class}} object of aligned contigs to a reference genome from two
#' and more de novo assemblies and finds positions where de novo assembly breaks in at least two independent assemblies.
#'
#' @param gr A \code{\link{GRanges-class}} object of contigs aligned to the reference.
#' @param chromosomes A user defined set of chromosomes to be analyzed.
#' @param binsize A binsize to count overlapping breaks in assemblies in.
#' @param stepsize A stepsize to smooth count signal and increase breakpoint resolution (default: binsize/2)
#' @param bsgenome A \code{\link{BSgenome-class}} object of reference genome used for binning step.
#' @param ID.col An unique identifier to split submitted genomic ranges that belong to different assembly and haplotype, 
#' recommended format: <sample>.<haplotype>.<assembly.ID>. 
#' @return A \code{\link{GRanges-class}} object of universal assembly breakpoints that appeat at least one in at least two 
#' independed assemblies based on <assembly.ID>.
#' @author David Porubsky
#' 
getUABs <- function(gr, chromosomes=NULL, binsize=500000, stepsize=500000/2, bsgenome=NULL, ID.col=NULL) {
  ## Extract start and end position of each contig as a unique breakpoint
  if (!is.null(ID.col) & is.character(ID.col)) {
    start.gr <- GRanges(seqnames=seqnames(gr), ranges=IRanges(start=start(gr), end=start(gr)), ID=unlist(mcols(gr)[eval(ID.col)]))
    end.gr <- GRanges(seqnames=seqnames(gr), ranges=IRanges(start=end(gr), end=end(gr)), ID=unlist(mcols(gr)[eval(ID.col)]))
    breaks.gr <- c(start.gr, end.gr)
    breaks.gr <- keepSeqlevels(breaks.gr, value = chromosomes, pruning.mode = 'coarse')
  } else {
    stop("Please specify column to use as an ID in 'ID.col', required!!!")
  }  
  
  ## Remove breaks at the very end of each chromosome
  if (!is.null(ID.col) & is.character(ID.col)) {
    grl <- split(gr, mcols(gr)[eval(ID.col)])
    mask.ranges <- endoapply(grl, range)   
    mask.ranges <- unlist(mask.ranges, use.names = FALSE)
    mask.start.gr <- GRanges(seqnames=seqnames(mask.ranges), ranges=IRanges(start=start(mask.ranges), end=start(mask.ranges)))
    mask.end.gr <- GRanges(seqnames=seqnames(mask.ranges), ranges=IRanges(start=end(mask.ranges), end=end(mask.ranges)))
    mask.gr <- c(mask.start.gr, mask.end.gr)
    breaks.gr <- subsetByOverlaps(breaks.gr, mask.gr, invert = TRUE)
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
  
  ## Bin genome into user defined bins
  if (!is.null(bsgenome)) {
    bins.gr <- makeBins(bsgenome = bsgenome, chromosomes = chromosomes, binsize = binsize, stepsize = stepsize)
  } else {
    stop("Please specify parameter 'bsgenome', required!!!")
  }
  
  ## In each bin calculate number contig breaks for tested assemblies
  mcols(breaks.gr) <- tidyr::separate(as.data.frame(mcols(breaks.gr)), col = 'ID', sep="\\.", into=c('sample','hap','tech'))
  breaks.grl <- split(breaks.gr, breaks.gr$tech)
  for (i in seq_along(breaks.grl)) {
    sub.breaks.gr <- breaks.grl[[i]]
    tech.id <- unique(sub.breaks.gr$tech)
    ## Break == 1, no break == 0
    b.counts <- countOverlaps(bins.gr, sub.breaks.gr)
    b.counts[b.counts > 0] <- 1
    mcols(bins.gr)[tech.id] <- b.counts
  }
  bins.gr$sums <- rowSums(as.data.frame(mcols(bins.gr)))
  bins.gr$uab <- bins.gr$sums > 1
  
  ## Keep only bins where bith canu and peregrine break in at least one haplotype
  roi.gr <- bins.gr[bins.gr$uab == TRUE]
  roi.gr <- reduce(roi.gr)
  ## Assign contig breaks to each region of interest
  hits <- findOverlaps(roi.gr, breaks.gr)
  breaks.per.roi <- split(breaks.gr[subjectHits(hits)], queryHits(hits))
  ## Get range between first and last contig break in each region of interest
  breaks.roi.ranges <- endoapply(breaks.per.roi, range)
  breaks.roi.ranges <- unlist(breaks.roi.ranges, use.names = FALSE)
  breaks.roi.ranges$midpoint <- start(breaks.roi.ranges) + (end(breaks.roi.ranges) - start(breaks.roi.ranges)) / 2 
  seqlengths(breaks.roi.ranges) <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)[seqlevels(breaks.roi.ranges)]
  
  return(breaks.roi.ranges)
}


#' Plot genome-wide visualisation of UAB's along with positions of aligned contigs to the reference.
#'
#' This function takes two \code{\link{GRanges-class}} objects; one that contains aligned contigs to a reference genome and second
#' that contains prediced universal assebly breaks and visualize them on the reference genome ideogram in the form of lollipops.
#'
#' @param gr A \code{\link{GRanges-class}} object of contigs aligned to the reference.
#' @param breaks.gr A \code{\link{GRanges-class}} object of predicted UAB's.
#' @param chromosomes A user defined set of chromosomes to be analyzed.
#' @param genome.ideo A reference genome id (e.g 'hg38') to plot the ideogram for. 
#' @param colors A set of colors used to distinguish contigs coming from different de novo assemblies.
#' @param annot.SD A data frame of positions of segmental duplication in the chosen reference genome (e.g. 'hg38') 
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object.
#' @author David Porubsky
#' 
plotUABlollipop <- function(gr=NULL, breaks.gr=NULL, ID.col=NULL, chromosomes=NULL, genome.ideo='hg38', colors=NULL, annot.SD=NULL) {
  ## Split Genomic Ranges by user defined identifier
  if (!is.null(ID.col) & is.character(ID.col)) {
    grl <- split(gr, mcols(gr)[eval(ID.col)])
  } else {
    stop("Please specify column to use as an ID in 'ID.col', required!!!")
  }  
  
  ## Add level information for plotting
  min.level <- 2
  plt.grl <- GRangesList()
  for (i in seq_along(grl)) {
    sub.gr <- grl[[i]]
    sub.gr <- sort(sub.gr)
    sub.gr$level <- rep(c(min.level, min.level + 1), length(sub.gr))[1:length(sub.gr)]
    plt.grl[[i]] <- sub.gr
    min.level <- min.level + 2
  }
  plt.gr <- unlist(plt.grl, use.names = FALSE)
  ## Prepare data for plotting
  plt.df <- as.data.frame(plt.gr)
  plt.df$ID <- factor(plt.df[,ID.col], levels=unique(plt.df[,ID.col]))
  
  ## Get coressponding ideogram from the database
  suppressMessages( ideogramCyto <- getIdeogram(genome.ideo, cytobands = TRUE) )
  ## Remove chromosomes that are not defined in chromosomes
  ideogramCyto <- keepSeqlevels(ideogramCyto, value = chromosomes, pruning.mode = 'coarse')
  seqlevels(ideogramCyto) <- chromosomes
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
                          panel.spacing.x=unit(1, "lines"))
  ## Convert ideogram bands stored in GRanges object into the data.frame
  ideo.df <- as.data.frame(ideogramCyto)
  
  ## Get chromosome breaks and labels
  max.len <- signif(max(ideo.df$end), digits = 2)
  breaks <- seq(from = 0, to = max.len, length.out = 6)
  labels <- breaks / 1000000
  labels <- paste0(labels, 'Mb')
  chr.num <- length(unique(ideo.df $seqnames))
  
  ## Plot the ideogram using ggplot2
  ideo <- ggplot() + 
    geom_rect(ideo.df, aes(ymin=start, ymax=end, xmin=0, xmax=0.9, fill=gieStain), color='black', show.legend=FALSE) + 
    scale_fill_giemsa() +
    facet_grid(. ~ seqnames, switch = 'x') + 
    scale_y_continuous(breaks = breaks, labels = labels) +
    theme_vertical +
    xlab("") +
    ylab("")
  
  ## Remove chromosomes that are not defined in chromosomes
  plt.df <- plt.df[plt.df$seqnames %in% chromosomes,]
  plt.df$seqnames <- factor(plt.df$seqnames, levels = chromosomes)
  ## Get the mid position of each genomic ranges
  plt <- ideo + 
    new_scale("fill") +
    geom_rect(data=plt.df, aes(ymin=start, ymax=end, xmin=level-1, xmax=level, fill=ID), inherit.aes=FALSE) +
    #scale_fill_manual(values = c('HG00733.H1'='cadetblue4','HG00733.H2'='cadetblue2','NA12878.H1'='darkgoldenrod4','NA12878.H2'='darkgoldenrod2')) +
    facet_grid(. ~ seqnames, switch = 'x')
  
  if (!is.null(colors) & is.character(colors)) {
    plt <- plt + scale_fill_manual(values = colors, name="")
  }
  
  if (!is.null(annot.SD) & is.data.frame(annot.SD)) {
    plt <- plt + new_scale("fill") +
      geom_rect(data=annot.SD, aes(ymin=start, ymax=end, xmin=0, xmax=-1), fill='red', inherit.aes=FALSE) +
      facet_grid(. ~ seqnames, switch = 'x')
  }  
  
  ## Add regions where contigs break in both assemblies
  if (!is.null(breaks.gr) & class(breaks.gr)[1] == "GRanges") {
    breaks.df <- as.data.frame(breaks.gr)
    breaks.df$annot <- factor(breaks.df$annot, levels=unique(breaks.df$annot))
    plt <- plt +
      geom_segment(data=breaks.df, aes(x=0, xend=12, y=midpoint, yend=midpoint), inherit.aes = FALSE, alpha=0.25) +
      geom_point(data=breaks.df, aes(x=12, y=midpoint, color=annot), inherit.aes = FALSE) +
      scale_color_manual(values = brewer.pal(n = 9, name = 'Set1'), name="") +
      facet_grid(. ~ seqnames, switch = 'x')
  }  
  return(plt)
}


#' Plot genome-wide positions of UAB's onto the chosen reference genome.
#'
#' This function takes a \code{\link{GRanges-class}} objects that contains prediced universal assebly breaks and 
#' visualize them on the reference genome as red squares.
#'
#' @param breaks.gr A \code{\link{GRanges-class}} object of predicted UAB's.
#' @param chromosomes A user defined set of chromosomes to be analyzed.
#' @param bsgenome A \code{\link{BSgenome-class}} object of reference genome to get the chromosome lengths from. 
#' @param centromeres.gr A \code{\link{GRanges-class}} object of centromere positions in the chosen reference genome.
#' @param gaps.gr A \code{\link{GRanges-class}} object of gap positions in the chosen reference genome.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object.
#' @author David Porubsky
#'
plotUABgenome <- function(breaks.gr, chromosomes=NULL, bsgenome=NULL, centromeres.gr=NULL, gaps.gr=NULL) {
  ## Load BSgenome object
  if (class(bsgenome) != 'BSgenome') {
    if (is.character(bsgenome)) {
      bsgenome <- tryCatch({
        suppressPackageStartupMessages(library(bsgenome, character.only=TRUE))
        bsgenome <- eval(parse(text=bsgenome)) ## replacing string by object
      }, error = function(e) {return(NULL)})
    }
  }
  
  ## Keep only user defined chromosomes
  chroms.in.data <- unique(seqnames(breaks.gr))
  if (is.null(chromosomes)) {
    chromosomes <- chroms.in.data
  }
  chroms2use <- intersect(chromosomes, chroms.in.data)
  breaks.gr <- keepSeqlevels(breaks.gr, value = chroms2use, pruning.mode = 'coarse')
  
  ## Get ideogram
  if (!is.null(bsgenome)) {
    seq.len <- seqlengths(bsgenome)[chroms2use]
    ideo.df <- data.frame(seqnames=names(seq.len), length=seq.len)
    ideo.df$seqnames <- factor(ideo.df$seqnames, levels=chroms2use)
    ideo.df$levels <- 1:length(seq.len)
  } else {
    stop("Please specify parameter 'bsgenome', required!!!")
  }
  
  ## Plot ideogram
  ideo <- ggplot() + geom_rect(data=ideo.df, aes(xmin=0, xmax=length, ymin=0, ymax=1), fill='black') 
  
  ## Add centromere and gap annot
  if (!is.null(centromeres.gr) & class(centromeres.gr)[1] == "GRanges") {
    cent.df <- as.data.frame(reduce(centromeres.gr))
    ideo <- ideo + geom_rect(data=cent.df, aes(xmin=start, xmax=end, ymin=0, ymax=1), fill='gray')
  }  
  if (!is.null(gaps.gr) & class(gaps.gr)[1] == "GRanges") {
    gaps.df <- as.data.frame(reduce(gaps.gr))
    ideo <- ideo + geom_rect(data=gaps.df, aes(xmin=start, xmax=end, ymin=0, ymax=1), fill='white')
  }
  
  ## Add assembly breakpoints
  plt.df <- as.data.frame(breaks.gr)
  plt <- ideo +
    geom_point(data=plt.df, aes(x = midpoint, y = 0.5), color='red', size=3, shape=18) +
    facet_grid(seqnames ~ ., switch = 'y') +
    scale_x_continuous(expand = c(0,0)) +
    theme_void() +
    theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank()) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(strip.text.y = element_text(angle = 180))
  
  return(plt)
}


#' Export FASTA sequences from a set of Genomic Ranges.
#'
#' This function takes a \code{\link{GRanges-class}} object and extracts a genomic sequence 
#' from these regions either from an original range or a range expanded on each side by
#' defined number of bases.
#' 
#' @param gr A \code{\link{GRanges-class}} object of genomic regions to extract genomic sequence from. 
#' @param bsgenome A \code{\link{BSgenome-class}} object of reference genome to get the genomic sequence from. 
#' @param asm.fasta An assembly FASTA file to extract DNA sequence determined by 'gr' parameter.
#' @param expand Expand ends of each genomic region by this length (in bp).
#' @param fasta.save A path to a filename where to store final FASTA file.
#' @author David Porubsky
#'
regions2FASTA <- function(gr, bsgenome=NULL, asm.fasta=NULL, expand=0, fasta.save=NULL) {
  ## Load BSgenome object
  if (class(bsgenome) != 'BSgenome') {
    if (is.character(bsgenome)) {
      bsgenome <- tryCatch({
        suppressPackageStartupMessages(library(bsgenome, character.only=TRUE))
        bsgenome <- eval(parse(text=bsgenome)) ## replacing string by object
      }, error = function(e) {return(NULL)})
    }
  }
  
  if (!is.null(asm.fasta)) {	
    ## Check if submitted fasta file is indexed
    asm.fasta.idx <- paste0(asm.fasta, ".fai")
    if (!file.exists(asm.fasta.idx)) {
      fa.idx <- Rsamtools::indexFa(file = asm.fasta)
    }
  }
  
  ## Expand regions of interest by certain size (downstream and upstream)
  if (expand > 0) {
    gr.seqLen <- seqlengths(gr)[as.character(seqnames(gr))]
    start(gr) <- pmax(1, start(gr) - expand)
    end(gr) <- pmin(gr.seqLen, end(gr) + expand)
  }
  
  if (!is.null(bsgenome)) {
    ## Extract FASTA from BSgenome object
    gr.seq <- BSgenome::getSeq(bsgenome, gr)
    names(gr.seq) <- as.character(gr)
  } else if (is.character(asm.fasta)) {
    ## Extract FASTA from user defined FASTA file
    fa.file <- open(Rsamtools::FaFile(asm.fasta))
    ## Remove sequences not present in the FASTA index
    fa.idx <- Rsamtools::scanFaIndex(fa.file)
    gr <- suppressWarnings( subsetByOverlaps(gr, fa.idx) )
    ## Read in contigs for a given cluster
    gr.seq <- Rsamtools::scanFa(file = fa.file, param = gr, as = "DNAStringSet")
    names(gr.seq) <- as.character(gr)
  } else {
    stop("Please set a 'bsgenome' or 'asm.fasta' parameter!!!")
  }

  ## Write final FASTA
  if (is.character(fasta.save)) {
    Biostrings::writeXStringSet(x = gr.seq, filepath = fasta.save, format = 'fasta')
  } else {
    warning("Please speficify 'fasta.save' if you want to export FASTA into a file!!!")
  }  
  return(gr.seq)
}


# concatAlign <- function(gr, max.gap=100000) {
#   print(gr)
#   gr <- keepSeqlevels(gr, value = as.character(unique(seqnames(gr))), pruning.mode = 'coarse')
#   gr$merge <- lead(width(gaps(gr)), default = Inf) < max.gap | lag(width(gaps(gr)), default = Inf) < max.gap
#   gr <- primatR::collapseBins(gr, id.field = ncol(mcols(gr)))
#   return(gr)
# }
