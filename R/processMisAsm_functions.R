#' Load minimap2 alignments in PAF format
#' 
#' @param paf.file A PAF alignment of query sequence to the target sequence.
#' @param min.mapq Required minimum mapping quality of each alignment.
#' @param seqnames A sequence names to be plotted only.
#' @return A \code{data.frame} object of all PAF alignments.
#' @author David Porubsky
#' 
loadPAFalign <- function(paf.file, min.mapq=10, seqnames=NULL) {
  ## Load PAF minimap output
  paf.align <- read.table(paf.file, stringsAsFactors = FALSE)
  paf.align <- paf.align[,c(1:12)]
  colnames(paf.align) <- c('seqnames1', 'seq.len1', 'start1', 'end1', 'dir','seqnames2','seq.len2','start2','end2','n.match','aln.len','mapq')
  
  ## Filter contigs by mapping quality
  if (min.mapq > 0) {
    paf.align <- paf.align[paf.align$mapq >= min.mapq,]
  }
  
  ## Keep only user defined seqnames
  if (!is.null(seqnames) & is.character(seqnames)) {
    paf.align <- paf.align[paf.align$seqnames1 %in% seqnames | paf.align$seqnames2 %in% seqnames,]
  }
  return(paf.align)
}


#' Plot PAF alignments as a circular links
#' 
#' @param paf A \code{data.frame} of loaded PAF alignments.
#' @param genome.ideo Banded genome ideogram to obtaing target sequences from (e.g. hg83).
#' @param track.ranges  A \code{data.frame} of genomic ranges to plot as a annotation track.
#' @return A circlize plot object.
#' @author David Porubsky
#' 
plotCircosLinks <- function(paf=NULL, genome.ideo=NULL, track.ranges=NULL) {
  
  ## Construct bed files
  bed1 <- paf[,c(1,3,4,5)]
  bed1$dir <- dplyr::recode(bed1$dir, "+"="chartreuse4", "-"="firebrick4")
  bed2 <- paf[,c(6,8,9)]
  
  ## Restric ideogram to a certain region
  if (!is.null(genome.ideo)) {
    target.ideo <- genome.ideo[genome.ideo$seqnames %in% unique(bed2$seqnames2),]
    query.ideo <- data.frame(seqnames=unique(paf$seqnames1), start=1, end=unique(paf$seq.len1), name='ctg', gieStain='gneg')
    ideo <- rbind(target.ideo, query.ideo)
    circos.initializeWithIdeogram(ideo)
  } else {
    target.ideo <- data.frame(seqnames=unique(paf$seqnames2), start=1, end=unique(paf$seq.len2))
    query.ideo <- data.frame(seqnames=unique(paf$seqnames1), start=1, end=unique(paf$seq.len1))
    ideo <- rbind(target.ideo, query.ideo)
    circos.genomicInitialize(ideo)
  } 
  
  ## Construct circlize plot
  circos.genomicLink(bed1, bed2, col = bed1$dir, border = NA)
  if (!is.null(track.ranges)) {
    if (length(unique(sd.df.sub$seqnames)) > 1) {
      sd.df.sub.l <- split(sd.df.sub, as.character(sd.df.sub$seqnames))
      for (i in seq_along(sd.df.sub.l)) {
        sd.df.sub.sub <- sd.df.sub.l[[i]]
        circos.genomicRect(sd.df.sub.sub, ytop = -1, ybottom = 0, col='orange', border='orange', sector.index = as.character(unique(sd.df.sub.sub$seqnames)))
      }
    } else {
      circos.genomicRect(sd.df.sub, ytop = -1, ybottom = 0, col='orange', border='orange', sector.index = unique(bed2$seqnames2))
    }  
  }
  circos.clear() 
}