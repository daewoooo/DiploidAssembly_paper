#' Find meiotic recombination breakpoints using circular binary segmentation (CBS)
#'
#' @param child A \code{\link{GRanges-class}} object of SNV positions found in a child of a given family.
#' @param par1 A \code{\link{GRanges-class}} object of SNV positions found in a parent1 (mother or father) of a given family.
#' @param par2 A \code{\link{GRanges-class}} object of SNV positions found in a parent2 (mother or father) of a given family.
#' @return A \code{list} object 
#' @author David Porubsky
#'
getRecombMap <- function(child=NULL, par1=NULL, par2=NULL) {
  ## Get only shared SNVs between datasets
  comparison.all <- child
  comparison.all$par1.H1 <- 'N'
  comparison.all$par1.H2 <- 'N'
  comparison.all$par2.H1 <- 'N'
  comparison.all$par2.H2 <- 'N'
  
  shared.par1 <- findOverlaps(par1, child)
  comparison.all$par1.H1[subjectHits(shared.par1)] <- par1$H1[queryHits(shared.par1)]
  comparison.all$par1.H2[subjectHits(shared.par1)] <- par1$H2[queryHits(shared.par1)]
  
  shared.par2 <- findOverlaps(par2, child)
  comparison.all$par2.H1[subjectHits(shared.par2)] <- par2$H1[queryHits(shared.par2)]
  comparison.all$par2.H2[subjectHits(shared.par2)] <- par2$H2[queryHits(shared.par2)]
  
  comparisons.grl <- GRangesList()
  hap.segm.grl <- GRangesList()
  recomb.break.grl <- GRangesList()
  for (i in seq_along(chromosomes)) {
    chr <- chromosomes[i]
    message("Working on chromosome: ", chr)
    comparison.obj <- comparison.all[seqnames(comparison.all) == chr]
    ## Compare child.H1 to par1
    comparison.obj$c1_to_par1 <- 'N'
    mask <- comparison.obj$H1 != 'N' & comparison.obj$par1.H1 != 'N' & comparison.obj$par1.H2 != 'N' & comparison.obj$par1.H1 != comparison.obj$par1.H2
    comparison.obj$c1_to_par1[mask] <- ifelse(comparison.obj$H1[mask] == comparison.obj$par1.H1[mask], 'H1', 'H2')
    c1_to_par1.breaks <- length(rle(comparison.obj$c1_to_par1[mask])$lengths)
    ## Compare child.H2 to par1
    comparison.obj$c2_to_par1 <- 'N'
    mask <- comparison.obj$H2 != 'N' & comparison.obj$par1.H1 != 'N' & comparison.obj$par1.H2 != 'N' & comparison.obj$par1.H1 != comparison.obj$par1.H2
    comparison.obj$c2_to_par1[mask] <- ifelse(comparison.obj$H2[mask] == comparison.obj$par1.H1[mask], 'H1', 'H2')
    c2_to_par1.breaks <- length(rle(comparison.obj$c2_to_par1[mask])$lengths)
    ## Compare child.H1 to par2
    comparison.obj$c1_to_par2 <- 'N'
    mask <- comparison.obj$H1 != 'N' & comparison.obj$par2.H1 != 'N' & comparison.obj$par2.H2 != 'N' & comparison.obj$par2.H1 != comparison.obj$par2.H2
    comparison.obj$c1_to_par2[mask] <- ifelse(comparison.obj$H1[mask] == comparison.obj$par2.H1[mask], 'H1', 'H2')
    c1_to_par2.breaks <- length(rle(comparison.obj$c1_to_par2[mask])$lengths)
    ## Compare child.H2 to par2
    comparison.obj$c2_to_par2 <- 'N'
    mask <- comparison.obj$H2 != 'N' & comparison.obj$par2.H1 != 'N' & comparison.obj$par2.H2 != 'N' & comparison.obj$par2.H1 != comparison.obj$par2.H2
    comparison.obj$c2_to_par2[mask] <- ifelse(comparison.obj$H2[mask] == comparison.obj$par2.H1[mask], 'H1', 'H2')
    c2_to_par2.breaks <- length(rle(comparison.obj$c2_to_par2[mask])$lengths)
    
    ## Get inherited parental homologs
    if ((c1_to_par1.breaks + c2_to_par2.breaks) < (c2_to_par1.breaks + c1_to_par2.breaks)) {
      comparisons <- comparison.obj[,c('c1_to_par1', 'c2_to_par2')]
      H1.inherited <- unique(as.character(par1$sampleName))
      H2.inherited <- unique(as.character(par2$sampleName))
    } else {
      comparisons <- comparison.obj[,c('c1_to_par2', 'c2_to_par1')]
      H1.inherited <- unique(as.character(par2$sampleName))
      H2.inherited <- unique(as.character(par1$sampleName))
    }
    names(mcols(comparisons)) <- c('H1.comp', 'H2.comp')
    H1.recomb <- findRecomb(comparison = comparisons[,'H1.comp'], minSeg = 1000, smooth = 2, collapse.amb = TRUE)
    H2.recomb <- findRecomb(comparison = comparisons[,'H2.comp'], minSeg = 1000, smooth = 2, collapse.amb = TRUE)
    
    ## Add field of inherited homologs
    H1.recomb$hap.segm$inherited <- H1.inherited
    H1.recomb$recomb.break$inherited <- H1.inherited
    H2.recomb$hap.segm$inherited <- H2.inherited
    H2.recomb$recomb.break$inherited <- H2.inherited
    
    comparisons$H1.comp[comparisons$H1.comp != 'N'] <- paste0(comparisons$H1.comp[comparisons$H1.comp != 'N'], '.', H1.inherited)
    comparisons$H2.comp[comparisons$H2.comp != 'N'] <- paste0(comparisons$H2.comp[comparisons$H2.comp != 'N'], '.', H2.inherited)
    
    comparisons.grl[[i]] <- comparisons
    hap.segm.grl[[i]] <- c(H1.recomb$hap.segm, H2.recomb$hap.segm)
    recomb.break.grl[[i]] <- c(H1.recomb$recomb.break, H2.recomb$recomb.break)
  }  
  ## Return final data object
  comparisons.gr <- unlist(comparisons.grl, use.names = FALSE)
  comparisons.gr <- comparisons.gr[comparisons.gr$H1.comp != 'N' | comparisons.gr$H2.comp != 'N']
  hap.segm.gr <- unlist(hap.segm.grl, use.names = FALSE)  
  hap.segm.gr$ID <- paste0(hap.segm.gr$match, '.', hap.segm.gr$inherited)
  recomb.break.gr <- unlist(recomb.break.grl, use.names = FALSE)
  recomb.break.gr <- recomb.break.gr[width(recomb.break.gr) > 1]
  return(list(parental.snvs = comparisons.gr, parental.segments = hap.segm.gr, recombination.breakpoints = recomb.break.gr))
}  


#' Load SNVs obtained from Peter Audano's pipeline comparing de novo assembly in respect to the reference genome 
#' 
#' @param snv.tab A text file with SNV calls in respect to the reference genome.
#' @param rmdup Set to \code{TRUE} if duplicated SNV positions should be removed.
#' @param chromosomes Chromosomes names (e.g. 'chr1') that will be kept.
#' @param family Membership in a family, either 'child' or 'parent'.
#' @param sampleName A unique sample identifier.
#' @return A \code{\link{GRanges-class}} object with haplotype specific alleles.
#' @author David Porubsky
#'
loadSNVtable <- function(snv.tab=NULL, rmdup=TRUE, chromosomes=NULL, family='child', sampleName='') {
  ## Load SNV data
  snv.df <- data.table::fread(snv.tab, stringsAsFactors = FALSE)
  snv.gr <- GenomicRanges::makeGRangesFromDataFrame(snv.df, seqnames.field='#CHROM', start.field='POS', end.field='POS', keep.extra.columns = TRUE)
  ## Keep only user defined chromosome names
  if (!is.null(chromosomes) & is.character(chromosomes)) {
    snv.gr <- GenomeInfoDb::keepSeqlevels(snv.gr, value = chromosomes, pruning.mode = 'coarse')
  }  
  
  ## Remove duplicated SNV positions
  if (rmdup) {
    snv.gr <- snv.gr[!duplicated(snv.gr)]
  }
  
  if (family == 'child') {
    ## for child keep only HET SNVs
    snv.gr.H1 <- snv.gr[grep(snv.gr$HAP_SAMPLES, pattern = "^h1$")]
    snv.gr1 <- snv.gr.H1[,0]
    snv.gr1$H1 <- snv.gr.H1$ALT
    snv.gr1$H2 <- snv.gr.H1$REF
    snv.gr.H2 <- snv.gr[grep(snv.gr$HAP_SAMPLES, pattern = '^h2$')]
    snv.gr2 <- snv.gr.H2[,0]
    snv.gr2$H1 <- snv.gr.H2$REF
    snv.gr2$H2 <- snv.gr.H2$ALT
    snv.gr <- sort(c(snv.gr1, snv.gr2))
    snv.gr$sampleName = sampleName
  } else if (family == 'parent') {
    ## For parents require at least one to HOM
    snv.gr.H1 <- snv.gr[grep(snv.gr$HAP_SAMPLES, pattern = "^h1$")]
    snv.gr1 <- snv.gr.H1[,0]
    snv.gr1$H1 <- snv.gr.H1$ALT
    snv.gr1$H2 <- snv.gr.H1$REF
    snv.gr.H2 <- snv.gr[grep(snv.gr$HAP_SAMPLES, pattern = '^h2$')]
    snv.gr2 <- snv.gr.H2[,0]
    snv.gr2$H1 <- snv.gr.H2$REF
    snv.gr2$H2 <- snv.gr.H2$ALT
    snv.gr.H1H2 <- snv.gr[grep(snv.gr$HAP_SAMPLES, pattern = '^h1,h2$')]
    snv.gr3 <- snv.gr.H1H2[,0]
    snv.gr3$H1 <- snv.gr.H1H2$ALT
    snv.gr3$H2 <- snv.gr.H1H2$ALT
    snv.gr <- sort(c(snv.gr1, snv.gr2, snv.gr3))
    snv.gr$sampleName = sampleName
  } else {
    message("Please set 'family' parameter to either 'child' or 'parent'.")
  }
  return(snv.gr)
}

#' Find meiotic recombination breakpoints using circular binary segmentation (CBS)
#' 
#' This function takes as an input a \code{\link{VRanges-class}} object with extra metacolumn that contains
#' for each variant parantal identity. P1 -> parent1 & P2 -> parent2.
#' 
#' @param comparison A \code{\link{VRanges-class}} object 
#' @param minSeg Minimal length (number of variants) being reported as haplotype block (\code{fastseg} parameter).
#' @param smooth Number of consecutive variants being considered as a random error and so being corrected (flipped).
#' @param collapse.amb Set to \code{TRUE} if segments with ambiguous haplotype assignments should be collapsed.
#' @importFrom fastseg fastseg
#' @importFrom dplyr recode
#' @author David Porubsky
#' 
findRecomb <- function(comparison=NULL, minSeg=100, smooth=3, collapse.amb=TRUE) {
  
  ## Helper function
  switchValue <- function(x) {
    if (x == 1) {
      x <- 0
    } else {
      x <- 1  
    } 
  }
  
  collapseBins <- function(gr, id.field=0) {
    ind.last <- cumsum(runLength(Rle(mcols(gr)[,id.field]))) ##get indices of last range in a consecutive(RLE) run of the same value
    ind.first <- c(1,cumsum(runLength(Rle(mcols(gr)[,id.field]))) + 1) ##get indices of first range in a consecutive(RLE) run of the same value
    ind.first <- ind.first[-length(ind.first)]  ##erase last index from first range indices 
    collapsed.gr <- GRanges(seqnames=seqnames(gr[ind.first]), ranges=IRanges(start=start(gr[ind.first]), end=end(gr[ind.last])), mcols=mcols(gr[ind.first]))
    names(mcols(collapsed.gr)) <- names(mcols(gr[ind.first]))
    return(collapsed.gr)
  }
  
  ## Remove missing values
  comparison.filt <- comparison[mcols(comparison)[,1] != 'N']
  
  if (length(comparison.filt) >= 2*minSeg) {
    ## Recode comparison into in 0/1 vector
    comparison.filt$comp.vector <- dplyr::recode(mcols(comparison.filt)[,1], 'H1' = 0, 'H2' = 1)
    
    ## Run CBS segmentation on 0/1 vector 
    segs <- fastseg(comparison.filt$comp.vector, minSeg=minSeg)
    
    while (any(segs$num.mark <= smooth)) {
      toSwitch <- which(segs$num.mark <= smooth)
      switch.segs <- segs[toSwitch]
      switch.pos <- mapply(function(x,y) {x:y}, x=switch.segs$startRow, y=switch.segs$endRow)
      switch.pos <- unlist(switch.pos)
      
      switched.vals <- sapply(comparison.filt$comp.vector[switch.pos], switchValue) #SWITCH
      #comparison$comp.vector <- comparison$comp.vector[-switch.pos]  #DELETE
      comparison.filt$comp.vector[switch.pos] <- switched.vals
      segs <- fastseg(comparison.filt$comp.vector, minSeg=minSeg)
    }
    gen.ranges <- IRanges(start=start(comparison.filt)[segs$startRow], end=end(comparison.filt)[segs$endRow])
    ranges(segs) <- gen.ranges
    ## Add chromosome name
    seqlevels(segs) <- unique(as.character(seqnames(comparison.filt)))
    #segs$index <- index
    
    segs$match[segs$seg.mean <= 0.25] <- 'hap1'
    segs$match[segs$seg.mean >= 0.75] <- 'hap2'
    segs$match[segs$seg.mean > 0.25 & segs$seg.mean < 0.75] <- 'amb'
    
    ## Remove segments with mixed H1 and H2 signal
    if (collapse.amb) {
      segs <- segs[segs$match != 'amb']
    }
    
  } else {
    message("    Low density of informative SNVs, skipping ...")
    segs <- GenomicRanges::GRanges(seqnames=seqlevels(comparison), 
                                   ranges=IRanges(start=min(start(comparison)), end=max(end(comparison))),
                                   ID=NA, num.mark=0, seg.mean=0, startRow=0, endRow=0, match=NA
    )
  }
  
  ## Get haplotype segments
  if (length(segs) > 1) {
    segm <- collapseBins(gr = segs, id.field = 6)
  } else {
    segm <- segs
  }
  
  ## Meiotic recombination breakpoints
  if (length(segm) > 1) {
    suppressWarnings( recomb.break <- gaps(segm, start = start(segm)) )
    start(recomb.break) <- start(recomb.break) - 1
    end(recomb.break) <- end(recomb.break) + 1
  } else {
    ## Create a dummy breakpoint
    recomb.break <- GenomicRanges::GRanges(seqnames = seqlevels(comparison), ranges = IRanges(start=1, end=1))
  } 
  return(list(hap.segm = segm, recomb.break = recomb.break))
}