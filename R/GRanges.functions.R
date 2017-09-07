

#' Get boundaries of GRanges.
#' 
#' This function returns a set of unique start and end points of the
#' input GRange object.
#'
#' @param gr \code{\link[GenomicRanges]{GRanges}} object.
#'
getBoundaries <- function(gr){
  
  startGR <- GRanges(seqnames(gr), IRanges(start(gr), start(gr)), seqinfo=seqinfo(gr), mcols=mcols(gr))
  endGR <- GRanges(seqnames(gr), IRanges(end(gr), end(gr)), seqinfo=seqinfo(gr), mcols=mcols(gr))
  
  boundaryGR <- c(startGR, endGR)
  
  boundaryGR <- unique(boundaryGR)
  
  return(boundaryGR)
}


#' Finding overlapping genomic ranges by reciprocal overlap cutoff.
#'
#' Given two ranges of length l1 and l2 with overlapping length ol, the
#' reciprocal overlap is defined as lo / max(l1, l2).
#'
#' @param query,subject A GRanges or GRangesList object. See findOverlaps() in
#'   GenomicRanges package.
#' @param fraction The minimal fraction of reciprocal overlap allowed between
#'   ranges.
#' @param ... Futher arguments passed to findOverlaps() function.
#' @value An \code{\link[S4Vectors]{Hits}} object with indices of overlappin
#'   ranges.
reciprocalOverlaps <- function(query, subject, fraction = 0.5, ...){
  
  # get any overlap hit object
  hits <- findOverlaps(query, subject, type = "any", ...)
  
  # compute width of the overlapping region
  overlaps <- width(pintersect(query[queryHits(hits)], subject[subjectHits(hits)]))
  
  # get max length of query and subject for each overlap
  maxLen <- apply(
    cbind(width(query[queryHits(hits)]), width(subject[subjectHits(hits)])),
    1,max)
  
  # compute fraction of overlap
  fracOverlap <- overlaps / maxLen
  
  # filter for overlap >= threshold
  return(hits[fracOverlap >= fraction])
  
}
