

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
