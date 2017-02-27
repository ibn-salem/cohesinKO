#' A function to scrren TADs for overlap with GRBs
#' 
#' Code from Nathan Harmston, Dimitris Polychronopoulos, and Elizabeth Ing-Simmons
screen.grbs = function(tads, grbs, top=0.8, bottom=0.2, plot=FALSE){
  
  overlaps_all <- findOverlaps(tads, grbs)
  overlaps_tab <- table(queryHits(overlaps_all))
  tads$id <- 1:length(tads)
  
  tads$overlap_count <- ifelse(tads$id %in% names(overlaps_tab), 
                               overlaps_tab[as.character(tads$id)], 
                               0)
  
  #get total width and fraction overlap
  tads$overlap <- 0
  tads$fractionoverlap <- 0
  overlap_idx <- queryHits(overlaps_all)
  for (i in overlap_idx){
    tads$overlap[i] <- sum(width(BiocGenerics::intersect(tads[i], grbs)))
  }
  tads$fractionoverlap <- tads$overlap/width(tads)
  
  if(plot){
    plot(density(tads$fractionoverlap), type = "l", 
         main = "Fraction overlap between TADs and GRBs")
    abline(v = c(bottom, top), lty = 2)
  }
  
  tads$class <- "nonGRB"
  tads$class[tads$fractionoverlap > bottom] <- "screened"
  tads$class[tads$fractionoverlap > top & tads$overlap_count == 1] <- "GRB"
  tads$class[tads$overlap_count > 1] <- "screened"
  return(tads)
}
