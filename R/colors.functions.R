


#' Add an alpha value to a colour
#' 
#' From: http://www.magesblog.com/2013/04/how-to-change-alpha-value-of-colours-in.html
#' 
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

#' Mixes two group of colors.
#' 
#' Adds transparancy of 50% to both groups of colors and than overlays them
#' pairwise.
#' 
#' @param col1 vector of colors
#' @param col2 vector of colors
#' @param ratio numeric between 0 and 1 as fraction of first color. Second color
#'   will be \code{1-ratio}.
mixColors <- function(col1, col2, ratio=.5){
  
  rgb1 <- col2rgb(col1, alpha=TRUE)
  rgb2 <- col2rgb(col2, alpha=TRUE)
  
  cols <- c()
  
  for(i in 1:length(col1)){
    
    for (j in 1:length(col2)){
      
      c1 <- rgb1[,i]
      c2 <- rgb2[,j]
      
      m <- ratio*c1 + (1-ratio)*c2
      
      mixedCol <- rgb(m[1], m[2], m[3], alpha=m[4], max=255)
      
      cols <- c(cols, mixedCol)
      
    }
  }
  
  return(cols)
}
