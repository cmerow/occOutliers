
#' omit outlying pres
#' @param xy data.frame with 2 columns
#' @param percent numeric on [0,100]
#' @param unin character; input units. option are "m" and "km"
#' @param unout character; input units. option are "ha","km2", "m2"
# @param ... additional functions to be passed to \code{\link[raster]{writeRaster}}
# @description  Divide raster by the sum of all cells.
# @export
# unin = c("m");  unout ='m2'
.presPercentile=function (xy, percent = 95, unin = c("m", "km"), unout = c("ha","km2", "m2")) {
  # for testing
  # xy=pres; percent = NULL; unin = c("m", "km"); unout = c("ha","km2", "m2")
  #if (!inherits(xy, "SpatialPoints")) # not needed; i already check it
    #stop("xy should be of class SpatialPoints")
  if (ncol(sp::coordinates(xy)) > 2)
    stop("xy should be defined in two dimensions")
  pfs <- sp::proj4string(xy)
  if(!is.null(percent)){
    if (length(percent) > 1)
      stop("only one value is required for percent")
    if (percent > 100) {
      warning("Using all relocations (percent>100)")
      percent <- 100
    }
  }
  unin <- match.arg(unin)
  unout <- match.arg(unout)
  if (inherits(xy, "SpatialPointsDataFrame")) {
    if (ncol(xy) != 1) {
      #warning("xy should contain only one column (the id of the animals), id ignored")
      id <- factor(rep("a", nrow(as.data.frame(xy))))
    } else {
      id <- xy[[1]]
    }
  } else {
    id <- factor(rep("a", nrow(as.data.frame(xy))))
  }
  
  if (min(table(id)) < 5) stop("At least 5 relocations are required to fit an home range")
  id <- factor(id)
  xy <- as.data.frame(sp::coordinates(xy))
  r <- split(xy, id)
  est.cdg <- function(xy) apply(xy, 2, mean)
  cdg <- lapply(r, est.cdg)
  levid <- levels(id)
  res =lapply(1:length(r), function(i) {
    k <- levid[i]
    df.t <- r[[levid[i]]]
    cdg.t <- cdg[[levid[i]]]
    dist.cdg <- function(xyt) {
      d <- sqrt(((xyt[1] - cdg.t[1])^2) + ((xyt[2] - cdg.t[2])^2))
      return(d)
    }
    di <- apply(df.t, 1, dist.cdg)
    key <- c(1:length(di))
    if(!is.null(percent)){
      acons <- key[di <= stats::quantile(di, percent/100)]
    } else { acons=key }
    xy.t <- df.t[acons, ]
    sp::coordinates(xy.t)=c(1,2)
    return(list(xy.t=xy.t,dist.from.centroid=di))
  })
  res
}

#' find records mor than 1.5 times the interquartile range beyond the upper quartile 
#' @param dists a numeric vector
.iqrOutlier=function(dists){
  
  #q1=stats::quantile(dists, .25)
  q3=stats::quantile(dists, .25)
  iqr=stats::IQR(dists)
  which(dists > (q3 + 1.5*iqr))
  #subset(df, df$A> (Q1 - 1.5*IQR) & df$A< (Q3 + 1.5*IQR))
}

#' Vectorized version of grep
#' @description vectorized version of grep
#' @param pattern character string containing a regular expression (or character string for ‘fixed = TRUE’) to be matched in the given character vector.  Coerced by ‘as.character’ to a character string if possible.  If a character vector of length 2 or more is supplied, the first element is used with  a warning.
#' @param x a character vector where matches are sought, or an object which can be coerced by ‘as.character’ to a character vector.
# @export
.vgrep=function(pattern,x){mapply(function(y){grep(y,pattern)},x)}
