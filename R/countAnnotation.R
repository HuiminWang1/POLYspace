#' Count annotations for neighborhoods
#'
#' Summarizes how often each annotation label appears for the given
#' neighborhoods of interest. Counts are transformed by \code{log2} (with a
#' pseudocount of 1).
#'
#' @param neighborhoods_of_interest A vector identifying neighborhoods of interest.
#' @param annotation A vector of labels.
#'
#' @return A named numeric vector of \code{log2}-transformed counts per
#'   neighborhood of interest (pseudocount = 1).
#'
#' @export
countAnnotation = function(neighborhoods_of_interest,
                           annotation){

  V1 = NULL

  frame = rep(0,length(neighborhoods_of_interest))
  names(frame) = neighborhoods_of_interest

  if(length(annotation)>0){
    anno = data.table::data.table(V1 = annotation)
    anno = anno[anno$V1 %in% names(frame),,drop=FALSE]
    ans = anno[,.N,by = V1]
    N = ans$N
    N = log2(N+1)
    names(N) = ans$V1

    frame[names(N)] = as.numeric(N)
  }

  return(frame)

}

