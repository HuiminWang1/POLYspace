#' Retrieve single-sample enrichment results
#'
#' Extracts per-sample enrichment results from a \code{POLYspace} object or from
#' a list of \code{POLYspace} objects.
#'
#' @param object A \code{POLYspace} object or a list of \code{POLYspace} objects.
#' @param targets Targets to include.
#' @param neighborhoods_of_interest Character vector (e.g., \code{c("union","intersect")}),
#'   or a user-specified set of neighborhoods to include.
#'
#' @return A list with components \code{results} (data.table of
#'   single-sample enrichment) and \code{targets} (the targets used).
#'
#' @export
retrieveSingleResult = function(object,
                                targets = 'default',
                                neighborhoods_of_interest = c('union','intersect')){


  if(is.list(object)){

    id_of_samples = sapply(object,function(x) x@parameters$id_of_samples)

    if(identical('default',targets)){
      targets = lapply(object,function(x) x@targets)
      targets = do.call('c',targets)
      targets = unique(targets)
    }

    if(!(identical(names(object),id_of_samples))){
      names(object) = id_of_samples
    }

    sr = lapply(id_of_samples, function(id) {
      re = vector("list", length(targets))
      x = object[[id]]
      re[match(x@targets,targets)] = x@results$enrichment[[id]]
      return(re)
    })
    names(sr) = id_of_samples

  }else{

    id_of_samples = object@parameters$id_of_samples
    targets = object@targets

    sr = object@results$enrichment

  }


  if(identical('union',neighborhoods_of_interest)){

    neighborhoods_of_interest = lapply(seq_along(targets), function(ta){

      noi = lapply(sr, function(x){
        x[[ta]][['neighborhood']]
      })
      noi = do.call('c',noi)
      noi = sort(unique(noi))
      return(noi)

    })

  }

  if(identical('intersect',neighborhoods_of_interest)){

    neighborhoods_of_interest = lapply(seq_along(targets), function(ta){

      noi = lapply(sr, function(x){
        x[[ta]][['neighborhood']]
      })
      noi = Reduce(intersect,noi)
      noi = sort(unique(noi))
      return(noi)

    })

  }



  results = lapply(seq_along(targets), function(ta){
    noi = neighborhoods_of_interest[[ta]]

    ma = lapply(sr, function(x){
      x[[ta]]
    })
    ma = data.table::as.data.table(do.call('rbind',ma))
    ma = ma[ma$neighborhood %in% noi]
    return(ma)
  })

  return(list(results = results,
              targets = targets))


}


