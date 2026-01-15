#' Retrieve single-sample enrichment results
#'
#' Extracts per-sample enrichment results from a \code{POLYspace} object or from
#' a list of \code{POLYspace} objects.
#'
#' @param object A \code{POLYspace} object or a list of \code{POLYspace} objects.
#' 
#' @param targets Targets to include. \code{"default"} uses the union of targets
#'   across all input objects.
#'   
#' @param neighborhoods_of_interest Neighborhoods to include in the feature matrix.
#'   This can be specified as:
#'   \itemize{
#'     \item \code{"union"}: for each neighborhood shape, use the union of
#'     neighborhoods observed across all samples;
#'     \item \code{"intersect"}: for each neighborhood shape, use only neighborhoods
#'     shared by all samples;
#'     \item a user-defined list: a list whose elements correspond to targets, where
#'     each element specifies the neighborhoods of interest for that target.
#'   }
#'   
#' @return A list with components:
#'   \itemize{
#'     \item \code{results}: a \code{data.table} containing single-sample
#'     enrichment results;
#'     \item \code{targets}: the targets used in the analysis.
#'   }
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


