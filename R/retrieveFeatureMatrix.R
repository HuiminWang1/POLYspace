#' Retrieve neighborhood feature matrices
#'
#' Retrieves neighborhood-level feature matrices and related information from
#' a \code{POLYspace} object or a list of \code{POLYspace} objects.
#'
#' @param object A \code{POLYspace} object or a list of \code{POLYspace} objects.
#'
#' @param region Logical; whether to use region-level neighborhoods when
#'   constructing the feature matrix, which determines the level of enrichment
#'   analysis to include.
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
#' @return A list containing:
#'   \itemize{
#'     \item \code{matrices}: a list of feature matrices, with ordering consistent
#'     with the targets in the output list. Rows correspond to neighborhoods and
#'     columns correspond to sample IDs;
#'     \item \code{neighborhoods_info}: metadata describing the neighborhoods
#'     included in the matrices, such as regional or cellular composition;
#'     \item \code{targets}: the targets used for feature matrix retrieval.
#'   }
#'
#' @export
retrieveFeatureMatrix = function(object,
                                 region = TRUE,
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

  matrices = vector("list", length(targets))
  neighborhoods_info = vector("list", length(targets))

  for(ta in seq_along(targets)){
    noi = neighborhoods_of_interest[[ta]]

    ma = lapply(sr, function(x){
      x[[ta]]

    })
    ma = data.table::as.data.table(do.call('rbind',ma))
    ma = ma[ma$neighborhood %in% noi]
    if(nrow(ma)>0){

      testdt = reshape2::acast(ma[,c('neighborhood','sample_id','enrichment')],
                               neighborhood~sample_id, value.var = 'enrichment')
      if(length(setdiff(id_of_samples,colnames(testdt)))>0){

        new_col_id = setdiff(id_of_samples,colnames(testdt))
        new_col = matrix(NA, nrow = nrow(testdt), ncol = length(new_col_id))
        colnames(new_col) = new_col_id

        testdt = cbind(testdt, new_col)

      }
      testdt = testdt[,id_of_samples]
      matrices[[ta]] = testdt

      if(region){
        testdt_info = ma[,c('neighborhood','celltype','region','composition')]
        testdt_info = testdt_info[!duplicated(testdt_info)]
        neighborhoods_info[[ta]] = testdt_info
      }else{

        testdt_info = ma[,c('neighborhood','composition')]
        testdt_info = testdt_info[!duplicated(testdt_info)]
        neighborhoods_info[[ta]] = testdt_info
      }

    }

  }


  return(list(matrices = matrices,
              neighborhoods_info = neighborhoods_info,
              targets = targets))

}



