#' Differential neighborhoods across samples
#'
#' Performs differential enrichment analysis of neighborhoods across samples.
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
#'   
#' @param des Design matrix specifying the experimental design.
#'
#' @param coef Coefficient name or index in \code{des} to be tested.
#'
#' @param p_adjust_method Character; method for p-value adjustment
#'   (e.g., \code{"fdr"}).
#'
#' @param return_object Logical; if \code{TRUE}, return the updated
#'   \code{POLYspace} object instead of a result table.

#' @return If \code{return_object = TRUE}, an updated \code{POLYspace} object with
#'   the differential result listattached; otherwise, a list of differential results.
#'
#'   Within the list, it contains the following components:
#'   \itemize{
#'     \item \code{results}: differential neighborhood results, with ordering
#'     consistent with the \code{targets} list;
#'     \item \code{targets}: the targets used in the analysis;
#'     \item \code{des}: the design matrix used for differential testing;
#'     \item \code{coef}: the coefficient tested in the design matrix.
#'   }
#'
#' @export
polyEnrichment = function(object,
                          region = TRUE,
                          targets = 'default',
                          neighborhoods_of_interest = c('union','intersect'),
                          des,
                          coef,
                          p_adjust_method = 'fdr',
                          return_object = FALSE){


  matrices_targets = retrieveFeatureMatrix(object = object,
                                           region = region,
                                           targets = targets,
                                           neighborhoods_of_interest = neighborhoods_of_interest)
  matrices = matrices_targets$matrices
  neighborhoods_info = matrices_targets$neighborhoods_info

  if(identical('default',targets)){

    targets = matrices_targets$targets

  }


  results = lapply(1:length(targets), function(ta){

    testdt = matrices[[ta]]
    if(!is.null(testdt)){

      testdt = testdt[,rownames(des)]
      testdt = testdt[rowSums(is.na(testdt)) != ncol(testdt),,drop=FALSE]
      if(nrow(testdt)>0){
        res = limma::topTable(limma::eBayes(limma::lmFit(testdt,design=des)),
                              n=nrow(testdt),
                              coef=coef,
                              adjust.method=p_adjust_method)
        res[,'target'] = ta
        res = na.omit(res)
        testdt_info = neighborhoods_info[[ta]]
        if(region){

          res[,c('celltype','region','composition')] =
            testdt_info[match(rownames(res),testdt_info$neighborhood),c('celltype','region','composition')]

        }else{
          res[,c('composition')] =
            testdt_info[match(rownames(res),testdt_info$neighborhood),c('composition')]

        }

        return(res)
      }

    }


  })
  if(return_object){

    POLYspace@results$multi = list(results = results,
                                   targets = targets,
                                   des = des,
                                   coef = coef)
    return(POLYspace)

  }else{
    out = list(results = results,
               targets = targets,
               des = des,
               coef = coef)

    return(out)

  }

}


