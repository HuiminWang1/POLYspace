#' Enrichment analysis across samples
#'
#' Tests enrichment of neighborhoods across samples.
#'
#' @param object A \code{POLYspace} object or a list of \code{POLYspace} objects.
#' @param region Logical; stratify by region.
#' @param targets Targets to test.
#' @param neighborhoods_of_interest Neighborhoods of interest to include in the test.
#' @param des Design matrix.
#' @param coef Coefficient name or index in \code{des} to be tested.
#' @param p_adjust_method Method for p-value adjustment (e.g., \code{"fdr"}).
#' @param return_object Logical; if \code{TRUE}, return the updated object instead of a table.
#'
#' @return If \code{return_object = TRUE}, an updated \code{POLYspace} object with
#'   enrichment results attached; otherwise a list of enrichment results.
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


