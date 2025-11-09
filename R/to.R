#' Convert segments to corresponding compositions
#'
#' Collapses segment labels into a composition string.
#'
#' @param segments Segment labels.
#' @param delimiter_node Delimiter used to join labels (default `"--"`).
#'
#' @return A string of corresponding compositions.
#'
#' @export
toComposition = function(segments,
                         delimiter_node = '--'){

  value = rn = NULL

  fun = do.call('rbind',segments)
  fun = data.table::as.data.table(fun)
  fun[, rn := .I]
  fun = data.table::melt(fun, id.vars = "rn")
  data.table::setorder(fun, rn, value)
  fun = fun[, .(composition = paste(value, collapse = delimiter_node)), by = rn]
  data.table::setorder(fun, rn)

  return(fun$composition)


}

#' Convert segments to cell-type neighborhoods
#'
#' Builds cell-type–level neighborhood labels from segment labels.
#'
#' @param segments Segment labels.
#' @param delimiter_node Delimiter between node labels (default "--").
#' @param delimiter_region_celltype Delimiter between region and cell type (default "_").
#' @param ta Target specification used to define neighborhood shape.
#' @param ag Isomorphism agent to canonicalize label order.
#'
#' @return A string of cell-type neighborhood labels.
#'
#' @export
toCelltypeNeighborhood = function(segments,
                                  delimiter_node = '--',
                                  delimiter_region_celltype = '_',
                                  ta,
                                  ag){


  cn = lapply(segments, function(p)
    sub(paste0(".*",delimiter_region_celltype), "", p))
  cn = data.table::as.data.table(do.call('rbind',cn))
  cn = isomorphismMapping(subs = cn,
                          ta = ta,
                          ag = ag,
                          delimiter_node = delimiter_node)
  return(cn)


}







