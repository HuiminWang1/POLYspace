#' Identify neighborhoods with unique compositions
#'
#' Drops duplicates by composition and returns the unique neighborhoods.
#'
#' @param neighborhood A list of neighborhoods.
#' @param delimiter_node A string used to join nodes. Default is "--".
#'
#' @return A list of neighborhoods with unique compositions.
#'
#' @export
findUniqueComposition = function(neighborhood,
                                 delimiter_node = '--'){

  value = rn = urn = composition = NULL

  fun = do.call('rbind',neighborhood)
  fun = data.table::as.data.table(fun)
  fun[, rn := .I]
  fun = data.table::melt(fun, id.vars = "rn")
  data.table::setorder(fun, rn, value)
  fun = fun[, .(composition = paste(value, collapse = delimiter_node)), by = rn]
  fun = fun[, .(urn = data.table::first(rn)), by = composition]

  return(neighborhood[fun$urn])

}
