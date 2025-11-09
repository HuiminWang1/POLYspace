#' @keywords internal
#' @import data.table
#' @importFrom methods new
#' @importFrom stats sd p.adjust na.omit setNames
#' @importFrom utils combn
"_PACKAGE"
utils::globalVariables(c(
  ".",
  "..dims",
  ".I", ".N",
  "di", "variable", "rank_identifier",
  "x1","x2","y1","y2"
))
