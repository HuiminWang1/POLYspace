#' POLYspace S4 class
#'
#' A container for \code{POLYspace} analyses, storing inputs,
#' intermediate objects, and analysis results.
#'
#' @slot samples \code{list}; per-sample tables. Each table typically contains
#'   columns such as \code{id}, \code{x}, \code{y} (and optional \code{z} for 3D samples),
#'   \code{celltype}, and \code{region} (which may be omitted if a global analysis is performed).
#'
#' @slot nets \code{list}; spatial adjacency networks for each sample.
#'
#' @slot targets \code{list}; specifications of neighborhood topologies used as targets
#'   for neighborhood mining and enrichment analyses.
#'
#' @slot neighborhoods \code{list}; discovered neighborhoods per sample and target.
#'   Each element of the list corresponds to a sample; within each sample,
#'   each sub-list element corresponds to a target, with ordering consistent
#'   with the \code{targets} slot.
#'
#' @slot annotations \code{list}; canonical labels assigned to discovered neighborhoods.
#'   Each element of the list corresponds to a sample; within each sample,
#'   each sub-list element corresponds to a target, with ordering consistent
#'   with the \code{targets} slot.
#'
#' @slot results \code{list}; analysis outputs, including per-sample enrichment results
#'   and group-level differential results.
#'
#' @slot parameters \code{list}; run-time parameters and analysis options.
#'
#' @exportClass POLYspace
setClass('POLYspace',slots = list(
  samples = 'list',
  nets = 'list',
  targets = 'list',
  neighborhoods = 'list',
  annotations = 'list',
  results = 'list',
  parameters = 'list'
))

#' Show a welcome message
#'
#' Prints a brief welcome for a POLYspace object.
#'
#' @param POLYspace A \code{POLYspace} object.
#'
#' @return No return value, called for side effects.
#'
#' @export
showWelcomeMessage = function(POLYspace)
{
  cat("Welcome to POLYspace!\n")
  cat("************************************\n")
  cat("  INPUT INFO:\n")
  cat("    - Number of samples:", POLYspace@parameters$number_of_samples, "\n")
  cat("    - ID of samples:", POLYspace@parameters$id_of_samples, "\n")
  cat("    - Number of cells:", POLYspace@parameters$number_of_cells, "\n")
  cat("************************************\n")
}

#' Create a POLYspace object
#'
#' Constructs a \code{POLYspace} container from input samples.
#'
#' @param samples A list of per-sample tables. Each table typically contains
#'   columns such as \code{id}, \code{x}, \code{y} (and optional \code{z}),
#'   \code{celltype}, and optionally \code{region} if region information is available.
#'
#' @param id_of_samples Optional character vector of sample IDs. If \code{NULL},
#'   sample IDs are automatically generated as \code{"1"}, \code{"2"}, \code{"3"}, \dots
#' @return A \code{POLYspace} object.
#'
#' @export
createPOLYspaceObject = function(
    samples,
    id_of_samples = NULL
)
{

  S = length(samples)
  stopifnot(is.list(samples))

  if(is.null(id_of_samples)){
    id_of_samples = as.character(seq_len(S))
  }

  for (s in seq_along(samples)) {
    if (!data.table::is.data.table(samples[[s]])) {
      data.table::setDT(samples[[s]])
      cat(sprintf("Sample data %s coerced to a data table\n", id_of_samples[s]))
    }
  }

  names(samples) = id_of_samples

  parameters = list(number_of_samples = S,
                    id_of_samples = id_of_samples,
                    number_of_cells = sapply(samples, nrow))
  POLYspace = new(Class = "POLYspace",
                  samples = samples,
                  nets = list(),
                  targets = list(),
                  neighborhoods = list(),
                  annotations = list(),
                  results = list(),
                  parameters = parameters)

  showWelcomeMessage(POLYspace)

  return(POLYspace)


}



