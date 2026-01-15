#' Construct spatial networks
#'
#' Builds the spatial networks for a \code{POLYspace} object.
#'
#' @param POLYspace A \code{POLYspace} object.
#'
#' @param filtering Logical; if \code{TRUE}, apply edge filtering to remove edges
#'   longer than the specified \code{cutoff}.
#'
#' @param cutoff Either \code{"default"} (per-sample mean edge length + 2 × SD),
#'   or a named numeric vector of per-sample thresholds (with names corresponding
#'   to sample IDs). This parameter is used only when \code{filtering = TRUE}.
#'
#' @param ncores Integer; number of CPU cores used for parallel computation.
#'   Samples are distributed across cores in a sample-wise manner (one sample per core).
#'
#' @param mc.preschedule Logical; passed to \code{parallel::mclapply()} to control
#'   task scheduling behavior.
#'
#' @return An updated \code{POLYspace} object.
#'
#' @export
constructSpatialNetwork = function(POLYspace,
                                   filtering = TRUE,
                                   cutoff = 'default',
                                   ncores = 1,
                                   mc.preschedule = TRUE){


  samples= POLYspace@samples
  id_of_samples = POLYspace@parameters$id_of_samples
  S = POLYspace@parameters$number_of_samples

  nets = parallel::mclapply(1:S,function(s){

    data = samples[[s]]
    dims  = intersect(c("x", "y", "z"), names(data))
    coords = as.matrix(data[, ..dims])

    tri = RCDT::delaunay(coords)
    e   = data.table::as.data.table(tri$edges[,1:2])
    data.table::setnames(e, c("u1", "u2"))

    for (d in dims) {
      e[,paste0(d, "1")] = coords[e$u1, d]
      e[,paste0(d, "2")] = coords[e$u2, d]
    }

    diffs = coords[e$u1, , drop = FALSE] - coords[e$u2, , drop = FALSE]
    e[, di := sqrt(rowSums(diffs^2))]

    return(e)

  },mc.cores = ncores,mc.preschedule = mc.preschedule)

  if(filtering){

    if(!(identical('default',cutoff))){
      cutoff = cutoff[id_of_samples]
    }else{
      cutoff = lapply(1:S,function(s){
        e = nets[[s]]
        return(mean(e$di)+2*sd(e$di))
      })
      cutoff = unlist(cutoff)
    }

    nets = lapply(1:S,function(s){
      e = nets[[s]]
      e = e[e$di<=cutoff[s],]
      return(e)
    })

  }else{

    cutoff = rep(NA, times = S)

  }

  names(nets) = id_of_samples
  names(cutoff) = id_of_samples
  POLYspace@nets = nets
  POLYspace@parameters$cutoff_of_networks = cutoff

  return(POLYspace)

}
