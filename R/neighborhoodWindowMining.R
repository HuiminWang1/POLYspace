#' Auto-select a window size from point density
#'
#' Chooses a square window side length so that each window contains, on average,
#' about `window_mean_points` points.
#'
#' @param x Numeric vector of x-coordinates.
#' @param y Numeric vector of y-coordinates.
#' @param window_mean_points Target mean number of points per window (default 10000).
#' @param epsilon Numeric tolerance when solving for the window size (default 1).
#'
#' @return A single numeric value: the window side length.
#'
#' @export
autoWindowSize = function(x,
                          y,
                          window_mean_points = 10000,
                          epsilon = 1) {

  area = (max(x) - min(x)) * (max(y) - min(y))
  density = length(x) / area
  window_size = sqrt(window_mean_points / density)
  window_size = window_size + epsilon

  return(window_size)

}

#' Generate sliding windows
#'
#' Creates sliding windows over a spatial area, optionally choosing
#' the window size automatically from point density.
#'
#' @param sa A spatial sample table with coordinates.
#' @param ne A spatial network for the sample.
#' @param id The sample ID.
#' @param ms Size of the target.
#' @param window_size Numeric side length, or `"auto"` to infer from density.
#' @param window_mean_points Target mean points per window when `window_size="auto"`.
#' @param epsilon Tolerance used when solving for the auto window size.
#'
#' @return A list of windows and the window size.
#'
#' @export
generateWindows = function(sa,
                           ne,
                           id,
                           ms,
                           window_size = 'auto',
                           window_mean_points = 10000,
                           epsilon = 1){


  if(identical("auto",window_size)){

    window_size = autoWindowSize(x = sa$x,
                                 y = sa$y,
                                 window_mean_points = window_mean_points,
                                 epsilon = epsilon)


  }


  # Step 1: compute nearest neighbor distances
  buffer = sum(sort(ne$di, decreasing = TRUE)[1:ms])
  step = window_size - (buffer+epsilon)
  if(step<=0){
    stop('Please try a larger window\n')
  }

  message(sprintf("ID = %s; auto window = %.2f; step = %.2f",
                  id, window_size,  step))

  # Step 2: compute ranges
  x_bounds = range(sa$x) + c(-epsilon, epsilon)
  y_bounds = range(sa$y) + c(-epsilon, epsilon)

  x_starts = seq(x_bounds[1], x_bounds[2], by = step)
  y_starts = seq(y_bounds[1], y_bounds[2], by = step)


  grid_coords = data.table::CJ(ix = x_starts, iy = y_starts)

  return(list(grid_coords = grid_coords,
              window_size = window_size))

}


#' Windowed neighborhood mining
#'
#' Mines specified neighborhood targets in sliding windows over a \code{POLYspace} object.
#'
#' @param POLYspace A \code{POLYspace} object.
#' @param targets List of targets to search.
#' @param window_size Numeric side length, or \code{"auto"} to infer from density.
#' @param window_mean_points Target mean points per window when \code{window_size = "auto"}.
#' @param epsilon Tolerance used for auto window size.
#' @param ncores Number of CPU cores.
#' @param mc.preschedule Logical; passed to \code{parallel::mclapply()}.
#'
#' @return An updated \code{POLYspace} object whose \code{neighborhoods} slot includes
#'   windowed target discoveries.
#'
#' @export
neighborhoodWindowMining = function(POLYspace,
                                    targets = list('singlet',
                                               'pair',
                                               matrix(c(0,1,0,1,0,1,0,1,0),3,3)),
                                    window_size = 'auto',
                                    window_mean_points = 10000,
                                    epsilon = 1,
                                    ncores = 1,
                                    mc.preschedule = TRUE){




  samples = POLYspace@samples
  nets= POLYspace@nets
  id_of_samples = POLYspace@parameters$id_of_samples

  neighborhoods = vector("list", length(targets))

  if('singlet' %in% targets){

    number_of_cells = POLYspace@parameters$number_of_cells
    subs = lapply(id_of_samples,function(id) as.list(1:number_of_cells[id]))
    names(subs) = id_of_samples
    neighborhoods[[which(targets == 'singlet')]] = subs

  }


  if('pair' %in% targets){

    subs = lapply(id_of_samples,function(id) {
      ne = nets[[id]]
      ps = Map(c, ne$u1, ne$u2)
      return(ps)
    })
    names(subs) = id_of_samples
    neighborhoods[[which(targets == 'pair')]] = subs

  }


  if(length(setdiff(unlist(targets),c('singlet','pair')))>0){

    id_of_targets = which(!(targets %in% c('singlet','pair')))
    ms = max(sapply(targets[id_of_targets],ncol))

    auto_window = identical('auto',window_size)
    if(auto_window){
      window_size = rep("auto", length(id_of_samples))
    }
    window_infos = lapply(id_of_samples,function(id){

      info = generateWindows(sa = samples[[id]],
                             ne = nets[[id]],
                             id = id,
                             ms = ms,
                             window_size = window_size[which(id_of_samples==id)],
                             window_mean_points = window_mean_points,
                             epsilon = epsilon)
      return(info)


    })
    names(window_infos) = id_of_samples
    window_counts = sapply(window_infos, function(x) nrow(x[['grid_coords']]))

    multi_core_index = do.call(rbind, lapply(id_of_targets, function(i) {

      do.call(rbind, lapply(id_of_samples, function(j) {
        data.frame(
          sample_index = j,
          target_index = i,
          window_index = 1:window_counts[[j]]
        )
      }))

    }))


    window_target_subs = parallel::mclapply(1:nrow(multi_core_index),function(i){

      sample_index = multi_core_index[i,'sample_index']
      target_index = multi_core_index[i,'target_index']
      window_index = multi_core_index[i,'window_index']

      ne = nets[[sample_index]]
      ta = targets[[target_index]]

      ix = window_infos[[sample_index]]$grid_coords$ix[window_index]
      iy = window_infos[[sample_index]]$grid_coords$iy[window_index]
      window_size = window_infos[[sample_index]]$window_size

      x_rng = ix + c(0, window_size)
      y_rng = iy + c(0, window_size)



      pts_center = ne[
        data.table::between(x1, x_rng[1], x_rng[2]) &
          data.table::between(x2, x_rng[1], x_rng[2]) &
          data.table::between(y1, y_rng[1], y_rng[2]) &
          data.table::between(y2, y_rng[1], y_rng[2])
      ]
      g = igraph::graph_from_data_frame(pts_center[,c('u1','u2')], directed = FALSE)
      gs = igraph::graph_from_adjacency_matrix(ta,mode='undirected')
      full_neighborhood = igraph::subgraph_isomorphisms(gs,
                                                        g,
                                                        method = 'lad',
                                                        induced = TRUE)
      if(length(full_neighborhood)>0){
        ##
        full_neighborhood = lapply(full_neighborhood,names)
        full_neighborhood = lapply(full_neighborhood,as.numeric)
        unique_neighborhood = findUniqueComposition(full_neighborhood)

        return(unique_neighborhood)

      }



    },mc.cores = ncores,mc.preschedule = mc.preschedule)


    multi_core_index2 = base::expand.grid(id_of_samples,id_of_targets,stringsAsFactors = FALSE)
    colnames(multi_core_index2) = c('sample_index','target_index')

    multi_core_index[,'comb'] = paste(multi_core_index$sample_index,multi_core_index$target_index,sep='_')
    multi_core_index2[,'comb'] = paste(multi_core_index2$sample_index,multi_core_index2$target_index,sep='_')

    target_subs = parallel::mclapply(multi_core_index2$comb,function(i){


      full_neighborhood = do.call('c',window_target_subs[multi_core_index$comb == i])
      if(length(full_neighborhood)>0){

        unique_neighborhood = findUniqueComposition(full_neighborhood)
        return(unique_neighborhood)

      }


    },mc.cores = ncores,mc.preschedule = mc.preschedule)

    for(i in id_of_targets){

      one_target_index = (multi_core_index2$target_index == i)
      neighborhoods[[i]] = target_subs[one_target_index]
      names(neighborhoods[[i]]) = multi_core_index2$sample_index[one_target_index]

    }


  }

  neighborhoods = lapply(id_of_samples, function(id) {

    sapply(neighborhoods, function(x) x[[id]])

  })
  names(neighborhoods) = id_of_samples

  POLYspace@neighborhoods = neighborhoods
  POLYspace@targets = targets

  return(POLYspace)

}
