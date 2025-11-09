#' Mine neighborhood
#'
#' Finds specified neighborhoods in a \code{POLYspace} object.
#'
#' @param POLYspace A \code{POLYspace} object.
#' @param targets A list of targets, e.g., \code{"singlet"}, \code{"pair"},
#'   or an adjacency matrix (numeric/logical) for custom targets.
#' @param ncores Number of CPU cores to use.
#' @param mc.preschedule Logical; passed to \code{parallel::mclapply()}.
#'
#' @return An updated \code{POLYspace} object whose \code{neighborhoods} slot includes
#'   the discovered neighborhoods.
#'
#' @export
neighborhoodMining = function(POLYspace,
                              targets = list('singlet',
                                             'pair',
                                             matrix(c(0,1,0,1,0,1,0,1,0),3,3)),
                              ncores = 1,
                              mc.preschedule = TRUE){


  nets = POLYspace@nets
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

    multi_core_index = base::expand.grid(id_of_samples,id_of_targets,stringsAsFactors = FALSE)
    colnames(multi_core_index) = c('sample_index','target_index')

    target_subs = parallel::mclapply(1:nrow(multi_core_index),function(i){

      sample_index = multi_core_index[i,'sample_index']
      target_index = multi_core_index[i,'target_index']

      ne = nets[[sample_index]]
      ta = targets[[target_index]]

      g = igraph::graph_from_data_frame(ne[,c('u1','u2')], directed = FALSE)
      gs = igraph::graph_from_adjacency_matrix(ta,mode='undirected')
      full_neighborhood = igraph::subgraph_isomorphisms(gs,
                                                        g,
                                                        method = 'lad',
                                                        induced = TRUE)
      ##
      if(length(full_neighborhood)>0){

        full_neighborhood = lapply(full_neighborhood,names)
        full_neighborhood = lapply(full_neighborhood,as.numeric)

        unique_neighborhood = findUniqueComposition(full_neighborhood)

        return(unique_neighborhood)
      }



    },mc.cores = ncores,mc.preschedule = mc.preschedule)

    for(i in id_of_targets){

      one_target_index = (multi_core_index$target_index == i)
      neighborhoods[[i]] = target_subs[one_target_index]
      names(neighborhoods[[i]]) = multi_core_index$sample_index[one_target_index]

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




