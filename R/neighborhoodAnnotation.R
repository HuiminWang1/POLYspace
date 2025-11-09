#' Annotate neighborhoods
#'
#' Adds annotations to neighborhoods in a \code{POLYspace} object.
#'
#' @param POLYspace A \code{POLYspace} object.
#' @param region Logical; include region-level information.
#' @param delimiter_node Delimiter between node labels (default \code{"--"}).
#' @param delimiter_region_celltype Delimiter between region and cell type (default \code{"_"}).
#' @param ncores Number of CPU cores to use.
#' @param mc.preschedule Logical; passed to \code{parallel::mclapply()}.
#'
#' @return An updated \code{POLYspace} object with neighborhood annotations.
#'
#' @export
neighborhoodAnnotation = function(POLYspace,
                                  region = TRUE,
                                  delimiter_node = '--',
                                  delimiter_region_celltype = '_',
                                  ncores = 1,
                                  mc.preschedule = TRUE){


  samples = POLYspace@samples
  targets = POLYspace@targets
  neighborhoods = POLYspace@neighborhoods
  id_of_samples = POLYspace@parameters$id_of_samples

  annotations = vector("list", length(targets))

  samples = lapply(id_of_samples,function(id){

    sa = samples[[id]]
    if(!('region' %in% colnames(sa)) | !region){
      sa[,'label'] = sa$celltype
      cat('Degenerate into global analysis in sample',id,'\n')
    }else{
      sa[,'label'] = paste(sa$region,sa$celltype,sep=delimiter_region_celltype)
    }
    return(sa)

  })
  names(samples) = id_of_samples


  neighborhood_size = rowSums(sapply(neighborhoods,lengths))
  agents = lapply(1:length(targets),function(x){

    if (!(targets[x] %in% c('singlet','pair')) & neighborhood_size[x]>0){
      ta = targets[[x]]
      ag = generateIsomorphismAgent(ta,
                                    delimiter_node = delimiter_node)
      return(ag)
    }

  })


  if('singlet' %in% targets){
    aot = lapply(id_of_samples,function(id) samples[[id]]$label)
    names(aot) = id_of_samples
    annotations[[which(targets == 'singlet')]] = aot
  }


  if(length(setdiff(unlist(targets),c('singlet')))>0){


    other_target_positions = lapply(neighborhoods, function(x) {
      x[(targets %in% c('singlet'))] = list(NULL)
      which(lengths(x) != 0)
    })
    multi_core_index = data.frame(
      sample_index = rep(id_of_samples, times = lengths(other_target_positions)),
      target_index = unlist(other_target_positions, use.names = FALSE)
    )

    annotations_of_non_singlet = parallel::mclapply(1:nrow(multi_core_index),function(i){

      sample_index = multi_core_index[i,'sample_index']
      target_index = multi_core_index[i,'target_index']

      sa = samples[[sample_index]]
      label = sa$label

      subs = neighborhoods[[sample_index]][[target_index]]
      subs = data.table::as.data.table(do.call('rbind',subs))
      subs = subs[, lapply(.SD, function(e) label[e])]

      ta = targets[[target_index]]
      ag = agents[[target_index]]

      canonical_annotation = isomorphismMapping(subs = subs,
                                                ta = ta,
                                                ag = ag)


      return(canonical_annotation)

    },mc.cores = ncores,mc.preschedule = mc.preschedule)

  }

  for(i in unique(multi_core_index$target_index)){
    one_target_index = (multi_core_index$target_index == i)
    annotations[[i]] = annotations_of_non_singlet[one_target_index]
    names(annotations[[i]]) = multi_core_index$sample_index[one_target_index]
  }

  annotations = lapply(id_of_samples, function(id) {

    sapply(annotations, function(x) x[[id]])

  })
  names(annotations) = id_of_samples

  POLYspace@annotations = annotations
  POLYspace@parameters$agents = agents

  return(POLYspace)

}
