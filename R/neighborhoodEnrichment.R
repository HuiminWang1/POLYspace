#' Neighborhood enrichment analysis
#'
#' Performs enrichment analysis of neighborhoods.
#'
#' @param POLYspace A \code{POLYspace} object.
#' @param region Logical; include region-level stratification.
#' @param delimiter_node Delimiter between node labels (default \code{"--"}).
#' @param delimiter_region_celltype Delimiter between region and cell type (default \code{"_"}).
#' @param neighborhoods_of_interest Character vector of neighborhoods of interest or \code{"default"}.
#' @param permutation_times Integer; number of permutations (default \code{100}).
#' @param p_adjust_method Method for p-value adjustment (e.g., \code{"fdr"}).
#' @param ncores Number of CPU cores.
#' @param mc.preschedule Logical; passed to \code{parallel::mclapply()}.
#' @param set_seed Logical; if \code{TRUE}, set the RNG seed to \code{seed}.
#' @param seed Integer RNG seed.
#'
#' @return An updated \code{POLYspace} object whose \code{results} slot includes
#'   tables of enrichment statistics.
#'
#' @export
neighborhoodEnrichment = function(POLYspace,
                                  region = TRUE,
                                  delimiter_node = '--',
                                  delimiter_region_celltype = '_',
                                  neighborhoods_of_interest = 'default',
                                  permutation_times = 100,
                                  p_adjust_method = 'fdr',
                                  ncores = 1,
                                  mc.preschedule = TRUE,
                                  set_seed = TRUE,
                                  seed = 123){

  samples = POLYspace@samples
  targets = POLYspace@targets
  neighborhoods = POLYspace@neighborhoods
  annotations = POLYspace@annotations
  id_of_samples = POLYspace@parameters$id_of_samples
  agents = POLYspace@parameters$agents

  samples = lapply(id_of_samples,function(x){

    sa = samples[[x]]
    if(!('region' %in% colnames(sa)) | !region){
      sa[,'label'] = sa$celltype
    }else{
      sa[,'label'] = paste(sa$region,sa$celltype,sep=delimiter_region_celltype)
    }
    return(sa)

  })
  names(samples) = id_of_samples


  results = lapply(id_of_samples,function(i){

    shuffled_samples = lapply(1:permutation_times, function(n){

      if(set_seed == TRUE) {
        set.seed(seed = seed+n)
      }

      sa = samples[[i]]
      sa$celltype = sample(sa$celltype)

      if(!('region' %in% colnames(sa)) | !region){
        sa[,'label'] = sa$celltype
      }else{
        sa[,'label'] = paste(sa$region,sa$celltype,sep=delimiter_region_celltype)
      }

      return(sa)

    })

    out = lapply(1:length(targets), function(j){

      ta = targets[[j]]
      ag = agents[[j]]
      nei = neighborhoods[[i]][[j]]
      ann = annotations[[i]][[j]]

      if (identical("full", neighborhoods_of_interest)) {
        
        noi = base::sort(unique(ann))
        
      } else if (identical("default", neighborhoods_of_interest)) {
        
        if (list(ta) %in% c("singlet", "pair")) {
          
          noi = base::sort(unique(ann))
          
        } else {
          
          noijd_frame = lapply(ann, function(x) {
            y = stringr::str_split(x, delimiter_node)[[1]]
            sapply(
              stringr::str_split(y, delimiter_region_celltype),
              function(z) z[2]
            )
          })
          
          jd = sapply(noijd_frame, function(x) {
            length(unique(x)) == length(x)
          })
          
          noi = ann[jd]
        }
        
      } else {
        
        noi = neighborhoods_of_interest[[j]]
        ann = ann[ann %in% noi]
        
      }
    
      if(length(noi)>0 & length(nei)>0){

        p.adjusted = enrichness = NULL

        ob = countAnnotation(neighborhoods_of_interest = noi,
                             annotation = ann)



        segments = strsplit(noi, delimiter_node)
        if(!identical('singlet',ta)){
          coi = toComposition(segments = segments,
                              delimiter_node = delimiter_node)
        }else{
          coi = noi
        }

        nu = parallel::mclapply(1:permutation_times, function(n){


          sa = shuffled_samples[[n]]
          label = sa$label

          subs = nei
          subs = data.table::as.data.table(do.call('rbind',subs))
          subs = subs[, lapply(.SD, function(e) label[e])]

          if(!identical('singlet',ta)){
            caan = isomorphismMapping(subs = subs,
                                      ta = ta,
                                      ag = ag,
                                      noi = noi,
                                      coi = coi)
          }else{

            caan = label[label %in% noi]

          }
          null_ob = countAnnotation(neighborhoods_of_interest = noi,
                                    annotation = caan)

          return(null_ob)

        },mc.cores = ncores,mc.preschedule = mc.preschedule)

        nu = do.call('cbind',nu)
        jd = apply(nu, 1, sd) ==0

        zscore = (ob-rowMeans(nu))/apply(nu, 1, sd)
        pvalue = rowMeans(abs(nu - rowMeans(nu)) >= abs(ob -rowMeans(nu)))
        zscore[jd] = pvalue[jd] = NA
        p.adjusted =  p.adjust(pvalue,method = p_adjust_method)


        noi_comp = sapply(segments, function(p) paste(unique(sort(sub(paste0(".*",delimiter_region_celltype), "", p))),collapse = delimiter_node))
        if(region){

          if(!identical('singlet',ta)){

            noi_no_region = toCelltypeNeighborhood(segments = segments,
                                                   delimiter_region_celltype = delimiter_region_celltype,
                                                   ta = ta,
                                                   ag = ag)

          }else{

            noi_no_region = sapply(segments, function(p)
              sub(paste0(".*", delimiter_region_celltype), "", p))

          }

          noi_region = sapply(segments, function(p) paste(unique(sort(sub(paste0(delimiter_region_celltype, ".*"), "", p))),collapse = delimiter_node))


          res = data.table::data.table(pvalue = pvalue,
                                       p.adjusted = p.adjusted,
                                       enrichness = zscore,
                                       neighborhood = noi,
                                       celltype = noi_no_region,
                                       region = noi_region,
                                       composition = noi_comp,
                                       target = j,
                                       sample_id = i)
        }else{

          res = data.table::data.table(pvalue = pvalue,
                                       p.adjusted = p.adjusted,
                                       enrichness = zscore,
                                       neighborhood = noi,
                                       composition = noi_comp,
                                       target = j,
                                       sample_id = i)
        }
        res = na.omit(res)
        data.table::setorder(res,p.adjusted,-enrichness)
        return(res)

      }

    })

    return(out)

  })

  names(results) = id_of_samples
  POLYspace@results$enrichment = results

  return(POLYspace)

}
