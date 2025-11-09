#' Generate the isomorphism agent for a target
#'
#' Builds isomorphism mappings for a target defined by an adjacency matrix.
#'
#' @param target An adjacency matrix defining the target.
#' @param delimiter_node A character string used to join node labels.
#'   Default is "--".
#'
#' @return A data.table containing the isomorphism mappings for
#'   \code{target}.
#'
#' @export
generateIsomorphismAgent = function(target,
                                    delimiter_node = '--'){

  rn = canonical = NULL

  n_cells = unique(dim(target))
  target_graph = igraph::graph_from_adjacency_matrix(target,mode='undirected')

  partition = partitions::parts(n_cells)
  color = apply(partition, 2, function(x){

    cop = unique(combinat::permn(rep(1:length(x),times = x)))
    cop_order = data.table::as.data.table(do.call('rbind',cop))
    cop = cop[do.call(base::order, cop_order)]

    return(cop)

  })

  membership = lapply(color, function(x){

    nx = length(x)

    if(nx ==1){
      return(setNames(1L, "1"))
    }else{

      iso_pair_index = data.frame(t(combn(nx, 2)))

      for(i in 1:nrow(iso_pair_index)){
        g1 = target_graph
        igraph::V(g1)$color = x[[iso_pair_index[i,1]]]
        g2 = target_graph
        igraph::V(g2)$color = x[[iso_pair_index[i,2]]]
        iso_pair_index[i,'iso'] = igraph::isomorphic(g1, g2, method = "vf2")
      }

      if(sum(iso_pair_index$iso)>0){

        g = igraph::graph_from_data_frame(iso_pair_index[iso_pair_index$iso == TRUE,,drop=FALSE], directed = FALSE)
        comp = igraph::components(g)
        comp = comp$membership

        if(nx != length(comp)){
          extra = setdiff(seq_len(nx),names(comp))
          extra_comp = (max(comp)+1):(max(comp)+length(extra))
          names(extra_comp) = extra
          comp = c(comp,extra_comp)
        }
      }else{
        comp = seq_len(nx)
        names(comp) = seq_len(nx)
      }

      comp = comp[as.character(seq_len(nx))]
      return(comp)
    }

  })

  canonical = lapply(seq_along(membership),function(x){
    return(paste(x,membership[[x]],sep='_'))
  })
  canonical = do.call('c',canonical)

  color = do.call('c',color)
  identifier = do.call('rbind',color)
  identifier = data.table::as.data.table(identifier)
  identifier =  identifier[, do.call(paste, c(.SD, sep = delimiter_node))]

  find_canonical_color = data.table::data.table(identifier = identifier,
                                                canonical = canonical)
  find_canonical_color[, rn := .I]
  urn_index = find_canonical_color[, .(urn = data.table::first(rn)), by = .(canonical)]
  canonical_color = color[urn_index$urn[data.table::chmatch(canonical,urn_index$canonical)]]

  agent = data.table::data.table(identifier = identifier)
  agent[, 'mapping'] = Map(match, canonical_color, color)

  return(agent)
}
