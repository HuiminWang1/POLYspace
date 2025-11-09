#' Get the canonical label for each neighborhood
#'
#' Computes the canonical label for each neighborhood.
#'
#' @param subs Neighborhoods identified by POLYspace.
#' @param ta   The shape of the neighborhoods.
#' @param ag   The isomorphism mappings.
#' @param delimiter_node A delimiter used when constructing labels. Default is \code{"--"}.
#' @param noi  Neighborhoods of interest.
#' @param coi Compositions of interest.
#'
#' @return Annotations to each neighborhood.
#'
#' @export
isomorphismMapping = function(subs,
                              ta,
                              ag,
                              delimiter_node = '--',
                              noi = NULL,
                              coi = NULL){

  rn = value = freq = rank = composition = NULL


  if(identical('pair',ta)){

    subs[, rn := .I]
    subs = data.table::melt(subs, id.vars = "rn")
    data.table::setorder(subs, rn, value)
    if(!is.null(coi)){
      subs[, composition := paste(value, collapse = delimiter_node), by = rn]
      subs = subs[subs$composition %in% coi,]
    }
    subs = subs[,c('rn','value')]


  }else{
    n_cells = ncol(ta)
    colnames(subs) = as.character(1:ncol(subs))

    subs[, rn := .I]
    subs = data.table::melt(subs, id.vars = "rn")
    subs$variable = as.numeric(subs$variable)
    data.table::setorder(subs, rn, value)


    subs[, composition := paste(value, collapse = delimiter_node), by = rn]
    if(!is.null(coi)){
      subs = subs[subs$composition %in% coi,]
    }

    subs[, freq := .N, by = .(composition, value)]
    subs[, rank := data.table::frank(list(-freq, value), ties.method = "dense"), by = composition]
    data.table::setorder(subs, rn, variable)
    subs[, rank_identifier := paste(rank, collapse = delimiter_node), by = rn]
    subs_frame = subs[,c('rn','rank_identifier')]

    subs = subs[,c('rn','value')]

    subs_frame = subs_frame[!duplicated(subs_frame)]
    data.table::setorder(subs_frame, rn)

    ag_mapping = ag$mapping[match(subs_frame$rank_identifier,ag$identifier)]
    ag_mapping = Map(`+`, ag_mapping, (seq_along(ag_mapping) - 1) * n_cells)
    ag_mapping = unlist(ag_mapping)
    subs$value = subs$value[ag_mapping]

  }


  canonical_annotation = subs[, paste(value, collapse = delimiter_node), by = rn][["V1"]]
  if(!is.null(noi)){
    canonical_annotation = canonical_annotation[canonical_annotation %in% noi]
  }
  if(length(canonical_annotation)>0){
    return(canonical_annotation)
  }


}

