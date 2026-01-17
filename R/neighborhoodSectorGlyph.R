#' Generate sector-level data for neighborhood glyph construction
#'
#' Parse a neighborhood string and adjacency matrix to generate structured
#' sector-level data used for neighborhood sector glyph visualization.
#'
#' @param string A character string encoding neighborhood structure
#' 
#' @param delimiter_node Character delimiter separating nodes in \code{string}.
#'   Default is \code{"--"}.
#'   
#' @param delimiter_region_celltype Character delimiter separating region and
#'   cell type labels within a node. Default is \code{"_"}.
#'   
#' @param adj_matrix #' @param adj_matrix Either
#'   \itemize{
#'     \item a square adjacency matrix defining relationships between nodes
#'       in the neighborhood, or
#'     \item a character (e.g. singlet or pair) encoding adjacency relationships.
#'   }
#'
#' @return A data structure containing sector-level information 
#' for downstream glyph construction.
#' 
#' @export
generateSectorData = function(string,
                              delimiter_node = '--',
                              delimiter_region_celltype = '_',
                              adj_matrix){
  
  rn = value = composition = NULL
  
  labels = unlist(strsplit(string, delimiter_node, fixed = TRUE))
  regions = ifelse(grepl(delimiter_region_celltype, labels), 
                   sub(paste0(delimiter_region_celltype, ".*"), "", labels), NA)
  celltypes = sub(paste0(".*",delimiter_region_celltype), "", labels)
  

  data = data.frame(
    region = regions,
    celltype = celltypes,
    label = labels,
    id = 1:length(labels),
    stringsAsFactors = FALSE
  )
  
  if(identical('singlet',adj_matrix)){
    adj_matrix = matrix(1)
  }
  if(identical('pair',adj_matrix)){
    adj_matrix = matrix(c(0,1,1,0),2,2)
  }

  edge_full = which(adj_matrix == 1, arr.ind = TRUE)
  find_unique_edge = reshape2::melt(edge_full)
  find_unique_edge = data.table::as.data.table(find_unique_edge)
  colnames(find_unique_edge) = c('rn','variable','value')
  data.table::setorder(find_unique_edge, rn, value)
  find_unique_edge[, composition := paste(value, collapse = "--"), by = rn]
  urn = find_unique_edge[, .(urn = data.table::first(rn)), by = composition]
  edge_data = edge_full[urn$urn,,drop=FALSE]
  
  #edges = matrix(data$label[edge_data],ncol = 2)
  
  return(list(data = data, edges = edge_data))
}

#' Generate point positions within a sector
#'
#' Compute polar coordinates for points evenly distributed within a
#' sector defined by start and end angles.
#'
#' @param start_angle Numeric scalar giving the starting angle (in radians)
#'   of the sector.
#' @param end_angle Numeric scalar giving the ending angle (in radians)
#'   of the sector.
#' @param n_points Integer specifying the number of points to generate.
#' @param radius Numeric scalar specifying the radial distance of points
#'   from the center.
#'
#' @return A data.frame containing the angular and radial coordinates of
#'   generated points.
#'   @export
generateSectorPointPosition = function(start_angle, end_angle, n_points, radius) {
  
  r = radius
  
  if (n_points == 1) {
    angle = (start_angle + end_angle) / 2
    return(data.frame(x = r * sin(angle), y = r * cos(angle)))
  }

  angles = seq(start_angle, end_angle, length.out = n_points + 2)[-c(1, n_points + 2)]
  points = data.frame(
    x = r * sin(angles),
    y = r * cos(angles)
  )
  
  return(points)
}


#' Draw a neighborhood sector glyph
#'
#' Visualize a spatial neighborhood as a sector-based glyph, where regions
#' are represented as angular sectors and cell types are positioned within
#' each sector according to adjacency relationships.
#'
#' @param string A character string encoding the neighborhood structure,
#'   including region and cell type information.
#' @param delimiter_node Character delimiter separating nodes in \code{string}.
#'   Default is \code{"--"}.
#' @param delimiter_region_celltype Character delimiter separating region and
#'   cell type labels within a node. Default is \code{"_"}.
#' @param adj_matrix Either a square adjacency matrix defining relationships
#'   between nodes in the neighborhood, or a character vector encoding
#'   singlet or pairwise adjacency relationships.
#' @param radius Numeric scalar specifying the outer radius of the glyph.
#' @param cell_radius Numeric scalar specifying the radius used for placing
#'   cell-type points within each sector.
#' @param region_color_plate Named character vector mapping regions to colors.
#' @param celltype_color_plate Named character vector mapping cell types to
#'   colors.
#' @param region_seq Optional character vector specifying the order of regions
#'   around the glyph. If \code{NULL}, regions are ordered automatically.
#' @param region_alpha Numeric scalar controlling the transparency of region
#'   sectors.
#' @param edge_color Color used to draw edges between adjacent nodes.
#' @param sector_boundary_color Color of sector boundaries.
#' @param sector_boundary_linewidth Numeric scalar specifying the line width
#'   of sector boundaries.
#' @param edge_linewidth Numeric scalar specifying the line width of edges.
#' @param cell_size Numeric scalar specifying the size of cell-type points.
#' @param text Logical; whether to draw text labels.
#' @param text_vjust Numeric value controlling vertical justification of text.
#' @param text_size Numeric scalar specifying the size of text labels.
#'
#' @return A ggplot object representing the neighborhood sector glyph.
#' @export
neighborhoodSectorGlyph = function(
    string,
    delimiter_node = "--",
    delimiter_region_celltype = "_",
    adj_matrix,
    radius = 1,
    cell_radius = 0.5,
    region_color_plate = NULL,
    celltype_color_plate,
    region_seq = NULL,
    region_alpha = 0.3,
    edge_color = "gray30",
    sector_boundary_color = "black",
    sector_boundary_linewidth = 0.5,
    edge_linewidth = 0.8,
    cell_size = 3,
    text = TRUE,
    text_vjust = 2,
    text_size = 4
) {
  sector_info = generateSectorData(string = string,
                                   delimiter_node = delimiter_node,
                                   delimiter_region_celltype = delimiter_region_celltype,
                                   adj_matrix = adj_matrix)
  data = sector_info$data
  edges = sector_info$edges
  
  if(is.null(region_seq)){
    regions = sort(unique(na.omit(data$region)))
  }else{
    regions = region_seq
  }
  
  n_regions = length(regions)
  
  points_df = data.frame()
  sector_df = data.frame()
  
  if (n_regions >= 1) {
    angles = if(n_regions == 1) c(0, 2*pi) else seq(0, 2*pi, length.out = n_regions + 1)
    
    ###
    angles = -angles
    sector_df = data.frame(
      region = regions,
      start = angles[-(n_regions + 1)],
      end = angles[-1],
      color = region_color_plate[regions]
    )
    
    
    
    
    for (region in regions) {
      region_data = data[data$region == region, ]
      n_points = nrow(region_data)
      region_idx = which(sector_df$region == region)
      
      if(n_regions == 1){
        temp_positions = generateSectorPointPosition(0, 2*pi, n_points, cell_radius)
      } else {
        temp_positions = generateSectorPointPosition(
          sector_df$start[region_idx],
          sector_df$end[region_idx],
          n_points,
          cell_radius
        )
      }
      
      temp_df = data.frame(
        x = temp_positions$x,
        y = temp_positions$y,
        label = region_data$label,
        id = region_data$id,
        celltype = region_data$celltype,
        region = region
      )
      
      points_df = rbind(points_df, temp_df)
    }
  } else {
    temp_positions = generateSectorPointPosition(0, 2*pi, nrow(data), cell_radius)
    points_df = data.frame(
      x = temp_positions$x,
      y = temp_positions$y,
      label = data$label,
      celltype = data$celltype,
      id = data$id,
      region = NA
    )
  }
  points_df = points_df[match(data$id,points_df$id),]
  edges_df = data.frame()
  if(nrow(edges) > 0){
    for (i in 1:nrow(edges)) {
      from_point = points_df[edges[i, 1], c("x", "y")]
      to_point = points_df[edges[i, 2], c("x", "y")]
      edges_df = rbind(edges_df, data.frame(
        x = from_point$x, y = from_point$y,
        xend = to_point$x, yend = to_point$y
      ))
    }
  }
  
  p = ggplot2::ggplot()
  
  if(n_regions >= 1){
    p = p + ggforce::geom_arc_bar(
      data = sector_df,
      ggplot2::aes(x0=0,y0=0,r0=0,r=radius,start=start,end=end,fill=region),
      color=if(n_regions==1) NA else sector_boundary_color, alpha=region_alpha
    ) + ggplot2::scale_fill_manual(values = region_color_plate)
    p = p + ggforce::geom_circle(ggplot2::aes(x0=0,y0=0,r=radius),color=sector_boundary_color,linewidth=sector_boundary_linewidth)
    
  } else {
    p = p + ggforce::geom_circle(ggplot2::aes(x0=0,y0=0,r=radius),color=sector_boundary_color,linewidth=sector_boundary_linewidth)
  }
  
  if(nrow(edges_df) > 0){
    p = p + ggplot2::geom_segment(
      data = edges_df, ggplot2::aes(x=x,y=y,xend=xend,yend=yend),
      color=edge_color,linewidth=edge_linewidth
    )
  }
  
  if(text){
    p = p + 
      ggplot2::geom_point(data=points_df,ggplot2::aes(x,y,color =celltype),size=cell_size) +
      ggplot2::scale_color_manual(values = celltype_color_plate)+
      ggplot2::geom_text(data=points_df,ggplot2::aes(x,y,label=celltype),vjust=text_vjust,size=text_size) +
      ggplot2::coord_fixed() + ggplot2::theme_void()+ ggplot2::theme(legend.position = "none")
  }else{
    p = p + 
      ggplot2::geom_point(data=points_df,ggplot2::aes(x,y,color =celltype),size=cell_size) +
      ggplot2::scale_color_manual(values = celltype_color_plate)+
      ggplot2::coord_fixed() + ggplot2::theme_void()+ ggplot2::theme(legend.position = "none")
  }
  
  
  return(p)
}




