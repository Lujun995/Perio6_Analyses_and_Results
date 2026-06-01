library(ggplot2)
library(cowplot)
library(stringr)

paint_tooth <- function(t_ids, colors, output_file = "temp.svg"){
  names(colors) <- t_ids
  tooth_id_to_path <-
    c("2994", "5171", "5187", "3063", "3876", "3892", "3900", "3908",
      "4364", "4356", "4348", "4332", "4320", "4306", "4290", "4278",
      "4368", "4382", "4400", "4416", "4426", "4436", "4452", "4454",
      "3899", "3897", "3881", "3870", "3962", "3946", "3928", "3912")
  names(tooth_id_to_path) <- sprintf("%02d", 1:32)
  svg <- xml2::read_xml("../tooth numbering/Universal_Numbering_System.svg")
  for(t_id in t_ids){
    tooth_path <- 
      xml2::xml_find_first(svg, str_glue("//*[@id='path{tooth_id_to_path[t_id]}']"))
    style <- xml2::xml_attr(tooth_path, "style")
    style <- sub("fill:[^;]*", str_glue("fill:{colors[t_id]}"), style)
    xml2::xml_set_attr(tooth_path, "style", style)
  }
  xml2::write_xml(svg, output_file) # have to write it to disk first
  # otherwise, the image will be displayed wrongly during format conversion?
}


plot_tooth_heatmap <- 
  function(values, option = "rocket", begin = 0.6, end = 1, direction = -1,
           n_colors = 10, output_svg = "heatmap.svg", legend_title = "Value") {
  # matrix. plot multiple heatmaps with shared legend
  if(is.matrix(values)) {
    return(plot_tooth_heatmap_matrix(values, option, begin, end, direction, 
                                     n_colors, output_svg, legend_title))
  }
  
  # vector
  if(is.null(names(values))) {
    stop("`values` must be a named numeric vector, with names as tooth IDs.")
  }
  return(plot_single_tooth_heatmap(values, option, begin, end, direction, 
                                   n_colors, output_svg, legend_title))
}

# basic function
plot_single_tooth_heatmap <- function(values, option, begin, end, direction, 
                                      n_colors, output_svg, legend_title, 
                                      plot_title = NULL, global_range = NULL) {
  # values should be a vector here
  if(is.matrix(values))
    stop("values should be a vector in plot_single_tooth_heatmap")
  # use global range
  val_range <- if(!is.null(global_range)) global_range else 
    range(values, na.rm = TRUE)
  
  
  
  # palette (should be consistent with the legend below)
  pal <- substr(
    viridisLite::viridis(n_colors, begin = begin, end = end, 
                         direction = direction, option = option),
    1, 7
  )
  breaks <- seq(val_range[1], val_range[2], length.out = n_colors + 1)
  # Map values to colors
  bin_indices <- cut(values, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  colors <- pal[bin_indices]
  
  # paint teeth
  paint_tooth(t_ids = names(values), colors = colors, output_file = output_svg)
  svg_magick <- magick::image_read(rsvg::rsvg(output_svg, width = 500)) # read back
  
  # create the plot with both SVG and legend
  legend_data <- data.frame(
    x = 0.5, # center the data
    y = seq(val_range[1], val_range[2], length.out = n_colors+5),
    value = seq(val_range[1], val_range[2], length.out = n_colors+5)# dummy values
  )
  p <- ggplot(legend_data, aes(x = x, y = value, fill = value)) +
    annotation_custom(
      grid::rasterGrob(svg_magick, interpolate = TRUE),
      xmin = -2, xmax = 2,  # the scale does not matter. It will be centered at (0,0)
      ymin = val_range[1], ymax = val_range[2]) +
    geom_tile(data = data.frame(x = 10, y = legend_data$y, 
                                value = legend_data$value),
              aes(x = x, y = y, fill = value), 
              width = 0.1, alpha = 0) +  # make it invisible
    scale_fill_viridis_b(
      option = option, 
      begin = begin, 
      end = end, 
      direction = direction, 
      limits = val_range, # the same limits
      name = legend_title,
      breaks = breaks) +  # the same breaks!
    coord_cartesian(xlim = c(-2.1, 2.1), ylim = val_range, expand = FALSE) +
    theme_void() + 
    theme(legend.key.height = unit(0.05, "npc"))
  if (!is.null(plot_title)) {
    p <- p + ggtitle(plot_title, ) +
      theme(plot.title = element_text(size = 11, hjust = 0.5))
  }
  return(p)
}

plot_tooth_heatmap_matrix <- function(values, option, begin, end, direction, 
                                      n_colors, output_svg, legend_title) {
  if((!is.matrix(values)) | (nrow(values) == 1))
    stop("`values` is not a matrix.")
  if(is.null(rownames(values)) || is.null(colnames(values))) {
    stop("`values` matrix must have both row names (features) and column names (tooth IDs).")
  }
  
  global_val_range <- range(values, na.rm = TRUE)
  plot_list <- list()

  # loop over rows
  for(i in 1:nrow(values)){
    feature_name <- rownames(values)[i]
    feature_values <- values[i, ]
    names(feature_values) <- colnames(values)
    feature_values <- feature_values[!is.na(feature_values)]
    if(length(feature_values) == 0){
      warning(paste("No valid values for feature:", feature_name))
      next
    }
    
    temp_svg <- gsub("\\.svg$", paste0("_", feature_name, ".svg"), output_svg)
    p <- plot_single_tooth_heatmap(
      feature_values, option, begin, end, direction, 
      n_colors, temp_svg, legend_title, 
      plot_title = feature_name,
      global_range = global_val_range # remember to set this
    )
    
    plot_list[[feature_name]] <- p
  }
  
  final_plot <- 
    do.call(what = lemon::grid_arrange_shared_legend,
            args = c(plot_list, ncol = length(plot_list), nrow = 1,
                     position = "right", plot = FALSE))
  
  return(final_plot)
}
