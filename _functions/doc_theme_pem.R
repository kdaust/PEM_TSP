# Copyright 2021 Province of British Columbia
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
# http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and limitations under the License.

require(extrafont)
require(hrbrthemes)
require(ggthemes)
library(ggplot2)
library(ggpubr)
require(colorspace)
 #font_import() ## use this once to load up the fonts to C://Windows/Fonts

## On Windows machines need to explicitly load specific fonts to use in the scripts
loadfonts(device = "win")
 windowsFonts(`Arial Narrow` = windowsFont("Arial Narrow"))
 windowsFonts(Times=windowsFont("TT Times New Roman"))
 windowsFonts(Helvetica=windowsFont("Helvetica")) 

 # Theming scripts to use in publications
theme_pem <- function(base_size=12, base_family="Helvetica") {
  thm <- theme_pem_foundation(base_size = base_size, base_family = base_family)
  thm
}

theme_pem_facet <- function(base_size = 12, base_family = "Helvetica") {
  
  theme_pem_foundation(base_size = base_size, base_family = base_family) + 
    theme(
      panel.spacing = unit(.6,"lines"),
      panel.border = element_rect(colour = "black", fill = NA),
      strip.background = element_rect(colour = "black", fill = "grey85"))
  
}

theme_pem_foundation <- function(base_size, base_family) {
  theme_few(
    base_size = base_size,
    base_family = base_family) + 
    theme(
      text = element_text(colour = "black"),
      line = element_line(colour = "black", size = 0.5,
                          linetype = 1, lineend = "butt"),
      rect = element_rect(fill = "white", colour = "black",
                          size = 0.5, linetype = 1),
      axis.line = element_line(colour = "black"),
      axis.line.x = element_line(colour = "black"),
      axis.line.y = element_line(colour = "black"),
      axis.text = element_text(colour = 'black'),
      axis.text.y = element_text(hjust = 1),
      axis.ticks = element_blank(),
      plot.title = element_text(vjust = 2),
      legend.title = element_text(face = "plain"),
      panel.background = element_blank(),
      panel.border = element_blank(),
      #panel.grid.minor = element_blank(),
      #panel.grid.major = element_line(colour = "grey80",size = 0.5),
      axis.title.y = element_text(vjust = 1, angle = 90),
      axis.title.x = element_text(vjust = 0),
      #panel.spacing = unit(0.25, "lines"),
      plot.background = element_blank(),
      legend.key = element_blank(),
      #complete = TRUE
    )
}

bgc_colours <- function(which = NULL) {
  cols <- c(BAFA = "#E5D8B1",
            SWB  = "#A3D1AB",
            BWBS = "#ABE7FF",
            ESSF = "#9E33D3",
            CMA  = "#E5C7C7",
            SBS  = "#2D8CBD",
            MH   = "#A599FF",
            CWH  = "#208500",
            ICH  = "#85A303",
            IMA  = "#B2B2B2",
            SBPS = "#36DEFC",
            MS   = "#FF46A3",
            IDF  = "#FFCF00",
            BG   = "#FF0000",
            PP   = "#DE7D00",
            CDF  = "#FFFF00")
  
  if (is.null(which)) {
    return(cols)
  } else {
    if (!all(which %in% names(cols))) stop("Unknown Biogeoclimatic Zone code(s) specified", call. = FALSE)
    return(cols[which])
  }
}

###from BBC package for saving plots to files

save_plot <- function (plot_grid, width, height, filenamepath) {
  grid::grid.draw(plot_grid)
  #save it
  ggplot2::ggsave(filename = filenamepath, 
                  plot=plot_grid, width=(width/72), height=(height/72),  bg="white", dpi = 600)
}

#Left align text
left_align <- function(plot_name, pieces){
  grob <- ggplot2::ggplotGrob(plot_name)
  n <- length(pieces)
  grob$layout$l[grob$layout$name %in% pieces] <- 2
  return(grob)
}

create_footer <- function (source_name, logo_image_path) {
  #Make the footer
  footer <- grid::grobTree(grid::linesGrob(x = grid::unit(c(0, 1), "npc"), y = grid::unit(1.1, "npc")),
                           grid::textGrob(source_name,
                                          x = 0.004, hjust = 0, gp = grid::gpar(fontsize=16)),
                           grid::rasterGrob(png::readPNG(logo_image_path), x = 0.944))
  return(footer)
  
}

#' Arrange alignment and save BBC ggplot chart
#'
#' Running this function will save your plot with the correct guidelines for publication for a BBC News graphic.
#' It will left align your title, subtitle and source, add the BBC blocks at the bottom right and save it to your specified location.
#' @param plot_name The variable name of the plot you have created that you want to format and save
#' @param source_name The text you want to come after the text 'Source:' in the bottom left hand side of your side
#' @param filenamepath Exact filepath that you want the plot to be saved to
#' @param width_pixels Width in pixels that you want to save your chart to - defaults to 640
#' @param height_pixels Height in pixels that you want to save your chart to - defaults to 450
#' @param logo_image_path File path for the logo image you want to use in the right hand side of your chart,
#'  which needs to be a PNG file - defaults to BBC blocks image that sits within the data folder of your package
#' @return (Invisibly) an updated ggplot object.

#' @keywords finalise_plot
#' @examples
#' finalise_plot(plot_name = myplot,
#' source = "The source for my data",
#' save_filepath = "filename_that_my_plot_should_be_saved_to-nc.png",
#' width_pixels = 640,
#' height_pixels = 450,
#' logo_image_path = "logo_image_filepath.png"
#' )
#'
#' @export
finalise_plot <- function(plot_name,
                          #save_filepath=file.path(Sys.getenv("TMPDIR"), "tmp-nc.pdf"),
                          filenamepath,
                          width_pixels=640,
                          height_pixels=450)#,
                          #logo_image_path = file.path(system.file("data", package = 'bbplot'),"placeholder.png")) 
  
  
  {
  
    #export_plot <- function(plot_name, source, save_filepath, width_pixels = 640, height_pixels = 450) {
  theme_pem() #+ 

  #Draw your left-aligned grid
  plot_left_aligned <- left_align(plot_name, c("subtitle", "title", "caption"))
  plot_grid <- plot_left_aligned
  #plot_grid <- ggpubr::ggarrange(plot_left_aligned, #footer,ncol = 1, nrow = 1)#,
                               
  #                                heights = c(1, 0.045/(height_pixels/450)))
  ## print(paste("Saving to", save_filepath))
  #save_plot(plot_grid, width_pixels, height_pixels, save_filepath)
save_plot(plot_left_aligned, width_pixels, height_pixels, filenamepath)
  ## Return (invisibly) a copy of the graph. Can be assigned to a
  ## variable or silently ignored.
  invisible(plot_grid)
}


