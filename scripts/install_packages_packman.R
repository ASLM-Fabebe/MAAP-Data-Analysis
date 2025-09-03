# First install pacman

if (!require("pacman", quietly = TRUE)){
  install.packages("pacman")
  library(pacman)
}


# Install the rest of the packages using pacman

pacman::p_load(
  AMR,
  parallel,
  tidyverse,
  readxl,
  remotes,
  ggplot2,
  colorRamp2,
  ggalt,
  lubridate,
  RColorBrewer,
  openxlsx,
  naniar,
  glue,
  readr,
  readxl,
  ggrepel,
  ggtext,
  BiocManager,
  ggrepel,
  rhandsontable,
  shiny,
  writexl,
  zoo,
  shinythemes,
  rhandsontable
)


# if (!require("ComplexHeatmap", quietly = TRUE)){
#   BiocManager::install("ComplexHeatmap")
#   library(ComplexHeatmap)
# }
#
#
# if (!require("tidyHeatmap", quietly = TRUE)){
#   install.packages("tidyHeatmap")
#   library(tidyHeatmap)
# }
#
# if (!require("ggtree", quietly = TRUE)){
#   BiocManager::install("ggtree")
#   library(ggtree)
# }
