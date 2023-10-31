# GNU General Public License v3.0 (https://github.com/IanevskiAleksandr/sc-type/blob/master/LICENSE)
# Written by Aleksandr Ianevski <aleksandr.ianevski@helsinki.fi>, June 2021
#
# modified by Oliver Burton <ob240@cam.ac.uk>, August 2023
# Functions on this page:
# prepare_marker_lists: prepare marker lists from input Cell Type excel file
#
# @params: path_to_db_file - DB file with cell types
# @tissue - source of cells (e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain)
#

prepare_marker_lists <- function(path_to_db_file, tissue = "Immune", parent.cell.type = 0){
  
  cell_markers <- read_xlsx(path_to_db_file)
  
  cell_markers <- cell_markers[cell_markers$Tissue.restricted %in% tissue,]
  # filter based on parent cell type, with option to skip if unselected
  if ( length(parent.cell.type >1 )){
    parent.cell.type <- paste0(parent.cell.type, collapse = "|")
    cell_markers <- cell_markers[grep(parent.cell.type, cell_markers$Family.tree),]
  } else if ( parent.cell.type !=0 ){
    cell_markers <- cell_markers[grep(parent.cell.type, cell_markers$Family.tree),]
  } else {
    cell_markers <- cell_markers
  }
  
  cell_markers$Positive.Markers <- gsub(" ","",cell_markers$Positive.Markers)
  cell_markers$Negative.Markers <- gsub(" ","",cell_markers$Negative.Markers)
  
  cell_markers$Positive.Markers <- gsub("///",",",cell_markers$Positive.Markers)
  cell_markers$Positive.Markers <- gsub(" ","",cell_markers$Positive.Markers)
  cell_markers$Negative.Markers <- gsub("///",",",cell_markers$Negative.Markers)
  cell_markers$Negative.Markers <- gsub(" ","",cell_markers$Negative.Markers)
   
  Pos.Markers <- lapply(1:nrow(cell_markers), function(x) gsub(" ","",unlist(strsplit(toString(cell_markers$Positive.Markers[x]),","))))
  names(Pos.Markers) <- cell_markers$Cell.type
  Neg.Markers <- lapply(1:nrow(cell_markers), function(x) gsub(" ","",unlist(strsplit(toString(cell_markers$Negative.Markers[x]),","))))
  names(Neg.Markers) <- cell_markers$Cell.type
  
  list(markers_positive = Pos.Markers, markers_negative = Neg.Markers)
}
