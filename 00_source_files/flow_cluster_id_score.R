# Adapted from scType
# GNU General Public License v3.0 (https://github.com/IanevskiAleksandr/sc-type/blob/master/LICENSE)
# Written by Aleksandr Ianevski <aleksandr.ianevski@helsinki.fi>, June 2021
#
# edited by Oliver Burton <ob240@cam.ac.uk>, August 2023
#
# Functions on this page:
# flow_cluster_id_score: calculate cluster ID scores and assign cell types
#
# @params: cluster.input.data - input transposed thresholded flow expression matrix (rownames - markers, column names - clusters),
#
# @params: marker_pos - list of markers positively expressed in the cell type 
# @params: marker_neg - list of markers that should not be expressed in the cell type (NULL if not applicable)

flow_cluster_id_score <- function(cluster.input.data, marker_pos, marker_neg = NULL, ...){
  
  # subset to markers found in the data
  names_mkp_cp <- names(marker_pos)
  names_mkn_cp <- names(marker_neg)
  
  marker_pos <- lapply(1:length(marker_pos), function(x){ 
    MarkerToKeep = rownames(cluster.input.data) %in% as.character(marker_pos[[x]])
    rownames(cluster.input.data)[MarkerToKeep]})
  
  marker_neg = lapply(1:length(marker_neg), function(x){ 
    MarkerToKeep = rownames(cluster.input.data) %in% as.character(marker_neg[[x]])
    rownames(cluster.input.data)[MarkerToKeep]})
  
  names(marker_pos) <- names_mkp_cp
  names(marker_neg) <- names_mkn_cp
  
  # subselect only with marker genes
  cluster.input.data = cluster.input.data[unique(c(unlist(marker_pos),unlist(marker_neg))), ]
  
  # deal with cases with only one marker
  if (is.null(ncol(cluster.input.data))){
    cluster.input.data <- data.frame(matrix(cluster.input.data,1))
    row.names(cluster.input.data) <- unique(c(unlist(marker_pos),unlist(marker_neg)))
  }
  
  # combine scores
  combined.score = do.call("rbind", lapply(names(marker_pos), function(mkrs){ 
    sapply(1:ncol(cluster.input.data), function(x) {
      mk_pos = cluster.input.data[marker_pos[[mkrs]], x]
      mk_neg = cluster.input.data[marker_neg[[mkrs]], x] * -1
      
      sum_t1 = (sum(mk_pos) / sqrt(length(mk_pos)))
      sum_t2 = sum(mk_neg) / sqrt(length(mk_neg))
      
      if(is.na(sum_t2)){
        sum_t2 = 0;
      }
      sum_t1 + sum_t2
    })
  })) 
  
  
  dimnames(combined.score) = list(names(marker_pos), colnames(cluster.input.data))
  combined.score.max <- combined.score[!apply(is.na(combined.score) | combined.score == "", 1, all),] # remove na rows
  
  combined.score.max
}
