# run clustering

# run clustering via FlowSOM or Phenograph----------
cat(
  "\n
As the first part of the analysis, the script will cluster your cells into groups.

For this, you can use either Phenograph or FlowSOM.

Phenograph is fast, and will automatically determine how many clusters to generate.
Sometimes it overclusters, particularly in cases with lots of homogeneous cells.

If you prefer to use FlowSOM, you'll need to decide how many clusters you want to find.

One strategy that may be helpful is first run Phenograph, assess approximately how
many clusters appear to really be distinct, then repeat the analysis using FlowSOM
with a defined number of clusters.
\n"
)

clustering.method <- menu( c("Phenograph", "FlowSOM"),
                           title = "Please choose your clustering approach.")

cat("\n")

if (clustering.method == 2){
  fcs.cluster.n <- readline( "Please enter the number of clusters to find using FlowSOM: ")
  fcs.cluster.n <- as.numeric(fcs.cluster.n)
  
  cat("\n")
  
  fcs.cluster <- sprintf( "%02d", 1 : fcs.cluster.n )
  
  fcs.cluster.label <- sprintf( "Cluster-%s", fcs.cluster )
  names( fcs.cluster.label ) <- fcs.cluster
  
  fcs.cluster.color <- rep( 
    fcs.color.pool, 
    ceiling( fcs.cluster.n / fcs.color.pool.n ) 
  )[ 1 : fcs.cluster.n ]
  names( fcs.cluster.color ) <- fcs.cluster
  
  fcs.cluster.line.type <- rep( 
    fcs.line.type.pool, 
    ceiling( fcs.cluster.n / fcs.line.type.pool.n ) 
  )[ 1 : fcs.cluster.n ]
  names( fcs.cluster.line.type ) <- fcs.cluster
  
  fcs.flow.som.dim <- 10
  
  # get FlowSOM clusters---------------
  cat("Clustering data with FlowSOM via EmbedSOM\n")
  set.seed.here( fcs.seed.base, "get flowsom clusters" )
  
  # obtain fast approximate knn
  non.self.knn <- 1
  knn.M <- 32
  while( length(non.self.knn) > 0) {
    print(knn.M)
    set.seed.here( fcs.seed.base, "get knn with Rcpp" )
    dmrd.knn <- RcppHNSW::hnsw_knn( dmrd.data, k= fcs.tsne.perplexity,
                                    distance= 'l2', n_threads= fcs.tsne.threads.n, 
                                    ef = nrow(dmrd.data)/10, M = knn.M, 
                                    ef_construction = nrow(dmrd.data)/10 )
    test.dat.self.idx <- sapply( 1 : nrow(dmrd.knn$idx), function( ri ) {
      ri.idx <- which( dmrd.knn$idx[ ri, ] == ri )
      ifelse( length( ri.idx ) == 1, ri.idx, NA )
    } )
    
    non.self.knn <- test.dat.self.idx[is.na(test.dat.self.idx)]
    
    knn.M <- knn.M*2
    
  }
  
  # build som objects--embedsom method
  flow.som <- EmbedSOM::SOM(dmrd.data, xdim = fcs.flow.som.dim, 
                            ydim = fcs.flow.som.dim, batch = TRUE,
                            parallel = TRUE, threads = fcs.tsne.threads.n )
  
  # get clusters
  flow.som.mapping <- flow.som$mapping[ , 1 ]
  flow.som.codes <- flow.som$codes
  
  # get clusters from som mapping
  consensus.cluster <- ConsensusClusterPlus( t( flow.som.codes ),
                                             maxK = fcs.cluster.n, reps = 100, pItem = 0.9, pFeature = 1,
                                             clusterAlg = "hc", innerLinkage = "average", finalLinkage = "average",
                                             distance = "euclidean", 
                                             seed = set.seed.here( fcs.seed.base, "get clusters from som mapping" ) )
  
  flow.som.event.cluster <- consensus.cluster[[ fcs.cluster.n ]]$
    consensusClass[ flow.som.mapping ]
  
  # reorder clusters from bigger to smaller
  flow.som.cluster.rank <- 1 + fcs.cluster.n - 
    rank( table( flow.som.event.cluster ), ties.method = "last" )
  flow.som.event.cluster <- flow.som.cluster.rank[ flow.som.event.cluster ]
  names( flow.som.event.cluster ) <- NULL
  
  # set clusters as a factor
  flow.som.event.cluster <- factor( flow.som.event.cluster, 
                                    levels = 1 : fcs.cluster.n )
  levels( flow.som.event.cluster ) <- fcs.cluster
  
  dmrd.event.cluster <- flow.som.event.cluster
  dmrd.event.cluster.n <- as.vector( table( dmrd.event.cluster ) )
  names( dmrd.event.cluster.n ) <- fcs.cluster
  
  length( dmrd.event.cluster )
  table( dmrd.event.cluster )
  
} else {
  # get phenograph clusters---------------
  cat("Clustering data with Phenograph\n")
  
  # obtain fast approximate knn
  non.self.knn <- 1
  knn.M <- 32
  while( length(non.self.knn) > 0) {
    print(knn.M)
    set.seed.here( fcs.seed.base, "get knn with Rcpp" )
    dmrd.knn <- RcppHNSW::hnsw_knn( dmrd.data, k= fcs.tsne.perplexity,
                                    distance= 'l2', n_threads= fcs.tsne.threads.n, 
                                    ef = nrow(dmrd.data)/10, M = knn.M, 
                                    ef_construction = nrow(dmrd.data)/10 )
    test.dat.self.idx <- sapply( 1 : nrow(dmrd.knn$idx), function( ri ) {
      ri.idx <- which( dmrd.knn$idx[ ri, ] == ri )
      ifelse( length( ri.idx ) == 1, ri.idx, NA )
    } )
    
    non.self.knn <- test.dat.self.idx[is.na(test.dat.self.idx)]
    
    knn.M <- knn.M*2
    
  }
  
  ind <- dmrd.knn$idx
  links <- FastPG::rcpp_parallel_jce( ind )
  links <- FastPG::dedup_links( links )
  clusters <- FastPG::parallel_louvain( links )
  
  phenograph.event.cluster <- factor(clusters$communities)
  
  fcs.cluster.n <- length(unique(phenograph.event.cluster))
  phenograph.cluster.rank <- 1 + fcs.cluster.n - 
    rank( table( phenograph.event.cluster ), ties.method = "last" )
  phenograph.event.cluster <- phenograph.cluster.rank[ phenograph.event.cluster ]
  names( phenograph.event.cluster ) <- NULL
  
  fcs.cluster <- sprintf( "%02d", 1 : fcs.cluster.n )
  
  fcs.cluster.label <- sprintf( "Cluster-%s", fcs.cluster )
  names( fcs.cluster.label ) <- fcs.cluster
  
  fcs.cluster.color <- rep( 
    fcs.color.pool, 
    ceiling( fcs.cluster.n / fcs.color.pool.n ) 
  )[ 1 : fcs.cluster.n ]
  names( fcs.cluster.color ) <- fcs.cluster
  
  fcs.cluster.line.type <- rep( 
    fcs.line.type.pool, 
    ceiling( fcs.cluster.n / fcs.line.type.pool.n ) 
  )[ 1 : fcs.cluster.n ]
  names( fcs.cluster.line.type ) <- fcs.cluster
  
  # set clusters as a factor
  phenograph.event.cluster <- factor( phenograph.event.cluster, 
                                      levels = 1 : fcs.cluster.n )
  levels( phenograph.event.cluster ) <- fcs.cluster
  
  # reorder events
  dmrd.event.cluster <- phenograph.event.cluster
  
  dmrd.event.cluster.n <- as.vector( table( dmrd.event.cluster ) )
  names( dmrd.event.cluster.n ) <- fcs.cluster
  
  length( dmrd.event.cluster )
  table( dmrd.event.cluster )
}



## cluster naming------------

cat(
"Next, the script will try to identify and name the cell types present in every cluster.
For this to work best, you'll need to enter three pieces of information:
1) Whether you're using human or mouse cells
2) Which tissues you're using
3) If you've pre-selected only certain cell types, which cells those are.
\n"
)

cat(
"Please consider the automated naming as a guide, and review the names by checking the heatmaps
and density plots for the clusters. If you don't see the correct cell types being identified,  
check the instructions for cluster naming and add your cell type definitions to the database spreadsheet.
\n"
)

# select correct databases based on species---------

species.options <- c("Mouse", "Human")

species.used <- menu( c("Mouse", "Human"), 
                      title = "Please select the species your cells come from. ")

cell.marker.filename <- paste0(species.options[species.used], "_marker_names.xlsx")
cell.type.filename <- paste0(species.options[species.used], "_celltype_database.xlsx")

cell.database <- read_xlsx( file.path(fcs.src.dir, cell.type.filename) )

# restrict cell types based on tissues used---------

cat("\n
To identify tissue-specific cell types, tell the script which tissue sources
you've used. Immune cells is the default, so include this in addition to any
extra sources you've used.")
tissue.options <- unique(cell.database$Tissue.restricted)
tissue.type <- multi.menu( tissue.options, title = "Select the tissue or tissues your cells come from. ")

tissue.type <- unique(cell.database$Tissue.restricted)[tissue.type]

if ( length(tissue.type) > 1 ){
  cell.database <- dplyr::filter(cell.database, Tissue.restricted %in% tissue.type)
} else if ( tissue.type == 0 ){
  cell.database <- cell.database
} else {
  cell.database <- dplyr::filter(cell.database, Tissue.restricted %in% tissue.type)
}


# restrict cell ID based on known input-------
selected.cell.type <- multi.menu( unique(cell.database$Parent.cell.type), 
                                  title = "If you have pre-gated on specific cell types, please them here.
                       If you're using all viable cells, select 1.")

selected.cell.type <- unique(cell.database$Parent.cell.type)[selected.cell.type]


# match markers to database----------
marker.synonyms <- read_xlsx( file.path(fcs.src.dir, cell.marker.filename) )

marker.alternates <- lapply(1:nrow(marker.synonyms), 
                            function(x) gsub(" ","",unlist(strsplit(toString(marker.synonyms$Alternates[x]),","))))
marker.alternates <- lapply(1:nrow(marker.synonyms), function(x)
                            gsub("\\*"," ", unlist(marker.alternates[x])))
names(marker.alternates) <- marker.synonyms$Marker.Name

new.marker.names <- lapply( 1:length(fcs.channel.label), function(x){
  names(marker.alternates)[grep( paste0( fcs.channel.label[x], "\\b"), marker.alternates,
                                 ignore.case = TRUE)]
} )


# scale data for cluster matching--------
flow.data.cluster.median <- apply( dmrd.data, 2, tapply, 
                                   dmrd.event.cluster, median )

scaled.fcs.cluster.median <- scale(flow.data.cluster.median, scale = FALSE)

colnames(scaled.fcs.cluster.median) <- new.marker.names

# prepare marker lists and match clusters-----------

source( file.path( fcs.src.dir, "prepare_marker_lists.r" ) )

source( file.path( fcs.src.dir, "flow_cluster_id_score.r" ) )

marker.list <- prepare_marker_lists( file.path( fcs.src.dir, cell.type.filename ), 
                                             tissue.type, selected.cell.type )

auto.fcs.cluster.label <- fcs.cluster.label


# match clusters to cell type database

scaled.id.score <- flow_cluster_id_score(t(scaled.fcs.cluster.median),
                                                 marker_pos = marker.list$markers_positive, 
                                                 marker_neg = marker.list$markers_negative )

for (cluster in 1:ncol(scaled.id.score)){
  auto.fcs.cluster.label[cluster] <- names( which.max(scaled.id.score[,cluster]) )
}

# plot heatmaps of cluster ID scores--------

jpeg( file.path(fcs.density.figure.dir, "cluster_id_heatmap.jpg"), 
      width = 1000, height = 1000 )
heatmap(scaled.id.score, Rowv = NA, Colv = NA, scale = "none",
        margins = c(5,10),
        xlab = "Clusters", ylab = "matching cell types")
dev.off()

jpeg( file.path(fcs.density.figure.dir, "cluster_id_heatmap_dendro.jpg"), 
      width = 1000, height = 1000 )
heatmap(scaled.id.score, scale = "none",
        margins = c(5,10),
        xlab = "Clusters", ylab = "matching cell types")
dev.off()



# differentiate similar clusters by variable markers----------------

rownames(flow.data.cluster.median) <- auto.fcs.cluster.label
colnames(flow.data.cluster.median) <- fcs.channel.label

# find (positions of) clusters with non-unique names
non.unique.cluster.names <- auto.fcs.cluster.label[duplicated(auto.fcs.cluster.label)]
# sort into groups
non.unique.groups <- unique(non.unique.cluster.names)

# append high variance marker labels to all clusters until unique

for (cluster.group in non.unique.groups) {
  
  # collect clusters needing renaming
  data.temp <- flow.data.cluster.median[ grepl( cluster.group, rownames(flow.data.cluster.median)),]
  
  # find most variable channels
  variances <- apply(data.temp, 2, var)
  sorted.variances <- sort(variances, decreasing = TRUE)
  
  # find min and max (neg and pos) for each channel
  max.expression <- apply(data.temp, 2, max)
  min.expression <- apply(data.temp, 2, min)
  
  for( cluster in 1:nrow(data.temp)){
    # get position of cluster in full list
    position.in.cluster.list <- rownames(flow.data.cluster.median)[
      grepl(cluster.group, rownames(flow.data.cluster.median))][cluster]
    
    # append markers to name for each positive variable up to the number of clusters in the group or max 4
    
    markers.to.append.n <- ifelse( nrow(data.temp) < 5, nrow(data.temp), 4 )
    
    for( marker in names(sorted.variances)[1:markers.to.append.n]){
      rownames(flow.data.cluster.median)[as.numeric(names(position.in.cluster.list))] <- 
        if ( abs(flow.data.cluster.median[as.numeric(names(position.in.cluster.list)),marker] - max.expression[marker] ) < abs(flow.data.cluster.median[as.numeric(names(position.in.cluster.list)),marker] - min.expression[marker] )){
          paste( rownames(flow.data.cluster.median)[as.numeric(names(position.in.cluster.list))], marker )
        } else{
          rownames(flow.data.cluster.median)[as.numeric(names(position.in.cluster.list))] <- rownames(flow.data.cluster.median)[as.numeric(names(position.in.cluster.list))]
        }
    }
    
  }
  
}

# set names as fcs.cluster.label
fcs.cluster.label <- rownames(flow.data.cluster.median)


cat("
    Plotting histograms for each marker for samples...\n")

source( file.path( fcs.src.dir, "plot_data_histograms.r") )

cat("
    Plotting histograms for each marker for clusters...\n")

source( file.path( fcs.src.dir, "plot_cluster_histograms.r") ) 

cat("\n")


# give option to rename clusters--------

happy.with.cluster.names <- 0

while( happy.with.cluster.names !=1){
  clusters.to.rename <- multi.menu( fcs.cluster.label, title = "Do you want to rename any clusters? ")
  clusters.to.rename <- as.numeric( clusters.to.rename )
  
  if (length(clusters.to.rename) > 1){
    for(cluster in clusters.to.rename){
      fcs.cluster.label[cluster] <- readline(paste0("Please enter a new name for the cluster currently called ", fcs.cluster.label[cluster], ": ") )
    }
    
    cat("\n")
    
    cat(
      "Your clusters will be named as follows: 
      \n")
    
    cat( paste0(fcs.cluster.label, collapse = "\n"))
    cat( "\n" )
    cat( "\n" )
    
    happy.with.cluster.names <- menu(c("Yes", "No"), title = "Are you happy with the cluster names? ")
    
  } else if (clusters.to.rename == 0){
      happy.with.cluster.names <- 1
      
      
      cat("\n")
      
      cat(
        "Your clusters will be named as follows: 
      \n")
      
      cat( paste0(fcs.cluster.label, collapse = "\n"))
      
    } else{
      for(cluster in clusters.to.rename){
      fcs.cluster.label[cluster] <- readline(paste0("Please enter a new name for the cluster currently called ", fcs.cluster.label[cluster], ": ") )
      }
      
      
      cat("\n")
      
      cat(
        "Your clusters will be named as follows: 
        \n")

      cat( paste0(fcs.cluster.label, collapse = "\n"))
      cat( "\n" )
      cat( "\n" )

      happy.with.cluster.names <- menu(c("Yes", "No"), title = "Are you happy with the cluster names? ")
    }

}



# save cluster counts---------------
cat("
    Exporting cluster counts and percentages as spreadsheets...\n")

flow.cluster.count <- sapply( fcs.cluster, function( fc )
  table( dmrd.event.sample[ dmrd.event.cluster == fc ] ) )

stopifnot( rownames( flow.cluster.count ) == flow.sample )
stopifnot( colnames( flow.cluster.count ) == fcs.cluster )

rownames( flow.cluster.count ) <- flow.sample.label
colnames( flow.cluster.count ) <- fcs.cluster.label

write.csv( flow.cluster.count, file = file.path( fcs.cluster.table.dir, 
                                                 sprintf( "%s.csv", fcs.cluster.table.counts ) ) )

flow.cluster.percent <- flow.cluster.count/rowSums(flow.cluster.count)*100

write.csv( flow.cluster.percent, file = file.path( fcs.cluster.table.dir, 
                                                   sprintf( "%s.csv", "cluster_proportions" ) ) )

# plot density distributions of clusters with new names---------------

# plot density distributions by cluster

if( length(clusters.to.rename) > 1){
  cat("
    Plotting histograms for each marker for clusters...\n")
  
  source( file.path( fcs.src.dir, "plot_cluster_histograms.r") ) 
} else if (clusters.to.rename != 0){
  cat("
    Plotting histograms for each marker for clusters...\n")
  
  source( file.path( fcs.src.dir, "plot_cluster_histograms.r") ) 
}


setup.end.time <- Sys.time()
setup.time <- round(difftime(setup.end.time, setup.start.time, units='mins'), 2)
setup.time <- as.numeric(setup.time)

cat("\n")

cat(
"Proceed to the next section."
)
