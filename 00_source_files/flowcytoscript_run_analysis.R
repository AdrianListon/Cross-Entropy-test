#Run the analysis

analysis.start.time <- Sys.time()

# plot heatmaps---------------
cat("Plotting heatmaps\n")

heatmap.type <- c( "by_condition", "by_sample", "by_cluster" )

for ( ht in heatmap.type )
{
  if ( ht == "by_condition" ) {
    flow.data.group.median <- apply( dmrd.data, 2, tapply, 
                                     dmrd.event.condition, median )
    margin.col <- fcs.heatmap.label.factor.col * 
      max( nchar( fcs.condition.label ) )
    group.label <- fcs.condition.label
    group.color <- fcs.condition.color
  }
  else if ( ht == "by_sample" ) {
    flow.data.group.median <- apply( dmrd.data, 2, tapply, 
                                     dmrd.event.sample, median )
    margin.col <- fcs.heatmap.label.factor.col * 
      max( nchar( flow.sample.label ) )
    group.label <- flow.sample.label
    group.color <- flow.sample.color
  }
  else if ( ht == "by_cluster" ) {
    flow.data.group.median <- apply( dmrd.data, 2, tapply, 
                                     dmrd.event.cluster, median )
    margin.col <- fcs.heatmap.label.factor.col * 
      max( nchar( fcs.cluster.label ) )
    group.label <- fcs.cluster.label
    group.color <- fcs.cluster.color
  }
  else
    stop( "wrong heatmap type" )
  
  margin.row <- 
    fcs.heatmap.label.factor.row * max( nchar( fcs.channel.label ) )
  
  if ( margin.row < 5 )
    margin.row <- 5
  if ( margin.col < 5 )
    margin.col <- 5
  
  png( filename = file.path( fcs.heatmap.figure.dir, 
                             sprintf( "%s_%s.png", fcs.heatmap.figure, ht ) ), 
       width = fcs.heatmap.width, height = fcs.heatmap.height )
  heatmap( t( flow.data.group.median ), scale = "row", 
           labRow = fcs.channel.label,  labCol = group.label, 
           col = fcs.heatmap.palette, ColSideColors = group.color, 
           cexRow = fcs.heatmap.font.size, cexCol = fcs.heatmap.font.size, 
           margins = c( margin.col, margin.row ) )
  dev.off()
}


# plot histograms---------------
cat("Plotting histograms of cluster frequency\n")

histogram.type <- c( "by_cluster", "by_sample", "by_assay" )

condition.sample.n <- as.vector( table( flow.sample.condition ) )

dmrd.event.sample.n <- as.numeric( table( dmrd.event.sample ) )

for ( ht in histogram.type )
{
  flow.data.cluster.fraction <- lapply( fcs.cluster, function( fc ) {
    sample.event.n <- as.vector( table( 
      dmrd.event.sample[ dmrd.event.cluster == fc ] ) )
    
    if ( ht == "by_cluster" )
    {
      sample.fraction <- log2( 
        ( sample.event.n / dmrd.event.cluster.n[ fc ] ) / 
          ( dmrd.event.sample.n / dmrd.event.n ) 
      )
      sample.fraction[ is.infinite( sample.fraction ) ] <- NA
      fraction <- tapply( sample.fraction, flow.sample.condition, mean, 
                          na.rm = TRUE )
      std.err <- tapply( sample.fraction, flow.sample.condition, sd, 
                         na.rm = TRUE ) / sqrt( condition.sample.n )
    }
    else if ( ht == "by_sample" )
    {
      sample.fraction <- 100 * sample.event.n / dmrd.event.sample.n
      fraction <- tapply( sample.fraction, flow.sample.condition, mean, 
                          na.rm = TRUE )
      std.err <- tapply( sample.fraction, flow.sample.condition, sd, 
                         na.rm = TRUE ) / sqrt( condition.sample.n )
    }
    else if ( ht == "by_assay" )
    {
      sample.fraction <- 100 * sample.event.n / dmrd.event.n
      fraction <- tapply( sample.fraction, flow.sample.condition, sum, 
                          na.rm = TRUE )
      std.err <- NA
    }
    else
      stop( "wrong histogram type" )
    
    data.frame( cluster = fc, condition = fcs.condition, fraction, std.err )
  } )
  
  flow.data.cluster.fraction <- do.call( rbind, flow.data.cluster.fraction )
  flow.data.cluster.fraction$condition <- factor( 
    flow.data.cluster.fraction$condition, levels = fcs.condition )
  
  histogram.plot <- ggplot( flow.data.cluster.fraction, 
                            aes( x = cluster, y = fraction, fill = condition ) ) + 
    geom_bar( stat = "identity", position = position_dodge2() ) + 
    geom_errorbar( aes( ymin = fraction - std.err, 
                        ymax = fraction + std.err ), 
                   linewidth = fcs.histogram.error.bar.size,
                   position = position_dodge2() ) + 
    scale_x_discrete( labels = fcs.cluster.label, name = "" ) + 
    scale_fill_manual( values = fcs.condition.color, 
                       limits = fcs.condition, breaks = fcs.condition, 
                       labels = fcs.condition.label ) + 
    theme_bw() + 
    theme( 
      axis.text = element_text( size = fcs.histogram.font.size, 
                                angle = 90 ), 
      axis.title = element_text( size = fcs.histogram.font.size + 1 ), 
      legend.title = element_blank(), 
      legend.text = element_text( size = fcs.histogram.font.size + 1 ), 
      legend.key.size = unit( fcs.histogram.legend.key.size, "lines" ), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank() 
    )
  
  if ( ht == "by_cluster" )
    histogram.plot <- histogram.plot + 
    ylab( "Average log2( frequency ratio )" )
  else if ( ht == "by_sample" )
    histogram.plot <- histogram.plot + ylab( "Average frequency" )
  else if ( ht == "by_assay" )
    histogram.plot <- histogram.plot + ylab( "Frequency" )
  
  ggsave( 
    file.path( fcs.histogram.figure.dir, 
               sprintf( "%s_%s.png", fcs.histogram.figure, ht ) ), 
    histogram.plot, 
    width = fcs.histogram.width, 
    height = fcs.histogram.height + 
      fcs.histogram.label.factor.height * 
      max( nchar( fcs.cluster.label ) ) 
  )
}


# run pca on median values for each marker for each sample------------
cat("Running PCA on MFIs\n")

flow.data.group.median <- data.frame(apply( dmrd.data, 2, tapply, 
                                            dmrd.event.sample, median ))

pca.mfi.result <- prcomp(flow.data.group.median)

# add labels and extract coordinates
flow.data.group.median$Sample <- rownames(flow.data.group.median)
flow.data.group.median$sample.color <- flow.sample.color[flow.data.group.median$Sample]
flow.data.group.median$Group <- flow.sample.condition[flow.data.group.median$Sample]
flow.data.group.median[ , c("PC1", "PC2")] <- pca.mfi.result$x[, 1:2]
pca.mfi.loadings <- data.frame(pca.mfi.result$rotation)
pca.mfi.loadings$var <- rownames(pca.mfi.result$rotation)
mfi.scaling.loadings <- max(abs( flow.data.group.median[,c("PC1", "PC2")] )) / max( abs(pca.mfi.loadings[, 1:2])) / 2
pca.mfi.loadings[, 1:2] <- pca.mfi.loadings[, 1:2] * mfi.scaling.loadings
pc.mfi.variance <- summary(pca.mfi.result)
pc.mfi.variance <- data.frame(pc.mfi.variance$importance)
pc1.mfi.proportion <- pc.mfi.variance[2,1]*100
pc2.mfi.proportion <- pc.mfi.variance[2,2]*100

names(new.group.labels) <- fcs.condition

# plot pca on markers with loadings
pca.marker.var.plot <- ggplot(data = flow.data.group.median, 
                              aes(PC1, PC2, col = Group )) +
  geom_point( size = pca.point.size ) +
  scale_color_manual( values = fcs.condition.color,
                      breaks = fcs.condition,
                      labels = new.group.labels)+
  geom_segment( inherit.aes = FALSE, data = pca.mfi.loadings, 
                aes(x = 0, y = 0, xend = PC1, yend = PC2 ),
                color = "black", linewidth = pca.loading.arrow.size, 
                arrow = arrow(length = unit(0.03, "npc")) ) +
  geom_label( inherit.aes = FALSE, data = pca.mfi.loadings, 
              aes(PC1 * 1.2, PC2 * 1.2, label = pca.mfi.loadings$var ), 
              size = pca.loading.text.size,
              label.size = pca.loading.label.size ) +
  xlab(paste( "PC1 (", pc1.mfi.proportion, "%)", sep = "" )) +
  ylab(paste( "PC2 (", pc2.mfi.proportion, "%)", sep = "" )) +
  theme_bw() +
  theme( panel.grid = element_blank(), legend.title = element_blank() )

ggsave(filename = file.path( fcs.pca.figure.dir, 
                             paste0( pca.figure.mfi, "_loadings.png" )),
       pca.marker.var.plot, width = fcs.pca.figure.width, height = fcs.pca.figure.height )

# plot pca on markers without loadings
pca.marker.sample.plot <- ggplot(data = flow.data.group.median, 
                                 aes(PC1, PC2, col = Group)) +
  geom_point(size = pca.point.size ) +
  scale_color_manual( values = fcs.condition.color,
                      breaks = fcs.condition,
                      labels = new.group.labels)+
  xlab(paste( "PC1 (", pc1.mfi.proportion, "%)", sep = "" )) +
  ylab(paste( "PC2 (", pc2.mfi.proportion, "%)", sep = "" ))+
  theme_bw() +
  theme(panel.grid = element_blank(), legend.title = element_blank())

ggsave(filename = file.path( fcs.pca.figure.dir, 
                             paste( pca.figure.mfi, "_samples.png" )),
       pca.marker.sample.plot, width = fcs.pca.figure.width, height = fcs.pca.figure.height )


# run pca on cluster frequencies------------
cat("Running PCA on cluster distributions\n")

pca.cluster.result <- prcomp(flow.cluster.percent)

# add labels and extract coordinates

flow.cluster.frame <- data.frame(flow.cluster.percent)

flow.cluster.frame$Sample <- names(flow.sample.color)

flow.cluster.frame$sample.color <- flow.sample.color[flow.cluster.frame$Sample]

flow.cluster.frame$Group <- flow.data.group.median$Group

flow.cluster.frame[ , c("PC1", "PC2")] <- pca.cluster.result$x[, 1:2]
pca.cluster.loadings <- data.frame(pca.cluster.result$rotation)
pca.cluster.loadings$var <- rownames(pca.cluster.result$rotation)
cluster.scaling.loadings <- max(abs( flow.cluster.frame[,c("PC1", "PC2")] )) / max( abs(pca.cluster.loadings[, 1:2])) / 2
pca.cluster.loadings[, 1:2] <- pca.cluster.loadings[, 1:2] * cluster.scaling.loadings
pc.cluster.variance <- summary(pca.cluster.result)
pc.cluster.variance <- data.frame(pc.cluster.variance$importance)
pc1.cluster.proportion <- pc.cluster.variance[2,1]*100
pc2.cluster.proportion <- pc.cluster.variance[2,2]*100

# plot pca on clusters with loadings
pca.cluster.var.plot <- ggplot(data = flow.cluster.frame, 
                              aes(PC1, PC2, col = Group )) +
  geom_point( size = pca.point.size ) +
  scale_color_manual( values = fcs.condition.color,
                      breaks = fcs.condition,
                      labels = new.group.labels)+
  geom_segment( inherit.aes = FALSE, data = pca.cluster.loadings, 
                aes(x = 0, y = 0, xend = PC1, yend = PC2 ),
                color = "black", linewidth = pca.loading.arrow.size, 
                arrow = arrow(length = unit(0.03, "npc")) ) +
  geom_label( inherit.aes = FALSE, data = pca.cluster.loadings, 
              aes(PC1 * 1.2, PC2 * 1.2, label = pca.cluster.loadings$var ), 
              size = pca.loading.text.size,
              label.size = pca.loading.label.size ) +
  xlab(paste( "PC1 (", pc1.cluster.proportion, "%)", sep = "" )) +
  ylab(paste( "PC2 (", pc2.cluster.proportion, "%)", sep = "" )) +
  theme_bw() +
  theme( panel.grid = element_blank(), legend.title = element_blank() )

ggsave(filename = file.path( fcs.pca.figure.dir, 
                             paste0( pca.figure.cluster, "_loadings.png" )),
       pca.cluster.var.plot, width = fcs.pca.figure.width, height = fcs.pca.figure.height )


# plot pca on clusters without loadings
pca.cluster.sample.plot <- ggplot(data = flow.cluster.frame, 
                                  aes(PC1, PC2, col = Group )) +
  geom_point(size = pca.point.size ) +
  scale_color_manual( values = fcs.condition.color,
                      breaks = fcs.condition,
                      labels = new.group.labels)+
  xlab(paste( "PC1 (", pc1.cluster.proportion, "%)", sep = "" )) +
  ylab(paste( "PC2 (", pc2.cluster.proportion, "%)", sep = "" ))+
  theme_bw() +
  theme(panel.grid = element_blank(), legend.title = element_blank())

ggsave(filename = file.path( fcs.pca.figure.dir, 
                             paste( pca.figure.cluster, "_samples.png" )),
       pca.cluster.sample.plot, width = fcs.pca.figure.width, height = fcs.pca.figure.height )


# run ANOVA with Tukey's post-test on marker MFIs and cluster frequencies-------

source( file.path( fcs.src.dir, "flowcytoscript_anova.r") )


# calculate tsne representation---------------

cat("\n")

{
  if ( fcs.use.cached.results && file.exists( fcs.tsne.cache.file.path ) )
  {
    cat( "Using cached results for tsne\n" )
    
    load( fcs.tsne.cache.file.path )
  }
  else
  {
    fcs.tsne.learning.rate <- ifelse( nrow(dmrd.data)/fcs.tsne.exaggeration.factor > 2000, 
                                      nrow(dmrd.data)/fcs.tsne.exaggeration.factor, 2000 )
    
    set.seed.here( fcs.seed.base, "calculate tsne representation" )
    
    cat( "Calculating tSNE\n" )
    
    tsne.result <- Rtsne_neighbors( dmrd.knn$idx, dmrd.knn$dist, 
                                    perplexity = fcs.tsne.perplexity, 
                                    exaggeration_factor = fcs.tsne.exaggeration.factor,
                                    eta = fcs.tsne.learning.rate,
                                    max_iter = fcs.tsne.iter.n, stop_lying_iter = fcs.tsne.early.iter.n,
                                    check_duplicates = FALSE, pca = FALSE, 
                                    num_threads = fcs.tsne.threads.n )
    
    save( tsne.result, file = fcs.tsne.cache.file.path )
  }
}

tsne.data <- tsne.result$Y


# plot tsne convergence---------------
cat( "Plotting tSNE convergence\n" )

tsne.iter <- 1 + 50 * ( 1 : length( tsne.result$itercosts ) - 1 )
tsne.cost <- tsne.result$itercosts

png( filename = file.path( fcs.tsne.figure.dir, 
                           sprintf( "%s.png", fcs.tsne.figure.convergence ) ), 
     width = fcs.tsne.figure.convergence.width, 
     height = fcs.tsne.figure.convergence.height )
par( mar = c( 5, 5.8, 2, 1.4 ) )
plot( tsne.iter, tsne.cost, log = "y", xlab = "Iteration", 
      ylab = "Total cost", xlim = c( 0, fcs.tsne.iter.n ), 
      ylim = c( 1, max( tsne.cost ) ), col = "blue3", 
      pch = 20, cex = 1, cex.lab = 2.5, cex.axis = 2 )
abline( h = tsne.cost[ length( tsne.cost ) ], lty = 2, lwd = 1.5, 
        col = "blue3" )
dev.off()


# plot tsne representations---------------
cat( "Plotting tSNE representations\n" )

tsne.data.max <- max( abs( tsne.data ) )

# plot tSNE with clusters colored

plot.cluster.drmd.figures(
  tsne.data, tsne.data.max, 
  fcs.tsne.figure.lims.factor, fcs.tsne.figure.point.size, 
  fcs.tsne.figure.dir, fcs.tsne.figure.plot, 
  dmrd.data, dmrd.event.cluster, dmrd.event.condition
)

# plot tSNE with everything

if (plot.everything == 1){
  plot.all.dmrd.figures( 
    tsne.data, tsne.data.max, 
    fcs.tsne.figure.lims.factor, fcs.tsne.figure.point.size, 
    fcs.tsne.figure.dir, fcs.tsne.figure.plot, 
    dmrd.data, dmrd.event.cluster, dmrd.event.condition
  )
}



# calculate umap representation---------------

{
  if ( fcs.use.cached.results && file.exists( fcs.umap.cache.file.path ) )
  {
    cat( "Using cached results for umap\n" )
    
    load( fcs.umap.cache.file.path )
  }
  else
  {
    set.seed.here( fcs.seed.base, "calculate umap representation" )
    
    cat( "Calculating UMAP\n" )
    
    umap.result <- uwot::umap(dmrd.data, n_epochs = fcs.umap.iter.n, 
                              nn_method = dmrd.knn, n_threads = fcs.tsne.threads.n,
                              n_sgd_threads = fcs.tsne.threads.n, 
                              batch = TRUE, verbose = TRUE, ret_model = TRUE,
                              ret_nn = TRUE, ret_extra = "sigma")
    
    save(umap.result, file = fcs.umap.cache.file.path )
    
  }
}

umap.data <- umap.result$embedding
dimnames( umap.data ) <- NULL


# plot umap representations---------------
cat( "Plotting UMAP representations\n" )

#umap.data.max <- max( abs( apply( umap.data, 2, function( x ) 
#  quantile( x, c( 0.25, 0.75 ) ) + c( -1, 1 ) * IQR( x ) ) ) )

umap.data.max <- max( abs(umap.data) )

# plot umap with clusters colored

plot.cluster.drmd.figures(
  umap.data, umap.data.max, 
  fcs.umap.figure.lims.factor, fcs.umap.figure.point.size, 
  fcs.umap.figure.dir, fcs.umap.figure.plot, 
  dmrd.data, dmrd.event.cluster, dmrd.event.condition 
)

# plot umap with everything

if(plot.everything == 1){
  plot.all.dmrd.figures( 
    umap.data, umap.data.max, 
    fcs.umap.figure.lims.factor, fcs.umap.figure.point.size, 
    fcs.umap.figure.dir, fcs.umap.figure.plot, 
    dmrd.data, dmrd.event.cluster, dmrd.event.condition 
  )
}



# trex: highlight changed regions---------------
cat( "Finding regions of over or under-representation (T-REX)\n" )

# subset data based on trex groups
trex.dmrd.data <- subset(dmrd.data, dmrd.event.condition %in% trex.condition)
trex.umap.data <- umap.data
rownames(trex.umap.data) <- rownames(dmrd.data)
trex.umap.data <- subset(trex.umap.data, dmrd.event.condition %in% trex.condition)
trex.tsne.data <- tsne.data
rownames(trex.tsne.data) <- rownames(dmrd.data)
trex.tsne.data <- subset(trex.tsne.data, dmrd.event.condition %in% trex.condition)

# calculate knn for the trex groups


if ( length(fcs.condition) != length(trex.condition) ) {
  non.self.knn <- 1
  knn.M <- knn.M/2
  while( length(non.self.knn) > 0) {
    print(knn.M)
    set.seed.here( fcs.seed.base, "get knn with Rcpp" )
    trex.knn <- RcppHNSW::hnsw_knn( trex.dmrd.data, k= fcs.tsne.perplexity,
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
} else {
  trex.knn <- dmrd.knn
}

trex.knn.index <- trex.knn$idx
first.condition.length <- nrow(subset(dmrd.data, dmrd.event.condition == trex.condition[1]))
trex.knn.index[trex.knn.index <= first.condition.length] <- 0
trex.knn.index[trex.knn.index > first.condition.length] <- 1

# calculate percent change in each KNN region
percent.change.knn <- (rowSums(trex.knn.index) / fcs.tsne.perplexity * 100)

# create trex plots
cat("Plotting changed regions on tSNE and UMAP\n")

plot.changed.regions(
  trex.tsne.data, percent.change.knn,
  trex.figure.point.size, trex.condition, 
  trex.figure.dir, fcs.tsne.figure.plot
)

plot.changed.regions(
  trex.umap.data, percent.change.knn,
  trex.figure.point.size, trex.condition,
  trex.figure.dir, fcs.umap.figure.plot
)


# run crossentropy if desired-----------

if(run.crossentropy == 1){
  # calculate cross-entropy test for tsne by condition---------------
  
  set.seed.here( fcs.seed.base, "calculate cross-entropy test for tsne" )
  
  ce.diff.test.tsne.res <- ce.diff.test.tsne( 
    dmrd.data, tsne.data, 
    dmrd.event.condition, 
    partition.label = fcs.condition.label, 
    partition.color = fcs.condition.color, 
    partition.line.type = fcs.condition.line.type, 
    base.test = fcs.ce.diff.base.test, 
    base.dist = fcs.ce.diff.base.dist, 
    prob.sample.n = fcs.ce.diff.prob.sample.n, 
    dendrogram.order.weight = fcs.ce.diff.figure.dendrogram.weight.condition, 
    result = file.path( fcs.ce.diff.tsne.figure.dir, 
                        sprintf( "%s_condition.txt", fcs.ce.diff.tsne.result ) ), 
    cdf.figure = file.path( fcs.ce.diff.tsne.figure.dir, 
                            sprintf( "%s_condition.png", fcs.ce.diff.tsne.figure.cdf ) ), 
    dendrogram.figure = file.path( fcs.ce.diff.tsne.figure.dir, 
                                   sprintf( "%s_condition.png", fcs.ce.diff.tsne.figure.dendrogram ) ) )
  
  
  # calculate cross-entropy test for tsne by sample---------------
  
  set.seed.here( fcs.seed.base, "calculate cross-entropy test for tsne" )
  
  ce.diff.test.tsne.res <- ce.diff.test.tsne( 
    dmrd.data, tsne.data, 
    dmrd.event.sample, 
    partition.label = flow.sample.label, 
    partition.color = flow.sample.color, 
    partition.line.type = flow.sample.line.type.single, 
    base.test = fcs.ce.diff.base.test, 
    base.dist = fcs.ce.diff.base.dist, 
    prob.sample.n = fcs.ce.diff.prob.sample.n, 
    dendrogram.order.weight = flow.ce.diff.figure.dendrogram.weight.sample, 
    result = file.path( fcs.ce.diff.tsne.figure.dir, 
                        sprintf( "%s_sample.txt", fcs.ce.diff.tsne.result ) ), 
    cdf.figure = file.path( fcs.ce.diff.tsne.figure.dir, 
                            sprintf( "%s_sample.png", fcs.ce.diff.tsne.figure.cdf ) ), 
    dendrogram.figure = file.path( fcs.ce.diff.tsne.figure.dir, 
                                   sprintf( "%s_sample.png", fcs.ce.diff.tsne.figure.dendrogram ) ) )
  
  
  # calculate cross-entropy test for tsne by cluster---------------
  
  fcs.ce.diff.figure.dendrogram.weight.cluster <- 1 : fcs.cluster.n
  names( fcs.ce.diff.figure.dendrogram.weight.cluster ) <- fcs.cluster
  
  set.seed.here( fcs.seed.base, "calculate cross-entropy test for tsne" )
  
  ce.diff.test.tsne.res <- ce.diff.test.tsne( 
    dmrd.data, tsne.data, 
    dmrd.event.cluster, 
    partition.label = fcs.cluster.label, 
    partition.color = fcs.cluster.color, 
    partition.line.type = fcs.cluster.line.type, 
    base.test = fcs.ce.diff.base.test, 
    base.dist = fcs.ce.diff.base.dist, 
    prob.sample.n = fcs.ce.diff.prob.sample.n, 
    dendrogram.order.weight = fcs.ce.diff.figure.dendrogram.weight.cluster, 
    result = file.path( fcs.ce.diff.tsne.figure.dir, 
                        sprintf( "%s_cluster.txt", fcs.ce.diff.tsne.result ) ), 
    cdf.figure = file.path( fcs.ce.diff.tsne.figure.dir, 
                            sprintf( "%s_cluster.png", fcs.ce.diff.tsne.figure.cdf ) ), 
    dendrogram.figure = file.path( fcs.ce.diff.tsne.figure.dir, 
                                   sprintf( "%s_cluster.png", fcs.ce.diff.tsne.figure.dendrogram ) ) )
  
  
  # calculate cross-entropy test for umap by condition---------------
  
  set.seed.here( fcs.seed.base, "calculate cross-entropy test for umap" )
  
  ce.diff.test.umap.res <- ce.diff.test.umap( 
    umap.result$nn$precomputed$dist, umap.result$nn$precomputed$idx, 
    umap.data, umap.result, 
    dmrd.event.condition, 
    partition.label = fcs.condition.label, 
    partition.color = fcs.condition.color, 
    partition.line.type = fcs.condition.line.type, 
    base.test = fcs.ce.diff.base.test, 
    base.dist = fcs.ce.diff.base.dist, 
    prob.sample.n = fcs.ce.diff.prob.sample.n, 
    dendrogram.order.weight = fcs.ce.diff.figure.dendrogram.weight.condition, 
    result = file.path( fcs.ce.diff.umap.figure.dir, 
                        sprintf( "%s_condition.txt", fcs.ce.diff.umap.result ) ), 
    cdf.figure = file.path( fcs.ce.diff.umap.figure.dir, 
                            sprintf( "%s_condition.png", fcs.ce.diff.umap.figure.cdf ) ), 
    dendrogram.figure = file.path( fcs.ce.diff.umap.figure.dir, 
                                   sprintf( "%s_condition.png", fcs.ce.diff.umap.figure.dendrogram ) ) )
  
  
  # calculate cross-entropy test for umap by sample---------------
  
  set.seed.here( fcs.seed.base, "calculate cross-entropy test for umap" )
  
  ce.diff.test.umap.res <- ce.diff.test.umap( 
    umap.result$nn$precomputed$dist, umap.result$nn$precomputed$idx, 
    umap.data, umap.result, 
    dmrd.event.sample, 
    partition.label = flow.sample.label, 
    partition.color = flow.sample.color, 
    partition.line.type = flow.sample.line.type.single, 
    base.test = fcs.ce.diff.base.test, 
    base.dist = fcs.ce.diff.base.dist, 
    prob.sample.n = fcs.ce.diff.prob.sample.n, 
    dendrogram.order.weight = flow.ce.diff.figure.dendrogram.weight.sample, 
    result = file.path( fcs.ce.diff.umap.figure.dir, 
                        sprintf( "%s_sample.txt", fcs.ce.diff.umap.result ) ), 
    cdf.figure = file.path( fcs.ce.diff.umap.figure.dir, 
                            sprintf( "%s_sample.png", fcs.ce.diff.umap.figure.cdf ) ), 
    dendrogram.figure = file.path( fcs.ce.diff.umap.figure.dir, 
                                   sprintf( "%s_sample.png", fcs.ce.diff.umap.figure.dendrogram ) ) )
  
  
  # calculate cross-entropy test for umap by cluster---------------
  
  set.seed.here( fcs.seed.base, "calculate cross-entropy test for umap" )
  
  ce.diff.test.umap.res <- ce.diff.test.umap( 
    umap.result$nn$precomputed$dist, umap.result$nn$precomputed$idx, 
    umap.data, umap.result, 
    dmrd.event.cluster, 
    partition.label = fcs.cluster.label, 
    partition.color = fcs.cluster.color, 
    partition.line.type = fcs.cluster.line.type, 
    base.test = fcs.ce.diff.base.test, 
    base.dist = fcs.ce.diff.base.dist, 
    prob.sample.n = fcs.ce.diff.prob.sample.n, 
    dendrogram.order.weight = fcs.ce.diff.figure.dendrogram.weight.cluster, 
    result = file.path( fcs.ce.diff.umap.figure.dir, 
                        sprintf( "%s_cluster.txt", fcs.ce.diff.umap.result ) ), 
    cdf.figure = file.path( fcs.ce.diff.umap.figure.dir, 
                            sprintf( "%s_cluster.png", fcs.ce.diff.umap.figure.cdf ) ), 
    dendrogram.figure = file.path( fcs.ce.diff.umap.figure.dir, 
                                   sprintf( "%s_cluster.png", fcs.ce.diff.umap.figure.dendrogram ) ) )
  
}

# calculate stats for T-REX groups-------------

trex.group.cluster.long <- flow.cluster.frame.long %>%
  dplyr::filter( Group %in% trex.condition)

p.values.all <- NULL
p.values.all <- list()

for (cluster in unique(trex.group.cluster.long$Cluster)) {
  temp <- trex.group.cluster.long %>% 
    dplyr::filter(Cluster == cluster) %>%
    group_by(Group) %>%
    suppressMessages(summarise(Percent = Percent))
  
  res_aov <- aov( Percent ~ Group, data = temp )
  
  fitted.em <- emmeans(res_aov, "Group", data = temp )
  
  p.values <- data.frame( pairs(fitted.em, adjust = "tukey" ) )
  p.values <- p.values[,c(1,6)]
  p.values$significant <- ifelse( p.values$p.value < 0.05, "Yes", "No" )
  p.values.all[cluster] <- p.values$p.value
  
  write.csv(p.values, file.path( trex.figure.dir, paste0( cluster, "_frequencies_anova.csv")))
  
}

#names(p.values.all) <- fcs.cluster.label
p.values.all.available <- p.values.all[!is.na(p.values.all)]
p.values.all.significant <- p.values.all.available[p.values.all.available < 0.05]

# format for report
if( length(tissue.type) >1 ){
  tissue.type <- tissue.type
} else if ( tissue.type == "Immune"){
  tissue.type <- "Immune system"
}

if( length(p.values.all.significant) !=0 ){
  p.value.message <- length( p.values.all.significant )
} else {
  p.value.message <- "none"
}

analysis.end.time <- Sys.time()
analysis.calc.time <- round(difftime(analysis.end.time, analysis.start.time, units='mins'), 2)
analysis.calc.time <- as.numeric(analysis.calc.time)


## inform user analysis is finished---------
cat(
  "\n
 Analysis complete!\n
  \n"
)


