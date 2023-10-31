# definitions of variables and graphics parameters


# graphics parameters---------------

fcs.color.pool <- c( 
  brewer.pal( 8, "Set1" )[ -6 ], 
  brewer.pal( 7, "Set2" )[ c( 1, 3, 6 ) ], 
  adjustcolor( brewer.pal( 8, "Set1" )[ -6 ], 
               red.f = 0.9, green.f = 0.8, blue.f = 0.7 ), 
  adjustcolor( brewer.pal( 7, "Set2" )[ c( 1, 3, 6 ) ], 
               red.f = 0.9, green.f = 0.8, blue.f = 0.7 ), 
  adjustcolor( brewer.pal( 8, "Set1" )[ -6 ], 
               red.f = 0.8, green.f = 0.6, blue.f = 0.5 ), 
  adjustcolor( brewer.pal( 7, "Set2" )[ c( 1, 3, 6 ) ], 
               red.f = 0.8, green.f = 0.6, blue.f = 0.5 ), 
  adjustcolor( brewer.pal( 8, "Set1" )[ -6 ], 
               red.f = 0.3, green.f = 0.3, blue.f = 0.3 ), 
  adjustcolor( brewer.pal( 7, "Set2" )[ c( 1, 3, 6 ) ], 
               red.f = 0.3, green.f = 0.3, blue.f = 0.3 ) )
fcs.color.pool.n <- length( fcs.color.pool )

fcs.line.type.pool <- 1:6
fcs.line.type.pool.n <- length( fcs.line.type.pool )

fcs.condition.color <- rep( 
  fcs.color.pool, 
  ceiling( fcs.condition.n / fcs.color.pool.n ) 
)[ 1 : fcs.condition.n ]
names( fcs.condition.color ) <- fcs.condition

fcs.condition.line.type <- rep( 
  fcs.line.type.pool, 
  ceiling( fcs.condition.n / fcs.line.type.pool.n ) 
)[ 1 : fcs.condition.n ]
names( fcs.condition.line.type ) <- fcs.condition


# density parameters---------------

fcs.density.data.sample.n <- 20000

fcs.density.partition.all <- "all"
fcs.density.partition.all.label <- c( "all" = "All" )
fcs.density.partition.all.color <- c( "all" = "grey" )

fcs.density.font.size <- 4.5

fcs.density.line.size <- 0.2
fcs.density.line.alpha <- 0.3

fcs.density.figure.width.base <- 0.4
fcs.density.figure.height.base <- 0.1

fcs.density.figure.dir <- "./figure_density"

fcs.density.figure.sample <- "density_sample"
fcs.density.figure.cluster <- "density_cluster"


# cluster parameters---------------

fcs.cluster.figure.dir <- "./figure_cluster"
fcs.cluster.table.dir <- "./table_cluster"
fcs.cluster.table.counts <- "cluster_counts"



# heatmap parameters---------------

fcs.heatmap.palette.n <- 100
fcs.heatmap.palette <- colorRampPalette( brewer.pal( 9, "YlOrRd" ) )( 
  fcs.heatmap.palette.n )

fcs.heatmap.font.size <- 2.5

fcs.heatmap.label.factor.row <- 1.4
fcs.heatmap.label.factor.col <- 1.4

fcs.heatmap.width <- 2000
fcs.heatmap.height <- 2000

fcs.heatmap.figure.dir <- "./figure_heatmap"

fcs.heatmap.figure <- "heatmap"


# histogram parameters---------------

fcs.histogram.font.size <- 7
fcs.histogram.error.bar.size <- 0.2
fcs.histogram.legend.key.size <- 0.7

fcs.histogram.label.factor.height <- 0.05

fcs.histogram.width <- 4.5
fcs.histogram.height <- 3

fcs.histogram.figure.dir <- "./figure_histogram"

fcs.histogram.figure <- "histogram"

# statistics parameters---------
fcs.mfi.stats.dir <- "./marker_stats/"
fcs.cluster.stats.dir <- "./cluster_stats/"


# dimensionality reduction parameters---------------
# Tip: Use the third option in most cases.
# More cells will take longer. 
# You might downsample to 50-100k cells at first, then run again with more cells if you like the result.

fcs.dmrd.data.sample.n <- NULL
fcs.dmrd.data.sample.n.per.condition <- NULL
fcs.dmrd.data.sample.n.per.sample <- 2000

# setting the appearance of the output plots
fcs.dmrd.gradient.color <- c( "black", "blue", "green", "yellow", "red" )
fcs.dmrd.gradient.palette.n <- 100
fcs.dmrd.density.palette <- colorRampPalette( fcs.dmrd.gradient.color )( 
  fcs.dmrd.gradient.palette.n )

fcs.dmrd.color.alpha <- 0.3

fcs.dmrd.group.title.size <- 8

fcs.dmrd.legend.title.size <- 7
fcs.dmrd.legend.label.size <- 7
fcs.dmrd.legend.point.size <- 3

fcs.dmrd.label.factor.width <- 0.1

# Tip: set the number of rows in the output figures here.
fcs.dmrd.figure.nrow <- ceiling( fcs.condition.n / 4 )
fcs.dmrd.figure.ncol <- ceiling( fcs.condition.n / fcs.dmrd.figure.nrow )

# Tip: set the figure size here.
fcs.dmrd.figure.width <- 3
fcs.dmrd.figure.height <- 3

# pca parameters-------------
fcs.pca.figure.dir <- "./figure_pca"
pca.point.size <- 2
fcs.pca.figure.width <- 4.5
fcs.pca.figure.height <- 3
pca.figure.mfi <- "pca_sample_mfi"
pca.figure.cluster <- "pca_cluster"
pca.figure.factor <- "pca_factor"
pca.loading.label.size <- 0.1
pca.sample.label.size <- 3
pca.loading.text.size <- 2
pca.loading.arrow.size <- 0.1

# trex parameters-------------

trex.figure.dir <- "figure_changed_regions"
trex.figure.width <- 4
trex.figure.height <- 4
trex.figure.point.size <- 0.5
trex.legend.text.size <- 2
trex.title.text.size <- 10

# trex plot color levels
trex.percent.15_85 <- "lightgray"
trex.percent.5_15 <- "lightskyblue"
trex.percent.0_5 <- "navyblue"
trex.percent.85_95 <- "lightcoral"
trex.percent.95_100 <- "darkred"

# tsne parameters---------------
# Tip: more iterations take longer.
# This version uses a higher learning rate, so fewer iterations are needed
# 1000 should be fine for most data sets, and corresponds to about 5000 with the original approach.
# 2000 should be plenty.
# Tip 2: Set fcs.tsne.threads.n to 0 for max processing speed.
# This is done via detectCores() automatically.
# To enable you to do something else meanwhile, 
# set it to one or two less than the number of threads on your processor (check specs online).

# set tSNE iter to automatic based on nrows dmrd.data
fcs.tsne.iter.n <- 750
fcs.tsne.early.iter.n <- fcs.tsne.iter.n/10
fcs.tsne.threads.n <- parallel::detectCores()

# Tip: don't change this unless you know why
# Values between 10 and 100 really don't change things much for most data sets.
# If you have a really large data set or want to preserve more global structure,
# you can try increasing the perplexity to 100 or greater.
# This will slow it down quite a bit and perhaps generate a more UMAP-like tSNE.
fcs.tsne.perplexity <- 30
fcs.tsne.exaggeration.factor <- 4
fcs.tsne.learning.rate <- 2000

fcs.tsne.figure.lims.factor <- 1.0
fcs.tsne.figure.point.size <- 0.5

fcs.tsne.figure.convergence.width <- 1200
fcs.tsne.figure.convergence.height <- 800

fcs.tsne.figure.dir <- "./figure_tsne"
fcs.tsne.cluster.figure.dir <- "./figure_tsne_cluster"
fcs.tsne.figure.convergence <- "tsne_convergence"
fcs.tsne.figure.plot <- "tsne_plot"

fcs.tsne.cache.file.path <- "./tsne_cache.dat"


# umap parameters---------------
# Tip: you don't generally need to change the UMAP iterations

fcs.umap.iter.n <- 500

fcs.umap.figure.lims.factor <- 1.1
fcs.umap.figure.point.size <- 0.3

fcs.umap.figure.dir <- "./figure_umap"
fcs.umap.cluster.figure.dir <- "./figure_umap_cluster"
fcs.umap.figure.plot <- "umap_plot"

fcs.umap.cache.file.path <- "./umap_cache.dat"


# cross-entropy test parameters---------------
# Tips: Set this depending on your RAM and number of groups.
# You won't be able to analyze more than about 100k cells unless you have >32GB RAM.
# The crossentropy test works best with at least 10k cells per group.
# Multiple hypothesis testing will greatly reduce your ability to distinguish statistical differences.


fcs.ce.diff.prob.sample.n <- 120000

# Tip: set to "ks" unless you have a good statistical reason for using rank testing.
# In that case, use "rank" and "median".

fcs.ce.diff.base.test <- "ks"
fcs.ce.diff.base.dist <- "ks"

fcs.ce.diff.test.alpha <- 0.05

fcs.ce.diff.figure.font.size <- 2
fcs.ce.diff.figure.line.width <- 3

fcs.ce.diff.figure.cdf.resolution <- 500
fcs.ce.diff.figure.cdf.all.color <- "black"
fcs.ce.diff.figure.cdf.all.label <- "All"

fcs.ce.diff.figure.dendrogram.weight.condition <- 1 : fcs.condition.n
names( fcs.ce.diff.figure.dendrogram.weight.condition ) <- fcs.condition

fcs.ce.diff.figure.cdf.width <- 1200
fcs.ce.diff.figure.cdf.height <- 800

fcs.ce.diff.figure.dendrogram.width <- 1000
fcs.ce.diff.figure.dendrogram.height <- 600


# cross-entropy test parameters for tsne--------------

fcs.ce.diff.tsne.perplexity.factor <- 3

fcs.ce.diff.tsne.figure.dir <- "./figure_tsne_ce_diff"

fcs.ce.diff.tsne.figure.cdf <- "tsne_ce_diff_cdf"
fcs.ce.diff.tsne.figure.dendrogram <- "tsne_ce_diff_dendrogram"
fcs.ce.diff.tsne.result <- "tsne_ce_diff_result"

fcs.ce.diff.tsne.cache.file.path <- "./tsne_ce_diff_cache.dat"


# cross-entropy test parameters for umap---------------

fcs.ce.diff.umap.figure.dir <- "./figure_umap_ce_diff"

fcs.ce.diff.umap.figure.cdf <- "umap_ce_diff_cdf"
fcs.ce.diff.umap.figure.dendrogram <- "umap_ce_diff_dendrogram"
fcs.ce.diff.umap.result <- "umap_ce_diff_result"

fcs.ce.diff.umap.cache.file.path <- "./umap_ce_diff_cache.dat"
