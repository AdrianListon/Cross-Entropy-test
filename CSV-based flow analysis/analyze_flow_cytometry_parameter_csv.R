# Analysis of 20210608 tSNE-diff data


# data parameters---------------

# Tip: create two folders inside your analysis directory called "Data" and "Output".
# Put the fcs files in the Data folder.
# Put the parameter file in the Output folder. Double click on it to start Rstudio.
# Opening in this way will set your working directory to the Output folder.
# Type getwd(). Copy the result into the line below, swapping "Data" for "Output"
fcs.data.dir <- "D:/tSNE-diff paper/Analysis/20210608 tSNE-diff flow/Data/csv 20k events"

fcs.condition <- c( "LN", "Spleen",
                    "siLPL" )

fcs.condition.label <- c( 
    "LN" = "LN", "Spleen" = "Spleen",
    "siLPL" = "siLPL"
)

# Run the analyze_flow_cytometry_csv script through line 50 and insert the output below.
# Remove unwanted channels.

fcs.channel <- c( 
  "RORgT","CD44","Gr.1","IgM","F4_80","Ki67","CD90.2","CCR9",
  "TCRgd","PDCA.1","CD11c","Ly.6C","CD103","IgD","NK1.1",
  "CTLA.4","c.Kit","CD62L","GITR","CD150","CXCR3","Siglec.F",
  "TCRb","PD.1","XCR1","CD127","CCR2","CD45","CD4","CD8","CD3",
  "CD19","CD11b","CD38","GATA.3","CD86","Foxp3","CD172a","CD64",
  "Helios","CD24",
  "MHCII","NKp46","CD69","B220","CD25","ICOS","KLRG1","T.bet"  
)

fcs.channel.label <- c( 
  "RORgT" = "RORgT",
  "CD44" = "CD44",
  "Gr.1" = "Gr.1",
  "IgM" = "IgM",
  "F4_80" = "F4_80",
  "Ki67" = "Ki67",
  "CD90.2" = "CD90.2",
  "CCR9" = "CCR9",
  "TCRgd" = "TCRgd",
  "PDCA.1" = "PDCA.1",
  "CD11c" = "CD11c",
  "Ly.6C" = "Ly.6C",
  "CD103" = "CD103",
  "IgD" = "IgD",
  "NK1.1" = "NK1.1",
  "CTLA.4" = "CTLA.4",
  "c.Kit" = "c.Kit",
  "CD62L" = "CD62L",
  "GITR" = "GITR",
  "CD150" = "CD150",
  "CXCR3" = "CXCR3",
  "Siglec.F" = "Siglec.F",
  "TCRb" = "TCRb",
  "PD.1" = "PD.1",
  "XCR1" = "XCR1",
  "CD127" = "CD127",
  "CCR2" = "CCR2",
  "CD45" = "CD45",
  "CD4" = "CD4",
  "CD8" = "CD8",
  "CD3" = "CD3",
  "CD19" = "CD19",
  "CD11b" = "CD11b",
  "CD38" = "CD38",
  "GATA.3" = "GATA.3",
  "CD86" = "CD86",
  "Foxp3" = "Foxp3",
  "CD172a" = "CD172a",
  "CD64" = "CD64",
  "Helios" = "Helios",
  "CD24" = "CD24",
  "MHCII" = "MHCII",
  "NKp46" = "NKp46",
  "CD69" = "CD69",
  "B220" = "B220",
  "CD25" = "CD25",
  "ICOS" = "ICOS",
  "KLRG1" = "KLRG1",
  "T.bet" = "T.bet"
)

fcs.condition.n <- length( fcs.condition )
fcs.channel.n <- length( fcs.channel )


# general parameters---------------
# Tip: edit the directory to match the location of the script files. 
fcs.src.dir <- "D:/Flowcytoscript/Modifications/CSV version/00_src"

# Tip: use today's date
fcs.seed.base <- 20220623


# Tip: Set to FALSE while optimizing the tSNE and FlowSOM. 
# Once you have a good run, set to TRUE and then change labels, colors, etc.
fcs.use.cached.results <- TRUE

fcs.sample.number.width <- 2
fcs.event.number.width <- 6


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
# Sampling largely irrelevant for this method

fcs.density.data.sample.n <- NULL

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
# Tip: use more clusters for more diverse cell collections. Use as few as possible to speed up processing.
# For example, for a broad immune phenotyping panel use 40-50 clusters.
# For pre-gated cell types (e.g., CD8s), use 10-12 clusters.

fcs.cluster.n <- 40

# Tip: naming clusters can be done automatically as clusters 1:n via the first command below.
# To rename your clusters based on marker expression, insert a # before the first command,
# remove the # before the second command,
# and rename the clusters appropriately.

fcs.cluster <- sprintf( "%02d", 1 : fcs.cluster.n )
#fcs.cluster <- c("B cells", "CD4 T cells", "CD8 T cells", "Monocytes",
#                "Neutrophils", "NK cells", "Eosinophils",
#                "cDC", "Macrophages")

fcs.cluster.label <- fcs.cluster
names( fcs.cluster.label ) <- fcs.cluster

# Tip: this controls the grouping of the fcs files that are generated based on the flowSOM clustering.
# If some of the clusters are related (e.g., CD14+ and CD16+ monocytes), 
# you might wish to group these into a single output folder.
# If so, you may use a command such as "a" = c(1, 5, 7),
# where the numbers are the numbered flowSOM clusters.
fcs.cluster.group <- as.list(1:fcs.cluster.n)
names(fcs.cluster.group) <- make.unique(rep(letters, length.out = fcs.cluster.n), sep='')

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

# 24 is recommended for EmbedSOM
fcs.flow.som.dim <- 24

fcs.cluster.table.dir <- "./table_cluster"
fcs.cluster.data.dir <- "./data_cluster"

fcs.cluster.table.counts <- "cluster_counts"

fcs.cluster.data <- "cluster_group"


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


# dimensionality reduction parameters---------------
# Tip: Use the third option in most cases.
# More cells will take longer. 
# You might downsample to 50-100k cells at first, then run again with more cells if you like the result.

fcs.dmrd.data.sample.n <- NULL
fcs.dmrd.data.sample.n.per.condition <- NULL
fcs.dmrd.data.sample.n.per.sample <- 2000

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
fcs.dmrd.figure.nrow <- 1
fcs.dmrd.figure.ncol <- ceiling( fcs.condition.n / fcs.dmrd.figure.nrow )

# Tip: set the figure size here.
fcs.dmrd.figure.width <- 2
fcs.dmrd.figure.height <- 2


# tsne parameters---------------
# Tip: more iterations take longer.
# More iterations are needed for more cells
# For a first look, try 1000.
# For a final figure on ~100k cells, use 5000.
# Tip 2: Set fcs.tsne.thread.n to 0 for max processing speed.
# To enable you to do something else meanwhile, 
# set it to one or two less than the number of threads on your processor (check specs online).

fcs.tsne.iter.n <- 5000
fcs.tsne.thread.n <- 0

# Tip: don't change this unless you know why
fcs.tsne.perplexity <- 30
fcs.tsne.theta <- 0.5
fcs.tsne.exaggeration.factor <- 12

fcs.tsne.figure.lims.factor <- 1.0
fcs.tsne.figure.point.size <- 1.2

fcs.tsne.figure.convergence.width <- 1200
fcs.tsne.figure.convergence.height <- 800

fcs.tsne.figure.dir <- "./figure_tsne"

fcs.tsne.figure.convergence <- "tsne_convergence"
fcs.tsne.figure.plot <- "tsne_plot"

fcs.tsne.cache.file.path <- "./tsne_cache.dat"


# umap parameters---------------
# Tip: you don't generally need to change the UMAP iterations

fcs.umap.iter.n <- 1000

fcs.umap.figure.lims.factor <- 0.8
fcs.umap.figure.point.size <- 1.2

fcs.umap.figure.dir <- "./figure_umap"

fcs.umap.figure.plot <- "umap_plot"

fcs.umap.cache.file.path <- "./umap_cache.dat"

# embedsom parameters---------------
fcs.embedsom.figure.lims.factor <- 1.0
fcs.embedsom.figure.point.size <- 0.05
fcs.ce.diff.embedsom.figure.dir <- "./figure_embedsom_ce_diff"
fcs.embedsom.figure.dir <- "./figure_embedsom"
fcs.ce.diff.embedsom.result <- "embedsom_ce_diff_result"
fcs.ce.diff.embedsom.cache.file.path <- "./embedsom_ce_diff_cache.dat"
fcs.ce.diff.embedsom.figure.cdf <- "embedsom_ce_diff_cdf"
fcs.ce.diff.embedsom.figure.dendrogram <- "embedsom_ce_diff_dendrogram"

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

fcs.ce.diff.figure.dendrogram.weight.cluster <- 1 : fcs.cluster.n
names( fcs.ce.diff.figure.dendrogram.weight.cluster ) <- fcs.cluster

fcs.ce.diff.figure.cdf.width <- 1200
fcs.ce.diff.figure.cdf.height <- 800

fcs.ce.diff.figure.dendrogram.width <- 2000
fcs.ce.diff.figure.dendrogram.height <- 800


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

