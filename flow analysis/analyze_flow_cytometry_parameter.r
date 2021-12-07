
# tSNE diff paper--20210608 tSNE-diff flow--Figure 5E SpleenPercentLymph LymphPercentSpleen


# data parameters

# Tip: create two folders inside your analysis directory called "Data" and "Output".
# Put the fcs files in the Data folder.
# Put the parameter file in the Output folder. Double click on it to start Rstudio.
# Opening in this way will set your working directory to the Output folder.
# Type getwd(). Copy the result into the line below, swapping "Data" for "Output"
fcs.data.dir <- "D://tSNE-diff paper/Analysis/20210608 tSNE-diff flow/Data/LymphPercentSpleen/"

fcs.condition <- c( "LymphpercentSpleen", "Sp_", "SpleenpercentLymph" )

fcs.condition.label <- c( 
    "LymphpercentSpleen" = "LymphpercentSpleen",
    "Sp_" = "Spleen",
    "SpleenpercentLymph" = "SpleenpercentLymph"
)

# Run the get_channels script and insert the output below.
fcs.channel <- c( 
    "FJComp-APC-A", 
    "FJComp-APC-Fire 750-A",
    "FJComp-APC-Fire810-A",
    "FJComp-Alexa Fluor 532-A", 
    "FJComp-Alexa Fluor 561-A", 
    "FJComp-Alexa Fluor 700-A",  
    "FJComp-Alexa Fluor 790-A", 
    "FJComp-BB515-A", 
    "FJComp-BB660-P2-A", 
    "FJComp-BB700-A", 
    "FJComp-BB755-P-A", 
    "FJComp-BB790-P-A", 
    "FJComp-BUV395-A", 
    "FJComp-BUV496-A",
    "FJComp-BUV563-A",
    "FJComp-BUV615-A",
    "FJComp-BUV661-A",
    "FJComp-BUV737-A",
    "FJComp-BUV805-A",
    "FJComp-BV421-A", 
    "FJComp-BV480-A", 
    "FJComp-BV510-A",
    "FJComp-BV570-A",
    "FJComp-BV605-A",
    "FJComp-BV650-A",
    "FJComp-BV711-A", 
    "FJComp-BV750-A",
    "FJComp-Nova Blue 530-A",
    "FJComp-Nova Blue 585-A",
    "FJComp-Nova Blue 610-A",
    "FJComp-Nova Red 685-A",
    "FJComp-Nova Yellow 690-A",
    "FJComp-PE Fire 640-A",
    "FJComp-PE Fire 810-A",
    "FJComp-PE-A",
    "FJComp-PE-Cy5-A",
    "FJComp-PE-Cy5.5-A",
    "FJComp-PE-Cy7-A",
    "FJComp-PE-Dazzle594-A",
    "FJComp-Pacific Blue-A",
    "FJComp-Pacific Orange-A",
    "FJComp-PerCP-A",
    "FJComp-PerCP-eFluor 710-A",
    "FJComp-Qdot 705-A",
    "FJComp-SBUV445-A",
    "FJComp-SBV515-A",
    "FJComp-Super Bright 436-A",
    "FJComp-Super Bright 780-A",
    "FJComp-eFluor 660-A"
)

fcs.channel.label <- c( 
    "FJComp-APC-A" = "RORgT", 
    "FJComp-APC-Fire 750-A" = "CD44",
    "FJComp-APC-Fire810-A" = "Gr-1",
    "FJComp-Alexa Fluor 532-A" = "IgM", 
    "FJComp-Alexa Fluor 561-A" = "F4-80", 
    "FJComp-Alexa Fluor 700-A" = "Ki67",  
    "FJComp-Alexa Fluor 790-A" = "CD90.2", 
    "FJComp-BB515-A" = "CCR9", 
    "FJComp-BB660-P2-A" = "TCRgd", 
    "FJComp-BB700-A" = "PDCA-1", 
    "FJComp-BB755-P-A" = "CD11c", 
    "FJComp-BB790-P-A" = "Ly-6C", 
    "FJComp-BUV395-A" = "CD103", 
    "FJComp-BUV496-A" = "IgD",
    "FJComp-BUV563-A" = "NK1.1",
    "FJComp-BUV615-A" = "CTLA-4",
    "FJComp-BUV661-A" = "c-Kit",
    "FJComp-BUV737-A" = "CD62L",
    "FJComp-BUV805-A" = "GITR",
    "FJComp-BV421-A" = "CD150", 
    "FJComp-BV480-A" = "CXCR3", 
    "FJComp-BV510-A" = "Siglec F",
    "FJComp-BV570-A" = "TCRb",
    "FJComp-BV605-A" = "PD-1",
    "FJComp-BV650-A" = "XCR1",
    "FJComp-BV711-A" = "CD127", 
    "FJComp-BV750-A" = "CCR2",
    "FJComp-Nova Blue 530-A" = "CD45",
    "FJComp-Nova Blue 585-A" = "CD4",
    "FJComp-Nova Blue 610-A" = "CD8",
    "FJComp-Nova Red 685-A" = "CD3",
    "FJComp-Nova Yellow 690-A" = "CD19",
    "FJComp-PE Fire 640-A" = "CD11b",
    "FJComp-PE Fire 810-A" = "CD38",
    "FJComp-PE-A" = "GATA-3",
    "FJComp-PE-Cy5-A" = "CD86",
    "FJComp-PE-Cy5.5-A" = "Foxp3",
    "FJComp-PE-Cy7-A" = "CD172a",
    "FJComp-PE-Dazzle594-A" = "CD64",
    "FJComp-Pacific Blue-A" = "Helios",
    "FJComp-Pacific Orange-A" = "CD24",
    "FJComp-PerCP-A" = "MHCII",
    "FJComp-PerCP-eFluor 710-A" = "NKp46",
    "FJComp-Qdot 705-A" = "CD69",
    "FJComp-SBUV445-A" = "B220",
    "FJComp-SBV515-A" = "CD25",
    "FJComp-Super Bright 436-A" = "ICOS",
    "FJComp-Super Bright 780-A" = "KLRG1",
    "FJComp-eFluor 660-A" = "T-bet"
)

# Tip: set each to 200 to start.
# Modify (increasing) until the density plots resemble the distributions you are familiar with.
fcs.channel.asinh.scale <- c( 
    "FJComp-APC-A" = 2000, 
    "FJComp-APC-Fire 750-A" = 4000, 
    "FJComp-APC-Fire810-A" = 5000, 
    "FJComp-Alexa Fluor 532-A" = 2000, 
    "FJComp-Alexa Fluor 561-A" = 8000, 
    "FJComp-Alexa Fluor 700-A" = 3000, 
    "FJComp-Alexa Fluor 790-A" = 2000, 
    "FJComp-BB515-A" = 4000,  
    "FJComp-BB660-P2-A" = 2000,  
    "FJComp-BB700-A" = 4000, 
    "FJComp-BB755-P-A" = 2000, 
    "FJComp-BB790-P-A" = 4000, 
    "FJComp-BUV395-A" = 1000, 
    "FJComp-BUV496-A" = 2000, 
    "FJComp-BUV563-A" = 3000, 
    "FJComp-BUV615-A" = 1000, 
    "FJComp-BUV661-A" = 2000, 
    "FJComp-BUV737-A" = 2000, 
    "FJComp-BUV805-A" = 2000, 
    "FJComp-BV421-A" = 6000, 
    "FJComp-BV480-A" = 3000, 
    "FJComp-BV510-A" = 20000, 
    "FJComp-BV570-A" = 3000, 
    "FJComp-BV605-A" = 2000,
    "FJComp-BV650-A" = 2000, 
    "FJComp-BV711-A" = 1000, 
    "FJComp-BV750-A" = 3000, 
    "FJComp-Nova Blue 530-A" = 1000, 
    "FJComp-Nova Blue 585-A" = 4000, 
    "FJComp-Nova Blue 610-A" = 3000, 
    "FJComp-Nova Red 685-A" = 2000, 
    "FJComp-Nova Yellow 690-A" = 5000, 
    "FJComp-PE Fire 640-A" = 10000, 
    "FJComp-PE Fire 810-A" = 2500, 
    "FJComp-PE-A" = 12000, 
    "FJComp-PE-Cy5-A" = 8000, 
    "FJComp-PE-Cy5.5-A" = 8000, 
    "FJComp-PE-Cy7-A" = 3000, 
    "FJComp-PE-Dazzle594-A" = 3000, 
    "FJComp-Pacific Blue-A" = 4000, 
    "FJComp-Pacific Orange-A" = 5000, 
    "FJComp-PerCP-A" = 6000,
    "FJComp-PerCP-eFluor 710-A" = 4000, 
    "FJComp-Qdot 705-A" = 5000, 
    "FJComp-SBUV445-A" = 1500,
    "FJComp-SBV515-A" = 2000, 
    "FJComp-Super Bright 436-A" = 8000, 
    "FJComp-Super Bright 780-A" = 8000,
    "FJComp-eFluor 660-A" = 3000
)

fcs.condition.n <- length( fcs.condition )
fcs.channel.n <- length( fcs.channel )


# general parameters
# Tip: edit the directory to match the location of the script files. 
fcs.src.dir <- "D:/Carlos tSNE script/00_src"

# Tip: use today's date
fcs.seed.base <- 20210608

# Tip: Set to FALSE while optimizing the tSNE and FlowSOM. 
# Once you have a good run, set to TRUE and then change labels, colors, etc.
fcs.use.cached.results <- FALSE

fcs.sample.number.width <- 2
fcs.event.number.width <- 6


# graphics parameters

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


# density parameters
# Tip: set to 50000 while you are setting the asinh scaling.
# Once done, set to NULL and re-run the script on all the data.

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


# cluster parameters
# Tip: In general use the same number for fcs.cluster.n and fcs.flow.som.dim unless you know why to do otherwise.
# Tip: use more clusters for more diverse cell collections. Use as few as possible to speed up processing.
# For example, for a broad immune phenotyping panel use 40-50 clusters.
# For pre-gated cell types (e.g., CD8s), use 10-12 clusters.
fcs.cluster.n <- 40

# Tip: naming clusters can be done automatically as clusters 1:n via the first command below.
# To rename your clusters based on marker expression, insert a # before the first command,
# remove the # before the second command,
# and rename the clusters appropriately.
fcs.cluster <- sprintf( "%02d", 1 : fcs.cluster.n )
#fcs.cluster <- c("B cells", "Naive CD4 T cells", "Activated CD4 T cells", 
#                 "Naive CD8 T cells", "Activated CD8 T cells", "Macrophages",
#                 "Monocytes", "Neutrophils", "NK cells", 
#                 "cDCs", "pDCs", "Eosinophils")


fcs.cluster.label <- sprintf( "Cluster-%s", fcs.cluster )
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

fcs.flow.som.dim <- 40

fcs.cluster.table.dir <- "./table_cluster"
fcs.cluster.data.dir <- "./data_cluster"

fcs.cluster.table.counts <- "cluster_counts"

fcs.cluster.data <- "cluster_group"


# heatmap parameters

fcs.heatmap.palette.n <- 100
fcs.heatmap.palette <- colorRampPalette( brewer.pal( 9, "YlOrRd" ) )( 
    fcs.heatmap.palette.n )

fcs.heatmap.font.size <- 2.5

fcs.heatmap.label.factor.row <- 1.4
fcs.heatmap.label.factor.col <- 1.4

fcs.heatmap.width <- 2500
fcs.heatmap.height <- 2500

fcs.heatmap.figure.dir <- "./figure_heatmap"

fcs.heatmap.figure <- "heatmap"


# histogram parameters

fcs.histogram.font.size <- 7
fcs.histogram.error.bar.size <- 0.2
fcs.histogram.legend.key.size <- 0.7

fcs.histogram.label.factor.height <- 0.05

fcs.histogram.width <- 4.5
fcs.histogram.height <- 3

fcs.histogram.figure.dir <- "./figure_histogram"

fcs.histogram.figure <- "histogram"


# dimensionality reduction parameters
# Tip: Use the third option in most cases.
# More cells will take longer. 
# You might downsample to 50-100k cells at first, then run again with more cells if you like the result.
fcs.dmrd.data.sample.n <- NULL
fcs.dmrd.data.sample.n.per.condition <- 8297
fcs.dmrd.data.sample.n.per.sample <- NULL

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
fcs.dmrd.figure.nrow <- 2
fcs.dmrd.figure.ncol <- ceiling( fcs.condition.n / fcs.dmrd.figure.nrow )

fcs.dmrd.figure.width <- 3
fcs.dmrd.figure.height <- 3


# tsne parameters
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


# umap parameters
# Tip: you don't generally need to change the UMAP iterations
fcs.umap.iter.n <- 1000

fcs.umap.figure.lims.factor <- 0.8
fcs.umap.figure.point.size <- 1.2

fcs.umap.figure.dir <- "./figure_umap"

fcs.umap.figure.plot <- "umap_plot"

fcs.umap.cache.file.path <- "./umap_cache.dat"


# cross-entropy test parameters
# Tips: Set this depending on your RAM and number of groups.
# You won't be able to analyze more than about 100k cells unless you have >32GB RAM.
# The crossentropy test works best with at least 10k cells per group.
# Multiple hypothesis testing will greatly reduce your ability to distinguish statistical differences.
fcs.ce.diff.prob.sample.n <- NULL

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


# cross-entropy test parameters for tsne

fcs.ce.diff.tsne.perplexity.factor <- 3

fcs.ce.diff.tsne.figure.dir <- "./figure_tsne_ce_diff"

fcs.ce.diff.tsne.figure.cdf <- "tsne_ce_diff_cdf"
fcs.ce.diff.tsne.figure.dendrogram <- "tsne_ce_diff_dendrogram"
fcs.ce.diff.tsne.result <- "tsne_ce_diff_result"

fcs.ce.diff.tsne.cache.file.path <- "./tsne_ce_diff_cache.dat"


# cross-entropy test parameters for umap

fcs.ce.diff.umap.figure.dir <- "./figure_umap_ce_diff"

fcs.ce.diff.umap.figure.cdf <- "umap_ce_diff_cdf"
fcs.ce.diff.umap.figure.dendrogram <- "umap_ce_diff_dendrogram"
fcs.ce.diff.umap.result <- "umap_ce_diff_result"

fcs.ce.diff.umap.cache.file.path <- "./umap_ce_diff_cache.dat"

