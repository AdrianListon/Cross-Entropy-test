
# script for the analysis of flow cytometry data---------------

library( digest )
require( dunn.test )
library( ggplot2 )
library( ggridges )
require( RANN )
library( RColorBrewer )
library( reshape2 )
library( Rtsne )
library( umap )
library( dplyr )
library( EmbedSOM )

# function to set random seed depending on base number and string---------------

set.seed.here <- function( seed.base, seed.char )
{
    seed.add <- strtoi( substr( digest( seed.char, "xxhash32" ), 2, 8 ), 16 )
    seed.new <- seed.base + seed.add    
    set.seed( seed.new )
    invisible( seed.new )
}

# source parameters---------------

param.filename <- "./analyze_flow_cytometry_parameter_csv.r"

source( param.filename )

# select channels for analysis---------------

channel.selection.data <- read.csv( list.files(fcs.data.dir, "\\.csv$", full.names = TRUE)[1] )
channels <- data.frame( desc = unname( colnames( channel.selection.data ) ))
for( i in 1:dim(channels)[1] ){
    channels$out[i] = paste0("\"", channels$desc[i], "\" = \"", channels$desc[i], "\",\n") 
}
descs.filtered <- channels$desc[!is.na(channels$desc) & channels$desc != '-'] 
channels.filtered <- filter(channels, desc %in% descs.filtered)

# paste the output below into fcs.channel in the parameter file, deleting unwanted channels
cat('"', paste0(channels.filtered$desc, collapse = '","'), '"', sep = "")

source( param.filename )

# add the output to fcs.channel.label in the parameter file
temp <- channels$out[channels$desc %in% fcs.channel]
temp[length(temp)] <- gsub(',', '', temp[length(temp)])
cat(paste(temp, collapse = ""))

source( param.filename )


# check consistency in parameters---------------

stopifnot( names( fcs.channel.label ) == fcs.channel )

stopifnot( names( fcs.condition.label ) == fcs.condition )
stopifnot( names( fcs.condition.color ) == fcs.condition )
stopifnot( names( fcs.condition.line.type ) == fcs.condition )

stopifnot( names( fcs.cluster.label ) == fcs.cluster )
stopifnot( names( fcs.cluster.color ) == fcs.cluster )
stopifnot( names( fcs.cluster.line.type ) == fcs.cluster )

stopifnot( sum( is.null( fcs.dmrd.data.sample.n ), 
    is.null( fcs.dmrd.data.sample.n.per.condition ), 
    is.null( fcs.dmrd.data.sample.n.per.sample ) ) >= 2 )

stopifnot( unlist( fcs.cluster.group ) >= 1 & 
    unlist( fcs.cluster.group ) <= fcs.cluster.n )


# source functions---------------

source( file.path( fcs.src.dir, "ce_diff_test.r" ) )
source( file.path( fcs.src.dir, "ce_diff_test_tsne.r" ) )
source( file.path( fcs.src.dir, "ce_diff_test_umap.r" ) )
source( file.path( fcs.src.dir, "plot_all_dmrd_figures.r" ) )
source( file.path( fcs.src.dir, "plot_all_embedsom_figures.r" ) )
source( file.path( fcs.src.dir, "plot_dimensionality_reduction.r" ) )
source( file.path( fcs.src.dir, "plot_dimensionality_reduction_embedSOM.r" ) )


# create dirs---------------

figure.dir <- c( 
    fcs.ce.diff.tsne.figure.dir, 
    fcs.ce.diff.umap.figure.dir, 
    fcs.density.figure.dir, 
    fcs.heatmap.figure.dir, 
    fcs.histogram.figure.dir, 
    fcs.tsne.figure.dir, 
    fcs.umap.figure.dir,
    fcs.embedsom.figure.dir,
    fcs.ce.diff.embedsom.figure.dir
)

table.dir <- fcs.cluster.table.dir

data.dir <- sapply( names( fcs.cluster.group ), function( fcg.name ) 
    sprintf( "%s/%s_%s", fcs.cluster.data.dir, fcs.cluster.data, fcg.name ) )

for ( the.dir in c( figure.dir, table.dir, data.dir ) )
    if ( ! file.exists( the.dir ) )
        dir.create( the.dir, recursive = TRUE )


# read csv data---------------

flow.data.filename.all <- list.files( fcs.data.dir, "\\.csv$" )

flow.data.filename <- grep( paste0( fcs.condition, collapse = "|" ), 
    flow.data.filename.all, value = TRUE )

sample.name.format <- paste0( "%s.%0", fcs.sample.number.width, "d" )
event.name.format <- paste0( "%s.%0", fcs.event.number.width, "d" )

flow.data.filename.sample <- rep( "", length( flow.data.filename ) )
names( flow.data.filename.sample ) <- flow.data.filename

sample.idx.next <- rep( 1, fcs.condition.n )
names( sample.idx.next ) <- fcs.condition

flow.data <- lapply( flow.data.filename, function( flow.data.fn ) {
 #   cat( flow.data.fn, "\n" )
    sample.data <- as.matrix(read.csv( file.path( fcs.data.dir, flow.data.fn ) ))
    
    condition <- fcs.condition[ sapply( fcs.condition, grepl, flow.data.fn ) ]
    stopifnot( length( condition ) == 1 )
    
    if ( ! all( fcs.channel %in% colnames( sample.data ) ) )
    {
        cat( sprintf( "File: %s\n", flow.data.fn ) )
        print( sort( fcs.channel[ 
            ! fcs.channel %in% colnames( sample.data ) ] ) )
        print( sort( colnames( sample.data ) ) )
        stop( "mismatch in names of channels" )
    }

    sample.name <- sprintf( sample.name.format, condition, 
        sample.idx.next[ condition ] )
    
    sample.data <- sample.data[ , fcs.channel, drop = FALSE ]
    
    event.n <- nrow( sample.data )
    if ( event.n > 0 ) {
        event.name <- sprintf( event.name.format, sample.name, 1 : event.n )
        rownames( sample.data ) <- event.name
    }
    
    flow.data.filename.sample[ flow.data.fn ] <<- sample.name
    sample.idx.next[ condition ] <<- sample.idx.next[ condition ] + 1
    
    sample.data
} )

flow.data <- do.call( rbind, flow.data )

# define samples---------------
flow.sample <- flow.data.filename.sample
names( flow.sample ) <- NULL

stopifnot( flow.sample == 
    unique( sub( "\\.[0-9]+$", "", rownames( flow.data ) ) ) )

flow.sample.n <- length ( flow.sample )

flow.sample.condition <- factor( sub( "\\.[0-9]+$", "", flow.sample ), 
    levels = fcs.condition )
names( flow.sample.condition ) <- flow.sample

# reorder samples to follow order of conditions
flow.sample <- flow.sample[ order( flow.sample.condition ) ]
flow.sample.condition <- flow.sample.condition[ flow.sample ]

flow.sample.label <- sapply( flow.sample, function( fs ) {
    sample.cond <- sub( "^(.*)\\.[0-9]+$", "\\1", fs )
    sample.num <- sub( "^.*\\.([0-9]+)$", "\\1", fs )
    sprintf( "%s-%s", fcs.condition.label[ sample.cond ], sample.num )
} )

flow.sample.filename <- sapply( flow.sample, function( fs  ) 
    names( which( flow.data.filename.sample == fs ) ) )

# define events
flow.event <- rownames( flow.data )
flow.event.n <- length( flow.event )

flow.event.sample <- factor( sub( "\\.[0-9]+$", "", flow.event ), 
    levels = flow.sample )
names( flow.event.sample ) <- flow.event

flow.event.condition <- factor( sub( "\\.[0-9]+$", "", flow.event.sample ), 
    levels = fcs.condition )
names( flow.event.condition ) <- flow.event

# reorder events to follow order of samples
flow.event.order <- order( flow.event.sample )

flow.data <- flow.data[ flow.event.order, ]
flow.event <- flow.event[ flow.event.order ]
flow.event.sample <- flow.event.sample[ flow.event.order ]
flow.event.condition <- flow.event.condition[ flow.event.order ]

flow.event.sample.n <- as.vector( table( flow.event.sample ) )
names( flow.event.sample.n ) <- flow.sample

flow.event.condition.n <- as.vector( table( flow.event.condition ) )
names( flow.event.condition.n ) <- fcs.condition

flow.data.filename

table( flow.sample.condition )

str( flow.data )
flow.event.condition.n
flow.event.sample.n


# define figure parameters for samples---------------

flow.sample.color <- fcs.condition.color[ flow.sample.condition ]
names( flow.sample.color ) <- flow.sample

flow.sample.color.single <- unlist( lapply( fcs.condition, function( fc ) {
    cond.sample.n <- sum( flow.sample.condition == fc )
    rep( 
        fcs.color.pool, 
        ceiling( cond.sample.n / fcs.color.pool.n ) 
    )[ 1 : cond.sample.n ]
} ) )
names( flow.sample.color.single ) <- flow.sample

flow.sample.line.type <- fcs.condition.line.type[ flow.sample.condition ]
names( flow.sample.line.type ) <- flow.sample

flow.sample.line.type.single <- unlist( lapply( fcs.condition, function( fc ) {
    cond.sample.n <- sum( flow.sample.condition == fc )
    rep( 
        fcs.line.type.pool, 
        ceiling( cond.sample.n / fcs.line.type.pool.n ) 
    )[ 1 : cond.sample.n ]
} ) )
names( flow.sample.line.type.single ) <- flow.sample

flow.ce.diff.figure.dendrogram.weight.sample <- 
    fcs.ce.diff.figure.dendrogram.weight.condition[ flow.sample.condition ]
names( flow.ce.diff.figure.dendrogram.weight.sample ) <- flow.sample


# plot density distributions of transformed data---------------

set.seed.here( fcs.seed.base, "plot density distributions of transformed data" )

{
if ( ! is.null( fcs.density.data.sample.n ) && 
    fcs.density.data.sample.n < flow.event.n )
    density.data.idx <- sort( sample( flow.event.n, 
        fcs.density.data.sample.n ) )
else
    density.data.idx <- 1 : flow.event.n
}

density.data <- flow.data[ density.data.idx, ]

density.data.ggdf.all <- melt( density.data, 
    varnames = c( "partition", "channel" ), value.name = "density.value" )
density.data.ggdf.all$partition <- fcs.density.partition.all
density.data.ggdf.all$channel <- as.character( density.data.ggdf.all$channel )

density.data.ggdf.condition <- lapply( fcs.condition, function( fc ) {
    ggdf.condition <- melt( 
        density.data[ flow.event.condition[ density.data.idx ] == fc, , 
            drop = FALSE], 
        varnames = c( "partition", "channel" ), 
        value.name = "density.value" )
    if ( nrow( ggdf.condition ) > 0 ) {
        ggdf.condition$partition <- fc
        ggdf.condition$channel <- as.character( ggdf.condition$channel )
    }
    ggdf.condition
} )
density.data.ggdf.condition <- do.call( rbind, density.data.ggdf.condition )

density.data.ggdf.sample <- lapply( flow.sample, function( fs ) {
    ggdf.sample <- melt( 
        density.data[ flow.event.sample[ density.data.idx ] == fs, , 
            drop = FALSE ], 
        varnames = c( "partition", "channel" ), 
        value.name = "density.value" )
    if ( nrow( ggdf.sample ) > 0 ) {
        ggdf.sample$partition <- fs
        ggdf.sample$channel <- as.character( ggdf.sample$channel )
    }
    ggdf.sample
} )
density.data.ggdf.sample <- do.call( rbind, density.data.ggdf.sample )

density.data.ggdf <- rbind( density.data.ggdf.all, density.data.ggdf.condition, 
    density.data.ggdf.sample )

density.data.partition <- c( fcs.density.partition.all, fcs.condition, 
    flow.sample )
density.data.partition.n <- length( density.data.partition )

density.data.partition.label <- c( fcs.density.partition.all.label, 
    fcs.condition.label, flow.sample.label )

density.data.ggdf$partition <- factor( density.data.ggdf$partition, 
    levels = density.data.partition )

density.data.ggdf$channel <- factor( density.data.ggdf$channel, 
    levels = fcs.channel )

density.plot.color <- c( fcs.density.partition.all.color, fcs.condition.color, 
    flow.sample.color )

density.plot <- ggplot( density.data.ggdf, aes( x = density.value, 
        y = partition, color = partition, fill = partition ) ) + 
    geom_density_ridges( size = fcs.density.line.size, 
        alpha = fcs.density.line.alpha, show.legend = FALSE ) + 
    labs( x = NULL, y = NULL ) + 
    scale_y_discrete( limits = rev( density.data.partition ), 
        breaks = rev( density.data.partition ), 
        labels = rev( density.data.partition.label ) ) + 
    scale_color_manual( values = density.plot.color ) + 
    scale_fill_manual( values = density.plot.color ) + 
    facet_wrap( vars( channel ), nrow = 1, labeller = labeller( 
        channel = fcs.channel.label ) ) + 
    theme_ridges( line_size = fcs.density.line.size, 
        font_size = fcs.density.font.size ) + 
    theme( strip.background = element_rect( fill = "white" ) )

ggsave( 
    file.path( fcs.density.figure.dir, 
        sprintf( "%s.png", fcs.density.figure.sample ) ), 
    density.plot, 
    width = fcs.density.figure.width.base * ( fcs.channel.n + 1 ), 
    height = fcs.density.figure.height.base * ( density.data.partition.n + 1 ) 
)


# select data for dimensionality reduction---------------

set.seed.here( fcs.seed.base, "select data for dimensionality reduction" )

{
if ( ! is.null( fcs.dmrd.data.sample.n ) )
{
    if ( fcs.dmrd.data.sample.n < flow.event.n )
        dmrd.data.idx <- sort( sample( flow.event.n, fcs.dmrd.data.sample.n ) )
    else
        dmrd.data.idx <- 1 : flow.event.n
}
else if ( ! is.null( fcs.dmrd.data.sample.n.per.condition ) )
{
    dmrd.data.idx <- unlist( sapply( fcs.condition, function( fc ) {
        fc.idx <- which( flow.event.condition == fc )
        if ( fcs.dmrd.data.sample.n.per.condition < length( fc.idx ) )
            sort( sample( fc.idx, fcs.dmrd.data.sample.n.per.condition ) )
        else
            fc.idx
    } ) )
    names( dmrd.data.idx ) <- NULL
}
else if ( ! is.null( fcs.dmrd.data.sample.n.per.sample ) )
{
    dmrd.data.idx <- unlist( sapply( flow.sample, function( fs ) {
        fs.idx <- which( flow.event.sample == fs )
        if ( fcs.dmrd.data.sample.n.per.sample < length( fs.idx ) )
            sort( sample( fs.idx, fcs.dmrd.data.sample.n.per.sample ) )
        else
            fs.idx
    } ) )
    names( dmrd.data.idx ) <- NULL
}
else
    dmrd.data.idx <- 1 : flow.event.n
}

dmrd.data <- flow.data[ dmrd.data.idx, ]

dmrd.event.sample <- flow.event.sample[ dmrd.data.idx ]
dmrd.event.condition <- flow.event.condition[ dmrd.data.idx ]

str( dmrd.data )
table( dmrd.event.condition )
table( dmrd.event.sample )


# get flowsom clusters---------------

set.seed.here( fcs.seed.base, "get flowsom clusters" )

# build som objects--embedsom method
flow.som <- EmbedSOM::SOM(flow.data, xdim = fcs.flow.som.dim, 
                          ydim = fcs.flow.som.dim, batch = TRUE,
                          parallel = TRUE, threads = fcs.tsne.thread.n )

# get clusters
embedsom.cluster <- hclust( dist( flow.som$codes ) )
flow.som.event.cluster <- cutree( embedsom.cluster, fcs.cluster.n )[flow.som$mapping[ , 1] ]

# check quality of flowSOM clustering (optional)---------------
embed.som <- EmbedSOM::EmbedSOM(data=flow.data, map=flow.som, 
                                parallel = TRUE, threads = fcs.tsne.thread.n)

# reorder clusters from bigger to smaller---------------
flow.som.cluster.rank <- 1 + fcs.cluster.n - 
    rank( table( flow.som.event.cluster ), ties.method = "last" )
flow.som.event.cluster <- flow.som.cluster.rank[ flow.som.event.cluster ]
names( flow.som.event.cluster ) <- NULL

# set clusters as a factor
flow.som.event.cluster <- factor( flow.som.event.cluster, 
    levels = 1 : fcs.cluster.n )
levels( flow.som.event.cluster ) <- fcs.cluster

# reorder events
flow.event.cluster <- flow.som.event.cluster
dmrd.event.cluster <- flow.event.cluster[ dmrd.data.idx ]

flow.event.cluster.n <- as.vector( table( flow.event.cluster ) )
names( flow.event.cluster.n ) <- fcs.cluster

dmrd.event.cluster.n <- as.vector( table( dmrd.event.cluster ) )
names( dmrd.event.cluster.n ) <- fcs.cluster

length( flow.event.cluster )
table( flow.event.cluster )

length( dmrd.event.cluster )
table( dmrd.event.cluster )


# save cluster counts---------------

flow.cluster.count <- sapply( fcs.cluster, function( fc )
    table( flow.event.sample[ flow.event.cluster == fc ] ) )

stopifnot( rownames( flow.cluster.count ) == flow.sample )
stopifnot( colnames( flow.cluster.count ) == fcs.cluster )

rownames( flow.cluster.count ) <- flow.sample.label
colnames( flow.cluster.count ) <- fcs.cluster.label

write.csv( flow.cluster.count, file = file.path( fcs.cluster.table.dir, 
    sprintf( "%s.csv", fcs.cluster.table.counts ) ) )

flow.cluster.percent <- flow.cluster.count/rowSums(flow.cluster.count)*100

write.csv( flow.cluster.percent, file = file.path( fcs.cluster.table.dir, 
    sprintf( "%s.csv", "cluster_proportions" ) ) )


# plot density distributions by cluster---------------

density.data.ggdf.cluster <- lapply( fcs.cluster, function( fc ) { 
    ggdf.cluster <- melt( 
        density.data[ flow.event.cluster[ density.data.idx ] == fc, , 
            drop = FALSE ], 
        varnames = c( "partition", "channel" ), 
        value.name = "density.value" 
    )
    if ( nrow( ggdf.cluster ) > 0 ) {
        ggdf.cluster$partition <- fc
        ggdf.cluster$channel <- as.character( ggdf.cluster$channel )
    }
    ggdf.cluster
} )
density.data.ggdf.cluster <- do.call( rbind, density.data.ggdf.cluster )

density.data.ggdf <- rbind( density.data.ggdf.all, density.data.ggdf.cluster )

density.data.partition <- c( fcs.density.partition.all, fcs.cluster )
density.data.partition.n <- length( density.data.partition )

density.data.partition.label <- c( fcs.density.partition.all.label, 
    fcs.cluster.label )

density.data.ggdf$partition <- factor( density.data.ggdf$partition, 
    levels = density.data.partition )

density.data.ggdf$channel <- factor( density.data.ggdf$channel, 
    levels = fcs.channel )

density.plot.color <- c( fcs.density.partition.all.color, fcs.cluster.color )

density.plot <- ggplot( density.data.ggdf, aes( x = density.value, 
        y = partition, color = partition, fill = partition ) ) + 
    geom_density_ridges( size = fcs.density.line.size, 
        alpha = fcs.density.line.alpha, show.legend = FALSE ) + 
    labs( x = NULL, y = NULL ) + 
    scale_y_discrete( limits = rev( density.data.partition ), 
        breaks = rev( density.data.partition ), 
        labels = rev( density.data.partition.label ) ) + 
    scale_color_manual( values = density.plot.color ) + 
    scale_fill_manual( values = density.plot.color ) + 
    facet_wrap( vars( channel ), nrow = 1, labeller = labeller( 
        channel = fcs.channel.label ) ) + 
    theme_ridges( line_size = fcs.density.line.size, 
        font_size = fcs.density.font.size ) + 
    theme( strip.background = element_rect( fill = "white" ) )

ggsave( 
    file.path( fcs.density.figure.dir, 
        sprintf( "%s.png", fcs.density.figure.cluster ) ), 
    density.plot, 
    width = fcs.density.figure.width.base * ( fcs.channel.n + 1 ), 
    height = fcs.density.figure.height.base * ( density.data.partition.n + 1 ) 
)


# plot embedsom figures---------------

embedsom.data.max <- max( embed.som )

plot.all.embedsom.figures( 
    embed.som, embedsom.data.max, 
    fcs.embedsom.figure.lims.factor, fcs.embedsom.figure.point.size, 
    fcs.embedsom.figure.dir, "embedSOM_plot", 
    flow.data, flow.event.cluster, flow.event.condition 
)


# plot heatmaps---------------

heatmap.type <- c( "by_condition", "by_sample", "by_cluster" )

for ( ht in heatmap.type )
{
    if ( ht == "by_condition" ) {
        flow.data.group.median <- apply( flow.data, 2, tapply, 
            flow.event.condition, median )
        margin.col <- fcs.heatmap.label.factor.col * 
            max( nchar( fcs.condition.label ) )
        group.label <- fcs.condition.label
        group.color <- fcs.condition.color
    }
    else if ( ht == "by_sample" ) {
        flow.data.group.median <- apply( flow.data, 2, tapply, 
            flow.event.sample, median )
        margin.col <- fcs.heatmap.label.factor.col * 
            max( nchar( flow.sample.label ) )
        group.label <- flow.sample.label
        group.color <- flow.sample.color
    }
    else if ( ht == "by_cluster" ) {
        flow.data.group.median <- apply( flow.data, 2, tapply, 
            flow.event.cluster, median )
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

histogram.type <- c( "by_cluster", "by_sample", "by_assay" )

condition.sample.n <- as.vector( table( flow.sample.condition ) )

for ( ht in histogram.type )
{
    flow.data.cluster.fraction <- lapply( fcs.cluster, function( fc ) {
        sample.event.n <- as.vector( table( 
            flow.event.sample[ flow.event.cluster == fc ] ) )
        
        if ( ht == "by_cluster" )
        {
            sample.fraction <- log2( 
                ( sample.event.n / flow.event.cluster.n[ fc ] ) / 
                ( flow.event.sample.n / flow.event.n ) 
            )
            sample.fraction[ is.infinite( sample.fraction ) ] <- NA
            fraction <- tapply( sample.fraction, flow.sample.condition, mean, 
                na.rm = TRUE )
            std.err <- tapply( sample.fraction, flow.sample.condition, sd, 
                na.rm = TRUE ) / sqrt( condition.sample.n )
        }
        else if ( ht == "by_sample" )
        {
            sample.fraction <- 100 * sample.event.n / flow.event.sample.n
            fraction <- tapply( sample.fraction, flow.sample.condition, mean, 
                na.rm = TRUE )
            std.err <- tapply( sample.fraction, flow.sample.condition, sd, 
                na.rm = TRUE ) / sqrt( condition.sample.n )
        }
        else if ( ht == "by_assay" )
        {
            sample.fraction <- 100 * sample.event.n / flow.event.n
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
            size = fcs.histogram.error.bar.size,
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


# calculate tsne representation---------------

{
    if ( fcs.use.cached.results && file.exists( fcs.tsne.cache.file.path ) )
    {
        cat( "Using cached results for tsne\n" )
        
        load( fcs.tsne.cache.file.path )
    }
    else
    {
        set.seed.here( fcs.seed.base, "calculate tsne representation" )
        
        cat( "Calculating tsne\n" )
        
        tsne.result <- Rtsne( dmrd.data, perplexity = fcs.tsne.perplexity, 
            theta = fcs.tsne.theta, 
            exaggeration_factor = fcs.tsne.exaggeration.factor, 
            max_iter = fcs.tsne.iter.n, check_duplicates = FALSE, pca = FALSE, 
            num_threads = fcs.tsne.thread.n )
        
        save( tsne.result, file = fcs.tsne.cache.file.path )
    }
}

tsne.data <- tsne.result$Y

str( tsne.result )


# plot tsne convergence---------------

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

tsne.data.max <- max( abs( tsne.data ) )

plot.all.dmrd.figures( 
    tsne.data, tsne.data.max, 
    fcs.tsne.figure.lims.factor, fcs.tsne.figure.point.size, 
    fcs.tsne.figure.dir, fcs.tsne.figure.plot, 
    dmrd.data, dmrd.event.cluster, dmrd.event.condition 
)


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
        
        cat( "Calculating umap\n" )
        
        umap.config <- umap.defaults
        umap.config$n_epochs <- fcs.umap.iter.n
        umap.config$verbose <- TRUE
        
        umap.result <- umap( dmrd.data, config = umap.config )
        
        save( umap.result, file = fcs.umap.cache.file.path )
    }
}

umap.data <- umap.result$layout
dimnames( umap.data ) <- NULL

str( umap.result )


# plot umap representations---------------

umap.data.max <- max( abs( apply( umap.data, 2, function( x ) 
    quantile( x, c( 0.25, 0.75 ) ) + c( -1, 1 ) * IQR( x ) ) ) )

plot.all.dmrd.figures( 
    umap.data, umap.data.max, 
    fcs.umap.figure.lims.factor, fcs.umap.figure.point.size, 
    fcs.umap.figure.dir, fcs.umap.figure.plot, 
    dmrd.data, dmrd.event.cluster, dmrd.event.condition 
)


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

print( ce.diff.test.tsne.res )


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

print( ce.diff.test.tsne.res )


# calculate cross-entropy test for tsne by cluster---------------

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

print( ce.diff.test.tsne.res )


# calculate cross-entropy test for umap by condition---------------

set.seed.here( fcs.seed.base, "calculate cross-entropy test for umap" )

ce.diff.test.umap.res <- ce.diff.test.umap( 
    umap.result$knn$distances, umap.result$knn$indexes, 
    umap.data, umap.result$config, 
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

print( ce.diff.test.umap.res )


# calculate cross-entropy test for umap by sample---------------

set.seed.here( fcs.seed.base, "calculate cross-entropy test for umap" )

ce.diff.test.umap.res <- ce.diff.test.umap( 
    umap.result$knn$distances, umap.result$knn$indexes, 
    umap.data, umap.result$config, 
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

print( ce.diff.test.umap.res )


# calculate cross-entropy test for umap by cluster---------------

set.seed.here( fcs.seed.base, "calculate cross-entropy test for umap" )

ce.diff.test.umap.res <- ce.diff.test.umap( 
    umap.result$knn$distances, umap.result$knn$indexes, 
    umap.data, umap.result$config, 
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

print( ce.diff.test.umap.res )


# calculate cross-entropy test for embedsom by condition---------------

set.seed.here( fcs.seed.base, "calculate cross-entropy test for embedSOM" )

ce.diff.test.embedsom.res <- ce.diff.test.tsne( 
    dmrd.data, embed.som, 
    dmrd.event.condition, 
    partition.label = fcs.condition.label, 
    partition.color = fcs.condition.color, 
    partition.line.type = fcs.condition.line.type, 
    base.test = fcs.ce.diff.base.test, 
    base.dist = fcs.ce.diff.base.dist, 
    prob.sample.n = fcs.ce.diff.prob.sample.n, 
    dendrogram.order.weight = fcs.ce.diff.figure.dendrogram.weight.condition, 
    result = file.path( fcs.ce.diff.embedsom.figure.dir, 
                        sprintf( "%s_condition.txt", fcs.ce.diff.embedsom.result ) ), 
    cdf.figure = file.path( fcs.ce.diff.embedsom.figure.dir, 
                            sprintf( "%s_condition.png", fcs.ce.diff.embedsom.figure.cdf ) ), 
    dendrogram.figure = file.path( fcs.ce.diff.embedsom.figure.dir, 
                                   sprintf( "%s_condition.png", fcs.ce.diff.embedsom.figure.dendrogram ) ) )

print( ce.diff.test.embedsom.res )


# output session info---------------

sessionInfo()

