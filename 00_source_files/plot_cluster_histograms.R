set.seed.here( fcs.seed.base, "plot density distributions of data" )

{
  if ( ! is.null( fcs.density.data.sample.n ) && 
       fcs.density.data.sample.n < dmrd.data.n )
    cluster.density.data.idx <- sort( sample( dmrd.data.n, 
                                      fcs.density.data.sample.n ) )
  else
    cluster.density.data.idx <- 1 : dmrd.data.n
}

cluster.density.data <- dmrd.data[ cluster.density.data.idx, ]

cluster.density.frame <- data.frame(cluster.density.data)
cluster.density.frame$Sample <- rownames(cluster.density.frame)


density.data.ggdf.cluster <- lapply( fcs.cluster, function( fc ) { 
  ggdf.cluster <- pivot_longer( 
    cluster.density.frame[ dmrd.event.cluster[ cluster.density.data.idx ] == fc, , 
                  drop = FALSE ], 
    !Sample, cols_vary = "slowest",
    names_to = "channel" , values_to = "density.value"
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

if( input.file.type == 2 ){
  channel.factor <- fcs.channel.label
  names(channel.factor) <- unique(density.data.ggdf$channel)
  #density.data.ggdf$channel <- channel.factor[density.data.ggdf$channel]
}else {
  channel.factor <- fcs.channel.label
}

#density.data.ggdf$channel <- factor( density.data.ggdf$channel, 
#                                     levels = channel.factor )

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
    channel = channel.factor ), scales = "free_x" ) + 
  theme_ridges( line_size = fcs.density.line.size, 
                font_size = fcs.density.font.size ) + 
  theme( strip.background = element_rect( fill = "white" ) )

ggsave( 
  file.path( fcs.density.figure.dir, 
             sprintf( "%s.jpeg", fcs.density.figure.cluster ) ), 
  density.plot, 
  width = fcs.density.figure.width.base * ( fcs.channel.n + 1 ), 
  height = fcs.density.figure.height.base * ( density.data.partition.n + 1 ) 
)


