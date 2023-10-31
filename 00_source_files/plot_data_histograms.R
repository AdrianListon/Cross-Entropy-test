set.seed.here( fcs.seed.base, "plot density distributions of data" )

{
  if ( ! is.null( fcs.density.data.sample.n ) && 
       fcs.density.data.sample.n < dmrd.event.n )
    density.data.idx <- sort( sample( dmrd.event.n, 
                                      fcs.density.data.sample.n ) )
  else
    density.data.idx <- 1 : dmrd.event.n
}

density.data <- dmrd.data[ density.data.idx, ]

density.frame <- data.frame(density.data)
density.frame$Sample <- rownames(density.frame)
density.data.ggdf.all <- pivot_longer( density.frame, !Sample, cols_vary = "slowest",
                                             names_to = "channel" , values_to = "density.value" )
density.data.ggdf.all$partition <- fcs.density.partition.all
density.data.ggdf.all$channel <- as.character( density.data.ggdf.all$channel )

density.data.ggdf.condition <- lapply( fcs.condition, function( fc ) {
  ggdf.condition <- pivot_longer( 
    density.frame[ dmrd.event.condition[ density.data.idx ] == fc, , 
                  drop = FALSE], !Sample, cols_vary = "slowest",
    names_to = "channel" , values_to = "density.value" )
  if ( nrow( ggdf.condition ) > 0 ) {
    ggdf.condition$partition <- fc
    ggdf.condition$channel <- as.character( ggdf.condition$channel )
  }
  ggdf.condition
} )
density.data.ggdf.condition <- do.call( rbind, density.data.ggdf.condition )

density.data.ggdf.sample <- lapply( flow.sample, function( fs ) {
  ggdf.sample <- pivot_longer( 
    density.frame[ dmrd.event.sample[ density.data.idx ] == fs, , 
                  drop = FALSE ], !Sample, cols_vary = "slowest",
    names_to = "channel" , values_to = "density.value" )
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

if( input.file.type == 2 ){
  channel.factor <- fcs.channel.label
  names(channel.factor) <- unique(density.data.ggdf$channel)
  #density.data.ggdf$channel <- channel.factor[density.data.ggdf$channel]
}else {
  channel.factor <- fcs.channel.label
}

#density.data.ggdf$channel <- factor( density.data.ggdf$channel, 
#                                     levels = channel.factor )

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
    channel = channel.factor ), scales = "free_x" ) + 
  theme_ridges( line_size = fcs.density.line.size, 
                font_size = fcs.density.font.size ) + 
  theme( strip.background = element_rect( fill = "white" ) )

ggsave( 
  file.path( fcs.density.figure.dir, 
             sprintf( "%s.jpeg", fcs.density.figure.sample ) ), 
  density.plot, 
  width = fcs.density.figure.width.base * ( fcs.channel.n + 1 ), 
  height = fcs.density.figure.height.base * ( density.data.partition.n + 1 ) 
)
