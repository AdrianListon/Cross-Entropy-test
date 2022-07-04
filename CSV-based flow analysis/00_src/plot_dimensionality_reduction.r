
# plots one dimensionality reduction figure


plot.dimensionality.reduction <- function( 
    dmrd.x, dmrd.y, 
    event.group = NULL, group.label = NULL, 
    event.partition = NULL, partition.label = NULL, partition.color = NULL, 
    event.level = NULL, 
    dmrd.lims = NULL, dmrd.nrow = NULL, dmrd.ncol = NULL, 
    point.size = NULL, show.guide = FALSE, guide.name = NULL 
)
{
    if ( ! is.null( event.partition ) )
    {
        if ( is.null( event.group ) )
            ggdf <- data.frame( dmrd.x, dmrd.y, event.partition )
        else
            ggdf <- data.frame( dmrd.x, dmrd.y, event.partition, event.group )
        
        if ( is.null( partition.label ) )
            plot.label <- waiver()
        else
            plot.label <- partition.label
        
        if ( show.guide )
            plot.guide <- guide_legend( keyheight = 0.8, 
                override.aes = list( size = fcs.dmrd.legend.point.size ), 
                label.theme = element_text( size = fcs.dmrd.legend.label.size ), 
                title = guide.name, title.position = "top", 
                title.theme = element_text( size = fcs.dmrd.legend.title.size ) )
        else
            plot.guide <- FALSE
        
        dmrd.plot <- ggplot( ggdf, aes( x = dmrd.x, 
                y = dmrd.y, color = event.partition ) ) + 
            scale_color_manual( values = partition.color, labels = plot.label, 
                guide = plot.guide )
    }
    else if ( ! is.null( event.level ) )
    {
        if ( is.null( event.group ) )
            ggdf <- data.frame( dmrd.x, dmrd.y, event.level )
        else
            ggdf <- data.frame( dmrd.x, dmrd.y, event.level, event.group )
        
        if ( show.guide )
            plot.guide <- guide_colorbar( barwidth = 0.8, barheight = 10, 
                title = guide.name, title.position = "top", 
                title.theme = element_text( size = fcs.dmrd.legend.title.size ) )
        else
            plot.guide <- FALSE
        
        dmrd.plot <- ggplot( ggdf, aes( x = dmrd.x, y = dmrd.y, 
                color = event.level ) ) + 
            scale_color_gradientn( colors = fcs.dmrd.density.palette, 
                labels = NULL, guide = plot.guide )
    }
    else
    {
        if ( is.null( event.group ) )
            ggdf <- data.frame( dmrd.x, dmrd.y )
        else
            ggdf <- data.frame( dmrd.x, dmrd.y, event.group )
        
        dmrd.plot <- ggplot( ggdf, aes( x = dmrd.x, y = dmrd.y ) )
    }
    
    if ( is.null( dmrd.lims ) ) {
        dmrd.xy.max <- max( abs( c( dmrd.x, dmrd.y ) ) )
        dmrd.lims <- dmrd.xy.max * c( -1, 1 )
    }
    
    if ( is.null( point.size ) )
        point.size <- 0.5
    
    dmrd.plot <- dmrd.plot + 
        coord_fixed() + 
        lims( x = dmrd.lims, y = dmrd.lims ) + 
        geom_point( shape = 20, stroke = 0, size = point.size ) + 
        theme_bw() + 
        theme( axis.title = element_blank(), 
            axis.text  = element_blank(), 
            axis.ticks = element_blank(), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank() )
    
    if ( ! is.null( event.group ) )
    {
        dmrd.plot <- dmrd.plot + 
            facet_wrap( vars( event.group ), 
                labeller = as_labeller( group.label ), 
                nrow = dmrd.nrow, ncol = dmrd.ncol ) + 
            theme( strip.background = element_rect( fill = "white" ), 
                strip.text = element_text( size = fcs.dmrd.group.title.size ) )
    }
    
    dmrd.plot
}

