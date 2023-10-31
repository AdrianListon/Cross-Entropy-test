
# plots all dimensionality reduction figures


plot.all.dmrd.figures <- function( 
    redu.data, redu.data.max, 
    redu.figure.lims.factor, redu.figure.point.size, 
    redu.figure.dir, redu.figure.plot, 
    dmrd.data, dmrd.event.cluster, dmrd.event.condition 
)
{
    redu.lims <- redu.figure.lims.factor * redu.data.max * c( -1, 1 )
    
    the.dmrd.figure.width <- fcs.dmrd.figure.width + 
        fcs.dmrd.label.factor.width * 
            max( nchar( fcs.condition.label ), nchar( flow.sample.label ), 
                nchar( fcs.cluster.label ) )
    
    the.dmrd.figure.width.multi <- 
        fcs.dmrd.figure.ncol * fcs.dmrd.figure.width + 
        fcs.dmrd.label.factor.width * 
            max( nchar( fcs.condition.label ), nchar( flow.sample.label ), 
                nchar( fcs.cluster.label ) )
    
    the.dmrd.figure.height <- fcs.dmrd.figure.height
    
    the.dmrd.figure.height.multi <- 
        fcs.dmrd.figure.nrow * fcs.dmrd.figure.height
    
    # plot all events colored by cluster
    
    redu.plot <- plot.dimensionality.reduction( 
        redu.data[ , 1 ], redu.data[ , 2 ], 
        event.partition = dmrd.event.cluster, 
        partition.label = fcs.cluster.label, 
        partition.color = fcs.cluster.color, 
        dmrd.lims = redu.lims, 
        point.size = redu.figure.point.size, 
        show.guide = TRUE 
    )
    
    ggsave( 
        file.path( redu.figure.dir, 
            sprintf( "%s_all_events__cluster.png", redu.figure.plot ) ), 
        redu.plot, 
        width = the.dmrd.figure.width, height = the.dmrd.figure.height 
    )
    
    # plot all events colored by condition
    
    redu.plot <- plot.dimensionality.reduction( 
        redu.data[ , 1 ], redu.data[ , 2 ], 
        event.partition = dmrd.event.condition, 
        partition.label = fcs.condition.label, 
        partition.color = adjustcolor( fcs.condition.color, 
            alpha.f = fcs.dmrd.color.alpha ), 
        dmrd.lims = redu.lims, 
        point.size = redu.figure.point.size, 
        show.guide = TRUE 
    )
    
    ggsave( 
        file.path( redu.figure.dir, 
            sprintf( "%s_all_events__condition.png", redu.figure.plot ) ), 
        redu.plot, 
        width = the.dmrd.figure.width, height = the.dmrd.figure.height 
    )
    
    # plot all events colored by each marker level
    
    for ( fch in names(fcs.channel.label) )
    {
        dmrd.event.level <- dmrd.data[ , fch ]
        
        redu.plot <- plot.dimensionality.reduction( 
            redu.data[ , 1 ], redu.data[ , 2 ], 
            event.level = dmrd.event.level, 
            dmrd.lims = redu.lims, 
            point.size = redu.figure.point.size, 
            show.guide = TRUE, guide.name = fcs.channel.label[ fch ] 
        )
        
        ggsave( 
            file.path( redu.figure.dir, 
                sprintf( "%s_all_events__%s.png", redu.figure.plot, 
                         fcs.channel.label[ fch ] ) ), 
            redu.plot, 
            width = the.dmrd.figure.width, height = the.dmrd.figure.height 
        )
    }
    
    # plot all conditions colored by cluster
    
    redu.plot <- plot.dimensionality.reduction( 
        redu.data[ , 1 ], redu.data[ , 2 ], 
        event.group = dmrd.event.condition, 
        group.label = fcs.condition.label, 
        event.partition = dmrd.event.cluster, 
        partition.label = fcs.cluster.label, 
        partition.color = fcs.cluster.color, 
        dmrd.lims = redu.lims, 
        dmrd.nrow = fcs.dmrd.figure.nrow, dmrd.ncol = fcs.dmrd.figure.ncol, 
        point.size = redu.figure.point.size, 
        show.guide = TRUE 
    )
    
    ggsave( 
        file.path( redu.figure.dir, 
            sprintf( "%s_all_conditions__cluster.png", redu.figure.plot ) ), 
        redu.plot, 
        width = the.dmrd.figure.width.multi, 
        height = the.dmrd.figure.height.multi 
    )
    
    # plot all conditions colored by each marker level
    
    for ( fch in names(fcs.channel.label) )
    {
        dmrd.event.level <- dmrd.data[ , fch ]
        
        redu.plot <- plot.dimensionality.reduction( 
            redu.data[ , 1 ], redu.data[ , 2 ], 
            event.group = dmrd.event.condition, 
            group.label = fcs.condition.label, 
            event.level = dmrd.event.level, 
            dmrd.lims = redu.lims, 
            dmrd.nrow = fcs.dmrd.figure.nrow, dmrd.ncol = fcs.dmrd.figure.ncol, 
            point.size = redu.figure.point.size, 
            show.guide = TRUE, guide.name = fcs.channel.label[ fch ] 
        )
        
        ggsave( 
            file.path( redu.figure.dir, 
                sprintf( "%s_all_conditions__%s.png", redu.figure.plot, 
                         fcs.channel.label[ fch ] ) ), 
            redu.plot, 
            width = the.dmrd.figure.width.multi, 
            height = the.dmrd.figure.height.multi 
        )
    }
    
    # plot each condition colored by cluster
    
    for ( cond in fcs.condition )
    {
        dmrd.cond.idx <- which( dmrd.event.condition == cond )
        
        redu.cond.data <- redu.data[ dmrd.cond.idx, ]
        dmrd.cond.event.cluster <- dmrd.event.cluster[ dmrd.cond.idx ]
        
        redu.plot <- plot.dimensionality.reduction( 
            redu.cond.data[ , 1 ], redu.cond.data[ , 2 ], 
            event.partition = dmrd.cond.event.cluster, 
            partition.label = fcs.cluster.label, 
            partition.color = fcs.cluster.color, 
            dmrd.lims = redu.lims, 
            point.size = redu.figure.point.size, 
            show.guide = TRUE 
        )
        
        ggsave( 
            file.path( redu.figure.dir, 
                sprintf( "%s_%s__cluster.png", redu.figure.plot, 
                    fcs.condition.label[ cond ] ) ), 
            redu.plot, 
            width = the.dmrd.figure.width, height = the.dmrd.figure.height 
        )
    }
    
    # plot each condition colored by sample
    
    for ( cond in fcs.condition )
    {
        dmrd.cond.idx <- which( dmrd.event.condition == cond )
        
        redu.cond.data <- redu.data[ dmrd.cond.idx, ]
        dmrd.cond.event.sample <- dmrd.event.sample[ dmrd.cond.idx ]
        
        redu.plot <- plot.dimensionality.reduction( 
            redu.cond.data[ , 1 ], redu.cond.data[ , 2 ], 
            event.partition = dmrd.cond.event.sample, 
            partition.label = flow.sample.label, 
            partition.color = adjustcolor( flow.sample.color.single, 
                alpha.f = fcs.dmrd.color.alpha ), 
            dmrd.lims = redu.lims, 
            point.size = redu.figure.point.size, 
            show.guide = TRUE 
        )
        
        ggsave( 
            file.path( redu.figure.dir, 
                sprintf( "%s_%s__sample.png", redu.figure.plot, 
                    fcs.condition.label[ cond ] ) ), 
            redu.plot, 
            width = the.dmrd.figure.width, height = the.dmrd.figure.height 
        )
    }
}

