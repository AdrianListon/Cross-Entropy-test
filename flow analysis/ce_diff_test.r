
# calculates test on differences of cross-entropy

ce.diff.test <- function( 
    cross.entropy, 
    event.partition, 
    partition.label, partition.color, partition.line.type, 
    base.test, base.dist, 
    dendrogram.order.weight, 
    result, cdf.figure, dendrogram.figure 
)
{
    partition <- levels( event.partition )
    partition.n <- length( partition )
    
    cross.entropy.split <- split( cross.entropy, event.partition )
    
    if ( base.test == "ks" )
    {
        if ( partition.n == 2 )
        {
            ks.test.res <- ks.test( cross.entropy.split[[ 1 ]], 
                cross.entropy.split[[ 2 ]] )
            
            test.res <- list( ks.single = ks.test.res )
        }
        else if ( partition.n > 2 )
        {
            comparison.n <- partition.n * ( partition.n - 1 ) / 2
            
            ks.pair <- vector( "list", comparison.n )
            comparison <- character( comparison.n )
            D.stat <- numeric( comparison.n )
            p.value <- numeric( comparison.n )

            k <- 1
            
            for ( i in 1 : ( partition.n - 1 ) )
                for ( j in (i+1) : partition.n )
                {
                    ks.pair[[ k ]] <- ks.test( cross.entropy.split[[ i ]], 
                        cross.entropy.split[[ j ]] )
                    
                    comparison[ k ] <- sprintf( "%s - %s", 
                        partition.label[ partition[ i ] ], 
                        partition.label[ partition[ j ] ] )
                    
                    D.stat[ k ] <- ks.pair[[ k ]]$statistic
                    
                    p.value[ k ] <- ks.pair[[ k ]]$p.value
                    
                    k <- k + 1
                }
            
            p.value.adj <- p.adjust( p.value, "holm" )
            
            test.res <- list( ks.multiple = list( ks.pair = ks.pair, 
                comparison = comparison, D.stat = D.stat, 
                p.value = p.value, p.value.adj = p.value.adj ) )
        }
        else
            stop( "no partitions for testing cross-entropy differences" )
    }
    else if ( base.test == "rank" )
    {
        if ( partition.n == 2 )
        {
            wilcox.test.res <- wilcox.test( cross.entropy.split[[ 1 ]], 
                cross.entropy.split[[ 2 ]] )
            
            test.res <- list( rank.single = wilcox.test.res )
        }
        else if ( partition.n > 2 )
        {
            kruskal.test.res <- kruskal.test( cross.entropy, 
                event.partition )
            
            dunn.test.res <- dunn.test( cross.entropy, event.partition, 
                method = "holm", alpha = fcs.ce.diff.test.alpha, altp = TRUE, 
                kw = FALSE, table = FALSE, list = TRUE )
            
            test.res <- list( rank.multiple = list( 
                kruskal = kruskal.test.res, dunn = dunn.test.res ) )
        }
        else
            stop( "no partitions for testing cross-entropy differences" )
    }
    else
        stop( "wrong base test for cross-entropy differences" )
    
    if ( ! is.null( result ) )
    {
        result.file <- file( result, "w" )
        sink( result.file )
        
        tr.name <- names( test.res )
        stopifnot( length( tr.name ) == 1 )
        
        if ( tr.name == "ks.single" )
        {
            cat( "\n**  Kolmogorov-Smirnov test\n")
            
            print( test.res$ks.single )
        }
        else if ( tr.name == "ks.multiple" )
        {
            cat( "\n**  Multiple Kolmogorov-Smirnov tests with Holm correction\n\n")
            
            comparison.width <- max( nchar( test.res$ks.multiple$comparison ) )
            
            for ( i in 1 : length( test.res$ks.multiple$comparison ) )
                cat( sprintf( "%-*s\t\tD = %g\t\tpv = %g\t\tadj-pv = %g\n", 
                    comparison.width, 
                    test.res$ks.multiple$comparison[ i ], 
                    test.res$ks.multiple$D.stat[ i ], 
                    test.res$ks.multiple$p.value[ i ], 
                    test.res$ks.multiple$p.value.adj[ i ] ) )
        }
        else if ( tr.name == "rank.single" )
        {
            cat( "\n**  Wilcoxon rank sum test\n")
            
            print( test.res$rank.single )
        }
        else if ( tr.name == "rank.multiple" )
        {
            cat( "\n**  Kruskal-Wallis rank sum test\n")
            
            print( test.res$rank.multiple$kruskal )
            
            cat( "\n**  Dunn post-hoc test with Holm correction\n\n")
            
            comparison.width <- max( nchar( 
                test.res$rank.multiple$dunn$comparison ) )
            
            for ( i in 1 : length( test.res$rank.multiple$dunn$comparisons ) )
                cat( sprintf( "%-*s\t\tZ = %g\t\tpv = %g\t\tadj-pv = %g\n", 
                    comparison.width, 
                    test.res$rank.multiple$dunn$comparisons[ i ], 
                    test.res$rank.multiple$dunn$Z[ i ], 
                    test.res$rank.multiple$dunn$altP[ i ], 
                    test.res$rank.multiple$dunn$altP.adjusted[ i ] ) )
        }
        else
        {
            sink()
            close( result.file )
            stop( "unknown test in ce-diff result" )
        }
        
        sink()
        close( result.file )
    }
    
    if ( ! is.null( cdf.figure ) )
    {
        if ( is.null( partition.label ) )
            partition.label = partition
        
        if ( is.null( partition.color ) )
            partition.color <- rainbow( partition.n )
        
        if ( is.null( partition.line.type ) )
            partition.line.type <- rep( 1, partition.n )
        
        png( filename = cdf.figure, width = fcs.ce.diff.figure.cdf.width, 
            height = fcs.ce.diff.figure.cdf.height )
        
        par( mar = c( 5.5, 6, 2, 1.5 ) )
        
        plot( ecdf( cross.entropy ), ylim = c( 0, 1 ), 
            xlab = "Cross-entropy", ylab = "CDF", main = "", 
            cex.lab = 3, cex.axis = 2.5, 
            col = fcs.ce.diff.figure.cdf.all.color, 
            lwd = fcs.ce.diff.figure.line.width - 1, do.points = FALSE )
        
        for ( pall in partition )
        {
            ces <- cross.entropy.split[[ pall ]]
            ces.n <- length( ces )
            
            if ( ces.n < fcs.ce.diff.figure.cdf.resolution ) {
                plot( ecdf( ces ), col = partition.color[ pall ], 
                    lty = partition.line.type[ pall ], 
                    lwd = fcs.ce.diff.figure.line.width, 
                    do.points = FALSE, add = TRUE )
            }
            else {
                ecdf.x <- sort( ces )
                ecdf.y <- 1 : ces.n / ces.n
                
                lines( ecdf.x, ecdf.y, col = partition.color[ pall ], 
                    lty = partition.line.type[ pall ], 
                    lwd = fcs.ce.diff.figure.line.width )
            }
        }
        
        legend( "bottomright", 
            legend = c( fcs.ce.diff.figure.cdf.all.label, partition.label ), 
            col = c( fcs.ce.diff.figure.cdf.all.color, partition.color ), 
            lty = partition.line.type, 
            lwd = fcs.ce.diff.figure.line.width, 
            cex = fcs.ce.diff.figure.font.size )
        
        dev.off()
    }
    
    if ( ! is.null( dendrogram.figure ) && partition.n > 2 )
    {
        cross.entropy.dist <- matrix( 0, nrow = partition.n, 
            ncol = partition.n )
        
        for ( i in 1 : ( partition.n - 1 ) )
            for ( j in (i+1) : partition.n )
            {
                if ( base.dist == "ks" )
                    cross.entropy.dist[ i, j ] <- ks.test( 
                        cross.entropy.split[[ i ]], 
                        cross.entropy.split[[ j ]] 
                    )$statistic
                else if ( base.dist == "median" )
                    cross.entropy.dist[ i, j ] <- abs( 
                        median( cross.entropy.split[[ i ]] ) - 
                        median( cross.entropy.split[[ j ]] )
                    )
                else
                    stop( "wrong base dist for cross-entropy differences" )
                
                cross.entropy.dist[ j, i ] <- cross.entropy.dist[ i, j ]
            }
        
        cross.entropy.hclust <- hclust( as.dist( cross.entropy.dist ) )
        
        if ( ! is.null( dendrogram.order.weight ) )
            cross.entropy.hclust <- as.hclust( reorder( 
                as.dendrogram( cross.entropy.hclust ), 
                dendrogram.order.weight, 
                agglo.FUN = mean
            ) )
        
        if ( is.null( partition.label ) )
            partition.label = partition
        
        png( filename = dendrogram.figure, 
            width = fcs.ce.diff.figure.dendrogram.width, 
            height = fcs.ce.diff.figure.dendrogram.height )
        
        par( mar = c( 5, 5.6, 4, 1.4 ) )
        
        plot( cross.entropy.hclust, 
            labels = partition.label, hang = -1, 
            xlab = "", ylab = "", main = "", sub = "", cex.axis = 3, 
            cex = fcs.ce.diff.figure.font.size )
        
        dev.off()
    }
    
    test.res
}

