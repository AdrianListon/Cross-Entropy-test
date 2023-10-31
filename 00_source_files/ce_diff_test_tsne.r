
# calculates cross-entropy for tsne plots and calls ce.diff.test


ce.diff.test.tsne <- function( 
    orig.data, tsne.data, 
    event.partition, 
    partition.label = NULL, partition.color = NULL, partition.line.type = NULL, 
    base.test = "ks", base.dist = "ks", 
    prob.sample.n = NULL, dendrogram.order.weight = NULL, 
    result = NULL, cdf.figure = NULL, dendrogram.figure = NULL 
)
{
    stopifnot( nrow( orig.data ) == nrow( tsne.data ) && 
        nrow( orig.data ) == length( event.partition ) )
    
    data.n <- nrow( orig.data )
    
    if ( ! is.null( prob.sample.n ) && prob.sample.n < data.n )
        prob.sample.idx <- sample( data.n, prob.sample.n )
    else
        prob.sample.idx <- 1 : data.n
    
    if ( fcs.use.cached.results && 
        file.exists( fcs.ce.diff.tsne.cache.file.path ) )
    {
        cat( "Using cached results for probability\n" )
        
        load( fcs.ce.diff.tsne.cache.file.path )
    }
    else
    {
        cat( "Calculating probability\n" )
        
        # sampling here temporary, until optimizing dist( tsne.dat ) below
        orig.tsne.prob <- calculate.probability.tsne( 
            orig.data[ prob.sample.idx, ], 
            tsne.data[ prob.sample.idx, ] 
        )
        
        save( orig.tsne.prob, file = fcs.ce.diff.tsne.cache.file.path )
    }
    
    cross.entropy.all <- calculate.cross.entropy( orig.tsne.prob$orig, 
        orig.tsne.prob$tsne )
    
    event.partition.all <- event.partition[ prob.sample.idx ]
    
    ce.diff.test( 
        cross.entropy.all, 
        event.partition.all, 
        partition.label, partition.color, partition.line.type, 
        base.test, base.dist, 
        dendrogram.order.weight, 
        result, cdf.figure, dendrogram.figure
    )
}


calculate.probability.tsne <- function( orig.dat, tsne.dat )
{
    orig.dat.n <- nrow( orig.dat )
    tsne.dat.n <- nrow( tsne.dat )
    
    stopifnot( orig.dat.n == tsne.dat.n )
    
    # find nearest neighbors in original space and their distances
    
    orig.dat.knn <- RcppHNSW::hnsw_knn( normalize_input( orig.dat ), 
        k = fcs.ce.diff.tsne.perplexity.factor * fcs.tsne.perplexity + 1,
        distance= 'l2', n_threads= fcs.tsne.threads.n, ef = nrow(orig.dat)/10, M = 64, ef_construction = nrow(orig.dat)/10 )
    
    orig.dat.self.idx <- sapply( 1 : orig.dat.n, function( ri ) {
        ri.idx <- which( orig.dat.knn$idx[ ri, ] == ri )
        ifelse( length( ri.idx ) == 1, ri.idx, NA )
    } )
    
    stopifnot( ! is.na( orig.dat.self.idx ) )
    
    orig.neigh <- t( sapply( 1 : orig.dat.n, function( ri ) 
        orig.dat.knn$idx[ ri, - orig.dat.self.idx[ ri ] ] ) )
    
    orig.dist2 <- t( sapply( 1 : orig.dat.n, function( ri ) 
        orig.dat.knn$dist[ ri, - orig.dat.self.idx[ ri ] ]^2 ) )
    
    # calculate probabilities associated to distances in original space
    
    orig.stdev <- apply( orig.dist2, 1, function( dd2 ) {
        tsne.perplexity.error <- function( ss, dd2 ) {
            p <- exp( - dd2 / (2*ss^2) )
            if ( sum( p ) < .Machine$double.eps )
                p <- 1
            p <- p / sum( p )
            p <- p[ p > 0 ]
            2^( - sum( p * log2( p ) ) ) - fcs.tsne.perplexity
        }
        
        dd2.min.idx <- 1
        dd2.ascen <- sort( dd2 )
        while( dd2.ascen[ dd2.min.idx ] == 0 )
            dd2.min.idx <- dd2.min.idx + 1
        ss.lower <- dd2.ascen[ dd2.min.idx ]
        
        dd2.max.idx <- 1
        dd2.descen <- sort( dd2, decreasing = TRUE )
        while( is.infinite( dd2.descen[ dd2.max.idx ] ) )
            dd2.max.idx <- dd2.max.idx + 1
        ss.upper <- dd2.descen[ dd2.max.idx ]
        
        while( tsne.perplexity.error( ss.upper, dd2 ) < 0 )
        {
            ss.lower <- ss.upper
            ss.upper <- 2 * ss.upper
        }
        
        while( tsne.perplexity.error( ss.lower, dd2 ) > 0 )
        {
            ss.upper <- ss.lower
            ss.lower <- ss.lower / 2
        }
        
        uniroot( tsne.perplexity.error, dd2, 
            interval = c( ss.lower, ss.upper ), 
            tol = ( ss.upper - ss.lower ) * .Machine$double.eps^0.25 )$root
    } )
    
    orig.prob <- t( sapply( 1 : orig.dat.n, function( i ) {
        p <- exp( - orig.dist2[ i, ] /  ( 2 * orig.stdev[ i ]^2 ) )
        p / sum( p )
    } ) )
    
    # symmetrize probabilities in original space
    
    for ( i in 1 : orig.dat.n )
        for ( j2 in 1 : length( orig.neigh[ i, ] ) )
        {
            j <- orig.neigh[ i, j2 ]
            
            i2 <- match( i, orig.neigh[ j, ] )
            
            if ( ! is.na( i2 ) )
            {
                if ( j > i )
                {
                    sym.prob <- ( orig.prob[ i, j2 ] + orig.prob[ j, i2 ] ) / 2
                    orig.prob[ i, j2 ] <- sym.prob
                    orig.prob[ j, i2 ] <- sym.prob
                }
            }
            else
                orig.prob[ i, j2 ] <- orig.prob[ i, j2 ] / 2
        }
    
    orig.prob <- sweep( orig.prob, 1, rowSums( orig.prob ), "/" )
    
    # get distances in tsne space for closest neighbors in original space
    
    tsne.dist2 <- t( sapply( 1 : tsne.dat.n, function( i )
        sapply( orig.neigh[ i, ], function( j )
            sum( ( tsne.dat[ i, ] - tsne.dat[ j, ] )^2 )
        )
    ) )
    
    # calculate probabilities associated to distances in tsne representation
    
    tsne.prob.factor <- tsne.dat.n / 
        ( 2 * sum( 1 / ( 1 + dist( tsne.dat )^2 ) ) )
    
    tsne.prob <- t( apply( tsne.dist2, 1, function( dd2 ) 
        p <- tsne.prob.factor / ( 1 + dd2 )
    ) )
    
    list( orig = orig.prob, tsne = tsne.prob )
}


calculate.cross.entropy <- function( prim.prob, secd.prob )
{
    prim.prob.n <- nrow( prim.prob )
    secd.prob.n <- nrow( secd.prob )
    
    prim.prob.m <- ncol( prim.prob )
    secd.prob.m <- ncol( secd.prob )
    
    stopifnot( prim.prob.n == secd.prob.n && prim.prob.m == secd.prob.m )
    
    sapply( 1 : prim.prob.n, function( i ) 
        - sum( prim.prob[ i, ] * log( secd.prob[ i, ] ) )
    )
}

