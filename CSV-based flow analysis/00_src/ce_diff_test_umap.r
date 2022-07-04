
# calculates cross-entropy for umap plots and calls ce.diff.test


ce.diff.test.umap <- function( 
    orig.dist, orig.knn, umap.data, umap.param, 
    event.partition, 
    partition.label = NULL, partition.color = NULL, partition.line.type = NULL, 
    base.test = "ks", base.dist = "ks", 
    prob.sample.n = NULL, dendrogram.order.weight = NULL, 
    result = NULL, cdf.figure = NULL, dendrogram.figure = NULL 
)
{
    stopifnot( nrow( orig.dist ) == nrow( umap.data ) && 
        nrow( orig.dist ) == length( event.partition ) )
    
    data.n <- nrow( orig.dist )
    
    if ( ! is.null( prob.sample.n ) && prob.sample.n < data.n )
        prob.sample.idx <- sample( data.n, prob.sample.n )
    else
        prob.sample.idx <- 1 : data.n
    
    if ( fcs.use.cached.results && 
            file.exists( fcs.ce.diff.umap.cache.file.path ) )
    {
        cat( "Using cached results for probability\n" )
        
        load( fcs.ce.diff.umap.cache.file.path )
    }
    else
    {
        cat( "Calculating probability\n" )
        
        orig.umap.prob <- calculate.probability.umap( orig.dist, orig.knn, 
            umap.data, umap.param )
        
        save( orig.umap.prob, file = fcs.ce.diff.umap.cache.file.path )
    }
    
    cross.entropy.all <- calculate.fuzzy.cross.entropy( 
        orig.umap.prob$orig[ prob.sample.idx, ], 
        orig.umap.prob$umap[ prob.sample.idx, ] 
    )
    
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


calculate.probability.umap <- function( orig.dis, orig.kn, umap.dat, 
    umap.param )
{
    orig.dat.n <- nrow( orig.dis )
    umap.dat.n <- nrow( umap.dat )
    
    stopifnot( orig.dat.n == umap.dat.n )
    
    # get nearest neighbors in original space and their distances
    
    orig.self.idx <- sapply( 1 : orig.dat.n, function( ri ) {
        ri.idx <- which( orig.kn[ ri, ] == ri )
        ifelse( length( ri.idx ) == 1, ri.idx, NA )
    } )
    
    stopifnot( ! is.na( orig.self.idx ) )
    
    orig.neigh <- t( sapply( 1 : orig.dat.n, function( ri ) 
        orig.kn[ ri, - orig.self.idx[ ri ] ] ) )
    
    orig.dis.reduc <- t( sapply( 1 : orig.dat.n, function( ri ) 
        orig.dis[ ri, - orig.self.idx[ ri ] ] ) )
    
    # calculate probabilities associated to distances in original space
    
    umap.sigma <- apply( orig.dis.reduc, 1, function( dd ) {
        umap.sigma.error <- function( ss, dd ) {
            p <- exp( - pmax( 0, dd - min( dd ) ) / ss )
            sum( p ) - log2( length( p ) )
        }
        
        dd.ascen <- sort( dd )
        dd.ascen <- dd.ascen[ dd.ascen > 0 ]
        ss.lower <- dd.ascen[ 1 ]
        
        dd.descen <- sort( dd, decreasing = TRUE )
        dd.descen <- dd.descen[ ! is.infinite( dd.descen ) ]
        ss.upper <- dd.descen[ 1 ]
        
        while( umap.sigma.error( ss.upper, dd ) < 0 )
        {
            ss.lower <- ss.upper
            ss.upper <- 2 * ss.upper
        }
        
        while( umap.sigma.error( ss.lower, dd ) > 0 )
        {
            ss.upper <- ss.lower
            ss.lower <- ss.lower / 2
        }
        
        uniroot( umap.sigma.error, dd, 
            interval = c( ss.lower, ss.upper ), 
            tol = ( ss.upper - ss.lower ) * .Machine$double.eps^0.25 )$root
    } )
    
    orig.prob <- t( sapply( 1 : orig.dat.n, function( i ) 
        exp( - pmax( 0, orig.dis.reduc[ i, ] - min( orig.dis.reduc[ i, ] ) ) / 
            umap.sigma[ i ] )
    ) )
    
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
                    sym.prob <- orig.prob[ i, j2 ] + orig.prob[ j, i2 ] - 
                        orig.prob[ i, j2 ] * orig.prob[ j, i2 ]
                    orig.prob[ i, j2 ] <- sym.prob
                    orig.prob[ j, i2 ] <- sym.prob
                }
            }
        }
    
    # get distances in umap space for closest neighbors in original space
    
    umap.dist2 <- t( sapply( 1 : umap.dat.n, function( i )
        sapply( orig.neigh[ i, ], function( j )
            sum( ( umap.dat[ i, ] - umap.dat[ j, ] )^2 )
        )
    ) )
    
    # calculate probabilities associated to distances in umap representation
    
    umap.a <- umap.param$a
    umap.b <- umap.param$b
    
    umap.prob <- t( apply( umap.dist2, 1, function( dd2 ) 
        p <- 1 / ( 1 + umap.a * dd2 ^ umap.b )
    ) )
    
    list( orig = orig.prob, umap = umap.prob )
}


calculate.fuzzy.cross.entropy <- function( prim.prob, secd.prob )
{
    prim.prob.n <- nrow( prim.prob )
    secd.prob.n <- nrow( secd.prob )
    
    prim.prob.m <- ncol( prim.prob )
    secd.prob.m <- ncol( secd.prob )
    
    stopifnot( prim.prob.n == secd.prob.n && prim.prob.m == secd.prob.m )
    
    sapply( 1 : prim.prob.n, function( i ) 
        - sum( prim.prob[ i, ] * log( secd.prob[ i, ] ) + 
            ( 1 - prim.prob[ i, ] ) * log( 1 - secd.prob[ i, ] ) )
    )
}

