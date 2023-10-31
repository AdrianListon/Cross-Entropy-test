# User-interactive portion of the script
# set group names
if (Be.Chatty == TRUE){
  cat( 
"This part of the workflow requires your input.\n
Let's start by defining the groups in your experiment.
For this to work, the data files need to be named with the identifying
label for the group, and those tags need to be unique to each group.
Please type the tags (names) exactly as they appear in the files.
In the following step, you'll get the opportunity to assign new names
for each group, which will appear on the plots in the end.
\n")
  
  Sys.sleep(message.delay.time)
}

number.of.groups <- readline( "How many groups do you have? ")
number.of.groups <- as.numeric(number.of.groups)

input.group.names <- character()

for (group.n in 1:number.of.groups) {
  input.group.names[group.n] <- readline( paste0("Please enter the identifier label for Group ", group.n, ": " ))
}

fcs.condition <- unlist(strsplit(input.group.names, split = ","))

# check group name entry, re-name groups as needed
groups.to.correct <- multi.menu( fcs.condition, title = "Do you need to correct any of the group names? ")
groups.to.correct <- as.numeric(groups.to.correct)

while (groups.to.correct != 0) {
  for (group.n in groups.to.correct) {
    input.group.names[group.n] <- readline( paste0("Please enter the identifier label for Group ", group.n, ": " ))
    fcs.condition <- unlist(strsplit(input.group.names, split = ","))
    groups.to.correct <- multi.menu( fcs.condition, title = "Do you need to correct any of the group names? ")
    groups.to.correct <- as.numeric(groups.to.correct)
  }
}

want.to.rename.groups <- menu(c("Yes", "No"), 
  title = paste("Do you want to enter different names for the groups? These labels will appear on the plots. "))


fcs.group.names <- NULL
fcs.group.names$orig.name  <- fcs.condition
fcs.group.names$new.name  <- fcs.condition
new.group.labels <- fcs.condition

if (want.to.rename.groups == 1) {
  
  for (cond in 1:length(fcs.condition)) {
    new.group.labels[cond] <- readline( paste("Please enter the label you want for ", 
                                              fcs.condition[cond], " ", sep = ""))
  }
  
  fcs.condition.label <- new.group.labels
  names(fcs.condition.label) <- fcs.condition
  cat("\n")
  cat("\nYour groups will be labeled as follows: \n")
  print(fcs.condition.label)
  cat("\n")
  cat("\n")
  Sys.sleep(message.delay.time)
  
} else {
  
  fcs.condition.label <- fcs.condition
  names(fcs.condition.label) <- fcs.condition
  
  cat("\n")
  cat("\nYour groups will be labeled as follows: \n")
  print(fcs.condition.label)
  cat("\n")
  cat("\n")
  Sys.sleep(message.delay.time)
  
}



# select whether files are csv or fcs----------

if (Be.Chatty == TRUE){
  cat(
"\n
Flowcytoscript will accept either FCS or CSV files as input.

If you use FCS files, a biexponential transformation will be applied automatically 
so that the data should be correctly scaled. This transformation will be adapted
to the type of cytometer you have used for acquiring the data.

To have more control over the transformation (scaling) of the data, you can export
data as CSV files with the biexponential transformation already embedded
in the data. To create these CSV-channel-value files, see the instructions.
\n"
  )
  
  Sys.sleep(message.delay.time)
}

input.file.type <- menu(c("CSV", "FCS"), title = "Please select whether you are using CSV or FCS files.")


if (input.file.type == 2){
  data.input.files <- data.directory.fcs.files
  file.extension <- "\\.fcs$"
  
  flowFrame <- read.FCS(list.files(fcs.data.dir, "\\.fcs$", full.names = TRUE)[1], 
                        which.lines = 1:100, truncate_max_range = FALSE)
  
  channels <- data.frame(name = unname(pData(parameters(flowFrame))$name), 
                         desc = unname(pData(parameters(flowFrame))$desc))
  # filter out height and width channels
  descs.filtered <- channels$desc[!is.na(channels$desc) & channels$desc != '-'] 
  channels.filtered <- dplyr::filter(channels, desc %in% descs.filtered)
  area.only <- c("Height", "height", "Width", "width")
  area.only <- channels$desc[]
  dplyr::filter(channels, desc %in% descs.filtered)
  
  # if data are from the ID7000, truncate channel names to marker only, remove Height and Width
  ## note: this will probably be needed for ZE5 as well
  if( flowFrame@description$`$CYT` == "ID7000" ){
    non.area <- c("Height", "height", "Width", "width")
    area.channels <- channels$desc[!grepl(paste(non.area, collapse = "|"), channels$desc)]
    channels.filtered <- dplyr::filter(channels, desc %in% area.channels)
    
    id7000.pattern <- "\\S"
    
    lapply( channels.filtered$desc, strsplit, "\\s" )
    split.channel.names <- strsplit(channels.filtered$desc, "\\s")
    
    for ( channel in 1:length(split.channel.names) ){
      channels.filtered$desc[channel] <- split.channel.names[[channel]][1]
    }
  }
  
  selected.channels <- NULL
  
}else {
  
  data.input.files <- data.directory.csv.files
  file.extension <- "\\.csv$"

  setDTthreads( threads = fcs.tsne.threads.n )
  
  channel.selection.data <- fread( list.files(fcs.data.dir, "\\.csv$", full.names = TRUE)[1], check.names = TRUE )
  channels <- data.frame( desc = unname( colnames( channel.selection.data ) ))
  descs.filtered <- channels$desc[!is.na(channels$desc) & channels$desc != '-'] 
  channels.filtered <- dplyr::filter(channels, desc %in% descs.filtered)

  selected.channels <- NULL
}

# select channels for analysis---------------

if (Be.Chatty == TRUE){
  cat(
"\n
Now you'll need to select the markers (channels) you want to use for your analysis.
Enter the numbers of the channels you want. You may need to expand the console window
in order to see everything.\n
\n"
  )
  
  Sys.sleep(message.delay.time)
}

reselect.channels <- 1
while( reselect.channels == 1 ){
  selected.channels$orig.name <- channels.filtered$desc[multi.menu(channels.filtered$desc, 
                                                                   title = "Please select channels for analysis:", header = c("Channel", "Marker"))]
  # set as list
  fcs.channel <- selected.channels$orig.name
  
  cat("\nThese are the channels you've selected:\n")
  cat( paste(fcs.channel, collapse = "\n") )
  cat("\n
      ")
  
  reselect.channels <- menu(c("Yes", "No"), 
                            title = "Do you need to change your channel selection? ")
}


if (Be.Chatty == TRUE){
  cat(
    "
Next you'll have the option to rename the marker labels.
\n"
  )
  Sys.sleep(message.delay.time)
  
}

# rename channels
renaming.positions <- multi.menu(selected.channels$orig.name, 
                                 title = "To rename any channels, select them now")
channels.to.rename <- selected.channels$orig.name[renaming.positions]

new.channel.names <- selected.channels$orig.name

if ( length(renaming.positions) > 1 ){
  for (ch in renaming.positions) {
    new.channel.names[ch] <- readline( paste0("Please enter the new channel name for ", 
                                             selected.channels$orig.name[ch], " : " ) )
  }
} else if (renaming.positions !=0){
  for (ch in renaming.positions) {
    new.channel.names[ch] <- readline( paste0("Please enter the new channel name for ", 
                                             selected.channels$orig.name[ch], " : " ) )
  }
}

fcs.channel.label <- new.channel.names
names(fcs.channel.label) <- selected.channels$orig.name

if (Be.Chatty == TRUE){
  cat("\nThis is how your channels will be labeled: \n")
  print(fcs.channel.label)
  cat("\n")
  
  Sys.sleep(message.delay.time)
  
  cat(
"\n
Now we'll match the data files to the group names you entered earlier.\n
\n"
  )
  
  Sys.sleep(message.delay.time)
}


fcs.condition.n <- length( fcs.condition )
fcs.channel.n <- length( fcs.channel )

# set required variables-----------

fcs.seed.base <- as.numeric(paste( fcs.condition.n, fcs.channel.n, fcs.tsne.threads.n, 
                                   length(data.input.files), sep = ""))

fcs.use.cached.results <- TRUE

fcs.sample.number.width <- 2
fcs.event.number.width <- 6

source( file.path( fcs.src.dir, "flowcytoscript_graphics_parameters.r" ) )



# read in data-----------

flow.data.filename.all <- list.files( fcs.data.dir, file.extension )

flow.data.filename <- grep( paste0( fcs.condition, collapse = "|" ), 
                            flow.data.filename.all, value = TRUE )

sample.name.format <- paste0( "%s.%0", fcs.sample.number.width, "d" )
event.name.format <- paste0( "%s.%0", fcs.event.number.width, "d" )

flow.data.filename.sample <- rep( "", length( flow.data.filename ) )
names( flow.data.filename.sample ) <- flow.data.filename

sample.idx.next <- rep( 1, fcs.condition.n )
names( sample.idx.next ) <- fcs.condition

if (input.file.type==2){
  fcs.channel <- channels.filtered$name[ channels.filtered$desc %in% fcs.channel ]
  names(fcs.channel) <- fcs.channel.label
  flow.data <- lapply( flow.data.filename, function( flow.data.fn ) {

    sample.flow.frame <- read.FCS( file.path( fcs.data.dir, flow.data.fn ), 
                                   transformation = NULL, truncate_max_range = FALSE )
    
    condition <- fcs.condition[ sapply( fcs.condition, grepl, flow.data.fn ) ]
    stopifnot( length( condition ) == 1 )
    
    sample.data <- exprs( sample.flow.frame )
    
    if ( ! all( fcs.channel %in% colnames( sample.data ) ) )
    {
      cat( sprintf( "File: %s\n", flow.data.fn ) )
      print( sort( fcs.channel[ 
        ! fcs.channel %in% colnames( sample.data ) ] ) )
      print( sort( colnames( sample.data ) ) )
      stop( "mismatch in names of fcs channels" )
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
  
}else{
  flow.data <- lapply( flow.data.filename, function( flow.data.fn ) {
    
    sample.data <- as.matrix(fread( file.path( fcs.data.dir, flow.data.fn ), check.names = TRUE ))
    
    condition <- fcs.condition[ sapply( fcs.condition, grepl, flow.data.fn ) ]
    stopifnot( length( condition ) == 1 )
    
    if ( ! all( fcs.channel %in% colnames( sample.data ) ) )
    {
      cat( sprintf( "File: %s\n", flow.data.fn ) )
      print( sort( fcs.channel[ 
        ! fcs.channel %in% colnames( sample.data ) ] ) )
      print( sort( colnames( sample.data ) ) )
      cat("Channel mismatch error\n
        Please check that your files were all run with the same flow panel and try again.\n")
      Sys.sleep(message.delay.time*2)
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
}

flow.data <- do.call( rbind, flow.data )

if (input.file.type==2){
  colnames(flow.data) <- fcs.channel.label
}

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
names( flow.event.condition.n ) <- fcs.condition.label

cat("\n
Files per group:
\n")
print(table( flow.sample.condition ))

cat("\n
Events per group:
\n")
print(flow.event.condition.n)
cat("\n
Events per sample:
\n")
print(flow.event.sample.n)

if (Be.Chatty == TRUE){
  cat(
"\n
We found these data files matching your groups.
If this doesn't meet your expectations, you should start over
and double-check your file names vis-a-vis your group names.\n"
  )
  
  Sys.sleep(message.delay.time)
  
  
  
  # set downsampling-------------
  
  cat(
"\n
Please set the number of cells (events) you'd like to use for the analysis.
This will be set as a maximum number per file, so if you set it at 2000
but you only have 500 in some samples, all 500 will be used.
The more data you analyze, the longer it will take. If you aren't sure,
maybe try for a total of no more than 100,000 (for example, 2 groups
with 5 samples per group and 10000 cells/sample gives 100000 total.)
Please enter the number without punctuation.\n"
  )
  
  Sys.sleep(message.delay.time)
  cat("\n")
  
}

cat("For your analysis, please enter a maximum number of cells you'd like to 
analyze per sample. For samples with fewer cells than this number, all 
cells will be used.
    ")

fcs.dmrd.data.sample.n.per.sample <- readline( "Set downsampling number. To run all the cells without downsampling, enter 0: ")
fcs.dmrd.data.sample.n.per.sample <- as.numeric(fcs.dmrd.data.sample.n.per.sample)

if (fcs.dmrd.data.sample.n.per.sample==0){
  fcs.dmrd.data.sample.n.per.sample <- NULL
}

# plotting options------
plot.everything <- menu(c("Yes", "No"), 
                        title = "\nDo you want overlays of every marker on your tSNE and UMAP projections as well as plotting clusters?
Generating many plots can be slow with lots of cells.\n")

# pick conditions for T-REX---------

if( length(fcs.condition) > 2 ){
  trex.selection <- multi.menu(fcs.condition.label, 
                               title = "Please select two groups for the T-REX analysis of under- and over-represented regions")
  trex.condition <- fcs.condition[trex.selection]
} else {
  trex.selection <- c(1,2)
  trex.condition <- fcs.condition
}


# run crossentropy test? ------------
run.crossentropy <- menu( c("Yes", "No"), 
                          title = "\nDo you want to run the crossentropy statistical test on your tSNE and UMAP projections?
It can be slow if you have lots of events, but is a powerful tool.\n")

# define model system-------

experimental.system <- readline(prompt = "Please enter a name for your experiment: ")

# define figure parameters for samples---------------

cat(
  "\n
  Setting color palette...\n
  \n"
)

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


# create dirs---------------
cat("\n
    Creating output folders...\n")

figure.dir <- c( 
  fcs.ce.diff.tsne.figure.dir, 
  fcs.ce.diff.umap.figure.dir, 
  fcs.density.figure.dir, 
  fcs.heatmap.figure.dir, 
  fcs.histogram.figure.dir, 
  fcs.tsne.figure.dir, 
  fcs.umap.figure.dir,
  trex.figure.dir,
  fcs.pca.figure.dir,
  fcs.mfi.stats.dir,
  fcs.cluster.stats.dir
)

table.dir <- fcs.cluster.table.dir

for ( the.dir in c( figure.dir, table.dir ) )
  if ( ! file.exists( the.dir ) )
    dir.create( the.dir, recursive = TRUE )


# select data for dimensionality reduction---------------

cat("\n
    Selecting data for analysis...\n")

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
dmrd.data.n <- nrow(dmrd.data)
dmrd.event.sample <- flow.event.sample[ dmrd.data.idx ]
dmrd.event.condition <- flow.event.condition[ dmrd.data.idx ]
dmrd.event.n <- nrow(dmrd.data)

# transform data if using fcs files-------
## TBD: include transforms for FACSDiscover, ZE5, Attune, Bigfoot
if ( input.file.type == 2 ){
  
  if( flowFrame@description$`$CYT` == "Aurora" ){
    width.basis <- -1000
    max.value <- 4194303
    log.decades <- 5.5
  } else if( flowFrame@description$`$CYT` == "ID7000" ){
    width.basis <- -500
    max.value <- 1000000
    log.decades <- 5
  } else {
    width.basis <- -100
    max.value <- 262144
    log.decades <- 4.5
  }
  
  extra.neg.decades <- 0
  
  biexp.transform <- flowjo_biexp(channelRange = 1250, maxValue = max.value, 
                                  pos = log.decades, neg = extra.neg.decades,
                                  widthBasis = width.basis)
  
  # transform data
  dmrd.data.untransformed <- dmrd.data
  dmrd.data <- apply( dmrd.data, 2, biexp.transform)
  dmrd.data <- apply( dmrd.data, 2, FUN = "-", 250 )
  dmrd.data.long <- pivot_longer(data.frame( dmrd.data ), cols = everything(), 
                                 names_to = "parameter", values_to = "value")
  rownames(dmrd.data) <- rownames(dmrd.data.untransformed)
  
  #plot histograms of channels
  transformation.plot <- ggplot( dmrd.data.long, 
                                 aes(x = value, y = after_stat(count) ))+
    geom_density(fill='black', alpha = 0.4) +
    theme_classic()+
    facet_wrap(~parameter, scales = "free")+
    coord_cartesian(xlim = c(-50,1000))+
    xlab("Channel")
  
  ggsave(
    file.path( fcs.density.figure.dir, 
               "biexponential_transform_fcs.jpg" ), 
    transformation.plot, 
    width = fcs.density.figure.width.base * ( fcs.channel.n + 1 )*1.2, 
    height = fcs.density.figure.height.base *100
  )
}


# end---------------
cat(
  "\n
  Move to the next section.\n
  \n"
)


