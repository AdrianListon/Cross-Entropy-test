# check for files

data.directory.csv.files <- list.files(fcs.data.dir, "\\.csv$", full.names = TRUE)
data.directory.fcs.files <- list.files(fcs.data.dir, "\\.fcs$", full.names = TRUE)


if ( length(c(data.directory.csv.files, data.directory.fcs.files)) < 1){
  cat("There don't seem to be any data files in your Data folder.\n
Please create a folder for your analysis.
In that folder, place a copy of this script, and also a folder called 
Data containing exported flow data in either csv-channel-values or fcs format,
close Rstudio, and re-start Rstudio by double clicking on the script file.\n
See the instructions for more detail.\n")
}

stopifnot("No data files detected" = length(c(data.directory.csv.files, data.directory.fcs.files)) > 0)

# start packages
cat(
  "\n
  Loading packages...\n
  \n"
)

library( digest )
require( dunn.test )
library( ggplot2 )
library( ggridges )
library( ggrepel )
library( ggtext )
library( RColorBrewer )
library( Rtsne )
library( uwot )
library( dplyr )
library( tidyr )
library( FastPG )
require( RcppHNSW )
library( parallel )
library( readxl )
library( coda )
library( emmeans )
library( ConsensusClusterPlus )
library( EmbedSOM )
library( flowCore )
library( flowWorkspace )
library( data.table )


# ncores check
fcs.tsne.threads.n <- parallel::detectCores()

if (fcs.tsne.threads.n < 4){
  cat("You don't seem to have much processing power.\n
You might consider using a faster machine for any large datasets.")
  Sys.sleep(message.delay.time)
}


# function to set random seed depending on base number and string---------------

set.seed.here <- function( seed.base, seed.char )
{
  seed.add <- strtoi( substr( digest( seed.char, "xxhash32" ), 2, 8 ), 16 )
  seed.new <- seed.base + seed.add    
  set.seed( seed.new )
  invisible( seed.new )
}

# source parameters and functions---------------

source( file.path( fcs.src.dir, "ce_diff_test.r" ) )
source( file.path( fcs.src.dir, "ce_diff_test_tsne.r" ) )
source( file.path( fcs.src.dir, "ce_diff_test_umap.r" ) )
source( file.path( fcs.src.dir, "plot_cluster_dmrd.r" ) )
source( file.path( fcs.src.dir, "plot_all_dmrd_figures.r" ) )
source( file.path( fcs.src.dir, "plot_dimensionality_reduction.r" ) )
source( file.path( fcs.src.dir, "plot_changed_regions.r" ) )

if (Be.Chatty == TRUE){
  cat(
"\n
If there are no error messages, then the packages are loaded and source code has been located.\n
Warnings about packages being built under a slight different R version are usually not a problem.\n
Proceed to the next step.\n
\n"
  )
} else {
  cat(
"Move to the next code chunk"
  )
}

