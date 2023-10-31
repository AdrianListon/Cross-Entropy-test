# set up for simplified flowcytoscript

# Welcome and check R status
cat("This script will try to help you update R\n
      and also install Rtools.\n
    \n")
Sys.sleep(2)
cat(
  "If you can do this yourself, that may work better, particularly\n
  for non-Windows users.\n
  \n"
)

install.packages('installr')
suppressPackageStartupMessages( library(installr) )
updateR()
install.Rtools()

cat(
  "If you see error messages, try to manually install Rtools (Windows)\n
  https://cran.r-project.org/bin/windows/Rtools/  \n
  or Command line tools (Mac)\n
  https://clanfear.github.io/CSSS508/docs/compiling.html  \n
  \n"
)
Sys.sleep(2)
cat(
  "You'll want to be using R version 4 or higher.\n
  \n"
)
Sys.sleep(2)
cat(
  "And for faster processing on a Mac, you'll want OpenMP.\n
  https://mac.r-project.org/openmp/ \n
  \n"
)

