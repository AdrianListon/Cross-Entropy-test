# get channels script
# run this to get the channel names for the flowcytoscript parameter file

library(flowCore)
library(dplyr)
library( RColorBrewer )

# source parameters

param.filename <- "./analyze_flow_cytometry_parameter.r"

source( param.filename )

flowFrame <- read.FCS(list.files(fcs.data.dir, "\\.fcs$", full.names = TRUE)[1], truncate_max_range = FALSE)

channels <- data.frame(name = unname(pData(parameters(flowFrame))$name), 
                      desc = unname(pData(parameters(flowFrame))$desc))

for(i in 1:dim(channels)[1] ){
  channels$out1[i] = paste0("\"", channels$name[i], "\", #", channels$desc[i], "\n")
  channels$out2[i] = paste0("\"", channels$name[i], "\" = \"", channels$desc[i], "\",\n") 
  channels$out3[i] = paste0("\"", channels$name[i], "\" = 200, #", channels$desc[i], "\n") 
}

descs.filtered <- channels$desc[!is.na(channels$desc) & channels$desc != '-'] 
channels.filtered <- filter(channels, desc %in% descs.filtered)
cat('"', paste0(channels.filtered$desc, collapse = '","'), '"', sep = "")

# copy the previous output into the channels.of.interest and delete unnecessary markers

channels.of.interest = c("Bcl-6","Ki67","CD95",
                         "Ly-6C","CD127","GATA-3",
                         "CD62L","CXCR5","CD25","CD44","ICOS","RORgT","PD-1","CXCR3",
                         "CXCR4","CD86",
                         "T-bet","IRF4","GL7","CD69")


# add the output to fcs.channel in the parameter file
temp <- channels$out1[channels$desc %in% channels.of.interest]
temp[length(temp)] <- gsub(',', '', temp[length(temp)])
cat(paste(temp, collapse = ""))

# add the output to fcs.channel.label in the parameter file
temp <- channels$out2[channels$desc %in% channels.of.interest]
temp[length(temp)] <- gsub(',', '', temp[length(temp)])
cat(paste(temp, collapse = ""))

# add the output to fcs.channel.asinh.scale in the parameter file
temp <- channels$out3[channels$desc %in% channels.of.interest]
temp[length(temp)] <- gsub(',', '', temp[length(temp)])
cat(paste(temp, collapse = ""))

# save the analyze_flow_cytometry_parameter.r file
# and run the analyze_flow_cytometry.r script