# Optional add-on to flowcytoscript analysis
# Use T-REX approach from Cytolab to determine changing knn regions

# first run tSNE, UMAP and/or EmbedSOM using flowcytoscript

library(FNN)

## Set things up-------------
# pick two groups for identifying differences from the output below
# copy into trex.condition <- c() between parenthesis
print(fcs.condition)
trex.condition <- c("Spleen", "siLPL")


# subset data based on those groups
trex.dmrd.data <- subset(dmrd.data, dmrd.event.condition %in% trex.condition)
trex.umap.data <- umap.data
rownames(trex.umap.data) <- rownames(dmrd.data)
trex.umap.data <- subset(trex.umap.data, dmrd.event.condition %in% trex.condition)
trex.tsne.data <- tsne.data
rownames(trex.tsne.data) <- rownames(dmrd.data)
trex.tsne.data <- subset(trex.tsne.data, dmrd.event.condition %in% trex.condition)
trex.embed.som <- embed.som
rownames(trex.embed.som) <- rownames(flow.data)
trex.embed.som <- subset(trex.embed.som, flow.event.condition %in% trex.condition)
kvalue <- 60


## Now you may generate plots showing regions that are disproportionately in one condition or the other.
## Use the three sections below to generate these with either tSNE, UMAP or EmbedSOM plots


## Using tSNE as the plot--------------------------
# KNN search per cell 
tsne.neighbor.index <- knnx.index(trex.tsne.data,trex.tsne.data,k=kvalue)
first.condition.length <- nrow(subset(dmrd.data, dmrd.event.condition == trex.condition[1]))
tsne.neighbor.index[tsne.neighbor.index <= first.condition.length] <- 0
tsne.neighbor.index[tsne.neighbor.index > first.condition.length] <- 1

# calculate percent change in each KNN region
percent.change.tsne <- (rowSums(tsne.neighbor.index) / kvalue * 100)

# binning and plot info
range <- apply(apply( trex.tsne.data, 2, range), 2, diff)
graphical.ratio <- (range[1] / range[2])
test.round <- round(percent.change.tsne)
trex.plot.tsne <-
  data.frame(x = trex.tsne.data[, 1], y = trex.tsne.data[, 2], col = test.round)
trex.plot.tsne$cuts <- cut(trex.plot.tsne$col, c(0, 5, 15, 85, 95, 100), include.lowest = TRUE, right = FALSE)
trex.plot.tsne$cuts <- factor(trex.plot.tsne$cuts,
                         levels = c("[15,85)", "[5,15)", "[0,5)", "[85,95)", "[95,100]"))
ordered.plot.tsne <- trex.plot.tsne[order(trex.plot.tsne$cuts), ]


# create T-REX plot

trex_comparison <- paste(trex.condition[1], " vs ", trex.condition[2], sep = "") 
dir.create("figure_trex_tsne")
png(
  paste(
    "./figure_trex_tsne/",
    strftime(Sys.time(), "%Y-%m-%d_%H%M%S"),
    " tSNE TREX plot.png",
    sep = ""
  ),
  res = 200,
  width = 1500,
  height = 1500
)
final.trex.plot.tsne <-
  ggplot(ordered.plot.tsne) + geom_point(aes(x = x, y = y, colour = cuts), cex = 1) +
  scale_color_manual(
    name = "Ratio",
    values = c(
      "[15,85)" = "lightgray",
      "[5,15)" = "lightskyblue",
      "[0,5)" = "navyblue",
      "[85,95)" = "lightcoral",
      "[95,100]" = "darkred"
    )
  ) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank()) +
  labs (x = "tSNE x", y = "tSNE y", title = paste( trex_comparison," - Percent Change",sep = "")) + 
  coord_fixed(ratio = graphical.ratio)+
  guides(color = guide_legend(override.aes = list(size=3)))
print(final.trex.plot.tsne)
dev.off()

print(final.trex.plot.tsne)


## Using UMAP as the plot--------------------------
# KNN search per cell 
umap.neighbor.index <- knnx.index(trex.umap.data,trex.umap.data,k=kvalue)
first.condition.length <- nrow(subset(dmrd.data, dmrd.event.condition == trex.condition[1]))
umap.neighbor.index[umap.neighbor.index <= first.condition.length] <- 0
umap.neighbor.index[umap.neighbor.index > first.condition.length] <- 1

# calculate percent change in each KNN region
percent.change.umap <- (rowSums(umap.neighbor.index) / kvalue * 100)

# binning and plot info
range <- apply(apply( trex.umap.data, 2, range), 2, diff)
graphical.ratio <- (range[1] / range[2])
test.round <- round(percent.change.umap)
trex.plot.umap <-
  data.frame(x = trex.umap.data[, 1], y = trex.umap.data[, 2], col = test.round)
trex.plot.umap$cuts <- cut(trex.plot.umap$col, c(0, 5, 15, 85, 95, 100), include.lowest = TRUE, right = FALSE)
trex.plot.umap$cuts <- factor(trex.plot.umap$cuts,
                              levels = c("[15,85)", "[5,15)", "[0,5)", "[85,95)", "[95,100]"))
ordered.plot.umap <- trex.plot.umap[order(trex.plot.umap$cuts), ]
range <- apply(apply(trex.umap.data, 2, range), 2, diff)
graphical.ratio <- (range[1] / range[2])

# create T-REX plot

trex_comparison <- paste(trex.condition[1], " vs ", trex.condition[2], sep = "") 
dir.create("figure_trex_umap")
png(
  paste(
    "./figure_trex_umap/",
    strftime(Sys.time(), "%Y-%m-%d_%H%M%S"),
    " UMAP TREX plot.png",
    sep = ""
  ),
  res = 200,
  width = 1500,
  height = 1500
)
final.trex.plot.umap <-
  ggplot(ordered.plot.umap) + geom_point(aes(x = x, y = y, colour = cuts), cex = 1) +
  scale_color_manual(
    name = "Ratio",
    values = c(
      "[15,85)" = "lightgray",
      "[5,15)" = "lightskyblue",
      "[0,5)" = "navyblue",
      "[85,95)" = "lightcoral",
      "[95,100]" = "darkred"
    )
  ) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank()) +
  labs (x = "UMAP1", y = "UMAP2", title = paste( trex_comparison," - Percent Change",sep = "")) + 
  coord_fixed(ratio = graphical.ratio)+
  guides(color = guide_legend(override.aes = list(size=3)))
print(final.trex.plot.umap)
dev.off()

print(final.trex.plot.umap)



## Using EmbedSOM as the plot--------------------------
# KNN search per cell 
som.neighbor.index <- knnx.index(trex.embed.som,trex.embed.som,k=kvalue)
first.condition.length <- nrow(subset(flow.data, flow.event.condition == trex.condition[1]))
som.neighbor.index[som.neighbor.index <= first.condition.length] <- 0
som.neighbor.index[som.neighbor.index > first.condition.length] <- 1

# calculate percent change in each KNN region
percent.change.som <- (rowSums(som.neighbor.index) / kvalue * 100)

# binning and plot info
range <- apply(apply( trex.embed.som, 2, range), 2, diff)
graphical.ratio <- (range[1] / range[2])
test.round <- round(percent.change.som)
trex.plot.som <-
  data.frame(x = trex.embed.som[, 1], y = trex.embed.som[, 2], col = test.round)
trex.plot.som$cuts <- cut(trex.plot.som$col, c(0, 5, 15, 85, 95, 100), include.lowest = TRUE, right = FALSE)
trex.plot.som$cuts <- factor(trex.plot.som$cuts,
                              levels = c("[15,85)", "[5,15)", "[0,5)", "[85,95)", "[95,100]"))
ordered.plot.som <- trex.plot.som[order(trex.plot.som$cuts), ]
range <- apply(apply(trex.embed.som, 2, range), 2, diff)
graphical.ratio <- (range[1] / range[2])

# create T-REX plot

trex_comparison <- paste(trex.condition[1], " vs ", trex.condition[2], sep = "") 
dir.create("figure_trex_embedsom")
png(
  paste(
    "./figure_trex_embedsom/",
    strftime(Sys.time(), "%Y-%m-%d_%H%M%S"),
    " EmbedSOM TREX plot.png",
    sep = ""
  ),
  res = 200,
  width = 1500,
  height = 1500
)
final.trex.plot.som <-
  ggplot(ordered.plot.som) + geom_point(aes(x = x, y = y, colour = cuts), cex = 1) +
  scale_color_manual(
    name = "Ratio",
    values = c(
      "[15,85)" = "lightgray",
      "[5,15)" = "lightskyblue",
      "[0,5)" = "navyblue",
      "[85,95)" = "lightcoral",
      "[95,100]" = "darkred"
    )
  ) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank()) +
  labs (x = "EmbedSOM1", y = "EmbedSOM2", title = paste( trex_comparison," - Percent Change",sep = "")) + 
  coord_fixed(ratio = graphical.ratio)+
  guides(color = guide_legend(override.aes = list(size=3)))
print(final.trex.plot.som)
dev.off()

print(final.trex.plot.som)