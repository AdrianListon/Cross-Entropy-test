# plotting changed regions

plot.changed.regions <- function(
    redu.data, percent.change.knn,
    redu.figure.point.size, trex.condition,
    redu.figure.dir, redu.figure.plot
)
{
    # binning and plot info
  range <- apply(apply( redu.data, 2, range), 2, diff)
  graphical.ratio <- (range[1] / range[2])
  test.round <- round(percent.change.knn)
  trex.plot.redu <-
    data.frame(x = redu.data[, 1], y = redu.data[, 2], col = test.round)
  trex.plot.redu$cuts <- cut(trex.plot.redu$col, c(0, 5, 15, 85, 95, 100), include.lowest = TRUE, right = FALSE)
  trex.plot.redu$cuts <- factor(trex.plot.redu$cuts,
                               levels = c("[15,85)", "[5,15)", "[0,5)", "[85,95)", "[95,100]"))
  ordered.plot.redu <- trex.plot.redu[order(trex.plot.redu$cuts), ]
  
  # create T-REX plot
  trex_comparison <- paste(trex.condition[1], " vs ", trex.condition[2], sep = "") 
  trex.title <- paste( "<span style ='color:navyblue;'>", fcs.condition.label[trex.selection[1]], "</span> vs ",
         "<span style ='color:darkred;'>", fcs.condition.label[trex.selection[2]], "</span>", sep = "")
  
  trex.plot <-
    ggplot(ordered.plot.redu)+
    geom_point(aes(x = x, y = y, colour = cuts), 
               shape = 20, stroke = 0,
               size = redu.figure.point.size) +
    scale_color_manual(
      name = "Percent",
      values = c(
        "[15,85)" = trex.percent.15_85,
        "[5,15)" = trex.percent.5_15,
        "[0,5)" = trex.percent.0_5,
        "[85,95)" = trex.percent.85_95,
        "[95,100]" = trex.percent.95_100
      )
    ) +
    labs( title = trex.title )+
    theme_bw() + theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       axis.title = element_blank(), 
                       axis.text  = element_blank(), 
                       axis.ticks = element_blank(),
                       plot.title = element_markdown( size = trex.title.text.size )) +
    coord_fixed(ratio = graphical.ratio)+
    guides(color = guide_legend(override.aes = list(size=trex.legend.text.size)))
  
  ggsave( 
    file.path( redu.figure.dir, 
               paste0( redu.figure.plot, ".png"  ) ), 
    trex.plot, 
    width = trex.figure.width, height = trex.figure.height 
  )
}

