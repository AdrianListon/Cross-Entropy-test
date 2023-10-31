cat("Performing statistical analysis on cluster and marker distributions")

# run statistical analysis on cluster distribution and marker expression

# run ANOVA on MFIs per sample, stats per group
flow.data.group.median$Group <- flow.sample.condition[flow.data.group.median$Sample]

flow.data.group.median.long <- pivot_longer(flow.data.group.median, !Group &!Sample &!PC1 &!PC2 &!sample.color,
                                    names_to = "Marker", values_to = "MFI")


for (fcs.marker in unique(flow.data.group.median.long$Marker)) {
  temp <- flow.data.group.median.long %>% 
    dplyr::filter(Marker == fcs.marker) %>%
    group_by(Group) %>%
    suppressMessages(summarise(MFI = MFI))
  
  res_aov <- aov( MFI ~ Group, data = temp )
  
  fitted.em <- emmeans(res_aov, "Group", data = temp )
  
  p.values <- data.frame( pairs(fitted.em, adjust = "tukey" ) )
  p.values <- p.values[,c(1,6)]
  p.values$significant <- ifelse( p.values$p.value < 0.05, "Yes", "No" )
  
  write.csv(p.values, file = paste0( fcs.mfi.stats.dir, fcs.marker, "_mfi_anova.csv"))
  
}

# run ANOVA on clusters, stats per group
flow.cluster.frame <- data.frame(flow.cluster.percent)
flow.cluster.frame$Sample <- rownames(flow.data.group.median)
flow.cluster.frame$Group <- flow.sample.condition[flow.cluster.frame$Sample]

flow.cluster.frame.long <- pivot_longer(flow.cluster.frame, !Group &!Sample,
                                            names_to = "Cluster", values_to = "Percent")

for (flow.cluster in unique(flow.cluster.frame.long$Cluster)) {
  temp <- flow.cluster.frame.long %>% 
    dplyr::filter(Cluster == flow.cluster) %>%
    group_by(Group) %>%
    suppressMessages(summarise(Percent = Percent))
  
  res_aov <- aov( Percent ~ Group, data = temp )
  
  fitted.em <- emmeans(res_aov, "Group", data = temp )
  
  p.values <- data.frame( pairs(fitted.em, adjust = "tukey" ) )
  p.values <- p.values[,c(1,6)]
  p.values$significant <- ifelse( p.values$p.value < 0.05, "Yes", "No" )
  
  write.csv(p.values, file = paste0( fcs.cluster.stats.dir, flow.cluster, "_frequencies_anova.csv"))
  
}


