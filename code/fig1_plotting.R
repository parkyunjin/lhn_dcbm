### Plot Figure 1 ###

library(gridExtra)
library(ggplot2)
library(showtext)
library(grid)
library(igraph)


## BALANCED SBM
# Load graph data
balanced_sbm_g <- read_graph("/result/balanced_sbm_curved.graphml", format = "graphml")
# Make edge dataframe
edge_list <- as.data.frame(get.edgelist(balanced_sbm_g))
results1 <- cbind(edge_list, bfc = E(balanced_sbm_g)$bfc,  frc = E(balanced_sbm_g)$frc,  jac_coef = E(balanced_sbm_g)$jac_coef, lrc = E(balanced_sbm_g)$lrc, nrc = E(balanced_sbm_g)$nrc, same_group = E(balanced_sbm_g)$membership)
colnames(results1) <- c("node_i", "node_j","bfc", "frc", "jac_coef", "lrc", "nrc", "same_group")

## UNBALANCED DCBM
# Load graph data
dcsbm_g <- read_graph("/result/dcsbm_curved.graphml", format = "graphml")
dcsbm_g <- cal_linear_curv(dcsbm_g)
# Make edge dataframe
edge_list <- as.data.frame(get.edgelist(dcsbm_g))
results3 <- cbind(edge_list,  bfc = E(dcsbm_g)$bfc, frc = E(dcsbm_g)$frc, jac_coef = E(dcsbm_g)$jac_coef, lrc = E(dcsbm_g)$lrc, nrc = E(dcsbm_g)$nrc, same_group = E(dcsbm_g)$membership)
colnames(results3) <- c("node_i", "node_j","bfc", "frc", "jac_coef", "lrc", "nrc", "same_group")


## PLOTTING

# Add subplot labels
label_theme <- theme(
  plot.title = element_text(
    hjust = -0.1, vjust = 1.5, size = 12
  )
)

# Plot for BALANCED SBM
 p1 <- ggplot(results1, aes(x = bfc, fill = same_group)) +
   geom_histogram(bins = 30, color = "black", alpha = 0.5, position = "identity") +
   labs(
     x = "BFC value",
     y = "Frequency") +
   theme_minimal()+
   scale_fill_manual(values = c("Same" = "blue", "Different" = "red"), guide = "none") +
   theme(
     axis.title = element_text(size = 32),    
     axis.text = element_text(size = 8.5)
   )


p2 <- ggplot(results1, aes(x = frc, fill = same_group)) +
  geom_histogram(bins = 30,  color = "black", alpha = 0.5, position = "identity") +
  labs(
    x = "FRC value",
    y = "Frequency") +
  theme_minimal()+
  scale_fill_manual(values = c("Same" = "blue", "Different" = "red"), guide = "none")+
  theme(
    axis.title = element_text(size = 32),    
    axis.text = element_text(size = 8.5)
  )


p25 <- ggplot(results1, aes(x = jac_coef, fill = same_group)) +
  geom_histogram(bins = 30, color = "black", alpha = 0.5, position = "identity") +
  labs(
    x = "Jaccard value",
    y = "Frequency") +
  theme_minimal()+
  scale_fill_manual(values = c("Same" = "blue", "Different" = "red"), guide = "none")+
  theme(
    axis.title = element_text(size = 32),    
    axis.text = element_text(size = 8.5)
  )


p3 <- ggplot(results1, aes(x = lrc, fill = same_group)) +
  geom_histogram(bins = 30, color = "black", alpha = 0.5, position = "identity") +
  labs(
    x = "LRC value",
    y = "Frequency") +
  theme_minimal()+
  scale_fill_manual(values = c("Same" = "blue", "Different" = "red"), guide = "none")+
  theme(
    axis.title = element_text(size = 32),    
    axis.text = element_text(size = 8.5)
  )


p4 <- ggplot(results1, aes(x = nrc, fill = same_group)) +
  geom_histogram(bins = 30, color = "black", alpha = 0.5, position = "identity") +
  labs(
    x = "DCRC value",
    y = "Frequency") +
  theme_minimal()+
  scale_fill_manual(values = c("Same" = "blue", "Different" = "red"), guide = "none")+
  theme(
    axis.title = element_text(size = 32),    
    axis.text = element_text(size = 8.5)
  )
  





# Plot for DCSBM
 p5 <- ggplot(results3, aes(x = bfc, fill = same_group)) +
   geom_histogram(bins = 30, color = "black", alpha = 0.5, position = "identity") +
   labs(
     x = "BFC value",
     y = "Frequency") +
   theme_minimal()+
   scale_fill_manual(values = c("Same" = "blue", "Different" = "red"), guide = "none")+
   theme(
     axis.title = element_text(size = 32),    
     axis.text = element_text(size = 8.5)
   )
   


p6 <- ggplot(results3, aes(x = frc, fill = same_group)) +
  geom_histogram(bins = 30,  color = "black", alpha = 0.5, position = "identity") +
  labs(
    x = "FRC value",
    y = "Frequency") +
  theme_minimal()+
  scale_fill_manual(values = c("Same" = "blue", "Different" = "red"), guide = "none")+
  theme(
    axis.title = element_text(size = 32),    
    axis.text = element_text(size = 8.5)
  )

p65 <- ggplot(results3, aes(x = jac_coef, fill = same_group)) +
  geom_histogram(bins = 30, color = "black", alpha = 0.5, position = "identity") +
  labs(
    x = "Jaccard value",
    y = "Frequency") +
  theme_minimal()+
  scale_fill_manual(values = c("Same" = "blue", "Different" = "red"), guide = "none")+
  theme(
    axis.title = element_text(size = 32),    
    axis.text = element_text(size = 8.5)
  )


p7 <- ggplot(results3, aes(x = lrc, fill = same_group)) +
  geom_histogram(bins = 30, color = "black", alpha = 0.5, position = "identity") +
  labs(
    x = "LRC value",
    y = "Frequency") +
  theme_minimal()+
  scale_fill_manual(values = c("Same" = "blue", "Different" = "red"), guide = "none")+
  theme(
    axis.title = element_text(size = 32),    
    axis.text = element_text(size = 8.5)
  )
  


p8 <- ggplot(results3, aes(x = nrc, fill = same_group)) +
  geom_histogram(bins = 30, color = "black", alpha = 0.5, position = "identity") +
  labs(
    x = "DCRC value",
    y = "Frequency") +
  theme_minimal()+
  scale_fill_manual(values = c("Same" = "blue", "Different" = "red"), guide = "none")+
  theme(
    axis.title = element_text(size = 32),   
    axis.text = element_text(size = 8.5)
  )



# Titles
title0 <- textGrob("SBM", gp = gpar(fontsize = 36, fontface = "bold"))
title1 <- textGrob("DCBM", gp = gpar(fontsize = 36, fontface = "bold"))

# Save figures in pdf
pdf("figure_1.pdf", width = 32, height = 4)
remove_y_axis_title <- function(p) {
  p + theme(axis.title.y = element_blank())
}
grid.arrange(
  arrangeGrob(  
    title0,
    arrangeGrob(p1, p2, p25, p3, p4, ncol = 5),
    nrow = 2,
    heights = c(0.1, 1)
  ),
  arrangeGrob(  
    title1,
    arrangeGrob(p5, p6, p65, p7, p8, ncol = 5),
    nrow = 2,
    heights = c(0.1, 1)
  ),
  ncol = 2
  
)
dev.off()
