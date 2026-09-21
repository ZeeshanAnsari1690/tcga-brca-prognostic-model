install.packages(c("timeROC", "survival", "readr"))
# load library
library(timeROC)
library(survival)
library(readr)


data <- read.csv("C:\\Cancer\\TCGA-BRCA\\results\\dataset_with_risk_scores.csv")

print(colnames(data))

head(data[, c("survival_time", "event", "Risk_Score")])

table(data$event)

summary(data$survival_time)
table(data$event)

# Create the variables

time <- data$survival_time
status <- data$event
marker <- data$Risk_Score

# Run the time-dependent ROC analysis

roc <- timeROC(
  T = time,
  delta = status,
  marker = marker,
  cause = 1,
  weighting = "marginal",
  times = c(365, 1095, 1825),
  iid = TRUE
)

# display the auc
print(roc$AUC)

# Save publication-quality ROC curve
# # Install packages (only once)
 
install.packages("ggplot2")
install.packages("tidyr")

# Load libraries
library(ggplot2)
library(tidyr)

# Create data frame
roc_df <- rbind(
  
  data.frame(
    FPR = roc$FP[,1],
    TPR = roc$TP[,1],
    Time = paste0("1-Year (AUC = ", round(roc$AUC[1],3),")")
  ),
  
  data.frame(
    FPR = roc$FP[,2],
    TPR = roc$TP[,2],
    Time = paste0("3-Year (AUC = ", round(roc$AUC[2],3),")")
  ),
  
  data.frame(
    FPR = roc$FP[,3],
    TPR = roc$TP[,3],
    Time = paste0("5-Year (AUC = ", round(roc$AUC[3],3),")")
  )
)
# Load package
library(ggplot2)
library(grid)

# Publication-quality ROC curve
p <- ggplot(
  roc_df,
  aes(x = FPR,
      y = TPR,
      colour = Time)
) +
  
  # ROC curves
  geom_line(
    linewidth = 1.8,
    lineend = "round"
  ) +
  
  # Reference diagonal
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    colour = "grey55",
    linewidth = 0.8
  ) +
  
  # Journal colours
  scale_color_manual(
    values = c(
      "#0072B2",
      "#D55E00",
      "#009E73"
    )
  ) +
  
  # Labels
  labs(
    title = "Time-dependent ROC Curve",
    x = "False Positive Rate",
    y = "True Positive Rate",
    colour = NULL
  ) +
  
  # Add these lines here
  scale_x_continuous(
    breaks = seq(0, 1, 0.25),
    limits = c(0, 1)
  ) +
  
  scale_y_continuous(
    breaks = seq(0, 1, 0.25),
    limits = c(0, 1)
  ) +
  
  coord_equal(
    xlim = c(0, 1),
    ylim = c(0, 1),
    expand = FALSE
  ) 
  
  
  theme(
    
    # Title
    plot.title = element_text(
      size = 16,
      face = "bold",
      hjust = 0.5,
      margin =margin(b=2)
    ),
    
    # Axis titles
    axis.title = element_text(
      size = 15,
      face = "bold",
      colour = "black"
    ),
    
    # Axis labels
    axis.text = element_text(
      size = 11,
      colour = "black"
    ),
    
    # Axis lines
    axis.line = element_line(
      colour = "black",
      linewidth = 0.8
    ),
    
    # Tick marks
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.7
    ),
    
    axis.ticks.length = unit(0.18, "cm"),
    
    # Border
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.8
    ),
    
    # Remove grid
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    # Legend
    legend.position = c(0.77, 0.16),
    
    legend.background = element_rect(
      fill = "white",
      colour = "black",
      linewidth = 0.5
    ),
    
    legend.key = element_blank(),
    
    legend.text = element_text(
      size = 10
    )
  )

# Display plot
print(p)

# Save figure
ggsave(
  filename = "C:/Cancer/TCGA-BRCA/results/Publication_ROC.png",
  plot = p,
  width = 7.5,
  height = 6.5,
  units = "in",
  dpi = 600,
  bg = "white"
)



