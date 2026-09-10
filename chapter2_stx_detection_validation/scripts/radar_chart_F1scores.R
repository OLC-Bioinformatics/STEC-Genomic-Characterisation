# radar plot of F1 score, sensitivity/recall and precision
# made using base R graphics - not ggplot2

#set working directory
setwd("/Users/noorshubair/Desktop/R_code_files/chapter2/")

# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(readxl)
#install.packages("fmsb")
library(fmsb)

radar <- read_excel("all_tools_subtype_data_forR.xlsx", sheet = "radar_data")

# Rename radar-chart axis labels for consistency with the thesis
colnames(radar)[colnames(radar) == "Recall"] <- "Sensitivity/Recall"
colnames(radar)[colnames(radar) == "Precision"] <- "Precision (PPV)"
radardf <- as.data.frame(radar, row.names = "Tool")
radardf
row.names(radardf) <- radardf$Tool
radardf<- radardf[,-1]


# Overlapping plot of tools -----------------------------------------------


# fmsb requires maximum (1) and minimum (0) rows as the first two rows
max_min <- data.frame(
  `Precision (PPV)` = c(1.0, 0.7),
  `Sensitivity/Recall` = c(1.0, 0.7),
  F1 = c(1.0, 0.7),
  check.names = FALSE
)
row.names(max_min) <- c("Max","Min")
max_min
# Combine the max/min with your tool data
plot_data <- rbind(max_min, radardf)
#plot_data <- rbind(max_min, radar)
plot_data
#rename it to F1 score
plot_data <- plot_data %>% rename(`F1 Score` = 'F1')

plot_data

combined_radar_plot <- function() {
  # Set graphical margins
  op <- par(
    mar = c(1.0, 2.6, 1.5, 2.6),
    xpd = NA,
    cex = 2.0
  ) #cex is increasing the graphics parameter size so our point shapes   are larger

  # Generate the radar chart with all together
  radarchart(
   plot_data,
   axistype = 1,       # 1 = center axis labels
   pty = 16, #point shape
   pcol = c("#66C2A5","#E76F51","#E69F00","#8DA0CB"), # Line colors for tools
   plwd = 5,           # Line width
   pfcol = scales::alpha(c("#66C2A5","#E76F51","#E69F00","#8DA0CB"), 0.05), # Transparent fills
    cglcol = "grey",    # Grid line color
    cglty = 1,          # Grid line type
    axislabcol = "black",# Axis label color
    seg = 3, #controls the number of segments
    caxislabels = c("0.7","0.8","0.9","1.0"), # Axis labels from 0 to 1
    calcex = 0.7,
    vlcex = 0.8         # Variable label size
  )

  # Add a legend for the tools
  legend(
    x = 0.5, y = 1.0,
   legend = c("KMA 70% ID, ConClave 1","Sipprverse","GeneSeekr","StxTyper"),#rownames  (radardf),
   bty = "n", pch = 16, col = c("#66C2A5","#E76F51","#E69F00","#8DA0CB"),
   cex = 0.6, pt.cex = 1
  )

  # Restore default graphical parameters
  par(op)

} 

dpi = 300
jpeg(
  filename = "F1_scores_radar_plot_combined.jpg", 
  width = 10*dpi, 
  height = 10*dpi, 
  res = dpi
  )

combined_radar_plot()
dev.off()


# Separate plots for each tool --------------------------------------------


#create separate plots for each tool
create_beautiful_radarchart <- function(data, color = "#00AFBB", 
                                        vlabels = colnames(data), vlcex = 0.9, #outer variables font size
                                        caxislabels = NULL, title = NULL, ...){
  radarchart(
    data, axistype = 1,
    # Customize the polygon
    pcol = color, pfcol = scales::alpha(color, 0.4), plwd = 2, plty = 1,
    # Customize the grid
    cglcol = "grey50", cglty = 1, cglwd = 0.8,
    # Customize the axis
    axislabcol = "grey30",
    # control number of segments
    seg = 3,
    # Variable labels
    vlcex = vlcex, vlabels = vlabels,
    calcex = 0.8,
    cex.main = 1.1, #title size
    caxislabels = caxislabels, title = title, ...
  )
}

# Define colors and titles
colors <- c("#66C2A5","#E76F51","#E69F00","#8DA0CB")
titles <- c("KMA 70% ID, ConClave 1","Sipprverse","GeneSeekr","StxTyper")


#wrapper function
separate_radar_plots <- function() {
  # Reduce plot margin using par()
  # Split the screen in 3 parts
  op <- par(
    mar = c(1.0, 2.8, 2.0, 2.8),
    mfrow = c(2, 2),
    xpd = NA
  )

  # Create the radar chart
  for(i in 1:4){
    create_beautiful_radarchart(
      data = plot_data[c(1, 2, i+2), ], #caxislabels = c("0.7","0.8","0.9","1.0"),
      caxislabels = c("70%","80%","90%","100%"),
      color = colors[i],
      title = titles[i],
      vlcex = 0.95
    )
  }

  par(op)
}

separate_radar_plots()

#lets save it
dpi = 300
jpeg("F1_scores_radar_plots_separate.jpg", width = 8*dpi, height = 8*dpi, res = dpi)
separate_radar_plots()
dev.off()





# if want to use ggplot ---------------------------------------------------

#ggplot method

#remotes::install_github("ricardo-bion/ggradar")
library(ggradar)
library(scales)

set.seed(4)
df <- data.frame(matrix(runif(30), ncol = 10))
df[, 1] <- paste0("G", 1:3)
colnames(df) <- c("Group", paste("Var", 1:9))
df

ggradar(df)

radar

ggradar(radar,
        background.circle.colour = "grey95",
        background.circle.transparency = 0.5,
        gridline.min.colour = "black",
        gridline.mid.colour = "black",
        gridline.max.colour = "black",
        point.alpha = 0.5, #transparency for all group points
        line.alpha = 0.8, #transparency for group lines
        group.colours = c("#66C2A5","#E76F51","#E69F00","#8DA0CB")
)















