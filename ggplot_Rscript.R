library(tidyverse)
setwd("~/")
#find datafile and make it a dataframe
figureA <- read.csv("~/BIOS-IN5410/BIOS-IN5410_H2021/data/data_file_1.csv")
figureA
figureA <- as_tibble(figureA)
figureA
pdf(file = "ggPlotFigA.pdf")
figureA %>% ggplot() + geom_point(aes(x = length, y = count))
dev.off() # And save the plot

