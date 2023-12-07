library(tidyverse)
setwd("~/")
args <- commandArgs(trailingOnly=TRUE) # This captures whatever you type on the command line
df <- read_csv("~/BIOS-IN5410/BIOS-IN5410_H2021/data/data_file_1.csv")
# Create an object that stores the file name
name <- args[1]
# Plot gene length vs. count
pdf(file = name) # This tells R that you want to make a pdf file, and what to call it
plot(df$length, df$count) # Make the plot
dev.off() # And save the plot
