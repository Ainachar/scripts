library(tidyverse)
args <- commandArgs(trailingOnly=TRUE) # This captures whatever you type on the command line
data <- args[1]

df <- read_csv(data)
# Create an object that stores the file name
name <- args[2]
# Plot gene length vs. count
pdf(file = name) # This tells R that you want to make a pdf file, and what to call it
plot(df$length, df$count) # Make the plot
dev.off() # And save the plot
