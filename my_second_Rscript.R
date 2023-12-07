library(tidyverse)
setwd("~/")
df <- read_csv("~/BIOS-IN5410/BIOS-IN5410_H2021/data/data_file_1.csv")
# Plot gene length vs. count
pdf(file = "gene_length_vs_count.pdf") # This tells R that you want to make a pdf file, and what to call it
plot(df$length, df$count) # Make the plot
dev.off() # And save the plot
