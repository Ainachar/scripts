
# Load the data
data <- read.csv('~/ec34/biosin5410/microbiome/data/town_data.csv')
data
townA <- data$townA
townB <- data$townB
townB

# Make a summary overview of dog breeds and their counts in each town
countsA <- table(townA)
countsB <- table(townB)
countsB
countsA

# Sort data in a descending order
countsA <- sort(countsA, decreasing=TRUE)
countsB <- sort(countsB, decreasing=TRUE)

barplot(countsA, width = 1)
barplot(countsB)
dim(countsA)
class(countsA)
str(countsA)
str(countsB)
max(countsA)
typeof(countsA)

# Observed species
obs_sp <- function(counts) {
  obs_ind <- length(counts)
  return(obs_ind)
}
# Shannon diversity index
shannon <- function(counts) {
  total_count <- sum(counts)
  prop <- counts / total_count
  shan_ind <- -sum(prop * log(prop))
  return(shan_ind)
}
# Inversed Simpson index
simpson <- function(counts) {
  total_count <- sum(counts)
  prop <- counts / total_count
  simp_ind <- 1 / sum(prop^2)
  return(simp_ind)
}
obs_indA <- obs_sp(countsA)
obs_indA
shan_indA <- shannon(countsA)
shan_indA
simp_indA <- simpson(countsA)
obs_indB <- obs_sp(countsB)
shan_indB <- shannon(countsB)
simp_indB <- simpson(countsB)

obs_indA
obs_indB
shan_indA
shan_indB
simp_indA
simp_indB

# Calculate beta diversity between samples
jaccard <- function(townA, townB) {
  intersection <- length(intersect(townA, townB))
  combined_community <- unique(c(townA, townB))
  jac <- 1 - intersection / length(combined_community)
  return(jac)
}
bray_curtis <- function(countsA, countsB) {
  countsA <- data.frame(countsA)
  countsB <- data.frame(countsB)
  counts <- merge(countsA, countsB, by.x = 'townA', by.y='townB', all = TRUE)
  colnames(counts) <- c('Dog breed','townA','townB')
  counts$Diff <- abs(counts$townA - counts$townB)
  counts$Total <- counts$townA + counts$townB
  bray <- sum(counts$Diff) / sum(counts$Total)
  return(bray)
}

jac_ind <- jaccard(townA, townB)
bray_ind <- bray_curtis(countsA, countsB)
jac_ind
bray_ind

subsample_sets <- function(sample, step, lbl) {
  div_est <- data.frame(Depth = integer(),
                        Observed_species = numeric(),
                        Shannon = numeric(),
                        Simpson = numeric(),
                        stringsAsFactors = FALSE)
  depths <- seq(step, length(sample), by = step)
  for (dep in depths) {
    subsample <- sample(sample, dep)
    sub_counts <- table(subsample)
    div_est <- rbind(div_est, c(dep, obs_sp(sub_counts), shannon(sub_counts), simpson(sub_counts)))
  }
  colnames(div_est) <- c('Depth', 'Observed_species', 'Shannon', 'Simpson'
  )
  #Plot rarefaction curves
  par(mfrow=c(3, 1))
  plot(div_est$Depth, div_est$Observed_species, type='o', col='lightcoral'
       , xlab='Depth', ylab='Observed species')
  title(paste0('Observed species; ', lbl))
  plot(div_est$Depth, div_est$Shannon, type='o', col=' lightcoral', xlab='
Depth', ylab='Shannon index')
  title(paste0('Shannon index; ', lbl))
  plot(div_est$Depth, div_est$Simpson, type='o', col=' lightcoral', xlab='
Depth', ylab='Inversed Simpson index')
  title(paste0('Inversed Simpson index; ', lbl))
  return(div_est)
}

step <- 10 # Define minimal number of reads (steps)
div_estA <- subsample_sets(townA, step, 'townA')
div_estB <- subsample_sets(townB, step, 'townB')

#around 200 dogs

  
  
