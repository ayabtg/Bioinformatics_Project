# plot_profiles.R
# Minimal R script to plot energy profiles for base-pair potentials.
# Usage: Rscript plot_profiles.R input_profile.csv output_plot.png

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop('Usage: plot_profiles.R <input_csv> <output_png>')
}

library(ggplot2)

input <- args[1]
output <- args[2]

# Expected CSV columns: position, score
df <- try(read.csv(input), silent = TRUE)
if (inherits(df, 'try-error')) {
  stop('Cannot read input file')
}

p <- ggplot(df, aes(x=position, y=score)) + geom_line() + theme_minimal()

ggsave(output, plot = p)
