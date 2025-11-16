# plot_profiles.R
# Script R minimal pour tracer des profils (placeholder)
# Usage: Rscript plot_profiles.R input_profile.csv output_plot.png

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop('Usage: plot_profiles.R <input_csv> <output_png>')
}

library(ggplot2)

input <- args[1]
output <- args[2]

# Placeholder: on suppose que le CSV a des colonnes 'position' et 'score'
df <- try(read.csv(input), silent = TRUE)
if (inherits(df, 'try-error')) {
  stop('Impossible de lire le fichier input')
}

p <- ggplot(df, aes(x=position, y=score)) + geom_line() + theme_minimal()

ggsave(output, plot = p)
