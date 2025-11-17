#!/usr/bin/env Rscript
library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)
input <- args[1]
output <- args[2]

df <- read.csv(input)

p <- ggplot(df, aes(x = position, y = score)) +
  geom_line(color="blue") +
  geom_point(size=1) +
  theme_minimal() +
  ggtitle("RNA Statistical Potential Profile") +
  xlab("Nucleotide position") +
  ylab("Energy score")

ggsave(output, p, width=7, height=4)
