## Generate only the conditional localization-error plots from saved results.

library(dplyr)
library(ggplot2)

source(file.path("simulations", "power_common.R"))

frameworks <- data.frame(
  folder = c("power_gaussian", "power_robust", "power_ar1"),
  n = c(1000L, 1000L, 600L),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(frameworks))) {
  root <- file.path("simulations", frameworks$folder[i])
  results <- readRDS(file.path(root, "results.rds"))
  results <- add_localization_error(results, frameworks$n[i])
  output <- file.path(root, "plots", "07_loc_error.pdf")
  ggplot2::ggsave(
    output, localization_error_plot(results),
    width = 12, height = 4.5
  )
}
