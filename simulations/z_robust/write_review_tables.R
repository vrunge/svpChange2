## Create compact LaTeX tables used by review_robust.Rmd and the replacement memo.

library(dplyr)

source(file.path("simulations", "power_robust", "run_power.R"))

ROOT <- normalizePath(".", mustWork = TRUE)
OUT <- file.path(ROOT, "simulations", "z_robust")
ROBUST <- file.path(ROOT, "simulations", "power_robust")
LEGACY <- file.path(ROOT, "simulations", "other_simus",
                    "power_robust_legacy_20260915")

calibration <- readRDS(file.path(ROBUST, "robust_calibration.rds"))
new <- readRDS(file.path(ROBUST, "results.rds"))
old <- readRDS(file.path(LEGACY, "results.rds"))
old <- old[old$algorithm %in% c(
  "PELT", "RFPOP paper", "SVP MedianMood", "SVP Wilcoxon"
), , drop = FALSE]
new$method <- as.character(new$algorithm)
new$method[new$method == robust_method_labels(calibration)[[3L]]] <-
  "MedianMood"
new$method[new$method == robust_method_labels(calibration)[[4L]]] <-
  "Wilcoxon"
old$method <- c(PELT = "PELT", `RFPOP paper` = "RFPOP",
                `SVP MedianMood` = "MedianMood", `SVP Wilcoxon` = "Wilcoxon")[
                  as.character(old$algorithm)]
metrics <- c("F1", "Precision", "Recall", "CorrectNumCP", "MSE", "NumSegments")
summarise_method <- function(x, prefix) {
  x |>
    group_by(pattern, jump, method) |>
    summarise(across(all_of(metrics), mean), .groups = "drop") |>
    rename_with(~ paste0(prefix, .x), all_of(metrics))
}
comparison <- full_join(summarise_method(new, "new_"),
                        summarise_method(old, "old_"),
                        by = c("pattern", "jump", "method")) |>
  mutate(across(starts_with("new_"), as.numeric),
         across(starts_with("old_"), as.numeric))
for (metric in metrics) {
  comparison[[paste0("delta_", metric)]] <-
    comparison[[paste0("new_", metric)]] -
    comparison[[paste0("old_", metric)]]
}
write.csv(comparison, file.path(OUT, "paired_comparison.csv"), row.names = FALSE)

compact <- comparison |>
  group_by(pattern, method) |>
  summarise(across(starts_with("delta_"), mean), .groups = "drop")

escape_latex <- function(x) gsub("_", "\\\\_", x, fixed = TRUE)
write_tex_table <- function(x, path, caption, columns) {
  con <- file(path, "w")
  on.exit(close(con), add = TRUE)
  writeLines(c(
    "\\begin{table}[ht]", "\\centering", paste0("\\caption{", caption, "}"),
    paste0("\\begin{tabular}{", paste(rep("l", length(columns)), collapse = ""), "}"),
    "\\toprule", paste(columns, collapse = " & "), " \\\\", "\\midrule"
  ), con)
  for (i in seq_len(nrow(x))) {
    values <- vapply(x[i, , drop = FALSE], function(v) {
      if (is.numeric(v)) sprintf("%.4f", v) else escape_latex(as.character(v))
    }, character(1))
    writeLines(paste0(paste(values, collapse = " & "), " \\\\"), con)
  }
  writeLines(c("\\bottomrule", "\\end{tabular}", "\\end{table}"), con)
}

selected <- calibration$table |>
  filter(stage == "validation", selected) |>
  select(method, candidate, null_f1, null_f1_lower, null_f1_upper,
         rfpop_target_f1, seed)
write_tex_table(
  selected, file.path(OUT, "robust_calibration_table.tex"),
  "Validation of selected robust null thresholds.",
  c("Method", "Candidate", "F1", "Lower", "Upper", "RFPOP target", "Seed")
)
write_tex_table(
  compact, file.path(OUT, "before_after_table.tex"),
  "Mean new-minus-legacy metric change by pattern and method.",
  c("Pattern", "Method", "F1", "Precision", "Recall", "Correct", "MSE", "Segments")
)
