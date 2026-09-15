## Write LaTeX tables for the corrected AR(1) review memo.

AR1_ROOT <- file.path("simulations", "power_ar1")
REVIEW_ROOT <- file.path("simulations", "z_AR1")
LEGACY_ROOT <- file.path(
  "simulations", "other_simus", "power_ar1_legacy_20260915"
)
source(file.path(AR1_ROOT, "run_power.R"))

calibration <- readRDS(file.path(AR1_ROOT, "ar1_calibration.rds"))
new <- readRDS(file.path(AR1_ROOT, "results.rds"))
old <- readRDS(file.path(LEGACY_ROOT, "results.rds"))

metrics <- c("F1", "Precision", "Recall", "CorrectNumCP", "MSE", "NumSegments")
new$method <- as.character(new$algorithm)
new$method[new$method == ar1_method_labels(calibration)[[4L]]] <-
  "SVP AR1Focus"
old$method <- as.character(old$algorithm)
old <- old[old$method %in% c(
  "PELT AR1 approximate", "PELT inflated", "DeCAFS AR1", "SVP AR1Focus"
), , drop = FALSE]

new_summary <- new |>
  dplyr::group_by(pattern, jump, method) |>
  dplyr::summarise(dplyr::across(dplyr::all_of(metrics), mean), .groups = "drop")
old_summary <- old |>
  dplyr::group_by(pattern, jump, method) |>
  dplyr::summarise(dplyr::across(dplyr::all_of(metrics), mean), .groups = "drop")
comparison <- dplyr::inner_join(
  new_summary, old_summary,
  by = c("pattern", "jump", "method"), suffix = c("_new", "_old")
)
for (metric in metrics) {
  comparison[[paste0("delta_", metric)]] <-
    comparison[[paste0(metric, "_new")]] -
    comparison[[paste0(metric, "_old")]]
}
utils::write.csv(
  comparison, file.path(REVIEW_ROOT, "paired_comparison.csv"), row.names = FALSE
)

selected <- calibration$validation[
  calibration$validation$method %in% c(
    "DeCAFS AR1", "PELT AR1 approximate", "SVP AR1Focus", "PELT inflated"
  ),
  c("method", "candidate", "null_f1", "null_f1_lower", "null_f1_upper")
]
selected$method[selected$method == "SVP AR1Focus"] <- "SVP AR1Focus / right"

write_tex <- function(data, path, caption, digits = 4) {
  lines <- capture.output(
    print(knitr::kable(
      data, format = "latex", booktabs = TRUE, longtable = FALSE,
      caption = caption, digits = digits, escape = TRUE
    )
  ))
  writeLines(lines, path)
}

write_tex(
  selected,
  file.path(REVIEW_ROOT, "ar1_calibration_table.tex"),
  "Validation null F1 for the selected AR(1) methods."
)

mean_delta <- comparison |>
  dplyr::group_by(pattern, method) |>
  dplyr::summarise(
    F1 = mean(delta_F1),
    Precision = mean(delta_Precision),
    Recall = mean(delta_Recall),
    Correct = mean(delta_CorrectNumCP),
    MSE = mean(delta_MSE),
    Segments = mean(delta_NumSegments),
    .groups = "drop"
  )
write_tex(
  mean_delta,
  file.path(REVIEW_ROOT, "ar1_before_after_table.tex"),
  "Mean refreshed-minus-legacy metric change by pattern and method."
)
