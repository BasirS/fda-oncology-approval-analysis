# phase 4b: expedited programs analysis
#
# analyzing usage of FDA expedited programs (Accelerated Approval, Fast Track,
# Orphan Drug) to explain differential review times between therapeutic areas

# sourcing configuration and utilities
source("config.R")
source("utils.R")

# loading required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

print_section_header("Phase 4b: Expedited Programs Analysis")

# 1. loading classified data
input_file = file.path(RESULTS_DIR, "fda_analysis_clean_classified.csv")
classified_data = load_csv(input_file)

# filtering to analysis-ready sample (exclude Uncertain)
analysis_data = classified_data %>%
  filter(therapeutic_area != "Uncertain")

cat(paste("Analysis-ready sample (excluding Uncertain): n =", nrow(analysis_data), "\n"))

# 2. analyzing Accelerated Approval usage
acc_approval_table = table(
  analysis_data$therapeutic_area,
  analysis_data$`Accelerated Approval`
)

cat("\nAccelerated Approval by Therapeutic Area:\n")
print(addmargins(acc_approval_table))

acc_approval_pct = prop.table(acc_approval_table, margin = 1) * 100
cat("Percentages (row-wise):\n")
print(round(acc_approval_pct, 2))

# extracting key statistics
acc_oncology_yes = acc_approval_table["Oncology", "Yes"]
acc_oncology_total = sum(acc_approval_table["Oncology", ])
acc_oncology_pct_val = (acc_oncology_yes / acc_oncology_total) * 100

acc_other_yes = acc_approval_table["Other", "Yes"]
acc_other_total = sum(acc_approval_table["Other", ])
acc_other_pct_val = (acc_other_yes / acc_other_total) * 100

cat(sprintf("\nOncology: %d/%d (%.1f%%) received Accelerated Approval\n",
            acc_oncology_yes, acc_oncology_total, acc_oncology_pct_val))
cat(sprintf("Other: %d/%d (%.1f%%) received Accelerated Approval\n",
            acc_other_yes, acc_other_total, acc_other_pct_val))
cat(sprintf("Ratio: Oncology is %.1fx more likely to receive Accelerated Approval\n",
            acc_oncology_pct_val / acc_other_pct_val))

# 3. analyzing Fast Track usage
fast_track_table = table(
  analysis_data$therapeutic_area,
  analysis_data$`Fast Track Designation`
)

cat("\nFast Track by Therapeutic Area:\n")
print(addmargins(fast_track_table))

fast_track_pct = prop.table(fast_track_table, margin = 1) * 100
cat("Percentages (row-wise):\n")
print(round(fast_track_pct, 2))

# extracting key statistics
ft_oncology_yes = fast_track_table["Oncology", "Yes"]
ft_oncology_total = sum(fast_track_table["Oncology", ])
ft_oncology_pct_val = (ft_oncology_yes / ft_oncology_total) * 100

ft_other_yes = fast_track_table["Other", "Yes"]
ft_other_total = sum(fast_track_table["Other", ])
ft_other_pct_val = (ft_other_yes / ft_other_total) * 100

cat(sprintf("\nOncology: %d/%d (%.1f%%) received Fast Track\n",
            ft_oncology_yes, ft_oncology_total, ft_oncology_pct_val))
cat(sprintf("Other: %d/%d (%.1f%%) received Fast Track\n",
            ft_other_yes, ft_other_total, ft_other_pct_val))
cat(sprintf("Ratio: Oncology is %.1fx more likely to receive Fast Track\n",
            ft_oncology_pct_val / ft_other_pct_val))

# 4. analyzing Orphan Drug Designation usage
analysis_data = analysis_data %>%
  mutate(
    orphan_binary = if_else(
      `Orphan Drug Designation` %in% c("Yes", "yes"),
      "Yes",
      "No"
    )
  )

orphan_table = table(
  analysis_data$therapeutic_area,
  analysis_data$orphan_binary
)

cat("\nOrphan Drug by Therapeutic Area:\n")
print(addmargins(orphan_table))

orphan_pct = prop.table(orphan_table, margin = 1) * 100
cat("Percentages (row-wise):\n")
print(round(orphan_pct, 2))

# extracting key statistics
orp_oncology_yes = orphan_table["Oncology", "Yes"]
orp_oncology_total = sum(orphan_table["Oncology", ])
orp_oncology_pct_val = (orp_oncology_yes / orp_oncology_total) * 100

orp_other_yes = orphan_table["Other", "Yes"]
orp_other_total = sum(orphan_table["Other", ])
orp_other_pct_val = (orp_other_yes / orp_other_total) * 100

cat(sprintf("\nOncology: %d/%d (%.1f%%) received Orphan Designation\n",
            orp_oncology_yes, orp_oncology_total, orp_oncology_pct_val))
cat(sprintf("Other: %d/%d (%.1f%%) received Orphan Designation\n",
            orp_other_yes, orp_other_total, orp_other_pct_val))
cat(sprintf("Ratio: Oncology is %.1fx more likely to receive Orphan Designation\n",
            orp_oncology_pct_val / orp_other_pct_val))

# 5. creating summary table
summary_table = data.frame(
  program = c("Accelerated Approval", "Fast Track", "Orphan Drug"),
  oncology_pct = c(acc_oncology_pct_val, ft_oncology_pct_val, orp_oncology_pct_val),
  other_pct = c(acc_other_pct_val, ft_other_pct_val, orp_other_pct_val)
) %>%
  mutate(
    ratio = oncology_pct / other_pct,
    diff = oncology_pct - other_pct
  )

cat("\nSummary of expedited program usage:\n")
print(summary_table)

max_ratio_program = summary_table$program[which.max(summary_table$ratio)]
max_ratio_value = max(summary_table$ratio)

cat(sprintf("\nLargest disparity: %s (%.1fx more likely for Oncology)\n",
            max_ratio_program, max_ratio_value))

# 6. saving results
output_file = file.path(RESULTS_DIR, "expedited_programs_summary.csv")
save_csv(summary_table, output_file)

cat("\nPhase 4b complete\n")
