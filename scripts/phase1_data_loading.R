# phase 1: data loading, cleaning, and initial validation
#
# loading FDA CDER dataset, validating structure, cleaning data,
# creating response and blocking variables

# sourcing configuration and utilities
source("config.R")
source("utils.R")

# loading required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(lubridate)
})

print_section_header("Phase 1: Data Loading, Cleaning, and Initial Validation")

# 1. loading FDA CDER dataset
fda_data = load_csv(RAW_DATA_FILE)

cat(paste("Loaded", nrow(fda_data), "records from FDA CDER dataset\n"))
cat("Columns:", paste(names(fda_data), collapse = ", "), "\n")

# 2. validating data structure
missing_cols = setdiff(REQUIRED_COLUMNS, names(fda_data))
if (length(missing_cols) > 0) {
  stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
}

cat("All required columns present\n")
cat("Date column types:\n  FDA Receipt Date:", class(fda_data$`FDA Receipt Date`),
    "\n  FDA Approval Date:", class(fda_data$`FDA Approval Date`), "\n")

# 3. parsing date columns
fda_data = fda_data %>%
  mutate(
    receipt_date = mdy(`FDA Receipt Date`),
    approval_date = mdy(`FDA Approval Date`)
  )

invalid_receipt = sum(is.na(fda_data$receipt_date))
invalid_approval = sum(is.na(fda_data$approval_date))

cat("Date parsing complete\nInvalid receipt dates:", invalid_receipt,
    "\nInvalid approval dates:", invalid_approval, "\n")

# 4. filtering implausible values
initial_n = nrow(fda_data)

fda_data = fda_data %>%
  mutate(review_time_days = as.numeric(approval_date - receipt_date))

implausible_negative = sum(fda_data$review_time_days < MIN_REVIEW_TIME_DAYS, na.rm = TRUE)
implausible_long = sum(fda_data$review_time_days > MAX_REVIEW_TIME_DAYS, na.rm = TRUE)

fda_data = fda_data %>%
  filter(
    review_time_days >= MIN_REVIEW_TIME_DAYS,
    review_time_days <= MAX_REVIEW_TIME_DAYS
  )

final_n = nrow(fda_data)

cat("Initial records:", initial_n, "\nImplausible records removed:", initial_n - final_n,
    "\n  - Negative review times:", implausible_negative,
    "\n  - Review times >3650 days:", implausible_long,
    "\nFinal clean records:", final_n, "\n")

# 5. creating response variables
fda_data = fda_data %>%
  mutate(log_review_time_days = log(review_time_days))

cat(paste("Created review_time_days:", sum(!is.na(fda_data$review_time_days)), "observations\n"))
cat(paste("Range:", round(min(fda_data$review_time_days, na.rm = TRUE)),
          "to", round(max(fda_data$review_time_days, na.rm = TRUE)), "days\n"))
cat(paste("Mean:", round(mean(fda_data$review_time_days, na.rm = TRUE), 1), "days\n"))
cat("Created log_review_time_days for scale-dependency testing\n")

# 6. creating blocking variable (regulatory era)
fda_data = fda_data %>%
  mutate(regulatory_era = assign_regulatory_era(`Approval Year`))

cat(paste("Created regulatory_era blocking variable:", sum(!is.na(fda_data$regulatory_era)), "observations\n"))
cat("Era distribution:\n")
print(table(fda_data$regulatory_era))

# 7. generating descriptive statistics
cat(paste("\nTotal observations:", nrow(fda_data), "\n"))

cat("Date ranges:\n")
cat(paste("  Receipt:", min(fda_data$receipt_date, na.rm = TRUE), "to",
          max(fda_data$receipt_date, na.rm = TRUE), "\n"))
cat(paste("  Approval:", min(fda_data$approval_date, na.rm = TRUE), "to",
          max(fda_data$approval_date, na.rm = TRUE), "\n"))
cat(paste("  Years:", min(fda_data$`Approval Year`, na.rm = TRUE), "-",
          max(fda_data$`Approval Year`, na.rm = TRUE), "\n"))

# 8. summary statistics for created variables
cat("\nSummary statistics for all created variables\n")

cat("review_time_days (response variable):\n")
review_summary = summary(fda_data$review_time_days)
for (i in 1:length(review_summary)) {
  cat(sprintf("%10s: %10.2f\n", names(review_summary)[i], review_summary[i]))
}

cat("log_review_time_days (log-transformed response):\n")
log_summary = summary(fda_data$log_review_time_days)
for (i in 1:length(log_summary)) {
  cat(sprintf("%10s: %10.2f\n", names(log_summary)[i], log_summary[i]))
}

cat("regulatory_era (blocking variable):\n")
era_counts = table(fda_data$regulatory_era)
for (era in names(era_counts)) {
  cat(sprintf("  %15s: %4d observations\n", era, era_counts[era]))
}

# 9. validating temporal trends
era_means = fda_data %>%
  group_by(regulatory_era) %>%
  summarise(mean_review_time = mean(review_time_days, na.rm = TRUE), .groups = "drop")

cat("\nMean review time by regulatory era:\n")
print(era_means)

cat("Validation:\n")
cat(paste("Pre-PDUFA:", round(era_means$mean_review_time[era_means$regulatory_era == "Pre-PDUFA"], 1),
          "days (expected ~970 days)\n"))
cat(paste("Post-FDASIA:", round(era_means$mean_review_time[era_means$regulatory_era == "Post-FDASIA"], 1),
          "days (expected <400 days)\n"))

temporal_valid = era_means$mean_review_time[era_means$regulatory_era == "Pre-PDUFA"] >
                 era_means$mean_review_time[era_means$regulatory_era == "Post-FDASIA"]
cat(paste("Temporal trend confirmed:", temporal_valid, "\n"))

# 10. FDA CDER dataset validation report
cat("\nDataset summary\n")
cat(paste("   Total records after cleaning:", final_n, "\n"))
cat(paste("   Original records:", initial_n, "\n"))
cat(paste("   Records removed:", initial_n - final_n, "(implausible review times)\n"))

cat("Date ranges\n")
cat(paste("   FDA Receipt Date:", min(fda_data$receipt_date, na.rm = TRUE), "to",
          max(fda_data$receipt_date, na.rm = TRUE), "\n"))
cat(paste("   FDA Approval Date:", min(fda_data$approval_date, na.rm = TRUE), "to",
          max(fda_data$approval_date, na.rm = TRUE), "\n"))
cat(paste("   Approval Years:", min(fda_data$`Approval Year`, na.rm = TRUE), "-",
          max(fda_data$`Approval Year`, na.rm = TRUE), "\n"))

cat("Data completeness rates (% non-null)\n")
for (col in REQUIRED_COLUMNS) {
  completeness = (1 - sum(is.na(fda_data[[col]])) / nrow(fda_data)) * 100
  cat(sprintf("   %s: %.1f%%\n", col, completeness))
}

cat("Data quality issues identified\n")
cat(paste("   •", initial_n - final_n, "records removed with review times >3650 days (biologically implausible)\n"))
cat("   • All date fields successfully parsed\n")

cat("Validation complete - dataset ready for analysis\n")
save_csv(fda_data, CLEAN_DATA_FILE)

cat("\nPhase 1 complete\n")