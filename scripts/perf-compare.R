#!/usr/bin/env Rscript

get_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  key <- paste0("--", name, "=")
  hit <- grep(paste0("^", key), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(key, "", hit[1], fixed = TRUE)
}

base_path <- get_arg("base", NULL)
new_path <- get_arg("new", NULL)
metric <- get_arg("metric", "elapsed_mean")
fail_pct <- as.numeric(get_arg("fail_pct", "NA"))

if (is.null(base_path) || is.null(new_path)) {
  stop("Usage: Rscript scripts/perf-compare.R --base=BASE.csv --new=NEW.csv [--metric=elapsed_mean] [--fail_pct=10]")
}

base <- utils::read.csv(base_path, stringsAsFactors = FALSE)
new <- utils::read.csv(new_path, stringsAsFactors = FALSE)

required_cols <- c("case", metric)
if (!all(required_cols %in% names(base))) {
  stop(sprintf("Base CSV missing required columns: %s", paste(required_cols, collapse = ", ")))
}
if (!all(required_cols %in% names(new))) {
  stop(sprintf("New CSV missing required columns: %s", paste(required_cols, collapse = ", ")))
}

cmp <- merge(
  base[, required_cols],
  new[, required_cols],
  by = "case",
  suffixes = c("_base", "_new")
)

if (nrow(cmp) == 0) {
  stop("No overlapping case labels between baseline and new CSV")
}

base_col <- paste0(metric, "_base")
new_col <- paste0(metric, "_new")
cmp$delta <- cmp[[new_col]] - cmp[[base_col]]
cmp$pct_change <- 100 * cmp$delta / cmp[[base_col]]

cat("\nPerformance comparison\n")
cat(sprintf("metric=%s\n\n", metric))
print(cmp, row.names = FALSE, digits = 4)

if (!is.na(fail_pct)) {
  bad <- cmp$pct_change > fail_pct
  if (any(bad)) {
    cat(sprintf(
      "\nRegression detected: %d case(s) slower than %.2f%% threshold\n",
      sum(bad), fail_pct
    ))
    quit(status = 2)
  }
  cat(sprintf("\nNo regression above %.2f%% threshold\n", fail_pct))
}
