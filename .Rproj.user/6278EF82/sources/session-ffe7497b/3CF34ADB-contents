library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(zoo)
library(stringr)
library(ggh4x)
library(mgcv)
library(here)
target_methods <- c("original", "common1", "spec1", "all3", "wrong_r")
method_table_labels <- c(
  "original" = "Proposed",
  "common1"  = "$\\boldsymbol{r} = (1,2,2,2,2)$",
  "spec1"    = "$\\boldsymbol{r} = (2,1,2,2,2)$",
  "all3"     = "$\\boldsymbol{r} = (3,3,3,3,3)$",
  "wrong_r"  = "$\\boldsymbol{r} = (1,1,1,1,1)$"
)
pp <- 1200
w_levels <- c("w=0", "w=0.05", "w=0.1", "w=0.15", "w=0.2", "w=0.25", "w=0.3", "w=0.35", "w=0.4")

w_cols <- c("w=0", "w=0.05", "w=0.1", "w=0.15", "w=0.2", "w=0.25", "w=0.3", "w=0.35")


print_latex_table_by_v <- function(df_data) {
  w_cols <- c("w=0", "w=0.05", "w=0.1", "w=0.15", "w=0.2", "w=0.25", "w=0.3", "w=0.35")

  cat("\\begin{table}[H]\n")
  cat("\\centering\n")
  cat("\\small\n")
  cat("\\renewcommand\\tabcolsep{3pt}\n")
  cat(paste0("\\caption{Empirical rejection rates (size and power) of the proposed test and its variants with pre-specified $\\boldsymbol{r}$  for $\\boldsymbol{\\beta}$ under Gaussian noise and loadings.}\n"))
  cat(paste0("\\label{tab:beta_proposed_", pp, "}\n"))
  cat("\\begingroup\n")
  cat("\\renewcommand{\\arraystretch}{0.75}\n")
  cat("\\begin{tabular}{lcccccccc}\n")
  cat("\\toprule\n")
  cat("Approach & $w=0$ & $w=0.05$ & $w=0.1$ & $w=0.15$ & $w=0.2$ & $w=0.25$ & $w=0.3$ & $w=0.35$ \\\\\n")
  cat("\\midrule\n")

  v_unique <- sort(unique(df_data$V))

  for (idx in seq_along(v_unique)) {
    cur_v <- v_unique[idx]

    cat(paste0("\\multicolumn{9}{c}{$\\boldsymbol{\\Sigma}_{\\boldsymbol{v}}$: V", cur_v, "} \\\\\n"))

    df_block <- df_data %>% filter(V == cur_v)

    for (i in 1:nrow(df_block)) {
      row <- df_block[i, ]
      meth_name <- as.character(row$group_method)

      label_str <- method_table_labels[meth_name]
      if (is.na(label_str)) label_str <- meth_name

      vals_str <- sapply(row[w_cols], function(val) {
        v <- as.numeric(val)
        if (is.na(v)) return("-") else sprintf("%.3f", v)
      })

      cat(paste0(
        label_str, " & ",
        paste(vals_str, collapse = " & "),
        " \\\\\n"
      ))
    }
    if (idx < length(v_unique)) {
      cat("\\midrule\n")
    }
  }

  cat("\\bottomrule\n")
  cat("\\end{tabular}\n")
  cat("\\endgroup\n")
  cat("\\end{table}\n")
}
df_table_sub <- read.csv(
  here::here("simulation", "Results", "betacsv", "betatest_overall_sensitivity.csv"),
  check.names = FALSE
)[, -1]
print_latex_table_by_v(df_table_sub)
