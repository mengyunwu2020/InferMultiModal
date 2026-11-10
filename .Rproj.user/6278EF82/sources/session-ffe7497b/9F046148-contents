library(readr)
library(dplyr)
library(tidyr)
library(here)
df_wide <- read.csv(here::here("simulation", "Results", "betacsv","betatest_onlycommon.csv"), check.names = FALSE)[, -1]

print_latex_onlycommon <- function(df_data) {
  w_cols <- c("0", "0.05", "0.1", "0.15", "0.2", "0.25", "0.3", "0.35")

  cat("\\begin{table}[H]\n")
  cat("\\centering\n")
  cat("\\small\n")
  cat("\\renewcommand\\tabcolsep{3pt}\n")
  cat("\\caption{Empirical rejection rates (size and power) of the proposed and FARM tests under the scenario where both the noise and the loading matrix are Gaussian-distributed with $r_k=0$ for $k\\in[K]$ based on 2,000 replications at significance level $\\alpha=0.05$.}\n")
  cat("\\label{tab:beta_1200_2}\n")
  cat("\\begingroup\n")
  cat("\\renewcommand{\\arraystretch}{1}\n")
  cat("\\begin{tabular}{lccccccccc}\n")
  cat("\\toprule\n")
  cat("$\\boldsymbol{\\Sigma}_{\\boldsymbol{v}}$ & Method & $w=0$ & $w=0.05$ & $w=0.1$ & $w=0.15$ & $w=0.2$ & $w=0.25$ & $w=0.3$ & $w=0.35$ \\\\\n")
  cat("\\midrule\n")

  v_unique <- unique(df_data$V)

  for (v_val in v_unique) {
    sub_v <- df_data %>% filter(V == v_val)
    n_methods <- nrow(sub_v)

    for (i in 1:n_methods) {
      row <- sub_v[i, ]

      cov_text <- ifelse(i == 1, paste0("\\multirow{", n_methods, "}{*}{", v_val, "}"), "")

      vals <- sapply(row[w_cols], function(x) sprintf("%.3f", as.numeric(x)))
      vals_str <- paste(vals, collapse = "  &  ")

      cat(paste0(
        sprintf("%-19s", cov_text), "  &  ",
        sprintf("%-8s", as.character(row$method)), "  &  ",
        vals_str, "  \\\\\n"
      ))
    }
  }

  cat("\\bottomrule\n")
  cat("\\end{tabular}\n")
  cat("\\endgroup\n")
  cat("\\end{table}\n")
}

print_latex_onlycommon(df_wide)
