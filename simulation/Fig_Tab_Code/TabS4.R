library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(zoo)
library(here)
##----Table S4----

df_result <- read.csv(here::here("simulation", "Results", "betacsv","betatest_overall_lowernoise.csv"))[, -1]

w_cols <- c("w=0", "w=0.05", "w=0.1", "w=0.15", "w=0.2", "w=0.25", "w=0.3", "w=0.35")
target_methods <- c("Proposed", "FARM", "Joint")

df_wide <- df_result %>%
  mutate(
    s = ifelse(grepl("^w=", as.character(s)), as.character(s), paste0("w=", s)),
    method = case_when(
      tolower(method) == "proposed" ~ "Proposed",
      toupper(method) == "FARM"     ~ "FARM",
      tolower(method) == "joint"    ~ "Joint",
      TRUE                          ~ as.character(method)
    ),
    method = factor(method, levels = target_methods),
    V = ifelse(grepl("^V", as.character(V)), as.character(V), paste0("V", V))
  ) %>%
  filter(s %in% w_cols, method %in% target_methods) %>%
  group_by(V, method, s) %>%
  summarise(pval = mean(pval, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    id_cols = c(V, method),
    names_from = s,
    values_from = pval
  ) %>%
  arrange(V, method)

print_standard_table <- function(data) {
  cat("\\begin{table}[H]\n")
  cat("\\centering\n")
  cat("\\small\n")
  cat("\\renewcommand\\tabcolsep{2.5pt}\n")
  cat("\\caption{Empirical rejection rates (size and power) of the proposed and alternative tests under scenarios with a higher signal-to-noise ratio.}\n")
  cat("\\label{tab:beta_1200}\n")
  cat("\\begingroup\n")
  cat("\\renewcommand{\\arraystretch}{0.65}\n")
  cat("\\begin{tabular}{llcccccccc}\n")
  cat("\\toprule\n")
  cat("$\\boldsymbol{\\Sigma}_{\\boldsymbol{v}}$ & Method & $w=0$ & $w=0.05$ & $w=0.1$ & $w=0.15$ & $w=0.2$ & $w=0.25$ & $w=0.3$ & $w=0.35$ \\\\\n")
  cat("\\midrule\n")

  v_list <- c("V1", "V2", "V3")

  for (v_name in v_list) {
    sub_v <- data %>% filter(V == v_name)
    if (nrow(sub_v) == 0) next

    for (m_idx in seq_along(target_methods)) {
      m_name <- target_methods[m_idx]
      row_data <- sub_v %>% filter(method == m_name)
      if (nrow(row_data) == 0) next

      vals <- sapply(row_data[w_cols], function(x) sprintf("%.3f", as.numeric(x)))
      vals_str <- paste(vals, collapse = " & ")

      if (m_idx == 1) {
        cat(sprintf("\\multirow{3}{*}{%s} & %-8s & %s \\\\\n", v_name, m_name, vals_str))
      } else {
        cat(sprintf("                    & %-8s & %s \\\\\n", m_name, vals_str))
      }
    }
  }

  cat("\\bottomrule\n")
  cat("\\end{tabular}\n")
  cat("\\endgroup\n")
  cat("\\end{table}\n")
}

print_standard_table(df_wide)
