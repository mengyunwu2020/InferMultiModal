library(readr)
library(dplyr)
library(tidyr)
library(here)
pp <- 1200
w_cols <- c("w=0", "w=0.05", "w=0.1", "w=0.15", "w=0.2", "w=0.25", "w=0.3", "w=0.35")

tmp <- read_csv(here::here("simulation", "Results", "betacsv","betatest_overall_weaksetting.csv"), show_col_types = FALSE)[, -1]

target_methods <- c("Proposed", "Joint", "oracle_based", "oracle_r_based")

df_wide <- tmp %>%
  mutate(
    method = case_when(
      tolower(method) == "proposed" ~ "Proposed",
      tolower(method) == "joint" ~ "Joint",
      tolower(method) == "oracle_based" ~ "oracle_based",
      tolower(method) == "oracle_r_based" ~ "oracle_r_based",
      TRUE ~ method
    ),
    method = factor(method, levels = target_methods),
    V = ifelse(grepl("^V", as.character(V)), as.character(V), paste0("V", V)),
    s = ifelse(grepl("^w=", as.character(s)), as.character(s), paste0("w=", s))
  ) %>%
  filter(s %in% w_cols, method %in% target_methods) %>%
  group_by(B, V, method, s) %>%
  summarise(pval = mean(pval, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    id_cols = c(B, V, method),
    names_from = s,
    values_from = pval
  ) %>%
  arrange(B, V, method)

print_exact_latex_table <- function(data) {
  cat("\\begin{table}[H]\n")
  cat("\\centering\n")
  cat("\\small\n")
  cat("\\renewcommand\\tabcolsep{2.5pt}\n")
  cat("\\caption{Empirical rejection rates (size and power) of the proposed and alternative tests.}\n")
  cat(paste0("\\label{tab:beta_", pp, "}\n"))
  cat("\\begingroup\n")
  cat("\\renewcommand{\\arraystretch}{0.65}\n")
  cat("\\begin{tabular}{llcccccccc}\n")
  cat("\\toprule\n")
  cat("$\\boldsymbol{\\Sigma}_{\\boldsymbol{v}}$ & Method & $w=0$ & $w=0.05$ & $w=0.1$ & $w=0.15$ & $w=0.2$ & $w=0.25$ & $w=0.3$ & $w=0.35$ \\\\\n")
  cat("\\midrule\n")

  b_settings <- list(
    list(id = 1, title = "Gaussian-distributed Loading Matrix"),
    list(id = 2, title = "Uniform-distributed Loading Matrix")
  )

  for (b_idx in seq_along(b_settings)) {
    b_info <- b_settings[[b_idx]]
    cat(sprintf("\\multicolumn{10}{c}{%s} \\\\\n", b_info$title))

    sub_b <- data %>% filter(B == b_info$id)

    for (v_name in c("V1", "V2", "V3")) {
      sub_v <- sub_b %>% filter(V == v_name)
      if (nrow(sub_v) == 0) next

      for (m_idx in seq_along(target_methods)) {
        m_name <- target_methods[m_idx]
        row_data <- sub_v %>% filter(method == m_name)
        if (nrow(row_data) == 0) next

        vals <- sapply(row_data[w_cols], function(x) sprintf("%.3f", as.numeric(x)))
        vals_str <- paste(vals, collapse = " & ")
        m_label <- gsub("_", "\\\\_", m_name)

        if (m_idx == 1) {
          cat(sprintf("\\multirow{4}{*}{%s} & %-15s & %s \\\\\n", v_name, m_label, vals_str))
        } else {
          cat(sprintf("                    & %-15s & %s \\\\\n", m_label, vals_str))
        }
      }
    }

    if (b_idx < length(b_settings)) {
      cat("\\midrule\n")
    }
  }

  cat("\\bottomrule\n")
  cat("\\end{tabular}\n")
  cat("\\endgroup\n")
  cat("\\end{table}\n")
}

print_exact_latex_table(df_wide)
