library(readr)
library(dplyr)
library(tidyr)
library(here)
df_filteredll <- read.csv(here::here("simulation", "Results", "betacsv", "betatest_overall_weaksetting.csv"))[, -1]

w_cols <- c("w=0", "w=0.05", "w=0.1", "w=0.15", "w=0.2", "w=0.25", "w=0.3", "w=0.35")
target_methods <- c("Proposed", "FARM")

df_wide <- df_filteredll %>%
  mutate(
    s = paste0("w=", s),
    method = case_when(
      tolower(method) == "proposed" ~ "Proposed",
      toupper(method) == "FARM" ~ "FARM",
      TRUE ~ as.character(method)
    ),
    method = factor(method, levels = target_methods),
    V = ifelse(grepl("^V", as.character(V)), as.character(V), paste0("V", V))
  ) %>%
  filter(s %in% w_cols, method %in% target_methods) %>%
  group_by(noise, B, V, method, s) %>%
  summarise(pval = mean(pval, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    id_cols = c(noise, B, V, method),
    names_from = s,
    values_from = pval
  ) %>%
  arrange(B, V, method)

print_standard_latex_table <- function(data) {
  cat("\\begin{table}[H]\n")
  cat("\\centering\n")
  cat("\\small\n")
  cat("\\renewcommand\\tabcolsep{2.5pt}\n")
  cat("\\caption{Empirical rejection rates (size and power) of the proposed and FARM tests under scenarios with weak factors.}\n")
  cat("\\label{tab:betaweak_1200}\n")
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
    v_list <- c("V1", "V2", "V3")

    for (v_name in v_list) {
      sub_v <- sub_b %>% filter(V == v_name)
      if (nrow(sub_v) == 0) next

      for (m_idx in seq_along(target_methods)) {
        m_name <- target_methods[m_idx]
        row_data <- sub_v %>% filter(method == m_name)
        if (nrow(row_data) == 0) next

        vals <- sapply(row_data[w_cols], function(x) sprintf("%.3f", as.numeric(x)))
        vals_str <- paste(vals, collapse = " & ")

        if (m_idx == 1) {
          cat(sprintf("\\multirow{2}{*}{%s} & %-8s & %s \\\\\n", v_name, m_name, vals_str))
        } else {
          cat(sprintf("                    & %-8s & %s \\\\\n", m_name, vals_str))
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

print_standard_latex_table(df_wide)
