library(readr)
library(dplyr)
library(tidyr)
library(here)
df <- read_csv(here::here("simulation", "Results", "betacsv", "betatest_overall_mainsetting.csv"),show_col_types = FALSE)[, -1]

for (pp in c(1200)) {
  w_levels <- c("0", "0.05", "0.1", "0.15", "0.2", "0.25", "0.3", "0.35")

  df_clean <- df %>%
    filter(
      p == pp,
      noise <= 2,
      V <= 3,
      s %in% w_levels
    ) %>%
    dplyr::select(-p, -F) %>%
    mutate(
      method = factor(ifelse(method == "proposed", "Proposed", "FARM"),
                      levels = c("Proposed", "FARM"))
    ) %>%
    arrange(noise, B, V, method) %>%
    pivot_wider(names_from = s, values_from = pval) %>%
    arrange(noise, B, V, method)

  print_exact_latex_table <- function(data) {
    w_cols <- c("0", "0.05", "0.1", "0.15", "0.2", "0.25", "0.3", "0.35")

    cat("  \\begin{table} \n")
    cat("\\centering\n")
    cat("\\small\n")
    cat("\\renewcommand\\tabcolsep{2.5pt}\n")
    cat("\\caption{Empirical rejection rates (size and power) of the proposed and FARM tests.}%, based on 2,000 replications at significance level $\\alpha=0.05$.}\n")
    cat(paste0("\\label{tab:beta_", pp, "}\n"))
    cat("\\begingroup\n")
    cat("\\renewcommand{\\arraystretch}{.65} \n")
    cat("\\begin{tabular}{clcccccccccc}  \n")
    cat("\\toprule\n")
    cat("Noise & $\\boldsymbol{\\Sigma}_{\\boldsymbol{v}}$ & Method & $w=0$ & $w=0.05$ & $w=0.1$ & $w=0.15$ & $w=0.2$ & $w=0.25$ & $w=0.3$ & $w=0.35$  \\\\\n")
    cat("\\midrule\n")

    noise_vals <- c(1, 2)
    noise_labels <- c("1" = "Normal", "2" = "Uniform")

    for (n_idx in seq_along(noise_vals)) {
      nv <- noise_vals[n_idx]
      n_label <- noise_labels[as.character(nv)]
      sub_noise <- data %>% filter(noise == nv)

      for (bv in c(1, 2)) {
        sub_b <- sub_noise %>% filter(B == bv)

        if (bv == 1) {
          cat("\\multicolumn{11}{c}{Gaussian-distributed Loading Matrix} \\\\  \n")
        } else {
          cat("\\multicolumn{11}{c}{Uniform-distributed Loading Matrix} \\\\")
        }

        for (vv in 1:3) {
          sub_v <- sub_b %>% filter(V == vv)
          v_label <- paste0("V", vv)

          row_prop <- sub_v %>% filter(method == "Proposed")
          vals_prop <- sapply(row_prop[w_cols], function(x) sprintf("%.3f", as.numeric(x)))

          row_farm <- sub_v %>% filter(method == "FARM")
          vals_farm <- sapply(row_farm[w_cols], function(x) sprintf("%.3f", as.numeric(x)))

           if (bv == 1 && vv == 1) {
            cat(sprintf("\\multirow{12}{*}{\\centering %s}  &  \\multirow{2}{*}{%s}  & Proposed & %s \\\\\n",
                        n_label, v_label, paste(vals_prop, collapse = " & ")))
          } else {
            cat(sprintf("   &  \\multirow{2}{*}{%s}  & Proposed & %s \\\\\n",
                        v_label, paste(vals_prop, collapse = " & ")))
          }

           cat(sprintf("                                    &                        & FARM     & %s \\\\",
                      paste(vals_farm, collapse = " & ")))
        }

        cat("\n")
        if (bv == 1) {
          cat("\\cmidrule(lr){2-11}\n")
        }
      }

      if (n_idx < length(noise_vals)) {
        cat("\\midrule\n")
      }
    }

    cat("\\bottomrule\n")
    cat("\\end{tabular}\n")
    cat("\\endgroup\n")
    cat("\\end{table}\n")
  }

  print_exact_latex_table(df_clean)
}
