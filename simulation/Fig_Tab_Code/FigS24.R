library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(zoo)
library(here)
##--Figure S24----


plot_noise_allB_from_table <- function(noise_val, B_vals = NULL, v_vals = 1:3, method_vals = NULL, out_png = NULL, all_data_df) {
  if ("w" %in% names(all_data_df) && !"v" %in% names(all_data_df)) {
    # all_data_df <- all_data_df %>% rename(v = w)
    names(all_data_df)[names(all_data_df) == "w"] <- "v"

  }

  sub <- all_data_df %>% filter(noise == noise_val)
  if (!is.null(B_vals)) {
    sub <- sub %>% filter(B %in% B_vals)
  } else {
    B_vals <- sort(unique(sub$B))
  }

  if (!is.null(method_vals)) {
    sub <- sub %>% filter(method %in% method_vals)
  }

  df <- sub %>% filter(V %in% v_vals)
  if (nrow(df) == 0) return(NULL)

  all_v <- sort(unique(df$v))

  df_sum <- df %>%
    mutate(pval = as.numeric(pval), v = as.numeric(v)) %>%
    group_by(B, V, method, k, v) %>%
    summarise(pval = mean(pval, na.rm = TRUE), .groups = "drop") %>%
    group_by(B, V, method, k) %>%
    complete(v = all_v) %>%
    ungroup() %>%
    mutate(
      pval = ifelse(
        tolower(method) == "ifrm" & (is.na(pval) | is.nan(pval)),
        1,
        pval
      )
    ) %>%
    filter(!is.na(pval) & !is.nan(pval)) %>%
    arrange(B, V, method, k, v) %>%
    group_by(B, V, method, k) %>%
    mutate(pval = cummax(pval)) %>%
    ungroup()

  v_labeller <- function(value) paste0("V", value)

  method_labels <- c(
    "proposed" = "Proposed",
    "Proposed" = "Proposed",
    "jointadjust_lambin" = "joint_lambin",
    "jointadjust" = "Joint",
    "Joint" = "Joint",
    "joint" = "Joint",
    "IFRM" = "IFRM",
    "ifrm" = "IFRM"
  )

  modality_labels <- c(
    "1" = "1", "2" = "2", "3" = "3", "4" = "4",
    "kk1" = "1", "kk2" = "2", "kk3" = "3", "kk4" = "4",
    "common" = "common"
  )

  p <- ggplot(df_sum, aes(x = v, y = pval,
                          color = method,
                          shape = factor(k),
                          linetype = factor(k),
                          group = interaction(method, k))) +
    geom_line(linewidth = 0.8, na.rm = TRUE) +
    geom_point(size = 4, na.rm = TRUE) +
    facet_grid(B ~ V, labeller = labeller(V = as_labeller(v_labeller))) +
    scale_x_continuous(
      breaks = seq(
        min(df_sum$v),
        max(df_sum$v),
        by = 0.1
      )
    ) +
    labs(
      x = "Signal",
      y = "Empirical rejection rate",
      shape = "Modality",
      linetype = "Modality",
      color = "Method"
    ) +
    theme_bw() +
    theme(
      strip.text = element_blank(),
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 16),
      legend.key.size = unit(1.2, "lines"),
      legend.key.width = unit(3, "lines"),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 16, face = "bold"),
      legend.position = "bottom",
      legend.spacing.x = unit(0.3, "cm"),
      legend.spacing.y = unit(0.1, "cm"),
      panel.grid.minor = element_blank(),
      legend.box.background = element_blank()
    ) +
    scale_color_discrete(labels = function(x) {
      mapped <- method_labels[x]
      mapped[is.na(mapped)] <- x[is.na(mapped)]
      parse(text = mapped)
    }) +
    scale_shape_discrete(labels = modality_labels) +
    scale_linetype_discrete(labels = modality_labels) +
    guides(
      shape = guide_legend(override.aes = list(size = 4, linewidth = 0.8)),
      color = guide_legend(override.aes = list(size = 4, linewidth = 0.8))
    )

  if(!is.null(out_png)){
    ggsave(out_png, p, width = 16, height = 9, dpi = 300)
  }

  return(p)
}

df_clean <- read_csv(here::here("simulation", "Results", "betacsv","betatest_lowernoise.csv"), show_col_types = FALSE)[, -1]

if ("w" %in% names(df_clean)) {
  df_sub2 <- df_clean %>% filter(round(as.numeric(w), 2) != 0.4)
} else if ("v" %in% names(df_clean)) {
  df_sub2 <- df_clean %>% filter(round(as.numeric(v), 2) != 0.4)
} else {
  df_sub2 <- df_clean
}

selected_methods <- c("Proposed", "Joint", "IFRM")

for(ii in 1) {
  dd <- plot_noise_allB_from_table(
    noise_val = ii,
    B_vals = c(1, 2),
    method_vals = selected_methods,
    out_png =here::here("simulation", "Results", "betacsv","FigS24.png"),
    all_data_df = df_sub2
  )
  print(dd)
}
