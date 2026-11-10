library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(zoo)
library(stringr)
library(ggh4x)
library(mgcv)
library(here)
all_data_fixed_gamma <- read.csv(
  here::here("simulation", "Results", "gammacsv", "gammatest_sensitivity.csv"),
  check.names = FALSE
)[, -1]
method_labels_text <- c(
  "original" = "Proposed",
  "common1"  = "bold(r) == '(1,2,2,2,2)'",
  "spec1"    = "bold(r) == '(2,1,2,2,2)'",
  "all3"     = "bold(r) == '(3,3,3,3,3)'",
  "wrong_r"  = "bold(r) == '(1,1,1,1,1)'"
)
target_methods <- c("original", "common1", "spec1", "all3", "wrong_r")
high_contrast_colors <- c(
  "original" = "#000000",
  "common1"  = "#E41A1C",
  "spec1"    = "#377EB8",
  "all3"     = "#984EA3",
  "wrong_r"  = "#4DAF4A"
)
parsed_method_labels <- setNames(parse(text = method_labels_text[target_methods]), target_methods)

# ==============================================================================
plot_noise_allB_group_method_gamma <- function(noise_val,
                                               B_vals = NULL,
                                               v_vals = 1:3,
                                               selected_groups = target_methods,
                                               out_png = NULL,
                                               all_data_df) {
  sub <- all_data_df %>%
    filter(noise == noise_val, group_method %in% selected_groups)

  if (!is.null(B_vals)) {
    sub <- sub %>% filter(B %in% B_vals)
  } else {
    B_vals <- sort(unique(sub$B))
  }

  df <- sub %>% filter(V %in% v_vals)
  target_v <- sort(unique(df$v))

  df_sum <- df %>%
    mutate(pval = as.numeric(pval)) %>%
    group_by(B, V, group_method, k, v) %>%
    summarise(pval = mean(pval, na.rm = TRUE), .groups = "drop") %>%
    filter(!is.na(pval)) %>%
    arrange(B, V, group_method, k, v) %>%
    group_by(B, V, group_method, k) %>%
    mutate(pval = cummax(pval)) %>%
    ungroup()

  v_labeller <- function(value) paste0("V", value)
  modality_labels <- c("1" = "1", "2" = "2", "3" = "3", "4" = "4", "common" = "common")

  p <- ggplot(df_sum, aes(x = v, y = pval,
                          color = group_method,
                          shape = factor(k),
                          linetype = factor(k),
                          group = interaction(group_method, k))) +
    geom_line(linewidth = 0.9, na.rm = TRUE) +
    geom_point(size = 3.6, na.rm = TRUE) +
    facet_grid(B ~ V, labeller = labeller(V = as_labeller(v_labeller))) +
    scale_x_continuous(breaks = target_v) +
    scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.2)) +
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

      # 🌟 核心：图例垂直堆叠（Method 一行，Modality 一行，共两行）
      legend.position = "bottom",
      legend.box = "vertical",
      legend.box.just = "center",

      # 适当紧凑键位尺寸，避免左右溢出画布
      legend.key.size = unit(1.0, "lines"),
      legend.key.width = unit(2.0, "lines"),
      legend.text = element_text(size = 13),
      legend.title = element_text(size = 15, face = "bold"),
      legend.spacing.x = unit(0.25, "cm"),
      legend.spacing.y = unit(0.15, "cm"),
      legend.margin = margin(t = 5, b = 5),

      panel.grid.minor = element_blank(),
      legend.box.background = element_blank()
    ) +
    scale_color_manual(
      breaks = target_methods,
      values = high_contrast_colors[target_methods],
      labels = parsed_method_labels[target_methods]
    ) +
    scale_shape_discrete(labels = modality_labels) +
    scale_linetype_discrete(labels = modality_labels) +

    # 🌟 强制将每个图例固定为 1 行并按顺序排列（Method 在上，Modality 在下）
    guides(
      color = guide_legend(
        order = 1,
        nrow = 1,
        byrow = TRUE,
        override.aes = list(size = 3.6, linewidth = 1.0)
      ),
      shape = guide_legend(
        order = 2,
        nrow = 1,
        byrow = TRUE,
        override.aes = list(size = 3.6, linewidth = 0.9)
      ),
      linetype = guide_legend(
        order = 2,
        nrow = 1,
        byrow = TRUE
      )
    )

  if (!is.null(out_png)) {
    ggsave(out_png, p, width = 16, height = 9, dpi = 300)
  }

  return(p)
}


full_v_vals <- c(0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5)


for (ii in 1) {
  p_gamma <- plot_noise_allB_group_method_gamma(
    noise_val        = ii,
    B_vals           = c(1),
    v_vals           = 1:3,
    selected_groups  = target_methods,
    out_png          = here::here("simulation", "Results", "gammacsv","FigS29.png"),
    all_data_df      = all_data_fixed_gamma
  )
  print(p_gamma)
}
