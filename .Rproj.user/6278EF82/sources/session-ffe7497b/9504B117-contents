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

method_labels_text <- c(
  "original" = "Proposed",
  "common1"  = "bold(r) == '(1,2,2,2,2)'",
  "spec1"    = "bold(r) == '(2,1,2,2,2)'",
  "all3"     = "bold(r) == '(3,3,3,3,3)'",
  "wrong_r"  = "bold(r) == '(1,1,1,1,1)'"
)

parsed_method_labels <- setNames(parse(text = method_labels_text[target_methods]), target_methods)

method_table_labels <- c(
  "original" = "Proposed",
  "common1"  = "$\\boldsymbol{r} = (1,2,2,2,2)$",
  "spec1"    = "$\\boldsymbol{r} = (2,1,2,2,2)$",
  "all3"     = "$\\boldsymbol{r} = (3,3,3,3,3)$",
  "wrong_r"  = "$\\boldsymbol{r} = (1,1,1,1,1)$"
)

high_contrast_colors <- c(
  "original" = "#000000",
  "common1"  = "#E41A1C",
  "spec1"    = "#377EB8",
  "all3"     = "#984EA3", #
  "wrong_r"  = "#4DAF4A"
)

my_shapes <- setNames(c(15, 16, 17, 18, 8), target_methods)

v_labels    <- c("1" = "V1", "2" = "V2", "3" = "V3")
label_value <- c("1" = "Gaussian noise", "2" = "Uniform noise")
target_grid <- seq(from = 0.25, to = 0.60, by = 0.0125)

plot_noise_allB_group_method <- function(noise_val,
                                         B_vals = NULL,
                                         v_vals = 1:3,
                                         group_method_vals = target_methods,
                                         out_png = NULL,
                                         all_data_df) {

  sub <- all_data_df %>%
    filter(noise == noise_val, group_method %in% group_method_vals)

  if (!is.null(B_vals)) {
    sub <- sub %>% filter(B %in% B_vals)
  } else {
    B_vals <- sort(unique(sub$B))
  }

  df <- sub %>% filter(V %in% v_vals)
  all_v <- sort(unique(df$v))

  df_sum <- df %>%
    mutate(pval = as.numeric(pval)) %>%
    group_by(B, V, group_method, k, v) %>%
    summarise(pval = mean(pval, na.rm = TRUE), .groups = "drop") %>%
    group_by(B, V, group_method, k) %>%
    complete(v = all_v) %>%
    ungroup() %>%
    filter(!is.na(pval) & !is.nan(pval)) %>%
    arrange(B, V, group_method, k, v) %>%
    group_by(B, V, group_method, k) %>%
    mutate(pval = cummax(pval)) %>%
    ungroup()

  v_labeller <- function(value) paste0("V", value)

  modality_labels <- c(
    "kk1" = "1", "kk2" = "2", "kk3" = "3", "kk4" = "4", "common" = "common",
    "1" = "1", "2" = "2", "3" = "3", "4" = "4"
  )

  p <- ggplot(df_sum, aes(x = v, y = pval,
                          color = group_method,
                          shape = factor(k),
                          linetype = factor(k),
                          group = interaction(group_method, k))) +
    geom_line(linewidth = 0.9, na.rm = TRUE) +
    geom_point(size = 3.6, na.rm = TRUE) +
    facet_grid(B ~ V, labeller = labeller(V = as_labeller(v_labeller))) +
    scale_x_continuous(breaks = seq(min(df_sum$v), max(df_sum$v), by = 0.05)) +
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
      axis.text = element_text(size = 15),
      axis.title = element_text(size = 16),
      strip.text = element_text(size = 14, face = "bold"),
      legend.key.size = unit(1.2, "lines"),
      legend.key.width = unit(2.5, "lines"),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 15, face = "bold"),
      legend.position = "bottom",
      legend.spacing.x = unit(0.4, "cm"),
      panel.grid.minor = element_blank()
    ) +
    scale_color_manual(
      breaks = target_methods,
      values = high_contrast_colors[target_methods],
      labels = parsed_method_labels[target_methods]
    ) +
    scale_shape_discrete(labels = modality_labels) +
    scale_linetype_discrete(labels = modality_labels) +
    guides(
      shape = guide_legend(override.aes = list(size = 3.6, linewidth = 0.9)),
      color = guide_legend(override.aes = list(size = 3.6, linewidth = 1.0))
    )

  if (!is.null(out_png)) {
    ggsave(out_png, p, width = 16, height = 7.5, dpi = 300)
  }

  return(p)
}


df_combined_plot_beta=read.csv(here::here("simulation","Results","betacsv","betatest_sensitivity.csv"))
p_beta <- plot_noise_allB_group_method(
  noise_val         = 1,
  B_vals            = c(1),
  v_vals            = 1:3,
  group_method_vals = target_methods,
  out_png           = here::here("simulation","Results","betacsv","FigS28.png"),
  all_data_df       = df_combined_plot_beta
)
print(p_beta)
