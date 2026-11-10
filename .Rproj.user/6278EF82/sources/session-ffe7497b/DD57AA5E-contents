library(here)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(zoo)
library(stringr)
library(ggh4x)
library(mgcv)
df_long_interp_est <- read.csv(
  here::here("simulation", "Results", "estimationcsv","estimation_sensitivity.csv"),check.names = FALSE
)[, -1]
target_methods <- c("original", "common1", "spec1", "all3", "wrong_r")

v_labels    <- c("1" = "V1", "2" = "V2", "3" = "V3")
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

my_shapes <- setNames(c(15, 16, 17, 18, 8), target_methods)
high_contrast_colors <- c(
  "original" = "#000000",
  "common1"  = "#E41A1C",
  "spec1"    = "#377EB8",
  "all3"     = "#984EA3", #
  "wrong_r"  = "#4DAF4A"
)

p_values_est <- unique(df_long_interp_est$p)
for (p_val in p_values_est) {
  df_p <- df_long_interp_est %>% filter(p == p_val)
  nb_combos <- df_p %>% distinct(noise)

  for (i in 1:nrow(nb_combos)) {
    noise_val <- nb_combos$noise[i]
    df_sub <- df_p %>% filter(noise == noise_val)

    # LOESS 平滑
    df_plot_data <- df_sub %>%
      group_by(B, V, group_method) %>%
      group_modify(~ {
        if (nrow(.x) > 4) {
          fit <- loess(Value ~ xval, data = .x, span = 0.6)
          .x$pred_value <- predict(fit, newdata = .x)
        } else {
          .x$pred_value <- .x$Value
        }
        return(.x)
      }) %>%
      ungroup()

    # 动态 Y 轴
    df_plot_data <- df_plot_data %>%
      group_by(B, V) %>%
      mutate(
        min_val = min(pred_value, na.rm = TRUE),
        max_val = max(pred_value, na.rm = TRUE),
        ymin_BV = ifelse(is.finite(min_val), 0.85 * min_val, NA_real_),
        ymax_BV = ifelse(is.finite(max_val), 1.15 * max_val, NA_real_)
      ) %>%
      ungroup() %>%
      dplyr::select(-min_val, -max_val)

    p_est <- ggplot(df_plot_data, aes(x = xval, y = pred_value,
                                      color = group_method,
                                      shape = group_method,
                                      group = group_method)) +
      geom_line(linewidth = 1.0) +
      geom_point(size = 3.2, alpha = 0.85) +
      geom_blank(aes(y = ymin_BV)) +
      geom_blank(aes(y = ymax_BV)) +
      scale_color_manual(
        name   = "Method",
        breaks = target_methods,
        values = high_contrast_colors[target_methods],
        labels = parsed_method_labels[target_methods]
      ) +
      scale_shape_manual(
        name   = "Method",
        breaks = target_methods,
        values = my_shapes[target_methods],
        labels = parsed_method_labels[target_methods]
      ) +
      facet_grid(
        B ~ V,
        scales = "free_y",
        labeller = labeller(V = v_labels, B = label_value)
      ) +
      labs(y = "Average L2 error") +
      theme_bw() +
      theme(
        legend.position = "bottom",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 15, face = "bold"),
        strip.text = element_text(size = 15),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        strip.text.y = element_blank(),
        strip.background.y = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.spacing.x = unit(0.4, "cm")
      ) +
      guides(
        color = guide_legend(
          title = "Method",
          override.aes = list(
            size      = 3.5,
            linewidth = 1.0,
            shape     = unname(my_shapes[target_methods])
          )
        ),
        shape = "none"
      ) +
      xlab(expression(2 * sqrt(s * log(p) / n))) +
      coord_cartesian(xlim = c(0.25, 0.6))

    out_file <- here::here("simulation", "Results", "estimationcsv","FigS26.png")
    ggsave(filename = out_file, plot = p_est, width = 16, height = 9, dpi = 300)
  }
}
