library(ggplot2)
library(dplyr)
library(tidyr)
library(mgcv)
library(ggh4x)
library(here)
v_labels <- c("1" = "V1", "2" = "V2", "3" = "V3")
label_value <- c("1" = "Gaussian noise", "2" = "Uniform noise")
my_shapes <- c(15, 16, 17, 18, 8, 4, 3)

df <- read.csv(here::here("simulation","Results","estimationcsv","prediction_mainsetting.csv"))[, -1]

df_long <- df %>%
  pivot_longer(
    cols = c(Proposed, MCP, FARM, IFRM),
    names_to = "Method",
    values_to = "Value"
  ) %>%
  filter(!is.na(Value))

for (p_val in 1200) {
  df_p <- df_long %>% filter(p == p_val)
  combos <- df_p %>% distinct(F, B)

  for (j in 1:nrow(combos)) {
    f_val <- combos$F[j]
    b_val <- combos$B[j]

    df_sub <- df_p %>% filter(F == f_val, B == b_val)
    df_prop <- df_sub %>% filter(Method == "Proposed")

    if (nrow(df_prop) > 0) {
      df_prop <- df_prop %>%
        group_by(V, noise) %>%
        group_modify(~{
          if (nrow(.x) > 3) {
            fit <- tryCatch(
              mgcv::gam(Value ~ s(xval, k = 3), data = .x),
              error = function(e) NULL
            )
            if (!is.null(fit)) {
              .x$pred_value <- predict(fit, newdata = .x)
            } else {
              .x$pred_value <- .x$Value
            }
          } else {
            .x$pred_value <- .x$Value
          }
          .x
        }) %>%
        ungroup()
    }

    df_other <- df_sub %>% filter(Method != "Proposed")

    if (nrow(df_other) > 0) {
      df_other <- df_other %>%
        group_by(Method, V, noise) %>%
        group_modify(~{
          if (nrow(.x) > 3) {
            fit <- tryCatch(
              loess(Value ~ xval, data = .x, span = 0.5),
              error = function(e) NULL
            )
            if (!is.null(fit)) {
              .x$pred_value <- predict(fit, newdata = .x)
            } else {
              .x$pred_value <- .x$Value
            }
          } else {
            .x$pred_value <- .x$Value
          }
          .x
        }) %>%
        ungroup()
    }

    df_plot_data <- bind_rows(df_prop, df_other)

    p_plot <- ggplot(
      df_plot_data,
      aes(x = xval, y = pred_value, color = Method, group = Method)
    ) +
      geom_line(aes(linetype = Method), linewidth = 1.2) +
      geom_point(aes(shape = Method), size = 3.5, alpha = 0.85) +
      scale_shape_manual(values = my_shapes) +
      facet_grid(
        noise ~ V,
        scales = "free_y",
        labeller = labeller(V = v_labels, noise = label_value)
      ) +
      ggh4x::facetted_pos_scales(
        y = list(
          scale_y_continuous(limits = c(1, 1.5), breaks = seq(1, 1.5, by = 0.1), expand = c(0, 0)),
          scale_y_continuous(limits = c(1, 2.10), breaks = seq(1, 2.10, by = 0.2), expand = c(0, 0))
        )
      ) +
      labs(
        y = "Prediction error",
        color = "Method",
        shape = "Method",
        linetype = "Method"
      ) +
      theme_bw() +
      theme(
        legend.position = "bottom",
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        strip.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        strip.text.y = element_blank(),
        strip.background.y = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10))
      ) +
      guides(
        color = guide_legend(override.aes = list(size = 4)),
        shape = guide_legend(override.aes = list(size = 4)),
        linetype = guide_legend(keywidth = 3)
      ) +
      xlab(expression(2 * sqrt(s * log(p) / n))) +
      xlim(0.25, 0.6)

    out_file <- if (b_val == 1) here::here("simulation", "Results", "estimationcsv","Fig2.png") else here::here("simulation", "Results", "estimationcsv","FigS2.png")
    ggsave(filename = out_file, plot = p_plot, width = 16, height = 9, dpi = 300)
  }
}
