library(ggplot2)
library(tidyr)
library(dplyr)
library(ggh4x)
library(here)
v_labels <- c("1" = "V1", "2" = "V2", "3" = "V3")

df <- read.csv(here::here("simulation","Results","estimationcsv","estimation_lowernoise.csv"))[,-1]

df_long <- df %>%
  pivot_longer(
    cols = c(Proposed, MCP, FARM, IFRM,Joint),
    names_to = "Method",
    values_to = "Value"
  ) %>%
  filter(!is.na(Value)) %>%
  mutate(setting_label = paste0("F", F, "_V", V, "_noise", noise, "_B", B)) %>%
  arrange(p, noise, B, V, Method, xval)
for (p_val in 1200) {
  df_p <- df_long %>% filter(p == p_val)
  nb_combos <- df_p %>% distinct(noise)

  for (i in 1:nrow(nb_combos)) {
    noise_val <- nb_combos$noise[i]
    df_sub <- df_p %>% filter(noise == noise_val)

    df_sub <- df_sub %>%
      group_by(B, V) %>%
      mutate(
        min_val = min(Value, na.rm = TRUE),
        max_val = max(Value, na.rm = TRUE),
        ymin_BV = ifelse(is.finite(min_val), 0.8 * min_val, NA_real_),
        ymax_BV = ifelse(is.finite(max_val), .4 * max_val, NA_real_)
      ) %>%
      ungroup() %>%
      dplyr::select(-min_val, -max_val)

    p_plot <- ggplot(df_sub, aes(x = xval, y = Value, color = Method, group = Method, shape = Method)) +

      geom_point(size = 3.5, alpha = 0.8) +

      geom_smooth(
        data = subset(df_sub, Method == "Proposed"),
        method = "gam",
        formula = y ~ s(x, k = 2),
        se = FALSE,
        size = 1.2
      ) +
      geom_smooth(
        data = subset(df_sub, Method != "Proposed"),
        method = "loess",
        span = 0.5,
        se = FALSE,
        size = 1.2
      ) +

    geom_blank(aes(y = ymin_BV)) +
      geom_blank(aes(y = ymax_BV)) +
      # =======================================================================

    scale_shape_manual(values = c(15, 16, 17, 18, 8, 4, 3)) +

      facet_grid(
        B ~ V,
        scales = "free_y",
        labeller = labeller(V = v_labels, B = label_value)
      ) +

      labs(y = "Average L2 error", color = "Method", shape = "Method") +
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
        color = guide_legend(override.aes = list(size = 5)),
        shape = guide_legend(override.aes = list(size = 5))
      ) +

      xlab(expression(2 * sqrt(s * log(p) / n))) +
      xlim(0.25, 0.6)

    ggsave(
      filename = here::here("simulation", "Results", "estimationcsv","FigS22.png"),
      plot = p_plot,
      width = 16,
      height = 9
    )
  }
}



#------plot prediction-------
df <- read.csv(here::here("simulation","Results","estimationcsv","prediction_lowernoise.csv"))[,-1]
my_shapes <- c(15, 16, 17, 18, 8, 4, 3)
df_long <- df %>%
  pivot_longer(
    cols = c(Proposed, MCP, FARM, IFRM,Joint),
    names_to = "Method",
    values_to = "Value"
  ) %>%
  filter(!is.na(Value)) %>%
  mutate(setting_label = paste0("F", F, "_V", V, "_noise", noise, "_B", B)) %>%
  arrange(p, noise, B, V, Method, xval)
for (p_val in 1200) {
  df_p <- df_long %>% filter(p == p_val)

  nb_combos <- df_p %>% distinct(noise)

  for (i in 1:nrow(nb_combos)) {
    noise_val <- nb_combos$noise[i]
    df_sub <- df_p %>% filter(noise == noise_val)

    y_min <- min(df_sub$Value, na.rm = TRUE)

    y_max <- ifelse(noise_val == 1, .3, .3)

    current_span <- ifelse(i == 1, 0.6, 0.45)

    df_prop <- df_sub %>%
      filter(Method == "Proposed") %>%
      group_by(B, V) %>%
      group_modify(~ {
        if(nrow(.x) > 3) {
          model_gam <- mgcv::gam(Value ~ s(xval, k = 2), data = .x)
          .x$pred_value <- predict(model_gam, newdata = .x)
        } else {
          .x$pred_value <- .x$Value
        }
        return(.x)
      }) %>%
      ungroup()

    df_others <- df_sub %>%
      filter(Method != "Proposed") %>%
      group_by(Method, B, V) %>%
      group_modify(~ {
        if(nrow(.x) > 3) {
          model_loess <- loess(Value ~ xval, data = .x, span = current_span)
          .x$pred_value <- predict(model_loess, newdata = .x)
        } else {
          .x$pred_value <- .x$Value
        }
        return(.x)
      }) %>%
      ungroup()

    df_plot_data <- bind_rows(df_prop, df_others)


    p_plot <- ggplot(df_plot_data, aes(x = xval, color = Method, group = Method)) +
      geom_line(aes(y = pred_value, linetype = Method), linewidth = 1.2) +
      geom_point(aes(y = pred_value, shape = Method), size = 3.5, alpha = 0.8) +
      scale_shape_manual(values = my_shapes) +
      facet_grid(B ~ V, scales = "free_y", labeller = labeller(V = v_labels, B = label_both)) +

      facetted_pos_scales(
        y = list(
          B == 1 ~ scale_y_continuous(limits = c(0, .26), n.breaks = 5),
          B == 2 ~ scale_y_continuous(limits = c(0, .26), n.breaks = 5)
        )
      ) +

      labs(y = "Prediction error", color = "Method", shape = "Method", linetype = "Method") +
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
      xlab(expression(2*sqrt(s*log(p)/n))) +
      coord_cartesian(xlim = c(0.25, 0.6))
    ggsave(
      filename = here::here("simulation", "Results", "estimationcsv","FigS23.png"),
      plot = p_plot,
      width = 16,
      height = 9
    )
  }
}
