library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(xtable)
library(ggplot2)
library(here)

all_data<-read.csv(here::here("simulation","Results","betacsv","betatest_mainsetting.csv"))[,-1]

plot_noise_allB <- function(noise_val, B_vals = NULL, v_vals = 1:3, F_val = NA, out_png = NULL) {

  sub <- all_data %>% filter(noise == noise_val)
  if (!is.null(B_vals)) {
    sub <- sub %>% filter(B %in% B_vals)
  } else {
    B_vals <- sort(unique(sub$B))
  }


  if (is.na(F_val)) {
    F_val <- unique(sub$F)[1]
  }

  df <- sub %>% filter(F == F_val, V %in% v_vals)

  df_sum <- df %>%
    group_by(B, V, method, k, w) %>%
    summarise(pval = mean(pval, na.rm = TRUE), .groups = "drop")

  v_labeller <- function(value) paste0("V", value)

  p <- ggplot(df_sum, aes(x = w, y = pval,
                          color = method,
                          shape = factor(k),
                          linetype = factor(k),
                          group = interaction(method, k))) +
    geom_line(size = 0.8, na.rm = TRUE) +
    geom_point(size = 4, na.rm = TRUE) +

    facet_grid(B ~ V, labeller = labeller(V = as_labeller(v_labeller))) +
    labs(
      x = expression(w),
      y = "Empirical rejection rate",
      shape = "Modality",
      linetype = "Modality",
      color = "Method"
    ) +
    theme_bw(base_size = 14) +
    theme(
      strip.text = element_blank(),
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 16),
      legend.key.size = unit(1.2, "lines"),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 16),
      legend.position = "bottom",
      legend.spacing.x = unit(0.3, "cm"),
      legend.spacing.y = unit(0.1, "cm"),
      panel.grid.minor = element_blank(),

      legend.box.background = element_blank()
    ) +
    guides(
      shape = guide_legend(override.aes = list(size = 4)),
      color = guide_legend(override.aes = list(size = 4))
    )

  if(!is.null(out_png)){
    ggsave(out_png, p, width = 16, height = 9, dpi = 300)
  }

  return(p)
}


  plot_noise_allB(noise_val = 1, B_vals = c(1,2), out_png = here::here("simulation", "Results", "betacsv","Fig3.png"))
  plot_noise_allB(noise_val = 2, B_vals = c(1,2), out_png = here::here("simulation", "Results", "betacsv","FigS3.png"))

