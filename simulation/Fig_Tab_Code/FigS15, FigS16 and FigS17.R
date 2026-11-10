library(dplyr)
library(ggplot2)
library(here)
reshape_method <- function(df, method_name, modality_name) {
  df2 <- df %>%
    group_by(V, noise, B) %>%
    mutate(
      setting = cur_group_id(),
      v = F,
      k = modality_name,
      method = method_name
    ) %>%
    ungroup() %>%
    dplyr::select(setting, F, V, noise, B, v, k, pval, method)

  return(df2)
}

plot_noise_allB_from_table <- function(noise_val, B_vals = NULL, v_vals = 1:3, method_vals = NULL, out_png = NULL, all_data_df) {

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

  if(nrow(df) == 0){
    stop("No Data")
  }

  df_sum <- df %>%
    mutate(
      pval = as.numeric(pval),
      v = as.numeric(v)
    ) %>%
    group_by(B, V, method, k, v) %>%
    summarise(pval = mean(pval, na.rm = TRUE), .groups = "drop") %>%
    mutate(pval = ifelse(is.nan(pval), NA, pval)) %>%

    filter(!is.na(pval)) %>%

    arrange(B, V, method, k, v) %>%
    group_by(B, V, method, k) %>%
    mutate(pval = cummax(pval)) %>%
    ungroup()

  if(nrow(df_sum) == 0){
    stop("No Data")
  }

  target_v <- sort(unique(df_sum$v))
  v_labeller <- function(value) paste0("V", value)

  method_labels <- c(
    "proposed" = "Proposed",
    "Joint" = "Joint",
    "jointadjust" = "Jointa"
  )

  modality_labels <- c("1" = "1", "2" = "2", "3" = "3", "4" = "4", "common" = "common")

  p <- ggplot(df_sum, aes(x = v, y = pval,
                          color = method,
                          shape = factor(k),
                          linetype = factor(k),
                          group = interaction(method, k))) +
    geom_line(linewidth = 0.8, na.rm = TRUE) +
    geom_point(size = 4, na.rm = TRUE) +

    facet_grid(B ~ V, labeller = labeller(V = as_labeller(v_labeller))) +
    scale_x_continuous(breaks = target_v) +

    labs(x = "Signal", y = "Empirical rejection rate", shape = "Modality", linetype = "Modality", color = "Method") +
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
      panel.grid.minor = element_blank()
    ) +

    scale_color_discrete(labels = function(x) {
      mapped <- method_labels[x]
      mapped[is.na(mapped)] <- x[is.na(mapped)]
      parse(text = mapped)
    }) +
    scale_shape_discrete(labels = modality_labels) +
    scale_linetype_discrete(labels = modality_labels)

  if(!is.null(out_png)){
    ggsave(out_png, p, width = 16, height = 9, dpi = 300)
  }

  return(p)
}


final_data= read.csv(here::here("simulation","Results","gammacsv","gammatest_joint_mainsetting.csv"))[,-1]

selected_methods <- c("Proposed","Proposed~(bold(gamma[c]))","Joint","Joint ~ (bold(gamma[c]))")
for(ii in 1) {
  dd = plot_noise_allB_from_table(
    noise_val = ii,
    B_vals = c(1,2),
    method_vals = selected_methods,
    out_png = here::here("simulation", "Results", "gammacsv","FigS15.png"),
    all_data_df = final_data
  )
  print(dd)
}
selected_methods <- c("Proposed","Proposed~(bold(gamma[c]))","oracle_r_based","oracle_r_based ~ (bold(gamma[c]))")
for(ii in 1) {
  dd = plot_noise_allB_from_table(
    noise_val = ii,
    B_vals = c(1,2),
    method_vals = selected_methods,
    out_png = here::here("simulation", "Results", "gammacsv","FigS16.png"),
    all_data_df = final_data
  )
  print(dd)
}
selected_methods <- c("Proposed","Proposed~(bold(gamma[c]))","oracle_based","oracle_based ~ (bold(gamma[c]))")

for(ii in 1) {
  dd = plot_noise_allB_from_table(
    noise_val = ii,
    B_vals = c(1,2),
    method_vals = selected_methods,
    out_png = here::here("simulation", "Results", "gammacsv","FigS17.png"),
    all_data_df = final_data
  )
  print(dd)
}
