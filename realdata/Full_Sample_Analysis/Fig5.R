library(openxlsx)
library(readxl)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(here)

dat <- read_excel(
  here::here("realdata", "Results", "Realdata.xlsx"),
  sheet = "beta_Methylation"
)

library(ggplot2)

methylation <- data.frame(
  gene = dat$Gene,
  coef = dat$coefficient,
  p = dat$p.value,
  sd = dat$sd.x
)[dat$p.value<0.05,]
library(ggplot2)


p1=ggplot(methylation, aes(x = reorder(gene, -abs(coef)), y = coef)) +

  geom_col(aes(fill = p), width = 0.7, color = "black", size = 0.3) +

  geom_errorbar(aes(ymin = coef - sd, ymax = coef + sd),
                width = 0.2, alpha = 0.6) +

  geom_text(aes(label = sprintf("%.2f", sd),
                vjust = ifelse(coef >= 0, -1.5, 1.5)),
            size = 3.5, color = "black") +

  scale_fill_gradient(low = "#E41A1C", high = "grey90", name = "P-value") +

  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "right"
  ) +

  labs(title = "Methylation Analysis",
       subtitle = "Bar Height = Coefficient | Color = P-value | Label = SD",
       x = "Gene",
       y = "Coefficient")



dat <- read_excel(
  here::here("realdata", "Results", "Realdata.xlsx"),
  sheet = "beta_CNA"
)

cna <- data.frame(
  gene = dat$Gene,
  coef = dat$coefficient,
  p = dat$p.value,
  sd = dat$sd.x
)[dat$p.value<0.05,]
p2=ggplot(cna, aes(x = reorder(gene, -abs(coef)), y = coef)) +

  geom_col(aes(fill = p), width = 0.6, color = "black", size = 0.3) +

  geom_errorbar(aes(ymin = coef - sd, ymax = coef + sd),
                width = 0.2, alpha = 0.6) +

  geom_text(aes(label = sprintf("%.2f", sd),
                vjust = ifelse(coef >= 0, -1.5, 1.5)),
            size = 4, color = "black") +

  scale_fill_gradient(low = "#E41A1C", high = "grey90", name = "P-value") +

  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "right"
  ) +

  labs(title = "CNA Analysis",
       subtitle = "Bar Height = Coefficient | Color = P-value | Label = SD",
       x = "Gene",
       y = "Coefficient")



grid.arrange(p1, p2, ncol = 2)



methylation$Type <- "DNA methylation"
cna$Type <- "CNA"
combined_data <- rbind(methylation, cna)

ggplot(combined_data, aes(x = reorder(gene, -abs(coef)), y = coef)) +

  geom_col(aes(fill = p), width = 0.7, color = "black", size = 0.3) +

  geom_errorbar(aes(ymin = coef - sd, ymax = coef + sd),
                width = 0.3, alpha = 0.6) +

   geom_text(aes(label = sprintf("%.2f", sd),
                vjust = ifelse(coef >= 0, -2, 2.5)),
            size = 3, color = "black") +

    facet_grid(. ~ Type, scales = "free_x", space = "free_x") +

   scale_y_continuous(expand = expansion(mult = 0.2)) +

   scale_fill_gradient(low = "#E41A1C", high = "grey90", name = "p-value") +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +

  theme_bw() +
  theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold", color = "black"),
      strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(size = 12, face = "bold"),
     panel.grid.major.x = element_blank()
  ) +

  labs(
    x =NULL,
    y = "Coefficient")

