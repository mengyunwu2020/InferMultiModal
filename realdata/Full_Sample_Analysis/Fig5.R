library(gridExtra)
library(ggplot2)
library(dplyr)
library(here)

load(here::here("realdata", "Results", "fit_allsample.Rdata"))
cat("pvals for betatest:", '\n')
print(fit_ori$pvals[1:5])
load(here::here("realdata", "Results", "gamma_allsample.Rdata"))
cat("pvals for gammatest:", '\n')
print(gres_ori[1:5])

raw_p_values <- fit_ori$pvals[-c(1:5)]
p_fdr <- p.adjust(raw_p_values, method = "fdr")
results <- data.frame(Raw = raw_p_values, FDR = p_fdr)

sig_fdr <- which(results$FDR[301:900] < 0.05)

load(here::here("realdata", "Results", "name_gene.Rdata"))

pd = fit_ori$fit[["FV"]]
rrhat = fit_ori$fit$rr_hat
Fhat = pd[, 1:(sum(rrhat))]
V = pd[, (sum(rrhat) + 1):ncol(pd)]
sdv = apply(V, 2, sd)
sdf = apply(Fhat, 2, sd)

betahat <- fit_ori$fit$beta
pm <- c(300, 300, 300, 194)
ss <- 1

sheetname = c('beta_mRNA', 'beta_Methylation', 'beta_CNA', 'beta_ProteinExpression')
beta_tables <- list()

for (kk in 1:4) {
  betatmp <- betahat[ss:(ss + pm[kk] - 1)]
  pvtmp = results$FDR[ss:(ss + pm[kk] - 1)]
  sdv1 = sdv[ss:(ss + pm[kk] - 1)]
  names(betatmp) <- name_gene[[kk]][1:pm[kk]]
  ss <- ss + pm[kk]

  df <- data.frame(Gene = names(betatmp), coefficient = as.numeric(betatmp), sd.x = sdv1, p.value = pvtmp)
  beta_tables[[sheetname[kk]]] <- df
}

# --- 直接从内存读取 Methylation 数据 ---
dat <- beta_tables[["beta_Methylation"]]

methylation <- data.frame(
  gene = dat$Gene,
  coef = dat$coefficient,
  p = dat$p.value,
  sd = dat$sd.x
)[dat$p.value < 0.05, ]

p1 = ggplot(methylation, aes(x = reorder(gene, -abs(coef)), y = coef)) +
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

# --- 直接从内存读取 CNA 数据 ---
dat <- beta_tables[["beta_CNA"]]

cna <- data.frame(
  gene = dat$Gene,
  coef = dat$coefficient,
  p = dat$p.value,
  sd = dat$sd.x
)[dat$p.value < 0.05, ]

p2 = ggplot(cna, aes(x = reorder(gene, -abs(coef)), y = coef)) +
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
    x = NULL,
    y = "Coefficient"
  )
