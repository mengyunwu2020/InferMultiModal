library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)
library(dplyr)
library(here)

load(here::here("realdata", "Results", "fit_allsample.Rdata"))
fit1 <- fit_ori$fit
load(here::here("realdata", "Results", "name_gene.RData"))

p_go <- c(30, 30, 30, 30)
ontologies <- c("BP", "CC", "MF")

compute_go_terms <- function(genes) {
  genes <- unique(trimws(unlist(strsplit(genes, "\\s+"))))

  gene_entrez <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop = TRUE)
  gene_list <- gene_entrez$ENTREZID

  if (length(gene_list) == 0) return(NULL)

  go_list <- list()
  for (ont in ontologies) {
    ego <- tryCatch({
      enrichGO(gene          = gene_list,
               OrgDb         = org.Hs.eg.db,
               keyType       = "ENTREZID",
               ont           = ont,
               pAdjustMethod = "BH",
               pvalueCutoff  = 0.05,
               qvalueCutoff  = 0.2,
               readable      = TRUE)
    }, error = function(e) NULL)

    if (!is.null(ego)) {
      ego_df <- as.data.frame(ego)
      if (nrow(ego_df) > 0) {
        ego_df$Ontology <- ont
        go_list[[ont]] <- ego_df
      }
    }
  }

  if (length(go_list) > 0) {
    return(do.call(rbind, go_list) %>% arrange(p.adjust))
  } else {
    return(NULL)
  }
}

list_to_compact_table <- function(my_list,
                                  file = here::here("realdata", "Results", "Tab2.tex"),
                                  table_caption = "Data analysis: top five significantly enriched GO terms for the identified modality-shared factor and each modality-specific factor and their adjusted $p$-values (p.adjust).",
                                  table_label = "tab:GO",
                                  row_select = 5) {

  fileConn <- file(file, open = "wt")
  on.exit(close(fileConn))

  writeLines(c("\\begin{table}[htbp]",
               "\\centering",
               paste0("\\caption{", table_caption, "}"),
               paste0(" \\label{", table_label, "} "),
               "\\renewcommand\\tabcolsep{1pt}",
               "\\scriptsize",
               "\\begingroup",
               "\\renewcommand{\\arraystretch}{.56}   ",
               "\\begin{tabular}{p{2.2cm}p{2cm}p{9.5cm}p{2cm}}  ",
               "\\toprule",
               "\\textbf{Factor} & \\textbf{GO term} & \\textbf{Description} & \\textbf{p.adjust} \\\\",
               "\\midrule "), fileConn)

  valid_keys <- names(my_list)

  for (k_idx in seq_along(valid_keys)) {
    nm <- valid_keys[k_idx]
    mat <- my_list[[nm]]
    if (is.null(mat) || nrow(mat) == 0) next

    mat <- mat[1:min(row_select, nrow(mat)), , drop = FALSE]
    n_rows <- nrow(mat)

    if (nm == "common") {
      factor_label <- "\\shortstack[l]{Modality-shared\\\\factor }"
    } else if (grepl("^mrna", nm)) {
      idx <- sub("^mrna", "", nm)
      factor_label <- paste0("\\shortstack[l]{mRNA-specific\\\\factor ", idx, "}")
    } else if (grepl("^cna", nm)) {
      idx <- sub("^cna", "", nm)
      factor_label <- paste0("\\shortstack[l]{CNA-specific\\\\factor ", idx, "}")
    } else {
      factor_label <- nm
    }

    for (i in 1:n_rows) {
      row_vals <- as.character(mat[i, ])

      row_vals[2] <- gsub("([%_&])", "\\\\\\1", row_vals[2])
      row_vals[2] <- gsub("(^[a-z])", "\\U\\1", row_vals[2], perl = TRUE)

      p_num <- suppressWarnings(as.numeric(row_vals[3]))
      if (!is.na(p_num) && p_num > 0) {
        exp_val  <- floor(log10(abs(p_num)))
        base_val <- round(p_num / 10^exp_val, 2)
        row_vals[3] <- paste0("$", sprintf("%.2f", base_val), "\\times10^{", exp_val, "}$")
      } else if (is.na(p_num)) {
        row_vals[3] <- "--"
      }

      lead_col <- if (i == 1) paste0("\\multirow{", n_rows, "}{*}{", factor_label, "}") else ""
      line_str <- paste(c(lead_col, row_vals), collapse = " & ")
      writeLines(paste0(line_str, " \\\\"), fileConn)
    }

    if (k_idx < length(valid_keys)) {
      writeLines("\\midrule", fileConn)
    } else {
      writeLines("\\bottomrule", fileConn)
    }
  }

  writeLines(c("\\end{tabular}",
               "\\endgroup",
               "\\end{table}"), fileConn)
}

GO <- list()

select_gene <- list()
for (ii in 1:4) {
  bm <- abs(fit1[["loading_c"]][[ii]])
  select_gene[[ii]] <- name_gene[[ii]][order(bm, decreasing = TRUE)[1:p_go[ii]]]
  if (ii == 4) {
    genes_raw <- sub("\\|.*", "", select_gene[[ii]])
    genes <- unlist(strsplit(genes_raw, "\\s+"))
    select_gene[[ii]] <- unique(genes)
  }
}
res_common <- unlist(select_gene)
res_df <- compute_go_terms(res_common)
if (!is.null(res_df) && nrow(res_df) > 0) GO[["common"]] <- res_df

n_mrna_factors <- ncol(fit1[["loading_s"]][[1]])
for (rr in 1:n_mrna_factors) {
  select_gene <- list()
  for (ii in 1:4) {
    bm <- abs(fit1[["loading_s"]][[1]][, rr])
    select_gene[[ii]] <- name_gene[[ii]][order(bm, decreasing = TRUE)[1:p_go[ii]]]
    if (ii == 4) {
      genes_raw <- sub("\\|.*", "", select_gene[[ii]])
      genes <- unlist(strsplit(genes_raw, "\\s+"))
      select_gene[[ii]] <- unique(genes)
    }
  }
  res_mrna <- unlist(select_gene)
  res_df <- compute_go_terms(res_mrna)
  if (!is.null(res_df) && nrow(res_df) > 0) GO[[paste0('mrna', rr)]] <- res_df
}

n_cna_factors <- ncol(fit1[["loading_s"]][[3]])
for (rr in 1:n_cna_factors) {
  select_gene <- list()
  for (ii in 1:4) {
    bm <- abs(fit1[["loading_s"]][[3]][, rr])
    select_gene[[ii]] <- name_gene[[ii]][order(bm, decreasing = TRUE)[1:p_go[ii]]]
    if (ii == 4) {
      genes_raw <- sub("\\|.*", "", select_gene[[ii]])
      genes <- unlist(strsplit(genes_raw, "\\s+"))
      select_gene[[ii]] <- unique(genes)
    }
  }
  res_cna <- unlist(select_gene)
  res_df <- compute_go_terms(res_cna)
  if (!is.null(res_df) && nrow(res_df) > 0) GO[[paste0('cna', rr)]] <- res_df
}

results_dir <- here::here("realdata", "Results")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

GO_all_combined <- bind_rows(lapply(names(GO), function(nm) {
  df <- GO[[nm]]
  if (is.null(df) || nrow(df) == 0) return(NULL)
  data.frame(Factor = nm, df, stringsAsFactors = FALSE)
}))
write.csv(GO_all_combined, file.path(results_dir, "GO_all_factors.csv"), row.names = FALSE)

output_tex_path <- file.path(results_dir, "TabS6.tex")
list_to_compact_table(lapply(GO, function(df) df[, c("ID", "Description", "p.adjust")]),
                      file = output_tex_path,
                      table_caption = "Data analysis: top five significantly enriched GO terms for the identified modality-shared factor and each modality-specific factor and their adjusted $p$-values (p.adjust).",
                      table_label = "tab:GO",
                      row_select = 5)
