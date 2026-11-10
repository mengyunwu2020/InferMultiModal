library(openxlsx)
library(readxl)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(here)
load(here::here("realdata", "Results", "fit_allsample.Rdata"))
cat("pvals for betatest:",'\n')
print(fit_ori$pvals[1:5])
load(here::here("realdata", "Results", "gamma_allsample.Rdata"))
cat("pvals for gammatest:",'\n')
print(gres_ori[1:5])


raw_p_values <- fit_ori$pvals[-c(1:5)]
p_fdr <- p.adjust(raw_p_values, method = "fdr")
results <- data.frame(Raw = raw_p_values,  FDR = p_fdr)

sig_fdr <- which(results$FDR[301:900] < 0.05)


load(here::here("realdata", "Results","name_gene.Rdata"))

pd=fit_ori$fit[["FV"]]
rrhat=fit_ori$fit$rr_hat
Fhat=pd[,1:(sum(rrhat))]
V=pd[,(sum(rrhat)+1):ncol(pd)]
sdv=apply(V,2,sd)
sdf=apply(Fhat,2,sd)


betahat <- fit_ori$fit$beta
pm <- c(300, 300, 300, 194)
ss <- 1

wb <- createWorkbook()
sheetname=c('beta_mRNA','beta_Methylation','beta_CNA','beta_ProteinExpression')
for (kk in 1:4) {
  betatmp <- betahat[ss:(ss + pm[kk] - 1)]
  pvtmp=results$FDR[ss:(ss + pm[kk] - 1)]
  sdv1=sdv[ss:(ss + pm[kk] - 1)]
  names(betatmp) <- name_gene[[kk]][1:pm[kk]]
  ss <- ss + pm[kk]

  df <- data.frame(Gene = names(betatmp), coefficient = as.numeric(betatmp),sd.x=sdv1,p.value=pvtmp)

  addWorksheet(wb, sheetname[kk])
  writeData(wb, sheet = kk, df)
}
gammahat=fit_ori$fit$gamma
rrhat=fit_ori$fit$rr_hat
df <- data.frame(factor = 1, coefficient = as.numeric(gammahat[1:rrhat[1]]))

addWorksheet(wb, 'gamma_commonFactor')
writeData(wb, sheet = 5, df)

ss=rrhat[1]+1
sheetname=c('gamma_s_mRNA','gamma_s_Methylation','gamma_s_CNA','gamma_s_ProteinExpression')
for (kk in 1:4) {
  betatmp <- gammahat[ss:(ss + rrhat[(kk+1)] - 1)]
  names(betatmp) <-1:rrhat[(kk+1)]
  ss <- ss + rrhat[(kk+1)]

  df <- data.frame(factor = names(betatmp), coefficient = as.numeric(betatmp))

  addWorksheet(wb, sheetname[kk])
  writeData(wb, sheet = kk+5, df)
}

saveWorkbook(wb, file = here::here("realdata","Results","Realdata.xlsx"), overwrite = TRUE)
