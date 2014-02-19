HW2_aengstrom28
===============
```{r}
# Getting data
> source("http://bioconductor.org/biocLite.R")
> library(Biobase)
> library(data.table)
> library(GEOquery)
> library(limma)
> gse_new <- getGEO(filename = "/Users/jjengstrom/Biostat-578/Data/GSE40812_series_matrix.txt.gz")
> setwd("~/aengstrom28")

#Cleaning-up data
> gds <- gse_new
> pd <- pData(gds)
> colnames(pd)
#Using monocyte-derived macrophage data from pData
> mm_pd <- pd[pd$source_name_ch1=="Monocyte-derived Macrophage",]
> mm_eset <- gds[,rownames(mm_pd)]
#Change "characteristics_ch1" to HCV status (viral load + or -)
> mm_pd$HCV <- gsub(".*: ", "", mm_pd$characteristics_ch1)
> mm_pd$HCV <- ifelse(mm_pd$HCV=="Neg", "-", "+")
#Change "characteristics_ch1.2" to treatment (Mock or Poly IC) and characteristics_ch1.1 to cell (cell-type is macrophage)
> mm_pd$treatment <- gsub(".*: ", "", mm_pd$characteristics_ch1.2)

#Set-up design matrix and linear model
> mm_matrix <-model.matrix(~treatment+HCV,mm_pd)
> colnames(mm_matrix) <- c("Neg", "Pos")
> fit_mm_matrix <- lmFit(mm_eset, mm_matrix)
> ebay_mm_matrix <- eBayes(fit_mm_matrix)
> TopHCV <- topTable(ebay_mm_matrix,adjust="BH")
> colnames(TopHCV)
> colnames(fit_mm_matrix$coef)
> TopHCV[,35]

#Set-up contrasts
> cont_matrix <- makeContrasts(Neg-Pos,levels=mm_matrix)
> rownames(cont_matrix)
> colnames(fit_mm_matrix$coef)
> colnames(fit_mm_matrix$coef) <- c("Neg", "Pos")
> fit2 <- contrasts.fit(fit_mm_matrix, cont_matrix)
> ebay_fit2 <- eBayes(fit2)
> cont_matrix
> colnames(ebay_fit2)

#
> model.matrix ~ treatment+VL
> set-up contrast from model.matrix for treatment between poly-mock treated to find probes that are differentially expressed in the poly treated group
> do a topTable on the fit for treatment+VL
colnames(TopHCV
> sum(TopHCV$adj.p < 0.05)
eset_small <- eSet[TopHCV$adj.p <0.05,]
> set-up contrast for HCV+ - HCV-
>may need to do linear model pulling from new eset
> then set-up new contrast matrix using VL- VL+ form old matrix)

> mm_pd <- pd[pd$source_name_ch1=="Monocyte-derived Macrophage",]
> mm_eset <- gds[,rownames(mm_pd)]
> mm_pd$HCV <- gsub(".*: ", "", mm_pd$characteristics_ch1)
> mm_pd$HCV <- ifelse(mm_pd$HCV=="Neg", "-", "+")
> mm_pd$treatment <- gsub(".*: ", "", mm_pd$characteristics_ch1.2)

> HCV_matrix <-model.matrix(~0+HCV,mm_pd)
> colnames(HCV_matrix) <- c("Neg", "Pos")
> fit_HCV_matrix <- lmFit(mm_eset, HCV_matrix)
> ebay_HCV_matrix <- eBayes(fit_HCV_matrix)
> TopHCV <- topTable(ebay_HCV_matrix,adjust="BH")
> colnames(TopHCV)


> contrast_HCV <- makeContrasts(Neg-Pos,levels=HCV_matrix)
> sum(TopHCV$adj.P.Val < 0.1)
[1] 10
> HCV_eset_small <- mm_pd[TopHCV$adj.P.Val <0.1,]
> treatment_matrix <-model.matrix(~0+treatment,HCV_eset_small)

> fit_treatment_matrix <- lmFit(mm_eset, treatment_matrix)
> ebay_treatment_matrix <- eBayes(fit_treatment_matrix)
> TopHCV2 <- topTable(ebay_treatment_matrix,adjust="BH")
> colnames(TopHCV2)
     
> mm_pd$treatment

> colnames(treatment_matrix) <- c("Mock", "Poly")
> contrast_treatment <- makeContrasts(Mock-Poly,levels=treatment_matrix)
> sum(TopHCV2$adj.P.Val < 0.1)
[1] 10

> library(pheatmap)
> pheatmap(contrast_treatment)
