# model fitting using step 4 script
# source("workflow/4_fit_model.R")

# test for specific comparisons
contrast <- makeContrasts(HCC827 - H1975, levels = design)

res <- glmQLFTest(fit, contrast = contrast)
# examine number of DE genes
summary(decideTests(res))
# test for logFC threshold
tr <- glmTreat(fit, contrast = contrast, 
               lfc = log2(1.5))
is_de <- decideTests(tr)
summary(is_de)
# top DE results
topTags(tr)
# MD plot
# pdf("results/figure/MD.pdf", height = 4, width = 8)
plotMD(tr, main = "HCC827 vs H1975")
# dev.off()