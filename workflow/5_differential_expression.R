# model fitting using step 4 script
# source("workflow/4_fit_model.R")

# test for specific comparisons
contrast <- makeContrasts(HCC827 - H1975, levels = design)

res <- glmQLFTest(fit, contrast = contrast)

# test for DE
is_de <- decideTests(res)
summary(is_de)

# top DE results
topTags(res)

# Figure 4a
# pdf("results/figure/MD.pdf", height = 5, width = 8)
plotMD(res, main = "HCC827 vs H1975", cex = 0.3)
# dev.off()