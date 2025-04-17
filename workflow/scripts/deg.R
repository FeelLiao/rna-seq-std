library(tidyverse)
library(edgeR)

count_matrix <-snakemake@input[[1]]
deg_table <- snakemake@output[["deg_table"]]
normalize_method <- snakemake@params[["normalization"]]
fdr_threshold <- snakemake@params[["fdr_threshold"]]
log_threshold <- snakemake@params[["log2fc_threshold"]]
# contrast_type <- snakemake@params[["contrast"]]

count <- read_csv(count_matrix)

# sample group
samples_p <- colnames(count)[-1]
condition_p <- str_split(
  string = samples_p,
  pattern = "-", simplify = TRUE
)[, 1]
group_p <- tibble(samples_p, condition_p)

# 位置效应 Position 差异基因分析
p_expr_edger <- DGEList(
  counts = count,
  group = condition_p
) ## 构建deglist对象
p_expr_edger <- p_expr_edger[
  filterByExpr(p_expr_edger,
    group = condition_p
  ), ,
  keep.lib.sizes = FALSE
] ## 过滤低表达基因
p_expr_deg_norm <- calcNormFactors(p_expr_edger, method = normalize_method)
design_p <- model.matrix(~ 0 + condition_p)
rownames(design_p) <- colnames(p_expr_deg_norm)
colnames(design_p) <- levels(factor(condition_p))
p_deg_disp <- estimateDisp(p_expr_deg_norm, design = design_p, robust = TRUE)
p_deg_fit <- glmQLFit(p_deg_disp, design_p)

contra <- combn(unique(condition_p),2,FUN = function(x) paste(x, collapse = "-"))
contrt_p <- makeContrasts(
  contrasts = contra,
  levels = colnames(design_p)
)
p_deg_lrt <- glmQLFTest(p_deg_fit, contrast = contrt_p)
p_deg_edger <- na.omit(topTags(p_deg_lrt, n = nrow(p_deg_lrt)))

p_deg_edger <- tibble(p_deg_edger$table)

output <- p_deg_edger |>
  filter(
    FDR < fdr_threshold,
    if_all(starts_with("logFC"), ~ abs(.x) > log_threshold)
  )

write_csv(output, deg_table)
