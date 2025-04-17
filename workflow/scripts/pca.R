library(tidyverse)
library(FactoMineR)
library(ggrepel)

sample_file <- snakemake@input[[1]]
pca_pdf <- snakemake@output[["pdf"]]
pca_png <- snakemake@output[["png"]]
pca_eig_s <- snakemake@output[["eig"]]
pca_results <- snakemake@output[["result"]]
ncps <- snakemake@params[["ncp"]]

gene_tpm <- read_csv(sample_file) |>
  mutate(valid_samples = rowSums(across(-gene_id, ~ .x > 1)) >= 5) |>
  filter(valid_samples) |>
  select(-valid_samples) |>
  column_to_rownames("gene_id") |>
  t()

# PCA Calculation
pca_out <- PCA(gene_tpm, ncp = ncps, scale.unit = TRUE, graph = FALSE)

pca_eig_out <- pca_out$eig |>
  as.data.frame() |>
  rownames_to_column("PC") |>
  rename(
    Variance = `percentage of variance`,
    Cumulative = `cumulative percentage of variance`
  )

pca_sample <- data.frame(pca_out$ind$coord) |>
  rownames_to_column("Sample")
pca_eig1 <- round(pca_out$eig[1, 2], 2)
pca_eig2 <- round(pca_out$eig[2, 2], 2)
pca_sample <- pca_sample |>
  mutate(group = str_split_i(Sample, "-", 1))


p <- ggplot(data = pca_sample, aes(x = Dim.1, y = Dim.2)) +
  geom_point(aes(color = group), size = 5, alpha = 0.4) +
	geom_text_repel(
      aes(label = Sample),
      box.padding = 0.5, # 标签与点间距
      segment.color = "grey50", # 连接线颜色
      max.overlaps = 20, # 最大允许重叠次数
      min.segment.length = 0.2 # 最短连接线长度
    ) +
  theme_bw() +
  labs(
    x = paste("PCA1:", pca_eig1, "%"),
    y = paste("PCA2:", pca_eig2, "%"), color = ""
  )

ggsave(
  filename = pca_pdf,
  plot = p,
  width = 8,
  height = 6,
  device = "pdf",
)

ggsave(
  filename = pca_png,
  plot = p,
  width = 8,
  height = 6,
  dpi = 600,
  device = "png",
)

write_csv(pca_eig_out, pca_eig_s, col_names = TRUE)
write_csv(pca_sample, pca_results, col_names = TRUE)
