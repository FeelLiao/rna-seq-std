
# Snakemake工作流：rna-seq-std 

本工作流旨在提供一个兼具可靠性和可重复性的标准RNA-Seq分析流程。若要使用此工作流进行分析，你需要在本地机器上安装 [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html) 软件和 [Conda](https://www.anaconda.com/) 环境管理工具。 
 
强烈建议某些地区的用户使用镜像源来下载Conda软件包。在中国，推荐使用 [清华大学源](https://mirrors.tuna.tsinghua.edu.cn/help/anaconda/) ，尤其是在教育网中的用户。本项目需要两个非官方的Conda源：conda-forge 和 bioconda 。 
 
该工作流分为三个部分： 
- [x] 从NCBI下载SRA文件 
- [x] RNA-Seq上游分析 
  - [x] 使用 `fastp` 对原始 reads 进行质量控制 
  - [x] 使用 `hisat2` 构建基因组索引并对 clean reads 进行比对 
  - [x] 使用 `featureCounts` 进行转录本定量 
- [ ] RNA-Seq下游分析 
  - [ ] 使用 `edgeR` 进行差异表达基因分析
 
⚠️注意：该工作流仍在开发中，如果遇到任何问题，欢迎提交 issues 
 
## 使用方法

1. 克隆此仓库

```bash 
git clone 
```

2. 根据你的需求修改 `config/config.yaml` 文件，更多信息请阅读 [config](config/README.md) 。 

3. 使用 `snakemake` 运行此工作流。 

```bash 
snakemake -c 20 --use-conda 
``` 

  - -c：该工作流将使用的线程数。 

4. 结果：分析结果将存储在 `out` 目录中。 
 
## 特性

- ❌ snakemake 报告 
- [ ] Rmarkdown 报告 
 
## 参考资料 

- [kevinrue/snakemake_rnaseq_hisat2](https://github.com/kevinrue/snakemake_rnaseq_hisat2)
- [snakemake-wrappers-hisat2-index](https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/bio/hisat2/index.html)
- [snakemake-wrappers-hisat2-align](https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/bio/hisat2/align.html)
- [snakemake-wrappers-samtools-flagstat](https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/bio/samtools/flagstat.html)
- [snakemake-wrappers-fastp](https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/bio/fastp.html)
- [c-kroeger/snakemake-hisat2-stringtie](https://github.com/c-kroeger/snakemake-hisat2-stringtie)