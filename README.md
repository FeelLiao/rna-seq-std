# Snakemake workflow: rna-seq-std

[English](./README.md) | [中文](./README_zh.md)

This workflow aims to provide a standard RNA-Seq analysis with both reliability and reproducibility. Using this workflow for your analysis, you need the [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html) software and the [Conda](https://www.anaconda.com/) environment management installed in your local machine.

Using a mirror for conda package download is highly recommended for some regions. In china, [Tsinghua mirror](https://mirrors.tuna.tsinghua.edu.cn/help/anaconda/) is suggested, especially in CERNET group. Two unofficial conda source is need in this project: conda-forge and bioconda. 

The workflow is divided into 3 sections:

- [x] Download SRA file from NCBI
- [x] RNA-Seq upstream analysis.
  - [x] `fastp` for raw reads quality control.
  - [x] `hisat2` for genome index and clean reads alignment.
  - [x] `featureCounts` for transcript quantification.
- [ ] RNA-Seq downstream analysis. 
  - [ ] `edgeR` for different expressed gene analysis.

⚠️Notes: This workflow is still under development. If you have any problem, please submit a issue. 

## Usage

1. Clone this repository

```bash
git clone 
```

2. Change `config/config.yaml` to your needs, read [config](config/README.md) for more information.

3. Use `snakemake` to run this workflow.

```bash
snakemake -c 20 --use-conda
```

  - -c: threads this workflow will use.

4. Result: results is in `out` directory.

## Features

- [ ] snakemake reports.
- [ ] Rmarkdown reports.

## Reference

- [kevinrue/snakemake_rnaseq_hisat2](https://github.com/kevinrue/snakemake_rnaseq_hisat2)
- [snakemake-wrappers-hisat2-index](https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/bio/hisat2/index.html)
- [snakemake-wrappers-hisat2-align](https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/bio/hisat2/align.html)
- [snakemake-wrappers-samtools-flagstat](https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/bio/samtools/flagstat.html)
- [snakemake-wrappers-fastp](https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/bio/fastp.html)
- [c-kroeger/snakemake-hisat2-stringtie](https://github.com/c-kroeger/snakemake-hisat2-stringtie)
