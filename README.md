[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15148691.svg)](https://doi.org/10.5281/zenodo.15148691)
[![GitHub license](https://img.shields.io/github/license/FeelLiao/rna-seq-std)](https://github.com/FeelLiao/rna-seq-std/blob/main/LICENSE)
[![Tests](https://github.com/FeelLiao/rna-seq-std/actions/workflows/test.yaml/badge.svg)](https://github.com/FeelLiao/rna-seq-std/actions/workflows/test.yaml)
![GitHub Release](https://img.shields.io/github/v/release/FeelLiao/rna-seq-std)
[![Snakemake](https://img.shields.io/badge/Snakemake->=8.25.3-green)](https://snakemake.readthedocs.io/en/stable/)

# RNA-Seq Analysis Pipeline: A Scalable Snakemake Workflow for Transcriptome Profiling

[English](./README.md) | [中文](./README_zh.md)

A Snakemake workflow for trimming, alignment, quantification and differential expressed gene analysis for RNA-Seq experiments, starting from paired-end reads with SRA or raw sequence data files, including a comprehensive report generated by Rmarkdown.

This workflow aims to provide a comprehensive, reproducible, and user-friendly workflow for RNA-Seq data analysis from raw reads to biological insights. Using this workflow for your analysis, you need the [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html) software and the [Conda](https://www.anaconda.com/) environment management installed in your local machine.

> [!NOTE]  
> Using a mirror for conda package download is highly recommended for some regions. In china, [Tsinghua mirror](https://mirrors.tuna.tsinghua.edu.cn/help/anaconda/) is suggested, especially in CERNET group. 

**If you use this workflow in a publication, please don't forget to give credit to the authors by citing it using this DOI [10.5281/zenodo.15148691](https://doi.org/10.5281/zenodo.15148691)**

# 🌟 Features

This workflow provides an end-to-end solution for RNA-Seq analysis with distinctive advantages:

1. **Comprehensive Analysis Pipeline**
- Raw data processing: Automatic SRA download (via NCBI tools) or direct FASTQ input  
- Quality control: fastp with multiQC integration
- Alignment: HISAT2 mapping with genomic index auto-generation  
- Quantification: FeatureCounts for transcript-level counts  
- Differential Expression: edgeR integration (under development)  
- Novel Analysis: De novo transcript assembly (StringTie) and lncRNA prediction pipeline
 
2. **Production-Grade Reproducibility**
- Conda-based environment management for version-controlled dependencies  
- Parameter centralization via YAML configuration  
 
3. **Advanced Reporting**
- Interactive HTML reports via RMarkdown
- Sample metrics visualization (PCA, sample clustering, MA plots)

> [!Important]
> This workflow is still under development. If you have any problem, please submit a issue. 

# 🔭 Roadmap  
 
Implemented Features  
- [x] Automated SRA → FASTQ conversion  
- [x] MultiQC-integrated QC pipeline  
- [x] Transcript quantification consensus  
- [x] Novel isoform detection (StringTie)  
- [x] RMarkdown report generation  
 
Upcoming Features  
- [x] edgeR differential expression module  
- [ ] WGCNA co-expression network analysis  
- [ ] LncRNA classification pipeline 

# 🚀 Quick Start  

1. Clone this repository

```bash
git clone https://github.com/FeelLiao/rna-seq-std.git
```

2. Change `config/config.yaml` to your needs, read [config](config/README.md) for more information.

3. Use `snakemake` to run this workflow.

```bash
snakemake <target> -c 20 --use-conda --conda-cleanup-pkgs
```

4. Result: results is in `out` directory.

# 📑 Reference

- [kevinrue/snakemake_rnaseq_hisat2](https://github.com/kevinrue/snakemake_rnaseq_hisat2)
- [snakemake-wrappers-hisat2-index](https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/bio/hisat2/index.html)
- [snakemake-wrappers-hisat2-align](https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/bio/hisat2/align.html)
- [snakemake-wrappers-samtools-flagstat](https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/bio/samtools/flagstat.html)
- [snakemake-wrappers-fastp](https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/bio/fastp.html)
- [c-kroeger/snakemake-hisat2-stringtie](https://github.com/c-kroeger/snakemake-hisat2-stringtie)
 
 

 


