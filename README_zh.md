[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15148691.svg)](https://doi.org/10.5281/zenodo.15148691)
[![GitHub license](https://img.shields.io/github/license/FeelLiao/rna-seq-std)](https://github.com/FeelLiao/rna-seq-std/blob/main/LICENSE)
[![Tests](https://github.com/FeelLiao/rna-seq-std/actions/workflows/test.yaml/badge.svg)](https://github.com/FeelLiao/rna-seq-std/actions/workflows/test.yaml)
![GitHub Release](https://img.shields.io/github/v/release/FeelLiao/rna-seq-std)
[![Snakemake](https://img.shields.io/badge/Snakemake->=8.25.3-green)](https://snakemake.readthedocs.io/en/stable/)
 
# RNA-Seq分析流程：基于Snakemake的可扩展转录组分析工作流 
 
 
一个用于RNA-Seq实验的Snakemake工作流，涵盖从SRA或双端测序数据的原始序列文件开始的质控、比对、定量及差异表达基因分析，并支持通过Rmarkdown生成综合性分析报告。

本工作流旨在为RNA-Seq数据分析提供完整、可复现且用户友好的解决方案。使用前需确保本地环境已安装[Snakemake](https://snakemake.readthedocs.io/en/stable/index.html)和[Conda](https://www.anaconda.com/)包管理系统。
 
> [!NOTE]  
> 建议部分地区用户配置Conda镜像加速下载。中国境内用户推荐使用[清华源](https://mirrors.tuna.tsinghua.edu.cn/help/anaconda/)，尤其是教育网用户。
 
**若在学术出版物中使用本工作流，请通过DOI引用以注明来源：[10.5281/zenodo.15148691](https://doi.org/10.5281/zenodo.15148691)**
 
# 🌟 核心特性 
 
本工作流提供端到端RNA-Seq分析解决方案，具有以下显著优势：
 
1. 全流程分析体系
- 原始数据处理：支持SRA自动下载（通过NCBI工具）或直接FASTQ输入  
- 质量控制：集成fastp的MultiQC的质控模块 
- 序列比对：基于HISAT2的基因组自动索引构建与比对  
- 定量分析：FeatureCounts转录本计数  
- 差异表达：整合edgeR分析模块（开发中）  
- 创新分析：de novo转录本组装（StringTie）及lncRNA预测流程 
 
2. 工业级可复现性
- Conda环境管理确保依赖版本可控  
- 通过YAML配置文件集中管理参数  
 
3. 智能化报告系统
- 基于RMarkdown的交互式HTML报告 
- 样本指标可视化（PCA分析、样本聚类、MA图）
 
> [!Important]
> 本项目仍处于持续开发阶段，如有问题请提交issue反馈。
 
# 🔭 发展路线  
 
已实现功能  
- [x] SRA→FASTQ自动化转换  
- [x] 集成MultiQC的质控流程  
- [x] 转录本定量
- [x] 新转录本检测（StringTie）  
- [x] RMarkdown报告生成  
 
规划功能  
- [ ] edgeR差异基因分析 
- [ ] WGCNA共表达网络分析  
- [ ] lncRNA鉴定流程
 
# 🚀 快速入门  
 
1. 克隆仓库 
 
```bash 
git clone https://github.com/FeelLiao/rna-seq-std.git 
```
 
2. 根据需求修改`config/config.yaml`配置文件，详见[config](config/README.md)
 
3. 运行Snakemake工作流 
 
```bash 
snakemake <target> -c 20 --use-conda --conda-cleanup-pkgs 
```
 
4. 结果文件将输出至`out`目录 
 
# 📑 参考资料 

- [kevinrue/snakemake_rnaseq_hisat2](https://github.com/kevinrue/snakemake_rnaseq_hisat2)
- [snakemake-wrappers-hisat2-index](https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/bio/hisat2/index.html)
- [snakemake-wrappers-hisat2-align](https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/bio/hisat2/align.html)
- [snakemake-wrappers-samtools-flagstat](https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/bio/samtools/flagstat.html)
- [snakemake-wrappers-fastp](https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/bio/fastp.html)
- [c-kroeger/snakemake-hisat2-stringtie](https://github.com/c-kroeger/snakemake-hisat2-stringtie)