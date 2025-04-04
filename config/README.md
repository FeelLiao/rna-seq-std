# General configuration

To configure this workflow, modify `config/config.yaml` according to your needs, following the explanations provided in the file.

# Sample setup (samples)

The sample file is specified via comma-separated files (`.csv`).

## samples

The default sample sheet is `config/samples.tsv` (as configured in `config/config.yaml`).

Each sample refers to an actual physical sample, and replicates (both biological and technical) are specified as separate samples.
For each sample, you will always have to specify a sample name by `sample`.

In addition, you need to specify the `group` of samples. When replicates are contained in an experiment, the file name of a sample usually contains the sample name and replicates number separated by "-". 

Two more columns are needed to be fulfilled, `read1` and `read2`, which represent the sample path. 

Finally, `extra` column is not necessary unless you have more complicated experiment design.

**Here, we propose that you have a folder containing raw reads data, and the filenames fit the formation `sample-replicates_R1/2.fastq`. Hence, we provide a python script `config/sample_pre.py` to process the sample information of rawdata in that situation. Before you use this workflow, just run `python config/sample_pre.py -h` to see how this script will help you in generating `sample_sheet.csv`.**

```bash
python config/sample_pre.py -h
# usage: sample_pre.py [-h] [-i INPUT] [-e EXTENSION] [-o OUTPUT]

# Automatic create sample sheet for this RNA-seq pipeline

# options:
#   -h, --help            show this help message and exit
#   -i INPUT, --input INPUT
#                         Path of your raw data [default: ../rawdata]
#   -e EXTENSION, --extension EXTENSION
#                         file extension of raw data [default: fastq]
#   -o OUTPUT, --output OUTPUT
#                         output path of sample sheet [default: .]

# version: 0.1.0
```

# Reference genome (ref)

The reference genome and genome annotation used for RNA-Seq

## genome

Define the reference genome to use. file can be `fasta` and `gz` format.

You could download this via NCBI or EMBL database.

## annotation

Which genome annotation to use. Usually, the more precise the annotation is, the more persuasive quantification results is. 

Both gff and gtf format are acceptable.

# Clean

If you want to clean the output, default is `true`. If you set this to `true`, the output will be removed after the workflow finished.
This is useful when you want to save disk space. Only the intermediate files will be removed, the final output will be kept.

# Reports

If you want to generate reports, default is `true`. If you set this to `true`, the report will be generated after the workflow finished.
This is useful when you want to check the quality of the workflow. The report will be generated in the `out/reports` directory.

# SRA download configuration (SRA)

To successfully run the SRA download , you need a stable connection to NCBI, or its a annoying time consuming task.

This workflow will download sra and decompress these files into fastq format automatically. This means that the output files are only the decompress fastq files, not contain sra original files.

After the sra workflow finished, you will get a sample file named `sample_sheet_sra.csv` in the same directory of `acc_list`. For downstream analysis, you need to fulfill the sample information yourself.

You can run this by targeting sra in snakemake commandline.

```bash
snakemake sra -c 30 --use-conda --conda-cleanup-pkgs
```

## activate

When you plan to use SRA download, you need to set this to `true`. If you set this to `false`, downstream analysis of these sra files will be unavailable.

**If your plan to use SRA for downstream analysis, please change `samples` in config file after you run sra target.**

## acc_list

Path to SRA accession list. Usually, you can download it using NCBI SRA run selector.

## output_dir

The directory that will store the processed sra files, usually in fastq format. 

# featureCounts configuration (featureCounts)

Specify the parameters that featureCounts use. For more information, see [official website](https://subread.sourceforge.net/featureCounts.html)

# New Gene

If you want to add new gene annotation, default is `false`. If you set this to `true`, the new gene annotation will be added to the existing gene annotation. 

## activate

When you plan to add new gene annotation, you need to set this to `true`.

## stringtie_params

The parameters that stringtie use. For more information, see [official website](https://ccb.jhu.edu/software/stringtie/index.shtml?t=example)